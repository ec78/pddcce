library(haven)
library(pwt10)

d <- read_dta("examples/penn_sample.dta")
data("pwt10.01")
pwt <- pwt10.01

our_years <- sort(unique(d$year))  # 1960-2007
ids <- sort(unique(d$id))

cat("PWT10.01 dims:", nrow(pwt), "x", ncol(pwt), "\n")

# Compute per-capita output GDP
pwt$log_rgdpo_pc <- log(pwt$rgdpo / pwt$pop)
pwt$log_hc <- log(pwt$hc)

# Find balanced countries for all our years
pwt_sub <- pwt[pwt$year %in% our_years, ]
balanced_iso <- NULL
for (iso in unique(as.character(pwt_sub$isocode))) {
  sub <- pwt_sub[as.character(pwt_sub$isocode)==iso, c("log_rgdpo_pc","log_hc")]
  if (nrow(sub) == length(our_years) && all(!is.na(sub))) {
    balanced_iso <- c(balanced_iso, iso)
  }
}
cat("PWT10.01 balanced (rgdpo_pc + hc) for 1960-2007:", length(balanced_iso), "\n")

pwt_full <- pwt_sub[as.character(pwt_sub$isocode) %in% balanced_iso, ]
pwt_iso <- sort(balanced_iso)
n_pwt <- length(pwt_iso)

# Build matrices
build_our <- function(varname) {
  mat <- matrix(NA_real_, 93, length(our_years))
  for (yi in seq_along(our_years)) {
    sub <- d[d$year == our_years[yi], c("id", varname)]
    for (ii in seq_along(ids)) {
      v <- sub[[varname]][sub$id == ids[ii]]
      if (length(v)==1) mat[ii, yi] <- v
    }
  }
  mat
}
build_pwt <- function(varname) {
  mat <- matrix(NA_real_, n_pwt, length(our_years))
  for (yi in seq_along(our_years)) {
    sub <- pwt_full[pwt_full$year == our_years[yi], c("isocode", varname)]
    for (ci in seq_along(pwt_iso)) {
      v <- sub[[varname]][as.character(sub$isocode) == pwt_iso[ci]]
      if (length(v)==1) mat[ci, yi] <- v
    }
  }
  mat
}

cat("Building matrices...\n")
our_gdp <- build_our("log_rgdpo")
our_hc  <- build_our("log_hc")
pwt_gdp <- build_pwt("log_rgdpo_pc")
pwt_hc  <- build_pwt("log_hc")

# Combined score
cat("Computing pairwise scores...\n")
all_scores <- matrix(Inf, 93, n_pwt)
for (ii in seq_along(ids)) {
  g <- our_gdp[ii,]; h <- our_hc[ii,]
  if (any(is.na(g)) || any(is.na(h))) next
  for (ci in seq_along(pwt_iso)) {
    gp <- pwt_gdp[ci,]; hp <- pwt_hc[ci,]
    if (any(is.na(gp))||any(is.na(hp))) next
    s_gdp <- var(g - gp)   # base-year offset -> use var of difference
    s_hc  <- mean((h - hp)^2)
    all_scores[ii, ci] <- s_gdp + 5 * s_hc
  }
}

# Greedy 1-to-1 assignment
flat <- data.frame(
  ii=rep(seq_along(ids), n_pwt),
  ci=rep(seq_len(n_pwt), each=length(ids)),
  score=as.vector(all_scores)
)
flat <- flat[order(flat$score),]
flat <- flat[is.finite(flat$score),]

assigned_id  <- rep(FALSE, length(ids))
assigned_pwt <- rep(FALSE, n_pwt)
matches <- data.frame(id=ids, isocode=NA_character_, dist=Inf, stringsAsFactors=FALSE)

for (row in seq_len(nrow(flat))) {
  ii <- flat$ii[row]; ci <- flat$ci[row]
  if (!assigned_id[ii] && !assigned_pwt[ci]) {
    matches$isocode[ii] <- pwt_iso[ci]
    matches$dist[ii]    <- flat$score[row]
    assigned_id[ii]     <- TRUE
    assigned_pwt[ci]    <- TRUE
    if (all(assigned_id)) break
  }
}

# Add country names
pwt_names <- unique(data.frame(
  isocode=as.character(pwt$isocode),
  country=as.character(pwt$country),
  stringsAsFactors=FALSE
))
matches <- merge(matches, pwt_names, by="isocode", all.x=TRUE)
matches <- matches[order(matches$id),]

cat("\nFull mapping (PWT10.01, GDP+HC greedy):\n")
print(matches[, c("id","isocode","country","dist")])
cat("\nMax:", max(matches$dist), "Mean:", mean(matches$dist), "\n")
cat("< 0.001:", sum(matches$dist < 0.001), "\n")
cat("< 0.01:", sum(matches$dist < 0.01), "\n")
cat("< 0.05:", sum(matches$dist < 0.05), "\n")
cat("> 0.1:", sum(matches$dist > 0.1), "\n")

# Save
write.csv(matches[order(matches$id), c("id","isocode","country")],
          "validation/id_country_map_final.csv", row.names=FALSE)
cat("Saved to id_country_map_final.csv\n")
