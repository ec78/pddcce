library(haven)
library(pwt10)

d <- read_dta("examples/penn_sample.dta")
data("pwt10.0")
pwt <- pwt10.0

our_years <- sort(unique(d$year))
ids <- sort(unique(d$id))

# Compute all variables for matching:
# 1. log_rgdpo_pc = log(rgdpo/pop)  [per capita GDP]
# 2. log_hc                          [human capital - same scale across versions]
# We'll match on variance of difference for log_rgdpo_pc,
# plus squared distance for log_hc (scale-invariant across versions)

pwt$log_rgdpo_pc <- log(pwt$rgdpo / pwt$pop)
pwt$log_hc <- log(pwt$hc)

# Balanced countries
check_cols <- c("log_rgdpo_pc", "log_hc")
pwt_sub <- pwt[pwt$year %in% our_years, ]
balanced_iso <- NULL
for (iso in unique(as.character(pwt_sub$isocode))) {
  sub <- pwt_sub[as.character(pwt_sub$isocode)==iso, check_cols]
  if (nrow(sub) == length(our_years) && all(!is.na(sub))) {
    balanced_iso <- c(balanced_iso, iso)
  }
}
cat("Balanced PWT countries:", length(balanced_iso), "\n")

pwt_full <- pwt_sub[as.character(pwt_sub$isocode) %in% balanced_iso, ]
pwt_iso <- sort(balanced_iso)
n_pwt <- length(pwt_iso)

# Build time series matrices
build_our_mat <- function(varname) {
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

build_pwt_mat <- function(varname) {
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

cat("Building our matrices...\n")
our_gdp  <- build_our_mat("log_rgdpo")
our_hc   <- build_our_mat("log_hc")

cat("Building PWT matrices...\n")
pwt_gdp  <- build_pwt_mat("log_rgdpo_pc")
pwt_hc   <- build_pwt_mat("log_hc")

# Combined score: var(our_gdp - pwt_gdp) + 10*mean((our_hc - pwt_hc)^2)
# HC gets higher weight since it's on a smaller scale and more discriminating
all_scores <- matrix(Inf, 93, n_pwt)
for (ii in seq_along(ids)) {
  g_our <- our_gdp[ii, ]
  h_our <- our_hc[ii, ]
  if (any(is.na(g_our)) || any(is.na(h_our))) next
  for (ci in seq_along(pwt_iso)) {
    g_pwt <- pwt_gdp[ci, ]
    h_pwt <- pwt_hc[ci, ]
    if (any(is.na(g_pwt)) || any(is.na(h_pwt))) next
    score_gdp <- var(g_our - g_pwt)
    score_hc  <- mean((h_our - h_pwt)^2)
    all_scores[ii, ci] <- score_gdp + 5 * score_hc
  }
}

# Greedy 1-to-1 assignment: process by ascending score
flat_scores <- data.frame(
  ii = rep(seq_along(ids), n_pwt),
  ci = rep(seq_len(n_pwt), each=length(ids)),
  score = as.vector(all_scores)
)
flat_scores <- flat_scores[order(flat_scores$score), ]
flat_scores <- flat_scores[is.finite(flat_scores$score), ]

assigned_id <- rep(FALSE, length(ids))
assigned_pwt <- rep(FALSE, n_pwt)
matches <- data.frame(id=ids, isocode=NA_character_, dist=Inf, stringsAsFactors=FALSE)

for (row in seq_len(nrow(flat_scores))) {
  ii <- flat_scores$ii[row]
  ci <- flat_scores$ci[row]
  if (!assigned_id[ii] && !assigned_pwt[ci]) {
    matches$isocode[ii] <- pwt_iso[ci]
    matches$dist[ii] <- flat_scores$score[row]
    assigned_id[ii] <- TRUE
    assigned_pwt[ci] <- TRUE
    if (all(assigned_id)) break
  }
}

pwt_names <- unique(data.frame(
  isocode=as.character(pwt$isocode),
  country=as.character(pwt$country),
  stringsAsFactors=FALSE
))
matches <- merge(matches, pwt_names, by="isocode", all.x=TRUE)
matches <- matches[order(matches$id),]

cat("\nFull mapping (combined GDP+HC, greedy 1-to-1):\n")
print(matches[, c("id","isocode","country","dist")])

cat("\nMax score:", max(matches$dist, na.rm=TRUE), "\n")
cat("Mean score:", mean(matches$dist, na.rm=TRUE), "\n")
cat("Near-zero (< 0.001):", sum(matches$dist < 0.001, na.rm=TRUE), "\n")
cat("Good (< 0.01):", sum(matches$dist < 0.01, na.rm=TRUE), "\n")
cat("High (> 0.1):", sum(matches$dist > 0.1, na.rm=TRUE), "\n")

# No duplicates by construction
write.csv(matches[order(matches$id), c("id","isocode","country")],
          "validation/id_country_map.csv", row.names=FALSE)
cat("Saved to id_country_map.csv\n")
