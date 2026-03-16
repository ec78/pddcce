library(haven)
library(pwt10)

d <- read_dta("examples/penn_sample.dta")
data("pwt10.01")
pwt <- pwt10.01

our_years <- sort(unique(d$year))
ids <- sort(unique(d$id))

pwt$log_rgdpo_pc <- log(pwt$rgdpo / pwt$pop)

# Find ALL PWT countries with complete rgdpo_pc for 1960-2007
pwt_sub <- pwt[pwt$year %in% our_years, ]
iso_counts <- tapply(!is.na(pwt_sub$log_rgdpo_pc), as.character(pwt_sub$isocode), sum)
balanced_iso <- names(iso_counts[iso_counts == length(our_years)])
cat("PWT10.01 countries with complete log_rgdpo_pc 1960-2007:", length(balanced_iso), "\n")

pwt_full <- pwt_sub[as.character(pwt_sub$isocode) %in% balanced_iso, ]
pwt_iso <- sort(balanced_iso)
n_pwt <- length(pwt_iso)

# Build matrices
our_gdp <- matrix(NA_real_, 93, length(our_years))
for (yi in seq_along(our_years)) {
  sub <- d[d$year == our_years[yi], c("id","log_rgdpo")]
  for (ii in seq_along(ids)) {
    v <- sub$log_rgdpo[sub$id == ids[ii]]
    if (length(v)==1) our_gdp[ii, yi] <- v
  }
}

pwt_gdp <- matrix(NA_real_, n_pwt, length(our_years))
for (yi in seq_along(our_years)) {
  sub <- pwt_full[pwt_full$year == our_years[yi], c("isocode","log_rgdpo_pc")]
  for (ci in seq_along(pwt_iso)) {
    v <- sub$log_rgdpo_pc[as.character(sub$isocode) == pwt_iso[ci]]
    if (length(v)==1) pwt_gdp[ci, yi] <- v
  }
}

# GDP-only scores (variance of difference)
cat("Computing pairwise GDP-only scores...\n")
all_scores <- matrix(Inf, 93, n_pwt)
for (ii in seq_along(ids)) {
  g <- our_gdp[ii, ]
  if (any(is.na(g))) next
  for (ci in seq_along(pwt_iso)) {
    gp <- pwt_gdp[ci, ]
    if (any(is.na(gp))) next
    all_scores[ii, ci] <- var(g - gp)
  }
}

# Greedy 1-to-1
flat <- data.frame(
  ii=rep(seq_along(ids), n_pwt),
  ci=rep(seq_len(n_pwt), each=length(ids)),
  score=as.vector(all_scores)
)
flat <- flat[order(flat$score),]
flat <- flat[is.finite(flat$score),]

a_id  <- rep(FALSE, length(ids))
a_pwt <- rep(FALSE, n_pwt)
matches <- data.frame(id=ids, isocode=NA_character_, dist=Inf, stringsAsFactors=FALSE)

for (row in seq_len(nrow(flat))) {
  ii <- flat$ii[row]; ci <- flat$ci[row]
  if (!a_id[ii] && !a_pwt[ci]) {
    matches$isocode[ii] <- pwt_iso[ci]
    matches$dist[ii]    <- flat$score[row]
    a_id[ii]  <- TRUE
    a_pwt[ci] <- TRUE
    if (all(a_id)) break
  }
}

pwt_names <- unique(data.frame(
  isocode=as.character(pwt$isocode),
  country=as.character(pwt$country),
  stringsAsFactors=FALSE
))
matches <- merge(matches, pwt_names, by="isocode", all.x=TRUE)
matches <- matches[order(matches$id),]

cat("\nFull mapping (GDP-only greedy, all PWT countries):\n")
print(matches[, c("id","isocode","country","dist")])
cat("\nMax:", max(matches$dist), "Mean:", mean(matches$dist), "\n")
cat("< 0.001:", sum(matches$dist < 0.001), "\n")
cat("< 0.01:", sum(matches$dist < 0.01), "\n")
cat("< 0.05:", sum(matches$dist < 0.05), "\n")
cat("> 0.1:", sum(matches$dist > 0.1), "\n")

write.csv(matches[order(matches$id), c("id","isocode","country")],
          "validation/id_country_map_gdponly.csv", row.names=FALSE)
cat("Saved to id_country_map_gdponly.csv\n")
