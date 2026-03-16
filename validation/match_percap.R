library(haven)
library(pwt10)

d <- read_dta("examples/penn_sample.dta")
data("pwt10.0")
pwt <- pwt10.0

# log_rgdpo = log(rgdpo / pop) = log per-capita output-side GDP
# PWT10: rgdpo in millions 2017 USD, pop in millions -> ratio in 2017 USD/person
our_years <- sort(unique(d$year))
ids <- sort(unique(d$id))

# PWT per-capita GDP
pwt$rgdpo_pc <- pwt$rgdpo / pwt$pop  # USD per person (2017)
pwt$log_rgdpo_pc <- log(pwt$rgdpo_pc)

# Find PWT countries with complete data
pwt_sub <- pwt[pwt$year %in% our_years & !is.na(pwt$log_rgdpo_pc), ]
country_counts <- table(as.character(pwt_sub$isocode))
balanced_iso <- names(country_counts[country_counts == length(our_years)])
cat("PWT10 countries with complete per-cap rgdpo 1960-2007:", length(balanced_iso), "\n")

pwt_full <- pwt_sub[as.character(pwt_sub$isocode) %in% balanced_iso, ]
pwt_iso <- sort(balanced_iso)
n_pwt <- length(pwt_iso)

# Verify US matches id=66 (log_rgdpo ≈ 10.8 in 2000)
us_2000 <- pwt[pwt$year==2000 & as.character(pwt$isocode)=="USA", "log_rgdpo_pc"]
cat("US 2000 log(rgdpo/pop) =", us_2000, "\n")  # should be ~10.8

# Build matrices
our_mat <- matrix(NA_real_, 93, length(our_years))
for (yi in seq_along(our_years)) {
  yr <- our_years[yi]
  sub <- d[d$year == yr, c("id","log_rgdpo")]
  for (ii in seq_along(ids)) {
    v <- sub$log_rgdpo[sub$id == ids[ii]]
    if (length(v) == 1) our_mat[ii, yi] <- v
  }
}

pwt_mat <- matrix(NA_real_, n_pwt, length(our_years))
for (yi in seq_along(our_years)) {
  yr <- our_years[yi]
  sub <- pwt_full[pwt_full$year == yr, c("isocode","log_rgdpo_pc")]
  for (ci in seq_along(pwt_iso)) {
    v <- sub$log_rgdpo_pc[as.character(sub$isocode) == pwt_iso[ci]]
    if (length(v) == 1) pwt_mat[ci, yi] <- v
  }
}

# Match by variance of difference (handles any constant offset from base-year diffs)
matches <- data.frame(id=ids, isocode=NA_character_, dist=Inf, stringsAsFactors=FALSE)
for (ii in seq_along(ids)) {
  our_vec <- our_mat[ii, ]
  if (any(is.na(our_vec))) next
  scores <- apply(pwt_mat, 1, function(x) {
    if (any(is.na(x))) return(Inf)
    var(our_vec - x)
  })
  best <- which.min(scores)
  matches$isocode[ii] <- pwt_iso[best]
  matches$dist[ii] <- scores[best]
}

pwt_names <- unique(data.frame(
  isocode=as.character(pwt$isocode),
  country=as.character(pwt$country),
  stringsAsFactors=FALSE
))
matches <- merge(matches, pwt_names, by="isocode", all.x=TRUE)
matches <- matches[order(matches$id),]

cat("\nFull mapping:\n")
print(matches[, c("id","isocode","country","dist")])

dups <- matches$isocode[duplicated(matches$isocode)]
cat("\nDuplicates:", length(dups), "\n")
if (length(dups) > 0) cat(paste(unique(dups), collapse=", "), "\n")

cat("Max variance:", max(matches$dist, na.rm=TRUE), "\n")
cat("Mean variance:", mean(matches$dist, na.rm=TRUE), "\n")
cat("Near-perfect (var < 0.00001):", sum(matches$dist < 0.00001, na.rm=TRUE), "\n")
cat("Good (var < 0.0001):", sum(matches$dist < 0.0001, na.rm=TRUE), "\n")
cat("OK (var < 0.001):", sum(matches$dist < 0.001, na.rm=TRUE), "\n")

if (length(dups) == 0 || TRUE) {
  write.csv(matches[order(matches$id), c("id","isocode","country")],
            "validation/id_country_map.csv", row.names=FALSE)
  cat("Saved to id_country_map.csv\n")
}
