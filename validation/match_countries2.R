library(haven)
library(pwt10)

d <- read_dta("examples/penn_sample.dta")
data("pwt10.0")
pwt <- pwt10.0

# Normalize log_rgdpo relative to max within each year (to remove base-year effect)
# This gives the relative GDP distribution which should match across PWT versions

our_years <- sort(unique(d$year))

# Build wide matrix: rows=country (by id), cols=years
# Each cell = log_rgdpo - mean(log_rgdpo) within that year
ids <- sort(unique(d$id))
n_ids <- length(ids)
n_yr  <- length(our_years)
our_mat <- matrix(NA_real_, n_ids, n_yr)

for (yi in seq_along(our_years)) {
  yr <- our_years[yi]
  sub <- d[d$year == yr, c("id", "log_rgdpo")]
  sub <- sub[order(sub$id),]
  yr_mean <- mean(sub$log_rgdpo, na.rm=TRUE)
  for (ii in seq_along(ids)) {
    v <- sub$log_rgdpo[sub$id == ids[ii]]
    if (length(v) > 0) our_mat[ii, yi] <- v - yr_mean
  }
}

# Build PWT wide matrix for countries with data in 1960-2007
pwt_sub <- pwt[pwt$year %in% our_years & !is.na(pwt$rgdpo),]
pwt_sub$log_rgdpo <- log(pwt_sub$rgdpo)
pwt_countries <- unique(pwt_sub$isocode)
n_pwt <- length(pwt_countries)
pwt_mat <- matrix(NA_real_, n_pwt, n_yr)

for (yi in seq_along(our_years)) {
  yr <- our_years[yi]
  sub <- pwt_sub[pwt_sub$year == yr, c("isocode", "log_rgdpo")]
  yr_mean <- mean(sub$log_rgdpo, na.rm=TRUE)
  for (ci in seq_along(pwt_countries)) {
    v <- sub$log_rgdpo[sub$isocode == pwt_countries[ci]]
    if (length(v) > 0) pwt_mat[ci, yi] <- v - yr_mean
  }
}

# Only use years where all countries have data
valid_cols_our <- which(colSums(!is.na(our_mat)) == n_ids)
cat("Years with complete our data:", length(valid_cols_our), "\n")

# For each of our 93 ids, find best matching PWT country (Euclidean on demeaned log GDP)
matches <- data.frame(id=ids, isocode=NA_character_, dist=Inf,
                      stringsAsFactors=FALSE)

for (ii in seq_along(ids)) {
  our_vec <- our_mat[ii, valid_cols_our]
  if (any(is.na(our_vec))) next

  pwt_sub_mat <- pwt_mat[, valid_cols_our, drop=FALSE]
  dists <- apply(pwt_sub_mat, 1, function(x) {
    if (any(is.na(x))) return(Inf)
    sqrt(sum((x - our_vec)^2))
  })
  best <- which.min(dists)
  matches$isocode[ii] <- pwt_countries[best]
  matches$dist[ii] <- dists[best]
}

# Add country names
pwt_names <- unique(pwt[, c("isocode","country")])
matches <- merge(matches, pwt_names, by="isocode", all.x=TRUE)
matches <- matches[order(matches$id),]

cat("\nMapping (id -> country):\n")
print(matches[, c("id","isocode","country","dist")])

cat("\nMax dist:", max(matches$dist, na.rm=TRUE))
cat("\nMean dist:", mean(matches$dist, na.rm=TRUE))

dup <- matches$isocode[duplicated(matches$isocode)]
if (length(dup) > 0) {
  cat("\n\nDUPLICATE matches:", paste(unique(dup), collapse=", "), "\n")
  cat("Affected rows:\n")
  print(matches[matches$isocode %in% dup, c("id","isocode","country","dist")])
}

write.csv(matches[order(matches$id), c("id","isocode","country")],
          "validation/id_country_map.csv", row.names=FALSE)
cat("\nSaved mapping to validation/id_country_map.csv\n")
