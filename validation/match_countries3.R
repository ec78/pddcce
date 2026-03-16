library(haven)
library(pwt10)

d <- read_dta("examples/penn_sample.dta")
data("pwt10.0")
pwt <- pwt10.0

# Check isocode type
cat("isocode class:", class(pwt$isocode), "\n")
cat("Sample isocodes:", paste(head(unique(as.character(pwt$isocode)), 5), collapse=", "), "\n")

# Strategy: match by rank within year
# Countries that are consistently rank X across years => same country
# Use multiple years and compute average rank

our_years <- 1965:2000
ids <- sort(unique(d$id))

# Compute year-by-year rank of log_rgdpo among our 93 countries
our_ranks <- matrix(NA_real_, length(ids), length(our_years))
for (yi in seq_along(our_years)) {
  yr <- our_years[yi]
  sub <- d[d$year == yr, c("id","log_rgdpo")]
  sub <- sub[order(sub$id),]
  r <- rank(sub$log_rgdpo)
  our_ranks[, yi] <- r
}

# For PWT, compute ranks among ALL countries with data in each year (not just our 93)
# Instead use only the countries that appear in ALL our_years
pwt_sub <- pwt[pwt$year %in% our_years & !is.na(pwt$rgdpo), ]
pwt_always <- names(which(table(pwt_sub$isocode) == length(our_years)))
cat("PWT countries present in all", length(our_years), "years:", length(pwt_always), "\n")

pwt_sub <- pwt_sub[as.character(pwt_sub$isocode) %in% pwt_always, ]
pwt_iso <- sort(unique(as.character(pwt_sub$isocode)))
n_pwt <- length(pwt_iso)

pwt_ranks <- matrix(NA_real_, n_pwt, length(our_years))
for (yi in seq_along(our_years)) {
  yr <- our_years[yi]
  sub <- pwt_sub[pwt_sub$year == yr, c("isocode","rgdpo")]
  sub <- sub[order(sub$isocode),]
  # Match the ordering of pwt_iso
  idx <- match(pwt_iso, as.character(sub$isocode))
  if (any(is.na(idx))) next
  r <- rank(sub$rgdpo[idx])
  pwt_ranks[, yi] <- r
}

# Scale our ranks to match PWT scale (our: 1-93, PWT: 1-N)
# Normalize to 0-1
our_ranks_n <- (our_ranks - 1) / (length(ids) - 1)
pwt_ranks_n <- (pwt_ranks - 1) / (n_pwt - 1)

# Match each our country to best PWT country
matches <- data.frame(id=ids, isocode=NA_character_, dist=Inf, stringsAsFactors=FALSE)

for (ii in seq_along(ids)) {
  our_vec <- our_ranks_n[ii, ]
  if (any(is.na(our_vec))) next

  dists <- apply(pwt_ranks_n, 1, function(x) {
    if (any(is.na(x))) return(Inf)
    sqrt(sum((x - our_vec)^2))
  })
  best <- which.min(dists)
  matches$isocode[ii] <- pwt_iso[best]
  matches$dist[ii] <- dists[best]
}

# Add country names
pwt_names <- unique(data.frame(
  isocode = as.character(pwt$isocode),
  country = as.character(pwt$country),
  stringsAsFactors = FALSE
))
matches <- merge(matches, pwt_names, by="isocode", all.x=TRUE)
matches <- matches[order(matches$id),]

cat("\nFull mapping:\n")
print(matches[, c("id","isocode","country","dist")])
cat("\nMax dist:", max(matches$dist, na.rm=TRUE), "\n")

dups <- matches$isocode[duplicated(matches$isocode)]
cat("Duplicate isocodes:", length(dups), "\n")
if (length(dups) > 0) {
  cat(paste(unique(dups), collapse=", "), "\n")
}

write.csv(matches[order(matches$id), c("id","isocode","country")],
          "validation/id_country_map.csv", row.names=FALSE)
cat("Saved to validation/id_country_map.csv\n")
