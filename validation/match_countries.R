library(haven)
library(pwt10)

d <- read_dta("examples/penn_sample.dta")
data("pwt10.0")
pwt <- pwt10.0

# PWT 7 and PWT 10 differ in base year (2005 vs 2017 USD)
# But the ratio between countries in the same year should be ~same
# Use multiple years for robust matching

# Get our data rgdpo for years with good coverage
years_match <- c(1970, 1980, 1990, 2000)
our_wide <- data.frame(id=sort(unique(d$id)))
for (yr in years_match) {
  sub <- d[d$year == yr, c("id","log_rgdpo")]
  sub$rgdpo_n <- exp(sub$log_rgdpo)
  # Normalize within year
  sub$rgdpo_n <- sub$rgdpo_n / sum(sub$rgdpo_n, na.rm=TRUE)
  names(sub)[3] <- paste0("r", yr)
  our_wide <- merge(our_wide, sub[,c("id",paste0("r",yr))], by="id", all.x=TRUE)
}

# PWT countries with data in all match years
pwt_wide <- data.frame(isocode=unique(pwt$isocode[pwt$year %in% years_match]))
for (yr in years_match) {
  sub <- pwt[pwt$year == yr & !is.na(pwt$rgdpo), c("isocode","rgdpo")]
  sub$rgdpo_n <- sub$rgdpo / sum(sub$rgdpo, na.rm=TRUE)
  names(sub)[2] <- paste0("r", yr)
  pwt_wide <- merge(pwt_wide, sub, by="isocode", all.x=TRUE)
}
pwt_wide <- pwt_wide[complete.cases(pwt_wide),]
cat("PWT countries with full data:", nrow(pwt_wide), "\n")

# For each of our 93 countries, find best match in PWT
our_ids <- sort(unique(d$id))
matches <- data.frame(id=our_ids, isocode=NA_character_, dist=Inf)

for (i in seq_along(our_ids)) {
  oid <- our_ids[i]
  our_row <- our_wide[our_wide$id == oid, paste0("r", years_match)]
  if (any(is.na(our_row))) next

  # Euclidean distance to each PWT country
  dists <- apply(pwt_wide[, paste0("r", years_match)], 1, function(x) {
    sqrt(sum((x - as.numeric(our_row))^2))
  })
  best_idx <- which.min(dists)
  matches$isocode[i] <- pwt_wide$isocode[best_idx]
  matches$dist[i] <- dists[best_idx]
}

# Get country names
pwt_names <- unique(pwt[, c("isocode","country")])
matches <- merge(matches, pwt_names, by="isocode", all.x=TRUE)
matches <- matches[order(matches$id),]

cat("\nTop of mapping (id -> country):\n")
print(head(matches[, c("id","isocode","country","dist")], 20))
cat("\nMax distance:", max(matches$dist, na.rm=TRUE), "\n")
cat("Mean distance:", mean(matches$dist, na.rm=TRUE), "\n")

# Check for duplicate matches
dup_iso <- matches$isocode[duplicated(matches$isocode)]
if (length(dup_iso) > 0) {
  cat("\nDuplicate matches:", paste(dup_iso, collapse=", "), "\n")
} else {
  cat("\nAll matches unique!\n")
}

# Save the mapping
write.csv(matches[, c("id","isocode","country")],
          "validation/id_country_map.csv", row.names=FALSE)
cat("\nSaved to validation/id_country_map.csv\n")
