library(haven)
library(pwt10)

d <- read_dta("examples/penn_sample.dta")
data("pwt10.0")
pwt <- pwt10.0

our_years <- sort(unique(d$year))  # 1960-2007
ids <- sort(unique(d$id))

# PWT10 rgdpo is output-side real GDP (millions 2017 USD)
# Our log_rgdpo should be log of this (possibly different base year -> uniform offset)

# Step 1: Find PWT10 countries with rgdpo in all 48 years (1960-2007)
pwt_sub <- pwt[pwt$year %in% our_years & !is.na(pwt$rgdpo) & pwt$rgdpo > 0, ]
country_counts <- table(as.character(pwt_sub$isocode))
balanced_iso <- names(country_counts[country_counts == length(our_years)])
cat("PWT10 countries with complete rgdpo 1960-2007:", length(balanced_iso), "\n")

# If we find exactly 93, that's our set!
# If not, also check with hc and ck conditions
pwt_full <- pwt_sub[as.character(pwt_sub$isocode) %in% balanced_iso, ]

# Also check which have complete hc (needed for x_csa)
hc_sub <- pwt[pwt$year %in% our_years & !is.na(pwt$hc), ]
hc_counts <- table(as.character(hc_sub$isocode))
hc_balanced <- names(hc_counts[hc_counts == length(our_years)])
cat("PWT10 countries with complete hc 1960-2007:", length(hc_balanced), "\n")

# Intersection: complete rgdpo AND hc
both_balanced <- intersect(balanced_iso, hc_balanced)
cat("Both rgdpo and hc complete:", length(both_balanced), "\n")

# If we have more than 93, need additional filter
# Also check ck (capital stock) and pop
ck_sub <- pwt[pwt$year %in% our_years & !is.na(pwt$rnna), ]
ck_counts <- table(as.character(ck_sub$isocode))
ck_balanced <- names(ck_counts[ck_counts == length(our_years)])
cat("Complete rnna (capital):", length(ck_balanced), "\n")

# Try cn (capital stock nominal) or ck
if ("ck" %in% names(pwt)) {
  ck_col <- "ck"
} else if ("rnna" %in% names(pwt)) {
  ck_col <- "rnna"
} else {
  cat("Capital cols:", paste(grep("cn|ck|rnna|cap", names(pwt), value=TRUE), collapse=", "), "\n")
  ck_col <- NULL
}

# Now: match our 93 countries to PWT10 using direct comparison
# Log-difference approach: log_rgdpo (our) vs log(rgdpo) (PWT10)
# Should differ only by a constant (base year difference) WITHIN each country

# For each of our 93 IDs, match to PWT10 by minimizing variance of the difference
# Var(log_rgdpo_ours - log_rgdpo_pwt10) should be near 0 for the correct match

pwt_iso <- sort(balanced_iso)
n_pwt <- length(pwt_iso)

# Build our time series matrix
our_mat <- matrix(NA_real_, 93, length(our_years))
for (yi in seq_along(our_years)) {
  yr <- our_years[yi]
  sub <- d[d$year == yr, c("id","log_rgdpo")]
  for (ii in seq_along(ids)) {
    v <- sub$log_rgdpo[sub$id == ids[ii]]
    if (length(v) == 1) our_mat[ii, yi] <- v
  }
}

# Build PWT10 time series matrix
pwt_mat <- matrix(NA_real_, n_pwt, length(our_years))
for (yi in seq_along(our_years)) {
  yr <- our_years[yi]
  sub <- pwt_full[pwt_full$year == yr, c("isocode","rgdpo")]
  for (ci in seq_along(pwt_iso)) {
    v <- sub$rgdpo[as.character(sub$isocode) == pwt_iso[ci]]
    if (length(v) == 1) pwt_mat[ci, yi] <- log(v)
  }
}

# For each our ID, find PWT country that minimizes variance of (our - pwt)
# This handles the base-year offset
matches <- data.frame(id=ids, isocode=NA_character_, dist=Inf, stringsAsFactors=FALSE)

for (ii in seq_along(ids)) {
  our_vec <- our_mat[ii, ]
  if (any(is.na(our_vec))) next

  scores <- apply(pwt_mat, 1, function(x) {
    if (any(is.na(x))) return(Inf)
    # Remove mean difference (handles base year scaling)
    diff <- our_vec - x
    var(diff)  # Should be ~0 if same country
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

cat("\nFull mapping (variance method):\n")
print(matches[, c("id","isocode","country","dist")])

dups <- matches$isocode[duplicated(matches$isocode)]
cat("\nDuplicates:", length(dups), "\n")
cat("Max variance:", max(matches$dist, na.rm=TRUE), "\n")
cat("Mean variance:", mean(matches$dist, na.rm=TRUE), "\n")
cat("Near-zero matches (var < 0.0001):", sum(matches$dist < 0.0001, na.rm=TRUE), "\n")
cat("Good matches (var < 0.001):", sum(matches$dist < 0.001, na.rm=TRUE), "\n")

if (length(dups) == 0) {
  write.csv(matches[order(matches$id), c("id","isocode","country")],
            "validation/id_country_map.csv", row.names=FALSE)
  cat("Perfect mapping saved!\n")
}
