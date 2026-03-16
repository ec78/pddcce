library(pwt10)

# Check which version of PWT is in the package
cat("Available datasets:\n")
print(data(package="pwt10")$results[,"Item"])

# Check pwt10.0 - look for version info
data("pwt10.0")
cat("\npwt10.0 dims:", nrow(pwt10.0), "x", ncol(pwt10.0), "\n")
cat("Year range:", min(pwt10.0$year), "-", max(pwt10.0$year), "\n")

# The key insight: our data has no NAs. Check which PWT countries have
# zero NAs in rgdpo for 1960-2007
our_years <- 1960:2007
pwt_sub <- pwt10.0[pwt10.0$year %in% our_years, ]

# Count complete obs per country for rgdpo
rgdpo_counts <- tapply(!is.na(pwt_sub$rgdpo), as.character(pwt_sub$isocode), sum)
full_rgdpo <- names(rgdpo_counts[rgdpo_counts == 48])
cat("PWT10 countries with no NA rgdpo 1960-2007:", length(full_rgdpo), "\n")

# Also need hc
hc_counts <- tapply(!is.na(pwt_sub$hc), as.character(pwt_sub$isocode), sum)
full_hc <- names(hc_counts[hc_counts == 48])
cat("PWT10 countries with no NA hc 1960-2007:", length(full_hc), "\n")

both <- intersect(full_rgdpo, full_hc)
cat("Both rgdpo and hc complete:", length(both), "\n")
cat("Countries:", paste(sort(both), collapse=", "), "\n")
