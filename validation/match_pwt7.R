library(haven)
library(pwt)

d <- read_dta("examples/penn_sample.dta")

# Check what datasets pwt package provides
cat("pwt datasets:", paste(data(package="pwt")$results[,"Item"], collapse=", "), "\n")
data("pwt7.1")
p7 <- pwt7.1
cat("PWT 7.1 cols:", paste(names(p7)[1:15], collapse=", "), "\n")
cat("PWT 7.1 dims:", nrow(p7), "x", ncol(p7), "\n")
cat("PWT 7.1 year range:", min(p7$year), "-", max(p7$year), "\n")

# Check country identifier
cat("isocode sample:", paste(head(unique(as.character(p7$isocode)), 5), collapse=", "), "\n")
cat("country sample:", paste(head(unique(as.character(p7$country)), 5), collapse=", "), "\n")

# Check rgdpch or rgdpl column (PWT7 uses rgdpl = real GDP per capita, laspeyres)
# Our log_rgdpo - is it rgdpch or rgdpl?
rgdp_cols <- grep("rgdp", names(p7), value=TRUE)
cat("rgdp columns in pwt7:", paste(rgdp_cols, collapse=", "), "\n")

# PWT 7.1 has "rgdpch" = chain-weighted real GDP per capita
# "rgdpl" = Laspeyres index real GDP per capita
# Our data: "rgdpo" might refer to output-based measure
# In PWT 7.1, "rgdpna" doesn't exist; try "rgdpl" or "rgdpch"

# Check years
cat("PWT 7.1 years:", min(p7$year), "to", max(p7$year), "\n")
our_years <- sort(unique(d$year))
cat("Our years:", min(our_years), "to", max(our_years), "\n")

# PWT 7 countries in 1960-2007
p7_sub <- p7[p7$year %in% our_years, ]
n_by_country <- table(as.character(p7_sub$isocode))
pwt7_all <- names(n_by_country[n_by_country == length(our_years)])
cat("PWT 7.1 countries with data in all years:", length(pwt7_all), "\n")

# Pick best rgdp column
# In PWT 7.1, "rgdpch" is chain-weighted GDP per capita
# Our "log_rgdpo" is likely log of this (or closely related)
# Let's check "rgdpch" first
if ("rgdpch" %in% names(p7)) {
  cat("Using rgdpch\n")
  val_col <- "rgdpch"
} else if ("rgdpl" %in% names(p7)) {
  cat("Using rgdpl\n")
  val_col <- "rgdpl"
} else {
  cat("Available rgdp cols:", paste(rgdp_cols, collapse=", "), "\n")
  val_col <- rgdp_cols[1]
}

# Match using log of rgdp column, demeaned by year
p7_sub2 <- p7_sub[as.character(p7_sub$isocode) %in% pwt7_all & !is.na(p7_sub[[val_col]]), ]
pwt7_iso <- sort(unique(as.character(p7_sub2$isocode)))
n_pwt <- length(pwt7_iso)

our_mat <- matrix(NA_real_, 93, length(our_years))
ids <- sort(unique(d$id))

for (yi in seq_along(our_years)) {
  yr <- our_years[yi]
  sub <- d[d$year == yr, c("id","log_rgdpo")]
  ymean <- mean(sub$log_rgdpo, na.rm=TRUE)
  for (ii in seq_along(ids)) {
    v <- sub$log_rgdpo[sub$id == ids[ii]]
    if (length(v) == 1) our_mat[ii, yi] <- v - ymean
  }
}

pwt_mat <- matrix(NA_real_, n_pwt, length(our_years))
for (yi in seq_along(our_years)) {
  yr <- our_years[yi]
  sub <- p7_sub2[p7_sub2$year == yr, c("isocode", val_col)]
  sub$log_v <- log(sub[[val_col]])
  ymean <- mean(sub$log_v, na.rm=TRUE)
  for (ci in seq_along(pwt7_iso)) {
    v <- sub$log_v[as.character(sub$isocode) == pwt7_iso[ci]]
    if (length(v) == 1) pwt_mat[ci, yi] <- v - ymean
  }
}

# Match
matches <- data.frame(id=ids, isocode=NA_character_, dist=Inf, stringsAsFactors=FALSE)
for (ii in seq_along(ids)) {
  our_vec <- our_mat[ii, ]
  valid <- !is.na(our_vec)
  if (sum(valid) < 20) next
  dists <- apply(pwt_mat, 1, function(x) {
    v <- valid & !is.na(x)
    if (sum(v) < 20) return(Inf)
    sqrt(mean((x[v] - our_vec[v])^2))
  })
  best <- which.min(dists)
  matches$isocode[ii] <- pwt7_iso[best]
  matches$dist[ii] <- dists[best]
}

pwt7_names <- unique(data.frame(
  isocode=as.character(p7$isocode),
  country=as.character(p7$country),
  stringsAsFactors=FALSE
))
matches <- merge(matches, pwt7_names, by="isocode", all.x=TRUE)
matches <- matches[order(matches$id),]

cat("\nFull mapping (PWT7):\n")
print(matches[, c("id","isocode","country","dist")])

dups <- matches$isocode[duplicated(matches$isocode)]
cat("\nDuplicates:", length(dups), "\n")
cat("Max dist:", max(matches$dist, na.rm=TRUE), "\n")
cat("Mean dist:", mean(matches$dist, na.rm=TRUE), "\n")
cat("N with dist < 0.01:", sum(matches$dist < 0.01, na.rm=TRUE), "\n")
cat("N with dist < 0.05:", sum(matches$dist < 0.05, na.rm=TRUE), "\n")
