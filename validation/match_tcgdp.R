library(haven)
library(pwt)

d <- read_dta("examples/penn_sample.dta")
data("pwt7.1")
p7 <- pwt7.1

# tcgdp = total current GDP (millions)
# rgdpl = real GDP per capita (Laspeyres)
# rgdpch = real GDP per capita (chain)
# To get total: rgdpch * pop = total real GDP

# Our log_rgdpo is likely log of total output-side GDP
# Try: log(rgdpch * pop) vs log_rgdpo

p7$total_gdp <- p7$rgdpch * p7$pop  # Total real GDP (chain)

our_years <- sort(unique(d$year))
ids <- sort(unique(d$id))

p7_sub <- p7[p7$year %in% our_years & !is.na(p7$total_gdp) & p7$total_gdp > 0, ]
pwt7_iso <- sort(unique(as.character(p7_sub$isocode)))
n_pwt <- length(pwt7_iso)

# Demean by year
our_mat <- matrix(NA_real_, 93, length(our_years))
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
  sub <- p7_sub[p7_sub$year == yr, c("isocode","total_gdp")]
  sub$log_v <- log(sub$total_gdp)
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

cat("Full mapping (PWT7 tcgdp):\n")
print(matches[, c("id","isocode","country","dist")])

dups <- matches$isocode[duplicated(matches$isocode)]
cat("\nDuplicates:", length(dups), "\n")
cat("Max dist:", max(matches$dist, na.rm=TRUE), "\n")
cat("Mean dist:", mean(matches$dist, na.rm=TRUE), "\n")
cat("N with dist < 0.005:", sum(matches$dist < 0.005, na.rm=TRUE), "\n")
