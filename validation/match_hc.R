library(haven)
library(pwt10)

d <- read_dta("examples/penn_sample.dta")
data("pwt10.0")
pwt <- pwt10.0

# The PWT Human Capital Index (hc) uses the same methodology across versions
# log_hc in our data should closely match log(hc) in PWT10
# hc is bounded 1-4+ range and country-specific pattern is version-independent

our_years <- sort(unique(d$year))
ids <- sort(unique(d$id))

# For PWT, find countries with hc data in all our years
pwt_sub <- pwt[pwt$year %in% our_years & !is.na(pwt$hc), ]
pwt_always <- names(which(table(pwt_sub$isocode) == length(our_years)))
cat("PWT countries with hc in all years:", length(pwt_always), "\n")

pwt_hc_sub <- pwt_sub[as.character(pwt_sub$isocode) %in% pwt_always, ]
pwt_iso <- sort(unique(as.character(pwt_hc_sub$isocode)))
n_pwt <- length(pwt_iso)

# Build matrices: rows=country, cols=year
build_mat <- function(df, iso_list, yr_list, val_col) {
  mat <- matrix(NA_real_, length(iso_list), length(yr_list))
  for (yi in seq_along(yr_list)) {
    sub <- df[df$year == yr_list[yi], c("isocode", val_col)]
    for (ci in seq_along(iso_list)) {
      v <- sub[[val_col]][as.character(sub$isocode) == iso_list[ci]]
      if (length(v) == 1 && !is.na(v)) mat[ci, yi] <- v
    }
  }
  mat
}

cat("Building PWT hc matrix...\n")
pwt_mat <- build_mat(pwt_hc_sub, pwt_iso, our_years, "hc")
# log transform
pwt_log_hc <- log(pwt_mat)

# Our data log_hc matrix
our_mat <- matrix(NA_real_, length(ids), length(our_years))
for (yi in seq_along(our_years)) {
  yr <- our_years[yi]
  sub <- d[d$year == yr, c("id","log_hc")]
  for (ii in seq_along(ids)) {
    v <- sub$log_hc[sub$id == ids[ii]]
    if (length(v) == 1) our_mat[ii, yi] <- v
  }
}

cat("Our data missing hc:", sum(is.na(our_mat)), "\n")
cat("PWT missing hc:", sum(is.na(pwt_log_hc)), "\n")

# Match each our id to PWT country via minimum L2 distance on log_hc
matches <- data.frame(id=ids, isocode=NA_character_, dist=Inf, stringsAsFactors=FALSE)

for (ii in seq_along(ids)) {
  our_vec <- our_mat[ii, ]
  if (sum(!is.na(our_vec)) < 10) next

  valid_cols <- !is.na(our_vec)
  dists <- apply(pwt_log_hc, 1, function(x) {
    v <- valid_cols & !is.na(x)
    if (sum(v) < 10) return(Inf)
    sqrt(mean((x[v] - our_vec[v])^2))
  })
  best <- which.min(dists)
  matches$isocode[ii] <- pwt_iso[best]
  matches$dist[ii] <- dists[best]
}

# Country names
pwt_names <- unique(data.frame(
  isocode = as.character(pwt$isocode),
  country = as.character(pwt$country),
  stringsAsFactors = FALSE
))
matches <- merge(matches, pwt_names, by="isocode", all.x=TRUE)
matches <- matches[order(matches$id),]

cat("\nFull mapping:\n")
print(matches[, c("id","isocode","country","dist")])

dups <- matches$isocode[duplicated(matches$isocode)]
cat("\nDuplicate isocodes:", length(dups), "\n")
if (length(dups) > 0) cat(paste(unique(dups), collapse=", "), "\n")

cat("\nMax dist:", max(matches$dist, na.rm=TRUE), "\n")
cat("Mean dist:", mean(matches$dist, na.rm=TRUE), "\n")
cat("Cases with dist > 0.05:", sum(matches$dist > 0.05, na.rm=TRUE), "\n")

write.csv(matches[order(matches$id), c("id","isocode","country")],
          "validation/id_country_map.csv", row.names=FALSE)
cat("Saved to validation/id_country_map.csv\n")
