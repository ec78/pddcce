library(haven)
library(pwt10)

d <- read_dta("examples/penn_sample.dta")
data("pwt10.0")
pwt <- pwt10.0

# Check the scale of log_rgdpo in our data
# Look at 2000 values - sort by size
d2000 <- d[d$year == 2000, c("id","log_rgdpo")]
d2000 <- d2000[order(d2000$log_rgdpo, decreasing=TRUE),]
cat("Top 10 by log_rgdpo in 2000:\n")
print(head(d2000, 10))
cat("exp of these:\n")
print(head(exp(d2000$log_rgdpo), 10))

# Check PWT10 for same year
pwt2000 <- pwt[pwt$year == 2000 & !is.na(pwt$rgdpo), c("country","isocode","rgdpo","pop")]
pwt2000$log_rgdpo_bn <- log(pwt2000$rgdpo/1000)  # if rgdpo in millions -> /1000 = billions
pwt2000$log_rgdpo_tr <- log(pwt2000$rgdpo/1e6)    # trillions
pwt2000$log_rgdpo_raw <- log(pwt2000$rgdpo)        # raw millions
cat("\nPWT10 2000: US, Japan, Germany (raw log):\n")
for (cc in c("USA","JPN","DEU","GBR","FRA","CHN","IND")) {
  row <- pwt2000[as.character(pwt2000$isocode)==cc,]
  if (nrow(row) > 0) {
    cat(cc, "raw_log=", row$log_rgdpo_raw, "bn_log=", row$log_rgdpo_bn,
        "pop=", row$pop, "\n")
  }
}

# Try: maybe log_rgdpo is log GDP per capita in thousands
# For US 2000: GDP ~$36000 per capita. log(36000) ≈ 10.49. Hmm close to some top values.
# Or maybe it's in 2005 international USD (common for older PWT)?

# Let's check PWT9
if (require(pwt9, quietly=TRUE)) {
  data("pwt9.1")
  p9 <- pwt9.1
  p9_2000 <- p9[p9$year == 2000 & !is.na(p9$rgdpo), c("country","isocode","rgdpo")]
  cat("\nPWT9.1 2000 US:\n")
  print(p9_2000[as.character(p9_2000$isocode)=="USA",])
} else {
  cat("pwt9 not available\n")
}

# Bottom 10 our data
d2000_sorted <- d2000[order(d2000$log_rgdpo),]
cat("\nBottom 5 by log_rgdpo in 2000:\n")
print(head(d2000_sorted, 5))

# Try matching: our top GDP country in 2000 should be USA
# id with max log_rgdpo in 2000:
top_id <- d2000[1, "id"]
cat("\nOur top GDP country id:", top_id, "\n")
cat("log_rgdpo:", d2000[1, "log_rgdpo"], "-> exp:", exp(d2000[1, "log_rgdpo"]), "\n")

# Check if this is US in PWT10 (US 2000 rgdpo = ?)
us_2000 <- pwt2000[as.character(pwt2000$isocode)=="USA",]
cat("US 2000 rgdpo (millions 2017):", us_2000$rgdpo, "\n")
cat("log(US rgdpo):", log(us_2000$rgdpo), "\n")

# Check ratio
ratio <- exp(d2000[1, "log_rgdpo"]) / us_2000$rgdpo
cat("Ratio (our/PWT10):", ratio, "\n")

# Maybe our data is from PWT 9.1 (base year 2011)
# Install pwt9 if needed
