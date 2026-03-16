library(haven)
library(pwt10)

# Load our sample data
d <- read_dta("examples/penn_sample.dta")

# Load PWT 10.x data which has country names and ISO codes
data("pwt10.0")
pwt <- pwt10.0

cat("PWT cols:", paste(names(pwt)[1:10], collapse=", "), "\n")

# Our data: log_rgdpo is log of rgdpo; reverse to get rgdpo
# Try to match by year range and rgdpo values
# Our data years: 1960-2007, so subset PWT to same range
pwt_sub <- pwt[pwt$year >= 1960 & pwt$year <= 2007, ]

# Get unique countries in PWT
pwt_countries <- unique(pwt_sub[, c("country", "isocode")])
cat("PWT countries in 1960-2007:", nrow(pwt_countries), "\n")

# Our data has exp(log_rgdpo) = rgdpo (in millions 2017 USD)
# Match by finding which PWT countries have data that aligns with our sample
# For each numeric id in our data, try to match rgdpo values to PWT

# Get one obs per group in our data
d_ids <- aggregate(list(log_rgdpo=d$log_rgdpo, year=d$year),
                   list(id=d$id),
                   function(x) x[1])
d_ids <- d_ids[order(d_ids$id), ]
cat("First few id/year/log_rgdpo:", "\n")
print(head(d_ids, 5))

# For matching: take 1960 values (or first available year)
d_1965 <- d[d$year == 1965, c("id", "log_rgdpo")]
d_1965$rgdpo <- exp(d_1965$log_rgdpo)
cat("\nOur data 1965 sample rgdpo (id 1-5):\n")
print(head(d_1965[order(d_1965$id),], 5))

# PWT 1965 rgdpo
pwt_1965 <- pwt_sub[pwt_sub$year == 1965, c("country", "isocode", "rgdpo")]
pwt_1965 <- pwt_1965[!is.na(pwt_1965$rgdpo), ]
cat("\nPWT 1965 sample rgdpo:\n")
print(head(pwt_1965[order(pwt_1965$rgdpo),], 5))
