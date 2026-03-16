library(haven)
library(pwt10)

d <- read_dta("examples/penn_sample.dta")
data("pwt10.01")
pwt <- pwt10.01

our_years <- sort(unique(d$year))
ids <- sort(unique(d$id))

# High-score IDs from previous run:
high_ids <- c(14, 19, 20, 23, 34, 38, 40, 42, 45, 51, 62, 66, 77, 90)
# Mapped to: CAN, DZA, DEU, SGP, SLV, TWN, PRT, IRN, JOR, AUS, COL, CMR, ZWE, TZA

pwt$log_rgdpo_pc <- log(pwt$rgdpo / pwt$pop)

# For each high-score ID, find top 5 best PWT10.01 matches from ALL countries
cat("Top matches for high-score IDs (all PWT countries, rgdpo_pc only):\n\n")
for (hid in high_ids) {
  our_vec <- d$log_rgdpo[d$id == hid]
  our_yr  <- d$year[d$id == hid]

  cat(sprintf("=== id=%d ===\n", hid))

  # Our 2000 value
  v2000 <- our_vec[our_yr == 2000]
  cat(sprintf("  log_rgdpo 2000: %.4f (exp=%.0f)\n", v2000, exp(v2000)))

  scores <- tapply(seq_len(nrow(pwt)), as.character(pwt$isocode), function(rows) {
    p <- pwt[rows, c("year","log_rgdpo_pc")]
    p <- p[p$year %in% our_yr, ]
    if (nrow(p) < length(our_yr)) return(Inf)
    p <- p[order(p$year),]
    var(our_vec - p$log_rgdpo_pc)
  })

  scores_sorted <- sort(scores[is.finite(scores)])
  cat("  Top 5 matches:\n")
  top5 <- head(scores_sorted, 5)
  for (nm in names(top5)) {
    cname <- unique(pwt$country[as.character(pwt$isocode)==nm])
    cat(sprintf("    %s (%s): var=%.6f\n", nm, cname, top5[nm]))
  }
  cat("\n")
}
