library(haven)
library(pwt10)

d <- read_dta("examples/penn_sample.dta")
data("pwt10.0")
pwt <- pwt10.0

our_years <- sort(unique(d$year))
ids <- sort(unique(d$id))

# Compute PWT10 variables matching our 4 series:
# 1. log_rgdpo_pc = log(rgdpo/pop)
# 2. log_hc = log(hc)
# 3. log_ck: capital per worker -- PWT has rnna (capital stock) and emp (employed)
#    or cn / pop for capital per person
#    Our log_ck might be log(capital/labor) or log(capital share)...
#    In the MRW/Solow framework ck = capital per efficiency unit of labor
#    but log_ck naming suggests capital per worker or capital stock in logs
# 4. log_ngd = log(n + g + d) where n = pop growth, g=0.02, d=0.03
#    Computed from PWT pop series: n_t = log(pop_t) - log(pop_{t-1})

# Compute log_ngd from PWT pop
pwt <- pwt[order(as.character(pwt$isocode), pwt$year), ]
pwt$n_pop <- c(NA, diff(log(pwt$pop)))
pwt$n_pop[!duplicated(as.character(pwt$isocode))] <- NA  # first obs per country
pwt$ngd <- pwt$n_pop + 0.05  # n + g + delta
pwt$log_ngd <- log(ifelse(pwt$ngd > 0, pwt$ngd, NA))

# Capital per worker (log)
pwt$log_ck <- log(pwt$rnna / pwt$emp)  # rnna: capital, emp: employment

pwt$log_rgdpo_pc <- log(pwt$rgdpo / pwt$pop)
pwt$log_hc <- log(pwt$hc)

# Find balanced countries (all 4 vars complete for 1960-2007)
# Note: log_ngd needs lag so starts from 1961
eff_years <- our_years[our_years >= 1961]
pwt_sub <- pwt[pwt$year %in% eff_years, ]

balanced_iso <- NULL
for (iso in unique(as.character(pwt_sub$isocode))) {
  sub <- pwt_sub[as.character(pwt_sub$isocode)==iso,
                 c("log_rgdpo_pc","log_hc","log_ck","log_ngd")]
  if (nrow(sub) == length(eff_years) && all(!is.na(sub))) {
    balanced_iso <- c(balanced_iso, iso)
  }
}
cat("PWT10 balanced (all 4 vars, 1961-2007):", length(balanced_iso), "\n")

# Fallback: use just rgdpo_pc + hc for full year range
pwt_sub2 <- pwt[pwt$year %in% our_years, ]
balanced_iso2 <- NULL
for (iso in unique(as.character(pwt_sub2$isocode))) {
  sub <- pwt_sub2[as.character(pwt_sub2$isocode)==iso,
                  c("log_rgdpo_pc","log_hc")]
  if (nrow(sub) == length(our_years) && all(!is.na(sub))) {
    balanced_iso2 <- c(balanced_iso2, iso)
  }
}
cat("PWT10 balanced (rgdpo+hc, 1960-2007):", length(balanced_iso2), "\n")

# Use the eff_years (1961-2007) with all 4 vars
pwt_full <- pwt[pwt$year %in% eff_years & as.character(pwt$isocode) %in% balanced_iso, ]
pwt_iso <- sort(balanced_iso)
n_pwt <- length(pwt_iso)
d_eff  <- d[d$year %in% eff_years, ]
eff_n  <- length(eff_years)

# Build matrices for our data (using 1961-2007)
build_our <- function(varname) {
  mat <- matrix(NA_real_, 93, eff_n)
  for (yi in seq_along(eff_years)) {
    sub <- d_eff[d_eff$year == eff_years[yi], c("id", varname)]
    for (ii in seq_along(ids)) {
      v <- sub[[varname]][sub$id == ids[ii]]
      if (length(v)==1) mat[ii, yi] <- v
    }
  }
  mat
}
build_pwt <- function(varname) {
  mat <- matrix(NA_real_, n_pwt, eff_n)
  for (yi in seq_along(eff_years)) {
    sub <- pwt_full[pwt_full$year == eff_years[yi], c("isocode", varname)]
    for (ci in seq_along(pwt_iso)) {
      v <- sub[[varname]][as.character(sub$isocode) == pwt_iso[ci]]
      if (length(v)==1) mat[ci, yi] <- v
    }
  }
  mat
}

cat("Building matrices...\n")
our_gdp <- build_our("log_rgdpo")
our_hc  <- build_our("log_hc")
our_ck  <- build_our("log_ck")
our_ngd <- build_our("log_ngd")
pwt_gdp <- build_pwt("log_rgdpo_pc")
pwt_hc  <- build_pwt("log_hc")
pwt_ck  <- build_pwt("log_ck")
pwt_ngd <- build_pwt("log_ngd")

# Combined score with weights
# log_ngd and log_ck are in same units as our data -> direct comparison (no var method)
# log_rgdpo has base-year offset -> use var(diff)
# log_hc same across versions (methodology consistent)
cat("Computing scores...\n")
all_scores <- matrix(Inf, 93, n_pwt)
for (ii in seq_along(ids)) {
  g <- our_gdp[ii,]; h <- our_hc[ii,]
  k <- our_ck[ii,];  n <- our_ngd[ii,]
  if (any(is.na(g)) || any(is.na(h)) || any(is.na(k)) || any(is.na(n))) next

  for (ci in seq_along(pwt_iso)) {
    gp <- pwt_gdp[ci,]; hp <- pwt_hc[ci,]
    kp <- pwt_ck[ci,];  np <- pwt_ngd[ci,]
    if (any(is.na(gp))||any(is.na(hp))||any(is.na(kp))||any(is.na(np))) next

    s_gdp <- var(g - gp)
    s_hc  <- mean((h - hp)^2)
    s_ck  <- var(k - kp)   # capital also has base-year offset
    s_ngd <- mean((n - np)^2)
    all_scores[ii, ci] <- s_gdp + 3*s_hc + s_ck + 5*s_ngd
  }
}

# Greedy 1-to-1 assignment
flat <- data.frame(
  ii=rep(seq_along(ids), n_pwt),
  ci=rep(seq_len(n_pwt), each=length(ids)),
  score=as.vector(all_scores)
)
flat <- flat[order(flat$score),]
flat <- flat[is.finite(flat$score),]

assigned_id  <- rep(FALSE, length(ids))
assigned_pwt <- rep(FALSE, n_pwt)
matches <- data.frame(id=ids, isocode=NA_character_, dist=Inf, stringsAsFactors=FALSE)

for (row in seq_len(nrow(flat))) {
  ii <- flat$ii[row]; ci <- flat$ci[row]
  if (!assigned_id[ii] && !assigned_pwt[ci]) {
    matches$isocode[ii] <- pwt_iso[ci]
    matches$dist[ii]    <- flat$score[row]
    assigned_id[ii]     <- TRUE
    assigned_pwt[ci]    <- TRUE
    if (all(assigned_id)) break
  }
}

pwt_names <- unique(data.frame(
  isocode=as.character(pwt$isocode),
  country=as.character(pwt$country),
  stringsAsFactors=FALSE
))
matches <- merge(matches, pwt_names, by="isocode", all.x=TRUE)
matches <- matches[order(matches$id),]

cat("\nFull mapping (4-variable greedy):\n")
print(matches[, c("id","isocode","country","dist")])
cat("\nMax:", max(matches$dist, na.rm=TRUE))
cat("  Mean:", mean(matches$dist, na.rm=TRUE), "\n")
cat("< 0.01:", sum(matches$dist < 0.01), "\n")
cat("< 0.05:", sum(matches$dist < 0.05), "\n")
cat("> 0.1:", sum(matches$dist > 0.1), "\n")

write.csv(matches[order(matches$id), c("id","isocode","country")],
          "validation/id_country_map_4var.csv", row.names=FALSE)
cat("Saved to id_country_map_4var.csv\n")
