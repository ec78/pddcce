library(haven)
library(pwt10)

d <- read_dta("examples/penn_sample.dta")
data("pwt10.01")
pwt <- pwt10.01

our_years <- sort(unique(d$year))
ids <- sort(unique(d$id))

# Try both rgdpo/pop and cgdpo/pop to see which fits better
pwt$log_rgdpo_pc <- log(pwt$rgdpo / pwt$pop)
pwt$log_cgdpo_pc <- log(pwt$cgdpo / pwt$pop)
pwt$log_hc <- log(pwt$hc)

# Check USA id=89 with both variants
us <- pwt[as.character(pwt$isocode)=="USA" & pwt$year %in% our_years,
          c("year","log_rgdpo_pc","log_cgdpo_pc")]
our_us <- d[d$id==89, c("year","log_rgdpo")]

df <- merge(us, our_us, by="year")
cat("USA comparison (first 5 rows):\n")
print(head(df))
cat("var(our - rgdpo):", var(df$log_rgdpo - df$log_rgdpo_pc), "\n")
cat("var(our - cgdpo):", var(df$log_rgdpo - df$log_cgdpo_pc), "\n")

# Check Germany (should be around id=20 or id=30)
# Try all Germany-size IDs
d2000 <- d[d$year==2000, c("id","log_rgdpo")]
deu_2000_rgdpo_pc <- log(pwt$rgdpo[as.character(pwt$isocode)=="DEU" & pwt$year==2000] /
                         pwt$pop[as.character(pwt$isocode)=="DEU" & pwt$year==2000])
deu_2000_cgdpo_pc <- log(pwt$cgdpo[as.character(pwt$isocode)=="DEU" & pwt$year==2000] /
                         pwt$pop[as.character(pwt$isocode)=="DEU" & pwt$year==2000])
cat("\nDEU 2000 log(rgdpo/pop):", deu_2000_rgdpo_pc, "\n")
cat("DEU 2000 log(cgdpo/pop):", deu_2000_cgdpo_pc, "\n")

# Which ids are close to DEU level?
d2000_near_deu <- d2000[abs(d2000$log_rgdpo - deu_2000_rgdpo_pc) < 0.5, ]
cat("IDs with log_rgdpo close to DEU 2000:", paste(d2000_near_deu$id, collapse=", "), "\n")

# Try cgdpo matching for all 93 countries
# Find balanced
pwt_sub <- pwt[pwt$year %in% our_years, ]
balanced_iso <- NULL
for (iso in unique(as.character(pwt_sub$isocode))) {
  sub <- pwt_sub[as.character(pwt_sub$isocode)==iso, c("log_cgdpo_pc","log_hc")]
  if (nrow(sub) == length(our_years) && all(!is.na(sub))) {
    balanced_iso <- c(balanced_iso, iso)
  }
}
cat("\nPWT10.01 balanced (cgdpo_pc + hc) for 1960-2007:", length(balanced_iso), "\n")

pwt_full <- pwt_sub[as.character(pwt_sub$isocode) %in% balanced_iso, ]
pwt_iso <- sort(balanced_iso)
n_pwt <- length(pwt_iso)

build_pwt <- function(varname) {
  mat <- matrix(NA_real_, n_pwt, length(our_years))
  for (yi in seq_along(our_years)) {
    sub <- pwt_full[pwt_full$year == our_years[yi], c("isocode", varname)]
    for (ci in seq_along(pwt_iso)) {
      v <- sub[[varname]][as.character(sub$isocode) == pwt_iso[ci]]
      if (length(v)==1) mat[ci, yi] <- v
    }
  }
  mat
}

build_our <- function(varname) {
  mat <- matrix(NA_real_, 93, length(our_years))
  for (yi in seq_along(our_years)) {
    sub <- d[d$year == our_years[yi], c("id", varname)]
    for (ii in seq_along(ids)) {
      v <- sub[[varname]][sub$id == ids[ii]]
      if (length(v)==1) mat[ii, yi] <- v
    }
  }
  mat
}

our_gdp <- build_our("log_rgdpo")
our_hc  <- build_our("log_hc")
pwt_cgdp <- build_pwt("log_cgdpo_pc")
pwt_hc   <- build_pwt("log_hc")
pwt_rgdp <- build_pwt("log_rgdpo_pc")

# Check which variant (rgdpo or cgdpo) gives lower variance for US=89
us_ii <- which(ids==89)
us_ci_rgdp <- which(pwt_iso == "USA")

cat("\nUS var(our - rgdpo_pc):", var(our_gdp[us_ii,] - pwt_rgdp[us_ci_rgdp,]), "\n")
cat("US var(our - cgdpo_pc):", var(our_gdp[us_ii,] - pwt_cgdp[us_ci_rgdp,]), "\n")

# Try matching with cgdpo
all_scores_c <- matrix(Inf, 93, n_pwt)
for (ii in seq_along(ids)) {
  g <- our_gdp[ii,]; h <- our_hc[ii,]
  if (any(is.na(g)) || any(is.na(h))) next
  for (ci in seq_along(pwt_iso)) {
    gp <- pwt_cgdp[ci,]; hp <- pwt_hc[ci,]
    if (any(is.na(gp))||any(is.na(hp))) next
    s_gdp <- var(g - gp)
    s_hc  <- mean((h - hp)^2)
    all_scores_c[ii, ci] <- s_gdp + 5 * s_hc
  }
}

# Greedy 1-to-1
flat <- data.frame(ii=rep(seq_along(ids),n_pwt), ci=rep(1:n_pwt,each=length(ids)),
                   score=as.vector(all_scores_c))
flat <- flat[order(flat$score),]; flat <- flat[is.finite(flat$score),]
a_id <- rep(F,length(ids)); a_pwt <- rep(F,n_pwt)
matches_c <- data.frame(id=ids, isocode=NA_character_, dist=Inf, stringsAsFactors=FALSE)
for (row in seq_len(nrow(flat))) {
  ii <- flat$ii[row]; ci <- flat$ci[row]
  if (!a_id[ii] && !a_pwt[ci]) {
    matches_c$isocode[ii] <- pwt_iso[ci]; matches_c$dist[ii] <- flat$score[row]
    a_id[ii] <- TRUE; a_pwt[ci] <- TRUE
    if (all(a_id)) break
  }
}
pwt_names <- unique(data.frame(isocode=as.character(pwt$isocode),country=as.character(pwt$country),stringsAsFactors=FALSE))
matches_c <- merge(matches_c, pwt_names, by="isocode", all.x=TRUE)
matches_c <- matches_c[order(matches_c$id),]
cat("\ncgdpo matching stats:\n")
cat("> 0.1:", sum(matches_c$dist > 0.1), "\n")
cat("> 0.05:", sum(matches_c$dist > 0.05), "\n")
cat("Mean:", mean(matches_c$dist), "\n")
cat("High score cases:\n")
print(matches_c[matches_c$dist > 0.1, c("id","isocode","country","dist")])
