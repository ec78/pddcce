# =============================================================================
# dccelib Benchmark — R (plm) vs GAUSS dccelib
#
# Measures wall-clock time for MG and CCE-MG estimation on:
#   (a) Penn World Tables sample   (N=93, T~50)
#   (b) Simulated medium panel     (N=100, T=50)
#   (c) Simulated large panel      (N=200, T=100)
#
# Output: console table + validation/benchmark_r_results.csv
#
# Run from repo root:
#   "C:/Program Files/R/R-4.5.0/bin/Rscript.exe" validation/benchmark_r.R
# =============================================================================

library(haven)
library(plm)
library(dplyr)

cat("=============================================================\n")
cat("dccelib Benchmark — R (plm) timing and accuracy\n")
cat("=============================================================\n")
cat("R version:", R.version$version.string, "\n")
cat("plm version:", as.character(packageVersion("plm")), "\n\n")

# ---- Utility: repeat timing B times and return median ----
time_median <- function(B, fn) {
  times <- numeric(B)
  for (i in seq_len(B)) {
    t0 <- proc.time()["elapsed"]
    fn()
    times[i] <- proc.time()["elapsed"] - t0
  }
  median(times)
}

# ---- Utility: simulate a balanced panel ----
sim_panel <- function(N, T, k = 2, seed = 42) {
  set.seed(seed)
  f  <- rnorm(T)                       # single common factor
  id <- rep(1:N, each = T)
  yr <- rep(1:T, times = N)
  x1 <- rnorm(N * T) + 0.5 * rep(f, N)
  x2 <- rnorm(N * T) + 0.3 * rep(f, N)
  y  <- 0.5 * x1 + 0.3 * x2 + rep(rnorm(N), each = T) +
        rep(f, N) + rnorm(N * T, sd = 0.5)
  data.frame(id = id, year = yr, y = y, x1 = x1, x2 = x2)
}

# =============================================================
# PART 1: Accuracy check (Penn World Tables)
# =============================================================

cat("-------------------------------------------------------------\n")
cat("PART 1: Accuracy (Penn World Tables, N=93, T~50)\n")
cat("-------------------------------------------------------------\n")

d <- read_dta("examples/penn_sample.dta")
d <- na.omit(d)
pdata <- pdata.frame(d, index = c("id", "year"))

mg_fit  <- pmg(log_rgdpo ~ log_ck + log_ngd, data = pdata, model = "mg")
cce_fit <- pmg(log_rgdpo ~ log_ck + log_ngd, data = pdata, model = "cmg")

cat("\nMG coefficients:\n")
print(round(coef(mg_fit), 6))
cat("\nCCE-MG coefficients:\n")
print(round(coef(cce_fit), 6))

cat("\nExpected (from GAUSS dccelib validation):\n")
cat("  MG:     log_ck=0.305300  log_ngd=0.279783  inpt=5.391778\n")
cat("  CCE-MG: log_ck=0.316743  log_ngd=0.089055  inpt=1.145539\n\n")

# =============================================================
# PART 2: Timing benchmarks
# =============================================================

B_reps <- 10   # repetitions per timing cell

results <- data.frame(
  dataset      = character(),
  N            = integer(),
  T            = integer(),
  estimator    = character(),
  tool         = character(),
  time_sec     = numeric(),
  stringsAsFactors = FALSE
)

# ---- Panel sizes ----
panel_specs <- list(
  list(name = "Penn WTables",  N = 93,  T = 50, data = NULL),
  list(name = "Simulated M",   N = 100, T = 50, data = NULL),
  list(name = "Simulated L",   N = 200, T = 100, data = NULL)
)

panel_specs[[1]]$data <- d
panel_specs[[2]]$data <- sim_panel(100, 50)
panel_specs[[3]]$data <- sim_panel(200, 100)

for (spec in panel_specs) {
  d_i     <- spec$data
  d_i     <- na.omit(d_i)
  pd_i    <- pdata.frame(d_i, index = c("id", "year"))
  yvar    <- names(d_i)[3]
  xvars   <- names(d_i)[4:5]
  fml     <- as.formula(paste(yvar, "~", paste(xvars, collapse = " + ")))

  cat(sprintf("Dataset: %s  (N=%d, T~%d)\n",
              spec$name,
              length(unique(d_i$id)),
              round(nrow(d_i) / length(unique(d_i$id)))))

  # --- MG timing ---
  t_mg <- time_median(B_reps, function() pmg(fml, data = pd_i, model = "mg"))
  cat(sprintf("  MG   (R/plm):     %.4f s\n", t_mg))
  results <- rbind(results, data.frame(
    dataset = spec$name, N = length(unique(d_i$id)),
    T = round(nrow(d_i) / length(unique(d_i$id))),
    estimator = "MG", tool = "R/plm", time_sec = t_mg
  ))

  # --- CCE-MG timing ---
  t_cce <- time_median(B_reps, function() pmg(fml, data = pd_i, model = "cmg"))
  cat(sprintf("  CCE-MG (R/plm):   %.4f s\n", t_cce))
  results <- rbind(results, data.frame(
    dataset = spec$name, N = length(unique(d_i$id)),
    T = round(nrow(d_i) / length(unique(d_i$id))),
    estimator = "CCE-MG", tool = "R/plm", time_sec = t_cce
  ))

  cat("\n")
}

# =============================================================
# PART 3: Save results for blog post table
# =============================================================

write.csv(results, "validation/benchmark_r_results.csv", row.names = FALSE)

cat("=============================================================\n")
cat("Results saved to validation/benchmark_r_results.csv\n")
cat("Run validation/benchmark_gauss.e to get GAUSS timings,\n")
cat("then merge with validation/compile_benchmark.R.\n")
cat("=============================================================\n")
