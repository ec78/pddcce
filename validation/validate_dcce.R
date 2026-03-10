# =============================================================================
# dccelib Validation Script
# Compares GAUSS dccelib output against R plm::pmg() for:
#   1. Mean Group (MG) Estimator
#   2. Common Correlated Effects MG (CCE-MG)
#
# Reference: Pesaran (2006), Econometrica 74(4), 967-1012.
# R package: plm >= 2.6 — pmg() with model="mg" and model="cmg"
# =============================================================================

library(haven)
library(plm)
library(dplyr)

# -----------------------------------------------------------------------------
# 1. Load and prepare data
# -----------------------------------------------------------------------------
d <- read_dta("examples/penn_sample.dta")

# Remove rows with any NA — equivalent to GAUSS packr()
d_clean <- na.omit(d)

cat("Obs after na.omit:", nrow(d_clean), "\n")
cat("Groups:", length(unique(d_clean$id)), "\n")
cat("Years:", min(d_clean$year), "-", max(d_clean$year), "\n\n")

# Create panel data frame (plm requires this)
pdata <- pdata.frame(d_clean, index = c("id", "year"))

# -----------------------------------------------------------------------------
# 2. Mean Group (MG) Estimator
#    GAUSS: mg(data[., "id" "year" "log_rgdpo" "log_ck" "log_ngd"])
#    R:     pmg(log_rgdpo ~ log_ck + log_ngd, data=pdata, model="mg")
# -----------------------------------------------------------------------------
cat("============================================================\n")
cat("MODEL 1: Mean Group (MG) Estimator\n")
cat("============================================================\n")

mg_fit <- pmg(log_rgdpo ~ log_ck + log_ngd, data = pdata, model = "mg")
mg_sum <- summary(mg_fit)
print(mg_sum)

cat("\n--- Extracted coefficients (for GAUSS comparison) ---\n")
mg_coef <- mg_sum$CoefTable
print(round(mg_coef, 6))

# -----------------------------------------------------------------------------
# 3. CCE Mean Group (CCE-MG) Estimator
#    GAUSS: cce_mg(data[., "id" "year" "log_rgdpo" "log_ck" "log_ngd"])
#    R:     pmg(log_rgdpo ~ log_ck + log_ngd, data=pdata, model="cmg")
#
#    Note: plm cmg adds cross-sectional averages of all model variables
#    (y, log_ck, log_ngd). The GAUSS default (no x_csa) does the same.
# -----------------------------------------------------------------------------
cat("\n============================================================\n")
cat("MODEL 2: CCE Mean Group (CCE-MG) Estimator\n")
cat("============================================================\n")

cce_fit <- pmg(log_rgdpo ~ log_ck + log_ngd, data = pdata, model = "cmg")
cce_sum <- summary(cce_fit)
print(cce_sum)

cat("\n--- Extracted coefficients (for GAUSS comparison) ---\n")
cce_coef <- cce_sum$CoefTable
print(round(cce_coef, 6))

# -----------------------------------------------------------------------------
# 4. DCCE-MG — manual construction
#    GAUSS: dcce_mg with y_lags=1, cr_lags=3, x_csa=log_hc
#
#    Steps:
#      a) Add lagged log_rgdpo as regressor
#      b) Compute cross-sectional averages of: log_rgdpo, log_ck, log_ngd, log_hc
#      c) Add PT=3 lags of each CSA series
#      d) Run pmg with model="mg" on the augmented dataset
#      e) Report only coefficients on log_rgdpo_lag, log_ck, log_ngd
# -----------------------------------------------------------------------------
cat("\n============================================================\n")
cat("MODEL 3: Dynamic CCE-MG (DCCE-MG) Estimator\n")
cat("  y_lags=1, cr_lags=3, x_csa=log_hc\n")
cat("============================================================\n")

# --- a) Lagged y (within each group) ---
d_dcce <- d_clean %>%
  arrange(id, year) %>%
  group_by(id) %>%
  mutate(log_rgdpo_l1 = lag(log_rgdpo, 1)) %>%
  ungroup() %>%
  filter(!is.na(log_rgdpo_l1))   # drop first obs per group (lag creates NA)

cat("Obs after lagging y:", nrow(d_dcce), "\n")

# --- b) Cross-sectional averages of log_rgdpo, log_ck, log_ngd, log_hc ---
#        (log_hc is the x_csa variable; others are direct regressors)
csa_vars <- c("log_rgdpo", "log_ck", "log_ngd", "log_hc")

d_csa <- d_dcce %>%
  group_by(year) %>%
  summarise(across(all_of(csa_vars), mean, .names = "{.col}_bar"), .groups = "drop")

d_dcce <- d_dcce %>%
  left_join(d_csa, by = "year")

# --- c) Three lags of each CSA variable ---
PT <- 3
csa_bar_names <- paste0(csa_vars, "_bar")

d_dcce <- d_dcce %>%
  arrange(id, year) %>%
  group_by(id) %>%
  mutate(across(all_of(csa_bar_names),
                list(l1 = ~lag(., 1), l2 = ~lag(., 2), l3 = ~lag(., 3)),
                .names = "{.col}_{.fn}")) %>%
  ungroup() %>%
  filter(!is.na(log_rgdpo_bar_l3))   # drop rows until all lags available

cat("Obs after CSA lags:", nrow(d_dcce), "\n\n")

# --- d) Run MG on augmented data ---
pdata_dcce <- pdata.frame(d_dcce, index = c("id", "year"))

# Build formula: y ~ y_l1 + x1 + x2 + all CSA levels + all CSA lags
csa_lag_terms <- paste(
  c(csa_bar_names,
    outer(csa_bar_names, paste0("_l", 1:PT), paste0) |> as.vector()),
  collapse = " + "
)
dcce_formula <- as.formula(
  paste("log_rgdpo ~ log_rgdpo_l1 + log_ck + log_ngd +", csa_lag_terms)
)

dcce_fit <- pmg(dcce_formula, data = pdata_dcce, model = "mg")
dcce_sum <- summary(dcce_fit)

cat("--- All coefficients ---\n")
print(round(coef(dcce_sum), 6))

cat("\n--- Key coefficients (log_rgdpo_l1, log_ck, log_ngd) ---\n")
key_vars <- c("log_rgdpo_l1", "log_ck", "log_ngd")
dcce_table <- dcce_sum$CoefTable
dcce_key <- dcce_table[key_vars, , drop = FALSE]
print(round(dcce_key, 6))

# -----------------------------------------------------------------------------
# 5. Summary table for direct comparison with GAUSS output
# -----------------------------------------------------------------------------
cat("\n============================================================\n")
cat("SUMMARY TABLE — for comparison with GAUSS dccelib output\n")
cat("============================================================\n")

fmt_coef_table <- function(label, cf_obj) {
  cat(sprintf("\n[ %s ]\n", label))
  cat(sprintf("  %-22s  %9s   %9s   %8s\n", "Variable", "Coef", "SE", "t-value"))
  cat(sprintf("  %s\n", strrep("-", 60)))
  if (is.matrix(cf_obj)) {
    tstat_col <- if ("t-value" %in% colnames(cf_obj)) "t-value" else "z-value"
    for (v in rownames(cf_obj)) {
      cat(sprintf("  %-22s  %9.5f   %9.5f   %8.4f\n",
                  v, cf_obj[v, "Estimate"], cf_obj[v, "Std. Error"], cf_obj[v, tstat_col]))
    }
  } else {
    for (v in names(cf_obj)) {
      cat(sprintf("  %-22s  %9.5f\n", v, cf_obj[v]))
    }
  }
}

fmt_coef_table("MG", mg_coef)
fmt_coef_table("CCE-MG", cce_coef)
fmt_coef_table("DCCE-MG (key vars: y_lag, log_ck, log_ngd)", dcce_key)

cat("\n============================================================\n")
cat("R session info:\n")
cat("plm version:", as.character(packageVersion("plm")), "\n")
cat("R version:", R.version$version.string, "\n")
