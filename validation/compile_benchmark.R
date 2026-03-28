# =============================================================================
# compile_benchmark.R
#
# Reads timing results from all three benchmark scripts and produces:
#   1. A console table showing GAUSS vs R speed ratios
#   2. validation/benchmark_summary.md  — blog-ready Markdown table
#   3. validation/coverage_table.md     — estimator coverage comparison
#
# Run after all three benchmark scripts have completed:
#   Rscript validation/compile_benchmark.R
# =============================================================================

library(dplyr)
library(tidyr)

# ---------------------------------------------------------------------------
# 1. Load results
# ---------------------------------------------------------------------------
gauss_file  <- "validation/benchmark_gauss_results.csv"
r_file      <- "validation/benchmark_r_results.csv"
python_file <- "validation/benchmark_python_results.csv"

load_if_exists <- function(f) {
  if (file.exists(f)) {
    read.csv(f, stringsAsFactors = FALSE)
  } else {
    cat(sprintf("  Warning: %s not found — run the corresponding benchmark first.\n", f))
    NULL
  }
}

df_gauss  <- load_if_exists(gauss_file)
df_r      <- load_if_exists(r_file)
df_python <- load_if_exists(python_file)

timing_parts <- Filter(Negate(is.null), list(df_gauss, df_r, df_python))
if (length(timing_parts) == 0) stop("No benchmark results found.")

timing <- bind_rows(timing_parts)

# ---------------------------------------------------------------------------
# 2. Timing comparison table
# ---------------------------------------------------------------------------
cat("=============================================================\n")
cat("dccelib Benchmark Results — Speed Comparison\n")
cat("=============================================================\n\n")

# Pivot: rows = dataset x estimator, cols = tool times
pivot <- timing %>%
  group_by(dataset, N, T, estimator, tool) %>%
  summarise(time_sec = median(time_sec), .groups = "drop") %>%
  pivot_wider(names_from = tool, values_from = time_sec)

print(as.data.frame(pivot))

# ---------------------------------------------------------------------------
# 3. Speed ratio table (GAUSS / R)
# ---------------------------------------------------------------------------
if (!is.null(df_gauss) && !is.null(df_r)) {
  cat("\n--- Speed ratio: GAUSS dccelib / R plm (< 1 = GAUSS is faster) ---\n\n")

  ratio_tbl <- inner_join(
    df_gauss %>% rename(gauss_sec = time_sec) %>% select(-tool),
    df_r     %>% rename(r_sec     = time_sec) %>% select(-tool),
    by = c("dataset", "N", "T", "estimator")
  ) %>%
    mutate(ratio = round(gauss_sec / r_sec, 3)) %>%
    select(dataset, N, T, estimator, gauss_sec, r_sec, ratio)

  print(as.data.frame(ratio_tbl))
}

# ---------------------------------------------------------------------------
# 4. Write blog-ready Markdown timing table
# ---------------------------------------------------------------------------
summary_md <- "validation/benchmark_summary.md"

sink(summary_md)
cat("## dccelib Performance Benchmarks\n\n")
cat("Wall-clock time (seconds, median of 10 runs).\n\n")

# Build wide table
wide <- timing %>%
  group_by(dataset, N, T, estimator, tool) %>%
  summarise(time_sec = round(median(time_sec), 4), .groups = "drop") %>%
  pivot_wider(names_from = tool, values_from = time_sec,
              values_fn = ~ round(., 4))

# Sort columns: put GAUSS first
tool_cols <- setdiff(names(wide), c("dataset", "N", "T", "estimator"))
gauss_cols  <- sort(grep("GAUSS", tool_cols, value = TRUE))
other_cols  <- sort(setdiff(tool_cols, gauss_cols))
wide        <- wide[, c("dataset", "N", "T", "estimator", gauss_cols, other_cols)]

# Add ratio column if both exist
if (length(gauss_cols) == 1 && length(other_cols) >= 1) {
  r_col <- grep("R/plm", other_cols, value = TRUE)
  if (length(r_col) == 1) {
    wide$`GAUSS/R ratio` <- round(wide[[gauss_cols]] / wide[[r_col]], 2)
  }
}

# Markdown table
header <- paste(names(wide), collapse = " | ")
sep    <- paste(rep("---", ncol(wide)), collapse = " | ")
cat(header, "\n", sep, "\n", sep = "")
for (i in seq_len(nrow(wide))) {
  row <- wide[i, ]
  row_chr <- as.data.frame(lapply(row, as.character), stringsAsFactors = FALSE)
  row_chr[is.na(row_chr)] <- "—"
  cat(paste(unlist(row_chr), collapse = " | "), "\n")
}

cat("\n_GAUSS/R ratio < 1 indicates GAUSS dccelib is faster._\n")
sink()

cat(sprintf("\nBlog timing table written to %s\n", summary_md))


# ---------------------------------------------------------------------------
# 5. Write estimator coverage table (Markdown)
# ---------------------------------------------------------------------------
coverage_md <- "validation/coverage_table.md"

sink(coverage_md)
cat("## Panel Estimator Coverage: dccelib vs R plm vs Python\n\n")
cat("| Estimator / Feature | GAUSS dccelib | R plm | Python linearmodels |\n")
cat("|---------------------|:-------------:|:-----:|:-------------------:|\n")

capabilities <- list(
  c("Mean Group (MG)",                    "Yes", "Yes (`pmg`)", "No"),
  c("CCE Mean Group (CCE-MG)",            "Yes", "Yes (`pmg`)", "No"),
  c("Dynamic CCE-MG (DCCE-MG)",           "Yes", "No",         "No"),
  c("Pooled CCE with NW SE",              "Yes", "Yes (`pcce`)", "No"),
  c("PC-CCE-MG (PCA augmentation)",       "Yes (GML req.)", "No", "No"),
  c("HPJ bias correction",                "Yes", "No",         "No"),
  c("Wild bootstrap SEs",                 "Yes", "No",         "No"),
  c("CIPS panel unit root (Pesaran 2007)","Yes", "No (CD.test only)", "No"),
  c("Pesaran-Yamagata slope homogeneity", "Yes", "No",         "No"),
  c("Westerlund cointegration test",      "Yes", "No (in `westerlund` pkg)", "No"),
  c("CCE rank condition test",            "Yes", "No",         "No"),
  c("Long-run multipliers (delta method)","Yes", "No",         "No"),
  c("Residual / coefficient plots",       "Yes", "No",         "No"),
  c("LaTeX table export",                 "Yes", "No",         "No"),
  c("I(1) extension (KPY 2011)",          "Yes", "No",         "No"),
  c("Two-way CCE (Bai 2009)",             "Yes", "No",         "No")
)

for (row in capabilities) {
  cat(sprintf("| %-41s | %-13s | %-5s | %-19s |\n",
              row[1], row[2], row[3], row[4]))
}

cat("\n_Table reflects package capabilities as of 2026. ")
cat("plm CCE-MG (`pmg` model=\"cmg\") uses the Pesaran (2006) estimator ")
cat("without dynamic extension or auxiliary diagnostics._\n")
sink()

cat(sprintf("Coverage table written to %s\n", coverage_md))
cat("\nAll done.\n")
