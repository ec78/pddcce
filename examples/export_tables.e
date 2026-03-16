/*
** export_tables.e
**
** Demonstrates publication-ready LaTeX table export from dccelib:
**
**   1. Single-model table — mgOutToLatex()
**      Exports one mgOut result to a .tex file with significance stars,
**      CD stat footer, and mean R-squared footer.
**
**   2. Multi-model comparison table — mgOutToLatexMulti()
**      Places MG, CCE-MG, and DCCE-MG side by side in one table.
**      Ideal for the main results table in empirical papers.
**
**   3. Numeric coefficient matrix — coeftable()
**      Extracts [coef, se, t-stat, p-value] as a plain matrix for
**      downstream computations (Wald tests, coefficient plots, etc.).
**
** Output files are written to the same directory as this script.
** Include in LaTeX with: \input{filename.tex}
*/

new;
library dccelib;

// -----------------------------------------------------------------------
// Load data and run three estimators
// -----------------------------------------------------------------------
fname = __FILE_DIR $+ "penn_world.dta";
data  = packr(loadd(fname, ". + date($year, '%Y')"));
data  = order(data, "id"$|"year");

reg_data = data[., "id" "year" "log_rgdpo" "log_ck" "log_ngd"];

// MG
print "Estimating MG...";
struct mgOut mgO;
mgO = mg(reg_data);

// CCE-MG
print "Estimating CCE-MG...";
struct mgControl ctl;
ctl = mgControlCreate();
ctl.x_csa  = data[., "log_hc"];
ctl.report = 0;

struct mgOut cceO;
cceO = cce_mg(reg_data, ctl);

// DCCE-MG
print "Estimating DCCE-MG...";
ctl.y_lags  = 1;
ctl.cr_lags = 3;

struct mgOut dcceO;
dcceO = dcce_mg(reg_data, ctl);
print;


// -----------------------------------------------------------------------
// Part 1: Single-model LaTeX table
// -----------------------------------------------------------------------
// mgOutToLatex(mgO, filename [, se_type, note])
//
// Arguments:
//   mgO       -- mgOut struct from any estimator
//   filename  -- output .tex path (absolute or relative to CWD)
//   se_type   -- optional: "np" (default, NP SE), "nw" (Newey-West SE)
//   note      -- optional: string appended as a table note
//
// The table includes:
//   - Variable names, coefficients, standard errors (in parentheses)
//   - Significance stars: *** p<0.01, ** p<0.05, * p<0.10
//   - Footer row: CD statistic and p-value
//   - Footer row: Mean R-squared across groups
// -----------------------------------------------------------------------
print "=================================================================";
print "PART 1: Single-model LaTeX table (CCE-MG)";
print "=================================================================";

out_file_cce = __FILE_DIR $+ "cce_table.tex";

mgOutToLatex(cceO, out_file_cce, "np",
             "Penn World Tables. N=93 countries, T\\approx 50 years. " $+
             "NP standard errors in parentheses.");

print "Written to: " $+ out_file_cce;
print;


// -----------------------------------------------------------------------
// Part 2: Multi-model comparison table
// -----------------------------------------------------------------------
// mgOutToLatexMulti(mgO_arr, labels, filename [, note])
//
// Arguments:
//   mgO_arr  -- array of mgOut structs (up to 6 models)
//   labels   -- string array of column headers
//   filename -- output .tex path
//   note     -- optional: table note appended at the bottom
//
// The table places each model in its own column with shared row labels.
// Rows from different models are aligned by variable name.
// -----------------------------------------------------------------------
print "=================================================================";
print "PART 2: Multi-model comparison table (MG vs CCE-MG vs DCCE-MG)";
print "=================================================================";

out_file_multi = __FILE_DIR $+ "comparison_table.tex";

labels = "MG"$|"CCE-MG"$|"DCCE-MG";

// Pass models as a vertical array using the | operator
mgOutToLatexMulti(mgO|cceO|dcceO, labels, out_file_multi,
                  "Dependent variable: $\\log(RGDP_o)$. " $+
                  "Penn World Tables (N=93, T\\approx 50). " $+
                  "NP standard errors in parentheses. " $+
                  "Significance: $^{***}$p$<$0.01, $^{**}$p$<$0.05, $^{*}$p$<$0.10.");

print "Written to: " $+ out_file_multi;
print;
print "To include in LaTeX:";
print "  \\input{comparison_table.tex}";
print;


// -----------------------------------------------------------------------
// Part 3: Numeric coefficient matrix (coeftable)
// -----------------------------------------------------------------------
// coeftable(mgO) returns a k x 4 numeric matrix:
//   col 1: coefficients
//   col 2: standard errors (NP)
//   col 3: t-statistics
//   col 4: p-values
//
// Use this for:
//   - Programmatic access to results (no string parsing)
//   - Computing Wald statistics: W = R*b'*inv(R*V*R')*R*b
//   - Constructing coefficient plots
//   - Comparing estimates across bootstrap iterations
// -----------------------------------------------------------------------
print "=================================================================";
print "PART 3: Numeric coefficient matrix (coeftable)";
print "=================================================================";

ct = coeftable(cceO);

print;
print "coeftable(cceO) returns a " $+ ntos(rows(ct), 1) $+
      " x 4 matrix: [coef, se, t-stat, p-value]";
print ct;

print;
print "Example: extract just the capital coefficient (first row):";
print "  coef = ct[1, 1];   // " $+ ntos(ct[1, 1], 6);
print "  se   = ct[1, 2];   // " $+ ntos(ct[1, 2], 6);
print "  tval = ct[1, 3];   // " $+ ntos(ct[1, 3], 6);
print "  pval = ct[1, 4];   // " $+ ntos(ct[1, 4], 6);

// -----------------------------------------------------------------------
// Example: Wald test of joint significance (H0: b_ck = b_ngd = 0)
// Using the MG covariance matrix
// -----------------------------------------------------------------------
print;
print "=================================================================";
print "BONUS: Joint Wald test using mgO.cov_mg";
print "=================================================================";

// H0: log_ck = 0 and log_ngd = 0
// Restriction matrix R selects the first two slope coefficients
// (excluding the intercept which is stored last)
k = rows(cceO.b_mg);

// b_mg order: [slope vars..., intercept]
// R picks out the two slope coefficients
R = eye(k-1) ~ zeros(k-1, 1);    // k-1 x k restriction matrix (exclude intercept)
b = cceO.b_mg;
V = cceO.cov_mg;

W = b' * R' * invpd(R * V * R') * R * b;
df_wald = k - 1;
pval_wald = cdfchic(W, df_wald);

print;
print "Wald statistic (H0: all slopes = 0): " $+ ntos(W, 4);
print "Degrees of freedom:                  " $+ ntos(df_wald, 1);
print "p-value:                             " $+ ntos(pval_wald, 4);
