/*
** full_workflow.e
**
** Complete empirical workflow for panel estimation with cross-sectional
** dependence, following the recommended dccelib analysis sequence:
**
**   1. Load and prepare data
**   2. Baseline MG detect cross-sectional dependence
**   3. Panel unit root tests (CIPS)
**   4. CCE-MG absorb common factors
**   5. Slope homogeneity test
**   6. Dynamic CCE-MG add persistence
**   7. HPJ bias correction
**   8. Export results to LaTeX
**
** Data: Penn World Tables (N=93 countries, T~50 years)
** Model: log_rgdpo ~ log_ck + log_ngd, extra CSA: log_hc
** Variable specification: formula string throughout
*/

new;
library dccelib;

// -----------------------------------------------------------------------
// 1. Load and prepare data
// -----------------------------------------------------------------------
fname = __FILE_DIR $+ "penn_world.dta";
data  = packr(loadd(fname, ". + date($year, '%Y')"));
data  = order(data, "id"$|"year");

// Base control struct — csa(log_hc) specifies log_hc as a CSA-only variable
// directly inside the formula string, no separate x_csa_names needed.
ctl = mgControlCreate();
ctl.formula = "log_rgdpo ~ log_ck + log_ngd + csa(log_hc)";

// -----------------------------------------------------------------------
// 2. Baseline MG does cross-sectional dependence exist?
// -----------------------------------------------------------------------
print "=================================================================";
print "STEP 1: Mean Group (baseline no CD correction)";
print "=================================================================";

struct mgOut mgO;
mgO = mg(data, ctl);

print "";
print "CD statistic on MG residuals: " $+ ntos(mgO.cd_stat, 4);
print "CD p-value:                   " $+ ntos(mgO.cd_pval, 4);
print "";
print "A large CD statistic (|CD| >> 1.96) indicates cross-sectional";
print "dependence. CCE augmentation is warranted.";
print "";

// Residual diagnostics on plain MG
// Look for: non-normal histogram, heavy Q-Q tails, and high per-group
// SD variation all signs that common factors contaminate the residuals.
plotResiduals(mgO);

// -----------------------------------------------------------------------
// 3. Panel unit root tests (CIPS)
// -----------------------------------------------------------------------
print "=================================================================";
print "STEP 2: CIPS Panel Unit Root Tests";
print "=================================================================";

print "--- log_rgdpo ---";
{ cips_y, cadf_y }   = cips(data, "log_rgdpo");

print "--- log_ck ---";
{ cips_ck, cadf_ck } = cips(data, "log_ck");

print "--- log_ngd ---";
{ cips_ngd, cadf_ngd } = cips(data, "log_ngd");

print "";
print "If series are I(1), use ctl.i1 = 1 (KPY 2011 extension).";
print "";

// -----------------------------------------------------------------------
// 4. CCE-MG absorb common factors
// -----------------------------------------------------------------------
print "=================================================================";
print "STEP 3: CCE Mean Group";
print "=================================================================";

struct mgOut cceO;
cceO = cce_mg(data, ctl);

print "";
print "CD statistic after CCE augmentation: " $+ ntos(cceO.cd_stat, 4);
print "CD p-value:                          " $+ ntos(cceO.cd_pval, 4);
print "";
print "The CD statistic should fall sharply relative to plain MG,";
print "confirming that CSA augmentation spanned the factor space.";
print "";

// Residual diagnostics after CCE augmentation
// Compare to the MG plot: residuals should be better centred, the Q-Q
// plot closer to the diagonal, and per-group SDs more uniform.
plotResiduals(cceO);

// Caterpillar plot of per-group slope estimates
// Wide spread of country-level coefficients around the MG mean confirms
// slope heterogeneity pooling across countries would be misleading.
plotCoefficients(cceO);

// Sample ACF of pooled CCE-MG residuals
// Significant autocorrelation at short lags suggests serial persistence
// in growth rates that the static model has not absorbed motivating
// the dynamic CCE-MG specification in the next step.
plotClearLayout();
plotResidualACF(cceO);

// -----------------------------------------------------------------------
// 5. Slope homogeneity test
// -----------------------------------------------------------------------
print "=================================================================";
print "STEP 4: Slope Homogeneity Test (Pesaran-Yamagata 2008)";
print "=================================================================";
{ delta, pval, delta_adj, pval_adj } = slopehomo(cceO);

print "";
if pval_adj < 0.05;
    print "Slope homogeneity rejected at 5%: MG estimator is appropriate.";
else;
    print "Cannot reject slope homogeneity: pooled CCE may be efficient (ctl.pooled=1).";
endif;
print "";

// -----------------------------------------------------------------------
// 6. Dynamic CCE-MG account for persistence in growth rates
// -----------------------------------------------------------------------
print "=================================================================";
print "STEP 5: Dynamic CCE Mean Group (y_lags=1)";
print "=================================================================";

ctl_d = mgControlCreate();
ctl_d.formula  = "log_rgdpo ~ log_ck + log_ngd + csa(log_hc)";
ctl_d.y_lags   = 1;
ctl_d.cr_lags  = 3;

struct mgOut dcceO;
dcceO = dcce_mg(data, ctl_d);

// Long-run multipliers from the dynamic model
print "=================================================================";
print "STEP 5b: Long-Run Multipliers";
print "=================================================================";

struct longRunOut lrO;
lrO = longRunMG(dcceO);

// -----------------------------------------------------------------------
// 7. HPJ bias correction (recommended when T is moderate)
// -----------------------------------------------------------------------
print "=================================================================";
print "STEP 6: Half-Panel Jackknife Bias Correction";
print "=================================================================";

struct mgOut hpjO;
hpjO = hpj(data, ctl_d, "dcce_mg");

print "";
print "HPJ-corrected vs. uncorrected DCCE-MG coefficients:";
printCoefCompare(dcceO.mg_vars[1:3],
    dcceO.b_mg[1:3] ~ hpjO.b_mg[1:3],
    "DCCE-MG" $| "HPJ");
print "";

// -----------------------------------------------------------------------
// 8. Coefficient comparison across all models
// -----------------------------------------------------------------------
print "=================================================================";
print "SUMMARY: Coefficient progression across estimators";
print "=================================================================";

printCoefCompare(cceO.mg_vars[1:2],
    mgO.b_mg[1:2] ~ cceO.b_mg[1:2] ~ dcceO.b_mg[1:2] ~ hpjO.b_mg[1:2],
    "MG" $| "CCE-MG" $| "DCCE-MG" $| "HPJ");

// -----------------------------------------------------------------------
// 9. Export results to LaTeX
// -----------------------------------------------------------------------
print "=================================================================";
print "STEP 7: Export to LaTeX";
print "=================================================================";

mgOutToLatexMulti(mgO|cceO|dcceO,
    "MG" $| "CCE-MG" $| "DCCE-MG",
    __FILE_DIR $+ "tables/workflow_results.tex",
    "Penn World Tables. NP standard errors in parentheses. " $+
    "*** p<0.01, ** p<0.05, * p<0.10.");

print "Results written to examples/tables/workflow_results.tex";
