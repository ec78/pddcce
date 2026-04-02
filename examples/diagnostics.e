/*
** diagnostics.e
**
** Demonstrates the diagnostic testing workflow recommended before CCE estimation:
**   1. Pesaran (2004) CD test on MG residuals (built into mg/cce_mg/dcce_mg output)
**   2. Pesaran (2007) CIPS panel unit root test for each series
**   3. Pesaran-Yamagata (2008) slope homogeneity test on CCE-MG estimates
**
** References:
**   Pesaran, M.H. (2004). General diagnostic tests for cross-section dependence in
**     panels. CESifo Working Paper No. 1229.
**   Pesaran, M.H. (2007). A simple panel unit root test in the presence of
**     cross-section dependence. Journal of Applied Econometrics, 22(2), 265-312.
**   Pesaran, M.H. and Yamagata, T. (2008). Testing slope homogeneity in large panels.
**     Journal of Econometrics, 142(1), 50-93.
**
** Data: Penn World Tables (N=93 countries, T~50 years)
**   log_rgdpo   Real GDP (output side, log)
**   log_ck      Physical capital share (log)
**   log_ngd     Population growth + break-even investment (log)
**   log_hc      Human capital index (log)
*/

new;
library dccelib;

// -----------------------------------------------------------------------
// Load data
// -----------------------------------------------------------------------
fname = __FILE_DIR $+ "penn_world.dta";
data  = packr(loadd(fname, ". + date($year, '%Y')"));
data  = order(data, "id"$|"year");

reg_data = data[., "id" "year" "log_rgdpo" "log_ck" "log_ngd"];


// -----------------------------------------------------------------------
// Step 1: Mean Group baseline — inspect the CD statistic
// -----------------------------------------------------------------------
// MG is consistent under slope heterogeneity but does NOT correct for
// cross-sectional dependence. Its CD stat tells us whether we need CCE.
// -----------------------------------------------------------------------
print "=================================================================";
print "STEP 1: MG Baseline — CD Diagnostic";
print "=================================================================";

mgO = mg(reg_data);

// CD stat is also printed in the output table, but access it directly:
print;
print "CD statistic on MG residuals: " $+ ntos(mgO.cd_stat, 4);
print "p-value:                       " $+ ntos(mgO.cd_pval, 4);
print;
print "Interpretation: if |CD| >> 2.0, cross-sectional dependence is";
print "present and MG standard errors are invalid. CCE correction needed.";
print;


// -----------------------------------------------------------------------
// Step 2: CIPS Panel Unit Root Test
// -----------------------------------------------------------------------
// Test each series for unit roots before choosing static vs. dynamic CCE.
// Pesaran (2007) CIPS extends IPS to the cross-sectionally dependent case.
//
// Pesaran (2007) critical values, N=100, T=50, no trend (Table 2b):
//   10%: -2.11,   5%: -2.20,   1%: -2.37
// -----------------------------------------------------------------------
print "=================================================================";
print "STEP 2: CIPS Panel Unit Root Tests (Pesaran 2007)";
print "=================================================================";

// Test each series — cips() prints results automatically (report=1 by default)
print "--- log_rgdpo (Real GDP, log) ---";
{ cips_stat, cadf_vec } = cips(data[., "id" "year" "log_rgdpo"], 1);

print "--- log_ck (Capital share, log) ---";
{ cips_stat, cadf_vec } = cips(data[., "id" "year" "log_ck"], 1);

print "--- log_ngd (Pop. growth + break-even, log) ---";
{ cips_stat, cadf_vec } = cips(data[., "id" "year" "log_ngd"], 1);

print "--- log_hc (Human capital, log) ---";
{ cips_stat, cadf_vec } = cips(data[., "id" "year" "log_hc"], 1);

print;
print "If series are I(1): use ctl.i1=1 for the KPY (2011) extension,";
print "which adds first-differenced CSAs alongside levels.";
print;


// -----------------------------------------------------------------------
// Step 3: Slope Homogeneity Test
// -----------------------------------------------------------------------
// Tests H0: beta_i = beta for all i (pooling is valid)
// vs.  H1: slopes are heterogeneous (use MG estimator)
//
// Run on CCE-MG estimates — slope homogeneity should be assessed
// after controlling for CD, not on raw FE residuals.
// -----------------------------------------------------------------------
print "=================================================================";
print "STEP 3: Slope Homogeneity Test (Pesaran-Yamagata 2008)";
print "=================================================================";

cceCtl = mgControlCreate();
cceCtl.x_csa   = data[., "log_hc"];
cceCtl.report  = 0;    // suppress estimator output for this step

cceO = cce_mg(reg_data, cceCtl);

// slopehomo() prints results automatically (report=1 by default)
{ delta, pval, delta_adj, pval_adj } = slopehomo(cceO);

print;
print "Interpretation:";
print "  Reject H0 => slope heterogeneity is present => use MG estimators";
print "  Fail to reject => pooled CCE may be efficient (use ctl.pooled=1)";
print;


// -----------------------------------------------------------------------
// Step 4: Verify CD is absorbed by CCE augmentation
// -----------------------------------------------------------------------
// Re-run CCE-MG with reporting and check the residual CD stat.
// A small CD stat after CCE confirms the CSAs span the factor space.
// -----------------------------------------------------------------------
print "=================================================================";
print "STEP 4: CCE-MG with CD check on CCE residuals";
print "=================================================================";

cceCtl.report = 1;
cceO = cce_mg(reg_data, cceCtl);

print;
print "CD on CCE-MG residuals: " $+ ntos(cceO.cd_stat, 4) $+
      "   (compare to CD on MG residuals: " $+ ntos(mgO.cd_stat, 4) $+ ")";
print;
print "If CD falls substantially after CCE, the augmentation is working.";
print "If CD remains large, try increasing cr_lags or adding to x_csa.";
