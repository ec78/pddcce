/*
** advanced_cce.e
**
** Demonstrates three advanced CCE options available through the mgControl struct:
**
**   1. Pooled CCE (ctl.pooled = 1)
**      Runs pooled CCE with Newey-West HAC SE alongside the MG estimator.
**      Use when slope homogeneity cannot be rejected; more efficient than MG.
**
**   2. I(1) Extension (ctl.i1 = 1)
**      Adds first-differenced cross-sectional averages alongside levels
**      per Kapetanios, Pesaran and Yamagata (2011). Use when CIPS tests
**      indicate the series are I(1).
**
**   3. Two-Way Factor Structure (ctl.two_way = 1)
**      Time-demeans all data before CCE augmentation per Bai (2009).
**      Addresses panels with both unit-specific and time-specific
**      factor structures (e.g., global shocks + country trends).
**
** References:
**   Pesaran, M.H. (2006). Estimation and inference in large heterogeneous
**     panels with a multifactor error structure. Econometrica, 74(4), 967-1012.
**   Kapetanios, G., Pesaran, M.H. and Yamagata, T. (2011). Panels with
**     non-stationary multifactor error structures. Journal of Econometrics,
**     160(2), 326-348.
**   Bai, J. (2009). Panel data models with interactive fixed effects.
**     Econometrica, 77(4), 1229-1279.
*/

new;
library dccelib;

// -----------------------------------------------------------------------
// Load data
// -----------------------------------------------------------------------
fname = __FILE_DIR $+ "penn_sample.dta";
data  = packr(loadd(fname, ". + date($year, '%Y')"));
data  = order(data, "id"$|"year");

reg_data = data[., "id" "year" "log_rgdpo" "log_ck" "log_ngd"];


// -----------------------------------------------------------------------
// Option 1: Pooled CCE (ctl.pooled = 1)
// -----------------------------------------------------------------------
// Setting pooled=1 estimates pooled CCE with Newey-West HAC standard errors
// in the same call as CCE-MG. Results are stored in mgO.pcce.
//
// Use pooled CCE when:
//   - Slope homogeneity is not rejected (slopehomo test)
//   - You want a more efficient estimator under homogeneity
//   - You want to compare MG and pooled estimates side by side
// -----------------------------------------------------------------------
print "=================================================================";
print "OPTION 1: Pooled CCE (NW SE) alongside CCE-MG";
print "=================================================================";

struct mgControl ctl1;
ctl1 = mgControlCreate();
ctl1.x_csa  = data[., "log_hc"];
ctl1.pooled = 1;    // also run pooled CCE in the same call

struct mgOut cceO_pooled;
cceO_pooled = cce_mg(reg_data, ctl1);

// MG estimates are in cceO_pooled (printed automatically above)
// Pooled CCE results are in the embedded pcce struct:
print;
print "--- Pooled CCE (NW) results ---";
print cceO_pooled.pcce.out_pcce;
print;
print "Note: pooled CCE is inconsistent if slope heterogeneity is present.";
print "Run slopehomo() first to test H0: slopes are equal across units.";
print;


// -----------------------------------------------------------------------
// Option 2: I(1) Extension — KPY (2011)
// -----------------------------------------------------------------------
// When the common factors are integrated of order one, cross-sectional
// averages in levels do not fully proxy for the I(1) factors. The KPY
// extension adds Dbar_y_t and Dbar_x_t (first differences of CSAs)
// alongside the levels to restore CCE consistency.
//
// Use this when:
//   - CIPS tests do not reject unit roots in the levels of y or x
//   - You want robustness to I(1) factor structure
// -----------------------------------------------------------------------
print "=================================================================";
print "OPTION 2: I(1) Extension — KPY (2011)  [ctl.i1 = 1]";
print "=================================================================";

struct mgControl ctl2;
ctl2 = mgControlCreate();
ctl2.x_csa = data[., "log_hc"];
ctl2.i1    = 1;    // add first-differenced CSAs to h matrix

struct mgOut cceO_i1;
cceO_i1 = cce_mg(reg_data, ctl2);

print;
print "The model description reports '[KPY I(1) extension]'.";
print "The CSA variable list (csa_vars) includes the differenced CSAs.";
print;


// -----------------------------------------------------------------------
// Option 3: Two-Way Factor Structure — Bai (2009)
// -----------------------------------------------------------------------
// Standard CCE handles cross-sectional factor structure via CSA augmentation.
// The two-way extension additionally time-demeans all variables before
// computing the cross-sectional averages, absorbing time-varying aggregate
// shocks on top of the cross-sectional factor structure.
//
// Use this when:
//   - The panel has both unit-specific factor loadings AND time-specific
//     aggregate shocks
//   - A two-way interactive fixed effects structure is suspected
// -----------------------------------------------------------------------
print "=================================================================";
print "OPTION 3: Two-Way CCE  [ctl.two_way = 1]";
print "=================================================================";

struct mgControl ctl3;
ctl3 = mgControlCreate();
ctl3.x_csa   = data[., "log_hc"];
ctl3.two_way = 1;    // time-demean y, x, x_csa before CCE augmentation

struct mgOut cceO_2way;
cceO_2way = cce_mg(reg_data, ctl3);

print;
print "The model description reports '[Two-Way CCE]'.";
print;


// -----------------------------------------------------------------------
// Combining options: DCCE-MG with I(1) + pooled
// -----------------------------------------------------------------------
// Options can be combined. This example runs dynamic CCE with the I(1)
// extension and also estimates pooled CCE in the same call.
// -----------------------------------------------------------------------
print "=================================================================";
print "COMBINED: DCCE-MG with I(1) extension + pooled CCE";
print "=================================================================";

struct mgControl ctl4;
ctl4 = mgControlCreate();
ctl4.y_lags  = 1;
ctl4.cr_lags = 3;
ctl4.x_csa   = data[., "log_hc"];
ctl4.i1      = 1;    // I(1) extension
ctl4.pooled  = 1;    // also estimate pooled CCE

struct mgOut dcceO_full;
dcceO_full = dcce_mg(reg_data, ctl4);

print;
print "MG coefficients:";
print dcceO_full.b_mg;
print;
print "Pooled CCE coefficients:";
print dcceO_full.pcce.b_pcce;
