/*
** cce_penn.e
**
** CCE Mean Group estimator on Penn World Tables data.
** Demonstrates three equivalent ways to specify variables:
**
**   1. Pre-select columns (traditional column-order approach)
**   2. ctl.y_var / ctl.x_vars  (select by column name via mgControl)
**   3. ctl.formula             (Wilkinson formula string via mgControl)
**
** In all three cases the estimator is identical.
**
** Data: Penn World Tables (N=93 countries, T~50 years)
** Model: log_rgdpo ~ log_ck + log_ngd, extra CSA: log_hc
*/

new;
library dccelib;

fname = __FILE_DIR $+ "penn_world.dta";
data  = packr(loadd(fname, ". + date($year, '%Y')"));

// -----------------------------------------------------------------------
// Approach 1: Pre-select columns (traditional)
// -----------------------------------------------------------------------
print "=================================================================";
print "Approach 1: Column pre-selection (traditional)";
print "=================================================================";

reg_data = data[., "id" "year" "log_rgdpo" "log_ck" "log_ngd"];

ctl1 = mgControlCreate();
ctl1.x_csa = data[., "log_hc"];    // extra CSA variable as matrix (traditional)

cceO1 = cce_mg(packr(reg_data), ctl1);

// -----------------------------------------------------------------------
// Approach 2: ctl.y_var / ctl.x_vars (select by column name)
// -----------------------------------------------------------------------
print "=================================================================";
print "Approach 2: ctl.y_var / ctl.x_vars (select by column name)";
print "=================================================================";

ctl2 = mgControlCreate();
ctl2.y_var       = "log_rgdpo";
ctl2.x_vars      = "log_ck" $| "log_ngd";
ctl2.x_csa_names = "log_hc";   // extra CSA variable by name (new)

cceO2 = cce_mg(packr(data), ctl2);

// -----------------------------------------------------------------------
// Approach 3: ctl.formula (Wilkinson formula string)
// -----------------------------------------------------------------------
print "=================================================================";
print "Approach 3: ctl.formula with separate x_csa_names";
print "=================================================================";

ctl3 = mgControlCreate();
ctl3.formula     = "log_rgdpo ~ log_ck + log_ngd";
ctl3.x_csa_names = "log_hc";

cceO3 = cce_mg(packr(data), ctl3);

// -----------------------------------------------------------------------
// Approach 4: csa() inline in formula string
// -----------------------------------------------------------------------
print "=================================================================";
print "Approach 4: csa() inline in formula string";
print "=================================================================";

// CSA variable specified directly inside the formula — no separate
// x_csa_names field needed. Equivalent to Approach 3.
cceO4 = cce_mg(packr(data), "log_rgdpo ~ log_ck + log_ngd + csa(log_hc)");

// -----------------------------------------------------------------------
// Confirm all four give the same coefficients
// -----------------------------------------------------------------------
print "=================================================================";
print "Coefficient check (all four should match)";
print "=================================================================";
printCoefCompare(cceO1.mg_vars[1:2],
    cceO1.b_mg[1:2] ~ cceO2.b_mg[1:2] ~ cceO3.b_mg[1:2] ~ cceO4.b_mg[1:2],
    "Approach 1" $| "Approach 2" $| "Approach 3" $| "Approach 4");
