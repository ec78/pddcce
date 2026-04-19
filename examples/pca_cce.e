/*
** pca_cce.e
**
** Demonstrates the Principal Component CCE Mean Group (PC-CCE-MG) estimator.
**
** Standard CCE projects out common factors using raw cross-sectional averages
** of y and x. PC-CCE replaces those CSAs with their leading principal
** components, which is preferable when:
**   - The CSA matrix is close to rank-deficient (many near-collinear averages)
**   - The number of latent factors is small relative to the number of variables
**   - Better power is desired by concentrating factor information
**
** The number of PCs can be chosen automatically via the eigenvalue-ratio
** criterion of Ahn & Horenstein (2013), or fixed by the user.
**
** Data: Penn World Tables (N=93 countries, T~50 years)
** Model: log_rgdpo ~ log_ck + log_ngd, extra CSA: log_hc
*/

new;
library dccelib;

// -----------------------------------------------------------------------
// Load data
// -----------------------------------------------------------------------
fname = __FILE_DIR $+ "penn_world.dta";
data  = packr(loadd(fname, ". + date($year, '%Y')"));
data  = order(data, "id"$|"year");

// reg_data: standard column order [id, year, y, x1, x2]
reg_data = data[., "id" "year" "log_rgdpo" "log_ck" "log_ngd"];

// -----------------------------------------------------------------------
// PART 1: Standard CCE-MG (baseline)
// Uses ctl.formula + ctl.x_csa_names on the full data (new API)
// -----------------------------------------------------------------------
print "=================================================================";
print "BASELINE: Standard CCE-MG";
print "=================================================================";

ctl = mgControlCreate();
ctl.formula = "log_rgdpo ~ log_ck + log_ngd + csa(log_hc)";

struct mgOut cceO;
cceO = cce_mg(data, ctl);

// -----------------------------------------------------------------------
// PART 2: PC-CCE with automatic PC selection (Ahn-Horenstein ER criterion)
// -----------------------------------------------------------------------
print "=================================================================";
print "PC-CCE: Automatic PC selection (Ahn-Horenstein 2013)";
print "=================================================================";

ctl_pc = mgControlCreate();
ctl_pc.formula = "log_rgdpo ~ log_ck + log_ngd + csa(log_hc)";
ctl_pc.report  = 1;

struct mgOut pcceO;
pcceO = pcce_mg(data, ctl_pc);    // num_pc=0 → automatic selection

print "Model: " $+ pcceO.model;
print "";

// -----------------------------------------------------------------------
// PART 3: PC-CCE with fixed number of PCs
// -----------------------------------------------------------------------
print "=================================================================";
print "PC-CCE: Fixed at 2 principal components";
print "=================================================================";

struct mgOut pcceO_2pc;
pcceO_2pc = pcce_mg(data, ctl_pc, 2);

print "Model: " $+ pcceO_2pc.model;
print "";

// -----------------------------------------------------------------------
// PART 4: Formula-string API via positional arg
// -----------------------------------------------------------------------
print "=================================================================";
print "PC-CCE via positional formula string";
print "=================================================================";

ctl_f = mgControlCreate();
ctl_f.report = 1;

struct mgOut pcceO_f;
pcceO_f = pcce_mg(data, "log_rgdpo ~ log_ck + log_ngd + csa(log_hc)", ctl_f);

print "Model: " $+ pcceO_f.model;
print "";

// -----------------------------------------------------------------------
// PART 5: Rank check — compare CSA rank before and after PCA
// -----------------------------------------------------------------------
print "=================================================================";
print "Rank condition check (raw CSA)";
print "=================================================================";

print "--- Raw CSA rank ---";
struct cceRankOut rankO;
rankO = cce_rank(reg_data);
print_cce_rank(rankO);

print "";
print "If rank < k+1, CCE is inconsistent; PC-CCE resolves rank deficiency.";
print "";

// -----------------------------------------------------------------------
// PART 6: Coefficient comparison
// -----------------------------------------------------------------------
print "=================================================================";
print "Coefficient comparison: CCE-MG vs. PC-CCE";
print "=================================================================";

printCoefCompare(cceO.mg_vars[1:2],
    cceO.b_mg[1:2] ~ pcceO.b_mg[1:2],
    "CCE-MG" $| "PC-CCE (auto)");
