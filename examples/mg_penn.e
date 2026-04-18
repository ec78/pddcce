/*
** mg_penn.e
**
** Mean Group estimator on Penn World Tables data.
** Demonstrates three equivalent ways to specify variables:
**
**   1. Pre-select columns (traditional column-order approach)
**   2. Formula string as second argument: mg(data, "y ~ x1 + x2")
**   3. ctl.formula + ctl.groupvar/timevar (works even when year is not col 2)
**
** Data: Penn World Tables (N=93 countries, T~50 years)
** Model: log_rgdpo ~ log_ck + log_ngd
*/

new;
library dccelib;

fname = __FILE_DIR $+ "penn_world.dta";
data  = packr(loadd(fname, ". + date($year, '%Y')"));

// -----------------------------------------------------------------------
// Approach 1: Pre-select columns in standard [group, time, y, x] order
// -----------------------------------------------------------------------
print "=================================================================";
print "Approach 1: Column pre-selection (traditional)";
print "=================================================================";

reg_data = data[., "id" "year" "log_rgdpo" "log_ck" "log_ngd"];
mgO = mg(reg_data);

// -----------------------------------------------------------------------
// Approach 2: Formula string as second argument
// Uses pre-selected data so group/time are in columns 1/2
// -----------------------------------------------------------------------
print "=================================================================";
print "Approach 2: Formula string as second argument";
print "=================================================================";

mgO2 = mg(reg_data, "log_rgdpo ~ log_ck + log_ngd");

// -----------------------------------------------------------------------
// Approach 3: ctl.formula with explicit groupvar/timevar
// Works on the full unsorted data (year is column 6 in penn_world.dta)
// -----------------------------------------------------------------------
print "=================================================================";
print "Approach 3: ctl.formula + ctl.groupvar/timevar";
print "=================================================================";

ctl = mgControlCreate();
ctl.formula  = "log_rgdpo ~ log_ck + log_ngd";
ctl.groupvar = "id";
ctl.timevar  = "year";

mgO3 = mg(data, ctl);

// -----------------------------------------------------------------------
// Confirm all three give the same coefficients
// -----------------------------------------------------------------------
print "=================================================================";
print "Coefficient check (all three should match)";
print "=================================================================";
print "            Traditional   Formula arg   ctl.formula";
for i(1, 2, 1);
    print mgO.mg_vars[i] $+
          "     " $+ ntos(mgO.b_mg[i], 6) $+
          "   "   $+ ntos(mgO2.b_mg[i], 6) $+
          "   "   $+ ntos(mgO3.b_mg[i], 6);
endfor;
