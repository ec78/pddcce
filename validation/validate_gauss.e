/*
** validate_gauss.e
**
** Runs MG, CCE-MG, and DCCE-MG on the Penn World Tables sample data
** and writes results to validation/gauss_results.txt for comparison
** against the R reference values in validate_dcce.R.
**
** Run from repo root with:
**   C:\gauss26\tgauss.exe -b -nj validation\validate_gauss.e
*/

// Load source directly — no library install needed
#include "src/cce.sdf"
#include "src/dcceutil.src"
#include "src/cce_mg.src"

// Redirect all output to file
output file = "validation/gauss_results.txt" RESET;

print "============================================================";
print "dccelib Validation — GAUSS 26";
print "Penn World Tables sample data";
print "============================================================";
print "";

// Load data
fname = __FILE_DIR $+ "../examples/penn_sample.dta";
data = packr(loadd(fname, ". + date($year, '%Y')"));
data = order(data, "id"$|"year");

local fmt1, fmt2, fmt3;
fmt1 = "%-20s %15.6f %15.6f %15.6f";
fmt2 = "%-20s %15s %15s %15s";

proc (0) = print_coef_table(mgO);
    local fmt_h, fmt_r, i, k;
    fmt_h = "%-20s %15s %15s %15s %15s";
    fmt_r = "%-20s %15.6f %15.6f %15.6f %15.6f";
    k = rows(mgO.b_mg);
    sprintf(fmt_h, "Variable", "Coef", "SE(NP)", "t(NP)", "p-value");
    print chrs(61*ones(72, 1));
    for i(1, k, 1);
        sprintf(fmt_r,
            mgO.mg_vars[i],
            mgO.b_mg[i],
            mgO.se_mg[i],
            mgO.tvalue[i],
            mgO.pval[i]);
    endfor;
    print chrs(61*ones(72, 1));
    sprintf("CD test stat: %10.4f   p-value: %10.4f", mgO.cd_stat, mgO.cd_pval);
    print "";
endp;

// -------------------------------------------------------------------
// MODEL 1: Mean Group (MG)
// -------------------------------------------------------------------
print "============================================================";
print "MODEL 1: Mean Group (MG) Estimator";
print "  Formula: log_rgdpo ~ log_ck + log_ngd";
print "============================================================";
print "";

reg_data = data[., "id" "year" "log_rgdpo" "log_ck" "log_ngd"];

mgCtl = mgControlCreate();
mgCtl.report = 0;    // suppress auto-print; we print manually below

mgO = mg(reg_data, mgCtl);

print_coef_table(mgO);

// -------------------------------------------------------------------
// MODEL 2: CCE Mean Group (CCE-MG)
// -------------------------------------------------------------------
print "============================================================";
print "MODEL 2: CCE Mean Group (CCE-MG) Estimator";
print "  Formula: log_rgdpo ~ log_ck + log_ngd";
print "  CSA: ybar, log_ck_bar, log_ngd_bar (default — no x_csa)";
print "============================================================";
print "";

cceCtl = mgControlCreate();
cceCtl.report = 0;

cceO = cce_mg(reg_data, cceCtl);

print_coef_table(cceO);

// -------------------------------------------------------------------
// MODEL 3: DCCE Mean Group (DCCE-MG)
// -------------------------------------------------------------------
print "============================================================";
print "MODEL 3: Dynamic CCE Mean Group (DCCE-MG) Estimator";
print "  Formula: log_rgdpo ~ log_ck + log_ngd";
print "  y_lags=1, cr_lags=3, x_csa=log_hc";
print "============================================================";
print "";

dcceCtl = mgControlCreate();
dcceCtl.y_lags  = 1;
dcceCtl.cr_lags = 3;
dcceCtl.x_csa   = data[., "log_hc"];
dcceCtl.report  = 0;

dcceO = dcce_mg(reg_data, dcceCtl);

print_coef_table(dcceO);

print "============================================================";
print "Validation complete.";
print "============================================================";

output off;
end;
