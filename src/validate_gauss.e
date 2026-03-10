/*
** validate_gauss.e
**
** Runs MG, CCE-MG, and DCCE-MG on the Penn World Tables sample data
** and writes results to validation/gauss_results.txt for comparison
** against the R reference values in validate_dcce.R.
**
** Must be run from the src/ directory (so #include cce.sdf resolves):
**   cd src
**   C:\gauss26\tgauss.exe -b -nj validate_gauss.e
*/

#include cce.sdf
#include dcceutil.src
#include cce_mg.src

output file = "../validation/gauss_results.txt" RESET;

print "============================================================";
print "dccelib Validation - GAUSS 26";
print "Penn World Tables sample data";
print "============================================================";
print "";

// Load data — path relative to src/
fname = "../examples/penn_sample.dta";
data = packr(loadd(fname, ". + date($year, '%Y')"));
data = order(data, "id"$|"year");

// Helper: print coefficient table from mgOut struct
proc (0) = print_coef_table(struct mgOut mgO);
    // out_mg is a GAUSS dataframe — print it directly
    print mgO.out_mg;
    sprintf("CD stat: %10.4f   CD p-value: %10.4f", mgO.cd_stat, mgO.cd_pval);
    print "";
endp;

// -------------------------------------------------------------------
// MODEL 1: Mean Group (MG)
// -------------------------------------------------------------------
print "============================================================";
print "MODEL 1: Mean Group (MG) Estimator";
print "  log_rgdpo ~ log_ck + log_ngd";
print "============================================================";

reg_data = data[., "id" "year" "log_rgdpo" "log_ck" "log_ngd"];

struct mgControl mgCtl;
mgCtl = mgControlCreate();
mgCtl.report = 0;

struct mgOut mgO;
mgO = mg(reg_data, mgCtl);
print_coef_table(mgO);

// -------------------------------------------------------------------
// MODEL 2: CCE Mean Group (CCE-MG)
// -------------------------------------------------------------------
print "============================================================";
print "MODEL 2: CCE-MG Estimator";
print "  log_rgdpo ~ log_ck + log_ngd";
print "  CSA: ybar, log_ck_bar, log_ngd_bar (default)";
print "============================================================";

struct mgControl cceCtl;
cceCtl = mgControlCreate();
cceCtl.report = 0;

struct mgOut cceO;
cceO = cce_mg(reg_data, cceCtl);
print_coef_table(cceO);

// -------------------------------------------------------------------
// MODEL 3: DCCE Mean Group (DCCE-MG)
// -------------------------------------------------------------------
print "============================================================";
print "MODEL 3: DCCE-MG Estimator";
print "  log_rgdpo ~ log_ck + log_ngd";
print "  y_lags=1, cr_lags=3, x_csa=log_hc";
print "============================================================";

struct mgControl dcceCtl;
dcceCtl = mgControlCreate();
dcceCtl.y_lags  = 1;
dcceCtl.cr_lags = 3;
dcceCtl.x_csa   = data[., "log_hc"];
dcceCtl.report  = 0;

struct mgOut dcceO;
dcceO = dcce_mg(reg_data, dcceCtl);
print_coef_table(dcceO);

print "============================================================";
print "Validation complete. Results saved to validation/gauss_results.txt";
print "============================================================";

output off;
end;
