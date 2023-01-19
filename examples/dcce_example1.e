new;

#include cce.sdf
#include cce_mg.src
#include dcceutil.src

/*
** The dataset contains the following:
**    log_rgdpo           Real GDP
**    log_hc              Human Capital
**    log_ck              Physical Capital
**    log_ngd             Population growth + break even investments of 5%
*/

// First load data
fname = __FILE_DIR $+ "penn_sample.dta";
data = loadd(fname, ". + date($year, '%Y')");


// Control structure
struct mgControl mgCtl;
mgCtl = mgControlCreate();

mgCtl.y_lags = 1;
mgCtl.cr_lags = 3;

// Set up reg data
reg_data = data[., "id" "year" "log_rgdpo" "log_ck" "log_ngd"];

// Call regression
struct mgOut mgO;
mgO = dcce_mg(packr(reg_data), mgctl);
