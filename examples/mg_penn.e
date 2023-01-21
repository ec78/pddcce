new;
library dccelib;

/*
** The dataset contains the following:
**    log_rgdpo           Real GDP
**    log_hc              Human Capital
**    log_ck              Physical Capital
**    log_ngd             Population growth + break even investments of 5%
*/

// First load data
fname = "examples/penn_sample.dta";
data = loadd(fname, ". + date($year, '%Y')");


// Set up reg data
reg_data = data[., "id" "year" "log_rgdpo" "log_ck" "log_ngd"];

// Call regression
struct mgOut mgO;
mgO = mg(packr(reg_data));
