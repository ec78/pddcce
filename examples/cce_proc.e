/*
** This GAUSS code demonstrates MG, CCE-MG, and DCCE-MG estimation using
** the Penn World Tables sample dataset.
**
** Reference:
**   Pesaran, M.H. (2006), "Estimation and Inference in Large Heterogeneous
**   Panels with a Multifactor Error Structure," Econometrica, 74(4), 967-1012.
**
** Originally based on code by Takashi Yamagata (2008), updated to use
** the mgControl structure API.
*/

new;
library dccelib;

/*
** The dataset contains the following variables:
**    log_rgdpo    Real GDP per output (log)
**    log_hc       Human Capital index (log)
**    log_ck       Physical Capital stock (log)
**    log_ngd      Population growth + break-even investment 5% (log)
*/

// Load data
fname = __FILE_DIR $+ "penn_world.dta";
data = packr(loadd(fname, ". + date($year, '%Y')"));
data = order(data, "id"$|"year");

// Regression data: group | time | y | x1 | x2
reg_data = data[., "id" "year" "log_rgdpo" "log_ck" "log_ngd"];

/*
** Model 1: Mean Group (MG) Estimator
** No cross-sectional averages included.
*/
mgO = mg(reg_data);

/*
** Model 2: Common Correlated Effects MG (CCE-MG) Estimator
** Cross-sectional averages of all regressors augment each group regression.
** log_hc is included only in the CSA, not as a direct regressor.
*/
cceCtl = mgControlCreate();
cceCtl.x_csa = data[., "log_hc"];

cceO = cce_mg(reg_data, cceCtl);

/*
** Model 3: Dynamic CCE-MG (DCCE-MG) Estimator
** Adds lagged dependent variable and CSA lags to the CCE-MG specification.
*/
dcceCtl = mgControlCreate();
dcceCtl.y_lags  = 1;
dcceCtl.cr_lags = 3;
dcceCtl.x_csa   = data[., "log_hc"];

dcceO = dcce_mg(reg_data, dcceCtl);
