/*
** This gauss code computes Common Correlated Effects (CCE) estimates proposed in:  
**         PESARAN, M.H., (2006), ESTIMATION AND INFERENCE IN LARGE HETEROGENEOUS PANELS
**         WITH A MULTIFACTOR ERROR STRUCTURE, ECONOMETRICA, VOL74, NO.4, PP.967-1012.
** Assumes yearly data. If not, year variables should be replaced sequence of numbers 
** Please read the CCEgauss6_Readme.pdf file which accompanied with this gauss code 
*/
new;

#include cce.sdf
#include cce_mg.src
#include dcceutil.src

// Based on original coce from Takashi Yamagata 22.08.2008, revised 12/2022
output file = CCE_test.txt RESET;

// Load test data 
data_file = __FILE_DIR $+ "CCEtestdata.xls";
data = loadd(data_file, ". + date($year, '%Y')");
data = order(data, "id"$|"year");

// Load cross-sectionally common variables
// File Control Options
struct LoadFileControl _ctl;
_ctl = LoadFileControlCreate();
_ctl.xls.sheet = 2;

// Perform import
x_common = loadd(data_file, _ctl);

// Specify the x which would not go 
// for cross section averages
no_xbar = 0;

// Specify the x and id 
// with "zeros" due to normalization
zero_x  = 5;
zero_id = 3;

/*
** Model one: Mean Group Estimation
** In this model no cross-sectional averages are included
** The model has an option to include lag of y (Default)
** or not
*/

struct mgOut mgO;
mgO = mg(data, 0, no_xbar, zero_x, zero_id);

//struct cdOut cdO;
//cdO = cdtest(mgO.e_mg, n, tvec, starttvec-miny+1, endtvec-maxy+maxt);
 
struct mgOut cceOut;
cceOut = cce_mg(data, 0, no_xbar, zero_x, zero_id, 1);

struct mgOut dcceOut;
dcceOut = dcce_mg(data, 0, no_xbar, zero_x, zero_id, 1);
