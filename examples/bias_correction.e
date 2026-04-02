/*
** bias_correction.e
**
** Demonstrates bias correction and robust inference for dynamic CCE models:
**
**   1. Half-Panel Jackknife (HPJ) bias correction — hpj()
**      Corrects the O(1/T) bias that accumulates in MG estimators applied
**      to dynamic panels with moderate T. Recommended for DCCE-MG when T < 40.
**
**   2. Wild Bootstrap Standard Errors — mgBootstrap()
**      Rademacher wild bootstrap for inference robust to non-normality
**      and heteroskedasticity. Use as a robustness check alongside NP SEs,
**      or after HPJ correction to get HPJ-specific confidence intervals.
**
** References:
**   Dhaene, G. and Jochmans, K. (2015). Split-panel jackknife estimation of
**     fixed-effect models. Review of Economic Studies, 82(3), 991-1030.
**
** Data: Penn World Tables (N=93 countries, T~50 years)
*/

new;
library dccelib;

// -----------------------------------------------------------------------
// Load data
// -----------------------------------------------------------------------
fname = __FILE_DIR $+ "penn_world.dta";
data  = packr(loadd(fname, ". + date($year, '%Y')"));
data  = order(data, "id"$|"year");

reg_data = data[., "id" "year" "log_rgdpo" "log_ck" "log_ngd"];


// -----------------------------------------------------------------------
// Part 1: HPJ Bias Correction
// -----------------------------------------------------------------------
// The HPJ formula: b_hpj = 2*b_full - 0.5*(b_h1 + b_h2)
// where b_full is estimated on the full panel and b_h1, b_h2 on each
// temporal half. This removes the leading O(1/T) bias term.
//
// hpj() accepts the same data and mgControl struct you would pass to the
// underlying estimator. The third argument selects the estimator:
//   "mg"      -- plain Mean Group
//   "cce_mg"  -- CCE Mean Group (default)
//   "dcce_mg" -- Dynamic CCE-MG
// -----------------------------------------------------------------------
print "=================================================================";
print "PART 1: Half-Panel Jackknife (HPJ) Bias Correction";
print "=================================================================";

ctl = mgControlCreate();
ctl.y_lags  = 1;
ctl.cr_lags = 3;
ctl.x_csa   = data[., "log_hc"];

// Step 1: full-sample DCCE-MG (uncorrected baseline)
print "--- Full-sample DCCE-MG (uncorrected) ---";
dcceO = dcce_mg(reg_data, ctl);

// Step 2: HPJ-corrected DCCE-MG
print;
print "--- HPJ bias-corrected DCCE-MG ---";
hpjO = hpj(reg_data, ctl, "dcce_mg");

// Compare full-sample vs. HPJ estimates
print;
print "Coefficient comparison (first k=numvars coefficients):";
k = hpjO.numvars;

print "               Uncorrected       HPJ-corrected";
print "               -----------       -------------";
for i(1, k, 1);
    print hpjO.mg_vars[i] $+
          "     " $+ ntos(dcceO.b_mg[i], 6) $+
          "     " $+ ntos(hpjO.b_mg[i], 6);
endfor;

print;
print "b_stats_hpj columns: [full-panel, first-half, second-half, HPJ]";
print hpjO.b_stats_hpj;

// -----------------------------------------------------------------------
// HPJ also works for CCE-MG
// -----------------------------------------------------------------------
print;
print "--- HPJ bias-corrected CCE-MG ---";
cce_ctl = mgControlCreate();
cce_ctl.x_csa = data[., "log_hc"];

hpjO_cce = hpj(reg_data, cce_ctl, "cce_mg");


// -----------------------------------------------------------------------
// Part 2: Wild Bootstrap Standard Errors
// -----------------------------------------------------------------------
// The Rademacher wild bootstrap resamples residuals using random weights
// w ~ {-1, +1} with equal probability, preserving the heteroskedasticity
// structure of the original errors.
//
// mgBootstrap() returns:
//   se_boot  -- k×1 bootstrap standard errors
//   b_boot   -- B×k matrix of bootstrap coefficient draws
//
// The full b_boot matrix can be used for:
//   - Non-parametric confidence intervals (e.g., percentile method)
//   - Joint hypothesis tests via bootstrap Wald statistics
//   - Visualizing the sampling distribution of each coefficient
// -----------------------------------------------------------------------
print;
print "=================================================================";
print "PART 2: Wild Bootstrap Standard Errors";
print "=================================================================";

// Bootstrap SE for DCCE-MG (B=199 for speed; use B=999 in practice)
B = 199;

print "Running " $+ ntos(B, 1) $+ " bootstrap replications for DCCE-MG...";
{ se_boot, b_boot } = mgBootstrap(reg_data, ctl, B, "dcce_mg");

print;
print "Bootstrap SE vs. NP SE comparison:";
print "               NP SE        Bootstrap SE";
print "               -----        ------------";
for i(1, k, 1);
    print hpjO.mg_vars[i] $+
          "     " $+ ntos(dcceO.se_mg[i], 6) $+
          "     " $+ ntos(se_boot[i], 6);
endfor;

print;
print "b_boot is a " $+ ntos(B, 1) $+ "x" $+ ntos(k, 1) $+
      " matrix of bootstrap draws.";
print "Use percentile CI: [quantile(b_boot[.,j], 0.025), quantile(b_boot[.,j], 0.975)]";

// 95% percentile confidence intervals from bootstrap
print;
print "95% Bootstrap percentile confidence intervals:";
for i(1, k, 1);
    ci_lower = quantile(b_boot[., i], 0.025);
    ci_upper = quantile(b_boot[., i], 0.975);
    print hpjO.mg_vars[i] $+ "  [" $+
          ntos(ci_lower, 6) $+ ", " $+ ntos(ci_upper, 6) $+ "]";
endfor;

// -----------------------------------------------------------------------
// Part 3: Combining HPJ + Bootstrap
// -----------------------------------------------------------------------
// After HPJ bias correction, the NP SE from the full-sample estimate is
// a conservative approximation. For HPJ-specific SEs, run mgBootstrap()
// with "dcce_mg" (hpj wraps the estimator internally).
//
// Note: bootstrapping after HPJ is computationally expensive because
// each bootstrap iteration internally calls dcce_mg three times.
// Consider reducing B (e.g., B=99) for exploratory analysis.
// -----------------------------------------------------------------------
print;
print "=================================================================";
print "PART 3: Bootstrap SE applied directly to CCE-MG (faster)";
print "=================================================================";

cce_ctl2 = mgControlCreate();
cce_ctl2.x_csa = data[., "log_hc"];

{ se_cce_boot, b_cce_boot } = mgBootstrap(reg_data, cce_ctl2, B, "cce_mg");

print;
print "CCE-MG NP SE vs. Bootstrap SE:";
k_cce = hpjO_cce.numvars;
for i(1, k_cce, 1);
    print hpjO_cce.mg_vars[i] $+
          "     NP: " $+ ntos(hpjO_cce.se_mg[i], 6) $+
          "     Boot: " $+ ntos(se_cce_boot[i], 6);
endfor;
