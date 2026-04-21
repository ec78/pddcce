# dccelib — GAUSS Panel Data Library for Cross-Sectional Dependence

**Dynamic Common-Correlated Effects Estimators for GAUSS**

[![Version](https://img.shields.io/badge/version-1.2.0-blue)](package.json)
[![License](https://img.shields.io/badge/license-GAUSS%20Standard-lightgrey)](https://www.aptech.com)
[![Validated](https://img.shields.io/badge/validated-plm%20%286dp%29-brightgreen)](#validation)

dccelib implements panel data estimators for large heterogeneous panels subject to **cross-sectional dependence** — the setting where unobserved common factors (business cycles, commodity shocks, contagion) generate correlated residuals across units. Standard panel estimators (FE, RE, pooled OLS) are invalid in this setting: standard errors are understated, t-statistics are inflated, and coefficients may themselves be inconsistent.

The library is based on the Common Correlated Effects (CCE) framework of Pesaran (2006) and provides a complete workflow from diagnostic testing through estimation, bias correction, and publication-ready output.

---

## Table of Contents

- [Installation](#installation)
- [Quick Start](#quick-start)
- [Data Format](#data-format)
- [Core Estimators](#core-estimators)
  - [Mean Group (MG)](#mean-group-mg)
  - [CCE Mean Group (CCE-MG)](#cce-mean-group-cce-mg)
  - [Dynamic CCE-MG (DCCE-MG)](#dynamic-cce-mg-dcce-mg)
- [mgControl Options](#mgcontrol-options)
- [Diagnostic Tests](#diagnostic-tests)
  - [CD Test](#pesaran-2004-cd-test)
  - [CIPS Panel Unit Root](#cips-panel-unit-root)
  - [Slope Homogeneity](#slope-homogeneity)
- [Post-Estimation Tools](#post-estimation-tools)
  - [HPJ Bias Correction](#hpj-bias-correction)
  - [Wild Bootstrap SE](#wild-bootstrap-se)
  - [LaTeX Export](#latex-export)
- [Visualization](#visualization)
  - [plotResiduals](#plotresiduals)
  - [plotCoefficients](#plotcoefficients)
  - [plotResidualACF](#plotresidualacf)
- [Advanced Options](#advanced-options)
- [Output Structure](#output-structure)
- [Examples](#examples)
- [Validation](#validation)
- [References](#references)

---

## Installation

Install via the **GAUSS Application Manager** (recommended):

1. Open GAUSS and navigate to **Tools > GAUSS Application Manager**
2. Search for `dccelib` and click **Install**
3. Load the library at the top of your program:

```gauss
library dccelib;
```

> **Note:** Do not install manually from source. The Application Manager handles all dependency resolution and path configuration.

---

## Quick Start

```gauss
new;
library dccelib;

// Load data: columns must be [group | time | y | x1 | x2 ...]
fname = __FILE_DIR $+ "examples/penn_sample.dta";
data  = packr(loadd(fname, ". + date($year, '%Y')"));
data  = order(data, "id"$|"year");

reg_data = data[., "id" "year" "log_rgdpo" "log_ck" "log_ngd"];

// 1. MG baseline (no CD correction)
struct mgOut mgO;
mgO = mg(reg_data);

// 2. CCE-MG (corrects for cross-sectional dependence)
struct mgControl ctl;
ctl = mgControlCreate();
ctl.x_csa = data[., "log_hc"];   // extra variable for CSA proxy

struct mgOut cceO;
cceO = cce_mg(reg_data, ctl);

// 3. DCCE-MG (adds dynamics: lagged y and CSA lags)
ctl.y_lags  = 1;
ctl.cr_lags = 3;

struct mgOut dcceO;
dcceO = dcce_mg(reg_data, ctl);
```

---

## Data Format

All estimator procedures expect a GAUSS dataframe with columns in this order:

| Column | Role |
|--------|------|
| 1 | Panel group ID |
| 2 | Time variable |
| 3 | Dependent variable (y) |
| 4+ | Independent variables (x) |

The panel must be **sorted by group, then time**. Use `packr()` to remove missing rows before calling any estimator.

```gauss
data     = order(data, "id"$|"year");   // sort
reg_data = packr(data[., "id" "year" "log_rgdpo" "log_ck" "log_ngd"]);
```

Unbalanced panels are supported. The library detects the time-series length for each group automatically.

---

## Core Estimators

### Mean Group (MG)

The Pesaran and Smith (1995) MG estimator runs OLS separately for each panel unit and averages the slope estimates. It is consistent under full slope heterogeneity but does **not** correct for cross-sectional dependence. Use as a baseline.

```gauss
struct mgOut mgO;
mgO = mg(reg_data);
```

The CD statistic on the MG residuals (`mgO.cd_stat`) tells you whether cross-sectional dependence is present. A large statistic (relative to the standard normal) indicates that CCE correction is needed.

### CCE Mean Group (CCE-MG)

The Pesaran (2006) CCE-MG estimator augments each unit's OLS regression with cross-sectional averages (CSAs) of y and all x variables. The CSAs act as observable proxies for the unobserved common factors, removing the factor structure from the residuals.

```gauss
struct mgControl ctl;
ctl = mgControlCreate();

// Optional: include extra variables only in the CSA (not as regressors)
ctl.x_csa = data[., "log_hc"];

struct mgOut cceO;
cceO = cce_mg(reg_data, ctl);
```

After CCE-MG, `cceO.cd_stat` should be small (close to zero). If it remains large, increase `cr_lags` or add more variables to `x_csa`.

### Dynamic CCE-MG (DCCE-MG)

The dynamic extension adds lags of y and lags of the cross-sectional averages as regressors. This is appropriate when the dependent variable is serially persistent (GDP, investment, consumption).

```gauss
struct mgControl ctl;
ctl = mgControlCreate();
ctl.y_lags  = 1;   // lags of dependent variable
ctl.cr_lags = 3;   // lags of cross-sectional averages
ctl.x_csa   = data[., "log_hc"];

struct mgOut dcceO;
dcceO = dcce_mg(reg_data, ctl);
```

`cr_lags = 0` activates automatic lag selection via the Andrews (1991) rule: `PT = floor(T^(1/3))`.

---

## mgControl Options

Create the control structure with defaults using `mgControlCreate()`, then set any options you need:

```gauss
struct mgControl ctl;
ctl = mgControlCreate();
```

| Member | Default | Description |
|--------|---------|-------------|
| `y_lags` | `0` | Number of lags of y to include as regressors (DCCE-MG) |
| `cr_lags` | `0` | CSA lag order; `0` = automatic (Andrews 1991 rule) |
| `x_csa` | `0` | Extra matrix of variables to include only in the CSA |
| `pooled` | `0` | `1` = also estimate pooled CCE (Newey-West SE) in the same call |
| `i1` | `0` | `1` = add first-differenced CSAs (KPY 2011; for I(1) data) |
| `two_way` | `0` | `1` = time-demean data before CCE (Bai 2009 two-way structure) |
| `report` | `1` | `1` = print results table; `0` = suppress output |
| `no_xbar` | `0` | Column indices to exclude from CSA computation |
| `x_common` | `0` | Regressors that are common across units |

---

## Diagnostic Tests

### Pesaran (2004) CD Test

The CD test is run automatically by all three estimators. Results are stored in the output struct:

```gauss
// Access CD results after estimation
cceO.cd_stat;   // CD statistic (standard normal under H0: no CD)
cceO.cd_pval;   // two-sided p-value
```

Under H₀ of no cross-sectional dependence, the CD statistic is asymptotically N(0,1). A statistic greater than ~3 (p < 0.01) indicates significant CD.

### CIPS Panel Unit Root

The Pesaran (2007) CIPS test extends the IPS panel unit root test to panels with cross-sectional dependence. Run this before choosing between static and dynamic specifications.

```gauss
// Test whether log_rgdpo has a unit root (p = 1 augmentation lag)
local cips_stat, cadf_vec;
{ cips_stat, cadf_vec } = cips(data[., "id" "year" "log_rgdpo"], 1);
print_cips(cips_stat, cadf_vec, 1, 0);
```

**Arguments:**

| Argument | Description |
|----------|-------------|
| `data` | Dataframe: [group, time, y] |
| `p` | (optional) Number of augmentation lags. Default: automatic |
| `demean` | (optional) `0` = no trend (default), `1` = demeaned, `2` = with trend |

Critical values (Pesaran 2007, Table 2b, N=100, T=50):

| Level | No trend | With trend |
|-------|----------|------------|
| 10% | −2.11 | −2.64 |
| 5% | −2.20 | −2.73 |
| 1% | −2.37 | −2.90 |

Reject H₀ of unit root if CIPS statistic is **below** the critical value (more negative).

### Slope Homogeneity

The Pesaran-Yamagata (2008) test evaluates H₀: all slope coefficients are equal across units. Use this to decide between MG (heterogeneous slopes) and pooled CCE (homogeneous slopes).

```gauss
// Run after cce_mg or dcce_mg
local delta, pval, delta_adj, pval_adj;
{ delta, pval, delta_adj, pval_adj } = slopehomo(cceO);
print_slopehomo(delta, pval, delta_adj, pval_adj);
```

`slopehomo()` takes the estimated `mgOut` struct directly — no need to re-extract matrices manually. It uses the stored per-group X'X matrices (`mgO.xxi_vec`) and residual SDs (`mgO.sig_vec`).

Both the Δ statistic and the bias-adjusted Δ_adj are reported with two-sided p-values. Prefer Δ_adj in samples where N and T are moderate.

---

## Post-Estimation Tools

### HPJ Bias Correction

In dynamic panels with moderate T (T < 40), the MG estimator accumulates an O(1/T) bias. The half-panel jackknife (Dhaene and Jochmans, 2015) corrects this by splitting the time dimension in half and applying the jackknife formula:

b̂_hpj = 2·b̂_full − (b̂_h1 + b̂_h2) / 2

```gauss
struct mgControl ctl;
ctl = mgControlCreate();
ctl.y_lags  = 1;
ctl.cr_lags = 3;
ctl.x_csa   = data[., "log_hc"];

// Apply HPJ to dcce_mg estimates
struct mgOut hpjO;
hpjO = hpj(reg_data, ctl, "dcce_mg");

// Access results
hpjO.b_mg;           // HPJ bias-corrected coefficients
hpjO.b_stats;        // [full, h1, h2, hpj] estimates side by side
```

The third argument selects the estimator: `"mg"`, `"cce_mg"` (default), or `"dcce_mg"`.

Standard errors in `hpjO.se_mg` are the NP standard errors from the full-sample estimate, serving as a conservative bound. For HPJ-specific standard errors, follow with `mgBootstrap()`.

### Wild Bootstrap SE

When sample sizes are small or the asymptotic normal approximation is unreliable, bootstrap standard errors provide inference that is robust to non-normality and heteroskedasticity. dccelib uses the Rademacher wild bootstrap.

```gauss
struct mgControl ctl;
ctl = mgControlCreate();
ctl.y_lags  = 1;
ctl.cr_lags = 3;
ctl.x_csa   = data[., "log_hc"];

// B = 499 bootstrap replications (default = 499)
local se_boot, b_boot;
{ se_boot, b_boot } = mgBootstrap(reg_data, ctl, 499, "dcce_mg");

// b_boot is B×k matrix of bootstrap coefficient draws
// se_boot is k×1 vector of bootstrap standard errors
```

`b_boot` contains the full bootstrap distribution and can be used to construct confidence intervals or perform joint hypothesis tests.

### LaTeX Export

#### Single model

```gauss
// Export CCE-MG results to a .tex file
mgOutToLatex(cceO, "results/cce_results.tex");

// With options: NP standard errors, custom table note
mgOutToLatex(cceO, "results/cce_results.tex", "np", "Penn World Tables, N=93");
```

#### Multiple models side by side

```gauss
// Compare MG, CCE-MG, DCCE-MG in one table
struct mgOut models_arr;
models_arr = mgO | cceO | dcceO;

local labels;
labels = "MG"$|"CCE-MG"$|"DCCE-MG";

mgOutToLatexMulti(models_arr, labels, "results/comparison_table.tex",
                  "Dependent variable: log RGDP. Standard errors in parentheses.");
```

The output is a complete LaTeX `tabular` environment with significance stars (*** p<0.01, ** p<0.05, * p<0.10), CD statistic footer, and mean R² footer. Include in your paper with `\input{results/comparison_table.tex}`.

#### Numeric coefficient matrix

For downstream use (Wald tests, coefficient plots), extract results as a plain matrix:

```gauss
// Returns k×4 matrix: [coef, se, t-stat, p-value]
local ct;
ct = coeftable(cceO);
```

---

## Visualization

Three plot functions provide post-estimation diagnostics. All accept an `mgOut` struct returned by any estimator.

### plotResiduals

Produces a 4-panel residual diagnostic figure:

- **Residuals over observations** — checks for outliers and variance drift
- **Histogram** — checks distributional shape; normality and symmetry improve after CCE
- **Normal Q-Q plot** — heavy-tailed S-curves indicate common factor contamination in plain MG residuals
- **Per-group residual SD (sorted)** — ranked line from smallest to largest SD; reveals whether misfit is concentrated in a few atypical economies or spread evenly

```gauss
plotResiduals(mgO);    // plain MG: look for heavy tails and skewness
plotResiduals(cceO);   // CCE-MG: compare — distribution should tighten
```

### plotCoefficients

Caterpillar plot of per-group slope estimates, sorted ascending, with 95% confidence intervals and a horizontal line at the MG mean. One panel per regressor (up to 6 in a grid). Use this to visualise slope heterogeneity and motivate the MG estimator over pooled alternatives.

```gauss
plotCoefficients(cceO);
```

### plotResidualACF

Bar chart of the sample autocorrelation function of the pooled residuals, from lag 0 to `maxlag` (default 20), with ±1.96/√N significance bands. Significant lag-1 autocorrelation after CCE-MG motivates the dynamic extension (`dcce_mg`).

```gauss
plotResidualACF(cceO);          // default: up to lag 20
plotResidualACF(cceO, 30);      // extend to lag 30
```

---

## Advanced Options

### Pooled CCE in One Call

Run pooled CCE (with Newey-West HAC SE) alongside the MG estimator by setting `pooled = 1`. Results are stored in the embedded `pcce` struct:

```gauss
ctl.pooled = 1;
cceO = cce_mg(reg_data, ctl);

// Access pooled results
cceO.pcce.b_pcce;         // pooled coefficient estimates
cceO.pcce.se_pcce;        // NW standard errors
cceO.pcce.out_pcce;       // formatted output dataframe
```

### I(1) Extension (KPY 2011)

When the regressors and common factors are integrated of order one, standard CCE augmentation in levels is not sufficient. Setting `i1 = 1` adds first differences of the cross-sectional averages (Δȳₜ, Δx̄ₜ) alongside the levels:

```gauss
ctl.i1 = 1;
cceO = cce_mg(reg_data, ctl);   // or dcce_mg
```

Use this when CIPS tests indicate the series are I(1).

### Two-Way Factor Structure (Bai 2009)

When the data contain both unit-specific factor loadings and time-specific aggregate shocks, time-demeaning before CCE augmentation provides additional robustness:

```gauss
ctl.two_way = 1;
cceO = cce_mg(reg_data, ctl);
```

This applies within-time demeaning to y, x, and any `x_csa` variables before computing cross-sectional averages.

---

## Output Structure

All estimators return an `mgOut` struct with the following key fields:

### Coefficient estimates

| Field | Dimensions | Description |
|-------|-----------|-------------|
| `b_mg` | k×1 | Mean group coefficient estimates |
| `se_mg` | k×1 | NP standard errors (Pesaran 2006 eq.58) |
| `tvalue` | k×1 | t-statistics |
| `pval` | k×1 | Two-sided p-values |
| `ci` | k×2 | 95% confidence intervals [lb, ub] |
| `b_vec` | n×k | Individual group estimates |
| `b_stats` | k×4 | Heterogeneity: [min, mean, max, sd] across groups |

### Model diagnostics

| Field | Description |
|-------|-------------|
| `cd_stat` | Pesaran (2004) CD statistic |
| `cd_pval` | CD test p-value |
| `R_sq` | Mean within-group R² (averaged across groups) |

### Model metadata

| Field | Description |
|-------|-------------|
| `panel_var` | Name of the panel group variable |
| `time_var` | Name of the time variable |
| `y_varname` | Name of the dependent variable |
| `mg_vars` | String array of regressor names |
| `model` | Model description string |
| `nobs` | Total observations |
| `ngroups` | Number of panel units (N) |
| `df`, `df_csa` | Degrees of freedom |

### For slopehomo / downstream use

| Field | Description |
|-------|-------------|
| `xxi_vec` | n×k² per-group X'X matrices (row-vectorised) |
| `sig_vec` | n×1 per-group residual standard deviations |
| `e_mg` | Stacked pooled residuals |

---

## Examples

After installation, example scripts are in the **`examples/`** folder of the package directory:

| File | Description |
|------|-------------|
| [mg_penn.e](examples/mg_penn.e) | MG estimator baseline using Penn World Tables |
| [cce_penn.e](examples/cce_penn.e) | CCE-MG estimation with extra CSA variable |
| [dcce_penn.e](examples/dcce_penn.e) | DCCE-MG with lagged y and CSA lags |
| [cce_proc.e](examples/cce_proc.e) | Full workflow: MG → CCE-MG → DCCE-MG |
| [diagnostics.e](examples/diagnostics.e) | CIPS unit root + slope homogeneity testing |
| [advanced_cce.e](examples/advanced_cce.e) | Pooled CCE, I(1) extension, two-way CCE |
| [bias_correction.e](examples/bias_correction.e) | HPJ bias correction + wild bootstrap SE |
| [export_tables.e](examples/export_tables.e) | LaTeX single and multi-model table export |
| [pca_cce.e](examples/pca_cce.e) | PC-CCE-MG with automatic and fixed principal component selection |
| [full_workflow.e](examples/full_workflow.e) | Complete recommended workflow: MG → CIPS → CCE-MG → slope homogeneity → DCCE-MG → HPJ → LaTeX |

All examples use `penn_world.dta` (Penn World Tables, N=93 countries, T≈50 years).

---

## Validation

All three core estimators are validated against R's `plm::pmg()` to **six decimal places** on Penn World Tables data.

| Estimator | Variable | GAUSS | R (plm) |
|-----------|----------|-------|---------|
| MG | log_ck | 0.305300 | 0.305300 |
| MG | log_ngd | 0.279783 | 0.279783 |
| MG | intercept | 5.391778 | 5.391778 |
| CCE-MG | log_ck | 0.316743 | 0.316743 |
| CCE-MG | log_ngd | 0.089055 | 0.089055 |
| CCE-MG | intercept | 1.145539 | 1.145539 |
| DCCE-MG | y_l | 0.456422 | 0.456422 |
| DCCE-MG | log_ck | 0.153173 | 0.153173 |
| DCCE-MG | log_ngd | 0.009159 | 0.009159 |

To reproduce:

```bash
# R validation
Rscript validation/validate_dcce.R

# GAUSS validation (from repo root)
tgauss.exe -b -nj validation/validate_gauss.e
```

---

## References

**Core methodology:**

- Pesaran, M.H. (2006). Estimation and inference in large heterogeneous panels with a multifactor error structure. *Econometrica*, 74(4), 967–1012.
- Pesaran, M.H. and Smith, R. (1995). Estimating long-run relationships from dynamic heterogeneous panels. *Journal of Econometrics*, 68(1), 79–113.
- Pesaran, M.H. (2004). General diagnostic tests for cross-section dependence in panels. *CESifo Working Paper* No. 1229.

**Extensions:**

- Kapetanios, G., Pesaran, M.H. and Yamagata, T. (2011). Panels with non-stationary multifactor error structures. *Journal of Econometrics*, 160(2), 326–348.
- Bai, J. (2009). Panel data models with interactive fixed effects. *Econometrica*, 77(4), 1229–1279.
- Dhaene, G. and Jochmans, K. (2015). Split-panel jackknife estimation of fixed-effect models. *Review of Economic Studies*, 82(3), 991–1030.

**Companion tests:**

- Pesaran, M.H. (2007). A simple panel unit root test in the presence of cross-section dependence. *Journal of Applied Econometrics*, 22(2), 265–312.
- Pesaran, M.H. and Yamagata, T. (2008). Testing slope homogeneity in large panels. *Journal of Econometrics*, 142(1), 50–93.

---

## License

Non-commercial public use only. See the [GAUSS Standard License Agreement](https://www.aptech.com).

## Author

Eric Clower — [eric@aptech.com](mailto:eric@aptech.com)
[Aptech Systems, Inc.](https://www.aptech.com)

For bugs and feature requests, please open an issue on [GitHub](https://github.com/aptech).
