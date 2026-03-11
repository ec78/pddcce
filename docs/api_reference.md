# dccelib API Reference

**Version 0.3.0** — Eric Clower, Aptech Systems

This document provides a complete reference for all public procedures in dccelib. Internal procedures (prefixed with `__`) are not documented here.

---

## Table of Contents

- [Data Structures](#data-structures)
  - [mgControl](#mgcontrol)
  - [mgOut](#mgout)
  - [pcceNWOut](#pccenwout)
- [Core Estimators](#core-estimators)
  - [mg](#mg)
  - [cce_mg](#cce_mg)
  - [dcce_mg](#dcce_mg)
- [Utility Procedures](#utility-procedures)
  - [mgControlCreate](#mgcontrolcreate)
  - [coeftable](#coeftable)
- [Diagnostic Tests](#diagnostic-tests)
  - [cdtest](#cdtest)
  - [cips](#cips)
  - [print_cips](#print_cips)
  - [slopehomo](#slopehomo)
  - [print_slopehomo](#print_slopehomo)
- [Bias Correction](#bias-correction)
  - [hpj](#hpj)
- [Bootstrap Inference](#bootstrap-inference)
  - [mgBootstrap](#mgbootstrap)
- [LaTeX Export](#latex-export)
  - [mgOutToLatex](#mgouttolatex)
  - [mgOutToLatexMulti](#mgouttolatexml)

---

## Data Structures

Structures are defined in `src/cce.sdf` and included automatically when the library is loaded.

### mgControl

Control structure for all estimator procedures. Initialize with `mgControlCreate()`.

```gauss
struct mgControl ctl;
ctl = mgControlCreate();
```

| Member | Type | Default | Description |
|--------|------|---------|-------------|
| `y_lags` | matrix | `0` | Number of lags of y to include as regressors. Used by `dcce_mg`. |
| `cr_lags` | matrix | `0` | CSA lag order. `0` = automatic selection using Andrews (1991) rule: `floor(T^(1/3))`. |
| `x_csa` | matrix | `0` | Additional matrix of variables whose cross-sectional averages should be included as factor proxies but not as direct regressors. Rows must match `data`. |
| `pooled` | matrix | `0` | If `1`, also estimate pooled CCE with Newey-West HAC SE. Results stored in `mgOut.pcce`. |
| `i1` | matrix | `0` | If `1`, apply the Kapetanios-Pesaran-Yamagata (2011) I(1) extension: add first-differenced CSAs (Δȳₜ, Δx̄ₜ) to the augmentation matrix. |
| `two_way` | matrix | `0` | If `1`, time-demean all data before CCE augmentation (Bai 2009). Removes time-specific aggregate shocks. |
| `report` | matrix | `1` | If `1`, print the formatted results table to stdout. Set to `0` to suppress output. |
| `no_xbar` | matrix | `0` | Column indices (in `reg_data`) to exclude from CSA computation. |
| `x_common` | matrix | `0` | Regressors that are common across all units (not unit-specific). |
| `zero_x` | matrix | `0` | Column index of variables to normalize to zero (internal normalization). |
| `zero_id` | matrix | `0` | Group ID for the normalization reference group. |

---

### mgOut

Output structure returned by `mg()`, `cce_mg()`, `dcce_mg()`, and `hpj()`.

#### Coefficient estimates

| Member | Dimensions | Description |
|--------|-----------|-------------|
| `b_mg` | k×1 | Mean Group coefficient estimates. For HPJ output: bias-corrected estimates. |
| `se_mg` | k×1 | Standard errors (NP estimator, Pesaran 2006 eq.58). |
| `tvalue` | k×1 | t-statistics: `b_mg ./ se_mg`. |
| `pval` | k×1 | Two-sided p-values from t-distribution. |
| `ci` | k×2 | 95% confidence intervals: `[b_mg - 1.96*se_mg, b_mg + 1.96*se_mg]`. |
| `cov_mg` | k×k | NP covariance matrix of the MG estimator. |
| `b_vec` | n×k | Individual group OLS estimates. Row i = group i's slope vector. |
| `b_stats` | k×4 | Heterogeneity of slopes: `[min, mean, max, sd]` across groups. For HPJ output: `[full, h1, h2, hpj]`. |
| `out_mg` | dataframe | Formatted results dataframe (variable names + [coef, se, t, p]). |

#### Diagnostics

| Member | Description |
|--------|-------------|
| `cd_stat` | Pesaran (2004) CD statistic. Asymptotically N(0,1) under H₀ of no CD. |
| `cd_pval` | Two-sided p-value for the CD test. |
| `R_sq` | Mean within-group R², averaged across all panel units. |

#### For downstream procedures

| Member | Dimensions | Description |
|--------|-----------|-------------|
| `e_mg` | nobs×1 | Stacked residuals from all group OLS regressions. |
| `xxi_vec` | n×k² | Per-group X'X matrices, row-vectorised. Used by `slopehomo()`. |
| `sig_vec` | n×1 | Per-group residual standard deviations. Used by `slopehomo()`. |
| `se` | n×k | Per-group OLS standard errors. |
| `xxi` | k×k | X'X for the last group (legacy field). |

#### Model metadata

| Member | Description |
|--------|-------------|
| `panel_var` | Name of the panel group variable (column 1 of `data`). |
| `time_var` | Name of the time variable (column 2 of `data`). |
| `y_varname` | Name of the dependent variable. |
| `mg_vars` | String array of regressor names (including CSA variables). |
| `csa_vars` | String array of CSA variable names. |
| `model` | Model description string (e.g., `"CCE Mean Group [CCE-MG]"`). |
| `nobs` | Total number of observations. |
| `ngroups` | Number of panel units (N). |
| `numvars` | Number of original (non-CSA) regressors. |
| `k_reg` | Total regressors in each group's OLS (includes CSA variables). |
| `obs_grp` | n×1 vector of observations per group. |
| `df` | Degrees of freedom for MG inference. |
| `df_csa` | Degrees of freedom accounting for CSA augmentation. |
| `csa_lags` | Number of CSA lags used. |

#### Pooled CCE (when `mgCtl.pooled = 1`)

| Member | Description |
|--------|-------------|
| `pcce` | Embedded `pcceNWOut` struct with pooled CCE results. |

---

### pcceNWOut

Embedded within `mgOut.pcce` when `mgCtl.pooled = 1`.

| Member | Description |
|--------|-------------|
| `b_pcce` | k×1 pooled CCE coefficient estimates. |
| `se_pcce` | k×1 Newey-West HAC standard errors. |
| `tvalue_pcce` | k×1 t-statistics. |
| `out_pcce` | Formatted output dataframe. |
| `cov_pcce_rbst_nw` | k×k NW HAC covariance matrix. |
| `cov_pcce_hs2` | k×k Pesaran hs2 covariance matrix. |
| `e_pcce` | Pooled CCE residuals. |
| `sig2_pcce` | Pooled residual variance. |
| `sig_i_vec` | Per-group residual variances. |

---

## Core Estimators

### mg

**Mean Group estimator** (Pesaran and Smith 1995).

```gauss
proc (1) = mg(data [, mgCtl]);
```

**Arguments:**

| Parameter | Required | Description |
|-----------|----------|-------------|
| `data` | Yes | Dataframe: [group, time, y, x₁, ..., xₖ]. Sorted by group then time. No missing values. |
| `mgCtl` | No | `mgControl` struct. If omitted, uses all defaults. |

**Returns:** `mgOut` struct.

**Notes:**
- MG does not augment with CSAs — it does not correct for cross-sectional dependence.
- Use the returned `cd_stat` to test whether CCE correction is needed.
- If `mgCtl.pooled = 1`, pooled OLS (without CSA augmentation) is estimated and stored in `mgOut.pcce`.

**Example:**
```gauss
struct mgOut mgO;
mgO = mg(reg_data);
```

---

### cce_mg

**CCE Mean Group estimator** (Pesaran 2006).

```gauss
proc (1) = cce_mg(data [, mgCtl]);
```

**Arguments:** Same as `mg()`.

**Returns:** `mgOut` struct.

**Notes:**
- Augments each unit's regression with cross-sectional averages of all y and x variables (and `mgCtl.x_csa` if provided).
- `mgCtl.cr_lags` controls the number of CSA lags; `0` uses automatic selection.
- If `mgCtl.i1 = 1`, adds first-differenced CSAs (Δȳ, Δx̄) for I(1) robustness.
- If `mgCtl.two_way = 1`, time-demeans all variables before computing CSAs.
- If `mgCtl.pooled = 1`, also runs pooled CCE with NW SE; results in `mgOut.pcce`.

**Example:**
```gauss
struct mgControl ctl;
ctl = mgControlCreate();
ctl.x_csa = data[., "log_hc"];

struct mgOut cceO;
cceO = cce_mg(reg_data, ctl);
```

---

### dcce_mg

**Dynamic CCE Mean Group estimator** (Chudik and Pesaran 2015).

```gauss
proc (1) = dcce_mg(data [, mgCtl]);
```

**Arguments:** Same as `mg()`.

**Returns:** `mgOut` struct.

**Notes:**
- Extends `cce_mg()` with `mgCtl.y_lags` lags of the dependent variable and `mgCtl.cr_lags` lags of the cross-sectional averages.
- Both `i1` and `two_way` flags are supported.
- For small-T samples, follow with `hpj()` to correct O(1/T) bias.

**Example:**
```gauss
struct mgControl ctl;
ctl = mgControlCreate();
ctl.y_lags  = 1;
ctl.cr_lags = 3;
ctl.x_csa   = data[., "log_hc"];

struct mgOut dcceO;
dcceO = dcce_mg(reg_data, ctl);
```

---

## Utility Procedures

### mgControlCreate

Returns an `mgControl` struct with all members set to their default values.

```gauss
proc (1) = mgControlCreate();
```

**Returns:** `mgControl` struct with defaults (see [mgControl table](#mgcontrol)).

**Example:**
```gauss
struct mgControl ctl;
ctl = mgControlCreate();
ctl.y_lags = 1;    // override specific defaults
```

---

### coeftable

Extracts the coefficient table from an `mgOut` struct as a plain numeric matrix.

```gauss
proc (1) = coeftable(mgO);
```

**Arguments:**

| Parameter | Description |
|-----------|-------------|
| `mgO` | `mgOut` struct from any estimator. |

**Returns:** k×4 numeric matrix with columns `[coef, se, t-stat, p-value]`.

**Example:**
```gauss
local ct;
ct = coeftable(cceO);

// Access the capital coefficient row
local coef_ck, se_ck, t_ck, p_ck;
coef_ck = ct[1, 1];
se_ck   = ct[1, 2];
t_ck    = ct[1, 3];
p_ck    = ct[1, 4];
```

---

## Diagnostic Tests

### cdtest

**Pesaran (2004) CD test** for cross-sectional dependence.

```gauss
proc (1) = cdtest(e, n, tivec, starttvec, endtvec);
```

This procedure is called internally by all estimators. Access results via `mgOut.cd_stat` and `mgOut.cd_pval` rather than calling directly.

**Returns:** `cdOut` struct with fields: `cd_nt`, `pvalue_cd_nt`, `meanrho`, `used`.

---

### cips

**Pesaran (2007) CIPS panel unit root test.**

```gauss
proc (2) = cips(data [, p, demean]);
```

**Arguments:**

| Parameter | Required | Default | Description |
|-----------|----------|---------|-------------|
| `data` | Yes | — | Dataframe with columns: [group, time, y]. |
| `p` | No | Automatic | Number of augmentation lags. Automatic uses `floor(T^(1/4))`. |
| `demean` | No | `0` | `0` = no trend, `1` = demeaned, `2` = with linear trend. |

**Returns:**

| Return value | Description |
|-------------|-------------|
| `cips_stat` | CIPS statistic (mean of group-level CADF t-ratios). |
| `cadf_vec` | n×1 vector of individual group CADF t-ratios. |

**Critical values** (Pesaran 2007, Table 2b, no trend):

| N\T | 20 | 30 | 50 | 100 |
|-----|----|----|----|----|
| 20 (5%) | −2.27 | −2.22 | −2.20 | −2.18 |
| 50 (5%) | −2.25 | −2.22 | −2.20 | −2.18 |
| 100 (5%) | −2.24 | −2.21 | −2.20 | −2.18 |

Reject H₀ of unit root if `cips_stat` is more negative than the critical value.

**Example:**
```gauss
local cips_stat, cadf_vec;
{ cips_stat, cadf_vec } = cips(data[., "id" "year" "log_rgdpo"], 1, 0);
print_cips(cips_stat, cadf_vec, 1, 0);
```

---

### print_cips

Prints CIPS test results in a formatted table.

```gauss
proc (0) = print_cips(cips_stat, cadf_vec, p, demean);
```

**Arguments:** Results from `cips()` and the same `p` and `demean` values passed to `cips()`.

---

### slopehomo

**Pesaran-Yamagata (2008) slope homogeneity test.**

```gauss
proc (4) = slopehomo(mgO);
```

**Arguments:**

| Parameter | Description |
|-----------|-------------|
| `mgO` | `mgOut` struct. Must contain `xxi_vec`, `sig_vec`, and `b_vec` (populated by all core estimators). |

**Returns:**

| Return value | Description |
|-------------|-------------|
| `delta` | Δ̂ statistic. Asymptotically N(0,1) under H₀. |
| `pval` | Two-sided p-value for Δ̂. |
| `delta_adj` | Bias-adjusted Δ̂_adj statistic (preferred in moderate samples). |
| `pval_adj` | Two-sided p-value for Δ̂_adj. |

**Notes:**
- H₀: all slope coefficients are equal across units (β_i = β for all i).
- H₁: slope heterogeneity.
- Uses the Swamy S statistic: S = Σᵢ (b̂ᵢ − b̄)' (σᵢ⁻² Xᵢ'Xᵢ) (b̂ᵢ − b̄).
- Reject H₀ at 5% if |Δ̂_adj| > 1.96.

**Example:**
```gauss
local delta, pval, delta_adj, pval_adj;
{ delta, pval, delta_adj, pval_adj } = slopehomo(cceO);
print_slopehomo(delta, pval, delta_adj, pval_adj);
```

---

### print_slopehomo

Prints slope homogeneity test results in a formatted table.

```gauss
proc (0) = print_slopehomo(delta, pval, delta_adj, pval_adj);
```

**Arguments:** Results from `slopehomo()`.

---

## Bias Correction

### hpj

**Half-Panel Jackknife bias correction** (Dhaene and Jochmans 2015).

```gauss
proc (1) = hpj(data, mgCtl [, estimator_type]);
```

**Arguments:**

| Parameter | Required | Default | Description |
|-----------|----------|---------|-------------|
| `data` | Yes | — | Panel dataframe in the same format as the core estimators. |
| `mgCtl` | Yes | — | `mgControl` struct with the same settings as the target estimator. |
| `estimator_type` | No | `"cce_mg"` | String: `"mg"`, `"cce_mg"`, or `"dcce_mg"`. |

**Returns:** `mgOut` struct with:
- `b_mg`: HPJ bias-corrected coefficients.
- `se_mg`: NP standard errors from the full-sample estimate (conservative).
- `tvalue`, `pval`, `ci`: Recomputed using HPJ `b_mg` and full-sample `se_mg`.
- `b_stats`: k×4 matrix with columns `[b_full, b_h1, b_h2, b_hpj]`.
- `model`: Full-sample model string with `"[HPJ Bias-Corrected]"` appended.

**Notes:**
- Corrects the O(1/T) bias in MG estimators applied to dynamic panels.
- Recommended for `dcce_mg` when T < 40.
- The panel is split at `floor(T_i/2)` for each group i independently, supporting unbalanced panels.
- If `mgCtl.x_csa` is set, it is automatically split with the same indices.
- For HPJ-specific standard errors (rather than the conservative full-sample NP SE), follow with `mgBootstrap()`.

**Example:**
```gauss
struct mgControl ctl;
ctl = mgControlCreate();
ctl.y_lags  = 1;
ctl.cr_lags = 3;
ctl.x_csa   = data[., "log_hc"];

struct mgOut hpjO;
hpjO = hpj(reg_data, ctl, "dcce_mg");

print "HPJ coefficients:";
print hpjO.b_mg;

print "Full vs. HPJ comparison:";
print hpjO.b_stats;    // [full, h1, h2, hpj] side by side
```

---

## Bootstrap Inference

### mgBootstrap

**Wild Rademacher bootstrap standard errors** for MG estimators.

```gauss
proc (2) = mgBootstrap(data, mgCtl [, B, estimator_type]);
```

**Arguments:**

| Parameter | Required | Default | Description |
|-----------|----------|---------|-------------|
| `data` | Yes | — | Panel dataframe in the same format as the core estimators. |
| `mgCtl` | Yes | — | `mgControl` struct. |
| `B` | No | `499` | Number of bootstrap replications. Use B ≥ 999 for publication. |
| `estimator_type` | No | `"cce_mg"` | String: `"mg"`, `"cce_mg"`, or `"dcce_mg"`. |

**Returns:**

| Return value | Description |
|-------------|-------------|
| `se_boot` | k×1 vector of bootstrap standard errors (standard deviation of `b_boot` columns). |
| `b_boot` | B×k matrix of bootstrap coefficient draws. |

**Algorithm:**
1. Estimate the target model on the original data; obtain residuals ê_i per group.
2. For each bootstrap replication b = 1, ..., B:
   - Draw Rademacher weights w_i ~ {−1, +1} with equal probability for each group.
   - Construct y*_{it} = ŷ_{it} + w_i · ê_{it}.
   - Re-estimate the model on the bootstrap panel.
3. Compute `se_boot = std(b_boot)` column-wise.

**Notes:**
- Use `b_boot` for percentile confidence intervals: `quantile(b_boot[., j], 0.025)` and `quantile(b_boot[., j], 0.975)`.
- Computationally expensive for large B or `dcce_mg`; consider B = 99–199 for exploratory work.
- For Wald tests, construct the bootstrap distribution of the test statistic rather than using asymptotic critical values.

**Example:**
```gauss
struct mgControl ctl;
ctl = mgControlCreate();
ctl.y_lags  = 1;
ctl.cr_lags = 3;
ctl.x_csa   = data[., "log_hc"];

local se_boot, b_boot;
{ se_boot, b_boot } = mgBootstrap(reg_data, ctl, 499, "dcce_mg");

// Percentile CI for each coefficient
local j;
for j(1, cols(b_boot), 1);
    print "CI for coeff " $+ ntos(j, 1) $+ ": [" $+
          ntos(quantile(b_boot[., j], 0.025), 4) $+ ", " $+
          ntos(quantile(b_boot[., j], 0.975), 4) $+ "]";
endfor;
```

---

## LaTeX Export

### mgOutToLatex

**Exports a single model's results to a LaTeX table file.**

```gauss
proc (0) = mgOutToLatex(mgO, filename [, se_type, note]);
```

**Arguments:**

| Parameter | Required | Default | Description |
|-----------|----------|---------|-------------|
| `mgO` | Yes | — | `mgOut` struct from any estimator. |
| `filename` | Yes | — | Output file path (absolute or relative to CWD). Extension `.tex` recommended. |
| `se_type` | No | `"np"` | Standard error type: `"np"` (NP estimator) or `"nw"` (Newey-West, if available). |
| `note` | No | `""` | String appended as a table note at the bottom. |

**Output format:**
- LaTeX `tabular` environment with `\toprule`, `\midrule`, `\bottomrule` (requires `booktabs` package).
- Rows: variable name | coefficient | standard error (parenthesised).
- Significance stars: `***` p<0.01, `**` p<0.05, `*` p<0.10.
- Footer rows: CD statistic with p-value; mean R².
- Include with: `\input{filename.tex}`.

**Example:**
```gauss
mgOutToLatex(cceO, "tables/cce_results.tex", "np",
             "Penn World Tables. NP standard errors in parentheses.");
```

---

### mgOutToLatexMulti

**Exports multiple models side by side in one LaTeX table.**

```gauss
proc (0) = mgOutToLatexMulti(mgO_arr, labels, filename [, note]);
```

**Arguments:**

| Parameter | Required | Description |
|-----------|----------|-------------|
| `mgO_arr` | Yes | Array of `mgOut` structs, created with the `|` operator: `mgO|cceO|dcceO`. Supports 2–6 models. |
| `labels` | Yes | String array of column headers, one per model: `"MG"$|"CCE-MG"$|"DCCE-MG"`. |
| `filename` | Yes | Output file path. |
| `note` | No | Table note string. |

**Output format:**
- One column per model; variable names in the leftmost column.
- Coefficients and parenthesised SEs in each model column.
- Significance stars applied per coefficient per model.
- Footer: observations, CD statistic, R² for each model.

**Example:**
```gauss
mgOutToLatexMulti(mgO|cceO|dcceO,
                  "MG"$|"CCE-MG"$|"DCCE-MG",
                  "tables/comparison.tex",
                  "Dependent variable: $\\log(RGDP_o)$.");
```

---

## Recommended Workflow

```
1. Load and prepare data
   packr() → order() → select columns [group, time, y, x...]

2. Test for cross-sectional dependence
   mg() → inspect cd_stat on MG residuals

3. Test for unit roots
   cips() per variable → print_cips()

4. Estimate CCE-MG (or DCCE-MG if series are persistent)
   cce_mg() / dcce_mg()
   → use i1=1 if series are I(1)
   → use two_way=1 if two-way factor structure suspected

5. Test slope homogeneity
   slopehomo() → print_slopehomo()
   → if H0 rejected: MG estimates are appropriate
   → if H0 not rejected: pooled CCE may be efficient (pooled=1)

6. Check residual CD
   cceO.cd_stat → should be small after CCE augmentation

7. Bias correction (if dynamic model with moderate T)
   hpj() → HPJ-corrected b_mg

8. Bootstrap SE (robustness check)
   mgBootstrap()

9. Export results
   mgOutToLatex() / mgOutToLatexMulti()
```
