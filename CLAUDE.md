# CLAUDE.md — GAUSS dccelib Project Context

## Project Overview

**dccelib** is a GAUSS econometrics library implementing panel data estimators based on:

> Pesaran, M.H. (2006), "Estimation and Inference in Large Heterogeneous Panels with a Multifactor Error Structure," *Econometrica*, Vol. 74, No. 4, pp. 967–1012.

The library provides three main estimators plus companion diagnostic and inference tools.

**Author:** Eric Clower, Aptech Systems, Inc. (`eric@aptech.com`)
**Target:** GAUSS 23+
**License:** Non-commercial public use only
**Package version:** 1.1.0

---

## Repository Structure

```
pddcce/
├── CLAUDE.md               # This file
├── README.md               # User-facing overview
├── package.json            # GAUSS package manifest (name: dccelib, version: 1.1.0)
├── src/
│   ├── cce.sdf             # Structure definitions (mgControl, mgOut, cdOut, pcceNWOut, cceRankOut)
│   ├── cce_mg.src          # MG, CCE-MG, DCCE-MG estimators + pcceNW + cdtest
│   ├── dcceutil.src        # Utility procedures, printing, coeftable(), __selectVars()
│   ├── slopehomo.src       # Pesaran-Yamagata (2008) slope homogeneity tests
│   ├── cips.src            # Pesaran (2007) CIPS panel unit root test
│   ├── latex_export.src    # LaTeX table export (single and multi-model)
│   ├── bootstrap.src       # Wild bootstrap standard errors for MG estimators
│   ├── bias_correct.src    # Half-panel jackknife (HPJ) bias correction
│   └── diagnostics.src     # plotResiduals(), cce_rank(), print_cce_rank()
├── docs/
│   ├── api_reference.md    # Full public API reference
│   ├── blog_outline.md     # Medium blog post outline
│   └── research_proposal.md # Academic paper/proposal
├── examples/
│   ├── mg_penn.e           # MG estimator example (Penn World Tables)
│   ├── cce_penn.e          # CCE-MG example
│   ├── dcce_penn.e         # DCCE-MG example
│   ├── cce_proc.e          # Combined MG/CCE-MG/DCCE-MG example
│   ├── diagnostics.e       # CIPS + slope homogeneity workflow
│   ├── advanced_cce.e      # pooled CCE, I(1), two-way CCE options
│   ├── bias_correction.e   # HPJ bias correction + wild bootstrap SE
│   ├── export_tables.e     # LaTeX single and multi-model export
│   ├── penn_sample.dta     # Penn World Tables sample data (Stata format)
│   └── jasa2.dta           # Additional dataset
└── validation/
    ├── validate_dcce.R     # R plm::pmg() reference values
    ├── validate_gauss.e    # GAUSS batch validation script
    ├── inspect_data.R      # Data inspection utility
    └── gauss_results.txt   # Output from validate_gauss.e (gitignored)
```

---

## GAUSS Language Reference

For writing and reviewing GAUSS code, use the LLM reference at:
**https://github.com/aptech/gauss-llm-reference**

Key GAUSS conventions used in this codebase:
- Procedures: `proc(N) = name(args);` ... `retp(...);` ... `endp;`
- Structs declared in `.sdf` files, instantiated with `struct TypeName varname;`
- Optional args via `dynargsGet(index, default)` or `dynargsGet(1|3, d1, d2, d3)`
- `do while` / `endo` loops; `for i(start, end, step)` / `endfor`
- String arrays: `$~` horizontal concat, `$|` vertical concat
- Element-wise ops: `.*`, `./`, `.^`, `.==`, `.!=`
- `pinv()` pseudoinverse; `invpd()` positive-definite inverse
- `packr()` removes rows with missing values; `miss(x, 0)` replaces 0 with missing
- `asDF()` converts matrix to dataframe; `setcolnames()` / `getcolnames()`
- `aggregate()` for group aggregation; `unique()` / `counts()` / `indnv()`
- `#include` resolves relative to CWD, not to the including file's directory
- Forward slashes required in paths (backslash sequences `\e`, `\p`, etc. are escape codes)
- `sprintf(fmt, ...)` as a bare statement prints to stdout in GAUSS 26

---

## Data Structures (`src/cce.sdf`)

### `mgControl` — Input control structure
| Member | Default | Description |
|--------|---------|-------------|
| `zero_x` | 0 | Column ID for vars with zeros (normalization) |
| `zero_id` | 0 | Group ID for normalized groups |
| `y_lags` | 0 | 1 = include lagged y as regressor |
| `cr_lags` | 0 | CSA lag order (0 = auto via `__getPT`) |
| `no_xbar` | 0 | Columns to exclude from CSA |
| `x_common` | 0 | Common (non-group-specific) regressors |
| `x_csa` | 0 | Extra vars for CSA not in main regression |
| `report` | 1 | 1 = print results table automatically |
| `pooled` | 0 | 1 = also run pooled CCE (`pcceNW`) |
| `i1` | 0 | 1 = add first-differenced CSA (KPY 2011; for I(1) regressors) |
| `two_way` | 0 | 1 = time-demean data before CCE (two-way CCE with time FE) |
| `y_var` | `""` | Name of dependent variable column. If set, columns are automatically reordered. |
| `x_vars` | `""` | Name(s) of regressor columns (`"x1" $\| "x2"`). Used with `y_var`. |

### `mgOut` — Results output structure
| Member | Description |
|--------|-------------|
| `b_mg` | k×1 MG coefficient estimates |
| `se_mg` | k×1 NP standard errors (Pesaran 2006 eq.58) |
| `tvalue` | k×1 t-values |
| `pval` | k×1 p-values |
| `ci` | k×2 95% CI [lb, ub] |
| `R_sq` | Mean within-group R² |
| `b_stats` | k×4 heterogeneity: [min, mean, max, sd] of b_i across groups |
| `cov_mg` | k×k NP covariance matrix |
| `e_mg` | Stacked pooled residuals |
| `b_vec` | n×k individual group estimates |
| `xxi_vec` | n×k² per-group X'X (row-vectorised); used by `slopehomo` |
| `sig_vec` | n×1 per-group residual SDs |
| `se_nw` | n×k per-group Newey-West SEs |
| `cd_stat`, `cd_pval` | Pesaran (2004) CD test |
| `out_mg` | Formatted output dataframe |
| `pcce` | `pcceNWOut` struct (populated if `mgCtl.pooled=1`) |
| `panel_var`, `time_var`, `y_varname`, `model`, `mg_vars`, `csa_vars` | Model metadata |
| `nobs`, `ngroups`, `obs_grp`, `k_reg`, `df`, `df_csa`, `csa_lags` | Dimensions |

### `pcceNWOut` — Pooled CCE results
`b_pcce`, `se_pcce`, `tvalue_pcce`, `out_pcce`, `cov_pcce_rbst_nw`, `cov_pcce_hs2`

---

## Main Public Procedures

### Core estimators (`src/cce_mg.src`)

| Procedure | Description |
|-----------|-------------|
| `mg(data, [mgCtl])` | Mean Group estimator |
| `cce_mg(data, [mgCtl])` | CCE Mean Group estimator |
| `dcce_mg(data, [mgCtl])` | Dynamic CCE-MG (y lags + CSA lags) |
| `pcceNW(varname, y, x, h, n, tivec, b_vec, [PT])` | Pooled CCE with NW SE |
| `cdtest(e, n, tivec, starttvec, endtvec)` | Pesaran (2004) CD test |
| `getCSA(data, groupvar, timevar, [x_csa])` | Compute cross-sectional averages |

### Utility (`src/dcceutil.src`)

| Procedure | Description |
|-----------|-------------|
| `mgControlCreate()` | Create default `mgControl` struct |
| `coeftable(mgO)` | Return k×4 numeric matrix [coef, se, t, p] |
| `__print_mg_output(mgO)` | Print formatted results table |

### Diagnostics & tests (`src/slopehomo.src`, `src/cips.src`, `src/diagnostics.src`)

| Procedure | Description |
|-----------|-------------|
| `slopehomo(mgO)` | Pesaran-Yamagata (2008) Δ and Δ_adj slope homogeneity tests |
| `print_slopehomo(delta, pval, delta_adj, pval_adj)` | Print slope homogeneity results |
| `cips(data, [p, demean])` | Pesaran (2007) CIPS panel unit root test |
| `print_cips(cips_stat, cadf_vec, p, demean)` | Print CIPS results |
| `plotResiduals(mgO)` | 4-panel residual diagnostic plot (scatter, histogram, Q-Q, per-group SD) |
| `cce_rank(data, [mgCtl])` | Rank condition test for CCE (De Vos et al. 2024): checks CSA matrix rank |
| `print_cce_rank(rankO)` | Print formatted rank condition test results |

### Bias correction (`src/bias_correct.src`)

| Procedure | Description |
|-----------|-------------|
| `hpj(data, mgCtl, [estimator_type])` | Half-panel jackknife bias correction (Dhaene & Jochmans 2015) |

### Export & inference (`src/latex_export.src`, `src/bootstrap.src`)

| Procedure | Description |
|-----------|-------------|
| `mgOutToLatex(mgO, filename, [se_type, note])` | Export single model to .tex |
| `mgOutToLatexMulti(mgO_arr, labels, filename, [note])` | Export 2–6 models side-by-side |
| `mgBootstrap(data, mgCtl, [B, estimator_type])` | Wild bootstrap SEs (Rademacher weights) |

### Extended estimator flags (built into `cce_mg`, `dcce_mg`)

| Flag | Ref. | Description |
|------|------|-------------|
| `ctl.i1 = 1` | Kapetanios, Pesaran & Yamagata (2011) | Adds Δȳₜ and Δx̄ₜ to h matrix; ensures CCE consistency with I(1) common factors |
| `ctl.two_way = 1` | Bai (2009) | Time-demeans y, x, and x_csa before CCE; handles both time-specific and cross-sectional common shocks |

---

## Data Format Convention

All estimator procedures expect `data` in this column order:
1. **Group variable** (panel ID) — column 1
2. **Time variable** — column 2
3. **Dependent variable (y)** — column 3
4. **Independent variables (x)** — columns 4+

Sort by group then time. Use `packr()` to remove missing rows before calling.

**Alternative: specify variables by name** using `mgControl.y_var` and `mgControl.x_vars`.
The library will reorder columns automatically before estimation:

```gauss
ctl = mgControlCreate();
ctl.y_var  = "log_rgdpo";
ctl.x_vars = "log_ck" $| "log_ngd";
cceO = cce_mg(mydata, ctl);   // columns can be in any order
```

---

## Bug History (Resolved)

All previously identified bugs are fixed. Summary:

| File | Bug | Fix |
|------|-----|-----|
| `dcceutil.src` `__getPT` | `T^1/3` was T/3 (precedence) | `T^(1/3)` |
| `dcceutil.src` `_getTimeInfo` | Hardcoded `+3+` lag offset | Use `+PT+` |
| `dcceutil.src` `_getTimeInfo` | Hardcoded `"year"` column | Use `timevar` |
| `dcceutil.src` `__getylag` | Loop `for i(1,1,rows)` ran once | `for i(1,rows,1)` |
| `dcceutil.src` `__getylag` | `tivec[i]:tivec[i+1]` off-by-one | `tivec[i]+1:tivec[i+1]` |
| `dcceutil.src` `__pdbalanced` | `rows(t)==1` lowercase `t` undefined | `rows(T)==1` |
| `cce_mg.src` `cdtest` | Bare debug expressions printed 4278 lines | Removed |
| `cce_mg.src` `pcceNW` | Self-assignment `cov=cov`; 3 loops | Fixed; reduced to 2 loops |
| `cce_mg.src` `cce_mg` | Wrong model label "Mean Group" | "Common Correlated Effects..." |
| Paths in `.e` files | Backslash escape sequences | Forward slashes everywhere |

---

## Validation

All three estimators verified against R `plm::pmg()` to 6 decimal places:
- **MG:** log_ck=0.305300, log_ngd=0.279783, inpt=5.391778 ✅
- **CCE-MG:** log_ck=0.316743, log_ngd=0.089055, inpt=1.145539 ✅
- **DCCE-MG:** y_l=0.456422, log_ck=0.153173, log_ngd=0.009159 ✅

Run validation:
```
# R
"C:/Program Files/R/R-4.5.0/bin/Rscript.exe" validation/validate_dcce.R

# GAUSS (from repo root)
C:\gauss26\tgauss.exe -b -nj validation\validate_gauss.e
```

---

## Penn World Tables Example Data (`examples/penn_sample.dta`)

Key variables: `id` (country), `year`, `log_rgdpo`, `log_ck`, `log_ngd`, `log_hc`

Typical model: `log_rgdpo ~ log_ck + log_ngd`, with `log_hc` as extra CSA variable.
