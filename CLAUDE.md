# CLAUDE.md — GAUSS dccelib Project Context

## Project Overview

**dccelib** is a GAUSS econometrics library implementing panel data estimators based on:

> Pesaran, M.H. (2006), "Estimation and Inference in Large Heterogeneous Panels with a Multifactor Error Structure," *Econometrica*, Vol. 74, No. 4, pp. 967–1012.

The library provides three main estimators:
- **MG** — Mean Group Estimator
- **CCE-MG** — Common Correlated Effects Mean Group Estimator
- **DCCE-MG** — Dynamic Common Correlated Effects Mean Group Estimator

**Author:** Eric Clower, Aptech Systems, Inc. (`eric@aptech.com`)
**Target:** GAUSS 23+
**License:** Non-commercial public use only

---

## Repository Structure

```
pddcce/
├── CLAUDE.md               # This file
├── README.md               # User-facing overview
├── package.json            # GAUSS package manifest (name: dccelib, version: 0.1.0)
├── ccep.src                # Legacy pooled CCE code (older style, not in package.json src list)
├── src/
│   ├── cce.sdf             # Structure definitions (mgControl, mgOut, cdOut, pcceNWOut)
│   ├── cce_mg.src          # Main estimator procedures
│   └── dcceutil.src        # Utility/helper procedures and output formatting
└── examples/
    ├── mg_penn.e           # MG estimator example (Penn World Tables data)
    ├── cce_penn.e          # CCE-MG estimator example
    ├── dcce_penn.e         # DCCE-MG estimator example
    ├── cce_proc.e          # Older-style example using direct #include
    ├── penn_sample.dta     # Penn World Tables sample data (Stata format)
    └── jasa2.dta           # Additional dataset
```

---

## GAUSS Language Reference

For writing and reviewing GAUSS code, use the LLM reference at:
**https://github.com/aptech/gauss-llm-reference**

Key GAUSS conventions used in this codebase:
- Procedures use `proc(N) = name(args);` ... `retp(...);` ... `endp;`
- Structs declared in `.sdf` files, instantiated with `struct TypeName varname;`
- Optional args via `dynargsGet(index, default)`
- `do while` / `endo` loops; `for i(start, end, step)` / `endfor`
- String arrays use `$~` (horizontal concat) and `$|` (vertical concat)
- Element-wise operations: `.*`, `./`, `.^`, `.==`, etc.
- `pinv()` for pseudoinverse, `invpd()` for positive-definite inverse, `inv()` for general inverse
- `packr()` removes rows with missing values
- `asDF()` converts matrix to dataframe; `setcolnames()` / `getcolnames()` manage column metadata
- `aggregate()` for group-level aggregation; `unique()` for unique values
- `indnv()` finds index of elements; `counts()` counts occurrences

---

## Data Structures (`src/cce.sdf`)

### `mgControl` — Input control structure
| Member | Type | Default | Description |
|--------|------|---------|-------------|
| `zero_x` | matrix | 0 | Column identifier for vars with zeros (normalization) |
| `zero_id` | matrix | 0 | Group identifier for normalized groups |
| `y_lags` | matrix | 0 | Include lagged y (1 = yes) |
| `cr_lags` | matrix | 0 | Number of CSA lags (0 = auto) |
| `no_xbar` | matrix | 0 | Column IDs to exclude from cross-sectional averages |
| `x_common` | matrix | 0 | Common regressors (not group-specific) |
| `x_csa` | matrix | 0 | Extra vars for CSA not in main regression |
| `report` | matrix | 1 | Print output (1 = yes) |

### `mgOut` — Results output structure
Key members: `b_mg`, `se_mg`, `tvalue`, `pval`, `cov_mg`, `e_mg`, `b_vec`, `se_nw`, `nobs`, `ngroups`, `obs_grp`, `k_reg`, `df`, `df_csa`, `csa_lags`, `cd_stat`, `cd_pval`, `out_mg`, `panel_var`, `time_var`, `y_varname`, `mg_vars`, `csa_vars`, `model`

### `pcceNWOut` — Pooled CCE results
Members: `b_pcce`, `cov_pcce`, `sig2_pcce`, `e_pcce`, `se_pcce`, `tvalue_pcce`, `out_pcce`, `ci_vec`, `sig_i_vec`, `cov_pcce_rbst_nw`, `cov_pcce_hs2`

---

## Main Public Procedures (`src/cce_mg.src`)

### `dcce_mg(data, [mgCtl])` → `mgOut`
Dynamic CCE-MG estimator. Includes lagged y and CSA lags. Sets `model = "Dynamic Common Correlated Effects Estimator"`.

### `mg(data, [mgCtl])` → `mgOut`
Plain Mean Group estimator. No cross-sectional averages.

### `cce_mg(data, [mgCtl])` → `mgOut`
CCE Mean Group estimator. Includes cross-sectional averages in the regression.
> **Bug:** Currently sets `model = "Mean Group Estimator"` — should be `"Common Correlated Effects Mean Group Estimator"`.

### `getCSA(data, groupvar, timevar, [x_csa])` → `data_bar`
Computes cross-sectional averages using `aggregate(..., "mean", timevar)`.

### `pcceNW(varname_cce, y, x, h, n, tivec, starttvec, endtvec, years, b_vec_mg)` → `pcceNWOut`
Pooled CCE estimator with Newey-West standard errors. Older interface — not currently wired into the main workflow.

### `cdtest(e, n, tivec, starttvec, endtvec)` → `{cd_nt, pvalue_cd_nt, meanrho, used}`
Pesaran (2004) CD test for cross-sectional dependence.

---

## Internal/Private Procedures (`src/dcceutil.src`)

| Procedure | Outputs | Description |
|-----------|---------|-------------|
| `mgControlCreate()` | `mgControl` | Creates default control structure |
| `_mg(yname, xname, y, x, n, tivec, zero_x, zero_id, [PT])` | `mgOut` | Core group-by-group OLS loop |
| `__olsIndividual(y_i, x_i, t_i, k, zero_id)` | 6 values | OLS for one group |
| `_getTimeInfo(data, groupvar, timevar, [PT])` | 10 values | Extract time/ID index vectors |
| `_getDataMats(data, x_common, data_bar, no_xbar, ylag)` | 7 values | Build y, x, ybar, xbar matrices |
| `__getHMat(data_d, ybar, xbar, ...)` | `{h, hname}` | Build H matrix (z, ybar, xbar, lags) |
| `__getdlags(ybar, xbar, crlags)` | `{ybarlags, xbarlags, PT}` | Generate CSA lag matrices |
| `__getylag(y, starttvec, endtvec, datevec, isbalanced, n)` | `y_l` | Lag dependent variable |
| `__getPT(T)` | scalar | Auto lag length: `int(T^(1/3))` |
| `__reshapeH(h_i, starttvec, endtvec, datevec, bigT)` | `h` | Align H matrix to panel structure |
| `__getnames(data, x_common, x_csa, ylag)` | 5 strings | Extract variable names |
| `__delXBars(xbar, no_xbar)` | `xbar` | Remove excluded vars from CSA |
| `__pdbalanced(grp)` | scalar | Check if panel is balanced |
| `__print_mg_output(mgO)` | — | Print full results table |
| `__print_mg_header(mgO)` | — | Print model header block |
| `__print_mg_footer(mgO)` | — | Print CD test and variable lists |

---

## Data Format Convention

All estimator procedures expect `data` in this column order:
1. **Group variable** (panel ID) — column 1
2. **Time variable** — column 2
3. **Dependent variable (y)** — column 3
4. **Independent variables (x)** — columns 4+

Data must be sorted by group, then time. Use `packr()` to remove missing rows before calling estimators.

---

## Bugs / Issues

### Remaining legacy issue
- **`ccep.src`** — stale file returning 11 values (old style); doc header incorrectly describes an ADF test; not in `package.json` src list. Do not use.

### Fixed — Phase 1 (bug fixes)

| File | Item | Fix applied |
|------|------|-------------|
| `dcceutil.src` `__getPT` | `T^1/3` → `T^(1/3)` | Was computing T/3 not cube root due to operator precedence |
| `dcceutil.src` `_getTimeInfo` | `+ 3 +` → `+ PT +` | Hardcoded lag offset instead of using variable |
| `dcceutil.src` `_getTimeInfo` | `"year"` → `timevar` | Hardcoded column name broke non-year time variables |
| `dcceutil.src` `__getylag` | `for i(1, 1, rows)` → `for i(1, rows, 1)`; `tivec[i]` → `tivec[i]+1` | Loop ran once only; 0-based index invalid in GAUSS |
| `cce_mg.src` `cce_mg` | Model label | Corrected to `"Common Correlated Effects Mean Group Estimator"` |
| `examples/mg_penn.e` | Path | `"examples/..."` → `__FILE_DIR $+ "..."` |
| `examples/cce_proc.e` | API | Rewrote to use `mgControl` struct; uses Penn data |

### Fixed — Phase 2 (efficiency) + validation bugs

| File | Change |
|------|--------|
| `dcceutil.src` `__olsIndividual` | Compute `pinv(X'X)` once, reuse for `b` and `se` (was computed twice) |
| `cce_mg.src` `pcceNW` | 3 loops → 2; pre-allocated `e_pcce`; fixed `Tivec` case bug; fixed `cov_pcce_rbst_nw` self-assignment; compute `pinv(mx_mx)` once and reuse |
| `cce_mg.src` `cdtest` | Removed bare debug expressions `cumtivec[j]+1;` and `cumtivec[j+1];` that printed N×(N-1)/2 lines of noise (4,278 lines for N=93) |

### Validation status
All coefficients verified against R `plm::pmg()` to 6 decimal places:
- MG: ✅ log_ck=0.305300, log_ngd=0.279783, inpt=5.391778
- CCE-MG: ✅ log_ck=0.316743, log_ngd=0.089055, inpt=1.145539
- DCCE-MG: ✅ y_l=0.456422, log_ck=0.153173, log_ngd=0.009159

---

## Objectives / Development Roadmap

1. **Accuracy & efficiency**: ✅ Phase 1 bugs fixed. ✅ Phase 2 efficiency applied. ✅ Validated against R `plm::pmg()` to 6 decimal places.
2. **Examples**: ✅ All `.e` files updated to current API. Next: add output commentary.
3. **pcceNW integration**: Wire pooled CCE estimator into `cce_mg`/`dcce_mg` workflow. Currently a standalone proc, not called from main estimators.
4. **Documentation**: Produce formatted documentation compatible with GAUSS package standards.
5. **Library expansion**: See expansion roadmap (TBD).

## Validation Scripts
- `validation/validate_dcce.R` — R `plm::pmg()` reference values (run with Rscript.exe)
- `src/validate_gauss.e` — GAUSS batch validation (run: `cd src && C:\gauss26\tgauss.exe -b -nj validate_gauss.e`)

---

## Penn World Tables Example Data (`examples/penn_sample.dta`)

Variables used:
- `id` — country identifier (panel group variable)
- `year` — year (time variable, loaded as date with `date($year, '%Y')`)
- `log_rgdpo` — log real GDP per output
- `log_hc` — log human capital index
- `log_ck` — log physical capital stock
- `log_ngd` — log population growth + break-even investment (5%)

Typical model structure: `log_rgdpo ~ log_ck + log_ngd` with `log_hc` used in CSA.
