# dccelib v1.2.0 — Initial Release

A GAUSS library for panel data estimation in the presence of cross-sectional dependence, implementing the Common Correlated Effects (CCE) framework of Pesaran (2006) and a suite of companion diagnostic and inference tools.

---

## Background

Standard panel estimators — pooled OLS, fixed effects, random effects — assume error terms are independent across units after controlling for observed covariates. In practice this assumption fails whenever unobserved common factors (global business cycles, commodity price shocks, financial contagion, technology diffusion) affect multiple units simultaneously. The consequences are severe: standard errors are understated, t-statistics are inflated, and coefficient estimates may themselves be inconsistent when the common factors are correlated with the regressors.

The CCE estimator (Pesaran 2006) resolves this by augmenting each unit's regression with cross-sectional averages of the dependent variable and regressors, which serve as observable proxies for the unobserved factors. No factor estimation or assumption about the number of factors is required.

dccelib brings this framework — and its extensions — to GAUSS with a validated, struct-based API and a complete suite of companion tools.

---

## Estimators

| Procedure | Description |
|-----------|-------------|
| `mg()` | Mean Group estimator (Pesaran & Smith 1995) |
| `cce_mg()` | CCE Mean Group estimator (Pesaran 2006) |
| `dcce_mg()` | Dynamic CCE-MG with lagged y and CSA lags (Chudik & Pesaran 2015) |
| `pcce_mg()` | PC-CCE-MG: PCA-based factor proxies via SVD (no external dependency) |

All estimators support:
- **Formula-string API**: `mg(data, "y ~ x1 + x2")` or `ctl.formula = "y ~ x1 + x2"` — no manual column selection required
- **Inline CSA specification**: `csa()` inside the formula string — `"y ~ x1 + x2 + csa(z)"` — eliminates the need for a separate `x_csa_names` field. Multiple variables: `csa(z1, z2)`.
- **Named variable selection**: `ctl.y_var`, `ctl.x_vars`, `ctl.x_csa_names`, `ctl.groupvar`, `ctl.timevar`
- **I(1) extension** (`ctl.i1 = 1`): adds differenced CSAs for integrated regressors (Kapetanios, Pesaran & Yamagata 2011)
- **Two-way CCE** (`ctl.two_way = 1`): time-demeaning for two-way factor structures (Bai 2009)
- **Pooled CCE** (`ctl.pooled = 1`): estimates pooled CCE with Newey-West SE in the same call

---

## Diagnostic Tests

| Procedure | Description |
|-----------|-------------|
| `cips()` / `cips_test()` | Pesaran (2007) CIPS panel unit root test |
| `slopehomo()` | Pesaran-Yamagata (2008) Δ and Δ_adj slope homogeneity tests |
| `cce_rank()` | De Vos, Everaert & Sarafidis (2024) CCE rank condition test |
| `westerlundTest()` | Westerlund (2007) ECM-based panel cointegration test (Gₜ, Gₐ, Pₜ, Pₐ) |

---

## Post-Estimation Tools

| Procedure | Description |
|-----------|-------------|
| `hpj()` | Half-panel jackknife bias correction (Dhaene & Jochmans 2015) |
| `mgBootstrap()` / `mgBootstrapSE()` | Wild Rademacher bootstrap standard errors |
| `longRunMG()` | Delta-method long-run multipliers from DCCE-MG |

---

## Visualization

| Procedure | Description |
|-----------|-------------|
| `plotResiduals()` | 4-panel residual diagnostic plot |
| `plotCoefficients()` | Caterpillar plot of per-group slope estimates |
| `plotResidualACF()` | Sample ACF of pooled residuals |

---

## Output and Export

| Procedure | Description |
|-----------|-------------|
| `printCoefCompare()` | Aligned side-by-side coefficient comparison table |
| `coeftable()` | Extract k×4 numeric matrix [coef, se, t, p] |
| `mgOutToLatex()` | Export single model to LaTeX `tabular` |
| `mgOutToLatexMulti()` | Export 2–6 models side by side in one LaTeX table |

---

## Validation

All three core estimators verified against R `plm::pmg()` to 6 decimal places on Penn World Tables data (N=93 countries, T≈50 years):

| Estimator | log_ck | log_ngd | Intercept / y_lag |
|-----------|--------|---------|-------------------|
| MG | 0.305300 | 0.279783 | 5.391778 |
| CCE-MG | 0.316743 | 0.089055 | 1.145539 |
| DCCE-MG | 0.153173 | 0.009159 | y_lag: 0.456422 |

---

## Examples

Nine example scripts in `examples/` cover the complete recommended workflow:

`mg_penn.e` · `cce_penn.e` · `dcce_penn.e` · `cce_proc.e` · `diagnostics.e` · `advanced_cce.e` · `bias_correction.e` · `export_tables.e` · `pca_cce.e`

---

## Requirements

GAUSS 26+. No external library dependencies.

---

## References

- Pesaran (2006), *Econometrica* 74(4): 967–1012
- Pesaran & Smith (1995), *Journal of Econometrics* 68(1): 79–113
- Chudik & Pesaran (2015), *Journal of Econometrics* 188(2): 393–420
- Pesaran (2007), *Journal of Applied Econometrics* 22(2): 265–312
- Pesaran & Yamagata (2008), *Journal of Econometrics* 142(1): 50–93
- Kapetanios, Pesaran & Yamagata (2011), *Journal of Econometrics* 160(2): 326–348
- Dhaene & Jochmans (2015), *Review of Economic Studies* 82(3): 991–1030
- Westerlund (2007), *Oxford Bulletin of Economics and Statistics* 69(6): 709–748
- De Vos, Everaert & Sarafidis (2024)
- Ahn & Horenstein (2013), *Econometrica* 81(3): 1203–1227
- Bai (2009), *Econometrica* 77(4): 1229–1279
