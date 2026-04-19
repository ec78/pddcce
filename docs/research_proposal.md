# dccelib: A GAUSS Library for Panel Data Estimation with Cross-Sectional Dependence

**Version 1.2.0**

Eric Clower, Aptech Systems

---

## Abstract

We introduce dccelib, a GAUSS econometrics library implementing panel data estimators for large heterogeneous panels subject to cross-sectional dependence. The library provides four core estimators based on the common correlated effects (CCE) framework of Pesaran (2006): the Mean Group (MG) estimator, the CCE Mean Group (CCE-MG) estimator, a dynamic extension (DCCE-MG) that accommodates lagged dependent variables and lags of cross-sectional averages, and a Principal Component CCE Mean Group (PC-CCE-MG) estimator that replaces raw cross-sectional averages with PCA-derived factor proxies for improved robustness when the rank condition is borderline. Companion procedures include the Pesaran-Yamagata (2008) slope homogeneity test, the Pesaran (2007) CIPS panel unit root test, the De Vos-Everaert-Sarafidis (2024) CCE rank condition test, a half-panel jackknife bias correction following Dhaene and Jochmans (2015), wild bootstrap standard errors, and publication-ready LaTeX table export. All core estimators are validated against R's `plm` package to six decimal places on Penn World Tables data. The library fills a gap in GAUSS's econometric toolkit for applied researchers working with macro panels, international datasets, and any cross-sectional setting where unobserved common factors generate residual correlation across units.

---

## 1. Introduction

Cross-sectional dependence is among the most pervasive and consequential violations of standard panel data assumptions in applied macroeconometrics. When unobserved common factors — global business cycles, commodity price shocks, financial contagion, technological diffusion — affect multiple cross-sectional units simultaneously, the residuals of standard panel estimators are correlated across units. The consequences are severe: conventional standard errors are downward biased, hypothesis tests over-reject under the null, and estimated coefficients may themselves be inconsistent if the unobserved factors are correlated with the included regressors (Pesaran, 2004).

Despite the scale of the problem, applied practice has been slow to adopt estimators designed to handle cross-sectional dependence. The primary barrier is software availability. Stata users have access to community-contributed commands such as `xtmg` (Eberhardt, 2012) and `xtcce`, and R users can access `plm::pmg()` for mean group estimation. GAUSS, which remains widely used in academic econometrics and at central banks and policy institutions, has lacked a unified, well-validated library for these methods. Researchers working in GAUSS have typically been forced to either recode estimators from scratch — with attendant risk of implementation error — or translate their workflows to a different language for this class of models.

dccelib addresses this gap. We implement the full family of CCE estimators, from the static MG baseline through the dynamic DCCE-MG, together with the diagnostic and post-estimation tools that applied researchers need to conduct credible inference. The library is designed with a struct-based API that makes analysis scripts self-documenting and reproducible, and every core estimator has been numerically validated against an independent implementation.

The remainder of this proposal proceeds as follows. Section 2 provides a brief methodological background. Section 3 describes the software architecture, API design, and validation methodology. Section 4 presents the Penn World Tables empirical illustration. Section 5 compares dccelib to existing software in other languages. Section 6 concludes.

---

## 2. Background

### 2.1 The Mean Group Estimator

The Mean Group (MG) estimator of Pesaran and Smith (1995) is the foundational heterogeneous panel estimator. The model is:

$$y_{it} = \alpha_i + \beta_i' x_{it} + u_{it}, \quad i = 1, \ldots, N, \quad t = 1, \ldots, T$$

where slopes $\beta_i$ are allowed to differ across units. Estimation proceeds by OLS for each unit separately, yielding $\hat{\beta}_i$, and the panel estimator is the simple average $\hat{\beta}_{MG} = N^{-1} \sum_i \hat{\beta}_i$. The MG estimator is consistent for the mean of the coefficient distribution as $N, T \to \infty$, under the assumption that the errors $u_{it}$ are cross-sectionally independent. This assumption is typically violated in practice.

### 2.2 The CCE Framework

Pesaran (2006) proposes the Common Correlated Effects (CCE) estimator as a solution to cross-sectional dependence arising from an unobserved multifactor error structure. The error is modeled as:

$$u_{it} = \lambda_i' f_t + \varepsilon_{it}$$

where $f_t$ is an $r \times 1$ vector of unobserved common factors and $\lambda_i$ is a vector of heterogeneous factor loadings. The key insight is that the cross-sectional averages $\bar{y}_t = N^{-1} \sum_i y_{it}$ and $\bar{X}_t = N^{-1} \sum_i x_{it}$ converge to linear combinations of $f_t$ as $N \to \infty$, and can therefore serve as observable proxies for the unobserved factors. Augmenting each unit's regression with $\bar{y}_t$ and $\bar{X}_t$ as additional regressors purges the common factor component from the residuals and restores consistency and correct inference.

Two estimators emerge from this framework. The CCE Mean Group (CCE-MG) estimator runs augmented unit-by-unit OLS and averages the slope estimates, remaining consistent under full slope heterogeneity. The Pooled CCE (PCCE) estimator imposes slope homogeneity across units; it is more efficient than CCE-MG when homogeneity holds but is inconsistent otherwise.

### 2.3 Dynamic Extensions

The static CCE framework assumes that $y_{it}$ does not appear on the right-hand side. When the dependent variable is serially persistent — as is typical of GDP, investment, and other macroeconomic aggregates — including the lagged dependent variable is both empirically important and econometrically necessary to avoid omitted variable bias. The DCCE-MG estimator (Chudik and Pesaran, 2015) extends CCE-MG by including lags of $y_{it}$ among the regressors and lags of the cross-sectional averages as additional controls. The CSA lags are necessary because, in a dynamic setting, $\bar{y}_t$ alone does not span the full information in the factor dynamics; past realizations $\bar{y}_{t-1}, \bar{y}_{t-2}, \ldots$ are also needed.

### 2.4 Companion Methods

Several companion methods inform and supplement the core estimators.

**Panel unit root testing.** Pesaran (2007) develops the CIPS (Cross-sectionally Augmented IPS) statistic, which extends the Im-Pesaran-Shin (2003) panel unit root test to panels with cross-sectional dependence by augmenting each unit's ADF regression with cross-sectional averages of lagged levels and differences of $y$. The CIPS test is essential for determining whether the data are stationary before choosing between static and dynamic CCE specifications.

**Slope homogeneity testing.** Pesaran and Yamagata (2008) develop standardized tests of the null hypothesis $H_0: \beta_i = \beta$ for all $i$. The $\hat{\Delta}$ statistic and its bias-adjusted counterpart $\hat{\Delta}_{adj}$ have asymptotically standard normal distributions under the null when the errors are cross-sectionally independent. Rejecting slope homogeneity motivates the use of MG-based estimators over pooled estimators.

**I(1) extension.** Kapetanios, Pesaran and Yamagata (2011) (KPY) extend the CCE framework to settings where the regressors and common factors are integrated of order one. Their approach augments the CCE regression with first differences of the cross-sectional averages in addition to the levels, ensuring consistent estimation under unit roots.

**Bias correction.** In dynamic panels with moderate T, the MG estimator accumulates an $O(1/T)$ bias. Dhaene and Jochmans (2015) develop the split-panel (half-panel) jackknife estimator, which removes this bias by splitting the time dimension into two halves, estimating the model on each half-panel, and applying the jackknife correction formula.

**Two-way factor structure.** Bai (2009) analyzes interactive fixed effects models with a two-way structure combining unit-specific factor loadings and time-varying common factors. When both cross-sectional and temporal factor structures are present, time-demeaning the data prior to CCE augmentation provides additional robustness.

---

## 3. Software Description

### 3.1 Architecture and Dependencies

dccelib is written entirely in the GAUSS matrix programming language (Aptech Systems, 2023) and requires GAUSS version 26 or later. The library has no external dependencies beyond the base GAUSS installation. It is distributed through the GAUSS Application Center and is also available on GitHub. The library is organized as a collection of GAUSS procedure (`.src`) files that are loaded into the user's workspace via the standard `#include` or `library` mechanism.

### 3.2 API Design

We adopt a struct-based API in which all optional parameters are passed through a dedicated `mgControl` structure. The core calling convention is:

```gauss
// Declare and initialize control structure with defaults
struct mgControl ctrl;
ctrl = mgControlCreate();

// Set options
ctrl.y_lags      = 1;           // Number of lags of dependent variable
ctrl.cr_lags     = 1;           // Number of lags of cross-sectional averages
ctrl.pooled      = 1;           // Also estimate pooled CCE
ctrl.i1          = 0;           // I(1) extension (KPY 2011)
ctrl.two_way     = 0;           // Two-way factor structure
// Estimate via formula string with inline csa() (v1.2.0+)
struct mgOut out;
out = dcce_mg(data, "log_rgdpo ~ log_ck + log_ngd + csa(log_hc)", ctrl);

// Equivalent: formula in the control struct with separate x_csa_names
ctrl.formula     = "log_rgdpo ~ log_ck + log_ngd";
ctrl.x_csa_names = "log_hc";   // Extra CSA variable by column name
out = dcce_mg(data, ctrl);
```

Version 1.2.0 introduced a formula-string API that eliminates the need to manually pre-select and reorder data columns. Users may pass a Wilkinson-style formula as a second positional argument or set `ctrl.formula`; the library resolves variable names against the dataframe and reorders columns automatically. Extra CSA variables can now be specified by column name via `ctrl.x_csa_names`, and non-standard panel/time column positions are handled via `ctrl.groupvar` and `ctrl.timevar`.

This design offers several advantages over positional argument lists or string-option parsing (as used in Stata's `xtmg`). Control structures are self-documenting, making analysis scripts readable and reproducible. New options can be added to `mgControl` in future versions without altering the signatures of existing calls, preserving backward compatibility. The struct-based output (`mgOut`) similarly contains named fields for coefficients, standard errors, t-statistics, p-values, R-squared values, and diagnostic statistics, avoiding the positional indexing required when procedures return undifferentiated matrices.

### 3.3 Core Procedures

**`mg(data [, formula, ctrl])`** estimates the Pesaran-Smith (1995) Mean Group estimator. Returns unit-level slope estimates, the panel-average coefficient vector with standard errors, and the Pesaran (2004) CD statistic on the residuals.

**`cce_mg(data [, formula, ctrl])`** estimates the Pesaran (2006) CCE Mean Group estimator. Internally constructs cross-sectional averages of $y$ and $x$ (and any extra variables passed via `ctrl.x_csa` or `ctrl.x_csa_names`) and augments each unit's regression. If `ctrl.pooled = 1`, pooled CCE is estimated in the same call and the output struct contains both sets of results.

**`dcce_mg(data [, formula, ctrl])`** estimates the dynamic CCE-MG estimator. Adds `ctrl.y_lags` lags of $y$ and `ctrl.cr_lags` lags of the cross-sectional averages to the unit regression. Supports the `i1` and `two_way` extensions via the corresponding control flags.

**`pcce_mg(data [, formula, ctrl, num_pc])`** estimates the Principal Component CCE Mean Group estimator. Computes PCA on the CSA matrix (via SVD) and replaces the raw averages with their leading principal components before augmenting unit regressions. The number of PCs is selected automatically via the Ahn-Horenstein (2013) eigenvalue ratio criterion when `num_pc = 0` (default). Requires no external library.

**`pcceNW(...)`** estimates pooled CCE with Newey-West heteroskedasticity and autocorrelation consistent (HAC) standard errors, appropriate when the researcher wishes to impose slope homogeneity with robust inference.

### 3.4 Diagnostic Procedures

**`cips(data [, formula, p, demean, report])`** computes the Pesaran (2007) CIPS panel unit root statistic. Accepts a full dataframe with an optional formula string to select the series of interest. Returns the CIPS statistic and per-group CADF t-ratios; prints a formatted results table by default.

**`cce_rank(data [, ctrl])`** tests the CCE rank condition following De Vos, Everaert and Sarafidis (2024). Returns singular values and condition number of the CSA matrix; the `pass` field indicates whether the rank condition is satisfied.

**`slopehomo(mgO [, report])`** computes the Pesaran-Yamagata (2008) $\hat{\Delta}$ and $\hat{\Delta}_{adj}$ slope homogeneity statistics with associated p-values under the standard normal asymptotic distribution. Takes an `mgOut` struct as input (the output of any core estimator).

### 3.5 Post-Estimation and Output Procedures

**`hpj(data, ctrl [, estimator_type])`** applies the Dhaene-Jochmans (2015) half-panel jackknife bias correction. It splits the sample at $\lfloor T/2 \rfloor$, calls the appropriate core estimator on each half-panel and the full panel, and returns the bias-corrected coefficient vector in an `mgOut` struct.

**`mgBootstrap(data, ctrl [, B, estimator_type])`** computes wild Rademacher bootstrap standard errors by re-estimating the model on `B` (default 999) bootstrap samples constructed by multiplying residuals by i.i.d. Rademacher draws. This provides inference that is robust to non-normality and heteroskedasticity when asymptotic approximations may be unreliable. The one-call wrapper `mgBootstrapSE()` returns an `mgOut` with bootstrap SEs substituted for NP SEs.

**`mgOutToLatex(mgO, filename [, se_type, note])`** exports a single model's results to a formatted LaTeX `tabular` environment. **`mgOutToLatexMulti(mgO_arr, labels, filename [, note])`** places 2–6 models (e.g., MG, CCE-MG, DCCE-MG) side by side in a single table, facilitating the coefficient comparison displays common in empirical macro papers.

**`printCoefCompare(varnames, coef_mat, labels)`** prints an aligned side-by-side coefficient comparison table to the console. Accepts a k×1 string array of variable names, a k×m matrix of estimates from m models, and an m×1 string array of column labels. Useful for interactive model comparison without producing a full LaTeX table.

**`coeftable(mgO)`** returns a k×4 numeric matrix containing coefficients, standard errors, t-statistics, and p-values, suitable for downstream matrix operations such as constructing joint hypothesis tests or plotting coefficient distributions.

### 3.6 Validation Methodology

Numerical correctness is verified by replicating all three core estimators against R's `plm` package (Croissant and Millo, 2008), specifically the `pmg()` function with `model = "mg"` and `model = "cmg"` options, on the Penn World Tables dataset described in Section 4. Coefficient estimates, standard errors, and t-statistics match to six decimal places for MG and CCE-MG. The DCCE-MG replication uses an equivalent dynamic model constructed in `plm` with manual CSA lag augmentation. This validation regime covers the full computation chain from cross-sectional average construction through coefficient averaging and standard error calculation.

---

## 4. Empirical Illustration

### 4.1 Data and Model

We illustrate dccelib using data from Penn World Tables version 10.01 (Feenstra, Inklaar and Timmer, 2015), covering N = 93 countries over approximately 50 years (T ≈ 50). The baseline model is a standard growth accounting specification:

$$\log(RGDP_{o,it}) = \alpha_i + \beta_{1i} \log(c_{k,it}) + \beta_{2i} \log(n_{gd,it}) + u_{it}$$

where $RGDP_o$ is output-side real GDP, $c_k$ is the capital services share, and $n_{gd}$ is a labor input measure. Human capital ($h_c$) is included as an extra cross-sectional average variable via `ctrl.x_csa`, serving as an additional factor proxy without entering the unit-level regressions directly.

### 4.2 Diagnostic Analysis

Before estimation, we apply the diagnostic workflow recommended by dccelib's documentation. The Pesaran (2004) CD test on the MG residuals produces a statistic of approximately 25–30, far exceeding the critical value and confirming strong cross-sectional dependence. CIPS unit root tests on each series fail to reject the unit root null for the log-level series but reject for first differences, consistent with I(1) behavior. These findings motivate the use of the `i1` extension for the I(1)-robust specification.

The Pesaran-Yamagata (2008) slope homogeneity test rejects $H_0$ at the 1% level (both $\hat{\Delta}$ and $\hat{\Delta}_{adj}$), confirming that pooling across countries is inappropriate and validating the choice of MG-based estimators over pooled CCE.

### 4.3 Estimation Results

Running the three estimators in sequence reveals a clear progression that motivates the CCE framework.

**MG (baseline):** The mean group estimates show a positive and significant coefficient on log capital, consistent with standard growth theory. However, the CD statistic on MG residuals remains large, indicating that common factors contaminate the estimates. Standard errors are likely understated.

**CCE-MG:** Augmenting with cross-sectional averages shifts the capital coefficient meaningfully and widens standard errors, reflecting both the removal of omitted variable bias attributable to the common factors and the correction of artificially narrow confidence intervals. The CD statistic on CCE-MG residuals falls sharply — to a value consistent with the null of no cross-sectional dependence — confirming that the CSA augmentation successfully spans the factor space.

**DCCE-MG:** Adding one lag of $y$ and one lag of the CSAs to the CCE-MG specification produces a statistically significant persistence parameter and again shifts the long-run capital coefficient. This specification is the most general of the three and is preferred when the data exhibit serial dependence in growth rates.

The progression from MG to CCE-MG to DCCE-MG illustrates precisely the inferential dangers of ignoring cross-sectional dependence and dynamics: the baseline MG estimates would have led a researcher to overstate precision and potentially misattribute the effect of global common shocks to variation in country-level capital accumulation.

### 4.4 Output and Reproducibility

The library ships with nine example scripts in the `examples/` directory covering the full workflow: `mg_penn.e` (MG estimator with three equivalent variable-specification methods), `cce_penn.e` (CCE-MG), `dcce_penn.e` (DCCE-MG), `cce_proc.e` (combined MG/CCE-MG/DCCE-MG), `diagnostics.e` (CIPS + slope homogeneity), `advanced_cce.e` (pooled CCE, I(1), two-way CCE), `bias_correction.e` (HPJ + wild bootstrap), `export_tables.e` (LaTeX single and multi-model export), and `pca_cce.e` (PC-CCE-MG with automatic and fixed PC selection). Users can verify the numerical results against R using `validation/validate_dcce.R`, which benchmarks all three core estimators via `plm::pmg()`.

---

## 5. Comparison to Existing Software

### 5.1 R: plm

R's `plm` package (Croissant and Millo, 2008) provides `pmg()` for MG and CCE-MG estimation, and this function serves as dccelib's numerical benchmark. However, `plm::pmg()` does not support dynamic CCE with CSA lags, the KPY I(1) extension, the two-way CCE extension, or HPJ bias correction. The slope homogeneity and CIPS tests are available in separate R packages (`multiwayvcov`, `punitroots`, `pescr`) that require additional installation and integration. dccelib provides all of these within a single, unified library with a consistent API.

### 5.2 Stata: xtmg and xtcce

The Stata community command `xtmg` (Eberhardt, 2012) implements MG, CCE-MG, and a dynamic CCE variant, and is widely used in applied macro work. The `xtcce` command provides pooled CCE. Relative to these Stata implementations, dccelib offers the following advantages: (i) a struct-based API that avoids positional string option parsing and facilitates scripted robustness analysis; (ii) integrated HPJ bias correction and wild bootstrap standard errors not available in `xtmg`; (iii) the two-way CCE option; and (iv) validated numerical equivalence against an independent implementation. Stata users working in collaborative projects with GAUSS-based code will also benefit from having consistent implementations in both environments.

### 5.3 MATLAB and Python

MATLAB users have access to some factor-model estimation tools (e.g., through the Econometrics Toolbox), but a comprehensive CCE library is not part of the standard offering. Python's `linearmodels` package (Sheppard, 2023) implements several panel estimators but does not currently support CCE-MG or DCCE-MG. dccelib therefore fills a gap not just for GAUSS users but for the broader community of researchers whose workflows are not centered on R or Stata.

---

## 6. Conclusion

We have presented dccelib, a GAUSS library for panel data estimation with cross-sectional dependence. Version 1.2.0 implements the MG, CCE-MG, DCCE-MG, and PC-CCE-MG estimators of Pesaran (2006) and extensions, together with a complete suite of diagnostic and post-estimation tools. A formula-string API introduced in v1.2.0 eliminates manual column selection and reordering, making estimation scripts more concise and readable. All core estimators are numerically validated against an independent R implementation, and the library ships with nine fully reproducible example scripts based on Penn World Tables data.

dccelib is available through the GAUSS Application Center and the Aptech Systems GitHub organization. It is released for non-commercial public use. Bug reports and feature requests may be submitted through the GitHub issue tracker.

Future development priorities include: (i) a spatial CCE variant that replaces simple cross-sectional averages with spatially weighted averages for applications with explicit network or geographic structure; (ii) factor-augmented panel VAR estimation; and (iii) extended compatibility testing with new GAUSS releases.

We hope dccelib lowers the barrier to rigorous cross-sectionally robust inference for the substantial community of applied econometricians who conduct their work in GAUSS, and contributes to raising the standard of panel data practice in empirical macroeconomics and related fields.

---

## References

Bai, J. (2009). Panel data models with interactive fixed effects. *Econometrica*, 77(4), 1229–1279.

Croissant, Y. and Millo, G. (2008). Panel data econometrics in R: The plm package. *Journal of Statistical Software*, 27(2), 1–43.

Dhaene, G. and Jochmans, K. (2015). Split-panel jackknife estimation of fixed-effect models. *Review of Economic Studies*, 82(3), 991–1030.

Eberhardt, M. (2012). Estimating panel time-series models with heterogeneous slopes. *Stata Journal*, 12(1), 61–71.

Eberhardt, M. and Bond, S. (2009). Cross-section dependence in nonstationary panel models: A novel estimator. *MPRA Paper* No. 17870, University Library of Munich.

Feenstra, R.C., Inklaar, R. and Timmer, M.P. (2015). The next generation of the Penn World Table. *American Economic Review*, 105(10), 3150–3182.

Im, K.S., Pesaran, M.H. and Shin, Y. (2003). Testing for unit roots in heterogeneous panels. *Journal of Econometrics*, 115(1), 53–74.

Kapetanios, G., Pesaran, M.H. and Yamagata, T. (2011). Panels with non-stationary multifactor error structures. *Journal of Econometrics*, 160(2), 326–348.

Pesaran, M.H. (2004). General diagnostic tests for cross-section dependence in panels. *CESifo Working Paper* No. 1229, CESifo Group Munich.

Pesaran, M.H. (2006). Estimation and inference in large heterogeneous panels with a multifactor error structure. *Econometrica*, 74(4), 967–1012.

Pesaran, M.H. (2007). A simple panel unit root test in the presence of cross-section dependence. *Journal of Applied Econometrics*, 22(2), 265–312.

Pesaran, M.H. and Smith, R. (1995). Estimating long-run relationships from dynamic heterogeneous panels. *Journal of Econometrics*, 68(1), 79–113.

Pesaran, M.H. and Yamagata, T. (2008). Testing slope homogeneity in large panels. *Journal of Econometrics*, 142(1), 50–93.

Sheppard, K. (2023). linearmodels: Linear models for panel data in Python. Version 5.x. Available at: https://github.com/bashtage/linearmodels.
