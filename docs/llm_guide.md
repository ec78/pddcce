# Using LLMs with dccelib

**dccelib** is a GAUSS library not covered by the training data of most LLMs. This guide
explains how to give any AI assistant enough context to answer questions, write code, and
debug issues for you effectively.

Works with: **Claude**, **ChatGPT**, **Gemini**, **GitHub Copilot**, and any other
assistant that accepts a system prompt or context block.

---

## Quick Start

The fastest path is to paste the [Context Block](#context-block) below at the start of
any conversation, then ask your question. Everything else in this guide makes the answers
more accurate for specific tasks.

---

## Context Block

Copy and paste this entire block at the start of a new chat session before asking your
question. It gives the LLM everything it needs to understand the library.

````
### dccelib context — paste this before asking questions ###

dccelib is a GAUSS 26+ library implementing panel data estimators based on Pesaran (2006).
Install with `library dccelib;`. Source: github.com/eclow/pddcce.

--- DATA FORMAT ---
All estimators expect data as a dataframe with columns in this order:
  [group_var, time_var, y, x1, x2, ...]
Sorted by group then time. Remove missing rows with packr() before calling.
Alternative: set ctl.y_var and ctl.x_vars (string names) and pass data in any column order.

--- THREE CORE ESTIMATORS ---
  mg(data [, ctl])       Mean Group (Pesaran & Smith 1995). No CD correction.
  cce_mg(data [, ctl])   CCE Mean Group (Pesaran 2006). Augments with CSAs.
  dcce_mg(data [, ctl])  Dynamic CCE-MG (Chudik & Pesaran 2015). Adds lagged y + CSA lags.
All return an mgOut struct.

--- mgControl STRUCT (initialize with mgControlCreate()) ---
  ctl.y_lags   = 1;              // lags of y (dcce_mg)
  ctl.cr_lags  = 3;              // CSA lag order (0 = auto via floor(T^(1/3)))
  ctl.x_csa    = data[.,"var"]; // extra vars for CSA (not direct regressors)
  ctl.no_xbar_names = "x1" $| "x2";  // exclude these from CSA
  ctl.pooled   = 1;              // also run pooled CCE → stored in mgOut.pcce
  ctl.i1       = 1;              // add first-differenced CSAs (for I(1) regressors)
  ctl.two_way  = 1;              // time-demean before CCE (two-way FE)
  ctl.report   = 0;              // suppress printed output

--- mgOut STRUCT (key fields) ---
  mgO.b_mg      // k×1 MG coefficient estimates
  mgO.se_mg     // k×1 NP standard errors (Pesaran 2006 eq.58)
  mgO.tvalue    // k×1 t-statistics
  mgO.pval      // k×1 p-values
  mgO.ci        // k×2 confidence intervals [lb, ub]
  mgO.R_sq      // mean within-group R²
  mgO.cd_stat   // Pesaran (2004) CD statistic
  mgO.cd_pval   // CD p-value
  mgO.b_vec     // n×k individual group estimates
  mgO.b_stats   // k×4 heterogeneity: [min, mean, max, sd] across groups
  mgO.out_mg    // formatted output dataframe
  mgO.loglik / mgO.aic / mgO.bic  // model fit statistics

--- COMPANION PROCEDURES ---
  // Bias correction
  hpj(data, ctl [, "mg"|"cce_mg"|"dcce_mg"])   // HPJ bias correction; returns mgOut
  hpjO.b_stats_hpj   // k×4: [b_full, b_h1, b_h2, b_hpj]

  // Bootstrap SEs
  { se_boot, b_boot } = mgBootstrap(data, ctl [, B, estimator_type])
  mgBootstrapSE(data, ctl [, B, estimator_type])  // returns mgOut with bootstrap SEs

  // Long-run multipliers
  lrO = longRunMG(mgO)   // requires dcce_mg with y_lags >= 1
  print_longRun(lrO)
  lrO.lr_coef / lrO.lr_se / lrO.adj_speed

  // Diagnostic tests
  { cips_stat, cadf_vec } = cips(data [, p, demean])  // panel unit root
  print_cips(cips_stat, cadf_vec, p, demean)
  { delta, pval, delta_adj, pval_adj } = slopehomo(mgO)  // slope homogeneity
  wO = westerlundTest(data [, p, demean])              // panel cointegration
  print_westerlund(wO)
  rankO = cce_rank(data [, ctl])                       // CSA rank condition
  coeftable(mgO)   // returns k×4 matrix [coef, se, t, p]

  // Visualization
  plotResiduals(mgO)      // 4-panel residual diagnostic
  plotCoefficients(mgO)   // caterpillar plot of per-group estimates
  plotResidualACF(mgO)    // ACF of pooled residuals

  // LaTeX export
  mgOutToLatex(mgO, "file.tex" [, "np"|"nw", note])
  mgOutToLatexMulti(mgO|cceO|dcceO, "MG"$|"CCE-MG"$|"DCCE-MG", "file.tex")

--- GAUSS LANGUAGE NOTES ---
  String array concatenation: "a" $| "b"   (vertical), "a" $~ "b"  (horizontal)
  String comparison: str1 $== str2  (not ==)
  Struct syntax:  struct mgOut mgO;   mgO = cce_mg(data, ctl);
  Optional args:  extra args accessed via dynargsGet(1, default)
  No any() function — use explicit comparisons or sumc(vec .== val)
  local declarations must come before executable code in a procedure
  Scalar struct fields cause type errors in GAUSS 26 — use matrix type

--- MINIMAL WORKING EXAMPLE ---
  library dccelib;
  fname = "examples/penn_world.dta";
  data  = packr(loadd(fname, ". + date($year, '%Y')"));
  data  = order(data, "id" $| "year");
  reg_data = data[., "id" "year" "log_rgdpo" "log_ck" "log_ngd"];

  struct mgControl ctl;
  ctl = mgControlCreate();
  ctl.x_csa = data[., "log_hc"];

  struct mgOut cceO;
  cceO = cce_mg(reg_data, ctl);

### end dccelib context ###
````

---

## Tips by Task

### Estimating a model

State the dependent variable, regressors, and which estimator you want. Paste the context
block first, then ask:

> *"Using dccelib, estimate log_rgdpo on log_ck and log_ngd with cce_mg. My data is in
> a Stata file called mydata.dta with columns id, year, and those three variables plus
> log_hc for extra CSA. Write the GAUSS code."*

### Interpreting output

Paste the context block, then quote the printed table or copy the relevant field values:

> *"Here is my dcce_mg output: b_mg = [0.456, 0.153, 0.009], se_mg = [0.089, 0.054,
> 0.028], cd_stat = 1.23, cd_pval = 0.22. The regressors are [y_lag, log_ck, log_ngd].
> Explain what these results mean for a policy-focused audience."*

### Choosing between estimators

> *"My panel has N=45 countries, T=30 years, and I suspect cross-sectional dependence.
> The CIPS test rejects unit roots for all variables. Which dccelib estimator should I
> use and why?"*

### Debugging a GAUSS error

Paste the context block plus the exact error message and the relevant code block:

> *"I get 'error G0071: Type mismatch' on the line mgO.loglik = ll_tot. Here is my
> procedure: [paste code]. What is wrong?"*

### Extending the library

> *"I want to write a GAUSS procedure that takes an mgOut struct and returns a
> LaTeX-formatted table showing per-group coefficients for the first regressor. Follow
> dccelib conventions."*

---

## Providing More Context

For complex questions, attach `docs/api_reference.md` from the repository as an additional
file or paste its contents into the chat. Every LLM can work from it directly.

For debugging src files, paste the specific `.src` file section that is failing along with
the error message.

---

## Tool-Specific Notes

### Claude (claude.ai / Claude Code)

Claude works well with long context. You can safely paste the entire `api_reference.md`
(~800 lines) in addition to the context block. Ask it to reason step by step for debugging
tasks.

Recommended approach for larger tasks: use **Claude Code** (the CLI tool) in the
repository directory. The `CLAUDE.md` file in the repo root is automatically loaded,
giving Claude full project context without any pasting required.

Effective prompt pattern:
```
[paste context block]

Task: [your question or task]
Constraints: use only dccelib public procedures, GAUSS 26+ syntax.
```

### ChatGPT (GPT-4o / o1)

ChatGPT performs well on GAUSS econometrics when given the context block. For code
generation, explicitly request that it follow the GAUSS conventions listed in the context
block, since it tends to default to MATLAB or R syntax for array operations.

If responses still use R or MATLAB patterns, add:
> *"Write strictly valid GAUSS 26+ code. Use `$|` for vertical string concatenation,
> `.==` for element-wise comparison, and `struct` declarations at the top of each
> procedure."*

For long sessions, re-paste the context block at the start of a new conversation — ChatGPT
does not retain memory across sessions by default.

### Gemini (gemini.google.com / Gemini Advanced)

Gemini Advanced handles long context well. Paste the context block once at the session
start. For code generation tasks, be explicit that you want GAUSS code, not Python or R:

> *"Write this in GAUSS 26+ language for the dccelib library."*

Gemini is particularly useful for explaining econometric methodology (the theory behind
CCE, CIPS, Westerlund tests) even without the context block.

### GitHub Copilot

Copilot works best when you have the dccelib source files open in your editor. With the
`.src` and `.sdf` files visible, Copilot can offer inline completions for struct field
access and procedure calls.

For chat mode, paste the context block and describe what you want to add or fix. Copilot
is most effective for:
- Autocompleting `mgO.` field accesses when editing code
- Writing boilerplate struct initialization blocks
- Adapting existing example files to new datasets

Copilot does **not** have dccelib in its training data, so always provide the context
block for non-trivial questions.

---

## Prompt Templates

Copy-paste ready prompts for common tasks. Fill in the `[ ]` placeholders.

**Model estimation:**
```
[paste context block]

Write GAUSS code to estimate [cce_mg / dcce_mg] with:
- Dependent variable: [var name]
- Regressors: [var names]
- Extra CSA variables: [var names, or "none"]
- y_lags: [0 or 1]
- Data file: [filename and format]
Print the results table and save coefficients to local variables.
```

**Bias correction + bootstrap:**
```
[paste context block]

I have a dcce_mg result stored in dcceO with y_lags=1, cr_lags=3.
Write GAUSS code to:
1. Apply HPJ bias correction with hpj()
2. Compute wild bootstrap SEs with mgBootstrap() using B=499
3. Print both sets of results side by side
```

**Full diagnostic workflow:**
```
[paste context block]

Write a GAUSS script that runs a complete pre-estimation diagnostic workflow:
1. CIPS unit root test on [list variables] (p=1, demean=1)
2. Slope homogeneity test on a cce_mg baseline
3. Westerlund cointegration test on [y variable] and [x variables]
4. CD test from the baseline estimate
Print a clear header before each test result.
Data is in reg_data with columns [id, year, y, x1, x2].
```

**Debugging:**
```
[paste context block]

I am getting this GAUSS error:
[paste full error message and line number]

Here is the relevant code:
[paste the procedure or code block]

What is causing this error and how do I fix it?
```

**Comparing specifications:**
```
[paste context block]

Write GAUSS code to estimate three models on the same data:
1. mg() baseline
2. cce_mg() with automatic CSA lags
3. dcce_mg() with y_lags=1, cr_lags=3

Then export all three side-by-side to a LaTeX table using
mgOutToLatexMulti() with the note: "[your note text]".
```

---

## Common Mistakes LLMs Make

These are patterns that LLMs often get wrong without the context block. Mention them
explicitly if you see errors.

| Mistake | Correct dccelib / GAUSS usage |
|---------|-------------------------------|
| Uses `any(vec)` | No `any()` in GAUSS — use `sumc(vec) > 0` or explicit comparison |
| Uses `==` for string comparison | Use `$==` for strings |
| Concatenates strings with `+` | Use `$+` for string concatenation |
| Declares `local` after executable code | All `local` declarations must come first in a procedure |
| Reads `hpjO.b_stats` for HPJ estimates | Use `hpjO.b_stats_hpj` — `b_stats` holds `[min, mean, max, sd]` |
| Uses `scalar` field type in structs | Use `matrix` type for all numeric fields in GAUSS 26 |
| Forgets `packr()` before calling estimators | Always `data = packr(data)` to remove missing rows |
| Wrong column order in data | dccelib requires `[group, time, y, x...]` — or set `ctl.y_var` and `ctl.x_vars` |
| Uses `\` path separators in GAUSS strings | Forward slashes only: `"C:/data/myfile.dta"` |

---

## Key References

For deeper econometric questions about the methods implemented in dccelib:

- **Pesaran (2006)** — CCE estimator: *Econometrica* 74(4), 967–1012
- **Chudik & Pesaran (2015)** — Dynamic CCE: *Journal of Econometrics* 188(2), 393–420
- **Pesaran (2007)** — CIPS unit root test: *Journal of Econometrics* 141(2), 654–684
- **Pesaran & Yamagata (2008)** — Slope homogeneity: *Journal of Econometrics* 142(1), 50–93
- **Westerlund (2007)** — Panel cointegration: *Oxford Bulletin of Economics and Statistics* 69(6), 709–748
- **Dhaene & Jochmans (2015)** — HPJ bias correction: *Review of Economic Studies* 82(3), 991–1030
- **GAUSS language reference for LLMs**: https://github.com/aptech/gauss-llm-reference
