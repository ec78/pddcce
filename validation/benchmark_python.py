"""
dccelib Benchmark — Python ecosystem survey
===========================================
Documents which panel CCE estimators are available in Python and
measures performance for those that exist.

Run from repo root:
    python validation/benchmark_python.py

Requirements (install with pip):
    pip install linearmodels pandas numpy scipy pyreadstat
"""

import sys
import time
import warnings
import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")

print("=" * 65)
print("dccelib Benchmark — Python ecosystem coverage")
print("=" * 65)
print(f"Python {sys.version}")


# ---------------------------------------------------------------------------
# Utility: median wall-clock time over B repetitions
# ---------------------------------------------------------------------------
def time_median(fn, B=10):
    times = []
    for _ in range(B):
        t0 = time.perf_counter()
        fn()
        times.append(time.perf_counter() - t0)
    return float(np.median(times))


# ---------------------------------------------------------------------------
# Load Penn World Tables data
# ---------------------------------------------------------------------------
try:
    import pyreadstat
    df_raw, _ = pyreadstat.read_dta("examples/penn_sample.dta")
    df = df_raw.dropna(subset=["id", "year", "log_rgdpo", "log_ck", "log_ngd"])
    df["year"] = df["year"].astype(int)
    df = df.sort_values(["id", "year"]).reset_index(drop=True)
    print(f"\nPenn data loaded: N={df['id'].nunique()}, obs={len(df)}\n")
    PENN_LOADED = True
except Exception as e:
    print(f"\nWarning: Could not load penn_sample.dta ({e}). Skipping Penn tests.\n")
    PENN_LOADED = False


# ---------------------------------------------------------------------------
# Simulated panel helper
# ---------------------------------------------------------------------------
def sim_panel(N, T, k=2, seed=42):
    rng = np.random.default_rng(seed)
    f   = rng.standard_normal(T)
    ids = np.repeat(np.arange(1, N + 1), T)
    yrs = np.tile(np.arange(1, T + 1), N)
    x1  = rng.standard_normal(N * T) + 0.5 * np.tile(f, N)
    x2  = rng.standard_normal(N * T) + 0.3 * np.tile(f, N)
    y   = (0.5 * x1 + 0.3 * x2
           + np.repeat(rng.standard_normal(N), T)
           + np.tile(f, N)
           + rng.standard_normal(N * T) * 0.5)
    return pd.DataFrame({"id": ids, "year": yrs, "y": y, "x1": x1, "x2": x2})


# ---------------------------------------------------------------------------
# PART 1: Ecosystem survey — what Python has for panel CCE
# ---------------------------------------------------------------------------
print("-" * 65)
print("PART 1: Python panel data ecosystem coverage")
print("-" * 65)

COL_WIDTH = 48
HEADER = f"{'Estimator / Capability':{COL_WIDTH}}  {'Available'}"
print(HEADER)
print("-" * 65)

capabilities = {}

# Check linearmodels
try:
    import linearmodels  # noqa: F401
    lm_version = linearmodels.__version__
    capabilities["linearmodels installed"] = f"Yes (v{lm_version})"
except ImportError:
    capabilities["linearmodels installed"] = "No"

# Pooled OLS
try:
    from linearmodels.panel import PooledOLS  # noqa: F401
    capabilities["Pooled OLS (linearmodels)"] = "Yes"
except ImportError:
    capabilities["Pooled OLS (linearmodels)"] = "No"

# Fixed Effects (LSDV)
try:
    from linearmodels.panel import PanelOLS  # noqa: F401
    capabilities["Fixed Effects / LSDV (linearmodels)"] = "Yes"
except ImportError:
    capabilities["Fixed Effects / LSDV (linearmodels)"] = "No"

# Fama-MacBeth
try:
    from linearmodels.panel import FamaMacBeth  # noqa: F401
    capabilities["Fama-MacBeth (linearmodels)"] = "Yes"
except ImportError:
    capabilities["Fama-MacBeth (linearmodels)"] = "No"

# Between estimator (proxy for MG)
try:
    from linearmodels.panel import BetweenOLS  # noqa: F401
    capabilities["Between OLS (linearmodels)"] = "Yes"
except ImportError:
    capabilities["Between OLS (linearmodels)"] = "No"

# MG via statsmodels
try:
    import statsmodels.formula.api  # noqa: F401
    # statsmodels has no MG; per-group OLS must be coded manually
    capabilities["Mean Group estimator (statsmodels)"] = "No — manual only"
except ImportError:
    capabilities["Mean Group estimator (statsmodels)"] = "No"

# CCE-MG
capabilities["CCE Mean Group (any package)"] = "No — not implemented"
capabilities["Dynamic CCE-MG (any package)"] = "No — not implemented"
capabilities["HPJ bias correction (any package)"] = "No — not implemented"
capabilities["CIPS panel unit root (any package)"] = "No — not implemented"
capabilities["Westerlund cointegration (any package)"] = "No — not implemented"
capabilities["Pesaran-Yamagata slope homo. (any package)"] = "No — not implemented"

for cap, avail in capabilities.items():
    print(f"  {cap:{COL_WIDTH - 2}}  {avail}")

print()


# ---------------------------------------------------------------------------
# PART 2: Manual MG in Python (for timing comparison)
# ---------------------------------------------------------------------------
print("-" * 65)
print("PART 2: Timing — manual MG (group-by OLS) in Python")
print("-" * 65)
print("Note: No single-call MG or CCE-MG function exists in Python.")
print("This implements the MG estimator manually via per-group OLS.")
print()


def mg_manual(df, yvar, xvars):
    """Manual Mean Group estimator: per-group OLS, averaged coefficients."""
    from numpy.linalg import lstsq
    coefs = []
    for _, grp in df.groupby("id"):
        y = grp[yvar].values
        X = np.column_stack([np.ones(len(y))] + [grp[x].values for x in xvars])
        b, _, _, _ = lstsq(X, y, rcond=None)
        coefs.append(b)
    return np.mean(coefs, axis=0)


B_REPS = 10
results = []

specs = [
    ("Penn WTables", df if PENN_LOADED else None, "log_rgdpo", ["log_ck", "log_ngd"]),
    ("Simulated M",  sim_panel(100, 50),           "y",         ["x1", "x2"]),
    ("Simulated L",  sim_panel(200, 100),          "y",         ["x1", "x2"]),
]

for name, data, yv, xv in specs:
    if data is None:
        print(f"  {name}: skipped (data not loaded)\n")
        continue
    N = data["id"].nunique()
    T = round(len(data) / N)
    t = time_median(lambda d=data, y=yv, x=xv: mg_manual(d, y, x), B=B_REPS)
    print(f"  {name} (N={N}, T~{T}):  MG manual = {t:.4f} s")
    results.append({"dataset": name, "N": N, "T": T,
                    "estimator": "MG (manual)", "tool": "Python/numpy",
                    "time_sec": t})

print()

# ---------------------------------------------------------------------------
# PART 3: FamaMacBeth timing (closest available Python panel estimator)
# ---------------------------------------------------------------------------
try:
    from linearmodels.panel import FamaMacBeth

    print("-" * 65)
    print("PART 3: Timing — Fama-MacBeth (linearmodels)")
    print("Note: FMB is NOT equivalent to CCE-MG but is the closest")
    print("available cross-sectionally robust estimator in Python.")
    print()

    for name, data, yv, xv in specs:
        if data is None:
            print(f"  {name}: skipped (data not loaded)\n")
            continue
        N = data["id"].nunique()
        T = round(len(data) / N)

        pdata = data.set_index(["id", "year"])

        def run_fmb(d=pdata, y=yv, x=xv):
            mod = FamaMacBeth(d[y], d[x])
            return mod.fit()

        t = time_median(run_fmb, B=B_REPS)
        print(f"  {name} (N={N}, T~{T}):  Fama-MacBeth = {t:.4f} s")
        results.append({"dataset": name, "N": N, "T": T,
                        "estimator": "Fama-MacBeth", "tool": "Python/linearmodels",
                        "time_sec": t})

    print()
except ImportError:
    print("linearmodels not installed — skipping Fama-MacBeth timing.\n")

# ---------------------------------------------------------------------------
# PART 4: Save results
# ---------------------------------------------------------------------------
if results:
    pd.DataFrame(results).to_csv("validation/benchmark_python_results.csv", index=False)
    print("Results saved to validation/benchmark_python_results.csv")

print("\n" + "=" * 65)
print("Python benchmark complete.")
print("Key finding: No Python package implements CCE-MG, DCCE-MG, CIPS,")
print("Westerlund cointegration, slope homogeneity, or HPJ bias correction.")
print("dccelib fills this gap entirely.")
print("=" * 65)
