## dccelib Performance Benchmarks

Wall-clock time (seconds, median of 10 runs).

dataset | N | T | estimator | GAUSS/dccelib | Python/numpy | R/plm | GAUSS/R ratio
--- | --- | --- | --- | --- | --- | --- | ---
Penn WTables | 93 | 47 | CCE-MG | 0.014 | — | 0.015 | 0.93 
Penn WTables | 93 | 47 | CIPS | 0.0035 | — | — | — 
Penn WTables | 93 | 47 | DCCE-MG | 0.021 | — | — | — 
Penn WTables | 93 | 47 | MG | 0.0125 | — | 0.02 | 0.62 
Penn WTables | 93 | 47 | Westerlund | 0.008 | — | — | — 
Simulated L | 200 | 100 | CCE-MG | 0.093 | — | 0.045 | 2.07 
Simulated L | 200 | 100 | CIPS | 0.01 | — | — | — 
Simulated L | 200 | 100 | DCCE-MG | 0.1115 | — | — | — 
Simulated L | 200 | 100 | MG | 0.0965 | — | 0.05 | 1.93 
Simulated L | 200 | 100 | MG (manual) | — | 0.019 | — | — 
Simulated L | 200 | 100 | Westerlund | 0.0215 | — | — | — 
Simulated M | 100 | 50 | CCE-MG | 0.0285 | — | 0.02 | 1.43 
Simulated M | 100 | 50 | CIPS | 0.004 | — | — | — 
Simulated M | 100 | 50 | DCCE-MG | 0.042 | — | — | — 
Simulated M | 100 | 50 | MG | 0.019 | — | 0.02 | 0.95 
Simulated M | 100 | 50 | MG (manual) | — | 0.0096 | — | — 
Simulated M | 100 | 50 | Westerlund | 0.009 | — | — | — 

_GAUSS/R ratio < 1 indicates GAUSS dccelib is faster._
