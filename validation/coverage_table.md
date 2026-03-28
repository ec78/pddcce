## Panel Estimator Coverage: dccelib vs R plm vs Python

| Estimator / Feature | GAUSS dccelib | R plm | Python linearmodels |
|---------------------|:-------------:|:-----:|:-------------------:|
| Mean Group (MG)                           | Yes           | Yes (`pmg`) | No                  |
| CCE Mean Group (CCE-MG)                   | Yes           | Yes (`pmg`) | No                  |
| Dynamic CCE-MG (DCCE-MG)                  | Yes           | No    | No                  |
| Pooled CCE with NW SE                     | Yes           | Yes (`pcce`) | No                  |
| PC-CCE-MG (PCA augmentation)              | Yes (GML req.) | No    | No                  |
| HPJ bias correction                       | Yes           | No    | No                  |
| Wild bootstrap SEs                        | Yes           | No    | No                  |
| CIPS panel unit root (Pesaran 2007)       | Yes           | No (CD.test only) | No                  |
| Pesaran-Yamagata slope homogeneity        | Yes           | No    | No                  |
| Westerlund cointegration test             | Yes           | No (in `westerlund` pkg) | No                  |
| CCE rank condition test                   | Yes           | No    | No                  |
| Long-run multipliers (delta method)       | Yes           | No    | No                  |
| Residual / coefficient plots              | Yes           | No    | No                  |
| LaTeX table export                        | Yes           | No    | No                  |
| I(1) extension (KPY 2011)                 | Yes           | No    | No                  |
| Two-way CCE (Bai 2009)                    | Yes           | No    | No                  |

_Table reflects package capabilities as of 2026. plm CCE-MG (`pmg` model="cmg") uses the Pesaran (2006) estimator without dynamic extension or auxiliary diagnostics._
