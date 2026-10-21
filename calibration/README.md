# Calibration design tables

Two waves of history matching (Appendix Text S-3). Each folder holds the design (one row per
parameter combination) and the per-design-point summary outputs read by `analysis/finalize_trend2.py`,
`nroy_marginals_trend2.py` and `nroy_structure_trend2.py`. Raw simulation outputs are not distributed.

| Folder | COMPS experiment | Design | Files |
|---|---|---|---|
| `wave1/` | `af707f22-16bc-f111-92f1-000d3af5294c` | 200 Latin-hypercube points x 10 replicates over the prior | `wave1_design_trend2.csv`, `design_outputs_trend.csv`, `bestfit_trend.json` |
| `wave2/` | `496ded25-37bc-f111-92f1-000d3af5294c` | 100 points inside the wave-1 not-ruled-out region x 25 replicates | `wave2_design_trend2.csv` (+ report), `design_outputs_trend.csv`, `bestfit_trend.json` |

`design_outputs_trend.csv` columns: the six burden targets (incidence 2000/2010/2015/2020, prevalence 2020,
mortality 2020), the held-out incidence years 2005 and 2023, the seven care-cascade benchmarks and the
2000-2020 incidence trend. The production parameter set is wave-2 design 4.
