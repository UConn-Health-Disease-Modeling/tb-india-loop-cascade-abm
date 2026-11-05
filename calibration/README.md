# Calibration design tables

Two waves of history matching (Appendix Text S-3). Each folder holds the design (one row per
parameter combination) and the per-design-point summary outputs read by `analysis/calib4_final_fit.py`,
`calib5_nroy_marginals.py` and `calib6_nroy_structure.py`. Raw simulation outputs are not distributed.

| Folder | COMPS experiment | Design | Design written by |
|---|---|---|---|
| `wave1/` | `af707f22-16bc-f111-92f1-000d3af5294c` | 200 Latin-hypercube points x 10 replicates over the prior | `design.csv`: `analysis/calib1_wave1_design.py` |
| `wave2/` | `496ded25-37bc-f111-92f1-000d3af5294c` | 100 points inside the wave-1 not-ruled-out region x 25 replicates | `design.csv`, `design_report.txt`: `analysis/calib3_wave2_design.py` |

In both folders, `design_outputs.csv` and `best_fit.json` are written by
`analysis/calib2_design_point_outputs.py` from that wave's simulations.

`design_outputs.csv` columns: the six burden targets (incidence 2000/2010/2015/2020, prevalence 2020,
mortality 2020), the held-out incidence years 2005 and 2023, the seven care-cascade benchmarks and the
2000-2020 incidence trend. The production parameter set is wave-2 design 4.
