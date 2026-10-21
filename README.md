# Care-cascade failures as a determinant of infectious burden in tuberculosis: an agent-based modelling study in India

Code, derived data and figures for the manuscript of the same title (Jo et al., 2026, under review).

The model is an EMOD-TBHIV individual-based simulation of TB transmission and care in India. It is calibrated
to national incidence, prevalence and mortality, and every active TB episode is followed through its diagnostic
and treatment events and assigned to a care-cascade group (no loop, diagnostic loop only, any treatment
failure). The manuscript asks how much infectious burden is created by repeated failures of care, and what
interventions at different stages of the cascade would change.

## What is in the repository

| Folder | Contents | Who needs it |
|---|---|---|
| `figures/` | The manuscript figures, named by manuscript number | Everyone |
| `analysis/` | The scripts that produce every figure and table, and in `analysis/derived/` the data they read | Anyone checking how a number was obtained |
| `calibration/` | Design tables and per-design-point outputs of the two history-matching waves (Appendix Text S-3) | Anyone checking the calibration |
| `simulation/` | The four notebooks that ran the simulations on IDM COMPS, with the model definition and EMOD input files | Anyone re-running the model |
| `docs/` | Development notes on the calibration, kept for the record | Optional |

Raw simulation outputs (8–45 GB per experiment) are not distributed. Everything in `analysis/derived/` was
extracted from them, and every figure and table can be rebuilt from those files in a few minutes without EMOD.

## Where each figure and table comes from

Run the scripts from inside `analysis/`. Outputs go to `figures/Figure_TxCure/`; the copies in `figures/`
are the versions in the manuscript.

| Manuscript item | Script | Reads (in `analysis/derived/` unless noted) |
|---|---|---|
| Figure 1A, 1B; Table 2 cohort rows | `fig1_plot_txc.py` | `fig1_values_txc.json` |
| Figure 2A (calibration trajectories) | `fig2_rebuild_txc.py` | `production_burden_txc.csv.gz` |
| Figure 2B (cascade benchmarks); Table 2 delay phases | `fig2b_and_table2_txc.py` | `baseline_episodes_year2_txc.csv.gz` → `table2_values_txc.json` |
| Figure 3 (loop-group burden) | `fig3_rebuild_txc.py --plot-only` | `baseline_episodes_year2_txc.csv.gz`, `baseline_popyears_2015_2024.json` |
| Figure 4 (three-lens intervention impact) | `fig5_rebuild_txc.py` | `scenario_episodes_txc.csv.gz` |
| Figure 5 (End TB progress, with parameter propagation) | `fig5_endtb_panels_incl.py` | `burden_by_rep_incl.csv.gz`, `burden_propagation.csv.gz` |
| Appendix Figure S-2 (loop-closure sensitivity) | `figS_loop_sensitivity.py` | `incidence_by_rep_txc.csv.gz`, `loopsweep_incidence_txc.csv.gz` |
| Appendix Tables S-1 to S-6 (targets, design, fit) | `finalize_trend2.py` | `calibration/wave1/`, `calibration/wave2/` → `finalize_trend2_report.txt` |
| Appendix Table S-7 (NROY marginals) | `nroy_marginals_trend2.py` | same → `nroy_marginals_trend2.csv` |
| Appendix S-3.8 (structure of the NROY region) | `nroy_structure_trend2.py` | same → `nroy_structure_trend2.txt`, `nroy_corr_trend2.csv` |
| Appendix Table S-8 (across-parameter-set intervals) | `propagation_contrasts.py` | `burden_propagation.csv.gz` → `propagation_contrasts.csv` |

Headline numbers for the 2015–2024 episode cohort are collected in `analysis/derived/cohort_2015_2024_numbers.json`.
Appendix Figure S-1 is a hand-drawn schematic.

Notes on Figure 2: panel labels are A1–A3 and B1–B3; the y-axes of A2 and A3 start at zero; panel B2 reports
tests per *tested* episode, which is the denominator of the published benchmark.

## Reproducing the figures and tables (minutes, no EMOD)

```bash
pip install -r analysis/requirements-analysis.txt
cd analysis
python fig1_plot_txc.py
python fig2_rebuild_txc.py
python fig2b_and_table2_txc.py
python fig3_rebuild_txc.py --plot-only
python fig5_rebuild_txc.py              # manuscript Figure 4
python fig5_endtb_panels_incl.py        # manuscript Figure 5
python figS_loop_sensitivity.py
python propagation_contrasts.py
python finalize_trend2.py               # Appendix Tables S-1 to S-6
python nroy_marginals_trend2.py         # Appendix Table S-7
python nroy_structure_trend2.py         # Appendix S-3.8
```

Each script's docstring states the experiment it reads and the analysis decisions behind it.

## Reproducing from scratch (EMOD-TBHIV on IDM COMPS)

1. Create the simulation environment (`simulation/requirements-simulation.txt`, Python 3.10) and obtain the
   EMOD-TBHIV binary and schema through the `emodpy` bootstrap described in `simulation/manifest.py`.
   A COMPS account and `idmtools` are required.
2. From `simulation/`, run `01_India_TB_v23_wave1.ipynb` and `01_India_TB_v24_wave2.ipynb`
   (`SCENARIO = 'wave1_cal'` / `'wave2_cal'`), then `analysis/compute_targets_trend.py` on each experiment and
   `analysis/finalize_trend2.py` to obtain the calibrated point (wave-2 design 4) and the not-ruled-out region.
3. Run `01_India_TB_v28_txgated.ipynb` with `SCENARIO = 'combined'` for the production experiment, and the cell
   described in `analysis/run_experiment_B.py` for the parameter-propagation experiment.
4. Run the extraction scripts (`fig1_rebuild_txc.py`, `scenarios_rebuild_txc.py`, `extract_burden_incl.py`,
   `extract_propagation.py`, `fig4_incidence_ci_txc.py`), then the plot scripts above.

| Experiment (IDM COMPS, CALCULON) | ID |
|---|---|
| Calibration wave 1 (200 design points × 10 replicates) | `af707f22-16bc-f111-92f1-000d3af5294c` |
| Calibration wave 2 (100 × 25) | `496ded25-37bc-f111-92f1-000d3af5294c` |
| Production (8 scenarios × 100 replicates) | `d892bdf7-18bd-f111-92f1-000d3af5294c` |
| Parameter propagation (20 parameter sets × 8 scenarios × 10 replicates) | `f8dfe706-d7c0-f111-92f1-000d3af5294c` |
| Loop-closure sweep (Appendix Figure S-2) | `941e6423` |

## Software

Simulations: EMOD-TBHIV (DTK v2.21) through emodpy 1.35.1, emod-api 1.35.1, emodpy-tbhiv 1.3.1 and
idmtools 3.0.8 under Python 3.10.20. Analysis: numpy 2.2.6, pandas 2.3.3, matplotlib 3.10.8,
scikit-learn 1.7.2 (Gaussian-process emulators for history matching).

## Licence

CC0 1.0 (see `LICENSE`), following the IDM `TB_India_DTK` repository this work builds on.
