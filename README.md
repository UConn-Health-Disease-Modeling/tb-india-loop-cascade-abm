# Care-cascade failures as a determinant of infectious burden in tuberculosis: an agent-based modelling study in India

Code, derived data and figures for the manuscript of the same title (Jo et al., 2026, in preparation).
The model is an EMOD-TBHIV individual-based simulation of TB transmission and care in India, calibrated to
national incidence, prevalence and mortality, in which every active TB episode is followed through its
diagnostic and treatment events and assigned to a care-cascade "loop group".

Three questions are answered here:

1. **What is in the repository** and what is not (raw simulation outputs are not; see below).
2. **How each figure and table in the manuscript was produced** (the output map).
3. **How to reproduce** the results at two levels: from the bundled derived data (minutes, no EMOD needed)
   or from scratch (requires EMOD-TBHIV and the IDM COMPS cluster).

## Repository layout

```
.
├── 01_India_TB_v23_wave1.ipynb      calibration wave 1  (200 LHS design points x 10 replicates)
├── 01_India_TB_v24_wave2.ipynb      calibration wave 2  (100 NROY-restricted points x 25 replicates)
├── 01_India_TB_v28_txgated.ipynb    production run: 8 scenarios x 100 replicates at the calibrated point,
│                                     and (through helpers/run_experiment_B.py) the parameter-propagation run
├── 01_India_TB_v30_loopsweep.ipynb  residual loop-closure sensitivity (Appendix Figure S-2)
├── manifest.py, params.py, set_config.py, add_tb_drug_type.py   model definition used by the notebooks
├── python_scripts/, EP4/            EMOD pre/post-processing hooks uploaded with each simulation
├── Assets/, data/, drug.json, twonodes.csv, dtk_centos.id        model inputs (demographic overlay,
│                                     India all-cause mortality series, prevalence-survey reference tables)
├── helpers/                         calibration, extraction and figure scripts + the derived data they read
├── <experiment-id>/calibration/     per-design-point summaries of the two calibration experiments
├── figures/                         the manuscript figures, named by manuscript number
├── docs/                            calibration notes
├── requirements-simulation.txt      package versions of the simulation environment (Python 3.10)
└── requirements-analysis.txt        minimal environment for the scripts in helpers/
```

## Output map

Manuscript item → script → data it reads → experiment. Scripts live in `helpers/`. "Derived" files are
bundled in `helpers/`, so every row marked *bundled* reproduces without the raw outputs.

| Manuscript item | Script | Reads | Writes | Raw outputs needed? |
|---|---|---|---|---|
| Figure 1A/1B, Table 2 (cohort rows) | `fig1_rebuild_txc.py` → `fig1_plot_txc.py` | ReportEventRecorder of d892bdf7 → `fig1_values_txc.json`, `baseline_episodes_year2_txc.csv.gz` | `Fig1A_cascade_rates.png`, `Fig1B_patient_trajectories.png` | rebuild: yes; plot: bundled |
| Figure 2A | `fig2_rebuild_txc.py` | InsetChart of d892bdf7 → `production_burden_txc.csv.gz` | `Fig2_calibration_trend.png` | bundled cache |
| Figure 2B (panel B2 = tests per *tested* episode, the benchmark denominator), Table 2 delay phases | `fig2b_and_table2_txc.py` | `baseline_episodes_year2_txc.csv.gz` | `Fig2B_cascade_validation.png`, `table2_values_txc.json` | bundled |
| Figure 3 | `fig3_rebuild_txc.py --plot-only` | `baseline_episodes_year2_txc.csv.gz`, `baseline_popyears_2015_2024.json` | `Fig3_loopgroup_burden.png`, `fig3_values_txc.json` | bundled (`--plot-only`); the denominator needs InsetCharts |
| Figure 4 (three-lens heatmap) | `scenarios_rebuild_txc.py` → `fig5_rebuild_txc.py` | ReportEventRecorder of all 8 scenarios → `scenario_episodes_txc.csv.gz` | `Fig5_triple_lens.png`, `fig5_values_txc.csv` | rebuild: yes; plot: bundled |
| Figure 5 | `extract_burden_incl.py`, `extract_propagation.py` → `fig5_endtb_panels_incl.py` | InsetCharts of d892bdf7 and f8dfe706 → `burden_by_rep_incl.csv.gz`, `burden_propagation.csv.gz` | `Fig5_endtb_panels.png`, `fig5_propagation_table.csv` | extract: yes; plot: bundled |
| Appendix Table S-8, propagation text (S-3.8) | `propagation_contrasts.py` | `burden_propagation.csv.gz` | `propagation_contrasts.csv` | bundled |
| Appendix Figure S-2 | `fig4_incidence_ci_txc.py` → `figS_loop_sensitivity.py` | InsetCharts of d892bdf7 and 941e6423 → `incidence_by_rep_txc.csv.gz`, `loopsweep_incidence_txc.csv.gz` | `FigS_loop_sensitivity.png` | extract: yes; plot: bundled |
| Appendix Tables S-1 to S-6 | `targets_trend.py`, `calibration_wave1_lhs_trend2.py`, `wave2_design_trend2.py`, `compute_targets_trend.py`, `finalize_trend2.py` | the two `<experiment-id>/calibration/design_outputs_trend.csv` (bundled) | `finalize_trend2_report.txt`, `production_params_trend2.json` | design tables bundled |
| Appendix Table S-7 (NROY marginals) | `nroy_marginals_trend2.py` | same | report to stdout | design tables bundled |
| Appendix S-3.8 text (NROY structure: parameter rank correlations, conditional spans, sensitivity of the cascade benchmarks across burden-plausible design points) | `nroy_structure_trend2.py` | same | `nroy_structure_trend2.txt`, `nroy_corr_trend2.csv` | design tables bundled |
| Propagation run (20 parameter sets x 8 scenarios x 10 replicates) | `run_experiment_B.py` (called from the v28 notebook), `nroy_subsamples_b.json` | — | experiment f8dfe706 | COMPS |
| Historical cascade trends and intervention switches | `campaign_trend.py` | — | campaign events used by every notebook | — |

Experiments (IDM COMPS, CALCULON): calibration wave 1 `af707f22-16bc-f111-92f1-000d3af5294c`, wave 2
`496ded25-37bc-f111-92f1-000d3af5294c`, production `d892bdf7-18bd-f111-92f1-000d3af5294c`, parameter
propagation `f8dfe706-d7c0-f111-92f1-000d3af5294c`, loop-closure sweep `941e6423`. Raw outputs are
8–45 GB per experiment and are not distributed; the derived files above are what the figures are drawn from.

The trajectory results use episodes activating in calendar years 2015–2024 (`COHORT` in
`fig1_rebuild_txc.py`), with the loop group assigned per episode; `cohort_2015_2024_numbers.json` collects
the headline numbers. Appendix Figure S-1 is a schematic drawn by hand.

## Reproducing from the bundled derived data

```bash
pip install -r requirements-analysis.txt
cd helpers
python fig1_plot_txc.py
python fig2_rebuild_txc.py
python fig2b_and_table2_txc.py
python fig3_rebuild_txc.py --plot-only
python fig5_rebuild_txc.py           # manuscript Figure 4
python fig5_endtb_panels_incl.py     # manuscript Figure 5
python figS_loop_sensitivity.py
python propagation_contrasts.py
python nroy_marginals_trend2.py        # Appendix Table S-7
python nroy_structure_trend2.py        # Appendix S-3.8 numbers
```

Figures are written to `figures/Figure_TxCure/` next to the bundled copies in `figures/`.

## Reproducing from scratch

1. Create the simulation environment (`requirements-simulation.txt`; Python 3.10) and obtain the EMOD-TBHIV
   binary and schema through `emodpy` bootstrap into `download/` (see `manifest.py`). Simulations were run
   on the IDM COMPS platform; `idmtools` and a COMPS account are required.
2. Run `01_India_TB_v23_wave1.ipynb` and `01_India_TB_v24_wave2.ipynb` (`SCENARIO = 'wave1_cal'` /
   `'wave2_cal'`), then `helpers/compute_targets_trend.py` on each experiment and
   `helpers/finalize_trend2.py` to obtain the calibrated point (wave-2 design 4) and the NROY.
3. Run `01_India_TB_v28_txgated.ipynb` with `SCENARIO = 'combined'` for the production experiment, and the
   cell described in `helpers/run_experiment_B.py` for the propagation experiment.
4. Run the extraction scripts in the output map (`fig1_rebuild_txc.py`, `scenarios_rebuild_txc.py`,
   `extract_burden_incl.py`, `extract_propagation.py`, `fig4_incidence_ci_txc.py`), then the plot scripts.

Each script's docstring records the experiment it reads and the analysis decisions behind it.

## Software

Simulations: EMOD-TBHIV (DTK v2.21) through emodpy 1.35.1, emod-api 1.35.1, emodpy-tbhiv 1.3.1 and
idmtools 3.0.8 under Python 3.10.20. Analysis: numpy 2.2.6, pandas 2.3.3, matplotlib 3.10.8,
scikit-learn 1.7.2 (Gaussian-process emulators for history matching).

## Licence

CC0 1.0 (see `LICENSE`), following the IDM `TB_India_DTK` repository this work builds on.
