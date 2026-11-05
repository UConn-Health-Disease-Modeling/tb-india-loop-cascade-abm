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
| `simulation/` | The notebooks that ran the simulations on IDM COMPS, numbered in run order, with the model definition and EMOD input files | Anyone re-running the model |
| `docs/` | Development notes on the calibration, kept for the record | Optional |

Raw simulation outputs (8–45 GB per experiment) are not distributed. Everything in `analysis/derived/` was
extracted from them, and every figure and table can be rebuilt from those files in a few minutes without EMOD.
With a COMPS account, `python simulation/download_outputs.py` downloads the raw outputs the extraction scripts
read and puts each experiment where they look for it; after re-running the extraction scripts,
`python analysis/check_derived.py` confirms that the re-extracted tables match the committed ones.

## How the pieces fit together

```
simulation/1_calibration_wave1.ipynb ──┐
simulation/2_calibration_wave2.ipynb ──┴─► analysis/calib1 … calib6 ─► calibration/wave1, wave2
                                                                       analysis/derived/production_parameters.json
                                                                                 │
                                       ┌─────────────────────────────────────────┘ (production parameter set)
                                       ▼
simulation/3_production_scenarios.ipynb   (+ 3b_parameter_propagation_run.py)
simulation/4_loop_closure_sweep.ipynb
        │  raw outputs on COMPS (not distributed)
        ▼
analysis/extract1 … extract6 ─► analysis/derived/ ─► analysis/figure*, tableS8 ─► figures/
```

Script names in `analysis/` say which stage they belong to:

| Prefix | Stage | Needs raw simulation outputs? |
|---|---|---|
| `calib1_` … `calib6_` | History-matching calibration, in run order | Only `calib2_`; the rest read `calibration/` |
| `extract1_` … `extract6_` | Reads the raw outputs of an experiment and writes a table to `analysis/derived/` | Yes |
| `figure…_`, `table…_` | Builds one manuscript figure or table, named by its manuscript number | No (reads `analysis/derived/`) |
| `calibration_targets.py`, `campaign_schedules.py` | Shared definitions: the calibration targets, and the time-varying cascade schedules used to build the EMOD campaign | — |

## Where each figure and table comes from

Run the scripts from inside `analysis/`. Outputs go to `figures/reproduced/` under the same file names as the
manuscript copies in `figures/`, so the two can be compared directly.

| Manuscript item | Script | Reads (in `analysis/derived/` unless noted) |
|---|---|---|
| Figure 1A, 1B; Table 2 cohort rows | `figure1_cascade_events_and_trajectories.py` | `figure1_values.json` |
| Figure 2A (calibration trajectories) | `figure2a_calibration.py` | `production_baseline_burden.csv.gz` |
| Figure 2B (cascade benchmarks); Table 2 delay phases | `figure2b_table2_cascade_benchmarks.py` | `baseline_episodes.csv.gz` → `table2_values.json` |
| Figure 3 (loop-group burden) | `figure3_loop_group_burden.py --plot-only` | `baseline_episodes.csv.gz`, `baseline_population_years_2015_2024.json` |
| Figure 4 (three-lens intervention impact) | `figure4_three_lens_intervention_impact.py` | `scenario_episodes.csv.gz` → `figure4_values.csv` |
| Figure 5 (End TB progress, with parameter propagation) | `figure5_end_tb_progress.py` | `burden_by_replicate_both_prevalence.csv.gz`, `burden_by_parameter_set.csv.gz` |
| Appendix Figure S-2 (loop-closure sensitivity) | `figureS2_loop_closure_sensitivity.py` | `burden_by_replicate.csv.gz`, `loop_closure_sweep_incidence.csv.gz` |
| Appendix Tables S-1 to S-6 (targets, design, fit) | `calib4_final_fit.py` | `calibration/wave1/`, `calibration/wave2/` → `calibration_final_fit_report.txt` |
| Appendix Table S-7 (NROY marginals) | `calib5_nroy_marginals.py` | same → `nroy_marginals.csv` |
| Appendix S-3.8 (structure of the NROY region) | `calib6_nroy_structure.py` | same → `nroy_structure_report.txt`, `nroy_rank_correlations.csv` |
| Appendix Table S-8 (across-parameter-set intervals) | `tableS8_propagation_contrasts.py` | `burden_by_parameter_set.csv.gz` → `propagation_contrasts.csv` |

Headline numbers for the 2015–2024 episode cohort are collected in `analysis/derived/cohort_2015_2024_numbers.json`.
Appendix Figure S-1 is a hand-drawn schematic.

Notes on Figure 2: panel labels are A1–A3 and B1–B3; the y-axes of A2 and A3 start at zero; panel B2 reports
tests per *tested* episode, which is the denominator of the published benchmark.

## Reproducing the figures and tables (minutes, no EMOD)

```bash
pip install -r analysis/requirements-analysis.txt
cd analysis
python figure1_cascade_events_and_trajectories.py
python figure2a_calibration.py
python figure2b_table2_cascade_benchmarks.py
python figure3_loop_group_burden.py --plot-only
python figure4_three_lens_intervention_impact.py
python figure5_end_tb_progress.py
python figureS2_loop_closure_sensitivity.py
python tableS8_propagation_contrasts.py
python calib4_final_fit.py          # Appendix Tables S-1 to S-6
python calib5_nroy_marginals.py     # Appendix Table S-7
python calib6_nroy_structure.py     # Appendix S-3.8
```

Each script's docstring states the experiment it reads and the analysis decisions behind it.

## Reproducing from scratch (EMOD-TBHIV on IDM COMPS)

1. Create the simulation environment (`simulation/requirements-simulation.txt`, Python 3.10) and obtain the
   EMOD-TBHIV binary and schema through the `emodpy` bootstrap described in `simulation/manifest.py`.
   A COMPS account and `idmtools` are required.
2. Calibration. `analysis/calib1_wave1_design.py` writes the wave-1 design. From `simulation/`, run
   `1_calibration_wave1.ipynb` (`SCENARIO = 'wave1_cal'`), then `analysis/calib2_design_point_outputs.py` on
   that experiment and copy the `design_outputs.csv` and `best_fit.json` it writes into `calibration/wave1/`.
   `analysis/calib3_wave2_design.py` then writes the wave-2 design. Run `2_calibration_wave2.ipynb`
   (`SCENARIO = 'wave2_cal'`) and `calib2_design_point_outputs.py` again (copy into `calibration/wave2/`), then
   `analysis/calib4_final_fit.py` to obtain the calibrated point (wave-2 design 4) and the not-ruled-out region.
3. Run `3_production_scenarios.ipynb` with `SCENARIO = 'combined'` for the production experiment, and the cell
   described in `simulation/3b_parameter_propagation_run.py` for the parameter-propagation experiment.
   `4_loop_closure_sweep.ipynb` runs the loop-closure sweep behind Appendix Figure S-2.
4. Run the extraction scripts `extract1_baseline_episodes.py` to `extract6_loop_closure_sweep.py`
   (`figure2a_calibration.py` and `figure3_loop_group_burden.py` without `--plot-only` also extract what they
   need), then the figure and table scripts above.

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
