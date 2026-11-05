# analysis/

Every script that turns simulation output into a number, table or figure in the manuscript. Run the scripts
from inside this folder. The top-level `README.md` maps each manuscript item to its script.

## Scripts, by stage

**Shared definitions** (imported, not run)

| Script | What it defines | Used by |
|---|---|---|
| `calibration_targets.py` | WHO incidence trajectory and the burden and cascade targets | `calib2`–`calib6`, `figure2a` |
| `campaign_schedules.py` | Time-varying care-seeking, treatment-initiation and infectivity schedules for the EMOD campaign | the `simulation/` notebooks |

**Calibration** (history matching, Appendix Text S-3), in run order

| Script | What it does | Writes |
|---|---|---|
| `calib1_wave1_design.py` | 200-point Latin hypercube over the 10 free parameters | `calibration/wave1/design.csv` |
| `calib2_design_point_outputs.py` | Target values for every design point of one calibration experiment (needs raw outputs) | `design_outputs.csv`, `best_fit.json` |
| `calib2_target_readers.py` | Readers for EMOD output files, imported by `calib2_design_point_outputs.py` | — |
| `calib3_wave2_design.py` | Wave-1 emulators and the wave-2 design inside the not-ruled-out region | `calibration/wave2/design.csv`, `design_report.txt` |
| `calib4_final_fit.py` | Pooled wave-1 + wave-2 fit; selects the production parameter set (Appendix Tables S-1 to S-6) | `derived/production_parameters.json`, `derived/calibration_final_fit_report.txt` |
| `calib5_nroy_marginals.py` | Width of the not-ruled-out region on each parameter (Appendix Table S-7) | `derived/nroy_marginals.csv` |
| `calib6_nroy_structure.py` | Trade-offs among parameters inside that region (Appendix S-3.8) | `derived/nroy_structure_report.txt`, `derived/nroy_rank_correlations.csv` |

**Extraction** from raw simulation outputs (not distributed; the results are in `derived/`)

| Script | Experiment | Writes to `derived/` |
|---|---|---|
| `extract1_baseline_episodes.py` | Production, Baseline scenario | `baseline_episodes.csv.gz`, `figure1_values.json` |
| `extract2_scenario_episodes.py` | Production, all 8 scenarios | `scenario_episodes.csv.gz`, `scenario_infectious_time.csv.gz` |
| `extract3_burden_by_replicate.py` | Production, all 8 scenarios | `burden_by_replicate.csv.gz` |
| `extract4_burden_both_prevalence.py` | Production, all 8 scenarios | `burden_by_replicate_both_prevalence.csv.gz` |
| `extract5_parameter_propagation.py` | Parameter propagation | `burden_by_parameter_set.csv.gz` |
| `extract6_loop_closure_sweep.py` | Loop-closure sweep | `loop_closure_sweep_incidence.csv.gz` |

`check_derived.py` compares re-extracted tables with the committed versions (git HEAD), ignoring row order.

**Figures and tables** (read `derived/` only, run in minutes)

| Script | Manuscript item |
|---|---|
| `figure1_cascade_events_and_trajectories.py` | Figure 1A, 1B |
| `figure2a_calibration.py` | Figure 2A |
| `figure2b_table2_cascade_benchmarks.py` | Figure 2B, Table 2 |
| `figure3_loop_group_burden.py` | Figure 3 (use `--plot-only` to skip the raw-output step) |
| `figure4_three_lens_intervention_impact.py` | Figure 4 |
| `figure5_end_tb_progress.py` | Figure 5 |
| `figureS2_loop_closure_sensitivity.py` | Appendix Figure S-2 |
| `tableS8_propagation_contrasts.py` | Appendix Table S-8 |

Figures are written to `figures/reproduced/` under the same file names as the manuscript copies in `figures/`.

## derived/

| File | Written by | Read by |
|---|---|---|
| `baseline_episodes.csv.gz` | `extract1` | `figure2b`, `figure3` |
| `baseline_population_years_2015_2024.json` | `figure3` (first step) | `figure3` |
| `burden_by_parameter_set.csv.gz` | `extract5` | `figure5`, `tableS8` |
| `burden_by_replicate.csv.gz` | `extract3` | `figureS2` |
| `burden_by_replicate_both_prevalence.csv.gz` | `extract4` | `figure5` |
| `calibration_final_fit_report.txt` | `calib4`, `calib5` | — |
| `cohort_2015_2024_numbers.json` | — | headline numbers for the 2015–2024 cohort |
| `figure1_values.json` | `extract1` | `figure1` |
| `figure3_values.json`, `figure4_values.csv`, `figure5_parameter_uncertainty.csv` | the figure script of the same number | — (numbers quoted in the text) |
| `loop_closure_sweep_incidence.csv.gz` | `extract6` | `figureS2` |
| `nroy_marginals.csv` | `calib5` | `calib6` |
| `nroy_parameter_sets.json` | the 20 parameter sets drawn from the not-ruled-out region | `simulation/3b_parameter_propagation_run.py` |
| `nroy_rank_correlations.csv`, `nroy_structure_report.txt` | `calib6` | — |
| `production_baseline_burden.csv.gz` | `figure2a` (first run) | `figure2a` |
| `production_parameters.json` | `calib4` | — |
| `propagation_contrasts.csv` | `tableS8` | — |
| `scenario_episodes.csv.gz` | `extract2` | `figure4` |
| `scenario_infectious_time.csv.gz` | `extract2` | — |
| `table2_values.json` | `figure2b` | — |
