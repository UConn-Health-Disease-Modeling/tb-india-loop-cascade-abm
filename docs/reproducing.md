# Reproducing the TB India ABM Analysis

This guide walks through reproducing all manuscript figures and supplementary
tables from the EMOD-TBHIV simulation outputs.

## Prerequisites

- Python 3.7+ with packages from `requirements.txt`
- Pre-generated EMOD simulation outputs (or access to COMPS to run new ones)
- Jupyter Lab or Jupyter Notebook

## Step 1 — Set up the environment

```bash
git clone <repo-url>
cd release/
pip install -r requirements.txt
```

## Step 2 — Generate or obtain EMOD simulation outputs

The analysis expects EMOD outputs in the parent directory of the notebooks, in
this structure:

```
<experiment_uid>/
├── PropertyReport.json
├── ReportEventRecorder.csv
└── ...
```

### Option A: Run new simulations
Open `notebooks/01_simulation/India_TB_v5_simulation.ipynb`, configure your
COMPS credentials at the top, and **Restart & Run All**. This submits 8
intervention scenarios; expected runtime is several hours on COMPS.

### Option B: Use pre-generated outputs
If pre-generated outputs are provided as a separate download, extract them into
the project root so the notebook can find them via the experiment-UID pattern
in cell 3.

## Step 3 — Calibrate (one-time only)

If you change parameter priors or add new calibration targets, re-run:

```
notebooks/02_calibration/history_matching.ipynb
```

This regenerates the calibration database. The best-fit parameter set used in
the manuscript (wave1_pt8) is hard-coded in `src/params.py`, so this step is
only needed if you change the calibration design.

## Step 4 — Generate main paper figures

Open and run:

```
notebooks/03_figures/India_TB_v5_figures_main.ipynb
```

This produces 13 main paper figures saved to `figures/`:

| Cell | Figure | Output |
|---|---|---|
| Cell 11 | Fig 1 — Calibration overlay | calibration_overlay.png |
| Cell 14 | Fig 2 — Delay decomposition by loop group | Fig2_delay_decomposition.png |
| Cell 17 | Fig 3 — Loop-group composition | Fig2a_events_and_groups.png |
| Cell 19 | Fig 4 — Inside the 85% | FigG_inside_the_85.png |
| Cell 21 | Fig 5 — Transmission tail | FigH_transmission_tail.png |
| Cell 23 | Fig 6 — Cascade efficiency | Fig4_cascade_metrics.png |
| Cell 27 | Fig 7 — Loop-group waterfall | Fig5b_loopgroup_waterfall.png |
| Cell 29 | Fig 8 — Delay waterfall | Fig5c_delay_waterfall.png |
| Cell 31 | Fig 9 — Triple-lens heatmap | FigDEF_combined_triple_lens.png |
| Cell 35 | Fig 10 — Path to End-TB | FigP1_path_to_endtb.png |
| Cell 37 | Fig 11 — Strategic priority | FigP2_strategic_matrix.png |
| Cell 39 | Fig 12 — Pareto frontier | FigP3_pareto_frontier.png |
| Cell 41 | Fig 13 — Loop sensitivity | FigP7_loop_sensitivity.png |

## Step 5 — Generate appendix tables

Open and run:

```
notebooks/03_figures/India_TB_v5_figures_appendix.ipynb
```

This produces:

| Cell | Output | Location |
|---|---|---|
| Cell 14 | Tables S-1, S-2, S-3, S-4 (4 sheets) | `tables/calibration_tables_S1_S4.xlsx` |
| Cell 16 | Tables S-7.1 … S-7.7 (7 sheets) | `tables/Table_A1_model_parameters.xlsx` |
| Cell 16 | Unified parameter listing | `tables/Table_A1_unified.csv` |
| Cell 18 | Consolidated outcomes | `tables/appendix_table.csv` + `.xlsx` |

## Troubleshooting

### Notebook can't find simulation outputs
Make sure the experiment UID directory is in the project root (one level above
the notebook). The notebook scans for a directory matching the UUID-like
pattern `xxxxxxxx-xxxx-xxxx-xxxx-xxxxxxxxxxxx`.

### `scenario_results` not defined (figures notebook)
Run the **Figure 5 — Intervention Scenarios** cell first (cell 27 in main).
This builds `scenario_results`, which is required by P1, P3, P7, and the
appendix outcomes table.

### Figures look stale after parameter change
Always do **Kernel → Restart & Run All** rather than running individual cells.
Calibration values are loaded once at the top; partial re-execution can leave
stale variables in memory.

## Computational requirements

- Calibration: ~4 hours on 8-core CPU (Latin hypercube + GP emulators)
- Single scenario simulation: ~30 min on COMPS (5 replicates)
- All 8 scenarios in parallel: ~1 hour on COMPS
- Figure generation: <5 minutes total (post-processing only)
