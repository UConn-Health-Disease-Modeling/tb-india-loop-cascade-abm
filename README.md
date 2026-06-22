# Cascade-of-Care Loops as Transmission Multipliers in Tuberculosis: India ABM

**A loop-group-resolved agent-based model of tuberculosis transmission in India**

This repository contains the simulation code, calibration framework, analysis
notebooks, manuscript figures, and supplementary appendix tables for the TB
India agent-based modelling (ABM) study using EMOD-TBHIV.

---

## Table of Contents

1. [Overview](#overview)
2. [Repository structure](#repository-structure)
3. [Installation](#installation)
4. [Reproducing the analysis](#reproducing-the-analysis)
5. [Manuscript and appendix](#manuscript-and-appendix)
6. [Key findings (at a glance)](#key-findings-at-a-glance)
7. [Citation](#citation)
8. [Contact](#contact)

---

## Overview

We extend an EMOD-TBHIV individual-based simulation calibrated to Indian TB
incidence, prevalence, and notification time-series (1997–2024) to track every
patient episode through diagnostic and treatment events. Active TB episodes are
partitioned into four **loop groups** based on calibrated counts of false-negative
diagnostic tests, pre-treatment defaults, and treatment failures: **NoLooponly,
TestLooponly, TestTxfail, Txfailonly**.

We simulate 8 intervention scenarios (Baseline, ACF, Dx, Tx, and all combinations)
starting in 2025 and project outcomes to 2035, quantifying:

- The hidden patient-level heterogeneity in WHO "treatment success rate"
- Pareto-style concentration of cascade-driven transmission (Gini ≈ 0.45)
- Mechanistic complementarity of intervention bundles (ACF, Dx, Tx target
  different loop populations)
- Path to End-TB 2035 (incidence reduction targets)

---

## Repository structure

```
release/
├── README.md                              This file
├── LICENSE                                MIT licence
├── requirements.txt                       Python dependencies
│
├── manuscript/                            Final manuscript drafts (.docx)
│   ├── India_TB_model_main.docx
│   └── Appendix_supplementary.docx        (structured: S-1 ... S-7)
│
├── notebooks/                             Analysis notebooks (Jupyter)
│   ├── 01_simulation/                     EMOD simulation runs
│   │   └── India_TB_v5_simulation.ipynb
│   ├── 02_calibration/                    History-matching calibration
│   │   ├── history_matching.ipynb
│   │   ├── calibration_framework.ipynb
│   │   └── calibration.ipynb
│   └── 03_figures/                        Figure & table generation
│       ├── India_TB_v5_figures_main.ipynb         (13 main paper figures)
│       └── India_TB_v5_figures_appendix.ipynb     (tables S-1..S-4, S-7.x)
│
├── src/                                   Python modules (EMOD config + helpers)
│   ├── params.py                          EMOD parameter dictionary
│   ├── manifest.py                        Path constants
│   ├── set_config.py                      EMOD config builder
│   ├── add_tb_drug_type.py
│   ├── run_base_simulation.py             Simulation entry point
│   ├── recalibration_plan.py
│   ├── spatial_plotter.py
│   └── Base_analysis.py
│
├── data/                                  Input data files
│   ├── demographics.json                  EMOD demographics
│   ├── drug.json                          Drug & treatment parameters
│   ├── default_config.json                EMOD base configuration
│   ├── BGW.xlsx                           Birth/growth/weight data
│   ├── regional_migration.bin             Spatial migration matrix
│   └── regional_migration.bin.json
│
├── figures/                               Final manuscript figures (PNG)
│   ├── calibration_overlay.png            Fig 1
│   ├── Fig2_delay_decomposition.png       Fig 2 (loop-group infectious time)
│   ├── Fig2a_events_and_groups.png        Fig 3 (composition)
│   ├── FigG_inside_the_85.png             Fig 4 (Pareto disparity narrative)
│   ├── FigH_transmission_tail.png         Fig 5 (Lorenz + distribution)
│   ├── Fig4_cascade_metrics.png           Fig 6 (test distribution)
│   ├── Fig5b_loopgroup_waterfall.png      Fig 7
│   ├── Fig5c_delay_waterfall.png          Fig 8
│   ├── FigDEF_combined_triple_lens.png    Fig 9 (3-panel heatmap)
│   ├── FigP1_path_to_endtb.png            Fig 10 (trajectory)
│   ├── FigP2_strategic_matrix.png         Fig 11
│   ├── FigP3_pareto_frontier.png          Fig 12
│   └── FigP7_loop_sensitivity.png         Fig 13
│
├── tables/                                Manuscript tables (Excel/CSV)
│   ├── Table_A1_model_parameters.xlsx     7-sheet parameter listing (S-7)
│   ├── Table_A1_unified.csv               Single-table version
│   ├── calibration_tables_S1_S4.xlsx      Tables S-1 ... S-4
│   ├── appendix_table.csv                 Consolidated outcomes (CSV)
│   └── appendix_table.xlsx                Same as Excel
│
└── docs/                                  Methods documentation
    ├── methods.md                         Loop-group definitions + calibration
    └── reproducing.md                     Step-by-step reproduction guide
```

---

## Installation

### Requirements
- **Python 3.7+** (tested on 3.7 and 3.10)
- **EMOD-DTK v2.21** (for running new simulations; not required if reading
  pre-generated outputs)
- **Jupyter Lab** or **Jupyter Notebook**

### Dependencies
```bash
pip install -r requirements.txt
```

Key packages: `numpy`, `pandas`, `matplotlib`, `scipy`, `pdfplumber`,
`python-docx`, `openpyxl`, `mogp_emulator` (for calibration).

---

## Reproducing the analysis

The analysis runs in three sequential stages:

### Stage 1 — Run EMOD simulation
```
notebooks/01_simulation/India_TB_v5_simulation.ipynb
```
Configures and submits 8 intervention scenarios to COMPS (or local EMOD).
Outputs: `ReportEventRecorder.csv`, `PropertyReport.json` per scenario.

### Stage 2 — Calibrate parameters (history matching)
```
notebooks/02_calibration/history_matching.ipynb
```
Generates wave-based Latin hypercube samples, fits Matérn 5/2 GP emulators,
identifies NROY region, selects best-fit parameter set (wave1_pt8).

### Stage 3 — Generate figures and tables
```
notebooks/03_figures/India_TB_v5_figures_main.ipynb        # 13 main figures
notebooks/03_figures/India_TB_v5_figures_appendix.ipynb    # Tables S-1..S-7
```

Both notebooks are self-contained: they load simulation outputs, build
`loop_df` / `delay_df` / `episode_df`, and run all analyses end-to-end via
**Kernel → Restart & Run All**.

See **docs/reproducing.md** for the full step-by-step guide.

---

## Manuscript and appendix

| Document | Section | Contents |
|---|---|---|
| **Main paper** | Introduction | TB burden, cascade framework limits, contributions |
| | Methods | EMOD model, calibration, loop-group classification |
| | Results | 8 scenarios, Pareto, intervention specificity, End-TB path |
| | Discussion | Hidden tail, find-vs-fix, loop multiplier, limitations |
| **Appendix** | S-1 Model overview | Population, natural history, cascade |
| | S-2 Calibration | Targets (S-1), parameters (S-2), GP emulator protocol, fit (S-3/S-4) |
| | S-3 Loop-group classification | Definitions (S-5), algorithm pseudocode |
| | S-4 Intervention scenarios | Eight scenarios (S-6) |
| | S-5 Effective infectious-time | Phase definitions, HS delay formula |
| | S-6 Complete parameter listing | 7 subsection tables (S-7.1 … S-7.7) |
| | S-7 Sensitivity analyses | R_AMP, time-cost constants, target distributions |
| | References | 11 cited works |

---

## Key findings (at a glance)

- Of 88,529 simulated active TB episodes, only **25% progressed via single-pass
  cascade** (NoLooponly). 45% experienced diagnostic loops (TestLooponly), and
  25% experienced both diagnostic loops and treatment failure (TestTxfail).
- **Per-patient effective infectious time** differed substantially:
  9.0 months (NoLoop) vs 14.0 (TestLoop) vs 22.2 (TestTxfail).
- The 25% TestTxfail minority generated **42% of cascade-driven transmission**
  (Pareto Gini ≈ 0.45).
- ACF, Dx, and Tx interventions are **mechanistically complementary** rather
  than substitutable: each acts on a different loop population.
- **Comprehensive ACF+Dx+Tx achieves 79% incidence reduction by 2035**; no
  single-component bundle reaches the End-TB 90% target.

---

## Citation

If you use this code or analysis in your work, please cite:

> [Author list]. Cascade-of-Care Loops as Transmission Multipliers: A Loop-Group
> Resolved Agent-Based Model of Tuberculosis in India. [Journal]. [Year];[Vol]:
> [pages]. doi:[DOI]

---

## Contact

[Corresponding author name and email]

For EMOD model questions: see <https://github.com/InstituteforDiseaseModeling/EMOD>.
