# simulation/

The EMOD-TBHIV model and the notebooks that ran it on IDM COMPS. Run the notebooks from inside this folder.

## Notebooks, in run order

The notebooks share the same model-building cells. Which experiment a notebook runs is set by `SCENARIO`
(and the design file `CAL_DESIGN_CSV`) in its run-settings cell.

| Notebook | `SCENARIO` | Experiment | COMPS ID |
|---|---|---|---|
| `1_calibration_wave1.ipynb` | `wave1_cal` | Calibration wave 1, design `calibration/wave1/design.csv` | `af707f22-16bc-f111-92f1-000d3af5294c` |
| `2_calibration_wave2.ipynb` | `wave2_cal` | Calibration wave 2, design `calibration/wave2/design.csv` | `496ded25-37bc-f111-92f1-000d3af5294c` |
| `3_production_scenarios.ipynb` | `combined` | Production: 8 intervention scenarios at the calibrated point | `d892bdf7-18bd-f111-92f1-000d3af5294c` |
| `3b_parameter_propagation_run.py` | — | Run from notebook 3: the production task over 20 parameter sets from the not-ruled-out region | `f8dfe706-d7c0-f111-92f1-000d3af5294c` |
| `4_loop_closure_sweep.ipynb` | `loopsweep` | Loop-closure sweep (Appendix Figure S-2) | `941e6423` |

The time-varying cascade schedules used to build the campaign are imported from
`analysis/campaign_schedules.py`.

## Model files

| File or folder | Contents |
|---|---|
| `params.py` | Model parameters and intervention timing |
| `set_config.py` | EMOD configuration |
| `add_tb_drug_type.py`, `drug.json` | TB drug regimen definitions |
| `manifest.py` | Paths to the EMOD binary, schema, assets and input data |
| `Assets/` | Demographics overlays uploaded with every simulation |
| `data/` | India mortality data (model input) and prevalence-survey data the notebooks compare the model against |
| `python_scripts/` | EMOD pre- and post-processing scripts run on COMPS |
| `dtk_centos.id` | Container image used on COMPS |
| `requirements-simulation.txt` | Python environment for running the notebooks |
| `requirements.txt` | Full package list of the original development environment (referenced by `manifest.py`) |
