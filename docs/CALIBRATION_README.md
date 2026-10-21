# Wave-2 Calibration Workflow — 10 Free Parameters (200+100 design)

Standard 2-wave history-matching protocol following Vesga et al. 2019.

| Wave | Design points | Reps | Total sims | Wall-clock |
|------|---------------|------|------------|------------|
| 1    | **200 (LHS)** | 10   | **2,000**  | ~6 h       |
| 2    | 100 (NROY-LHS)| 10   | 1,000      | ~3 h       |
| Validation | 1 (best-fit) | 100 | 100     | ~30 min    |
| **Total** | 301          | —    | **3,100**  | **~9.5 h** |

Asymmetric 200 + 100 design (Vesga 2019 style) chosen to:
- Maximise wave-1 emulator R² (200 LHS in 10D gives discrepancy 0.019 vs
  0.035 for 100 LHS), improving the likelihood that all 10 emulators meet
  the > 0.9 LOOCV R² convergence criterion in a single wave-1 pass.
- Reduce the probability of needing a third wave by giving the wave-1 NROY
  identification more design support across the full 10D prior.
- Match the design protocol of the most directly comparable published TB-
  cascade calibration (Vesga et al. 2019; 200 wave-1 + 100 wave-2 over a
  13-parameter prior for India TB control modelling).

## Files

| File | Purpose |
|------|---------|
| `helpers/calibration_wave1_lhs.py`     | Generates `wave1_design.csv` (100 LHS points) |
| `helpers/wave1_design.csv`             | ✅ Generated (100 × 10 params) |
| `01_India_TB_v12_simulation_calibration.ipynb` | Submits sims to COMPS in `SCENARIO='wave1_cal'` or `'wave2_cal'` mode |
| `helpers/calibration_emulator_nroy.py` | Fits GP emulator, computes NROY |
| `helpers/calibration_wave2_lhs.py`     | Generates `wave2_design.csv` from wave-1 NROY |

## Workflow

### Step 1 — Wave-1 LHS design (✅ done)

```bash
python helpers/calibration_wave1_lhs.py
```

Produces `helpers/wave1_design.csv` (100 design points × 10 params).

### Step 2 — Wave-1 sweep on COMPS

In `01_India_TB_v12_simulation_calibration.ipynb`:

```python
SCENARIO       = 'wave1_cal'
CAL_DESIGN_CSV = 'helpers/wave1_design.csv'
replicates     = 10
```

→ Kernel → Restart & Run All. Submits 100 × 10 = 1,000 sims to COMPS.

### Step 3 — Wave-1 emulator + NROY

After wave-1 experiment completes:

```bash
python helpers/calibration_emulator_nroy.py <wave1_exp_uid> --wave 1
```

Produces:
- `<exp>/calibration/design_outputs.csv`
- `<exp>/calibration/emulator_metrics.csv` (per-target CV R²)
- `<exp>/calibration/nroy_samples.csv` (20,000 candidates with NROY flag)
- `<exp>/calibration/nroy_summary.txt` (volume reduction, best-fit point)

### Step 4 — Wave-2 LHS within NROY

```bash
python helpers/calibration_wave2_lhs.py <wave1_exp_uid> --n 100
```

Produces `helpers/wave2_design.csv` (100 design points constrained to NROY).

### Step 5 — Wave-2 sweep

In v12 notebook:

```python
SCENARIO       = 'wave2_cal'
CAL_DESIGN_CSV = 'helpers/wave2_design.csv'
replicates     = 10
```

→ Run All. Submits 100 × 10 = 1,000 sims.

### Step 6 — Wave-2 emulator + final NROY

```bash
python helpers/calibration_emulator_nroy.py <wave2_exp_uid> --wave 2
```

Final NROY volume should be < 5% of prior. Best-fit point becomes the
calibrated parameter set.

### Step 7 — Validation

Take the best-fit point from wave-2, plug into v10/v11 production subsample
dict, run 100 replicates × 8 scenarios as final validation.

## Calibration Targets (Appendix Table S-1)

| # | Target                  | z      | σ_obs | σ_disc |
|---|-------------------------|--------|-------|--------|
| 1 | Incidence/100k 2020     | 188    | 25    | 20     |
| 2 | Prevalence/100k 2020    | 304    | 20    | 30     |
| 3 | Mortality/100k 2020     | 33     | 4.5   | 4.0    |
| 4 | NoLoop fraction         | 0.20   | 0.05  | 0.05   |
| 5 | TestLoop fraction       | 0.50   | 0.08  | 0.08   |
| 6 | TestTxfail fraction     | 0.30   | 0.05  | 0.05   |
| 7 | Mean tests/case overall | 3.5    | 0.6   | 0.5    |
| 8 | Mean tests/case NoLoop  | 1.0    | 0.3   | 0.3    |
| 9 | Mean tests/case TestLoop| 3.0    | 0.7   | 0.6    |
|10 | Mean tests/case TestTxf | 4.5    | 1.0   | 0.8    |

Note: cascade fraction targets (4-6) now use the unified post-hoc target
distribution (0.20 / 0.50 / 0.30) per user decision to remove the post-hoc
reclassification step.

## Free Parameters (Wave-2 priors)

| # | Parameter                                | Min   | Max    | Transform |
|---|------------------------------------------|-------|--------|-----------|
| 1 | Base_Infectivity_Constant                | 0.005 | 0.100  | log       |
| 2 | TB Duration                              | 1.5   | 4.0    | linear    |
| 3 | TB Presymptomatic Duration Years         | 0.15  | 0.5    | linear    |
| 4 | Care Seeking Fast Duration Days          | 10    | 30     | log       |
| 5 | Care Seeking Slow Duration Days          | 60    | 180    | log       |
| 6 | sens_smear_pos_GXL                       | 0.45  | 0.85   | linear    |
| 7 | treatment_fraction_LQ                    | 0.30  | 0.85   | linear    |
| 8 | TB_Slow_Progressor_Rate                  | 5e-6  | 5e-5   | log       |
| 9 | Post_Infection_Acquisition_Multiplier    | 0.20  | 0.80   | linear    |
|10 | TB_Fast_Progressor_Fraction_Adult        | 0.05  | 0.20   | linear    |

## Dependencies

- scipy >= 1.7 (LatinHypercube, qmc)
- scikit-learn (Matérn-5/2 GP fallback)
- mogp_emulator (preferred, faster)

Install: `pip install mogp_emulator` (optional but recommended).
