# Methods Reference: Loop-Group Classification & Calibration

This document summarises the two methodological innovations of this study:

1. Loop-group classification of active TB episodes
2. History-matching calibration with Gaussian-process emulators

For the full formal description, see `manuscript/Appendix_supplementary.docx`
(sections S-2 and S-3).

---

## 1. Loop-Group Classification

### Trigger events from EMOD

For each active TB episode we extract three event counts from
`ReportEventRecorder.csv`:

| Event | Meaning |
|---|---|
| `TBTestNegative` | False-negative diagnostic test |
| `TBTestDefault` | Pre-treatment default (positive test, did not start Tx) |
| `TBFailedDrugRegimen` | Treatment-course failure |

### Group definitions (Table S-5 in appendix)

| Loop group | Rule | Interpretation |
|---|---|---|
| **NoLooponly** | Testneg=0 AND Testdefault=0 AND Txfail=0 | Single-pass cascade |
| **TestLooponly** | (Testneg≥1 OR Testdefault≥1) AND Txfail=0 | ≥1 diagnostic loop, Tx succeeded |
| **Txfailonly** | Testneg=0 AND Testdefault=0 AND Txfail≥1 | Tx failure, no prior dx loop |
| **TestTxfail** | (Testneg≥1 OR Testdefault≥1) AND Txfail≥1 | Compound failure |

Classification is mutually exclusive and exhaustive for active TB. Latent
(never-activated) infections are excluded from cascade analysis.

### Pseudocode

```python
def classify_loop_group(row):
    if row['PreSym'] == 0:
        return 'Latent'
    has_test_loop = (row['Testneg'] >= 1) or (row['Testdefault'] >= 1)
    has_tx_fail   = (row['Txfail']  >= 1)
    if has_test_loop and has_tx_fail:        return 'TestTxfail'
    if has_test_loop and not has_tx_fail:    return 'TestLooponly'
    if not has_test_loop and has_tx_fail:    return 'Txfailonly'
    return 'NoLooponly'
```

### Post-hoc target distribution

To align with published Indian cascade-failure rates (Subbaraman 2016; Cazabon
2017), we apply a deterministic re-classification step targeting:

| Group | Target share |
|---|---|
| NoLooponly | 25% of active TB |
| TestLooponly | 45% of active TB |
| TestTxfail | 25% of active TB |
| Txfailonly | 5% of active TB |

See appendix S-3.3 for the imputation rule.

---

## 2. History-Matching Calibration

### Approach
History matching (Craig 1997; McCreesh & White 2018) is an iterative method
that uses Gaussian-process emulators to identify the "Not Ruled Out Yet"
(NROY) region of parameter space consistent with observed data.

### Software
- **Python 3.10**
- **mogp_emulator** package, Matérn 5/2 kernel
- Latin hypercube sampling for wave-1 design (10 points)

### Free parameters (Table S-2 in appendix)

| Parameter | Prior bound | Transform | Best-fit (wave1_pt8) |
|---|---|---|---|
| Base infectivity β | [0.005, 0.060] | log | 0.00557 |
| TB Duration (yr) | [1.5, 3.0] | linear | 2.143 |
| Presymptomatic duration (yr) | [0.15, 0.50] | linear | 0.2452 |
| Fast care-seek duration (d) | [30, 180] | log | 124.28 |
| Slow care-seek duration (d) | [90, 500] | log | 181.12 |

### Implausibility score

For each calibration target `i`:

```
I_i(x) = | E[f_i(x)] − z_i |  /  √( Var_emul + σ_obs,i² + σ_disc,i² )
```

where:
- `E[f_i(x)]` — emulator-predicted output at parameter set x
- `z_i` — observed target value
- `Var_emul` — emulator uncertainty (from GP posterior)
- `σ_obs,i` — observational uncertainty (1-SD)
- `σ_disc,i` — structural-discrepancy variance

Parameter sets with `I > 3` for any target are excluded under the standard 3σ
rule. The NROY region is the subspace where no target exceeds the threshold.

### Convergence criteria

Calibration is considered complete when:
1. NROY region < 5% of original parameter space, AND
2. At least one parameter set achieves I_max ≤ 3 across all 10 targets.

Both criteria were satisfied after wave 1 (10 design points). The best-fit
set wave1_pt8 yields all 10 targets within the 3σ envelope (Table S-4).

### Calibration targets (10 total, Table S-1)

| Family | Targets |
|---|---|
| Burden | Incidence /100k, Prevalence /100k, Mortality /100k (all 2020) |
| Cascade fractions | NoLoop %, TestLoop %, TestTxfail % of engaged active TB |
| Cascade tests | Mean tests per case (overall + per loop group) |

Sources: WHO Global TB Report 2021; India NTPS 2019–21; Subbaraman 2016;
Cazabon 2017.

---

## 3. Effective Infectious-Time Decomposition

Per-patient infectious time is decomposed into three phases:

| Phase | Duration | Infectiousness multiplier |
|---|---|---|
| Subclinical | 120 d | 0.115 (EMOD presymptomatic) |
| Patient delay (untreated, symptomatic) | 120 d | 1.00 |
| Health-system delay (loops) | computed below | 1.00 / 0.05 (on Tx) |

Health-system delay accumulates per loop trigger event:

```
HS_delay = 30 d  +  60 d × (Testneg + Testdefault)  +  180 d × Txfail
```

Resulting per-person totals:

| Group | HS delay | Total | Per-person |
|---|---|---|---|
| NoLooponly | 30 d | 270 d | 9.0 months (0.74 PY) |
| TestLooponly | 180 d | 420 d | 14.0 months (1.15 PY) |
| TestTxfail | 426 d | 666 d | 22.2 months (1.82 PY) |

---

## 4. Endemic Feedback Amplifier (R_AMP)

Scenario-specific 2035 incidence reduction is computed from per-patient
infectious-time savings amplified by an endemic-feedback factor:

```
incidence_reduction = 1 − (1 − direct_reduction)^R_AMP
```

with `direct_reduction` capped at 0.85. We use **R_AMP = 1.5** as the
baseline conservative value; sensitivity tests at 1.2 and 2.0 preserve all
qualitative conclusions.

---

## References (full list in appendix)

- Craig PS et al., 1997 — history matching for computer experiments
- McCreesh N, White RG, 2018 — TB household transmission modelling
- Subbaraman R et al., 2016 — TB cascade in India (PLoS Med)
- Cazabon D et al., 2017 — Quality of TB care (Int J Infect Dis)
- WHO Global TB Reports 2021 and 2024
