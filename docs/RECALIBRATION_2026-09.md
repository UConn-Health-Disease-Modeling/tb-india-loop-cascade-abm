# Trend re-calibration and 2025 intervention start — run plan (2026-09-09)

Decisions taken: recalibrate against the WHO incidence trend (option b), re-run the
scenarios with a 2025 intervention start, and adopt the recomputed episode-level
numbers for the manuscript.

## 1. Why the published baseline rises

Cascade performance is not constant over time in the model — it gets **worse** at
GeneXpert introduction. The treatment-initiation fraction is hard-coded per epoch:

| Epoch | treatment_fraction |
|---|---|
| Pre-DOTS (to 1997) | 0.60 |
| DOTS (1997–2012) | 0.70 |
| GeneXpert (2012–) | **0.55** |

The 0.55 was hand-tuned to hit the tests-per-case cascade target
(notebook comment: *"v9: lower for normal-like tests/case distribution"*). The
epidemiological cost is visible in the simulated trajectory: incidence falls
195 → 158 per 100,000 between 2010 and 2015 as GeneXpert sensitivity ramps in,
then climbs back to 198 by 2020 and keeps rising as the lower treatment fraction
takes over. Net trend 2000–2020: **+0.05%/yr**, against WHO's **−2.43%/yr**.

## 2. What changed

### `params.py` (backed up to `params_backup_pre0909.py`)
- `intervention_day` = `dots_start + 28*365` → **2025** (was +22 y = 2019, which
  contradicted the manuscript's stated 2025 scenario start).
- Two new calibration parameters appended to `params.parameters`:
  `tf_gain_per_year` (0–0.012, guess 0.006) and
  `cs_gain_per_year` (0–0.050, guess 0.020).

### `helpers/campaign_trend.py` (new)
Replaces the two hand-set epochs with a monotone cascade-improvement schedule
1997 → 2024, in 3-year sub-epochs:
- `cascade_schedule(params, tf_gain, cs_gain, dx_boost, tx_boost, sens_pos_GXL)`
  returns non-overlapping diagnostic blocks with their own treatment fraction and
  sensitivities. The GeneXpert sensitivity ramp is preserved exactly (each
  sub-epoch gets the interpolated start/end sensitivities, so splitting the era
  does not restart the ramp).
- `careseeking_schedule(params, cs_gain)` returns the per-sub-epoch multiplier on
  the care-seeking probability.
- Zero is inside both priors, so "no secular improvement" — the published
  behaviour — remains reachable and the calibration can reject the mechanism.

**Scenario definition change (required):** the treatment-support scenario can no
longer be a fixed 0.85, because the calibrated 2024 baseline may already exceed
it (with `tf_gain` = 0.006/yr the 2024 value is 0.86, and a fixed 0.85 target
would be a *worsening*). It is now defined relative to the calibrated baseline:
`tf_Tx = 1 − 0.5 × (1 − tf_baseline)`, i.e. treatment support halves the residual
pre-treatment loss. `TX_BOOST_LOSS_REDUCTION` in `campaign_trend.py`.
Table 3 and the Methods text will need this wording.

### `helpers/targets_trend.py` (new)
Calibration targets change from three cross-sectional 2020 anchors to a
four-point WHO incidence trajectory plus the 2020 prevalence and mortality
anchors:

| target | value | sigma |
|---|---|---|
| incidence_2000 | 322 | 54.7 |
| incidence_2010 | 276 | 46.9 |
| incidence_2015 | 237 | 40.3 |
| incidence_2020 | 197 | 33.5 |
| prevalence_2020 | 304 | 36.1 |
| mortality_2020 | 33 | 6.0 |

Source: WHO Global TB Programme estimates for India (World Bank mirror,
indicator `SH.TBS.INCD`). 2005 (311) and 2023 (192) are held out as a check on
unused years. sigma = 17% of the point estimate, the same ratio as the published
2020 target (32/188).

*Vintage note:* the manuscript's published 2020 target was 188 (Global TB Report
2021); the current series puts 2020 at 197. The whole trajectory here comes from
one vintage so the trend is internally consistent, and the difference is well
inside one sigma.

Split is unchanged in spirit: **6 burden targets are emulator-gated calibration
targets; the 7 cascade targets stay validation-only** (3-sigma check at the
best fit).

### `helpers/calibration_wave1_lhs_trend.py` (new)
200-point Latin hypercube over the 10 free parameters (the previous 8 plus the
two improvement parameters); `sens_smear_pos_GXL` and `treatment_fraction_LQ`
are written as constant columns so the notebook's sweep shape is unchanged.
Output: `helpers/wave1_design_trend.csv` (already generated).
LHS is implemented in numpy — scipy is not installed in this environment.

### `helpers/compute_targets_trend.py` (new)
Same pipeline as `compute_wave1_targets.py`, extracting incidence at 2000, 2005,
2010, 2015, 2020 and 2023 instead of 2020 alone, and reporting the implied
annualised trend per design point. Checkpointed and resumable.

## 3. Notebook: patched — `01_India_TB_v20_simulation_trend.ipynb`

v19 is left untouched; v20 is the run notebook. Applied edits:

- **Cell 5** — `CAL_DESIGN_CSV = 'helpers/wave1_design_trend.csv'`.
- **Cell 14 `build_camp`** — signature gains `tf_gain_per_year=0.0,
  cs_gain_per_year=0.0`; imports `cascade_schedule` / `careseeking_schedule`.
  - the two DOTS-era-onward HSB blocks are replaced by one loop over
    `careseeking_schedule(...)` (per sub-epoch, both quality streams);
  - the Epoch-2 `DiagnosticTreatNeg` pair and the Epoch-3 `RampDTN` pair are
    replaced by one loop over `cascade_schedule(...)`.
  - **The GeneXpert era stays on `RampDTN`.** `RampDTN` is what emits
    `ProviderOrdersTBTest` / `ProviderTestNoR`; the resistance diagnostics are
    triggered on those events and `cascade_from_rer` counts them. Swapping in
    `DiagnosticTreatNeg` there would have silently removed both. Each sub-epoch
    receives interpolated start/end sensitivities (`sens_*_start` → `sens_*_end`,
    passed as RampDTN's `*2` and plain arguments) so the 8-year GeneXpert ramp
    continues across sub-epochs instead of restarting: 0.700 → 0.719 → 0.738 →
    0.750 for high-quality smear-positive sensitivity.
  - pre-DOTS epoch, all treatment blocks, resistance diagnostics, ACF and the
    cash tracker are unchanged.
- **Cell 19 `map_sample_to_model_input`** — the four competing
  `create_campaign_from_callback` calls are replaced by a `camp_kwargs` dict and
  a single registration at the end of the function (see §4), with the two new
  parameters wired in.

Checks run after patching: all three cells compile; no orphan references to the
removed `_gxl_tf` / `_base_tf` / `_gxl_sens_*` locals; resistance diagnostics,
treatment blocks, ACF and HSB blocks all still present; both schedules partition
1997–2050 with no gaps or overlaps; `from helpers.campaign_trend import ...`
resolves from the repo root.

## 4. The campaign-callback bug (fixed in v20, still needs verification)

`map_sample_to_model_input` used to call
`simulation.task.create_campaign_from_callback(partial(build_camp, X=...))`
**separately** for `Care Seeking Fast Duration Days`, `Care Seeking Slow Duration
Days`, `sens_smear_pos_GXL` and `treatment_fraction_LQ`.

Registration is not additive, so depending on how this emodpy build behaves,
either

- **last call wins** — the two calibrated care-seeking parameters were never
  applied in the 8-parameter waves or in the production run, and the campaign was
  built with the defaults 0.0208 / 0.05247 (48 d and 19 d), or
- **the callbacks accumulate on the module-level `camp` object** — every
  intervention was added up to four times.

I could not settle which: emodpy is not installed on this machine and the symptom
is not distinguishable from the stored outputs. v20 collects every
campaign-affecting parameter into `camp_kwargs` and registers the callback once,
which is correct under either behaviour.

**Verify before spending a wave.** Run 2–3 single simulations with deliberately
extreme care-seeking values (e.g. `Care Seeking Fast Duration Days` 10 vs 30) and
confirm the patient-delay distribution in `ReportEventRecorder` actually moves.
If it does not move under v19 but does under v20, the bug was real and the
published calibration inherits it — which would need saying in the paper.

## 5. Run sequence

1. Open `01_India_TB_v20_simulation_trend.ipynb` (§3 edits are already applied).
2. Smoke test: `SCENARIO='base'`, 2 replicates. Check that
   incidence declines over 2000–2020 for a design point with `tf_gain` near the
   top of its prior, and is flat for `tf_gain = 0`.
3. Wave 1: `SCENARIO='wave1_cal'`, `CAL_DESIGN_CSV='helpers/wave1_design_trend.csv'`,
   200 design points x 10 replicates = 2,000 sims.
4. `python helpers/compute_targets_trend.py <wave1_uid>` → burden emulator, NROY,
   wave-2 design.
5. Wave 2: 100 points in the NROY x 25 replicates = 2,500 sims, then the same
   target script.
6. Production: `SCENARIO='combined'`, 8 scenarios x 100 replicates, with the
   best-fit parameters as the subsample dict.
7. Regenerate everything downstream with the existing 2026-09 scripts, which all
   read the experiment UID: `fig1_rebuild_0909.py`, `fig1_plot_0909.py`,
   `scenarios_rebuild_0909.py`, `fig3_rebuild_0909.py`, `fig4_incidence_ci_0909.py`,
   `fig4_plot_0909.py`, `fig5_rebuild_0909.py`.

## 6. What to expect

- **Every downstream number changes.** The loop-group distribution, delay
  medians, Table 2, the abstract percentages and all four rebuilt figures are
  functions of the calibrated parameters. The 2026-09 recomputed values
  (26/39/4/31, 2.4 → 17.4 months, 31% → 62% of transmission, all-three 58%)
  should be treated as the *current* baseline and re-derived after the re-run.
- **The cascade validation targets will move, probably away from Subbaraman.** A
  higher treatment fraction means fewer pre-treatment defaults, so
  `tests_per_case` and the loop-group fractions will fall. They are validation
  targets, so this is reportable rather than fatal — but if the trend can only be
  reproduced at the cost of the cascade structure, that tension is itself a
  finding and belongs in the Discussion.
- **The improvement mechanism may be insufficient.** Reproducing −2.4%/yr from a
  cascade-only driver is demanding. If wave 1 shows the NROY empty or the best
  fit still far from the trend, the next lever to add is a declining contact rate
  / transmission over time. EMOD's `Infectivity_Scale_Type` cannot do this here —
  it only supports exponential *growth* and only for `GENERIC_SIM` /
  `ENVIRONMENTAL_SIM`, not `TBHIV_SIM` — so it would have to be a stepped
  `Base_Infectivity` change via config over sub-epochs, which requires a
  different sweep structure.


---

# Wave 1 results — experiment `fe725f4e-39b8-f111-92f1-000d3af5294c` (2026-09-24)

200 design points x 10 replicates. Processed with `helpers/compute_targets_trend.py`.

## The re-calibration works

Best fit is design 136, **max I = 1.48 across all six burden targets**, with an
implied 2000-2020 incidence trend of **-2.93%/yr** against WHO's -2.43%/yr. The
published model ran at +0.05%/yr, so the direction is fixed.

| target | model | target | I |
|---|---|---|---|
| incidence 2000 | 403 | 322 +/- 55 | 1.48 |
| incidence 2010 | 304 | 276 +/- 47 | 0.60 |
| incidence 2015 | 276 | 237 +/- 40 | 0.98 |
| incidence 2020 | 222 | 197 +/- 34 | 0.76 |
| prevalence 2020 | 281 | 304 +/- 36 | 0.64 |
| mortality 2020 | 27.4 | 33 +/- 6 | 0.94 |

Held-out years: 2023 model 195 vs WHO 192; 2005 model 365 vs WHO 311. The model
starts a little high and declines a little fast, landing accurately at the end.

All seven cascade validation targets sit within 3 sigma (max 2.11). The loop-group
composition shifts relative to the published run: 0.31 / 0.29 / 0.40 (no loop /
diagnostic loop / compound) versus 0.26 / 0.39 / 0.35.

## The trend is driven by care-seeking, not treatment initiation

Rank correlation with the 2000-2020 trend across the design:

| parameter | rho |
|---|---|
| **cs_gain_per_year** | **-0.42** |
| Post_Infection_Acquisition_Multiplier | -0.40 |
| Base_Infectivity_Constant | -0.31 |
| tf_gain_per_year | -0.07 |

By quartile of `cs_gain_per_year` the median trend is -1.26, -1.92, -2.51,
-2.93%/yr (monotone); by quartile of `tf_gain_per_year` it is flat
(-2.04, -2.32, -2.18, -1.85).

Reading: removing the 2012 downward step in the treatment fraction is what
removed the rising baseline, but the *slope* of the historical decline is set by
shortening patient delay. Raising treatment initiation further acts only on
patients already diagnosed and barely moves incidence. This belongs in the
Discussion — it is a cascade-mechanism result, not a technical detail.

## Wave-2 design — `helpers/wave2_design_trend.py`

The naive "refined box" recipe does not apply here: only 7% of wave-1 points
(14/200) pass, yet the per-parameter marginal ranges shrink to just 77-99% of the
prior. The NROY is a thin joint region, not a box.

A first emulator attempt on the raw outputs also failed — 2020 incidence spans
0 to 11,489 per 100,000 across the design (77 collapsed, 67 exploded), so the GP
was dominated by the exploded points and returned a 77% NROY.

The working recipe is two-stage:

1. **Viability gate** — random-forest classifier for "2020 incidence in
   [50, 800]", 5-fold CV AUC 0.93; keeps 24.6% of the prior.
2. **Log-scale emulation** — one GP per target on log10(output), fitted to the
   viable subset, implausibility evaluated in log space. CV R2 0.75-0.85
   (raw scale gave 0.52-0.82).

Emulator NROY: **6.6% of the prior**, which matches the 7% actual wave-1 pass
rate — the check that the first attempt failed.

Validation against the wave-1 simulations: **all 14 actual NROY points are
retained** (no known-good region is excluded), 177/186 non-NROY points are
correctly rejected, and best-fit design 136 is retained.

Output: `helpers/wave2_design_trend.csv` — 100 points chosen as k-means
centroids of the accepted candidates (median emulator implausibility 1.96).

## Running wave 2

In `01_India_TB_v20_simulation_trend.ipynb`, cell 5:

    SCENARIO       = 'wave2_cal'
    CAL_DESIGN_CSV = 'helpers/wave2_design_trend.csv'
    replicates     = 25

100 x 25 = 2,500 simulations. 25 replicates (not 10) because the earlier waves
showed that is what lifts the burden emulators past R2 = 0.9 and cuts the
cascade Monte-Carlo noise.

Then:

    python helpers/compute_targets_trend.py <wave2_uid>

and the combined wave-1 + wave-2 fit selects the production parameter set.
