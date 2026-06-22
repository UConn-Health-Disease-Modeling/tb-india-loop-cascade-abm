"""
recalibration_plan.py
=====================

Recalibration plan for the India TB EMOD model so that **raw** simulator
outputs (no post-hoc loop-group resampling, no calibrated delay override)
land near empirical India cascade-of-care data.

Why this file exists
--------------------
Running `delay_decomposition_raw.py` against the current simulator output
exposes two problems that the old override formula was masking:

1. Loop-group proportions are off (engagement-confounded):
       Current raw : NoLoop 52%, TestLoop 34%, TestTxfail 9%
       Subbaraman  : NoLoop 50%, TestLoop 35%, TestTxfail 13%, TxFailonly 2%
   BUT >75% of the current "NoLooponly" pool has health_system_delay = 0
   (median IQR 0-0). Those are NEVER-ENGAGED patients (no test events at
   all), not true single-pass successes. After excluding them, the
   engaged-cohort distribution is much further from target than the raw
   numbers suggest. Use `reclassify_loop_group_engagement_aware()` below.

2. Delay phases are longer than plausible empirical range for India.
   IMPORTANT: published empirical delays are systematically biased low
   (recall bias, survivor bias, only counts diagnosed cases). A realistic
   plausible range for the simulator is 60-120 d per phase, with up to
   180 d acceptable for patient/HS delay in low-quality-care subsets.

                       OBSERVED raw   Plausible range   Mid-target   Ratio
                                       (India)                       raw/target
       Subclinical    : 181 d         60-90 d            60 d        3.0x
       Patient delay  : 262 d         60-180 d           90 d        2.9x
       HS delay (mean): 280 d         60-180 d           90 d        3.1x
       Total          : ~700 d        180-390 d         240 d        2.9x

   By loop_group (raw means):
       NoLooponly      subclinical=181  patient=262  HS=18    total=461
       TestLooponly    subclinical=181  patient=258  HS=429   total=867
       TestTxfail      subclinical=178  patient=262  HS=1085  total=1525

This module:
- Documents which EMOD knob controls which delay phase / cascade outcome.
- Provides target values grounded in the published India cascade literature.
- Exposes `patch_config()` and `patch_campaign_kwargs()` so the next COMPS
  calibration run can start from these revised priors.
- Provides `reclassify_loop_group_engagement_aware()` for the analysis side,
  so that "never-engaged" patients are separated from true NoLoop.

This does NOT run a calibration. Calibration is a COMPS / CALCULON job.
This file gives you the starting parameter set for that job and the
analysis-side classifier to evaluate its output.

References
----------
- Subbaraman R et al. 2016. The TB cascade of care in India's public sector.
  PLoS Med 13(10): e1002149.
- Sreeramareddy CT et al. 2014. Time delays in diagnosis of pulmonary TB:
  a systematic review of literature. BMC Infect Dis 14:194 (India studies).
- Cazabon D et al. 2017. Quality of TB care in India: a systematic review.
  Int J Tuberc Lung Dis 22(8):816-826.
- Frascella B et al. 2021. Subclinical TB: review of prevalence and natural
  history. Clin Infect Dis 73(3): e830-e841.
"""

from __future__ import annotations

import numpy as np
import pandas as pd


# =============================================================================
# 1. Empirical targets (India)
# =============================================================================
EMPIRICAL_TARGETS = {
    # Delay phases: target MEAN days per episode.
    # Targets are set to the mid-upper end of the plausible India range,
    # NOT to the published empirical means -- those are biased low
    # (recall, survivor, diagnosed-only artifacts). Acceptable range:
    #   subclinical 60-90 d, patient 60-180 d, HS 60-180 d.
    "subclinical_mean_days":       60.0,   # Frascella 2021 lower bound; OK up to 90
    "patient_delay_mean_days":     90.0,   # mid of plausible 60-120 (180 OK in LQ care)
    "hs_delay_mean_days":          90.0,   # mid of plausible 60-120 (180 OK with loops)
    "total_delay_mean_days":      240.0,   # sum of the three above

    # Loop-group fractions among ACTIVE-TB engaged patients (Subbaraman 2016)
    "frac_NoLoop":      0.50,   # single-pass: 1 pos test + successful tx
    "frac_TestLoop":    0.35,   # >=1 test failure or pre-tx default, but eventual cure
    "frac_TestTxfail":  0.13,   # combined test loop + treatment failure
    "frac_TxFailonly":  0.02,   # rare: no test loop, but tx failed

    # Cascade outcome rates (per-event probabilities)
    "p_first_test_positive":   0.78,   # = test sensitivity weighted by smear-pos/neg mix
    "p_start_tx_given_pos":    0.85,   # treatment_fraction (current = 0.60, too low)
    "p_tx_failure_first_line": 0.15,   # NTEP failure rate first-line first attempt

    # Per-event delay intervals
    "days_per_retest":         45.0,   # patient who tested neg/defaulted, retest interval
    "days_per_tx_course":     180.0,   # 6-month first-line regimen
}


# =============================================================================
# 2. EMOD parameter recommendations
# =============================================================================
# Each entry: (param_name_in_emod, current_value, recommended_value,
#              affects_phase, rationale)
PARAMETER_RECOMMENDATIONS: list[dict] = [
    # ── Subclinical period ─────────────────────────────────────────────
    {
        "param":      "TB Presymptomatic Duration Years  (params.py sample dict)",
        "current":    0.5,            # OBSERVED produces mean 181 d
        "recommended": 0.165,         # 60 d
        "phase":      "subclinical",
        "rationale":  "Observed mean = 181 d; target 60 d (Frascella 2021). "
                      "Reduce calibration guess from 0.5 yr (182 d) -> 0.165 yr (60 d). "
                      "The Min should also drop from 10/365 to 30/365 so the optimizer "
                      "doesn't push it absurdly low.",
    },

    # ── Care-seeking (patient delay) ────────────────────────────────────
    # Observed mean = 262 d; target = 90 d (mid of plausible 60-120 range).
    # Required multiplicative boost on HSB probabilities = 262/90 = 2.9x.
    # MUCH milder than the 8.7x implied by literature-mean targets, because
    # published empirical delays are biased low (recall, survivor, diagnosed-
    # only) -- the model should land in the upper end of plausible.
    {
        "param":      "probability_per_step_slow  (HSB)",
        "current":    0.0208,         # OBSERVED produces ~262 d mean patient delay
        "recommended": 0.06,          # 0.0208 * 2.9 -- target ~90 d mean
        "phase":      "patient_delay",
        "rationale":  "Boost slow rate ~2.9x so OBSERVED mean patient delay shrinks "
                      "from 262 d -> ~90 d (mid of plausible 60-120 range; "
                      "up to 180 d acceptable for low-quality-care subset). "
                      "Validate after next sim run; iterate if needed.",
    },
    {
        "param":      "probability_per_step_fast  (HSB)",
        "current":    0.05247,        # currently 2x slow
        "recommended": 0.15,          # keep 2-2.5x slow ratio
        "phase":      "patient_delay",
        "rationale":  "Boost fast rate by same ~2.9x factor so high-quality-care "
                      "patients land near ~30-45 d (proper fast-seeker tail).",
    },
    {
        "param":      "fixed_initial_delay  (HSB)",
        "current":    7,
        "recommended": 3,
        "phase":      "patient_delay",
        "rationale":  "Median time from symptom onset to recognition of symptoms "
                      "as needing care is closer to 2-4 days in survey data.",
    },

    # ── Test sensitivity (drives false-negative loop) ──────────────────
    {
        "param":      "sens_smear_pos_GXH  (params.py, GeneXpert era, high quality)",
        "current":    0.75,
        "recommended": 0.88,
        "phase":      "hs_delay / loop_group",
        "rationale":  "GeneXpert real-world sensitivity for smear-pos pulmonary TB is "
                      "~88% (WHO 2020 Xpert handbook). 75% under-detects, inflates "
                      "TBTestNegative loop. Set to ~0.85-0.90.",
    },
    {
        "param":      "sens_smear_neg_GXH  (params.py, GeneXpert era, high quality)",
        "current":    0.65,
        "recommended": 0.75,
        "phase":      "hs_delay / loop_group",
        "rationale":  "Xpert sens for smear-neg pulmonary TB ~70-75%. 65% is "
                      "conservative-low.",
    },
    {
        "param":      "sens_smear_pos_GXL  (params.py, GeneXpert era, low quality)",
        "current":    0.65,
        "recommended": 0.78,
        "phase":      "hs_delay / loop_group",
        "rationale":  "Private sector with X-ray+clinical dx: roughly 75-80% sens "
                      "for smear-pos in India private-sector audits.",
    },
    {
        "param":      "sens_smear_neg_GXL  (params.py, GeneXpert era, low quality)",
        "current":    0.60,
        "recommended": 0.65,
        "phase":      "hs_delay / loop_group",
        "rationale":  "Small bump to keep public > private gradient.",
    },

    # ── Pre-treatment loss (drives TBTestDefault loop) ─────────────────
    {
        "param":      "treatment_fraction  (DiagnosticTreatNeg)",
        "current":    0.60,
        "recommended": 0.85,
        "phase":      "hs_delay / loop_group",
        "rationale":  "Subbaraman 2016: ~85% of diagnosed start treatment in NTEP. "
                      "0.60 implies 40% pre-tx loss, which over-fires TBTestDefault "
                      "and is one of the main drivers of inflated TestLoop fraction.",
    },

    # ── Treatment failure ──────────────────────────────────────────────
    {
        "param":      "TB_Drug_Efficacy_Multiplier_Failed",
        "current":    1.0,
        "recommended": 1.0,
        "phase":      "hs_delay / loop_group",
        "rationale":  "Keep at 1.0; underlying first-line failure rate should target "
                      "~15% per course. Tune drug.json regimen cure_rate "
                      "(currently 0.55 implied for v10 by header comment in "
                      "set_config.py line 81 'tf=0.55').",
    },
    {
        "param":      "TB_Drug_Cure_Rate  (drug.json: first-line)",
        "current":    "tf=0.55 implied",
        "recommended": "tf=0.85 cure on first course",
        "phase":      "hs_delay / loop_group",
        "rationale":  "v10 comment implies 55% cure rate; NTEP first-line cohort cure "
                      "is ~85% (TB India report 2022). Bumping cure rate reduces "
                      "TBFailedDrugRegimen events and shrinks TestTxfail group.",
    },

    # ── Active TB duration (controls untreated trajectory length) ─────
    {
        "param":      "set_duration_no_hiv(duration=...)",
        "current":    2.5,            # years
        "recommended": 2.0,           # years
        "phase":      "hs_delay (tail)",
        "rationale":  "Untreated active TB natural duration ~2 yr (Tiemersma 2011). "
                      "2.5 yr is on the high side and lengthens the HS-delay tail for "
                      "patients who fail multiple cascade attempts.",
    },
]


# =============================================================================
# 3. Apply recommendations to a live EMOD config object
# =============================================================================
def patch_config(config) -> None:
    """
    Apply the recommended values to an EMOD `config.parameters` namespace.

    Call this from your `set_param_fn` (in 01_India_TB_v5_full_cascade.py)
    AFTER `set_config.set_config(config)` so it overrides the defaults.

        from recalibration_plan import patch_config
        def set_param_fn(config):
            config = set_config.set_config(config)
            patch_config(config)        # <-- apply recalibrated values
            ...
    """
    p = config.parameters
    # Subclinical: observed 181 d -> target 60 d
    p.TB_Presymptomatic_Rate = 1.0 / 60.0
    # Active TB duration: keep death_fraction = 0.7, duration = 2.0 yr
    duration_days = 2.0 * 365.0
    death_frac = 0.7
    p.TB_Active_Cure_Rate     = (1.0 - death_frac) / duration_days
    p.TB_Active_Mortality_Rate = death_frac / (1.0 - death_frac) * p.TB_Active_Cure_Rate
    # Drug efficacy multipliers stay 1.0; cure rate is tuned in drug.json directly.


def patch_campaign_kwargs() -> dict:
    """
    Return the recommended kwargs for `build_camp()` in the campaign builder.

    Multipliers are calibrated to bring OBSERVED mean patient delay from
    262 d -> ~30 d (Subbaraman 2016 target). Iterate after running.

        from recalibration_plan import patch_campaign_kwargs
        camp_kwargs = patch_campaign_kwargs()
        build_camp(**camp_kwargs)
    """
    return {
        "probability_per_step_slow": 0.06,   # was 0.0208; bumped 2.9x  -> target ~90d
        "probability_per_step_fast": 0.15,   # was 0.05247; bumped 2.9x -> target ~30-45d
        # other kwargs (acf_coverage, personal_hc_budget, cascade_efficiency)
        # keep at user's existing scenario defaults.
    }


# =============================================================================
# 4b. Empirically-driven adjustment suggester
# =============================================================================
def suggest_adjustments_from_observed(delay_df: pd.DataFrame,
                                      base_slow_prob: float = 0.0208,
                                      base_fast_prob: float = 0.05247,
                                      base_presymp_years: float = 0.5) -> pd.DataFrame:
    """
    Given the actual raw delay_df from a sim run, compute the multiplicative
    parameter adjustments needed to land at empirical India targets.

    Reports the observed mean for each phase, the target, the ratio, and
    the recommended new parameter value. Caps probabilities below 0.5.

    Parameters
    ----------
    delay_df : output of compute_delay_decomposition(), filtered to engaged
               patients if you want -- otherwise the never-engaged tail will
               distort the patient_delay mean.
    base_*   : the current parameter values (defaults match your `params.py`
               / campaign builder).
    """
    df = reclassify_loop_group_engagement_aware(delay_df)
    engaged = df[df["loop_group_engaged"] != "NeverEngaged"]
    if len(engaged) == 0:
        print("No engaged episodes; cannot suggest adjustments.")
        return pd.DataFrame()

    obs_subclin = engaged["subclinical_days"].mean()
    obs_patient = engaged["patient_delay_days"].mean()
    obs_hs      = engaged["health_system_delay_days"].mean()

    tgt_subclin = EMPIRICAL_TARGETS["subclinical_mean_days"]
    tgt_patient = EMPIRICAL_TARGETS["patient_delay_mean_days"]
    tgt_hs      = EMPIRICAL_TARGETS["hs_delay_mean_days"]

    ratio_subclin = obs_subclin / tgt_subclin
    ratio_patient = obs_patient / tgt_patient
    ratio_hs      = obs_hs      / tgt_hs

    # Subclinical: rate is inversely proportional to mean, so duration_years
    # scales by target/observed
    new_presymp_years = base_presymp_years * (tgt_subclin / obs_subclin)

    # Patient delay: probability_per_step is inversely proportional to mean,
    # so probability scales by ratio (multiplied by observed/target)
    new_slow = min(0.49, base_slow_prob * ratio_patient)
    new_fast = min(0.49, base_fast_prob * ratio_patient)

    rows = [
        {"phase": "subclinical",
         "observed_mean_d": obs_subclin,  "target_mean_d": tgt_subclin,
         "obs/target": ratio_subclin,
         "param": "TB Presymptomatic Duration Years",
         "current_value": base_presymp_years,
         "recommended_value": round(new_presymp_years, 3)},
        {"phase": "patient_delay (slow)",
         "observed_mean_d": obs_patient,  "target_mean_d": tgt_patient,
         "obs/target": ratio_patient,
         "param": "probability_per_step_slow",
         "current_value": base_slow_prob,
         "recommended_value": round(new_slow, 4)},
        {"phase": "patient_delay (fast)",
         "observed_mean_d": obs_patient,  "target_mean_d": tgt_patient,
         "obs/target": ratio_patient,
         "param": "probability_per_step_fast",
         "current_value": base_fast_prob,
         "recommended_value": round(new_fast, 4)},
        {"phase": "hs_delay",
         "observed_mean_d": obs_hs,       "target_mean_d": tgt_hs,
         "obs/target": ratio_hs,
         "param": "(driven by test sens / treatment_fraction / tx cure)",
         "current_value": "see PARAMETER_RECOMMENDATIONS",
         "recommended_value": "see PARAMETER_RECOMMENDATIONS"},
    ]
    out = pd.DataFrame(rows)
    print("\n=== Empirically-driven recalibration suggestions ===")
    print(out.to_string(index=False))
    return out


def patch_params_module_attrs() -> dict:
    """
    Return the recommended attribute overrides for `params.py`. Set them
    by:

        import params
        from recalibration_plan import patch_params_module_attrs
        for k, v in patch_params_module_attrs().items():
            setattr(params, k, v)

    BEFORE the campaign-builder import path picks up `params.sens_smear_*`.
    """
    return {
        "sens_smear_pos_GXH": 0.88,
        "sens_smear_neg_GXH": 0.75,
        "sens_smear_pos_GXL": 0.78,
        "sens_smear_neg_GXL": 0.65,
        # Pre-DOTS sensitivities are usually left alone for historical fidelity,
        # but bump them slightly if the burn-in delay tail is also too long.
    }


def diagnostic_treat_neg_kwargs() -> dict:
    """
    Return the recommended kwargs for `DiagnosticTreatNeg()` calls in the
    campaign builder so pre-treatment loss matches NTEP (~15%).
    """
    return {"treatment_fraction": 0.85}


# =============================================================================
# 4. Analysis-side: engagement-aware loop-group classifier
# =============================================================================
def reclassify_loop_group_engagement_aware(delay_df: pd.DataFrame) -> pd.DataFrame:
    """
    Replace the simple loop_group label with a 5-level taxonomy that
    distinguishes patients who never engaged the cascade from true
    single-pass successes.

    Categories (computed from per-episode event counts):
        NeverEngaged   : never had a test event (n_tests_total == 0 inferred
                         from n_testneg == 0 AND n_testdef == 0 AND
                         health_system_delay_days == 0). These should NOT
                         be counted toward NoLoop in cascade reporting.
        NoLooponly     : true one-pass success (had a test, no failures).
        TestLooponly   : >=1 test failure (TBTestNegative or TBTestDefault),
                         no treatment failure.
        TxFailonly     : no test failure, >=1 treatment failure.
        TestTxfail     : both.

    Requires that `delay_df` have columns from `compute_delay_decomposition()`:
    n_testneg, n_testdef, n_txfail, health_system_delay_days.
    """
    df = delay_df.copy()
    # If the input is empty, return an empty DF with all expected columns so
    # downstream code (which assumes the schema) works without special-casing.
    if len(df) == 0:
        empty_cols = list(df.columns) + [
            'n_testneg', 'n_testdef', 'n_txfail',
            'subclinical_days', 'patient_delay_days', 'health_system_delay_days',
            'total_days', 'loop_group', 'loop_group_engaged', 'replicate',
            'Individual_ID',
        ]
        out = pd.DataFrame(columns=list(dict.fromkeys(empty_cols)))
        return out
    required = {"n_testneg", "n_testdef", "n_txfail", "health_system_delay_days"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"delay_df missing columns: {sorted(missing)}")

    n_tn = df["n_testneg"].fillna(0).astype(int)
    n_td = df["n_testdef"].fillna(0).astype(int)
    n_tf = df["n_txfail"].fillna(0).astype(int)
    hs   = df["health_system_delay_days"].fillna(0.0)

    # NeverEngaged = no test events AND no HS delay accrued
    never_engaged = (n_tn == 0) & (n_td == 0) & (n_tf == 0) & (hs <= 0.0)
    has_test_loop = (n_tn > 0) | (n_td > 0)
    has_tx_loop   = (n_tf > 0)

    cat = np.where(
        never_engaged, "NeverEngaged",
        np.where(has_test_loop & has_tx_loop, "TestTxfail",
        np.where(has_tx_loop,                 "TxFailonly",
        np.where(has_test_loop,               "TestLooponly",
                                              "NoLooponly"))))
    df["loop_group_engaged"] = cat
    return df


def summarize_target_vs_observed(delay_df: pd.DataFrame) -> pd.DataFrame:
    """
    Build a side-by-side table of current raw vs empirical-target values
    for the headline metrics, restricted to ENGAGED patients only.
    """
    df = reclassify_loop_group_engagement_aware(delay_df)
    engaged = df[df["loop_group_engaged"] != "NeverEngaged"]

    if len(engaged) == 0:
        print("No engaged patients in delay_df.")
        return pd.DataFrame()

    rows = []

    # Delay phase means
    rows.append({"metric": "Subclinical mean (d)",
                 "current": engaged["subclinical_days"].mean(),
                 "target":  EMPIRICAL_TARGETS["subclinical_mean_days"]})
    rows.append({"metric": "Patient delay mean (d)",
                 "current": engaged["patient_delay_days"].mean(),
                 "target":  EMPIRICAL_TARGETS["patient_delay_mean_days"]})
    rows.append({"metric": "HS delay mean (d)",
                 "current": engaged["health_system_delay_days"].mean(),
                 "target":  EMPIRICAL_TARGETS["hs_delay_mean_days"]})
    rows.append({"metric": "Total delay mean (d)",
                 "current": engaged["total_days"].mean(),
                 "target":  EMPIRICAL_TARGETS["total_delay_mean_days"]})

    # Loop-group fractions (engaged base)
    n_eng = len(engaged)
    for grp, key in [("NoLooponly",   "frac_NoLoop"),
                     ("TestLooponly", "frac_TestLoop"),
                     ("TestTxfail",   "frac_TestTxfail"),
                     ("TxFailonly",   "frac_TxFailonly")]:
        rows.append({"metric": f"Fraction {grp} (engaged)",
                     "current": (engaged["loop_group_engaged"] == grp).sum() / n_eng,
                     "target":  EMPIRICAL_TARGETS[key]})

    out = pd.DataFrame(rows)
    out["delta"]    = out["current"] - out["target"]
    out["pct_off"]  = 100.0 * (out["current"] - out["target"]) / out["target"]
    print(f"\n=== Current (RAW, engaged-only) vs empirical India targets "
          f"(n_engaged={n_eng:,}, n_neverengaged="
          f"{int((df['loop_group_engaged'] == 'NeverEngaged').sum()):,}) ===")
    print(out.to_string(index=False,
                        float_format=lambda x: f"{x:8.2f}" if abs(x) > 1
                                                else f"{x:8.4f}"))
    return out


# =============================================================================
# 5. Pretty-print recommendations
# =============================================================================
def print_recommendations() -> None:
    """Print the full parameter-recommendation table."""
    print("=" * 78)
    print("EMOD PARAMETER RECALIBRATION RECOMMENDATIONS")
    print("=" * 78)
    by_phase: dict[str, list[dict]] = {}
    for r in PARAMETER_RECOMMENDATIONS:
        by_phase.setdefault(r["phase"], []).append(r)

    for phase, recs in by_phase.items():
        print(f"\n--- {phase} ---")
        for r in recs:
            print(f"  param      : {r['param']}")
            print(f"  current    : {r['current']}")
            print(f"  recommended: {r['recommended']}")
            print(f"  rationale  : {r['rationale']}")
            print()


if __name__ == "__main__":
    print(__doc__)
    print()
    print_recommendations()
