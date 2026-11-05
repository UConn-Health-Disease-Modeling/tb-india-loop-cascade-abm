"""Time-varying cascade schedules used to build the EMOD campaign.

WHY THIS EXISTS
---------------
An earlier version of the model held cascade performance constant after 2012
and, worse, stepped the treatment-initiation fraction DOWN from 0.70 (DOTS era) to 0.55 at
GeneXpert introduction -- a value that was hand-tuned to hit the tests-per-case
cascade target.  The consequence was a simulated incidence trajectory that was flat
over 2000-2020 (+0.05%/yr) and rising thereafter (+0.59%/yr), whereas WHO
estimates India's incidence fell from 322 to 197 per 100,000 over 2000-2020
(-2.4%/yr).

This module replaces the two hand-set epochs with a single monotone cascade
improvement schedule from DOTS introduction (1997) to the present (2024), driven
by two new calibrated parameters:

  tf_gain_per_year  absolute annual increase in the treatment-initiation
                    fraction, applied from 1997; the 1997 value stays at the
                    DOTS-era 0.70 and the value is capped at TF_MAX.
  cs_gain_per_year  relative annual increase in the care-seeking probability
                    per step, applied from 1997 (shortens patient delay).

Both are zero-inclusive, so the previous "no secular improvement" behaviour
remains inside the prior and the calibration can reject the mechanism if the
data do not support it.

The GeneXpert sensitivity ramp is preserved exactly: the schedule interpolates
the pre-GeneXpert -> GeneXpert sensitivities over params.length_gene_xpert_ramp
and hands each sub-epoch its own start/end sensitivities, so splitting the era
into sub-epochs does not restart the ramp.

USAGE (in the simulation notebooks, `build_camp`)
-------------------------------------------------
The diagnostic blocks from DOTS introduction onward are built as:

    for blk in cascade_schedule(params, tf_gain_per_year, cs_gain_per_year,
                                dx_boost=dx_boost, tx_boost=tx_boost,
                                sens_pos_GXL=_base_sens_pos):
        if blk.use_ramp:           # GeneXpert era -- must stay RampDTN
            camp.add(ramp.RampDTN(..., ramp_time=blk.duration,
                base_sensitivity_smearpos=blk.sens_pos_end,
                base_sensitivity_smearpos2=blk.sens_pos_start, ...))
        else:
            camp.add(dtn.DiagnosticTreatNeg(..., 
                base_sensitivity_smearpos=blk.sens_pos_end, ...))

(the treatment and resistance-diagnostic blocks are unchanged; they already run
with duration=-1 from their epoch start)

The care-seeking multiplier is applied separately, because the HSB blocks are
built once:

    for blk in careseeking_schedule(params, cs_gain_per_year):
        camp.add(hsb.HSB(..., probability_per_step=adj_slow * blk.multiplier,
                         start_day=blk.start_day, duration=blk.duration, ...))
"""
from __future__ import annotations
from dataclasses import dataclass

YEAR = 365.0
TF_MAX = 0.90              # ceiling on treatment-initiation fraction
TF_DOTS_1997 = 0.70        # unchanged DOTS-era starting value
STEP_YEARS = 3.0           # sub-epoch length for the stepped ramp
IMPROVEMENT_END_YEAR = 2024.0   # improvement stops at the present; 2025+ is the
                                # scenario baseline and is held at its 2024 value
TX_BOOST_LOSS_REDUCTION = 0.5   # treatment-support scenario halves the residual
                                # pre-treatment loss (relative, not a fixed 0.85)


@dataclass
class CascadeBlock:
    """One diagnostic block covering [year_start, year_end).

    use_ramp says which intervention the notebook must instantiate.  The
    GeneXpert era MUST stay on ramp.RampDTN: it is what emits
    ProviderOrdersTBTest / ProviderTestNoR, which the resistance diagnostics are
    triggered on and which the cascade classification counts.  Substituting
    DiagnosticTreatNeg there would silently remove both.

    sens_*_start / sens_*_end are the sensitivities at the two ends of the block,
    so a RampDTN sub-epoch continues the overall GeneXpert ramp instead of
    restarting it (RampDTN's *2 arguments are the pre-ramp values).
    """
    qoc: str                 # 'High' or 'Low'
    trigger: str             # 'TBTestDOTSHigh' or 'TBTestDOTSLow'
    start_day: float
    duration: float
    treatment_fraction: float
    sens_pos_start: float
    sens_neg_start: float
    sens_pos_end: float
    sens_neg_end: float
    use_ramp: bool
    year_start: float
    year_end: float


@dataclass
class SeekingBlock:
    start_day: float
    duration: float
    multiplier: float
    year_start: float
    year_end: float


def _year_of(params, day: float) -> float:
    return params.base_year + day / YEAR


def _day_of(params, year: float) -> float:
    return (year - params.base_year) * YEAR


def _tf_at(year: float, tf_gain: float) -> float:
    """Treatment-initiation fraction at a calendar year."""
    y = min(year, IMPROVEMENT_END_YEAR)
    t = max(0.0, y - 1997.0)
    return min(TF_MAX, TF_DOTS_1997 + tf_gain * t)


def _sens_at(year: float, params, qoc: str, smear: str) -> float:
    """Sensitivity at a calendar year, interpolating the GeneXpert ramp."""
    pre = getattr(params, f"sens_smear_{smear}_pre_G{qoc[0]}")
    post = getattr(params, f"sens_smear_{smear}_GX{qoc[0]}")
    gx_start = _year_of(params, params.genexpert_introduction)
    gx_end = gx_start + params.length_gene_xpert_ramp / YEAR
    if year <= gx_start:
        return pre
    if year >= gx_end:
        return post
    f = (year - gx_start) / (gx_end - gx_start)
    return pre + f * (post - pre)


def cascade_schedule(params, tf_gain_per_year: float, cs_gain_per_year: float = 0.0,
                     dx_boost: bool = False, tx_boost: bool = False,
                     sens_pos_GXL: float | None = None,
                     loop_reduction: float = 0.0,
                     end_year: float = 2050.0):
    """Non-overlapping diagnostic/treatment blocks from DOTS start to end_year.

    dx_boost / tx_boost apply only from the intervention start date onward, so
    the historical trajectory is identical across scenarios.
    """
    blocks = []
    dots_year = _year_of(params, params.dots_start)
    interv_year = _year_of(params, params.intervention_day)

    edges = []
    y = dots_year
    while y < end_year - 1e-9:
        nxt = min(y + STEP_YEARS, end_year)
        # never let a block straddle the intervention start
        if y < interv_year < nxt:
            edges.append((y, interv_year))
            y = interv_year
            continue
        edges.append((y, nxt))
        y = nxt

    gx_start = _year_of(params, params.genexpert_introduction)

    def _pos_with_override(year, qoc):
        """Low-QoC GeneXpert sensitivity, honouring the calibration override."""
        if qoc != "Low" or sens_pos_GXL is None or year <= gx_start:
            return _sens_at(year, params, qoc, "pos")
        gx_end = gx_start + params.length_gene_xpert_ramp / YEAR
        pre = params.sens_smear_pos_pre_GL
        f = min(1.0, (year - gx_start) / (gx_end - gx_start))
        return pre + f * (sens_pos_GXL - pre)

    for y0, y1 in edges:
        mid = 0.5 * (y0 + y1)
        use_ramp = mid >= gx_start
        for qoc, trigger in (("High", "TBTestDOTSHigh"), ("Low", "TBTestDOTSLow")):
            s_pos0 = _pos_with_override(y0, qoc)
            s_pos1 = _pos_with_override(y1, qoc)
            s_neg0 = _sens_at(y0, params, qoc, "neg")
            s_neg1 = _sens_at(y1, params, qoc, "neg")
            tf = _tf_at(mid, tf_gain_per_year)

            if mid >= interv_year:
                # scenario levers switch on only from the intervention start
                if dx_boost:
                    s_pos0 = _sens_at(y0, params, "High", "pos")
                    s_pos1 = _sens_at(y1, params, "High", "pos")
                    s_neg0 = _sens_at(y0, params, "High", "neg")
                    s_neg1 = _sens_at(y1, params, "High", "neg")
                if loop_reduction > 0:
                    # remove a fraction of the residual loss at the two
                    # diagnostic-side loop drivers; the treatment-failure driver
                    # is handled by the drug definition (loop_drug)
                    _k = 1.0 - loop_reduction
                    s_pos0 = 1.0 - (1.0 - s_pos0) * _k
                    s_pos1 = 1.0 - (1.0 - s_pos1) * _k
                    s_neg0 = 1.0 - (1.0 - s_neg0) * _k
                    s_neg1 = 1.0 - (1.0 - s_neg1) * _k
                    tf = 1.0 - (1.0 - tf) * _k
                if tx_boost:
                    # Treatment support is defined RELATIVE to the calibrated
                    # 2024 baseline, not as a fixed 0.85.  With a calibrated
                    # secular improvement the 2024 baseline may already exceed
                    # 0.85, which would make a fixed target a no-op (or a
                    # worsening).  Halving the residual pre-treatment loss keeps
                    # the scenario meaningful whatever the calibration returns.
                    tf = 1.0 - (1.0 - tf) * (1.0 - TX_BOOST_LOSS_REDUCTION)

            blocks.append(CascadeBlock(
                qoc=qoc, trigger=trigger,
                start_day=_day_of(params, y0),
                duration=_day_of(params, y1) - _day_of(params, y0),
                treatment_fraction=round(tf, 4),
                sens_pos_start=round(s_pos0, 4), sens_neg_start=round(s_neg0, 4),
                sens_pos_end=round(s_pos1, 4), sens_neg_end=round(s_neg1, 4),
                use_ramp=use_ramp,
                year_start=y0, year_end=y1))
    return blocks


def careseeking_schedule(params, cs_gain_per_year: float, end_year: float = 2050.0):
    """Multiplier on the per-step care-seeking probability, by sub-epoch."""
    out = []
    dots_year = _year_of(params, params.dots_start)
    y = dots_year
    while y < end_year - 1e-9:
        nxt = min(y + STEP_YEARS, end_year)
        mid = min(0.5 * (y + nxt), IMPROVEMENT_END_YEAR)
        mult = (1.0 + cs_gain_per_year) ** max(0.0, mid - 1997.0)
        out.append(SeekingBlock(start_day=_day_of(params, y),
                                duration=_day_of(params, nxt) - _day_of(params, y),
                                multiplier=round(mult, 4),
                                year_start=y, year_end=nxt))
        y = nxt
    return out

@dataclass
class InfectivityPoint:
    year: float
    days_since_start: float
    multiplier: float


def infectivity_multiplier_map(params, decline_per_year: float,
                               start_year: float = 1997.0,
                               end_year: float = IMPROVEMENT_END_YEAR,
                               step_years: float = 3.0):
    """Multiplier_By_Duration map for a NodeInfectivityMult intervention.

    Represents the part of India's historical TB decline that is NOT produced by
    the care cascade -- improving living conditions, nutrition, housing density,
    demographic change.  Without such a lever the calibration has no choice but
    to attribute the entire observed decline to cascade improvement, which forces
    an implausibly efficient present-day cascade and removes the headroom that
    diagnostic and treatment interventions would otherwise act on.

    The multiplier falls geometrically at `decline_per_year` from `start_year`,
    then holds flat from `end_year` onward (the projection period inherits the
    level reached in 2024; it is a counterfactual, not a continuing trend).

    decline_per_year = 0 reproduces the previous behaviour exactly, so the prior
    contains the "cascade explains everything" hypothesis.
    """
    times, values, pts = [], [], []
    y = start_year
    while y <= end_year + 1e-9:
        mult = (1.0 - decline_per_year) ** (y - start_year)
        days = (y - start_year) * YEAR
        times.append(round(days, 1))
        values.append(round(mult, 5))
        pts.append(InfectivityPoint(y, days, mult))
        y += step_years
    # hold the final level out to the end of the simulation
    times.append(round((2050.0 - start_year) * YEAR, 1))
    values.append(values[-1])
    return {"Times": times, "Values": values}, pts


if __name__ == "__main__":
    import sys
    from pathlib import Path
    sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "simulation"))
    import params

    print("cascade schedule preview (tf_gain=0.006/yr, cs_gain=0.02/yr)\n")
    print(f"{'years':>16} {'QoC':>5} {'tf':>6} {'smear+ sens':>14} {'kind':>5}")
    for b in cascade_schedule(params, 0.006, 0.02):
        if b.year_start > 2032:
            break
        print(f"{b.year_start:7.1f}-{b.year_end:<7.1f} {b.qoc:>5} "
              f"{b.treatment_fraction:6.3f} {b.sens_pos_start:6.3f}->{b.sens_pos_end:<6.3f} "
              f"{'ramp' if b.use_ramp else 'dtn'}")
    print("\ncare-seeking multipliers")
    for s in careseeking_schedule(params, 0.02):
        if s.year_start > 2032:
            break
        print(f"{s.year_start:7.1f}-{s.year_end:<7.1f} x{s.multiplier}")


# ── loop-reduction drug naming ─────────────────────────────────────────
#  A cure proportion lives in the config, which carries no start day, so a
#  scenario that changes it has to distribute a differently named drug from the
#  intervention date (see the v27 note in the simulation notebook).  Both the
#  campaign builder and the sweep callback derive the name from here so they
#  cannot disagree.
BASE_CURE = {"HQ": 0.80, "LQ": 0.50}
#  The stock regimens also carry death 0.03, inactivation 0.10 and relapse 0.02.
#  add_tb_drug_type converts the SUM of the four proportions with -log(1 - sum),
#  so the cure proportion cannot simply be scaled up: what is actually reduced
#  is the residual failure mass, 1 - cure - these three, and the mass removed is
#  moved into cure.  The sum is held below 1 so the rate stays finite.
DRUG_OTHER = 0.03 + 0.10 + 0.02
DRUG_SUM_CAP = 0.995


def loop_drug(qoc: str, tx_boost: bool = False, loop_reduction: float = 0.0):
    """Return (drug_name, cure_proportion) for a care stream after 2025.

    qoc is "HQ" or "LQ".  With no treatment support and no loop reduction the
    stock drug is returned unchanged, so unaffected arms stay bit-identical.
    """
    cure = BASE_CURE[qoc]
    if tx_boost and qoc == "LQ":
        cure = BASE_CURE["HQ"]
    if loop_reduction > 0:
        fail = max(0.0, 1.0 - cure - DRUG_OTHER)
        cure = min(cure + fail * loop_reduction, DRUG_SUM_CAP - DRUG_OTHER)
    if abs(cure - BASE_CURE[qoc]) < 1e-9:
        return f"DOTS{qoc}", cure
    return f"DOTS{qoc}C{int(round(cure * 1000)):04d}", cure
