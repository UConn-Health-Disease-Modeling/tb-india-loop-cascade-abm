"""[2026-10-05 cohort revision]  The trajectory cohort is every active TB episode
ACTIVATING in calendar years 2015-2024 (2015.0 <= act_year < 2025.0) that resolves
within the simulation, and the loop group is assigned PER EPISODE from that
episode's own TBTestNegative / TBTestDefault / TBFailedDrugRegimen counts, as the
Methods and Table 2 define it.  The published 1850-2050 pool was 52.5% pre-1997
and 20.6% projected-future episodes, and the group had been assigned from each
individual's whole event history (97% agreement with the per-episode rule).
Panel A rates are per-episode column sums (treatment starts = first initiation
per episode), not recorder-wide event totals, which counted every event in the
run regardless of window.  The full per-episode table (all years, with act_year)
is dumped to analysis/derived/baseline_episodes.csv.gz for figure3_loop_group_burden.py.
One replicate (rep 40) fades out before 1960 and contributes no cohort episodes;
it is dropped from the per-replicate rate means (97 reps) but kept in the
population-year denominator of Figure 3D, as in every other figure.

[2026-09-30 round-2 production run — experiment d892bdf7-18bd-f111-92f1-000d3af5294c, loose-cascade calibration with the secular infectivity decline]

[2026-09-25 production re-run — experiment d892bdf7-18bd-f111-92f1-000d3af5294c, trend-calibrated parameters,
interventions from 2025]

Rebuild Figure 1A/1B for the 2026-09 revision (Sourya comments #59, #60).

#59  Fig 1A showed raw per-replicate event/person counts accumulated over the
     WHOLE 1850-2050 simulation in a 5%-sampled, 20x-scaled synthetic
     population -- uninterpretable, and it mixed a latent-infection stock with
     active-TB flows.  Rebuilt as percentages / rates on one denominator:
     active TB episodes (same definition as Table 2).
#60  Fig 1B was a hardcoded schematic with illustrative durations.  Rebuilt
     with the third (Test-loop only) trajectory and real model medians, with
     transition times annotated.

Episode definition (the one Table 2 uses): a
TBActivation followed by a resolution (SymptomaticCleared | DiseaseDeaths)
before the next TBActivation.  Loop group is assigned per episode from the
episode's counts of TBTestNegative / TBTestDefault / TBFailedDrugRegimen.

Writes analysis/derived/figure1_values.json
"""
from __future__ import annotations
import io, json, sqlite3, sys
from pathlib import Path
import numpy as np
import pandas as pd

if sys.stdout.encoding and sys.stdout.encoding.lower() != "utf-8":
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")

HERE = Path(__file__).resolve().parent
DATA = HERE / "derived"
ROOT = HERE.parent
EXP = ROOT / "d892bdf7-18bd-f111-92f1-000d3af5294c"

PRESYMP = "TBActivationPresymptomatic"
TESTS = {"TBTestDOTSLow", "TBTestDOTSHigh"}
RESOLUTION = {"SymptomaticCleared", "DiseaseDeaths"}
CASCADE = ({"TBActivation", PRESYMP} | TESTS | RESOLUTION |
           {"TBTestDefault", "TBTestNegative", "TBTestPositive",
            "TBStartDrugRegimen", "TBFailedDrugRegimen"})

GROUPS = ["NoLooponly", "TestLooponly", "TxFailonly", "TestTxfail"]
COHORT = (2015.0, 2025.0)          # activation in calendar years 2015-2024
DUMP = DATA / "baseline_episodes.csv.gz"


def process_rep(rf: Path):
    """Return (per-episode DataFrame, per-rep event totals dict)."""
    df = pd.read_csv(rf, usecols=["Time", "Event_Name", "Individual_ID"])

    n_infected = df.loc[df.Event_Name == "NewInfection", "Individual_ID"].nunique()
    n_activated = df.loc[df.Event_Name == PRESYMP, "Individual_ID"].nunique()
    ev_tot = df.Event_Name.value_counts().to_dict()

    d = df[df.Event_Name.isin(CASCADE)].sort_values(["Individual_ID", "Time"],
                                                    kind="mergesort")
    ind = d.Individual_ID.values
    ev = d.Event_Name.values
    t = d.Time.values.astype(float)

    # episode index: number of TBActivation events seen so far within the
    # individual, counting the row itself.  Rows carrying index k lie between
    # activation k and activation k+1, so the tests/resolution of episode k
    # carry index k and the PRESYMP of episode k+1 carries index k.
    is_act = (ev == "TBActivation")
    new_ind = np.empty(len(ind), bool)
    new_ind[0] = True
    new_ind[1:] = ind[1:] != ind[:-1]
    cum = np.cumsum(is_act)
    offset = np.where(new_ind, cum - is_act.astype(int), np.nan)
    offset = pd.Series(offset).ffill().values
    epi = (cum - offset).astype(int)

    w = pd.DataFrame({"ind": ind, "ev": ev, "t": t, "epi": epi})

    acts = w[w.ev == "TBActivation"][["ind", "t", "epi"]].rename(columns={"t": "t_act"})
    pres = (w[w.ev == PRESYMP].groupby(["ind", "epi"])["t"].max()
            .rename("t_presym").reset_index())
    pres["epi"] = pres["epi"] + 1
    firsttest = (w[w.ev.isin(TESTS)].groupby(["ind", "epi"])["t"].min()
                 .rename("t_test").reset_index())
    res = (w[w.ev.isin(RESOLUTION)].groupby(["ind", "epi"])["t"].min()
           .rename("t_res").reset_index())
    death = (w[w.ev == "DiseaseDeaths"].groupby(["ind", "epi"])["t"].min()
             .rename("t_death").reset_index())
    txstart = (w[w.ev == "TBStartDrugRegimen"].groupby(["ind", "epi"])["t"].min()
               .rename("t_tx").reset_index())
    lastfail = (w[w.ev == "TBFailedDrugRegimen"].groupby(["ind", "epi"])["t"].max()
                .rename("t_lastfail").reset_index())
    ecnt = (w[w.ev.isin(TESTS | {"TBTestNegative", "TBTestDefault",
                                 "TBFailedDrugRegimen"})]
            .assign(k=lambda x: x.ev.map(lambda v: "n_tests" if v in TESTS else
                                         {"TBTestNegative": "n_testneg",
                                          "TBTestDefault": "n_testdef",
                                          "TBFailedDrugRegimen": "n_txfail"}[v]))
            .groupby(["ind", "epi", "k"]).size().unstack(fill_value=0).reset_index())
    for c in ("n_tests", "n_testneg", "n_testdef", "n_txfail"):
        if c not in ecnt:
            ecnt[c] = 0

    e = (acts.merge(pres, on=["ind", "epi"], how="left")
             .merge(firsttest, on=["ind", "epi"], how="left")
             .merge(res, on=["ind", "epi"], how="left")
             .merge(death, on=["ind", "epi"], how="left")
             .merge(txstart, on=["ind", "epi"], how="left")
             .merge(lastfail, on=["ind", "epi"], how="left")
             .merge(ecnt, on=["ind", "epi"], how="left"))
    e = e[e.t_res.notna()]                      # episode must resolve (Table 2 rule)
    for c in ("n_tests", "n_testneg", "n_testdef", "n_txfail"):
        e[c] = e[c].fillna(0)

    e["subclin"] = np.where(e.t_presym.notna() & (e.t_presym < e.t_act),
                            e.t_act - e.t_presym, 0.0)
    has_t = e.t_test.notna() & (e.t_test <= e.t_res)
    e["patient"] = np.where(has_t, np.maximum(0.0, e.t_test - e.t_act),
                            np.maximum(0.0, e.t_res - e.t_act))
    e["hs"] = np.where(has_t, np.maximum(0.0, e.t_res - e.t_test), 0.0)
    e["tx_to_res"] = np.where(e.t_tx.notna(), np.maximum(0.0, e.t_res - e.t_tx), np.nan)
    e["test_to_tx"] = np.where(e.t_tx.notna() & has_t,
                               np.maximum(0.0, e.t_tx - e.t_test), np.nan)
    e["died"] = e.t_death.notna()
    e["ever_tested"] = has_t
    e["tx_started"] = e.t_tx.notna()            # first treatment initiation
    # loop group PER EPISODE (Methods / Table 2 definition)
    has_loop = (e.n_testneg >= 1) | (e.n_testdef >= 1)
    has_txf = e.n_txfail >= 1
    e["group"] = np.select([has_loop & has_txf, has_loop & ~has_txf, ~has_loop & has_txf],
                           ["TestTxfail", "TestLooponly", "TxFailonly"],
                           default="NoLooponly")
    e["act_year"] = 1850 + e.t_act / 365.0

    # Additive trajectory segments (sum = subclinical + [activation -> resolution]),
    # used to draw Figure 1B with real durations:
    #   subclin | patient (activation -> first test) | dxloop (first test -> Tx start)
    #   | txfail (Tx start -> last failed regimen) | ontx (-> resolution)
    has_tx = e.t_tx.notna() & has_t
    has_fail = e.t_lastfail.notna() & has_tx
    e["seg_subclin"] = e["subclin"]
    e["seg_patient"] = e["patient"]
    e["seg_dxloop"] = np.where(has_tx, e.t_tx - e.t_test,
                               np.where(has_t, e.t_res - e.t_test, 0.0))
    e["seg_txfail"] = np.where(has_fail, e.t_lastfail - e.t_tx, 0.0)
    last_tx_point = np.where(has_fail, e.t_lastfail, e.t_tx)
    e["seg_ontx"] = np.where(has_tx, e.t_res - last_tx_point, 0.0)
    for c in ("seg_dxloop", "seg_txfail", "seg_ontx"):
        e[c] = np.maximum(0.0, e[c].astype(float))
    # effective infectious days = duration x infectiousness weight
    # Weights are relative to an average UNTREATED ACTIVE case (= 1.0), which is
    # not the same as EMOD's smear-positive reference.  Adults are 65% smear-
    # positive (1.0), 10% extrapulmonary (0) and 25% smear-negative (0.34604),
    # so an average active case is 0.73651 as infectious as a smear-positive one.
    #   subclinical : TB_Active_Presymptomatic_Infectivity_Multiplier = 0.166
    #                 (scale factor on Base_Infectivity) / 0.73651 = 0.2254
    #   on treatment: TB_Reduced_Transmit = 1 is a 100% reduction, so 0.0
    W_SUBCLINICAL, W_ON_TREATMENT = 0.2254, 0.0
    e["eff_days"] = (W_SUBCLINICAL * e.seg_subclin + 1.0 * e.seg_patient
                     + 1.0 * e.seg_dxloop + 1.0 * e.seg_txfail
                     + W_ON_TREATMENT * e.seg_ontx)

    ev_tot.update({"_n_infected": n_infected, "_n_activated": n_activated})
    return e[["group", "act_year", "eff_days", "subclin", "patient", "hs",
              "tx_to_res", "test_to_tx", "ever_tested", "died", "tx_started",
              "n_tests", "n_testneg", "n_testdef", "n_txfail",
              "seg_subclin", "seg_patient", "seg_dxloop", "seg_txfail",
              "seg_ontx"]], ev_tot


def main():
    tags = pd.read_sql_query("SELECT SIM_ID, Scenario_Code, Run_Number FROM results",
                             sqlite3.connect(str(EXP / "results.db")))
    sims = (tags.loc[tags.Scenario_Code == 0].sort_values("Run_Number")
            .SIM_ID.astype(str).tolist())
    print(f"Baseline sims: {len(sims)}")

    epi_all, prog_rows = [], []
    for i, s in enumerate(sims):
        rf = EXP / s / "ReportEventRecorder.csv"
        if not rf.exists():
            continue
        e, ev = process_rep(rf)
        e["rep"] = i
        epi_all.append(e)
        prog_rows.append({"rep": i, "n_infected": ev["_n_infected"],
                          "n_activated": ev["_n_activated"]})
        if (i + 1) % 20 == 0:
            print(f"  {i+1}/{len(sims)} reps", flush=True)

    E_all = pd.concat(epi_all, ignore_index=True)
    E_all.to_csv(DUMP, index=False, compression="gzip")
    print(f"\nResolved episodes 1850-2050: {len(E_all):,} in {E_all.rep.nunique()} reps"
          f"  -> wrote {DUMP.name}")
    E = E_all[(E_all.act_year >= COHORT[0]) & (E_all.act_year < COHORT[1])].copy()
    print(f"Cohort {COHORT[0]:.0f}-{COHORT[1] - 1:.0f}: {len(E):,} episodes in "
          f"{E.rep.nunique()} reps (reps without cohort episodes: "
          f"{sorted(set(E_all.rep) - set(E.rep))})")

    # ---- Panel A: events per 100 cohort episodes, from per-episode columns ----
    E["self_clear"] = ~E.died
    acols = {"tests_any": "n_tests", "test_neg": "n_testneg",
             "test_default": "n_testdef", "tx_start": "tx_started",
             "tx_fail": "n_txfail", "self_clear": "self_clear", "tb_death": "died"}
    R = E.groupby("rep").agg(n_episodes=("group", "size"),
                             **{k: (v, "sum") for k, v in acols.items()})
    rate_cols = list(acols)
    rates = R[rate_cols].div(R.n_episodes, axis=0) * 100.0
    t_crit = 1.985  # t(0.975, df=96)
    panelA = pd.DataFrame({"mean": rates.mean(), "ci": rates.sem() * t_crit})

    # ---- Panel B: loop-group composition, % of cohort episodes ----
    fr = (E.groupby(["rep", "group"]).size().unstack(fill_value=0)
          .reindex(columns=GROUPS, fill_value=0))
    fr = fr.div(fr.sum(axis=1), axis=0) * 100.0
    panelB = pd.DataFrame({"mean": fr.mean().values, "ci": (fr.sem() * t_crit).values,
                           "q25": fr.quantile(0.25).values,
                           "q75": fr.quantile(0.75).values}, index=GROUPS)

    P = pd.DataFrame(prog_rows)
    prog = P.n_activated / P.n_infected * 100.0     # whole run, as before

    def med(s):
        return float(np.median(s)) if len(s) else float("nan")

    dl = {}
    for g in GROUPS + ["ALL"]:
        sub = E if g == "ALL" else E[E.group == g]
        dl[g] = {"n": int(len(sub)),
                 "subclin": med(sub.subclin), "patient": med(sub.patient),
                 "hs": med(sub.hs), "total": med(sub.subclin + sub.patient + sub.hs),
                 "test_to_tx": med(sub.test_to_tx.dropna()),
                 "tx_to_res": med(sub.tx_to_res.dropna()),
                 "pct_tested": float(sub.ever_tested.mean() * 100),
                 "pct_died": float(sub.died.mean() * 100),
                 # additive mean segments for the Fig 1B trajectories
                 "m_subclin": float(sub.seg_subclin.mean()),
                 "m_patient": float(sub.seg_patient.mean()),
                 "m_dxloop": float(sub.seg_dxloop.mean()),
                 "m_txfail": float(sub.seg_txfail.mean()),
                 "m_ontx": float(sub.seg_ontx.mean()),
                 "m_eff_days": float(sub.eff_days.mean()),
                 "med_eff_days": med(sub.eff_days),
                 "n_tests": float(sub.n_tests.mean()),
                 "n_testneg": float(sub.n_testneg.mean()),
                 "n_testdef": float(sub.n_testdef.mean()),
                 "n_txfail": float(sub.n_txfail.mean())}

    out = {"cohort": {"window": [int(COHORT[0]), int(COHORT[1]) - 1],
                      "rule": "2015.0 <= act_year < 2025.0 and resolves within the run",
                      "grouping": "per episode (n_testneg / n_testdef / n_txfail)"},
           "n_episodes": int(len(E)), "n_reps": int(E.rep.nunique()),
           "n_episodes_all_years": int(len(E_all)), "n_reps_all": int(E_all.rep.nunique()),
           "panelA_per100": {k: {"mean": float(panelA.loc[k, "mean"]),
                                 "ci": float(panelA.loc[k, "ci"])} for k in rate_cols},
           "panelB_pct": {g: {"mean": float(panelB.loc[g, "mean"]),
                              "ci": float(panelB.loc[g, "ci"]),
                              "iqr": [float(panelB.loc[g, "q25"]),
                                      float(panelB.loc[g, "q75"])]} for g in GROUPS},
           "progression_pct": {"mean": float(prog.mean()), "ci": float(prog.sem() * t_crit)},
           "group_delays_days": dl}
    (DATA / "figure1_values.json").write_text(json.dumps(out, indent=1), encoding="utf-8")

    print("\nPanel A (events per 100 cohort episodes):")
    print(panelA.round(1).to_string())
    print("\nPanel B (% of cohort episodes):")
    print(panelB.round(2).to_string())
    print(f"\nProgression to active TB (whole run): {prog.mean():.1f}% of infected individuals")
    print("\nPer-group median delays (days):")
    print(pd.DataFrame(dl).T.round(1).to_string())
    print("\nWrote analysis/derived/figure1_values.json")


if __name__ == "__main__":
    main()
