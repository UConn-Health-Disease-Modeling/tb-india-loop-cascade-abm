"""Compute all 10 calibration targets per design point for a wave-1 experiment.

Burden targets (1-3) come from PropertyReportTB.json (consistent with the
calibration-overlay figure); cascade targets (4-10) come from
ReportEventRecorder.csv. Results are checkpointed every few design points so an
interrupted run resumes where it left off.

Outputs (in <exp>/calibration/):
    design_outputs.csv   — design_id + 10 target values (mean across reps)
    wave1_bestfit.json   — best-fit design_id, per-target implausibility, summary

Usage:
    python analysis/compute_wave1_targets.py <exp_uid>
"""
from __future__ import annotations
import glob, io, json, os, sqlite3, sys
from collections import defaultdict
from pathlib import Path

if sys.stdout.encoding and sys.stdout.encoding.lower() != "utf-8":
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")

import numpy as np
import pandas as pd

# z, sigma_total  (matches Appendix Table S-1 / emulator TARGETS)
TARGETS = {
    "incidence_2020":            (188, 32.0),
    "prevalence_2020":           (304, 36.1),
    "mortality_2020":            (33,  6.0),
    "noloop_frac":               (0.20, 0.06),
    "testloop_frac":             (0.50, 0.10),
    "testtxfail_frac":           (0.30, 0.07),
    "tests_per_case_overall":    (3.5, 0.78),
    "tests_per_case_noloop":     (1.0, 0.42),
    "tests_per_case_testloop":   (3.0, 0.92),
    "tests_per_case_testtxfail": (4.5, 1.28),
}
BURDEN = ["incidence_2020", "prevalence_2020", "mortality_2020"]
CASCADE = [k for k in TARGETS if k not in BURDEN]

SCALE = 1e5
BASE_YEAR = 1850


def sum_pr(channels, base, n):
    tot = np.zeros(n)
    for ch, val in channels.items():
        if ch.startswith(base + ":") or ch == base:
            tot += np.array(val["Data"])
    return tot


def burden_from_pr(pf):
    """2020 incidence / prevalence / mortality per 100k from one PropertyReport."""
    with open(pf, encoding="utf-8") as fp:
        pr = json.load(fp)
    tstep = pr["Header"]["Simulation_Timestep"]
    n = pr["Header"]["Timesteps"]
    years = (BASE_YEAR + np.arange(n) * tstep / 365.0).astype(int)
    pop = np.maximum(sum_pr(pr["Channels"], "Statistical Population", n), 1)
    active = sum_pr(pr["Channels"], "Active TB Infections", n)
    new_active = sum_pr(pr["Channels"], "New Active Infections", n)
    dd = np.diff(sum_pr(pr["Channels"], "Disease Deaths", n), prepend=0)
    m = years == 2020
    if not m.any():
        return None
    inc = (new_active[m] / pop[m]).sum() * SCALE
    prev = (active[m] / pop[m]).mean() * SCALE
    mort = (dd[m] / pop[m]).sum() * SCALE
    return inc, prev, mort


def cascade_from_rer(csv):
    """Per-loop-group composition + tests/case from one ReportEventRecorder."""
    df = pd.read_csv(csv, usecols=["Time", "Event_Name", "Individual_ID"])
    df["year"] = BASE_YEAR + df["Time"] / 365.0
    df = df[df["year"] >= 2010]
    act = df[df.Event_Name == "TBActivation"]["Individual_ID"].unique()
    if len(act) == 0:
        return None
    sub = df[df.Individual_ID.isin(act)]
    # vectorised per-individual event counts
    piv = sub.pivot_table(index="Individual_ID", columns="Event_Name",
                          aggfunc="size", fill_value=0)
    zero = pd.Series(0, index=piv.index)
    def col(c): return piv[c] if c in piv.columns else zero
    n_test = col("ProviderTestNoR") + col("TBTestDefault")
    n_txfail = col("TBFailedDrugRegimen")
    n_tests = col("TBTestDOTSLow") + col("TBTestDOTSHigh") + col("ProviderOrdersTBTest")
    has_test = (n_test >= 1).values
    has_tx = (n_txfail >= 1).values
    grp = np.where(has_test & has_tx, "testtxfail",
          np.where(has_test, "testloop",
          np.where(has_tx, "testtxfail", "noloop")))
    tests = np.asarray(n_tests.values, dtype=float)
    return pd.DataFrame({"grp": grp, "n_tests": tests})


def cascade_targets(group_df):
    total = len(group_df)
    if total == 0:
        return {k: np.nan for k in CASCADE}
    g = group_df.groupby("grp")
    counts = group_df.grp.value_counts()
    def frac(name): return counts.get(name, 0) / total
    def mean_tests(name):
        s = group_df[group_df.grp == name].n_tests
        return s.mean() if len(s) else np.nan
    return {
        "noloop_frac": frac("noloop"),
        "testloop_frac": frac("testloop"),
        "testtxfail_frac": frac("testtxfail"),
        "tests_per_case_overall": group_df.n_tests.mean(),
        "tests_per_case_noloop": mean_tests("noloop"),
        "tests_per_case_testloop": mean_tests("testloop"),
        "tests_per_case_testtxfail": mean_tests("testtxfail"),
    }


def main():
    exp = sys.argv[1]
    root = Path(__file__).resolve().parent.parent
    exp_dir = root / exp
    out_dir = exp_dir / "calibration"
    out_dir.mkdir(exist_ok=True)
    ckpt = out_dir / "design_outputs.csv"

    conn = sqlite3.connect(str(exp_dir / "results.db"))
    tags = pd.read_sql_query("SELECT * FROM results", conn)
    sim2did = dict(zip(tags.SIM_ID, tags.design_id))
    n_dp = tags.design_id.nunique()
    print(f"results.db: {len(tags)} sims, {n_dp} design points", flush=True)

    done = {}
    if ckpt.exists():
        prev = pd.read_csv(ckpt)
        for _, r in prev.iterrows():
            done[int(r["design_id"])] = r.to_dict()
        print(f"Resuming: {len(done)}/{n_dp} design points already computed", flush=True)

    records = list(done.values())
    for k, (did, grp) in enumerate(tags.groupby("design_id")):
        did = int(did)
        if did in done:
            continue
        sims = grp.SIM_ID.tolist()
        incs, prevs, morts, casc_frames = [], [], [], []
        for s in sims:
            pf = exp_dir / str(s) / "PropertyReportTB.json"
            if pf.exists():
                b = burden_from_pr(pf)
                if b:
                    incs.append(b[0]); prevs.append(b[1]); morts.append(b[2])
            rf = exp_dir / str(s) / "ReportEventRecorder.csv"
            if rf.exists():
                cf = cascade_from_rer(rf)
                if cf is not None:
                    casc_frames.append(cf)
        rec = {"design_id": did,
               "incidence_2020": np.mean(incs) if incs else np.nan,
               "prevalence_2020": np.mean(prevs) if prevs else np.nan,
               "mortality_2020": np.mean(morts) if morts else np.nan}
        casc = pd.concat(casc_frames, ignore_index=True) if casc_frames else pd.DataFrame(columns=["grp", "n_tests"])
        rec.update(cascade_targets(casc))
        records.append(rec)
        if len(records) % 5 == 0 or len(records) == n_dp:
            pd.DataFrame(records).to_csv(ckpt, index=False)
            print(f"  [{len(records)}/{n_dp}] design points done (saved)", flush=True)

    df = pd.DataFrame(records).sort_values("design_id").reset_index(drop=True)
    df.to_csv(ckpt, index=False)
    print(f"\nSaved design outputs: {ckpt}", flush=True)

    # ── implausibility per target, per design point ──
    for t, (z, sig) in TARGETS.items():
        df[f"I_{t}"] = (df[t] - z).abs() / sig
    Icols = [f"I_{t}" for t in TARGETS]
    df["I_max"] = df[Icols].max(axis=1)
    df["I_burden"] = df[[f"I_{t}" for t in BURDEN]].max(axis=1)
    df["I_cascade"] = df[[f"I_{t}" for t in CASCADE]].max(axis=1)

    valid = df[df.I_max <= 3.0]
    best = df.sort_values("I_max").iloc[0]
    bid = int(best.design_id)

    print("\n" + "=" * 66)
    print("WAVE-1 SUMMARY (200 design points, empirical)")
    print("=" * 66)
    print(f"  n valid (I_max <= 3):   {len(valid)} / {len(df)}")
    print(f"  min I_max   (best-fit): {df.I_max.min():.2f}")
    print(f"  min I_burden:           {df.I_burden.min():.2f}")
    print(f"  min I_cascade:          {df.I_cascade.min():.2f}")
    print(f"  best-fit design_id:     {bid}")

    print("\n  Per-target at wave-1 best-fit:")
    print(f"  {'#':>2} {'target':28s} {'z':>8} {'sigma':>7} {'model':>9} {'I_i':>6}")
    out_targets = {}
    for i, (t, (z, sig)) in enumerate(TARGETS.items(), 1):
        mv = best[t]; Ii = best[f"I_{t}"]
        out_targets[t] = {"z": z, "sigma": sig, "model": float(mv), "I": float(Ii)}
        print(f"  {i:>2} {t:28s} {z:>8} {sig:>7} {mv:>9.3f} {Ii:>6.2f}")

    result = {
        "experiment": exp,
        "n_design": int(len(df)),
        "n_valid": int(len(valid)),
        "min_I_max": float(df.I_max.min()),
        "min_I_burden": float(df.I_burden.min()),
        "min_I_cascade": float(df.I_cascade.min()),
        "best_fit_design_id": bid,
        "targets": out_targets,
    }
    with open(out_dir / "wave1_bestfit.json", "w") as fp:
        json.dump(result, fp, indent=2)
    print(f"\nSaved: {out_dir/'wave1_bestfit.json'}")


if __name__ == "__main__":
    main()
