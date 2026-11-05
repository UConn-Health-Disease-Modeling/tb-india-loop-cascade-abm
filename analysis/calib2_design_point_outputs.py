"""Compute the trend-calibration targets per design point (2026-09 re-run).

Same pipeline as analysis/calib2_target_readers.py, with one change: incidence is
extracted at 2000, 2010, 2015 and 2020 instead of 2020 only, so the calibration
is gated on the shape of the epidemic as well as its level (Sourya comment #64).
2005 and 2023 are extracted too but held out of the implausibility calculation
as a check on unused years.

Usage:  python analysis/calib2_design_point_outputs.py <exp_uid>
Writes: <exp>/calibration/design_outputs.csv   (checkpointed, resumable)
        <exp>/calibration/best_fit.json
"""
from __future__ import annotations
import io, json, sqlite3, sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from calib2_target_readers import cascade_from_rer, cascade_targets, sum_pr, SCALE, BASE_YEAR  # noqa: E402
from calibration_targets import TARGETS, BURDEN, CASCADE, WHO_INCIDENCE, HOLDOUT_YEARS  # noqa: E402

if sys.stdout.encoding and sys.stdout.encoding.lower() != "utf-8":
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")

INC_YEARS = [2000, 2005, 2010, 2015, 2020, 2023]


def burden_series_from_pr(pf):
    """Incidence for every year in INC_YEARS, plus 2020 prevalence and mortality."""
    with open(pf, encoding="utf-8") as fp:
        pr = json.load(fp)
    tstep = pr["Header"]["Simulation_Timestep"]
    n = pr["Header"]["Timesteps"]
    years = (BASE_YEAR + np.arange(n) * tstep / 365.0).astype(int)
    pop = np.maximum(sum_pr(pr["Channels"], "Statistical Population", n), 1)
    active = sum_pr(pr["Channels"], "Active TB Infections", n)
    new_active = sum_pr(pr["Channels"], "New Active Infections", n)
    dd = np.diff(sum_pr(pr["Channels"], "Disease Deaths", n), prepend=0)

    out = {}
    for y in INC_YEARS:
        m = years == y
        out[f"incidence_{y}"] = (new_active[m] / pop[m]).sum() * SCALE if m.any() else np.nan
    m20 = years == 2020
    if m20.any():
        out["prevalence_2020"] = (active[m20] / pop[m20]).mean() * SCALE
        out["mortality_2020"] = (dd[m20] / pop[m20]).sum() * SCALE
    else:
        out["prevalence_2020"] = out["mortality_2020"] = np.nan
    return out


def main():
    exp = sys.argv[1]
    root = Path(__file__).resolve().parent.parent
    exp_dir = root / exp
    out_dir = exp_dir / "calibration"
    out_dir.mkdir(exist_ok=True)
    ckpt = out_dir / "design_outputs.csv"

    tags = pd.read_sql_query("SELECT * FROM results",
                             sqlite3.connect(str(exp_dir / "results.db")))
    n_dp = tags.design_id.nunique()
    print(f"results.db: {len(tags)} sims, {n_dp} design points", flush=True)

    done = {}
    if ckpt.exists():
        for _, r in pd.read_csv(ckpt).iterrows():
            done[int(r["design_id"])] = r.to_dict()
        print(f"Resuming: {len(done)}/{n_dp} design points already computed", flush=True)

    records = list(done.values())
    for did, grp in tags.groupby("design_id"):
        did = int(did)
        if did in done:
            continue
        burdens, casc_frames = [], []
        for s in grp.SIM_ID.tolist():
            pf = exp_dir / str(s) / "PropertyReportTB.json"
            if pf.exists():
                burdens.append(burden_series_from_pr(pf))
            rf = exp_dir / str(s) / "ReportEventRecorder.csv"
            if rf.exists():
                cf = cascade_from_rer(rf)
                if cf is not None:
                    casc_frames.append(cf)
        rec = {"design_id": did}
        if burdens:
            bdf = pd.DataFrame(burdens)
            rec.update(bdf.mean().to_dict())
        casc = (pd.concat(casc_frames, ignore_index=True) if casc_frames
                else pd.DataFrame(columns=["grp", "n_tests"]))
        rec.update(cascade_targets(casc))
        # annualised incidence change implied by the run
        if rec.get("incidence_2000", 0) and rec.get("incidence_2020", 0):
            rec["inc_change_pct_per_yr"] = (
                (rec["incidence_2020"] / rec["incidence_2000"]) ** (1 / 20) - 1) * 100
        records.append(rec)
        if len(records) % 5 == 0 or len(records) == n_dp:
            pd.DataFrame(records).to_csv(ckpt, index=False)
            print(f"  [{len(records)}/{n_dp}] design points done (saved)", flush=True)

    df = pd.DataFrame(records).sort_values("design_id").reset_index(drop=True)
    df.to_csv(ckpt, index=False)
    print(f"\nSaved: {ckpt}", flush=True)

    for t, (z, sig) in TARGETS.items():
        if t in df:
            df[f"I_{t}"] = (df[t] - z).abs() / sig
    Icols = [f"I_{t}" for t in TARGETS if f"I_{t}" in df]
    Iburden = [f"I_{t}" for t in BURDEN if f"I_{t}" in df]
    df["maxI_burden"] = df[Iburden].max(axis=1)
    df["maxI_all"] = df[Icols].max(axis=1)

    best = df.loc[df.maxI_burden.idxmin()]
    print("\nBest fit on the 6 calibration (burden) targets: design_id "
          f"{int(best.design_id)}  max I = {best.maxI_burden:.2f}")
    print(f"  implied incidence trend {best.get('inc_change_pct_per_yr', float('nan')):+.2f}%/yr "
          f"(WHO -2.43%/yr)")
    print(f"\n  {'target':22s} {'z':>7} {'sigma':>7} {'model':>9} {'I':>6}")
    for t, (z, sig) in TARGETS.items():
        if t in df:
            print(f"  {t:22s} {z:>7} {sig:>7} {best[t]:>9.3f} {best[f'I_{t}']:>6.2f}")
    for y in HOLDOUT_YEARS:
        k = f"incidence_{y}"
        if k in df:
            print(f"  [held out] {k}: model {best[k]:.1f} vs WHO {WHO_INCIDENCE[y]}")

    json.dump({"design_id": int(best.design_id),
               "maxI_burden": float(best.maxI_burden),
               "maxI_all": float(best.maxI_all),
               "targets": {t: {"z": TARGETS[t][0], "sigma": TARGETS[t][1],
                               "model": float(best[t]), "I": float(best[f"I_{t}"])}
                           for t in TARGETS if t in df},
               "holdout": {f"incidence_{y}": {"who": WHO_INCIDENCE[y],
                                              "model": float(best[f"incidence_{y}"])}
                           for y in HOLDOUT_YEARS if f"incidence_{y}" in df}},
              open(out_dir / "best_fit.json", "w"), indent=1)
    print(f"\nWrote {out_dir / 'best_fit.json'}")


if __name__ == "__main__":
    main()
