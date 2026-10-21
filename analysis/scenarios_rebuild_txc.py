"""[2026-09-30 round-2 production run — experiment d892bdf7-18bd-f111-92f1-000d3af5294c, loose-cascade calibration with the secular infectivity decline]

[2026-09-25 production re-run — experiment d892bdf7-18bd-f111-92f1-000d3af5294c, trend-calibrated parameters,
interventions from 2025]

Process ALL 8 production scenarios at episode level for the 2026-09 revision.

Feeds the rebuilt Figure 3 (comment #65: density curves instead of overlapping
histograms), Figure 5 (comment #67: why ACF+Tx can look worse than Tx inside a
loop group) and replicate-based uncertainty (comment #66).

Scenario_Code -> (ACF, Dx, Tx):
  0 Baseline | 1 Tx | 2 Dx | 3 Dx+Tx | 4 ACF | 5 ACF+Tx | 6 ACF+Dx | 7 ACF+Dx+Tx

Output: analysis/derived/scenario_episodes_txc.csv.gz  (one row per scenario x replicate
x loop group x activation window) and analysis/derived/scenario_dist_txc.csv.gz (per
episode infectious time, baseline scenario only, for the density panel).
"""
from __future__ import annotations
import io, sqlite3, sys
from pathlib import Path
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from fig1_rebuild_txc import process_rep, EXP, GROUPS  # noqa: E402

if sys.stdout.encoding and sys.stdout.encoding.lower() != "utf-8":
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")

HERE = Path(__file__).resolve().parent
DATA = HERE / "derived"
SCEN = {0: "Baseline", 1: "Tx", 2: "Dx", 3: "Dx+Tx",
        4: "ACF", 5: "ACF+Tx", 6: "ACF+Dx", 7: "ACF+Dx+Tx"}
WINDOWS = {"all": (0, 9999), "post2019": (2019, 9999), "post2025": (2025, 9999)}


def main():
    tags = pd.read_sql_query("SELECT SIM_ID, Scenario_Code, Run_Number FROM results",
                             sqlite3.connect(str(EXP / "results.db")))
    rows, dist = [], []
    for sc in sorted(SCEN):
        sims = tags.loc[tags.Scenario_Code == sc].sort_values("Run_Number")
        print(f"scenario {sc} ({SCEN[sc]}): {len(sims)} sims", flush=True)
        for rep, sim in enumerate(sims.SIM_ID.astype(str)):
            rf = EXP / sim / "ReportEventRecorder.csv"
            if not rf.exists():
                continue
            e, _ = process_rep(rf)
            e["infectious_days"] = e.seg_patient + e.seg_dxloop + e.seg_txfail
            e["total_days"] = (e.seg_subclin + e.seg_patient + e.seg_dxloop
                               + e.seg_txfail + e.seg_ontx)
            if sc == 0 and rep < 25:
                dist.append(e.loc[e.act_year >= 2019,
                                  ["group", "total_days", "eff_days"]].assign(rep=rep))
            for wname, (lo, hi) in WINDOWS.items():
                sub = e[(e.act_year >= lo) & (e.act_year < hi)]
                if not len(sub):
                    continue
                g = sub.groupby("group").agg(
                    n_epi=("total_days", "size"),
                    subclin_d=("seg_subclin", "sum"),
                    patient_d=("seg_patient", "sum"),
                    dxloop_d=("seg_dxloop", "sum"),
                    txfail_d=("seg_txfail", "sum"),
                    ontx_d=("seg_ontx", "sum"),
                    eff_d=("eff_days", "sum"),
                    total_d=("total_days", "sum"),
                    tests=("n_tests", "sum"),
                    testneg=("n_testneg", "sum"),
                    testdef=("n_testdef", "sum"),
                    txfail_n=("n_txfail", "sum"),
                    deaths=("died", "sum")).reset_index()
                g["scenario"] = SCEN[sc]
                g["scenario_code"] = sc
                g["rep"] = rep
                g["window"] = wname
                rows.append(g)
            if (rep + 1) % 25 == 0:
                print(f"   {SCEN[sc]} {rep+1}/{len(sims)}", flush=True)

    out = pd.concat(rows, ignore_index=True)
    out.to_csv(DATA / "scenario_episodes_txc.csv.gz", index=False,
               compression="gzip")
    print("wrote scenario_episodes_txc.csv.gz", out.shape)

    d = pd.concat(dist, ignore_index=True)
    d.to_csv(DATA / "scenario_dist_txc.csv.gz", index=False, compression="gzip")
    print("wrote scenario_dist_txc.csv.gz", d.shape)


if __name__ == "__main__":
    main()
