"""[2026-09-30 round-2 production run — experiment d892bdf7-18bd-f111-92f1-000d3af5294c, loose-cascade calibration with the secular infectivity decline]

[2026-09-25 production re-run — experiment d892bdf7-18bd-f111-92f1-000d3af5294c, trend-calibrated parameters,
interventions from 2025]

Extract per-replicate annual TB burden trajectories for all 8 scenarios.

Feeds Figure 4 with replicate-based uncertainty intervals (Sourya comment #66)
and quantifies the baseline trend question (comment #64).

Writes analysis/derived/burden_by_replicate.csv.gz with one row per
scenario x replicate x year: incidence per 100k, active symptomatic prevalence
per 100k, and TB deaths per 100k.
"""
from __future__ import annotations
import io, json, sqlite3, sys
from pathlib import Path
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from extract1_baseline_episodes import EXP  # noqa: E402

if sys.stdout.encoding and sys.stdout.encoding.lower() != "utf-8":
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")

HERE = Path(__file__).resolve().parent
DATA = HERE / "derived"
SCEN = {0: "Baseline", 1: "Tx", 2: "Dx", 3: "Dx+Tx",
        4: "ACF", 5: "ACF+Tx", 6: "ACF+Dx", 7: "ACF+Dx+Tx"}
BASE_YEAR = 1850
DT = 14.0


def main():
    tags = pd.read_sql_query("SELECT SIM_ID, Scenario_Code, Run_Number FROM results",
                             sqlite3.connect(str(EXP / "results.db")))
    rows = []
    for sc in sorted(SCEN):
        sims = tags.loc[tags.Scenario_Code == sc].sort_values("Run_Number")
        print(f"scenario {sc} ({SCEN[sc]})", flush=True)
        for rep, sim in enumerate(sims.SIM_ID.astype(str)):
            f = EXP / sim / "InsetChart.json"
            if not f.exists():
                continue
            ch = json.load(open(f))["Channels"]
            n = len(ch["New Active TB Infections"]["Data"])
            yr = (BASE_YEAR + np.arange(n) * DT / 365.0).astype(int)
            df = pd.DataFrame({
                "yr": yr,
                "inc": np.asarray(ch["New Active TB Infections"]["Data"], float),
                "pop": np.asarray(ch["Statistical Population"]["Data"], float),
                "prev": np.asarray(ch["Active Symptomatic Prevalence"]["Data"], float),
                # Disease Deaths is a CUMULATIVE channel, so it has to be
                # differenced over the WHOLE series before the year window is
                # applied; summing it raw gave a mort100k of ~78,000 per 100k.
                "dth": np.diff(np.asarray(ch["Disease Deaths"]["Data"], float),
                               prepend=0.0)})
            df = df[(df.yr >= 1990) & (df.yr <= 2050)]
            g = df.groupby("yr").agg(inc=("inc", "sum"), dth=("dth", "sum"),
                                     pop=("pop", "mean"), prev=("prev", "mean"))
            g["inc100k"] = g["inc"] / g["pop"] * 1e5
            g["mort100k"] = g["dth"] / g["pop"] * 1e5
            g["prev100k"] = g["prev"] * 1e5
            g = g[["inc100k", "mort100k", "prev100k"]].reset_index()
            g["scenario"] = SCEN[sc]
            g["rep"] = rep
            rows.append(g)
            if (rep + 1) % 25 == 0:
                print(f"   {SCEN[sc]} {rep+1}/{len(sims)}", flush=True)
    out = pd.concat(rows, ignore_index=True)
    out.to_csv(DATA / "burden_by_replicate.csv.gz", index=False, compression="gzip")
    print("wrote burden_by_replicate.csv.gz", out.shape)


if __name__ == "__main__":
    main()
