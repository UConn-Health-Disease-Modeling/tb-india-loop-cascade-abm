"""Extract burden trajectories from the NROY propagation run.

Same channels as extract_burden_incl.py, with the parameter set carried through
as a third key, so the analysis can separate what varies because the model is
stochastic from what varies because the calibration did not pin the parameters
down.

The propagation experiment tags each simulation with Scenario_Code, Run_Number
and Param_Set; the last is what this adds.  Prevalence is taken on both
channels, as in the main analysis: Active TB Prevalence is the one the
calibration targets and the one Figure 5 reports, and the symptomatic channel
is kept for continuity with the earlier extract.

Writes analysis/derived/burden_propagation.csv.gz
"""
from __future__ import annotations
import io
import json
import sqlite3
import sys
from pathlib import Path

import numpy as np
import pandas as pd

if sys.stdout.encoding and sys.stdout.encoding.lower() != "utf-8":
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")

HERE = Path(__file__).resolve().parent
DATA = HERE / "derived"
ROOT = HERE.parent
EXP = ROOT / "f8dfe706-d7c0-f111-92f1-000d3af5294c"
SCEN = {0: "Baseline", 1: "Tx", 2: "Dx", 3: "Dx+Tx",
        4: "ACF", 5: "ACF+Tx", 6: "ACF+Dx", 7: "ACF+Dx+Tx"}
BASE_YEAR, DT = 1850, 14.0
OUT = DATA / "burden_propagation.csv.gz"


def load_tags():
    """Scenario, replicate and parameter set for every simulation."""
    con = sqlite3.connect(str(EXP / "results.db"))
    cols = {r[1] for r in con.execute("PRAGMA table_info(results)")}
    want = ["SIM_ID", "Scenario_Code", "Run_Number"]
    for extra in ("Param_Set", "Is_Production"):
        if extra in cols:
            want.append(extra)
    df = pd.read_sql_query(f"SELECT {', '.join(want)} FROM results", con)
    con.close()
    missing = [c for c in ("Param_Set",) if c not in df.columns]
    if missing:
        raise SystemExit(
            "results.db has no Param_Set column, so the parameter dimension was "
            "not tagged. Columns present: " + ", ".join(sorted(cols)))
    return df


def main():
    tags = load_tags()
    tags["Param_Set"] = tags["Param_Set"].astype(float).round().astype(int)
    tags["Scenario_Code"] = tags["Scenario_Code"].astype(float).round().astype(int)
    print(f"{len(tags)} simulations: {tags.Param_Set.nunique()} parameter sets x "
          f"{tags.Scenario_Code.nunique()} scenarios x "
          f"{tags.Run_Number.nunique()} replicates")

    rows, missing = [], 0
    for (ps, sc), grp in tags.groupby(["Param_Set", "Scenario_Code"], sort=True):
        for rep, sim in enumerate(grp.sort_values("Run_Number").SIM_ID.astype(str)):
            f = EXP / sim / "InsetChart.json"
            if not f.exists():
                missing += 1
                continue
            ch = json.load(open(f))["Channels"]
            pop = np.asarray(ch["Statistical Population"]["Data"], float)
            n = len(pop)
            yr = (BASE_YEAR + np.arange(n) * DT / 365.0).astype(int)
            df = pd.DataFrame({
                "yr": yr, "pop": pop,
                "inc": np.asarray(ch["New Active TB Infections"]["Data"], float),
                "prev": np.asarray(ch["Active Symptomatic Prevalence"]["Data"], float),
                "prevall": np.asarray(ch["Active TB Prevalence"]["Data"], float),
                "presym": np.asarray(ch["Active Presymptomatic Prevalence"]["Data"], float),
                # cumulative channel: difference the whole series before windowing
                "dth": np.diff(np.asarray(ch["Disease Deaths"]["Data"], float),
                               prepend=0.0)})
            df = df[(df.yr >= 1990) & (df.yr <= 2050)]
            g = df.groupby("yr").agg(inc=("inc", "sum"), dth=("dth", "sum"),
                                     pop=("pop", "mean"), prev=("prev", "mean"),
                                     prevall=("prevall", "mean"),
                                     presym=("presym", "mean"))
            g["inc100k"] = g["inc"] / g["pop"] * 1e5
            g["mort100k"] = g["dth"] / g["pop"] * 1e5
            g["prev100k"] = g["prev"] * 1e5
            g["prevall100k"] = g["prevall"] * 1e5
            g["presym100k"] = g["presym"] * 1e5
            g = g[["inc100k", "mort100k", "prev100k", "prevall100k",
                   "presym100k"]].reset_index()
            g["scenario"] = SCEN[sc]
            g["param_set"] = ps
            g["rep"] = rep
            rows.append(g)
        print(f"  param {ps:>5}  {SCEN[sc]:<10} {len(grp)} sims", flush=True)

    out = pd.concat(rows, ignore_index=True)
    out.to_csv(OUT, index=False, compression="gzip")
    print(f"\nmissing InsetChart files: {missing}")
    print(f"wrote {OUT}  {out.shape}")
    print(f"  parameter sets with data: {out.param_set.nunique()}")
    print(f"  scenarios with data     : {out.scenario.nunique()}")


if __name__ == "__main__":
    main()
