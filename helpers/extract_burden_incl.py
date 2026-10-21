"""Re-extract burden trajectories carrying BOTH prevalence channels.

Figure 2 calibrates to the national prevalence survey using Active TB
Prevalence, which includes the presymptomatic phase that a survey's
X-ray-plus-sputum screen also detects.  Figure 5 previously projected the
narrower Active Symptomatic Prevalence, so the figure the paper calibrated and
the figure it reported were different quantities.  This writes both, leaving
the original extract untouched for the scripts that still read it.
"""
from __future__ import annotations
import io, json, sqlite3, sys
from pathlib import Path
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from fig1_rebuild_txc import EXP  # noqa: E402

if sys.stdout.encoding and sys.stdout.encoding.lower() != "utf-8":
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")

HERE = Path(__file__).resolve().parent
SCEN = {0: "Baseline", 1: "Tx", 2: "Dx", 3: "Dx+Tx",
        4: "ACF", 5: "ACF+Tx", 6: "ACF+Dx", 7: "ACF+Dx+Tx"}
BASE_YEAR, DT = 1850, 14.0


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
                "prevall": np.asarray(ch["Active TB Prevalence"]["Data"], float),
                "presym": np.asarray(ch["Active Presymptomatic Prevalence"]["Data"], float),
                # cumulative channel: difference the whole series first
                "dth": np.diff(np.asarray(ch["Disease Deaths"]["Data"], float), prepend=0.0)})
            df = df[(df.yr >= 1990) & (df.yr <= 2050)]
            g = df.groupby("yr").agg(inc=("inc", "sum"), dth=("dth", "sum"),
                                     pop=("pop", "mean"), prev=("prev", "mean"),
                                     prevall=("prevall", "mean"), presym=("presym", "mean"))
            g["inc100k"] = g["inc"] / g["pop"] * 1e5
            g["mort100k"] = g["dth"] / g["pop"] * 1e5
            g["prev100k"] = g["prev"] * 1e5
            g["prevall100k"] = g["prevall"] * 1e5
            g["presym100k"] = g["presym"] * 1e5
            g = g[["inc100k", "mort100k", "prev100k", "prevall100k",
                   "presym100k"]].reset_index()
            g["scenario"], g["rep"] = SCEN[sc], rep
            rows.append(g)
    out = pd.concat(rows, ignore_index=True)
    out.to_csv(HERE / "burden_by_rep_incl.csv.gz", index=False, compression="gzip")
    print("wrote burden_by_rep_incl.csv.gz", out.shape)


if __name__ == "__main__":
    main()
