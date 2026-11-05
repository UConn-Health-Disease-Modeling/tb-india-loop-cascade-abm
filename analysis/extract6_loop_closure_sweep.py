"""Extract annual TB incidence per replicate from the loop-closure sweep.

Reads experiment 941e6423 (simulation/4_loop_closure_sweep.ipynb): each of the
8 intervention scenarios, with a further 25%, 50% or 75% of the residual loss
removed at each of the three loop drivers from 2025, at 100 replicates each.
Read by figureS2_loop_closure_sensitivity.py; the zero point of that figure is
the production run (extract3_burden_by_replicate.py).

Incidence is computed exactly as in extract3_burden_by_replicate.py: new active
TB infections summed over each calendar year, divided by the mean statistical
population of that year, per 100,000.  Replicates are numbered by Run_Number
order within each scenario x loop-reduction cell, as in extract3, so the same
seed carries the same replicate number in both tables.

Writes analysis/derived/loop_closure_sweep_incidence.csv.gz with one row per
simulation x year (1997-2035): yr, inc100k, sim, scenario_code, loop_reduction, rep.
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
EXP_PREFIX = "941e6423"
BASE_YEAR = 1850
DT = 14.0
YEARS = (1997, 2035)


def find_experiment(root: Path = ROOT) -> Path:
    hits = sorted(p for p in root.glob(EXP_PREFIX + "*") if p.is_dir())
    if len(hits) != 1:
        sys.exit(f"expected one downloaded experiment folder {EXP_PREFIX}* under {root}, "
                 f"found {len(hits)}")
    return hits[0]


def incidence(insetchart: Path) -> pd.DataFrame:
    ch = json.load(open(insetchart))["Channels"]
    n = len(ch["New Active TB Infections"]["Data"])
    yr = (BASE_YEAR + np.arange(n) * DT / 365.0).astype(int)
    df = pd.DataFrame({
        "yr": yr,
        "inc": np.asarray(ch["New Active TB Infections"]["Data"], float),
        "pop": np.asarray(ch["Statistical Population"]["Data"], float)})
    df = df[(df.yr >= YEARS[0]) & (df.yr <= YEARS[1])]
    g = df.groupby("yr").agg(inc=("inc", "sum"), pop=("pop", "mean"))
    g["inc100k"] = g["inc"] / g["pop"] * 1e5
    return g[["inc100k"]].reset_index()


def main(exp: Path | None = None, out: Path | None = None):
    exp = exp or find_experiment()
    out = out or DATA / "loop_closure_sweep_incidence.csv.gz"
    tags = pd.read_sql_query(
        "SELECT SIM_ID, Scenario_Code, Loop_Reduction, Run_Number FROM results",
        sqlite3.connect(str(exp / "results.db")))
    tags = tags.sort_values("Run_Number", kind="stable").reset_index(drop=True)
    tags["rep"] = tags.groupby(["Scenario_Code", "Loop_Reduction"]).cumcount()
    print(f"{exp.name}: {len(tags)} simulations", flush=True)

    rows = []
    for k, t in enumerate(tags.itertuples(index=False)):
        f = exp / str(t.SIM_ID) / "InsetChart.json"
        if not f.exists():
            continue
        g = incidence(f)
        g["sim"] = str(t.SIM_ID)
        g["scenario_code"] = int(t.Scenario_Code)
        g["loop_reduction"] = float(t.Loop_Reduction)
        g["rep"] = int(t.rep)
        rows.append(g)
        if (k + 1) % 200 == 0:
            print(f"   {k + 1}/{len(tags)}", flush=True)
    d = pd.concat(rows, ignore_index=True)
    d = d[["yr", "inc100k", "sim", "scenario_code", "loop_reduction", "rep"]]
    d.to_csv(out, index=False, compression="gzip")
    print(f"wrote {out.name}", d.shape)
    return d


if __name__ == "__main__":
    main()
