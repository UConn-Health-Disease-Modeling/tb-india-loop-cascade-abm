"""Download from IDM COMPS the raw simulation outputs that the analysis scripts read.

For each experiment, every simulation's listed output files are downloaded into
<dest>/<experiment ID>/<simulation ID>/, and the simulation tags are written to
<dest>/<experiment ID>/results.db.  The extraction scripts in analysis/ expect
these folders at the top of the repository, which is the default <dest>.

Needs the simulation environment (requirements-simulation.txt) and a COMPS
account; you are asked to log in on the first run.

Usage, from the top folder of the repository:
    python simulation/download_outputs.py                  # all five experiments
    python simulation/download_outputs.py loopsweep        # one or more by name
    python simulation/download_outputs.py --dest D:/tb_raw # somewhere else

Experiments, in the order they are downloaded:
    production   d892bdf7  InsetChart.json, ReportEventRecorder.csv
    loopsweep    941e6423  InsetChart.json
    propagation  f8dfe706  InsetChart.json
    wave1        af707f22  PropertyReportTB.json, ReportEventRecorder.csv
    wave2        496ded25  PropertyReportTB.json, ReportEventRecorder.csv
"""
from __future__ import annotations
import argparse
import os
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]

EXPERIMENTS = {
    "production":  ("d892bdf7-18bd-f111-92f1-000d3af5294c", ["InsetChart.json", "ReportEventRecorder.csv"]),
    "loopsweep":   (None,                                   ["InsetChart.json"]),
    "propagation": ("f8dfe706-d7c0-f111-92f1-000d3af5294c", ["InsetChart.json"]),
    "wave1":       ("af707f22-16bc-f111-92f1-000d3af5294c", ["PropertyReportTB.json", "ReportEventRecorder.csv"]),
    "wave2":       ("496ded25-37bc-f111-92f1-000d3af5294c", ["PropertyReportTB.json", "ReportEventRecorder.csv"]),
}


def loopsweep_id(platform) -> str:
    """Full ID of the loop-closure sweep, looked up from one of its simulations.

    Only the short ID 941e6423 is recorded in the README; the derived table
    carries the full simulation IDs, and COMPS knows each one's experiment.
    """
    from idmtools.core import ItemType
    sim_id = pd.read_csv(ROOT / "analysis" / "derived" / "loop_closure_sweep_incidence.csv.gz",
                         usecols=["sim"], nrows=1).sim.iloc[0]
    sim = platform.get_item(sim_id, item_type=ItemType.SIMULATION)
    exp_id = str(sim.parent_id)
    if not exp_id.startswith("941e6423"):
        raise SystemExit(f"simulation {sim_id} belongs to {exp_id}, not 941e6423")
    return exp_id


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("names", nargs="*",
                    help=f"experiments to download, from {', '.join(EXPERIMENTS)} (default: all)")
    ap.add_argument("--dest", default=str(ROOT), help="folder to download into")
    args = ap.parse_args()
    names = args.names or list(EXPERIMENTS)
    unknown = [n for n in names if n not in EXPERIMENTS]
    if unknown:
        ap.error(f"unknown experiment(s): {', '.join(unknown)}")

    from idmtools.core.platform_factory import Platform
    from emodpy.emod_task import EMODTask
    platform = Platform("CALCULON")

    dest = Path(args.dest).resolve()
    dest.mkdir(parents=True, exist_ok=True)
    os.chdir(dest)  # the COMPS helpers write into <cwd>/<experiment ID>/
    for name in names:
        exp_id, files = EXPERIMENTS[name]
        exp_id = exp_id or loopsweep_id(platform)
        print(f"\n== {name}: {exp_id}  ({', '.join(files)})", flush=True)
        EMODTask.get_file_from_comps(exp_id, files)
        EMODTask.cache_experiment_metadata_in_sql(exp_id)
        n = {f: len(list((dest / exp_id).glob(f"**/{f}"))) for f in files}
        print(f"   downloaded: {n}; results.db: {(dest / exp_id / 'results.db').exists()}")


if __name__ == "__main__":
    main()
