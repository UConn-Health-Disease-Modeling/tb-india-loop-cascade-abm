"""Wave-1 Latin-hypercube design over the 10 free parameters.

Two choices here differ from an earlier calibration of the same model:

1. `infectivity_decline_per_year` is added.  It drives a NodeInfectivityMult
   campaign event that reduces node infectivity geometrically from 1997 to 2024
   and then holds flat.  Without it the cascade is the only mechanism the model
   has to reproduce the WHO decline, which forced a 2024 treatment-initiation
   fraction of 0.80 and a 2.5x improvement in care-seeking -- leaving no headroom
   for the diagnostic and treatment scenarios to act on.  Checked against that
   earlier calibration: no design point with a 2024 treatment fraction below 0.78
   AND a care-seeking improvement below 1.6x came within 3 sigma of the burden
   targets (best max I = 3.44).  Zero is inside the prior, so the
   "cascade explains the whole decline" solution remains reachable.

2. `TB Duration` is dropped from the design.  It has no MapTo and is applied only
   inside the `if 'Death Fraction' in sample:` branch of
   map_sample_to_model_input; no calibration design has ever carried
   'Death Fraction', so the swept range never reached the model.  It stays fixed
   at the production value (3.5 y with a 0.55 case-fatality ratio, set by the
   subsample) and is reported as a fixed parameter.

Free parameters: 10.  Output: calibration/wave1/design.csv
"""
from __future__ import annotations
import argparse, io, os, sys

if sys.stdout.encoding and sys.stdout.encoding.lower() != "utf-8":
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")

import numpy as np
import pandas as pd

PRIORS = [
    ("Base_Infectivity_Constant",             0.005,  0.100,  "log"),
    ("TB Presymptomatic Duration Years",       0.15,   0.5,    "linear"),
    ("Care Seeking Fast Duration Days",        10,     30,     "log"),
    ("Care Seeking Slow Duration Days",        60,     180,    "log"),
    ("TB_Slow_Progressor_Rate",                5e-6,   5e-5,   "log"),
    ("Post_Infection_Acquisition_Multiplier",  0.20,   0.80,   "linear"),
    ("TB_Fast_Progressor_Fraction_Adult",      0.05,   0.20,   "linear"),
    ("tf_gain_per_year",                       0.000,  0.012,  "linear"),
    ("cs_gain_per_year",                       0.000,  0.050,  "linear"),
    ("infectivity_decline_per_year",           0.000,  0.020,  "linear"),
]
FIXED = {"sens_smear_pos_GXL": 0.65, "treatment_fraction_LQ": 0.55}
COL_ORDER = [n for n, *_ in PRIORS[:7]] + ["sens_smear_pos_GXL", "treatment_fraction_LQ",
                                           "tf_gain_per_year", "cs_gain_per_year",
                                           "infectivity_decline_per_year"]


def latin_hypercube(n, d, rng):
    cut = (np.arange(n)[:, None] + rng.random((n, d))) / n
    for j in range(d):
        rng.shuffle(cut[:, j])
    return cut


def transform(u, lo, hi, mode):
    if mode == "linear":
        return lo + (hi - lo) * u
    return np.exp(np.log(lo) + (np.log(hi) - np.log(lo)) * u)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, default=200)
    ap.add_argument("--seed", type=int, default=929)
    args = ap.parse_args()
    here = os.path.dirname(os.path.abspath(__file__))
    out = os.path.join(os.path.dirname(here), "calibration", "wave1", "design.csv")

    rng = np.random.default_rng(args.seed)
    u = latin_hypercube(args.n, len(PRIORS), rng)
    df = pd.DataFrame({"design_id": np.arange(args.n)})
    for j, (name, lo, hi, mode) in enumerate(PRIORS):
        df[name] = transform(u[:, j], lo, hi, mode)
    for k, v in FIXED.items():
        df[k] = v
    df = df[["design_id"] + COL_ORDER]
    df.to_csv(out, index=False)
    print(f"wrote {out}  ({len(df)} points, {len(PRIORS)} free parameters)")
    print(df[[n for n, *_ in PRIORS]].describe().T[["min", "50%", "max"]].round(5).to_string())
    mult = (1 - df["infectivity_decline_per_year"]) ** 27
    print(f"\n2024 infectivity multiplier spanned by the design: "
          f"{mult.min():.2f}x - {mult.max():.2f}x (median {mult.median():.2f}x)")


if __name__ == "__main__":
    main()
