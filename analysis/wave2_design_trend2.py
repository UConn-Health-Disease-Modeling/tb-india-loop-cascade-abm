"""Wave-2 design for the second trend calibration, weighted toward a loose cascade.

Wave 1 of this round (experiment af707f22) showed that adding
`infectivity_decline_per_year` makes a loose present-day cascade *possible*:
one design point (164) reached max I = 1.53 with a 2024 treatment-initiation
fraction of 0.75 and only a 1.17x care-seeking improvement, against zero such
points in the first round.  But it is still not where the NROY concentrates --
the median passing point sits at a 2024 treatment fraction of 0.86.

A uniform wave-2 design would therefore land where round 1 already landed and
reproduce its conclusion.  This design instead allocates most of its points to
the region that matters for the manuscript's question:

    LOOSE = { 2024 treatment fraction < 0.82  AND  care-seeking improvement < 2.0x }

`LOOSE_SHARE` of the points are drawn from NROY n LOOSE and the rest from the
remainder of the NROY, so the two strata can be compared rather than one being
assumed.

Pipeline is otherwise the same as analysis/wave2_design_trend.py: a random-forest
viability gate, then log-scale GP emulators, then implausibility in log space.

Run with the notebook environment (needs sklearn/scipy):
    C:/Users/young/anaconda3/envs/tb_india/python.exe analysis/wave2_design_trend2.py

Output: analysis/wave2_design_trend2.csv, analysis/wave2_design_trend2_report.txt
"""
from __future__ import annotations
import io, sys
from pathlib import Path
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from targets_trend import TARGETS, BURDEN  # noqa: E402
from calibration_wave1_lhs_trend2 import PRIORS, FIXED, COL_ORDER  # noqa: E402

if sys.stdout.encoding and sys.stdout.encoding.lower() != "utf-8":
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
WAVE1 = "af707f22-16bc-f111-92f1-000d3af5294c"
N_WAVE2 = 100
LOOSE_SHARE = 0.70
N_CANDIDATES = 500_000
I_THRESHOLD = 3.0
SEED = 2609
VIABLE_LO, VIABLE_HI = 50.0, 800.0
TF_LOOSE, CS_LOOSE = 0.82, 2.0
NAMES = [p[0] for p in PRIORS]


def to_unit(df):
    U = np.empty((len(df), len(PRIORS)))
    for j, (n, lo, hi, mode) in enumerate(PRIORS):
        v = df[n].to_numpy(float)
        U[:, j] = ((np.log(v) - np.log(lo)) / (np.log(hi) - np.log(lo))
                   if mode == "log" else (v - lo) / (hi - lo))
    return U


def from_unit(U):
    out = {}
    for j, (n, lo, hi, mode) in enumerate(PRIORS):
        u = U[:, j]
        out[n] = (np.exp(np.log(lo) + u * (np.log(hi) - np.log(lo)))
                  if mode == "log" else lo + u * (hi - lo))
    return pd.DataFrame(out)


def main():
    from sklearn.gaussian_process import GaussianProcessRegressor
    from sklearn.gaussian_process.kernels import Matern, ConstantKernel as C
    from sklearn.model_selection import KFold, cross_val_score
    from sklearn.ensemble import RandomForestClassifier
    from sklearn.cluster import KMeans

    rng = np.random.default_rng(SEED)
    d = pd.read_csv(ROOT / "calibration" / "wave1" / "wave1_design_trend2.csv").merge(
        pd.read_csv(ROOT / "calibration" / "wave1" / "design_outputs_trend.csv"),
        on="design_id").dropna(subset=BURDEN)
    X = to_unit(d)
    viable = d["incidence_2020"].between(VIABLE_LO, VIABLE_HI).to_numpy()
    rep = [f"Wave-2 design (round 2), weighted toward a loose cascade",
           f"wave-1 experiment: {WAVE1}",
           f"design points: {len(d)}, viable: {int(viable.sum())}", ""]
    print("\n".join(rep))

    clf = RandomForestClassifier(n_estimators=500, min_samples_leaf=2,
                                 random_state=SEED).fit(X, viable)
    auc = cross_val_score(RandomForestClassifier(n_estimators=500, min_samples_leaf=2,
                                                 random_state=SEED),
                          X, viable, cv=5, scoring="roc_auc").mean()
    line = f"viability classifier 5-fold CV AUC = {auc:.3f}"
    print(line); rep.append(line)

    Xv, emus = X[viable], {}
    rep.append("\nburden emulators (log10):")
    for t in BURDEN:
        y = np.log10(np.maximum(d.loc[viable, t].to_numpy(float), 1e-3))
        kern = C(1.0) * Matern(length_scale=np.ones(Xv.shape[1]), nu=2.5)
        gp = GaussianProcessRegressor(kernel=kern, alpha=0.02, normalize_y=True,
                                      n_restarts_optimizer=5, random_state=42).fit(Xv, y)
        pred = np.empty_like(y)
        for tr, te in KFold(5, shuffle=True, random_state=1).split(Xv):
            g = GaussianProcessRegressor(kernel=kern, alpha=0.02, normalize_y=True,
                                         n_restarts_optimizer=2,
                                         random_state=42).fit(Xv[tr], y[tr])
            pred[te] = g.predict(Xv[te])
        r2 = 1 - ((y - pred) ** 2).sum() / ((y - y.mean()) ** 2).sum()
        emus[t] = gp
        line = f"  {t:18s} CV R2 = {r2:6.3f}"
        print(line); rep.append(line)

    U = rng.random((N_CANDIDATES, len(PRIORS)))
    gate = clf.predict_proba(U)[:, 1] >= 0.5
    Imax = np.full(len(U), np.inf); Imax[gate] = 0.0
    for t in BURDEN:
        z, s = TARGETS[t]
        s_log = s / (z * np.log(10))
        mu, sd = emus[t].predict(U[gate], return_std=True)
        Imax[gate] = np.maximum(Imax[gate],
                                np.abs(mu - np.log10(z)) / np.sqrt(sd ** 2 + s_log ** 2))
    keep = Imax <= I_THRESHOLD
    cand = from_unit(U[keep])
    cand_u = U[keep]
    tf2024 = 0.70 + cand["tf_gain_per_year"] * 27
    cs_mult = (1 + cand["cs_gain_per_year"]) ** 27
    loose = ((tf2024 < TF_LOOSE) & (cs_mult < CS_LOOSE)).to_numpy()
    line = (f"\nemulator NROY: {keep.mean()*100:.2f}% of the prior "
            f"({keep.sum():,} candidates); of those {loose.mean()*100:.1f}% are LOOSE "
            f"(2024 treatment fraction < {TF_LOOSE}, care-seeking < {CS_LOOSE}x)")
    print(line); rep.append(line)

    n_loose = min(int(round(N_WAVE2 * LOOSE_SHARE)), int(loose.sum()))
    n_rest = N_WAVE2 - n_loose
    line = f"allocating {n_loose} points to the loose stratum and {n_rest} to the rest"
    print(line); rep.append(line)

    picked = []
    for mask, k in ((loose, n_loose), (~loose, n_rest)):
        if k <= 0 or mask.sum() == 0:
            continue
        sub = cand_u[mask]
        km = KMeans(n_clusters=min(k, len(sub)), n_init=10, random_state=SEED).fit(sub)
        idx = np.where(mask)[0]
        chosen = {idx[np.argmin(((sub - c) ** 2).sum(axis=1))] for c in km.cluster_centers_}
        while len(chosen) < k:
            chosen.add(int(rng.choice(idx)))
        picked += sorted(chosen)[:k]

    W = from_unit(cand_u[np.array(picked)])
    W.insert(0, "design_id", np.arange(len(W)))
    for k, v in FIXED.items():
        W[k] = v
    W = W[["design_id"] + COL_ORDER]
    W.to_csv(ROOT / "calibration" / "wave2" / "wave2_design_trend2.csv", index=False)

    tf_w = 0.70 + W["tf_gain_per_year"] * 27
    cs_w = (1 + W["cs_gain_per_year"]) ** 27
    inf_w = (1 - W["infectivity_decline_per_year"]) ** 27
    rep += [f"\nwrote wave2_design_trend2.csv ({len(W)} points)",
            f"  2024 treatment fraction : {tf_w.min():.2f} - {tf_w.max():.2f} "
            f"(median {tf_w.median():.2f}); {(tf_w < TF_LOOSE).sum()} below {TF_LOOSE}",
            f"  care-seeking multiplier : x{cs_w.min():.2f} - x{cs_w.max():.2f} "
            f"(median x{cs_w.median():.2f}); {(cs_w < CS_LOOSE).sum()} below {CS_LOOSE}x",
            f"  infectivity multiplier  : x{inf_w.min():.2f} - x{inf_w.max():.2f} "
            f"(median x{inf_w.median():.2f})"]
    (ROOT / "calibration" / "wave2" / "wave2_design_trend2_report.txt").write_text("\n".join(rep), encoding="utf-8")
    print("\n".join(rep[-4:]))


if __name__ == "__main__":
    main()
