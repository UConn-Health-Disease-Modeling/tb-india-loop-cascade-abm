"""Combined wave-1 + wave-2 fit for the SECOND trend calibration (round 2).

Round 2 adds infectivity_decline_per_year and drops the inert TB Duration
column.  Wave 2 was weighted 70/30 toward the loose-cascade stratum, so this
script also reports the two strata separately: the question is whether a 2024
cascade loose enough to leave headroom for the diagnostic and treatment
scenarios is compatible with the WHO trajectory.

Pools the two calibration waves (200 x 10 reps + 100 x 25 reps = 300 design
points), refits the burden emulators, recomputes the NROY, and selects the
production parameter set.

Two details that matter:

* The waves carry different Monte-Carlo noise (10 vs 25 replicates), so the GP
  is given a per-point noise level alpha proportional to 1/reps rather than one
  homoscedastic WhiteKernel.  Without this the wave-1 points would drag the fit.
* Emulation is on log10 of the output, on the viable subset only, for the reason
  established in wave 2: the raw-scale spread (0 to 11,489 per 100,000) makes a
  raw GP useless.  See analysis/wave2_design_trend.py.

The production point is chosen as the SIMULATED design point with the lowest
max implausibility across the six burden targets — not an emulator optimum — so
the production run reproduces a parameter set that has actually been run.

Run with the notebook environment (needs sklearn/scipy):
    C:/Users/young/anaconda3/envs/tb_india/python.exe analysis/finalize_trend.py

Writes analysis/production_params_trend.json and analysis/finalize_trend_report.txt
"""
from __future__ import annotations
import io, json, sys
from pathlib import Path
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from targets_trend import TARGETS, BURDEN, CASCADE, WHO_INCIDENCE, HOLDOUT_YEARS  # noqa: E402
from wave2_design_trend2 import to_unit, PRIORS, VIABLE_LO, VIABLE_HI, SEED  # noqa: E402

if sys.stdout.encoding and sys.stdout.encoding.lower() != "utf-8":
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")

HERE = Path(__file__).resolve().parent
DATA = HERE / "derived"
ROOT = HERE.parent
WAVES = [
    ("wave1", "af707f22-16bc-f111-92f1-000d3af5294c", "wave1_design_trend2.csv", 10),
    ("wave2", "496ded25-37bc-f111-92f1-000d3af5294c", "wave2_design_trend2.csv", 25),
]
TF_LOOSE, CS_LOOSE = 0.82, 2.0
N_CANDIDATES = 400_000
I_THRESHOLD = 3.0
FREE = [p[0] for p in PRIORS]


def load():
    frames = []
    for name, uid, design, reps in WAVES:
        out = pd.read_csv(ROOT / "calibration" / name / "design_outputs_trend.csv")
        des = pd.read_csv(ROOT / "calibration" / name / design)
        d = des.merge(out, on="design_id")
        d["wave"] = name
        d["reps"] = reps
        d["uid"] = uid
        frames.append(d)
        print(f"{name}: {len(d)} design points x {reps} reps")
    return pd.concat(frames, ignore_index=True)


def main():
    from sklearn.gaussian_process import GaussianProcessRegressor
    from sklearn.gaussian_process.kernels import Matern, ConstantKernel as C
    from sklearn.model_selection import KFold
    from sklearn.ensemble import RandomForestClassifier

    rng = np.random.default_rng(SEED)
    d = load()
    for t, (z, s) in TARGETS.items():
        if t in d:
            d[f"I_{t}"] = (d[t] - z).abs() / s
    d["maxI_burden"] = d[[f"I_{t}" for t in BURDEN]].max(axis=1)
    d["maxI_cascade"] = d[[f"I_{t}" for t in CASCADE if f"I_{t}" in d]].max(axis=1)

    rep = [f"Combined wave-1 + wave-2 fit — trend re-calibration",
           f"design points: {len(d)} ({(d.wave=='wave1').sum()} x 10 reps + "
           f"{(d.wave=='wave2').sum()} x 25 reps)", ""]

    X = to_unit(d)
    viable = d["incidence_2020"].between(VIABLE_LO, VIABLE_HI).to_numpy()
    line = (f"viable points (2020 incidence in [{VIABLE_LO:.0f},{VIABLE_HI:.0f}]): "
            f"{viable.sum()}/{len(d)}  "
            f"(wave1 {(viable & (d.wave=='wave1')).sum()}, "
            f"wave2 {(viable & (d.wave=='wave2')).sum()})")
    print(line); rep.append(line)

    # ── viability gate (all points) ──────────────────────────────────────
    clf = RandomForestClassifier(n_estimators=500, min_samples_leaf=2,
                                 random_state=SEED).fit(X, viable)

    # ── burden emulators on log10, heteroscedastic by replicate count ────
    Xv = X[viable]
    reps_v = d.loc[viable, "reps"].to_numpy(float)
    emus, r2 = {}, {}
    rep.append("\nburden emulators (log10 scale, alpha ~ 1/reps):")
    for t in BURDEN:
        y = np.log10(np.maximum(d.loc[viable, t].to_numpy(float), 1e-3))
        alpha = 0.02 * (10.0 / reps_v)          # wave-2 points get 2.5x less noise
        kern = C(1.0) * Matern(length_scale=np.ones(Xv.shape[1]), nu=2.5)
        gp = GaussianProcessRegressor(kernel=kern, alpha=alpha, normalize_y=True,
                                      n_restarts_optimizer=5, random_state=42).fit(Xv, y)
        pred = np.empty_like(y)
        for tr, te in KFold(5, shuffle=True, random_state=1).split(Xv):
            g = GaussianProcessRegressor(kernel=kern, alpha=alpha[tr], normalize_y=True,
                                         n_restarts_optimizer=2,
                                         random_state=42).fit(Xv[tr], y[tr])
            pred[te] = g.predict(Xv[te])
        r2[t] = 1 - ((y - pred) ** 2).sum() / ((y - y.mean()) ** 2).sum()
        emus[t] = gp
        line = f"  {t:18s} CV R2 = {r2[t]:6.3f}"
        print(line); rep.append(line)

    # ── NROY over the prior ──────────────────────────────────────────────
    U = rng.random((N_CANDIDATES, len(PRIORS)))
    gate = clf.predict_proba(U)[:, 1] >= 0.5
    Imax = np.full(len(U), np.inf)
    Imax[gate] = 0.0
    for t in BURDEN:
        z, s = TARGETS[t]
        s_log = s / (z * np.log(10))
        mu, sd = emus[t].predict(U[gate], return_std=True)
        Imax[gate] = np.maximum(Imax[gate],
                                np.abs(mu - np.log10(z)) / np.sqrt(sd ** 2 + s_log ** 2))
    nroy = (Imax <= I_THRESHOLD).mean()
    line = (f"\nNROY volume: {nroy*100:.2f}% of the prior "
            f"(wave-1 only was 6.6%; viability gate keeps {gate.mean()*100:.1f}%)")
    print(line); rep.append(line)

    # ── loose vs tight cascade: the question this round was run to answer ──
    d["tf2024"] = 0.70 + d["tf_gain_per_year"] * 27
    d["cs_mult"] = (1 + d["cs_gain_per_year"]) ** 27
    d["inf_mult"] = (1 - d["infectivity_decline_per_year"]) ** 27
    d["loose"] = (d.tf2024 < TF_LOOSE) & (d.cs_mult < CS_LOOSE)
    rep.append("")
    rep.append("loose vs tight cascade (loose = 2024 treatment fraction < "
               f"{TF_LOOSE} and care-seeking < {CS_LOOSE}x):")
    for lab, sub in (("loose", d[d.loose]), ("tight", d[~d.loose])):
        if not len(sub):
            continue
        n3 = int((sub.maxI_burden <= 3).sum())
        line = (f"  {lab:5s} n={len(sub):3d}  passing 3 sigma: {n3:3d} "
                f"({n3/len(sub)*100:4.1f}%)  best max I = {sub.maxI_burden.min():.2f}")
        print(line); rep.append(line)
    lb = d[d.loose]
    if len(lb):
        b2 = lb.loc[lb.maxI_burden.idxmin()]
        line = (f"  best LOOSE point: {b2.wave} design {int(b2.design_id)}  "
                f"max I {b2.maxI_burden:.2f}  trend {b2.inc_change_pct_per_yr:+.2f}%/yr  "
                f"tf2024 {b2.tf2024:.3f}  care-seeking x{b2.cs_mult:.2f}  "
                f"infectivity x{b2.inf_mult:.2f}")
        print(line); rep.append(line)

    # ── production point: best SIMULATED design point ────────────────────
    best = d.loc[d.maxI_burden.idxmin()]
    line = (f"\nPRODUCTION POINT: {best.wave} design_id {int(best.design_id)}  "
            f"(experiment {best.uid})\n"
            f"  max I on the 6 burden targets = {best.maxI_burden:.2f}\n"
            f"  max I on the 7 cascade validation targets = {best.maxI_cascade:.2f}\n"
            f"  implied 2000-2020 trend {best.inc_change_pct_per_yr:+.2f}%/yr "
            f"(WHO {((WHO_INCIDENCE[2020]/WHO_INCIDENCE[2000])**(1/20)-1)*100:+.2f}%/yr)")
    print(line); rep.append(line)

    rep.append("\n  calibration targets:")
    for t in BURDEN:
        z, s = TARGETS[t]
        rep.append(f"    {t:18s} model {best[t]:8.1f}   target {z:>5} +/- {s:<5}  I={best[f'I_{t}']:.2f}")
    rep.append("\n  cascade validation targets (3-sigma check only):")
    for t in CASCADE:
        z, s = TARGETS[t]
        rep.append(f"    {t:26s} model {best[t]:7.3f}   target {z:>5} +/- {s:<5}  I={best[f'I_{t}']:.2f}")
    rep.append("\n  held-out years:")
    for y in HOLDOUT_YEARS:
        rep.append(f"    incidence_{y}: model {best[f'incidence_{y}']:.1f} vs WHO {WHO_INCIDENCE[y]}")
    rep.append("\n  production parameters:")
    for p in FREE:
        rep.append(f"    {p:38s} {best[p]:.6g}")
    rep.append(f"    {'sens_smear_pos_GXL':38s} {best['sens_smear_pos_GXL']:.6g}  (fixed)")
    rep.append(f"    {'treatment_fraction_LQ':38s} {best['treatment_fraction_LQ']:.6g}  (fixed)")

    top = d.nsmallest(8, "maxI_burden")[
        ["wave", "design_id", "maxI_burden", "maxI_cascade", "inc_change_pct_per_yr",
         "incidence_2020", "tf_gain_per_year", "cs_gain_per_year"]]
    rep.append("\n  top 8 design points by burden fit:")
    rep.append(top.round(4).to_string(index=False))

    json.dump({"production": {p: float(best[p]) for p in FREE}
               | {"sens_smear_pos_GXL": float(best["sens_smear_pos_GXL"]),
                  "treatment_fraction_LQ": float(best["treatment_fraction_LQ"])},
               "source": {"wave": str(best.wave), "design_id": int(best.design_id),
                          "experiment": str(best.uid)},
               "maxI_burden": float(best.maxI_burden),
               "maxI_cascade": float(best.maxI_cascade),
               "trend_pct_per_yr": float(best.inc_change_pct_per_yr),
               "nroy_fraction": float(nroy),
               "emulator_cv_r2": {k: float(v) for k, v in r2.items()},
               "targets": {t: {"z": TARGETS[t][0], "sigma": TARGETS[t][1],
                               "model": float(best[t]), "I": float(best[f"I_{t}"])}
                           for t in TARGETS if t in d}},
              open(DATA / "production_params_trend2.json", "w"), indent=1)
    (DATA / "finalize_trend2_report.txt").write_text("\n".join(rep), encoding="utf-8")
    print("\n".join(rep[rep.index("\n  production parameters:"):]))
    print(f"\nwrote analysis/production_params_trend.json")


if __name__ == "__main__":
    main()
