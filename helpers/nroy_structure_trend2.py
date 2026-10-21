"""Structure of the NROY region from the round-2 trend calibration.

Table S-7 reports that eight of ten parameters keep ~90% of their prior range
inside the NROY region.  This script asks what shape that region has: which
parameters trade off against one another (rank correlations among retained
draws), how much a parameter's retained range shrinks once another is fixed
(conditional spans), and how little the care-cascade benchmarks vary across
burden-plausible design points (why they could not break the trade-offs).

Same emulator pipeline as nroy_marginals_trend2.py (viability gate + six
log10 burden emulators + 400,000 prior draws, 3-sigma threshold).  Read-only.

Run with the notebook environment:
    C:/Users/young/anaconda3/envs/tb_india/python.exe helpers/nroy_structure_trend2.py
Writes helpers/nroy_structure_trend2.txt and helpers/nroy_corr_trend2.csv
"""
from __future__ import annotations
import io, sys
from pathlib import Path
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from targets_trend import TARGETS, BURDEN, CASCADE  # noqa: E402
from wave2_design_trend2 import to_unit, PRIORS, VIABLE_LO, VIABLE_HI, SEED  # noqa: E402

if sys.stdout.encoding and sys.stdout.encoding.lower() != "utf-8":
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
WAVES = [
    ("wave1", "af707f22-16bc-f111-92f1-000d3af5294c", "wave1_design_trend2.csv", 10),
    ("wave2", "496ded25-37bc-f111-92f1-000d3af5294c", "wave2_design_trend2.csv", 25),
]
N_CANDIDATES = 400_000
I_THRESHOLD = 3.0
SHORT = {
    "Base_Infectivity_Constant": "infectivity",
    "TB Presymptomatic Duration Years": "presymp_dur",
    "Care Seeking Fast Duration Days": "cs_delay_fast",
    "Care Seeking Slow Duration Days": "cs_delay_slow",
    "TB_Slow_Progressor_Rate": "slow_prog_rate",
    "Post_Infection_Acquisition_Multiplier": "reinfect_mult",
    "TB_Fast_Progressor_Fraction_Adult": "fast_frac",
    "tf_gain_per_year": "tx_init_gain",
    "cs_gain_per_year": "cs_improve",
    "infectivity_decline_per_year": "inf_decline",
}


def load():
    frames = []
    for name, uid, design, reps in WAVES:
        out = pd.read_csv(ROOT / uid / "calibration" / "design_outputs_trend.csv")
        des = pd.read_csv(HERE / design)
        d = des.merge(out, on="design_id")
        d["wave"], d["reps"] = name, reps
        frames.append(d)
    return pd.concat(frames, ignore_index=True)


def main():
    from sklearn.gaussian_process import GaussianProcessRegressor
    from sklearn.gaussian_process.kernels import Matern, ConstantKernel as C
    from sklearn.ensemble import RandomForestClassifier

    rng = np.random.default_rng(SEED)
    d = load()
    rep = []

    def say(s=""):
        print(s)
        rep.append(s)

    for t, (z, s) in TARGETS.items():
        if t in d:
            d[f"I_{t}"] = (d[t] - z).abs() / s
    d["maxI_burden"] = d[[f"I_{t}" for t in BURDEN]].max(axis=1)

    X = to_unit(d)
    viable = d["incidence_2020"].between(VIABLE_LO, VIABLE_HI).to_numpy()
    clf = RandomForestClassifier(n_estimators=500, min_samples_leaf=2,
                                 random_state=SEED).fit(X, viable)
    Xv = X[viable]
    reps_v = d.loc[viable, "reps"].to_numpy(float)
    emus = {}
    for t in BURDEN:
        y = np.log10(np.maximum(d.loc[viable, t].to_numpy(float), 1e-3))
        alpha = 0.02 * (10.0 / reps_v)
        kern = C(1.0) * Matern(length_scale=np.ones(Xv.shape[1]), nu=2.5)
        emus[t] = GaussianProcessRegressor(kernel=kern, alpha=alpha, normalize_y=True,
                                           n_restarts_optimizer=5, random_state=42).fit(Xv, y)

    U = rng.random((N_CANDIDATES, len(PRIORS)))
    gate = clf.predict_proba(U)[:, 1] >= 0.5
    Imax = np.full(len(U), np.inf)
    Imax[gate] = 0.0
    for t in BURDEN:
        z, s = TARGETS[t]
        s_log = s / (z * np.log(10))
        mu, sd = emus[t].predict(U[gate], return_std=True)
        Imax[gate] = np.maximum(Imax[gate], np.abs(mu - np.log10(z)) / np.sqrt(sd ** 2 + s_log ** 2))
    keep = Imax <= I_THRESHOLD
    say(f"NROY: {keep.mean()*100:.2f}% of prior ({keep.sum():,} of {N_CANDIDATES:,} draws); "
        f"viability gate keeps {gate.mean()*100:.1f}%")

    names = [SHORT.get(p[0], p[0]) for p in PRIORS]
    Un = pd.DataFrame(U[keep], columns=names)          # unit sampling scale
    Up = pd.DataFrame(U, columns=names)

    # ---- 1. rank correlations among retained draws --------------------------
    corr = Un.corr(method="spearman")
    corr.to_csv(HERE / "nroy_corr_trend2.csv")
    pairs = sorted(((a, b, corr.loc[a, b]) for i, a in enumerate(names) for b in names[i + 1:]),
                   key=lambda x: -abs(x[2]))
    say("\n1. Spearman correlations among NROY draws (prior draws are independent, so any "
        "correlation is induced by the targets):")
    for a, b, r in pairs[:10]:
        say(f"   {r:+.2f}  {a:15s} ~ {b}")
    say(f"   (remaining {len(pairs)-10} pairs all |r| <= {abs(pairs[10][2]):.2f})")

    # ---- 2. conditional spans: fix one parameter to a tercile, re-measure others ----
    def span(u):
        q5, q95 = np.percentile(u, [5, 95])
        return 100 * (q95 - q5)

    say("\n2. Retained 5-95% span of each parameter (percent of prior), overall and within "
        "terciles of base infectivity:")
    terc = pd.qcut(Un["infectivity"], 3, labels=["low", "mid", "high"])
    say("   %-16s%9s%9s%9s%9s" % ("parameter", "overall", "inf low", "inf mid", "inf high"))
    for p in names:
        row = [span(Un[p])] + [span(Un.loc[terc == k, p]) for k in ("low", "mid", "high")]
        say("   %-16s" % p + "".join("%8.0f%%" % v for v in row))
    say("   median of each parameter within infectivity terciles (unit scale, 0-1):")
    for p in ("slow_prog_rate", "fast_frac", "reinfect_mult", "cs_delay_slow"):
        meds = [Un.loc[terc == k, p].median() for k in ("low", "mid", "high")]
        say("   %-16s" % p + "".join("%8.2f " % v for v in meds))

    # ---- 3. how much the targets moved each marginal at all ---------------------
    say("\n3. Marginal shift: median of retained draws minus 0.5 (prior median), unit scale:")
    for p in names:
        say(f"   {p:16s} {Un[p].median()-0.5:+.3f}   span {span(Un[p]):.0f}%")

    # ---- 4. cascade benchmarks across burden-plausible SIMULATED points ----------
    ok = d.maxI_burden <= I_THRESHOLD
    say(f"\n4. Care-cascade benchmarks across the {ok.sum()} simulated design points with "
        f"burden max I <= 3 (of {len(d)}):")
    say("   %-26s%9s%9s%8s%9s%18s" % ("benchmark", "SD(all)", "SD(ok)", "sigma", "SD/sigma", "range (ok)"))
    for c in CASCADE:
        z, s = TARGETS[c]
        a = d.loc[ok, c]
        say("   %-26s%9.3f%9.3f%8.2f%9.2f   %6.2f - %6.2f" % (c, d[c].std(), a.std(), s, a.std() / s,
                                                                a.min(), a.max()))
    fixed = [c for c in ("sens_smear_pos_GXL", "treatment_fraction_LQ") if c in d]
    for c in fixed:
        say(f"   {c}: {'constant' if d[c].nunique() == 1 else 'VARIES'} across design "
            f"(value {d[c].iloc[0]:.4g})")

    # which free parameters do the cascade outputs respond to, among ok points?
    say("\n   Spearman correlation of each cascade benchmark with the free parameters "
        "(burden-plausible simulated points):")
    free_cols = [p[0] for p in PRIORS]
    Xo = d.loc[ok, free_cols].rename(columns=SHORT)
    for c in CASCADE:
        r = Xo.corrwith(d.loc[ok, c], method="spearman").sort_values(key=np.abs, ascending=False)
        say(f"   {c:26s} " + "  ".join(f"{k} {v:+.2f}" for k, v in r.head(3).items()))

    (HERE / "nroy_structure_trend2.txt").write_text("\n".join(rep), encoding="utf-8")
    print(f"\nwrote {HERE / 'nroy_structure_trend2.txt'}")


if __name__ == "__main__":
    main()
