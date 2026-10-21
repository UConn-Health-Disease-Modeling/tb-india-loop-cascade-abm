"""[2026-09-30 round-2 production run — experiment d892bdf7-18bd-f111-92f1-000d3af5294c, loose-cascade calibration with the secular infectivity decline]

Figure 2 for the manuscript — production run, 2026-09-25.

This is the figure that answers the reviewer comment about the baseline trend.
It plots the calibrated model against the full WHO incidence series instead of a
single 2020 anchor, using the 25 replicates of the combined wave-1 + wave-2 best
fit (wave-2 design 8, experiment fcc9da93).

Unlike figures/Figure_0909/Fig2_calibration.png — which showed the PUBLISHED
model before re-calibration — the model line here is fitted to the four WHO
incidence years, so the diamonds are targets it was actually fitted to.

Note this is the calibration point run in Baseline configuration only; the
manuscript's other figures need the production run (8 scenarios x 100 reps).

Writes figures/Figure_TxCure/Fig2_calibration_trend.png
"""
from __future__ import annotations
import io, json, sqlite3, sys
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).resolve().parent))
from targets_trend import WHO_INCIDENCE, TREND_YEARS, TARGETS  # noqa: E402

if sys.stdout.encoding and sys.stdout.encoding.lower() != "utf-8":
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")

mpl.rcParams["font.family"] = "DejaVu Sans"
mpl.rcParams["axes.unicode_minus"] = False
# 2026-10-01: these figures are drawn wide and printed at column width, so the
# defaults that matplotlib applies to anything without an explicit fontsize
# (tick labels above all) have to be raised along with the explicit sizes.
mpl.rcParams.update({"font.size": 14, "axes.titlesize": 17, "axes.labelsize": 15,
                     "xtick.labelsize": 13.5, "ytick.labelsize": 13.5,
                     "legend.fontsize": 13.5, "figure.titlesize": 18})

HERE = Path(__file__).resolve().parent
DATA = HERE / "derived"
ROOT = HERE.parent
OUT = ROOT / "figures" / "Figure_TxCure"
OUT.mkdir(parents=True, exist_ok=True)
EXP = ROOT / "d892bdf7-18bd-f111-92f1-000d3af5294c"
CACHE = DATA / "production_burden_txc.csv.gz"
BASE_YEAR, DT = 1850, 14.0


def extract():
    sims = pd.read_sql_query(
        "SELECT SIM_ID FROM results WHERE Scenario_Code=0",
        sqlite3.connect(str(EXP / "results.db"))).SIM_ID.astype(str)
    rows = []
    for rep, s in enumerate(sims):
        f = EXP / s / "InsetChart.json"
        if not f.exists():
            continue
        ch = json.load(open(f))["Channels"]
        n = len(ch["New Active TB Infections"]["Data"])
        yr = (BASE_YEAR + np.arange(n) * DT / 365.0).astype(int)
        df = pd.DataFrame({
            "yr": yr,
            "inc": np.asarray(ch["New Active TB Infections"]["Data"], float),
            "pop": np.asarray(ch["Statistical Population"]["Data"], float),
            "prev": np.asarray(ch["Active TB Prevalence"]["Data"], float),
            # Disease Deaths is a CUMULATIVE channel
            "dth": np.diff(np.asarray(ch["Disease Deaths"]["Data"], float), prepend=0.0)})
        df = df[(df.yr >= 1995) & (df.yr <= 2035)]
        g = df.groupby("yr").agg(inc=("inc", "sum"), dth=("dth", "sum"),
                                 pop=("pop", "mean"), prev=("prev", "mean"))
        rows.append(pd.DataFrame({"yr": g.index,
                                  "inc100k": g["inc"] / g["pop"] * 1e5,
                                  "mort100k": g["dth"] / g["pop"] * 1e5,
                                  "prev100k": g["prev"] * 1e5, "rep": rep}))
    d = pd.concat(rows, ignore_index=True)
    d.to_csv(CACHE, index=False, compression="gzip")
    print(f"extracted {d.rep.nunique()} replicates")
    return d


def main():
    d = pd.read_csv(CACHE) if CACHE.exists() else extract()
    n_rep = d.rep.nunique()

    fig, axes = plt.subplots(1, 3, figsize=(17, 5.4))
    panels = [("inc100k", "A1. TB incidence", "Incidence (per 100,000/yr)", "#2171B5"),
              ("prev100k", "A2. Active TB prevalence", "Prevalence (per 100,000)", "#238B45"),
              ("mort100k", "A3. TB mortality", "Mortality (per 100,000/yr)", "#CB181D")]
    for ax, (col, title, ylab, colr) in zip(axes, panels):
        g = d.groupby("yr")[col]
        m, half = g.mean(), g.sem() * 2.064          # t(0.975, df=24)
        lo, hi = g.quantile(0.25), g.quantile(0.75)
        ax.fill_between(m.index, lo.values, hi.values, color=colr, alpha=0.12, lw=0,
                        label="Model IQR across replicates")
        ax.fill_between(m.index, (m - half).values, (m + half).values, color=colr,
                        alpha=0.30, lw=0, label="Model mean 95% CI")
        ax.plot(m.index, m.values, color=colr, lw=2.4, label="Model (calibrated)")
        ax.set_title(title, fontsize=17.6, fontweight="bold", loc="left")
        ax.set_xlabel("Year", fontsize=16.2)
        ax.set_ylabel(ylab, fontsize=16.2)
        ax.grid(alpha=0.25)
        for sp in ("top", "right"):
            ax.spines[sp].set_visible(False)

    ax = axes[0]
    ys = sorted(WHO_INCIDENCE)
    ax.plot(ys, [WHO_INCIDENCE[y] for y in ys], "o-", color="#B30000", ms=6, lw=1.8,
            label="WHO estimate (India)")
    tv = [WHO_INCIDENCE[y] for y in TREND_YEARS]
    ax.errorbar(TREND_YEARS, tv, yerr=[0.17 * v for v in tv], fmt="D", color="#B30000",
                ms=8, capsize=5, lw=0, elinewidth=1.4,
                label="Calibration targets (±1σ)")
    mA = d.groupby("yr").inc100k.mean()
    mod = ((mA.loc[2020] / mA.loc[2000]) ** (1 / 20) - 1) * 100
    who = ((WHO_INCIDENCE[2020] / WHO_INCIDENCE[2000]) ** (1 / 20) - 1) * 100
    txt = ("2000-2020 trend" + chr(10) + f"  model {mod:+.2f}%/yr"
           + chr(10) + f"  WHO   {who:+.2f}%/yr")
    ax.text(0.03, 0.06, txt, transform=ax.transAxes, fontsize=14.9, va="bottom",
            bbox=dict(boxstyle="round,pad=0.4", fc="#F2F7FB", ec="#2171B5", lw=1.2))
    ax.set_ylim(0, max(max(WHO_INCIDENCE.values()), mA.max()) * 1.15)
    # below the axes: inside the panel the entries sat on top of the trajectory
    ax.legend(fontsize=12.8, frameon=False, loc="upper center",
              bbox_to_anchor=(0.5, -0.22), ncol=1)

    axes[1].errorbar([2020], [304], yerr=[36], fmt="D", color="#00441B", ms=8, capsize=5,
                     label="National Prevalence Survey 2019-21 (304)")
    axes[1].errorbar([2011], [383], yerr=[60], fmt="^", color="#41AB5D", ms=8, capsize=5,
                     label="Gujarat survey (383)")
    axes[1].legend(fontsize=12.8, frameon=False, loc="upper center",
                   bbox_to_anchor=(0.5, -0.22), ncol=1)
    axes[2].errorbar([2020], [33], yerr=[6], fmt="D", color="#67000D", ms=8, capsize=5,
                     label="WHO 2020 (33)")
    axes[2].legend(fontsize=12.8, frameon=False, loc="upper center",
                   bbox_to_anchor=(0.5, -0.22), ncol=1)

    # all three panels start at zero: a truncated axis exaggerates the decline
    # (reviewer comment).  The top is left where autoscaling put it, so the
    # survey markers and their error bars stay inside the axes.
    for ax in axes[1:]:
        ax.set_ylim(bottom=0)

    plt.tight_layout()
    p = OUT / "Fig2_calibration_trend.png"
    plt.savefig(p, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close()
    print("saved", p)

    mm = d.groupby("yr")[["inc100k", "prev100k", "mort100k"]].mean()
    for y in (2000, 2010, 2015, 2020, 2035):
        who_v = WHO_INCIDENCE.get(y, "-")
        print(f"  {y}: incidence {mm.loc[y,'inc100k']:6.1f} (WHO {who_v})"
              f"  prevalence {mm.loc[y,'prev100k']:6.1f}  mortality {mm.loc[y,'mort100k']:5.1f}")
    print(f"model trend {mod:+.2f}%/yr vs WHO {who:+.2f}%/yr")


if __name__ == "__main__":
    main()
