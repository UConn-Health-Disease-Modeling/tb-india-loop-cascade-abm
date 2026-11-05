"""Sensitivity of intervention impact to further cascade-loop closure.

Previously Panel C of Figure 4; separated out because Figure 4 now carries the
three End-TB burden measures.  Experiment 941e6423 removes a further fraction
of the residual loss at each of the three loop drivers from 2025, on top of
each scenario (24 design points x 100 replicates); the zero point is the
production run.

The y axis is broken: every curve lies below 35% while the End-TB target sits
at 85% on the vs-Baseline scale, so an unbroken axis leaves five sixths of the
panel empty and crushes the eight curves into the bottom quarter.

Writes figures/reproduced/FigureS2_loop_closure_sensitivity.png
"""
from __future__ import annotations
import io, json, sys
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt

if sys.stdout.encoding and sys.stdout.encoding.lower() != "utf-8":
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")

mpl.rcParams["font.family"] = "DejaVu Sans"
mpl.rcParams["axes.unicode_minus"] = False
mpl.rcParams.update({"font.size": 14, "axes.titlesize": 17, "axes.labelsize": 15,
                     "xtick.labelsize": 13.5, "ytick.labelsize": 13.5,
                     "legend.fontsize": 13.5})

HERE = Path(__file__).resolve().parent
DATA = HERE / "derived"
OUT = HERE.parent / "figures" / "reproduced"
OUT.mkdir(parents=True, exist_ok=True)

ORDER = ["Baseline", "ACF", "Dx", "Tx", "ACF+Dx", "ACF+Tx", "Dx+Tx", "ACF+Dx+Tx"]
LAB = {"Baseline": "Loop closure alone", "ACF": "Active case finding (ACF)",
       "Dx": "Diagnostic improvement (Dx)", "Tx": "Treatment support (Tx)",
       "ACF+Dx": "ACF + Dx", "ACF+Tx": "ACF + Tx",
       "Dx+Tx": "Dx + Tx", "ACF+Dx+Tx": "ACF + Dx + Tx"}
COL = {"Baseline": "#333333", "ACF": "#4292C6", "Dx": "#41AB5D", "Tx": "#F16913",
       "ACF+Dx": "#6BAED6", "ACF+Tx": "#FD8D3C", "Dx+Tx": "#74C476",
       "ACF+Dx+Tx": "#A50F15"}
SCEN = {0: "Baseline", 1: "Tx", 2: "Dx", 3: "Dx+Tx",
        4: "ACF", 5: "ACF+Tx", 6: "ACF+Dx", 7: "ACF+Dx+Tx"}


def main():
    d = pd.read_csv(DATA / "burden_by_replicate.csv.gz")
    d = d[(d.yr >= 2000) & (d.yr <= 2035)]
    sweep = pd.read_csv(DATA / "loop_closure_sweep_incidence.csv.gz")
    for name, df in (("production", d), ("loop sweep", sweep)):
        alive = df.groupby("rep").inc100k.max()
        dead = sorted(alive.index[alive == 0])
        if dead:
            print(f"{name}: dropping degenerate replicate(s) {dead} of {alive.size}")
    d = d[~d.rep.isin(d.groupby("rep").inc100k.max().pipe(lambda s: s.index[s == 0]))]
    sweep = sweep[~sweep.rep.isin(
        sweep.groupby("rep").inc100k.max().pipe(lambda s: s.index[s == 0]))]
    sweep["scenario"] = sweep.scenario_code.map(SCEN)

    # End-TB targets are reductions against 2015; this panel measures reduction
    # against the 2035 baseline, so the target lines have to be converted.
    y2015 = d[(d.scenario == "Baseline") & (d.yr == 2015)].inc100k.mean()
    base35 = d[(d.scenario == "Baseline") & (d.yr == 2035)].inc100k.mean()
    RATIO = base35 / y2015

    def vs_base(target):
        return (1 - (1 - target) / RATIO) * 100

    def vs_2015(r):
        return (1 - RATIO * (1 - r / 100)) * 100

    pts = []
    for sc in ORDER:                                   # f = 0: the production run
        v = d[(d.scenario == sc) & (d.yr == 2035)].inc100k
        r = (1 - v / base35) * 100
        pts.append(dict(scenario=sc, f=0.0, mean=r.mean(), ci=r.sem() * 1.984))
    s35 = sweep[sweep.yr == 2035]
    for (sc, f), g in s35.groupby(["scenario", "loop_reduction"]):
        r = (1 - g.inc100k / base35) * 100
        pts.append(dict(scenario=sc, f=f, mean=r.mean(), ci=r.sem() * 1.984))
    pts = pd.DataFrame(pts)

    fig = plt.figure(figsize=(12.5, 8.2))
    gs = fig.add_gridspec(2, 1, height_ratios=[1, 9], hspace=0.10)
    axt = fig.add_subplot(gs[0])                 # the End-TB line only
    ax = fig.add_subplot(gs[1], sharex=axt)

    for a in (ax, axt):
        for sc in ORDER:
            g = pts[pts.scenario == sc].sort_values("f")
            x = g.f * 100
            lw = 2.8 if sc in ("Baseline", "ACF+Dx+Tx") else 1.8
            a.plot(x, g["mean"], "-o", color=COL[sc], lw=lw, ms=5.5, label=LAB[sc],
                   zorder=3 if sc in ("Baseline", "ACF+Dx+Tx") else 2)
            a.fill_between(x, g["mean"] - g.ci, g["mean"] + g.ci,
                           color=COL[sc], alpha=0.12, lw=0)
        a.grid(alpha=0.25)
        a.set_xlim(-2.5, 77.5)

    hi = float((pts["mean"] + pts.ci).max())
    C_TOP = 5 * np.ceil(hi * 1.03 / 5)
    end = vs_base(0.90)
    ax.set_ylim(0, C_TOP)
    ax.set_yticks(np.arange(0, C_TOP, 5))        # no tick on the broken edge
    axt.set_ylim(end - 1.8, end + 2.2)
    axt.set_yticks([round(end)])

    ax.axhline(vs_base(0.50), color="#C8821A", ls="--", lw=1.7)
    ax.text(0.4, vs_base(0.50) + 0.7, "WHO interim target (-50% vs 2015)",
            color="#C8821A", fontsize=13.5, ha="left")
    axt.axhline(end, color="#B30000", ls="--", lw=1.7)
    axt.text(0.4, end + 0.7, "End-TB 2035 target (-90% vs 2015)",
             color="#B30000", fontsize=13.5, ha="left")

    ax.set_xlabel("Additional cascade-loop closure from 2025\n"
                  "(% of the residual loss removed at each loop driver)", fontsize=16.5)
    ax.set_ylabel("2035 incidence reduction vs Baseline (%)", fontsize=16.5)
    ax.yaxis.set_label_coords(-0.075, 0.62)
    ax.set_xticks([0, 25, 50, 75])
    ax.spines["top"].set_visible(False)
    axt.spines["bottom"].set_visible(False)
    axt.tick_params(bottom=False, labelbottom=False)

    ax2 = ax.twinx()
    ax2.set_ylim(vs_2015(0), vs_2015(C_TOP))
    ax2.set_ylabel("2035 incidence reduction vs 2015 (%)", fontsize=16.5)
    ax2.yaxis.set_label_coords(1.075, 0.62)
    ax2.spines["top"].set_visible(False)
    axt2 = axt.twinx()
    axt2.set_ylim(vs_2015(end - 1.8), vs_2015(end + 2.2))
    axt2.set_yticks([90])
    axt2.spines["bottom"].set_visible(False)

    h, l = ax.get_legend_handles_labels()
    fig.subplots_adjust(left=0.095, right=0.90, top=0.965, bottom=0.265)
    fig.legend(h, l, fontsize=13, ncol=3, frameon=False, loc="upper left",
               bbox_to_anchor=(0.095, 0.145), columnspacing=1.6, handlelength=1.8)

    # a squiggle on each spine for the interval the break leaves out
    bt, bb = axt.get_position(), ax.get_position()
    ymid = (bt.y0 + bb.y1) / 2
    t = np.linspace(0, 1, 200)
    for xc in (bb.x0, bb.x1):
        fig.add_artist(mpl.lines.Line2D(
            xc - 0.011 + t * 0.022, ymid + 0.0055 * np.sin(2 * np.pi * 2.5 * t),
            transform=fig.transFigure, color="#333333", lw=1.5,
            solid_capstyle="round", zorder=20))

    p = OUT / "FigureS2_loop_closure_sensitivity.png"
    plt.savefig(p, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close()
    print("saved", p)
    print("\n2035 incidence reduction vs Baseline (%), by additional loop closure")
    print(pts.pivot(index="scenario", columns="f", values="mean")
          .loc[ORDER].round(1).to_string())


if __name__ == "__main__":
    main()
