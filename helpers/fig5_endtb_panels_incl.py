"""Figure 5: progress against the End-TB targets on all three burden measures.

[2026-10-05 parameter uncertainty]  Lines and bars are the production run
(calibrated best fit, 98 replicates) as before.  The BANDS now carry
parameter uncertainty: experiment f8dfe706 propagated 20 parameter sets drawn
from the retained (NROY) region through all 8 scenarios at 10 replicates each
(helpers/burden_propagation.csv.gz, from extract_propagation.py).  For each
set the 10-replicate mean is formed first; the band is the 2.5th-97.5th
percentile of those set means across the 20 sets.  Panels A-C keep the
REPLICATE bands (reverted 2026-10-05 after a first version drew parameter
bands there): the across-set spread of absolute level is the 3-sigma
calibration tolerance, common to every arm and widened further by the
stratified draw that covers the secular-trend extremes, so it says nothing
about the contrasts the figure is about.  In panel D the bar keeps the
replicate bootstrap CI as the inner, thick error bar and adds the across-set
range as an outer, thin whisker; every set's reduction is measured against
that set's own 2015 baseline.  helpers/fig5_propagation_table.csv lists the
numbers (production value, replicate CI, set median and range) for the
appendix table.

The End-TB Strategy sets targets on TB deaths (-95% by 2035) and TB incidence
(-90%), both against 2015; prevalence carries no target but is shown because it
is what a cascade intervention moves most directly - prevalence is incidence
multiplied by duration, and duration is what the cascade controls.

  A  incidence trajectory      B  prevalence trajectory
  C  mortality trajectory      D  2035 reduction against 2015, all three

Panel D is drawn against the 2015 denominator so the End-TB targets sit at
fixed positions.  It uses the ratio of means with a replicate bootstrap rather
than the mean of per-replicate ratios: TB deaths are a low-count channel
(CV 0.50 against 0.17-0.19 for the other two) and the ratio estimator breaks
down on it, returning -36% for an arm whose median replicate shows no change
at all.

Prevalence is Active TB Prevalence, the channel Figure 2 calibrates to the
national prevalence survey: a survey screens on X-ray as well as symptoms
and so detects the presymptomatic phase too.  Reporting the narrower
symptomatic channel here would score the interventions on a different
quantity from the one the model was fitted to.

Writes figures/Figure_TxCure/Fig5_endtb_panels.png
"""
from __future__ import annotations
import io, sys
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt

if sys.stdout.encoding and sys.stdout.encoding.lower() != "utf-8":
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")

mpl.rcParams["font.family"] = "DejaVu Sans"
mpl.rcParams["axes.unicode_minus"] = False
mpl.rcParams.update({"font.size": 18, "axes.titlesize": 22, "axes.labelsize": 20,
                     "xtick.labelsize": 18.5, "ytick.labelsize": 18.5,
                     "legend.fontsize": 18})

HERE = Path(__file__).resolve().parent
OUT = HERE.parent / "figures" / "Figure_TxCure"
OUT.mkdir(parents=True, exist_ok=True)

ORDER = ["Baseline", "ACF", "Dx", "Tx", "ACF+Dx", "ACF+Tx", "Dx+Tx", "ACF+Dx+Tx"]
LAB = {"Baseline": "Baseline (status quo)", "ACF": "Active case finding (ACF)",
       "Dx": "Diagnostic improvement (Dx)", "Tx": "Treatment support (Tx)",
       "ACF+Dx": "ACF + Dx", "ACF+Tx": "ACF + Tx",
       "Dx+Tx": "Dx + Tx", "ACF+Dx+Tx": "ACF + Dx + Tx"}
COL = {"Baseline": "#333333", "ACF": "#4292C6", "Dx": "#41AB5D", "Tx": "#F16913",
       "ACF+Dx": "#6BAED6", "ACF+Tx": "#FD8D3C", "Dx+Tx": "#74C476",
       "ACF+Dx+Tx": "#A50F15"}
MEASURE = [("inc100k", "Incidence", "#2171B5"),
           ("prevall100k", "Prevalence", "#6A51A3"),
           ("mort100k", "Mortality", "#B30000")]
TARGET = {"inc100k": 90.0, "mort100k": 95.0}        # End-TB 2035, vs 2015
RNG = np.random.default_rng(20261004)
NBOOT = 2000
PROP = HERE / "burden_propagation.csv.gz"
BAND_ARMS = ()     # parameter bands in A-C are off; D carries the across-set range
Q_LO, Q_HI = 2.5, 97.5


def load_propagation():
    """Set-level means by scenario x set x year, degenerate replicates dropped."""
    P = pd.read_csv(PROP)
    P = P[(P.yr >= 2000) & (P.yr <= 2035)]
    dead = P.groupby(["param_set", "rep"]).inc100k.max().pipe(lambda s: s[s == 0].index)
    if len(dead):
        print(f"propagation: dropping {len(dead)} degenerate replicate(s)")
        P = P[~pd.MultiIndex.from_frame(P[["param_set", "rep"]]).isin(dead)]
    chans = [c for c, _, _ in MEASURE]
    S = P.groupby(["scenario", "param_set", "yr"])[chans].mean().reset_index()
    return S


def set_reductions(S, sc, ch, base_yr=2015):
    """Within each set: 1 - mean_2035(sc) / mean_<base_yr>(Baseline)."""
    a = S[(S.scenario == sc) & (S.yr == 2035)].set_index("param_set")[ch]
    b = S[(S.scenario == "Baseline") & (S.yr == base_yr)].set_index("param_set")[ch]
    k = a.index.intersection(b.index)
    r = (1 - a[k] / b[k]) * 100
    return r[np.isfinite(r)]


def boot_ci(num, den, nboot=NBOOT):
    """Percentile CI for 1 - mean(num)/mean(den), resampling replicates."""
    n = len(num)
    idx = RNG.integers(0, n, size=(nboot, n))
    r = (1 - num.values[idx].mean(axis=1) / den.values[idx].mean(axis=1)) * 100
    return np.percentile(r, 2.5), np.percentile(r, 97.5)


def main():
    d = pd.read_csv(HERE / "burden_by_rep_incl.csv.gz")
    d = d[(d.yr >= 2000) & (d.yr <= 2035)]
    alive = d.groupby("rep").inc100k.max()
    dead = sorted(alive.index[alive == 0])
    if dead:
        print(f"dropping degenerate replicate(s) {dead} of {alive.size}")
        d = d[~d.rep.isin(dead)]
    n_rep = d.rep.nunique()
    S = load_propagation()
    n_set = S.param_set.nunique()

    fig, axes = plt.subplots(2, 2, figsize=(18.5, 18.5))

    # ── A-C: one trajectory panel per measure ─────────────────────────
    traj = [(axes[0, 0], "inc100k", "A. TB incidence", "Incidence (per 100,000)"),
            (axes[0, 1], "prevall100k", "B. TB prevalence",
             "Active TB prevalence (per 100,000)"),
            (axes[1, 0], "mort100k", "C. TB mortality", "TB deaths (per 100,000/yr)")]
    for ax, ch, title, ylab in traj:
        for sc in ORDER:
            g = d[d.scenario == sc].groupby("yr")[ch]
            m, half = g.mean(), g.sem() * 1.984
            wide = sc in ("Baseline", "ACF+Dx+Tx")
            # each value summarises a whole year, so it is drawn at the midpoint
            # of that year; plotting at the year start makes the 2024-2025
            # segment leave the baseline a year before the interventions do
            xm = m.index.values + 0.5
            ax.plot(xm, m.values, color=COL[sc], lw=2.6 if wide else 1.7,
                    label=LAB[sc])
            ax.fill_between(xm, (m - half).values, (m + half).values,
                            color=COL[sc], alpha=0.20 if wide else 0.09, lw=0)
            if sc in BAND_ARMS:
                # parameter band: percentile across set means, by year
                # the package band is identical to the baseline band before
                # the interventions start, so it is drawn from 2024 onward only
                yrs = m.index if sc == "Baseline" else m.index[m.index >= 2024]
                gs = S[S.scenario == sc].groupby("yr")[ch]
                lo = gs.quantile(Q_LO / 100).reindex(yrs)
                hi = gs.quantile(Q_HI / 100).reindex(yrs)
                ax.fill_between(yrs.values + 0.5, lo.values, hi.values,
                                color=COL[sc], alpha=0.16, lw=0)
        if ch == "inc100k":
            ax.errorbar([2020.5], [197], yerr=[[34], [34]], fmt="D", color="#B30000",
                        ms=8, capsize=5, label="WHO 2020 estimate (197/100k)")
        ax.axvline(2025, color="#999999", ls=":", lw=1.4)
        ax.annotate("interventions\nstart (2025)", xy=(2024.4, 0.06),
                    xycoords=("data", "axes fraction"), fontsize=18,
                    color="#666666", ha="right")
        ax.set_xlabel("Year", fontsize=22)
        ax.set_ylabel(ylab, fontsize=22)
        ax.set_title(title, fontsize=23, fontweight="bold", loc="left")
        ax.set_ylim(bottom=0)
        ax.grid(alpha=0.25)
        for sp in ("top", "right"):
            ax.spines[sp].set_visible(False)
    handles, labels = axes[0, 0].get_legend_handles_labels()

    # ── D: 2035 against 2015, where the End-TB targets are defined ────
    ax = axes[1, 1]
    y15 = {ch: d[(d.scenario == "Baseline") & (d.yr == 2015)].set_index("rep")[ch]
           for ch, _, _ in MEASURE}
    y = np.arange(len(ORDER))
    h = 0.26
    for k, (ch, name, colr) in enumerate(MEASURE):
        off = (1 - k) * h
        for i, sc in enumerate(ORDER):
            v = d[(d.scenario == sc) & (d.yr == 2035)].set_index("rep")[ch]
            both = v.index.intersection(y15[ch].index)
            m = (1 - v[both].mean() / y15[ch][both].mean()) * 100
            lo, hi = boot_ci(v[both], y15[ch][both])
            r = set_reductions(S, sc, ch)
            plo, phi = np.percentile(r, Q_LO), np.percentile(r, Q_HI)
            ax.barh(y[i] + off, m, h, color=colr, edgecolor="white",
                    label=name if i == 0 else None)
            # outer thin whisker: range across parameter sets
            ax.errorbar([m], [y[i] + off], xerr=[[max(0.0, m - plo)], [max(0.0, phi - m)]],
                        fmt="none", ecolor="#444444", elinewidth=0.9, capsize=3.0,
                        capthick=0.9, zorder=3)
            # inner thick bar: replicate bootstrap CI at the production point
            ax.errorbar([m], [y[i] + off], xerr=[[m - lo], [hi - m]], fmt="none",
                        ecolor="#111111", elinewidth=3.2, capsize=0, zorder=4)
            ax.text(max(phi, hi) + 1.2, y[i] + off, f"{m:.0f}", va="center", fontsize=17)
    for ch, t in TARGET.items():
        colr = dict((c, k) for c, _, k in MEASURE)[ch]
        ax.axvline(t, color=colr, ls="--", lw=1.6)
        lab, dy, ha = (("End-TB incidence -90%", 1.005, "right")
                       if ch == "inc100k" else ("End-TB deaths -95%", 1.005, "left"))
        ax.annotate(lab, xy=(t, dy), xycoords=("data", "axes fraction"),
                    color=colr, fontsize=18, ha=ha, va="bottom")
    ax.set_yticks(y)
    ax.set_yticklabels([LAB[s] for s in ORDER], fontsize=19)
    ax.invert_yaxis()
    ax.set_xlabel("2035 reduction vs 2015 (%)", fontsize=22)
    ax.set_title("D. Progress against the End-TB denominator", fontsize=23,
                 fontweight="bold", loc="left", pad=30)
    ax.set_xlim(0, 112)
    ax.grid(axis="x", alpha=0.25)
    ax.text(0.5, -0.215, f"thick bar: 95% CI across replicates;  "
            f"thin whisker: {Q_LO:g}-{Q_HI:g} percentile across {n_set} calibrated parameter sets",
            transform=ax.transAxes, ha="center", va="top", fontsize=15,
            color="#444444")
    ax.legend(frameon=False, loc="upper center", bbox_to_anchor=(0.5, -0.125),
              ncol=3, fontsize=18, columnspacing=2.4, handlelength=1.6)
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)

    # panel D carries long scenario names on its y axis, so the columns need
    # more room between them than the trajectory panels alone would require
    fig.subplots_adjust(left=0.060, right=0.985, top=0.965, bottom=0.165,
                        wspace=0.45, hspace=0.40)
    fig.legend(handles, labels, fontsize=17.5, ncol=2, frameon=False,
               loc="upper left", bbox_to_anchor=(0.060, 0.105),
               columnspacing=1.8, handlelength=1.8, labelspacing=0.55)


    p = OUT / "Fig5_endtb_panels.png"
    plt.savefig(p, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close()
    print("saved", p)

    print(f"\n2035 vs 2015 (%): production ratio of means ({n_rep} reps) [bootstrap CI]"
          f" | across {n_set} sets: median [{Q_LO:g}-{Q_HI:g} pct]")
    recs = []
    for sc in ORDER:
        row = f"{LAB[sc]:26s}"
        for ch, name, _ in MEASURE:
            v = d[(d.scenario == sc) & (d.yr == 2035)].set_index("rep")[ch]
            both = v.index.intersection(y15[ch].index)
            m = (1 - v[both].mean() / y15[ch][both].mean()) * 100
            lo, hi = boot_ci(v[both], y15[ch][both])
            r = set_reductions(S, sc, ch)
            rb = set_reductions(S, sc, ch, base_yr=2035)       # vs the 2035 baseline
            recs.append({"scenario": sc, "measure": name, "comparison": "vs 2015",
                         "production": m, "boot_lo": lo, "boot_hi": hi,
                         "set_median": np.median(r), "set_lo": np.percentile(r, Q_LO),
                         "set_hi": np.percentile(r, Q_HI), "n_sets": len(r)})
            recs.append({"scenario": sc, "measure": name, "comparison": "vs 2035 baseline",
                         "production": np.nan, "boot_lo": np.nan, "boot_hi": np.nan,
                         "set_median": np.median(rb), "set_lo": np.percentile(rb, Q_LO),
                         "set_hi": np.percentile(rb, Q_HI), "n_sets": len(rb)})
            row += (f"  {m:5.1f} [{lo:3.0f}-{hi:3.0f}] | {np.median(r):5.1f} "
                    f"[{np.percentile(r, Q_LO):3.0f}-{np.percentile(r, Q_HI):3.0f}]")
        print(row)
    tab = pd.DataFrame(recs)
    tab.to_csv(HERE / "fig5_propagation_table.csv", index=False)
    print("wrote helpers/fig5_propagation_table.csv")


if __name__ == "__main__":
    main()
