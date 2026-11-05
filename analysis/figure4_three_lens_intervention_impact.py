"""Figure 4: three-lens intervention impact by care-cascade loop group.

Reads scenario_episodes.csv.gz (extract2_scenario_episodes.py, production run,
experiment d892bdf7-18bd-f111-92f1-000d3af5294c).

Episode durations come from the model's own event timestamps, not from an
assumption-based reconstruction with fixed times per test or treatment
failure.  The figure draws the loop-event lens (D) and the delay-phase lens
(B).  The loop-group lens (A) and the resource lens (C) are computed into
figure4_values.csv but not drawn: lens A moves with the reclassification of
patients between loop groups rather than with any per-patient effect, so a
group's burden can rise under an intervention purely because patients move
into it.  Each cell carries a paired-replicate 95% CI (same seeds across
scenarios), reported in the caption.

Uses episodes activating from 2025 (when the interventions switch on in the
simulations) so that pre-intervention episodes do not dilute the contrast.

Writes figures/reproduced/Figure4_three_lens_intervention_impact.png,
analysis/derived/figure4_values.csv and analysis/derived/figure4_loop_interval.csv
(untreated diagnostic-loop days per loop event, quoted in the text).
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
# these figures are drawn wide and printed at column width, so the
# defaults that matplotlib applies to anything without an explicit fontsize
# (tick labels above all) have to be raised along with the explicit sizes.
mpl.rcParams.update({"font.size": 14, "axes.titlesize": 17, "axes.labelsize": 15,
                     "xtick.labelsize": 13.5, "ytick.labelsize": 13.5,
                     "legend.fontsize": 13.5, "figure.titlesize": 18})

HERE = Path(__file__).resolve().parent
DATA = HERE / "derived"
OUT = HERE.parent / "figures" / "reproduced"
OUT.mkdir(parents=True, exist_ok=True)

WINDOW = "post2025"
ORDER = ["Baseline", "ACF", "Dx", "Tx", "ACF+Dx", "ACF+Tx", "Dx+Tx", "ACF+Dx+Tx"]
SLAB = {"Baseline": "Baseline", "ACF": "Active case\nfinding\n(ACF)", "Dx": "Diagnostic\nimprovement\n(Dx)",
        "Tx": "Treatment\nsupport\n(Tx)", "ACF+Dx": "ACF +\nDx", "ACF+Tx": "ACF +\nTx",
        "Dx+Tx": "Dx +\nTx", "ACF+Dx+Tx": "ACF + Dx\n+ Tx"}
GLAB = {"NoLooponly": "No loop", "TestLooponly": "Diagnostic loop only",
        "TxFailonly": "Treatment failure only",
        "TestTxfail": "Diagnostic loop + treatment failure"}
NOLOOP_MIN_DAYS = 222.93   # baseline mean total infectious days of a single-pass episode


NREP = None   # set by build(); the production run is 98 replicates per scenario


def build():
    global NREP
    d = pd.read_csv(DATA / "scenario_episodes.csv.gz")
    NREP = int(d.rep.nunique())
    d = d[d.window == WINDOW].copy()
    d["hs_d"] = d.dxloop_d + d.txfail_d + d.ontx_d
    d["excess_tests"] = d.testneg + d.testdef
    d["excess_drugmonths"] = d.txfail_n * 6.0
    d["excess_persondays"] = np.maximum(0.0, d.total_d - d.n_epi * NOLOOP_MIN_DAYS)

    # ---- lens A: transmission burden by loop group ------------------
    a = d.pivot_table(index=["rep", "scenario"], columns="group", values="eff_d",
                      aggfunc="sum").reset_index()
    # ---- lens B: infectious person-time by delay phase ---------------
    b = d.groupby(["rep", "scenario"])[["subclin_d", "patient_d", "hs_d"]].sum().reset_index()
    # ---- lens C: excess resources -----------------------------------
    c = d.groupby(["rep", "scenario"])[["excess_tests", "excess_drugmonths",
                                        "excess_persondays"]].sum().reset_index()
    # ---- lens D: what actually generates the loops -------------------
    # Loops account for most of the transmission-weighted infectious time, but
    # they are produced by three separate failures and no single intervention
    # touches more than one of them.  This lens is what shows that.
    e = d.groupby(["rep", "scenario"])[["testneg", "testdef", "txfail_n",
                                        "dxloop_d", "txfail_d", "ontx_d"]].sum().reset_index()

    rows = []
    for lens, tab, cols, labels in [
            ("A", a, ["NoLooponly", "TestLooponly", "TxFailonly", "TestTxfail"],
             [GLAB[g] for g in ["NoLooponly", "TestLooponly", "TxFailonly", "TestTxfail"]]),
            ("B", b, ["subclin_d", "patient_d", "hs_d"],
             ["Subclinical", "Patient delay", "Health-system delay"]),
            ("C", c, ["excess_tests", "excess_drugmonths", "excess_persondays"],
             ["Excess diagnostic tests", "Excess drug-months", "Excess infectious person-time"]),
            ("D", e, ["testneg", "testdef", "txfail_n", "dxloop_d", "txfail_d", "ontx_d"],
             ["False-negative tests", "Pre-treatment defaults", "Treatment failures",
              "Diagnostic loop, untreated", "Failed treatment, unresolved",
              "On effective treatment"])]:
        base = tab[tab.scenario == "Baseline"].set_index("rep")
        for sc in ORDER:
            cur = tab[tab.scenario == sc].set_index("rep")
            for col, lab in zip(cols, labels):
                paired = (1 - cur[col] / base[col]) * 100        # paired by seed
                paired = paired.replace([np.inf, -np.inf], np.nan).dropna()
                m = paired.mean()
                ci = paired.sem() * 1.984
                rows.append(dict(lens=lens, row=lab, scenario=sc, mean=m, ci=ci))
    return pd.DataFrame(rows)


def plot(v: pd.DataFrame):
    # lenses A (loop group) and C (resource use) are computed but no longer drawn:
    # A moved with reclassification rather than with any per-patient effect, and
    # C's drug-months row is identical to D's treatment-failure row.
    _events = ["False-negative tests", "Pre-treatment defaults", "Treatment failures"]
    _time = ["Diagnostic loop, untreated", "Failed treatment, unresolved",
             "On effective treatment"]
    lenses = [("D", "(A) Which loop driver is removed — number of loop events", _events),
              ("B", "(B) Where in the timeline — unweighted active TB person-time by delay phase", None),
              ("D", "(C) Health-system delay of (B) by state; time on treatment included",
               _time)]
    heights = [len(r) if r else v[v.lens == l].row.nunique() for l, _, r in lenses]
    fig, axes = plt.subplots(3, 1, figsize=(14.5, 10.0),
                             gridspec_kw={"height_ratios": heights, "hspace": 0.38})
    cmap = mpl.colors.LinearSegmentedColormap.from_list(
        "RdYlGn_upper", plt.get_cmap("RdYlGn")(np.linspace(0.5, 1.0, 256)))
    cmap.set_under("#F4A582")          # a value below 0 would show up orange
    norm = mpl.colors.Normalize(vmin=0, vmax=100)

    for k, (ax, (lens, title, keep)) in enumerate(zip(axes, lenses)):
        sub = v[v.lens == lens]
        rows = list(keep) if keep else list(dict.fromkeys(sub.row))
        M = np.array([[sub[(sub.row == r) & (sub.scenario == s)]["mean"].iloc[0]
                       for s in ORDER] for r in rows])
        ax.imshow(np.where(np.abs(M) < 0.5, 0.0, M), cmap=cmap, norm=norm,
                  aspect="auto")
        for i in range(len(rows)):
            for j in range(len(ORDER)):
                val = M[i, j]
                # a leading "+" on a reduction reads as an increase, so the
                # direction is spelled out on the cells that go the other way
                _txt = "0%" if abs(val) < 0.5 else f"{val:.0f}%"
                ax.text(j, i, _txt, ha="center", va="center",
                        fontsize=15.5, fontweight="bold",
                        color="#111111" if abs(val) < 70 else "white")
        ax.set_yticks(range(len(rows)))
        ax.set_yticklabels(rows, fontsize=15.5)
        ax.set_xticks(range(len(ORDER)))
        ax.set_xticklabels([SLAB[s] for s in ORDER] if k == len(lenses) - 1
                           else [""] * len(ORDER),
                           fontsize=14.9)
        ax.set_title(title, fontsize=17.6, fontweight="bold", loc="left")
        ax.set_xticks(np.arange(-.5, len(ORDER), 1), minor=True)
        ax.set_yticks(np.arange(-.5, len(rows), 1), minor=True)
        ax.grid(which="minor", color="white", lw=2)
        ax.tick_params(which="minor", length=0)

    sm = mpl.cm.ScalarMappable(cmap=cmap, norm=norm)
    cb = fig.colorbar(sm, ax=axes, fraction=0.025, pad=0.02)
    cb.set_label("Reduction vs Baseline (%)", fontsize=16.2)

    p = OUT / "Figure4_three_lens_intervention_impact.png"
    plt.savefig(p, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close()
    print("saved", p)


def loop_interval() -> pd.DataFrame:
    """Untreated diagnostic-loop days per loop event, by scenario.

    Diagnostic-loop time is summed over the episodes activating in WINDOW and over
    replicates, then divided by the number of loop events (false-negative tests
    plus test defaults) in the same episodes.  The Results and Discussion quote it
    as the interval to repeat testing (Baseline 63 days, ACF 17, Dx 61).
    """
    d = pd.read_csv(DATA / "scenario_episodes.csv.gz")
    s = d[d.window == WINDOW].groupby("scenario")[["dxloop_d", "testneg", "testdef"]].sum()
    out = (s.dxloop_d / (s.testneg + s.testdef)).rename("days_per_loop_event")
    return out.reindex(ORDER).reset_index()


if __name__ == "__main__":
    v = build()
    v.to_csv(DATA / "figure4_values.csv", index=False)
    gap = loop_interval()
    gap.to_csv(DATA / "figure4_loop_interval.csv", index=False)
    print(gap.round(1).to_string(index=False))
    plot(v)
    # the figure draws lenses B and D only, so the figure note quotes their maximum
    print("max paired 95pct CI half-width in the drawn cells (lenses B, D): {:.2f} pp"
          .format(v[v.lens.isin(["B", "D"])].ci.max()))
    print(v[v.lens == "D"].pivot(index="row", columns="scenario", values="mean")
          [ORDER].round(0).to_string())
