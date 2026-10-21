"""[2026-09-30 round-2 production run — experiment d892bdf7-18bd-f111-92f1-000d3af5294c, loose-cascade calibration with the secular infectivity decline]

[2026-09-25 production re-run — experiment d892bdf7-18bd-f111-92f1-000d3af5294c, trend-calibrated parameters,
interventions from 2025]

Draw Figure 1A and Figure 1B for the 2026-09 revision.

Reads helpers/fig1_values_txc.json (produced by fig1_rebuild_0909.py) and
writes figures/Figure_TxCure/Fig1A_cascade_rates.png and
figures/Figure_TxCure/Fig1B_patient_trajectories.png
"""
from __future__ import annotations
import json, sys
from pathlib import Path
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

mpl.rcParams["font.family"] = "DejaVu Sans"
mpl.rcParams["axes.unicode_minus"] = False
# 2026-10-01: these figures are drawn wide and printed at column width, so the
# defaults that matplotlib applies to anything without an explicit fontsize
# (tick labels above all) have to be raised along with the explicit sizes.
mpl.rcParams.update({"font.size": 16, "axes.titlesize": 19.5, "axes.labelsize": 17.5,
                     "xtick.labelsize": 16, "ytick.labelsize": 16,
                     "legend.fontsize": 16, "figure.titlesize": 20})

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
OUT = ROOT / "figures" / "Figure_TxCure"
OUT.mkdir(parents=True, exist_ok=True)
V = json.loads((HERE / "fig1_values_txc.json").read_text(encoding="utf-8"))

REDS = {"NoLooponly": "#FCAE91", "TestLooponly": "#FB6A4A",
        "TxFailonly": "#DE2D26", "TestTxfail": "#A50F15"}
GLAB = {"NoLooponly": "No loop",
        "TestLooponly": "Diagnostic\nloop only",
        "TxFailonly": "Treatment\nfailure\nonly",
        "TestTxfail": "Diagnostic loop\n+ treatment\nfailure"}

# 2026-10-05: raised to match Figures 3-5 at column width
FS_T, FS_L, FS_A = 19.5, 17.5, 16


# ======================================================================
# FIGURE 1A -- cascade event rates and loop-group composition,
#              both on one denominator: active TB episodes
# ======================================================================
def fig1a():
    ev_order = [("tests_any", "Diagnostic\ntests", "#6BAED6"),
                ("test_neg", "False-negative\nresults", "#4292C6"),
                ("test_default", "Pre-treatment\ndefaults", "#2171B5"),
                ("tx_start", "Treatment\nstarts", "#BDD7E7"),
                ("tx_fail", "Treatment\nfailures", "#08519C"),
                ("self_clear", "Resolved\n(cleared/cured)", "#BDBDBD"),
                ("tb_death", "TB\ndeaths", "#404040")]

    fig, axes = plt.subplots(1, 2, figsize=(19, 7.4),
                             gridspec_kw={"width_ratios": [1.45, 1]})

    ax = axes[0]
    xs = np.arange(len(ev_order))
    means = [V["panelA_per100"][k]["mean"] for k, _, _ in ev_order]
    cis = [V["panelA_per100"][k]["ci"] for k, _, _ in ev_order]
    cols = [c for _, _, c in ev_order]
    ax.bar(xs, means, yerr=cis, color=cols, edgecolor="white", width=0.72,
           error_kw=dict(ecolor="#333333", capsize=4, capthick=1.1, lw=1.1))
    for x, m in zip(xs, means):
        ax.text(x, m + 12, f"{m:.0f}", ha="center", va="bottom",
                fontsize=FS_A, fontweight="bold")
    ax.axvline(4.5, color="#999999", ls=":", lw=1.2)
    ax.text(2.25, max(means) * 1.02, "cascade events", ha="center",
            fontsize=FS_A, style="italic", color="#555555")
    ax.text(5.5, max(means) * 1.02, "episode outcome", ha="center",
            fontsize=FS_A, style="italic", color="#555555")
    ax.set_xticks(xs)
    ax.set_xticklabels([l for _, l, _ in ev_order], fontsize=FS_A,
                       rotation=32, ha="right", rotation_mode="anchor")
    ax.set_ylabel("Events per 100 active TB episodes", fontsize=FS_L)
    ax.set_title("A. Cascade activity per 100 active TB episodes",
                 fontsize=FS_T, fontweight="bold", loc="left")
    ax.set_ylim(0, max(means) * 1.13)
    ax.grid(axis="y", alpha=0.3)
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)

    ax = axes[1]
    groups = ["NoLooponly", "TestLooponly", "TxFailonly", "TestTxfail"]
    m = [V["panelB_pct"][g]["mean"] for g in groups]
    c = [V["panelB_pct"][g]["ci"] for g in groups]
    xs = np.arange(len(groups))
    ax.bar(xs, m, yerr=c, color=[REDS[g] for g in groups], edgecolor="white",
           width=0.68, error_kw=dict(ecolor="#333333", capsize=4, capthick=1.1, lw=1.1))
    for x, v in zip(xs, m):
        ax.text(x, v + 0.9, f"{v:.1f}%", ha="center", va="bottom",
                fontsize=FS_A, fontweight="bold")
    ax.set_xticks(xs)
    ax.set_xticklabels([GLAB[g] for g in groups], fontsize=FS_A,
                       rotation=0, ha="center")
    ax.set_ylabel("% of active TB episodes", fontsize=FS_L)
    ax.set_title("B. Care-cascade trajectory composition",
                 fontsize=FS_T, fontweight="bold", loc="left")
    ax.set_ylim(0, max(m) * 1.32)
    ax.grid(axis="y", alpha=0.3)
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
    anyfail = V["panelB_pct"]["TxFailonly"]["mean"] + V["panelB_pct"]["TestTxfail"]["mean"]
    anyloop = 100 - V["panelB_pct"]["NoLooponly"]["mean"]
    ax.text(0.02, 0.95,
            f"any cascade loop {anyloop:.0f}%    |    any treatment failure {anyfail:.0f}%",
            transform=ax.transAxes, fontsize=FS_A, color="#A50F15",
            fontweight="bold", va="top")

    plt.tight_layout()
    p = OUT / "Fig1A_cascade_rates.png"
    plt.savefig(p, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close()
    print("saved", p)


# ======================================================================
# FIGURE 1B -- three real patient trajectories with transition times
# ======================================================================
def fig1b():
    D = V["group_delays_days"]
    C = {"subclin": "#F5C97A", "patient": "#ED7D31", "dxloop": "#F4A6A6",
         "txfail": "#C00000", "ontx": "#BFBFBF"}
    W = {"subclin": 0.12, "patient": 1.0, "dxloop": 1.0, "txfail": 1.0, "ontx": 0.05}
    MO = 30.44

    tracks = [
        ("NoLooponly", "Patient A\nno loop", "#B45309"),
        ("TestLooponly", "Patient B\ndiagnostic loop", "#C2410C"),
        ("TestTxfail", "Patient C\nloop + Tx failure", "#C00000"),
    ]
    seg_keys = [("m_subclin", "subclin", "Subclinical"),
                ("m_patient", "patient", "Patient delay"),
                ("m_dxloop", "dxloop", "Diagnostic loop"),
                ("m_txfail", "txfail", "Failed treatment"),
                ("m_ontx", "ontx", "On treatment → cure")]

    fig, ax = plt.subplots(figsize=(16.5, 8.2))
    BAR_H, GAP = 0.58, 1.85
    ys = [2 * GAP, GAP, 0.0]

    for (g, label, lab_col), y in zip(tracks, ys):
        d = D[g]
        x = 0.0
        marks = []
        for key, state, _ in seg_keys:
            m_days = d[key] / MO
            if m_days <= 0.01:
                continue
            ax.barh(y, m_days, left=x, height=BAR_H, color=C[state],
                    edgecolor="white", lw=1.4)
            if m_days > 1.6:
                ax.text(x + m_days / 2, y, f"{d[key]:.0f} d",
                        ha="center", va="center", fontsize=15.5, fontweight="bold",
                        color="white" if state == "txfail" else "#222222")
            x += m_days
            marks.append((state, x))
        # transition markers with cumulative time (staggered so labels never collide)
        transitions, cum = [], 0.0
        for key, state, _ in seg_keys:
            cum += d[key] / MO
            if state == "subclin":
                transitions.append(("symptoms", cum))
            elif state == "patient":
                transitions.append(("1st test", cum))
            elif state == "dxloop" and d[key] > 1:
                transitions.append(("Tx start", cum))
        for k, (name, cm) in enumerate(transitions):
            off = (0.16, 0.42, 0.68)[k % 3]
            ax.plot([cm], [y - BAR_H / 2 - 0.04], marker="v", color="#333333", ms=6)
            ax.plot([cm, cm], [y - BAR_H / 2 - 0.06, y - BAR_H / 2 - off + 0.02],
                    color="#999999", lw=0.8)
            ax.text(cm, y - BAR_H / 2 - off, f"{name} {cm:.1f}",
                    ha="center", va="top", fontsize=13.5, color="#333333")
        total_mo = sum(d[k] for k, _, _ in seg_keys) / MO
        ax.plot([total_mo], [y], marker="|", color="#333333", ms=16, mew=2)
        ax.text(-0.6, y, label, ha="right", va="center", fontsize=17.6,
                fontweight="bold", color=lab_col)
        eff = d["m_eff_days"] / MO
        ax.text(total_mo + 0.45, y + 0.20,
                f"total {total_mo:.1f} mo", ha="left", va="center", fontsize=16.2,
                color="#333333")
        ax.text(total_mo + 0.45, y - 0.16,
                f"effective {eff:.1f} mo", ha="left", va="center",
                fontsize=16.9, fontweight="bold", color="#C00000")
        # composition note
        note = {"NoLooponly": f"{d['pct_tested']:.0f}% ever tested",
                "TestLooponly": (f"{d['n_tests']:.1f} tests · {d['n_testneg']:.1f} false-neg · "
                                 f"{d['n_testdef']:.1f} defaults"),
                "TestTxfail": (f"{d['n_tests']:.1f} tests · {d['n_txfail']:.1f} failed regimens")}[g]
        ax.text(0.15, y + BAR_H / 2 + 0.10, note, ha="left", va="bottom",
                fontsize=14.2, style="italic", color="#555555")

    ax.set_xlim(-6.2, 36.5)
    ax.set_ylim(-1.30, 2 * GAP + 0.80)
    ax.set_yticks([])
    ax.set_xticks(range(0, 31, 3))
    ax.set_xlabel("Months since TB activation (presymptomatic onset)",
                  fontsize=17.6)
    ax.tick_params(axis="x", labelsize=16.2)
    for sp in ("top", "right", "left"):
        ax.spines[sp].set_visible(False)
    ax.grid(axis="x", alpha=0.25)

    leg = [Patch(facecolor=C["subclin"], label="Subclinical (~23%)"),
           Patch(facecolor=C["patient"], label="Untreated — patient delay (100%)"),
           Patch(facecolor=C["dxloop"], label="Diagnostic loop (100%)"),
           Patch(facecolor=C["txfail"], label="Failed treatment (100%)"),
           Patch(facecolor=C["ontx"], label="On effective treatment (0%)")]
    # above the axes: inside the panel the box covered Patient A's duration labels
    ax.legend(handles=leg, loc="lower center", bbox_to_anchor=(0.5, 1.01), ncol=3,
              fontsize=15.5, framealpha=0.95, title="Infectiousness state",
              title_fontsize=16.2)

    plt.tight_layout()
    p = OUT / "Fig1B_patient_trajectories.png"
    plt.savefig(p, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close()
    print("saved", p)


if __name__ == "__main__":
    fig1a()
    fig1b()
