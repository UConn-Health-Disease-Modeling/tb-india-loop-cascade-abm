"""[2026-09-30 round-2 production run — experiment d892bdf7-18bd-f111-92f1-000d3af5294c, loose-cascade calibration with the secular infectivity decline]

Figure 2B (cascade benchmarks, 2015-2024 cohort) and all Table-2 values from the 2026-09 production run.

One pass over the 100 baseline replicates of experiment 3cd8bdf6 produces both:

  figures/reproduced/Figure2B_cascade_benchmarks.png
  analysis/derived/table2_values.json

Panels mirror the published Figure 2B: the distribution of diagnostic tests per
episode, mean tests per episode by loop group against the Subbaraman targets,
and the loop-group composition against the same targets.

Table-2 quantities, as defined in the manuscript: loop-group fractions
(mean + IQR across replicates), the three delay phases, and the excess-resource
metrics (renamed from "wasted" per reviewer comment #62).
"""
from __future__ import annotations
import io, json, sqlite3, sys
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).resolve().parent))
from extract1_baseline_episodes import process_rep, EXP, COHORT, DUMP as EPISODES  # noqa: E402

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
OUT = HERE.parent / "figures" / "reproduced"
OUT.mkdir(parents=True, exist_ok=True)
CACHE = HERE / "baseline_episodes_full_cache.csv.gz"

GROUPS = ["NoLooponly", "TestLooponly", "TxFailonly", "TestTxfail"]
SHOW = ["NoLooponly", "TestLooponly", "TestTxfail"]
GLAB = {"NoLooponly": "No loop", "TestLooponly": "Diagnostic loop only",
        "TestTxfail": "Any treatment failure"}
TARGET_TESTS = {"NoLooponly": 1.0, "TestLooponly": 3.0, "TestTxfail": 4.5}
TARGET_FRAC = {"NoLooponly": 0.20, "TestLooponly": 0.50, "TestTxfail": 0.30}
COLS = {"NoLooponly": "#DEEBF7", "TestLooponly": "#9ECAE1", "TestTxfail": "#3182BD"}


def collect():
    sims = pd.read_sql_query("SELECT SIM_ID FROM results WHERE Scenario_Code=0",
                             sqlite3.connect(str(EXP / "results.db"))).SIM_ID.astype(str)
    keep = []
    for rep, s in enumerate(sims):
        rf = EXP / s / "ReportEventRecorder.csv"
        if not rf.exists():
            continue
        e, _ = process_rep(rf)
        e["rep"] = rep
        keep.append(e[["rep", "group", "subclin", "patient", "hs", "n_tests",
                       "n_testneg", "n_testdef", "n_txfail", "eff_days",
                       "seg_subclin", "seg_patient", "seg_dxloop", "seg_txfail",
                       "seg_ontx"]])
        if (rep + 1) % 25 == 0:
            print(f"  {rep+1} reps", flush=True)
    d = pd.concat(keep, ignore_index=True)
    d.to_csv(CACHE, index=False, compression="gzip")
    return d


def q(s, p):
    return float(np.percentile(s, p)) if len(s) else float("nan")


def main():
    # 2026-10-05: benchmarks are evaluated on the same cohort as every other
    # trajectory result - episodes activating in 2015-2024, loop group per
    # episode - with treatment-failure-only episodes merged into the
    # any-treatment-failure stratum for the three displayed groups (the
    # Subbaraman stage "treatment failure" does not distinguish them).
    d_all = pd.read_csv(EPISODES)
    d = d_all[(d_all.act_year >= COHORT[0]) & (d_all.act_year < COHORT[1])].copy()
    d["group4"] = d["group"]
    d["group"] = d["group"].replace({"TxFailonly": "TestTxfail"})
    n_rep = d.rep.nunique()
    n_epi = len(d)
    print(f"{n_epi:,} cohort episodes ({COHORT[0]:.0f}-{COHORT[1]-1:.0f}) across {n_rep} replicates")

    # ── Figure 2B ────────────────────────────────────────────────────────
    fig, axes = plt.subplots(1, 3, figsize=(16.5, 5))
    ax = axes[0]
    v = d.n_tests.values
    ax.hist(v, bins=np.arange(-0.5, 16.5, 1), density=True, color="#6BAED6",
            edgecolor="white", label="Simulated")
    ax.axvline(v.mean(), color="#B30000", ls="--", lw=1.6)
    ax.axvline(np.median(v), color="#B30000", ls=":", lw=1.6)
    ax.text(0.97, 0.95, f"Mean = {v.mean():.1f}\nMedian = {np.median(v):.0f}\n"
                        f"SD = {v.std():.1f}", transform=ax.transAxes, ha="right",
            va="top", fontsize=14.9,
            bbox=dict(boxstyle="round,pad=0.4", fc="white", ec="#BBBBBB"))
    ax.set_xlabel("Number of diagnostic tests", fontsize=16.2)
    ax.set_ylabel("Density", fontsize=16.2)
    ax.set_title("B1. Tests per active TB episode", fontsize=17.6, fontweight="bold", loc="left")
    ax.set_xlim(-0.5, 16)

    ax = axes[1]
    # the published benchmark counts tests per confirmed case, so episodes
    # that never reached a test are outside its denominator; including them
    # put the no-loop group below one test per case, which the benchmark
    # cannot mean (reviewer comment)
    means = [d.loc[(d.group == g) & (d.n_tests > 0), "n_tests"].mean() for g in SHOW]
    ax.bar(range(len(SHOW)), means, color=[COLS[g] for g in SHOW], edgecolor="#2171B5")
    for i, (g, m) in enumerate(zip(SHOW, means)):
        ax.text(i, m * 0.5, f"{m:.1f}", ha="center", va="center", fontsize=16.2,
                fontweight="bold", color="#1A1A1A" if m < 1.5 else "white")
        ax.plot([i], [TARGET_TESTS[g]], "D", ms=10, mfc="white", mec="#111111", mew=1.4)
    ax.set_xticks(range(len(SHOW)))
    ax.set_xticklabels([GLAB[g] for g in SHOW], fontsize=14.9, rotation=18, ha="right")
    ax.set_ylabel("Tests per tested episode", fontsize=16.2)
    ax.set_title("B2. Mean tests per tested episode", fontsize=17.6,
                 fontweight="bold", loc="left")
    ax.plot([], [], "D", ms=9, mfc="white", mec="#111111", label="Published benchmark [3]")
    ax.legend(fontsize=13.5, frameon=False, loc="upper left")

    ax = axes[2]
    fr = [(d.group == g).mean() for g in SHOW]
    ax.bar(range(len(SHOW)), fr, color=[COLS[g] for g in SHOW], edgecolor="#2171B5")
    for i, (g, f) in enumerate(zip(SHOW, fr)):
        ax.text(i, f + 0.015, f"{f:.2f}", ha="center", va="bottom", fontsize=16.2,
                fontweight="bold", color="#1A1A1A")
        ax.plot([i], [TARGET_FRAC[g]], "D", ms=10, mfc="white", mec="#111111", mew=1.4)
    ax.set_xticks(range(len(SHOW)))
    ax.set_xticklabels([GLAB[g] for g in SHOW], fontsize=14.9, rotation=18, ha="right")
    ax.set_ylabel("Fraction of active TB episodes", fontsize=16.2)
    ax.set_title("B3. Loop-group composition", fontsize=17.6, fontweight="bold", loc="left")
    ax.set_ylim(0, 0.62)
    ax.plot([], [], "D", ms=9, mfc="white", mec="#111111", label="Published benchmark [3]")
    ax.legend(fontsize=13.5, frameon=False, loc="upper left")

    for a in axes:
        a.grid(axis="y", alpha=0.25)
        for sp in ("top", "right"):
            a.spines[sp].set_visible(False)
    fig.suptitle(f"Care-cascade benchmarks (n = {n_epi:,} active TB episodes activating "
                 f"{COHORT[0]:.0f}-{COHORT[1]-1:.0f})", fontsize=18.9, fontweight="bold", y=1.04)
    plt.tight_layout()
    # B2's title is wider than its axes and ran into B3's y-axis label
    fig.subplots_adjust(wspace=0.34)
    p = OUT / "Figure2B_cascade_benchmarks.png"
    plt.savefig(p, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close()
    print("saved", p)

    # ── Table 2 values ───────────────────────────────────────────────────
    d["group"] = d["group4"]
    per_rep = d.groupby("rep").group.value_counts(normalize=True).unstack(fill_value=0)
    out = {"n_episodes": int(n_epi), "n_reps": int(n_rep), "loop_fractions": {}}
    lab = {"NoLooponly": "NoLoop only", "TestLooponly": "Test loop only",
           "TxFailonly": "Tx fail only", "TestTxfail": "Test + Tx fail"}
    for g in GROUPS:
        s = per_rep.get(g, pd.Series(dtype=float))
        out["loop_fractions"][lab[g]] = {"mean": float(s.mean()),
                                         "iqr_lo": q(s, 25), "iqr_hi": q(s, 75)}
    anyfail = per_rep.get("TxFailonly", 0) + per_rep.get("TestTxfail", 0)
    out["loop_fractions"]["_any_treatment_failure"] = {"mean": float(anyfail.mean())}

    def block(s):
        s = pd.Series(s).dropna()
        return {"mean": float(s.mean()), "median": float(s.median()),
                "iqr_lo": q(s, 25), "iqr_hi": q(s, 75), "n": int(len(s))}

    hs_dx = d.loc[d.group.isin(["TestLooponly", "TestTxfail"]), "hs"]
    hs_tx = d.loc[d.group.isin(["TxFailonly", "TestTxfail"]), "hs"]
    out["delays"] = {"subclinical": block(d.subclin), "patient": block(d.patient),
                     "hs_all": block(d.hs), "hs_falseneg_loop": block(hs_dx),
                     "hs_txfail_loop": block(hs_tx),
                     "total_infectious": block(d.subclin + d.patient + d.hs)}
    out["effective_infectious_days"] = {lab[g]: block(d.loc[d.group == g, "eff_days"])
                                        for g in GROUPS}
    out["effective_infectious_days"]["ALL"] = block(d.eff_days)

    noloop_mean_total = float((d.loc[d.group == "NoLooponly",
                                     ["seg_subclin", "seg_patient", "seg_dxloop",
                                      "seg_txfail", "seg_ontx"]].sum(axis=1)).mean())
    total_days = d[["seg_subclin", "seg_patient", "seg_dxloop", "seg_txfail",
                    "seg_ontx"]].sum(axis=1)
    excess_pd = np.maximum(0.0, total_days - noloop_mean_total)
    out["excess"] = {
        "excess_dx_tests_per_ep": block(d.n_testneg + d.n_testdef),
        "excess_dx_tests_cohort_total": float((d.n_testneg + d.n_testdef).sum()),
        "excess_drugmonths_per_ep": block(6.0 * d.n_txfail),
        "excess_drugmonths_cohort_total": float((6.0 * d.n_txfail).sum()),
        "excess_persondays_per_ep": block(excess_pd),
        "excess_personyears_cohort_total": float(excess_pd.sum() / 365.0),
        "noloop_mean_total_days": noloop_mean_total,
    }
    (DATA / "table2_values.json").write_text(json.dumps(out, indent=1), encoding="utf-8")
    print("wrote analysis/derived/table2_values.json")

    print("\nloop fractions:", {k: round(v["mean"], 3)
                                for k, v in out["loop_fractions"].items()})
    for k, v in out["delays"].items():
        print(f"  {k:18s} mean {v['mean']:6.1f}  median {v['median']:6.1f} "
              f"[{v['iqr_lo']:.0f}-{v['iqr_hi']:.0f}]  n={v['n']:,}")
    for k, v in out["excess"].items():
        print(f"  {k:32s} {v if not isinstance(v, dict) else round(v['mean'], 2)}")


if __name__ == "__main__":
    main()
