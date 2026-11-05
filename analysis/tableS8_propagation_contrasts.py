"""Do the paper's contrasts survive parameter uncertainty?

Every estimate is formed twice over: within a parameter set across its
replicates, using the estimator the main text uses, and then across the twenty
parameter sets.  The second spread is the one the main text does not currently
report.

Contrasts between scenarios are computed WITHIN each parameter set before being
summarised across sets.  That matters: a difference can be stable even when the
two quantities being differenced are not, because the parameter set is common
to both arms, and the comparison the paper makes is of that kind.  Reading two
marginal intervals against each other would hide it.

Writes analysis/derived/propagation_contrasts.csv
"""
from __future__ import annotations
import io
import sys
from pathlib import Path

import numpy as np
import pandas as pd

if sys.stdout.encoding and sys.stdout.encoding.lower() != "utf-8":
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")

HERE = Path(__file__).resolve().parent
DATA = HERE / "derived"
SRC = DATA / "burden_by_parameter_set.csv.gz"
OUT = DATA / "propagation_contrasts.csv"
MEASURES = [("inc100k", "incidence"), ("prevall100k", "prevalence"),
            ("mort100k", "mortality")]
ARMS = ["ACF", "Dx", "Tx", "ACF+Dx", "ACF+Tx", "Dx+Tx", "ACF+Dx+Tx"]


def q(v):
    v = np.asarray(v, float)
    v = v[~np.isnan(v)]
    if not len(v):
        return np.nan, np.nan, np.nan
    return np.median(v), np.percentile(v, 2.5), np.percentile(v, 97.5)


def main():
    if not SRC.exists():
        raise SystemExit(f"{SRC} not found - run extract5_parameter_propagation.py first")
    d = pd.read_csv(SRC)
    # a replicate in which TB died out carries no information about an epidemic
    dead = d.groupby(["param_set", "rep"]).inc100k.max().pipe(lambda s: s[s == 0].index)
    if len(dead):
        d = d[~pd.MultiIndex.from_frame(d[["param_set", "rep"]]).isin(dead)]
        print(f"dropped {len(dead)} degenerate replicate(s)")
    sets = sorted(d.param_set.unique())
    print(f"{len(sets)} parameter sets, {d.rep.nunique()} replicates each\n")

    recs = []

    def add(name, per_set):
        m, lo, hi = q(list(per_set.values()))
        recs.append({"quantity": name, "median": m, "lo2.5": lo, "hi97.5": hi,
                     "n_sets": len(per_set)})
        return m, lo, hi

    # ---- reduction vs 2015, ratio of means within a set -----------------
    print("2035 reduction against 2015 (%), median [2.5-97.5] across parameter sets")
    print("  %-14s%-12s%10s%22s" % ("scenario", "measure", "median", "across sets"))
    for arm in ["Baseline", "ACF", "Tx", "ACF+Dx+Tx"]:
        for ch, lab in MEASURES:
            per = {}
            for p in sets:
                a = d[(d.param_set == p) & (d.scenario == arm) & (d.yr == 2035)][ch]
                b = d[(d.param_set == p) & (d.scenario == "Baseline") & (d.yr == 2015)][ch]
                if len(a) and len(b) and b.mean() > 0:
                    per[p] = (1 - a.mean() / b.mean()) * 100
            m, lo, hi = add(f"{arm} vs2015 {lab}", per)
            print("  %-14s%-12s%9.1f%%%22s" % (arm, lab, m, "%.1f - %.1f" % (lo, hi)))

    # ---- reduction vs the 2035 baseline, paired by seed within a set ----
    print("\n2035 reduction against the baseline (%), paired within a parameter set")
    print("  %-14s%-12s%10s%22s" % ("scenario", "measure", "median", "across sets"))
    for arm in ARMS:
        for ch, lab in MEASURES:
            per = {}
            for p in sets:
                a = d[(d.param_set == p) & (d.scenario == arm) & (d.yr == 2035)] \
                    .set_index("rep")[ch]
                b = d[(d.param_set == p) & (d.scenario == "Baseline") & (d.yr == 2035)] \
                    .set_index("rep")[ch]
                k = a.index.intersection(b.index)
                if len(k):
                    per[p] = float(((1 - a[k] / b[k]) * 100).mean())
            m, lo, hi = add(f"{arm} vsbase {lab}", per)
            print("  %-14s%-12s%9.1f%%%22s" % (arm, lab, m, "%.1f - %.1f" % (lo, hi)))

    # ---- Tx minus ACF, formed inside each parameter set -----------------
    print("\nTreatment support minus active case finding (percentage points),")
    print("differenced within each parameter set")
    print("  %-12s%10s%22s%14s" % ("measure", "median", "across sets", "sets with Tx>ACF"))
    for ch, lab in MEASURES:
        per, wins = {}, 0
        for p in sets:
            vals = {}
            for arm in ("Tx", "ACF"):
                a = d[(d.param_set == p) & (d.scenario == arm) & (d.yr == 2035)] \
                    .set_index("rep")[ch]
                b = d[(d.param_set == p) & (d.scenario == "Baseline") & (d.yr == 2035)] \
                    .set_index("rep")[ch]
                k = a.index.intersection(b.index)
                vals[arm] = float(((1 - a[k] / b[k]) * 100).mean()) if len(k) else np.nan
            diff = vals["Tx"] - vals["ACF"]
            if not np.isnan(diff):
                per[p] = diff
                wins += diff > 0
        m, lo, hi = add(f"Tx-ACF {lab}", per)
        print("  %-12s%9.1f%s%22s%10d/%d" % (
            lab, m, "pp", "%.1f - %.1f" % (lo, hi), wins, len(per)))

    # ---- duration proxy and what is left behind ------------------------
    print("\nprevalence-to-incidence ratio (years)")
    for arm, yr in (("Baseline", 2015), ("Baseline", 2035), ("ACF+Dx+Tx", 2035)):
        per = {}
        for p in sets:
            g = d[(d.param_set == p) & (d.scenario == arm) & (d.yr == yr)]
            if len(g) and g.inc100k.mean() > 0:
                per[p] = g.prevall100k.mean() / g.inc100k.mean()
        m, lo, hi = add(f"P/I {arm} {yr}", per)
        print("  %-14s%-6d%9.2f%22s" % (arm, yr, m, "%.2f - %.2f" % (lo, hi)))

    print("\npresymptomatic share of active TB prevalence in 2035 (%)")
    for arm in ("Baseline", "ACF+Dx+Tx"):
        per = {}
        for p in sets:
            g = d[(d.param_set == p) & (d.scenario == arm) & (d.yr == 2035)]
            if len(g) and g.prevall100k.mean() > 0:
                per[p] = g.presym100k.mean() / g.prevall100k.mean() * 100
        m, lo, hi = add(f"presym share {arm}", per)
        print("  %-14s%9.0f%%%22s" % (arm, m, "%.0f - %.0f" % (lo, hi)))

    pd.DataFrame(recs).to_csv(OUT, index=False)
    print(f"\nwrote {OUT}  ({len(recs)} quantities)")


if __name__ == "__main__":
    main()
