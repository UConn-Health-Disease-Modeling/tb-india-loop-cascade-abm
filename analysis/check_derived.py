"""Check that tables re-extracted from the raw outputs match the committed ones.

After downloading the raw outputs (simulation/download_outputs.py) and re-running
the extraction scripts, run from the top folder of the repository:

    python analysis/check_derived.py

Every file in analysis/derived/ is compared with its committed version (git HEAD).
A file that was not re-written is reported as "not re-run".  Tables are compared
by content, ignoring row order, with numbers agreeing to a relative tolerance of
1e-9; text reports must match exactly.  The calibration outputs that
calib2_design_point_outputs.py writes into each downloaded experiment are compared
with calibration/wave1/ and calibration/wave2/ in the same way.

Exit status is 1 if any compared file differs.
"""
from __future__ import annotations
import io, json, subprocess, sys
from pathlib import Path
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
RTOL = 1e-9
WAVES = {"wave1": "af707f22-16bc-f111-92f1-000d3af5294c",
         "wave2": "496ded25-37bc-f111-92f1-000d3af5294c"}


def committed(rel: str) -> bytes:
    return subprocess.run(["git", "show", f"HEAD:{rel}"], cwd=ROOT,
                          capture_output=True, check=True).stdout


def read_table(raw: bytes, name: str) -> pd.DataFrame:
    comp = "gzip" if name.endswith(".gz") else None
    return pd.read_csv(io.BytesIO(raw), compression=comp)


def compare_tables(a: pd.DataFrame, b: pd.DataFrame) -> str:
    if list(a.columns) != list(b.columns):
        return f"columns differ: {list(a.columns)} vs {list(b.columns)}"
    if len(a) != len(b):
        return f"row counts differ: {len(a)} vs {len(b)}"
    cols = list(a.columns)
    a = a.sort_values(cols, kind="stable").reset_index(drop=True)
    b = b.sort_values(cols, kind="stable").reset_index(drop=True)
    num = [c for c in cols if pd.api.types.is_numeric_dtype(a[c]) and pd.api.types.is_numeric_dtype(b[c])]
    txt = [c for c in cols if c not in num]
    if txt and not a[txt].astype(str).equals(b[txt].astype(str)):
        return f"text columns differ ({', '.join(txt)})"
    if num:
        x, y = a[num].to_numpy(float), b[num].to_numpy(float)
        if not np.allclose(x, y, rtol=RTOL, atol=0, equal_nan=True):
            bad = ~np.isclose(x, y, rtol=RTOL, atol=0, equal_nan=True)
            return (f"numbers differ in {bad.sum()} cell(s), "
                    f"largest absolute difference {np.nanmax(np.abs(x - y)[bad]):.3g}")
    return ""


def compare_json(a, b, path="") -> str:
    if isinstance(a, dict) and isinstance(b, dict):
        if set(a) != set(b):
            return f"keys differ at {path or '/'}"
        for k in a:
            msg = compare_json(a[k], b[k], f"{path}/{k}")
            if msg:
                return msg
        return ""
    if isinstance(a, list) and isinstance(b, list):
        if len(a) != len(b):
            return f"list lengths differ at {path}"
        for i, (u, v) in enumerate(zip(a, b)):
            msg = compare_json(u, v, f"{path}[{i}]")
            if msg:
                return msg
        return ""
    if isinstance(a, (int, float)) and isinstance(b, (int, float)) and not isinstance(a, bool):
        return "" if np.isclose(a, b, rtol=RTOL, atol=0, equal_nan=True) else f"value differs at {path}: {a} vs {b}"
    return "" if a == b else f"value differs at {path}: {a!r} vs {b!r}"


def compare(new: bytes, old: bytes, name: str) -> str:
    if new == old:
        return ""
    if name.endswith((".csv", ".csv.gz")):
        return compare_tables(read_table(new, name), read_table(old, name))
    if name.endswith(".json"):
        return compare_json(json.loads(new), json.loads(old))
    return "text differs"


def main() -> int:
    files = subprocess.run(["git", "ls-files", "analysis/derived"], cwd=ROOT,
                           capture_output=True, text=True, check=True).stdout.split()
    checks = [(rel, ROOT / rel, rel) for rel in files]
    for wave, exp in WAVES.items():
        for name in ("design_outputs.csv", "best_fit.json"):
            checks.append((f"{exp[:8]}/calibration/{name}",
                           ROOT / exp / "calibration" / name, f"calibration/{wave}/{name}"))

    failed = 0
    for label, path, ref in checks:
        old = committed(ref)
        if not path.exists():
            status = "not re-run"
        else:
            new = path.read_bytes()
            if new == old and path == ROOT / ref:
                status = "not re-run"
            else:
                msg = compare(new, old, ref)
                status = "SAME" if not msg else f"DIFFERENT: {msg}"
                failed += bool(msg)
        print(f"{status:12.12s}  {label}" if status in ("SAME", "not re-run") else f"{status}  {label}")
    print(f"\n{failed} file(s) differ" if failed else "\nno differences in the files that were re-run")
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
