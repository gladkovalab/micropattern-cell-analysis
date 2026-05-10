"""Print nested-ANOVA + Šídák pairwise stats per sheet, reusing plot_metrics.py
helpers. Same numbers that appear as brackets in the per-sheet figures.
"""
from __future__ import annotations

import argparse
import pathlib
import sys

import numpy as np
import polars as pl

REPO = pathlib.Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / "analysis"))
from plot_metrics import (  # noqa: E402
    SHEET_CONFIG, _test_pair, _collect, _nested_oneway_anova,
    load_template_matching, load_comparisons_table, join_with_metadata,
)

METRICS = [
    "peripheral_5um_percent_total",
    "perinuclear_5um_percent_total",
    "wedge_r_ks_vs_uniform",
    "wedge_r_ks_vs_60merNoTRAK",
]


def fmt_p(p: float) -> str:
    if not np.isfinite(p):
        return "ns"
    if p < 0.001:
        return f"*** {p:.1e}"
    if p < 0.01:
        return f"**  {p:.4f}"
    if p < 0.05:
        return f"*   {p:.4f}"
    return f"ns  {p:.4f}"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--template-matching",
                    default="analysis/wedge_r_ks_out_all_denoised/by_well")
    ap.add_argument("--comparisons-xlsx",
                    default=str(REPO / "config/Comparisons_table_v3.xlsx"))
    ap.add_argument("--sheets", nargs="*",
                    help="Sheets to summarize. Default: all known except mito.")
    args = ap.parse_args()

    df = load_template_matching(pathlib.Path(args.template_matching))
    df = join_with_metadata(df, pathlib.Path(args.comparisons_xlsx))

    sheets = args.sheets or [
        "TRAK isoform (peroxisome)",
        "TRAK isoform (60mer)",
        "TRAK1 helix muts",
        "TRAK2 helix muts",
        "MAPK9 siRNA",
    ]
    for sheet in sheets:
        if sheet not in SHEET_CONFIG:
            print(f"\n=== {sheet}: not configured, skipping ===")
            continue
        cfg = SHEET_CONFIG[sheet]
        sub = df.filter(pl.col("sheet") == sheet)
        if sub.height == 0:
            print(f"\n=== {sheet}: no rows ===")
            continue
        print(f"\n=== {sheet}  (n={sub.height} cells, family m={cfg['family_m']}) ===")
        for metric in METRICS:
            if metric not in sub.columns:
                continue
            print(f"\n  {metric}")
            # Per-condition mean ± SEM
            for cond in cfg["conditions"]:
                vals = (sub.filter(pl.col("condition") == cond)[metric]
                        .to_numpy().astype(float))
                vals = vals[np.isfinite(vals)]
                if vals.size == 0:
                    print(f"    {cond:25s}: (no data)")
                    continue
                sem = vals.std(ddof=1) / max(np.sqrt(vals.size), 1)
                print(f"    {cond:25s}: n={vals.size:3d}  "
                      f"mean={vals.mean():8.4f}  sem={sem:8.4f}")
            # Pairwise (nested ANOVA + Šídák)
            for pair in cfg["pairs"]:
                p = _test_pair(sub, metric, pair, cfg["family_m"])
                print(f"    {pair[0]:>18s}  vs  {pair[1]:<18s}  →  {fmt_p(p)}")


if __name__ == "__main__":
    main()
