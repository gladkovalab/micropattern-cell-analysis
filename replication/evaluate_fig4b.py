"""Evaluate the Fig 4B extended-metric rerun.

Loads replication/overnight_fig4b_out/combined_raw.csv (raw per-cell metrics),
merges plate/well/condition metadata from the original overnight_out/combined.csv
(by path), then runs nested ANOVA + Šídák (family m=3) for the three Fig 4B
pairs and prints a ranked report focused on no-TRAK vs TRAK2.
"""
from __future__ import annotations

import pathlib
import sys
import warnings

import numpy as np
import polars as pl

REPO = pathlib.Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / "replication"))

from replicate_stats import ConditionData, nested_oneway_anova, sidak_pairwise  # noqa: E402
from evaluate_metrics import cohens_d  # noqa: E402

warnings.filterwarnings("ignore")

NEW_CSV = REPO / "replication" / "overnight_fig4b_out" / "combined_raw.csv"
OLD_CSV = REPO / "replication" / "overnight_out" / "combined.csv"
OUT_DIR = REPO / "replication" / "overnight_fig4b_out"

SHEET = "TRAK isoform (mito)"
PAIRS = [("no TRAK", "TRAK1"), ("no TRAK", "TRAK2"), ("TRAK1", "TRAK2")]
FAMILY = 3

META_COLS = {"path", "plate", "well", "sheet", "condition",
             "template_matching_score", "lateral_pixel_pitch_um",
             "zsum_bg_threshold", "maxip_bg_threshold",
             "nuc_area_um2", "nuc_solidity", "nuc_eccentricity",
             "nuc_euler_number", "nuc_n_components", "nuc_largest_area_frac",
             "zsum_total_signal", "maxip_total_signal"}


def collect(df: pl.DataFrame, metric: str) -> list[ConditionData]:
    sub = df.filter(pl.col(metric).is_not_null() & pl.col(metric).is_not_nan())
    conds: list[ConditionData] = []
    for cond_name in sorted(sub["condition"].unique().to_list()):
        g = sub.filter(pl.col("condition") == cond_name)
        plate_cells: dict[str, np.ndarray] = {}
        for plate, grp in g.group_by("plate"):
            key = plate[0] if isinstance(plate, tuple) else plate
            plate_cells[key] = grp[metric].to_numpy().astype(float)
        conds.append(ConditionData(name=cond_name, plate_cells=plate_cells))
    return conds


def test_pair(df, metric, pair, family):
    conds = collect(df, metric)
    name_to_idx = {c.name: i for i, c in enumerate(conds)}
    if pair[0] not in name_to_idx or pair[1] not in name_to_idx:
        return None
    i, j = name_to_idx[pair[0]], name_to_idx[pair[1]]
    try:
        a = nested_oneway_anova(conds)
        r = sidak_pairwise(conds, a, pairs=[(i, j)])[0]
        p_raw = r["p_sidak"]
        p_sidak = 1 - (1 - p_raw) ** family if np.isfinite(p_raw) else np.nan
    except Exception:
        return None
    d = cohens_d(conds[i].all_cells, conds[j].all_cells)
    return {"d": d, "p_sidak": p_sidak, "n_i": int(conds[i].all_cells.size),
            "n_j": int(conds[j].all_cells.size)}


def main():
    if not NEW_CSV.exists():
        print(f"Not found: {NEW_CSV}"); return 1
    # Merge path-based metadata
    new = pl.read_csv(NEW_CSV)
    old = pl.read_csv(OLD_CSV).filter(pl.col("sheet") == SHEET).select(
        ["path", "plate", "well", "sheet", "condition"])
    # Build a lookup from path-basename (SMB paths may differ) — match by path directly
    df = new.join(old, on="path", how="left")
    missing = df.filter(pl.col("condition").is_null())
    if missing.height > 0:
        print(f"WARNING: {missing.height} cells have no condition match")
        print(missing.select("path").head(3))
    df = df.filter(pl.col("condition").is_not_null())
    print(f"Evaluating {df.height} cells across conditions: "
          f"{sorted(df['condition'].unique().to_list())}")
    print(f"Plates: {sorted(df['plate'].unique().to_list())}")

    metric_cols = [c for c in df.columns if c not in META_COLS and c not in ("sheet",)]
    metric_cols = [c for c in metric_cols if df[c].dtype in (pl.Float64, pl.Int64, pl.Float32)]
    print(f"Testing {len(metric_cols)} metrics × {len(PAIRS)} pairs\n")

    rows = []
    for pair in PAIRS:
        for m in metric_cols:
            r = test_pair(df, m, pair, FAMILY)
            if r is not None and np.isfinite(r["p_sidak"]):
                rows.append({"pair": f"{pair[0]} vs {pair[1]}", "metric": m, **r})
    summary = pl.from_dicts(rows)
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    summary.write_csv(OUT_DIR / "evaluation_summary.csv")

    for pair in PAIRS:
        label = f"{pair[0]} vs {pair[1]}"
        sub = summary.filter(pl.col("pair") == label).sort("p_sidak")
        print(f"\n{'=' * 80}\n{label}  (top 25 by p)\n{'=' * 80}")
        print(f"{'Metric':<55}  {'d':>6}  {'p (Šídák m=3)':>14}  sig")
        print("-" * 85)
        for r in sub.head(25).iter_rows(named=True):
            p = r["p_sidak"]
            sig = "***" if p < 0.001 else "**" if p < 0.01 else "*" if p < 0.05 else "ns"
            print(f"{r['metric']:<55}  {r['d']:+6.2f}  {p:14.5f}  {sig}")

    # Highlight the "newly significant" metrics on no-TRAK vs TRAK2 that the existing
    # overnight_out combined.csv does NOT have
    print(f"\n\n{'=' * 80}\nNewly significant metrics on no TRAK vs TRAK2 (not in original pipeline)\n{'=' * 80}")
    old_cols_lower = {c.lower() for c in old.columns}
    old_full = pl.read_csv(OLD_CSV)
    old_metric_cols = {c for c in old_full.columns if c not in ("path", "plate", "well", "sheet", "condition",
                                                                "template_matching_score", "lateral_pixel_pitch_um",
                                                                "zsum_bg_threshold", "maxip_bg_threshold")}
    ntr = summary.filter((pl.col("pair") == "no TRAK vs TRAK2") &
                         (pl.col("p_sidak") < 0.05)).sort("p_sidak")
    new_only = ntr.filter(~pl.col("metric").is_in(list(old_metric_cols)))
    print(f"{'Metric':<55}  {'d':>6}  {'p':>14}  sig")
    print("-" * 85)
    for r in new_only.iter_rows(named=True):
        p = r["p_sidak"]
        sig = "***" if p < 0.001 else "**" if p < 0.01 else "*"
        print(f"{r['metric']:<55}  {r['d']:+6.2f}  {p:14.5f}  {sig}")
    if new_only.height == 0:
        print("  (none) — best on no TRAK vs TRAK2:")
        best = summary.filter(pl.col("pair") == "no TRAK vs TRAK2").sort("p_sidak").head(10)
        for r in best.iter_rows(named=True):
            p = r["p_sidak"]
            sig = "***" if p < 0.001 else "**" if p < 0.01 else "*" if p < 0.05 else "ns"
            print(f"  {r['metric']:<55}  {r['d']:+6.2f}  {p:14.5f}  {sig}")


if __name__ == "__main__":
    main()
