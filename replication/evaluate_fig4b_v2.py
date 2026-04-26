"""Evaluate Fig 4B v2 with the Scheme 1 wedge-restricted polar metrics.

Loads replication/overnight_fig4b_v2_out/combined_raw.csv, joins plate/well/
condition metadata from overnight_out/combined.csv, runs nested ANOVA +
Šídák (family m=3) on the three Fig 4B pairs, and reports:
  1. Headline table for no-TRAK vs TRAK2 top-25 metrics
  2. Comparison of Y-axis (v1) vs wedge-r (v2) metrics for the same pair
  3. Ranked wedge-r-specific results for all three pairs
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

NEW_CSV = REPO / "replication" / "overnight_fig4b_v2_out" / "combined_raw.csv"
OLD_CSV = REPO / "replication" / "overnight_out" / "combined.csv"
OUT_DIR = REPO / "replication" / "overnight_fig4b_v2_out"

SHEET = "TRAK isoform (mito)"
PAIRS = [("no TRAK", "TRAK1"), ("no TRAK", "TRAK2"), ("TRAK1", "TRAK2")]
FAMILY = 3

META_COLS = {"path", "plate", "well", "sheet", "condition",
             "template_matching_score", "lateral_pixel_pitch_um",
             "zsum_bg_threshold", "maxip_bg_threshold",
             "nuc_area_um2", "nuc_solidity", "nuc_eccentricity",
             "nuc_euler_number", "nuc_n_components", "nuc_largest_area_frac",
             "zsum_total_signal", "maxip_total_signal",
             "wedge_opening_deg", "wedge_px_fraction",
             "pattern_bottom_dy_um_from_nuc", "pattern_bottom_dx_um_from_nuc",
             "pattern_top_dy_um_from_nuc", "pattern_top_dx_um_from_nuc",
             "pattern_left_dy_um_from_nuc", "pattern_left_dx_um_from_nuc",
             "pattern_right_dy_um_from_nuc", "pattern_right_dx_um_from_nuc"}


def collect(df, metric):
    sub = df.filter(pl.col(metric).is_not_null() & pl.col(metric).is_not_nan())
    conds = []
    for cn in sorted(sub["condition"].unique().to_list()):
        g = sub.filter(pl.col("condition") == cn)
        pc = {}
        for plate, grp in g.group_by("plate"):
            key = plate[0] if isinstance(plate, tuple) else plate
            pc[key] = grp[metric].to_numpy().astype(float)
        conds.append(ConditionData(name=cn, plate_cells=pc))
    return conds


def test_pair(df, metric, pair, family=3):
    conds = collect(df, metric)
    n2i = {c.name: i for i, c in enumerate(conds)}
    if pair[0] not in n2i or pair[1] not in n2i:
        return None
    i, j = n2i[pair[0]], n2i[pair[1]]
    try:
        a = nested_oneway_anova(conds)
        r = sidak_pairwise(conds, a, pairs=[(i, j)])[0]
        p = 1 - (1 - r["p_sidak"]) ** family if np.isfinite(r["p_sidak"]) else np.nan
    except Exception:
        return None
    return {"d": cohens_d(conds[i].all_cells, conds[j].all_cells), "p": p}


def main():
    if not NEW_CSV.exists():
        print(f"Not found: {NEW_CSV}"); return 1
    new = pl.read_csv(NEW_CSV)
    old = pl.read_csv(OLD_CSV).filter(pl.col("sheet") == SHEET).select(
        ["path", "plate", "well", "sheet", "condition"])
    df = new.join(old, on="path", how="left").filter(pl.col("condition").is_not_null())
    print(f"Evaluating {df.height} cells across {sorted(df['condition'].unique().to_list())}")

    metric_cols = [c for c in df.columns
                   if c not in META_COLS and c != "sheet"
                   and df[c].dtype in (pl.Float64, pl.Int64, pl.Float32)]
    print(f"Testing {len(metric_cols)} metrics × {len(PAIRS)} pairs\n")

    rows = []
    for pair in PAIRS:
        for m in metric_cols:
            r = test_pair(df, m, pair, FAMILY)
            if r and np.isfinite(r["p"]):
                rows.append({"pair": f"{pair[0]} vs {pair[1]}", "metric": m, **r})
    summary = pl.from_dicts(rows)
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    summary.write_csv(OUT_DIR / "evaluation_summary.csv")

    for pair in PAIRS:
        label = f"{pair[0]} vs {pair[1]}"
        sub = summary.filter(pl.col("pair") == label).sort("p")
        print(f"\n{'=' * 80}\n{label}  (top 25 by p)\n{'=' * 80}")
        print(f"{'Metric':<52}  {'d':>6}  {'p (Šídák m=3)':>14}  sig")
        print("-" * 85)
        for r in sub.head(25).iter_rows(named=True):
            sig = "***" if r["p"] < 0.001 else "**" if r["p"] < 0.01 else "*" if r["p"] < 0.05 else "ns"
            print(f"{r['metric']:<52}  {r['d']:+6.2f}  {r['p']:14.5f}  {sig}")

    # Focus on wedge-r metrics for no-TRAK vs TRAK2
    print(f"\n\n{'=' * 80}\nWedge-r metrics (Scheme 1) · no TRAK vs TRAK2, ranked\n{'=' * 80}")
    print(f"{'Metric':<52}  {'d':>6}  {'p':>14}  sig")
    print("-" * 85)
    wsub = summary.filter(
        (pl.col("pair") == "no TRAK vs TRAK2") &
        (pl.col("metric").str.contains("wedge_r_") | pl.col("metric").str.contains("wedge_mt_apex"))
    ).sort("p")
    for r in wsub.iter_rows(named=True):
        sig = "***" if r["p"] < 0.001 else "**" if r["p"] < 0.01 else "*" if r["p"] < 0.05 else "ns"
        print(f"{r['metric']:<52}  {r['d']:+6.2f}  {r['p']:14.5f}  {sig}")

    # Side-by-side: best Y-axis vs best wedge-r vs Mark
    print(f"\n\n{'=' * 80}\nSide-by-side: Mark / Y-axis (v1) / Wedge-r (v2) · no TRAK vs TRAK2\n{'=' * 80}")
    focus = [
        ("Mark perinuclear (zsum)", "zsum_perinuclear_5um_pct"),
        ("Mark peripheral (zsum)", "zsum_peripheral_5um_pct"),
        ("Y-Gini (maxip)", "maxip_y_gini"),
        ("Y-Gini (zsum)", "zsum_y_gini"),
        ("Y-entropy (maxip)", "maxip_y_entropy"),
        ("Wedge-r-Gini (maxip)", "maxip_wedge_r_gini"),
        ("Wedge-r-Gini (zsum)", "zsum_wedge_r_gini"),
        ("Wedge-r-entropy (maxip)", "maxip_wedge_r_entropy"),
        ("Wedge-r mean (zsum)", "zsum_wedge_r_mean_um"),
        ("Wedge-r mean (maxip)", "maxip_wedge_r_mean_um"),
        ("Wedge-r-0-20µm band (zsum)", "zsum_wedge_r_0_20um_frac_pct"),
        ("Wedge-r-20-35µm band (zsum)", "zsum_wedge_r_20_35um_frac_pct"),
        ("Wedge-r-35-55µm band (zsum)", "zsum_wedge_r_35_55um_frac_pct"),
        ("Wedge MT elongation (maxip)", "maxip_wedge_mt_apex_elongation"),
    ]
    print(f"{'Metric':<32} {'d':>6}  {'p':>14}  sig")
    print("-" * 70)
    for nice, m in focus:
        r = summary.filter((pl.col("pair") == "no TRAK vs TRAK2") &
                           (pl.col("metric") == m))
        if r.height == 0:
            print(f"{nice:<32} (not in output)")
            continue
        d = r["d"].item(); p = r["p"].item()
        sig = "***" if p < 0.001 else "**" if p < 0.01 else "*" if p < 0.05 else "ns"
        print(f"{nice:<32} {d:+6.2f}  {p:14.5f}  {sig}")


if __name__ == "__main__":
    main()
