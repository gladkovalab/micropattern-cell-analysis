"""Nested ANOVA + Šídák (m=3) for pattern-masked single-zone metrics on Fig 4B.

Does the pattern-mask fix carry through the plate-aware statistical framework
Mark uses in his Prism analysis?
"""
from __future__ import annotations

import pathlib
import sys
import warnings

import numpy as np
import polars as pl

REPO = pathlib.Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / "replication"))

from test_single_zone_metrics import add_radial_derived, test_metric, SHEET, PAIRS, FAMILY  # noqa: E402
from evaluate_metrics import cohens_d  # noqa: E402

warnings.filterwarnings("ignore")

COMBINED = REPO / "replication" / "overnight_out" / "combined.csv"


def main():
    df = pl.read_csv(COMBINED).filter(pl.col("sheet") == SHEET)
    for proj in ("zsum", "maxip"):
        for mask in ("crop", "pattern"):
            df = add_radial_derived(df, proj, mask)

    candidates = [
        # Mark baselines (crop)
        "zsum_crop_perinuclear_5um_pct",
        "zsum_crop_peripheral_5um_pct",
        # Mark-style metrics under pattern mask
        "zsum_pattern_perinuclear_5um_pct",
        "zsum_pattern_peripheral_5um_pct",
        "maxip_pattern_perinuclear_5um_pct",
        "maxip_pattern_peripheral_5um_pct",
        # Single-zone (pattern mask)
        "zsum_pattern_mean_dist_to_nucleus_um",
        "zsum_pattern_median_dist_to_nucleus_um",
        "zsum_pattern_q90_dist_to_nucleus_um",
        "zsum_pattern_apical_fraction_pct",
        "zsum_pattern_com_vs_pattern_offset_um",
        "zsum_pattern_com_offset_um",
        "zsum_pattern_radial_gini",
        "zsum_pattern_radial_entropy",
        "zsum_pattern_radial_cov",
        "zsum_pattern_radial_sd_r_um",
        "maxip_pattern_mean_dist_to_nucleus_um",
        "maxip_pattern_apical_fraction_pct",
        "maxip_pattern_radial_gini",
        "maxip_pattern_radial_cov",
        # Apical fraction on crop (best raw MW hit)
        "zsum_crop_apical_fraction_pct",
        "maxip_crop_apical_fraction_pct",
    ]

    rows = []
    for m in candidates:
        if m not in df.columns:
            continue
        r = test_metric(df, m)
        if r is not None:
            rows.append(r)

    for pair_label in ["no TRAK vs TRAK1", "no TRAK vs TRAK2", "TRAK1 vs TRAK2"]:
        print(f"\n\n=== Fig 4B: {pair_label} (nested ANOVA, Šídák m={FAMILY}) ===\n")
        print(f"{'Metric':<50}  {'d':>6}  {'p (Šídák)':>10}  {'sig':>4}")
        print("-" * 78)
        rows_sorted = sorted(
            rows,
            key=lambda r: r.get(f"{pair_label} p", 1.0) or 1.0,
        )
        for r in rows_sorted:
            d = r.get(f"{pair_label} d")
            p = r.get(f"{pair_label} p")
            if d is None or p is None:
                continue
            sig = "***" if p < 0.001 else "**" if p < 0.01 else "*" if p < 0.05 else "ns"
            print(f"{r['metric']:<50}  {d:+6.2f}  {p:10.4f}  {sig:>4}")


if __name__ == "__main__":
    main()
