"""Diagnose Fig 4B no-TRAK vs TRAK2: what is the phenotype actually?

1. Per-condition, per-plate radial profile means (to see where signal lives).
2. Re-run candidate metrics with pattern-mask variants (strip off-pattern signal).
3. Mann-Whitney as normality-free sensitivity check.
4. Per-plate Cohen's d to see if the effect is consistent across plates
   (informs whether the nested ANOVA penalty is over-harsh).
"""
from __future__ import annotations

import pathlib
import sys
import warnings

import numpy as np
import polars as pl
from scipy import stats

REPO = pathlib.Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / "replication"))

from test_single_zone_metrics import add_radial_derived  # noqa: E402
from evaluate_metrics import cohens_d  # noqa: E402

warnings.filterwarnings("ignore")

COMBINED = REPO / "replication" / "overnight_out" / "combined.csv"
SHEET = "TRAK isoform (mito)"


def main():
    df = pl.read_csv(COMBINED).filter(pl.col("sheet") == SHEET)
    for proj in ("zsum", "maxip"):
        for mask in ("crop", "pattern"):
            df = add_radial_derived(df, proj, mask)

    # -------- 1. Radial profile means per condition --------
    print("=" * 78)
    print("Radial profile means per condition (zsum, crop mask)")
    print("=" * 78)
    cols = ["zsum_crop_radial_0_2um_pct", "zsum_crop_radial_2_5um_pct",
            "zsum_crop_radial_5_10um_pct", "zsum_crop_radial_10_15um_pct",
            "zsum_crop_radial_ge15um_pct"]
    for cond in ("no TRAK", "TRAK1", "TRAK2"):
        sub = df.filter(pl.col("condition") == cond)
        print(f"\n{cond:>10s} (n={sub.height})")
        for c in cols:
            vals = sub[c].to_numpy()
            vals = vals[~np.isnan(vals)]
            print(f"  {c:<40s} {vals.mean():6.2f} ± {vals.std():5.2f}   "
                  f"(median {np.median(vals):5.2f})")

    # Same for pattern mask
    print("\n" + "=" * 78)
    print("Radial profile means per condition (zsum, pattern mask)")
    print("=" * 78)
    cols_p = ["zsum_pattern_radial_0_2um_pct", "zsum_pattern_radial_2_5um_pct",
              "zsum_pattern_radial_5_10um_pct", "zsum_pattern_radial_10_15um_pct",
              "zsum_pattern_radial_ge15um_pct"]
    for cond in ("no TRAK", "TRAK1", "TRAK2"):
        sub = df.filter(pl.col("condition") == cond)
        print(f"\n{cond:>10s} (n={sub.height})")
        for c in cols_p:
            vals = sub[c].to_numpy()
            vals = vals[~np.isnan(vals)]
            print(f"  {c:<40s} {vals.mean():6.2f} ± {vals.std():5.2f}   "
                  f"(median {np.median(vals):5.2f})")

    # -------- 2. Mann-Whitney & Welch t on top candidates, BOTH masks --------
    print("\n\n" + "=" * 78)
    print("Mann-Whitney & Welch t · no TRAK vs TRAK2 · both masks (zsum + maxip)")
    print("=" * 78)
    candidates = [
        "zsum_crop_perinuclear_5um_pct",
        "zsum_crop_peripheral_5um_pct",
        "zsum_crop_mean_dist_to_nucleus_um",
        "zsum_crop_radial_gini",
        "zsum_crop_radial_entropy",
        "zsum_crop_radial_cov",
        "zsum_crop_radial_sd_r_um",
        "zsum_crop_com_vs_pattern_offset_um",
        "zsum_crop_apical_fraction_pct",
        "zsum_pattern_perinuclear_5um_pct",
        "zsum_pattern_peripheral_5um_pct",
        "zsum_pattern_mean_dist_to_nucleus_um",
        "zsum_pattern_radial_gini",
        "zsum_pattern_radial_entropy",
        "zsum_pattern_radial_cov",
        "zsum_pattern_apical_fraction_pct",
        "zsum_pattern_com_vs_pattern_offset_um",
        "maxip_crop_radial_gini",
        "maxip_crop_radial_cov",
        "maxip_crop_apical_fraction_pct",
        "maxip_pattern_apical_fraction_pct",
        "maxip_pattern_radial_gini",
    ]
    print(f"{'Metric':<50}  {'d':>6}  {'MW p':>8}  {'Welch p':>8}  {'MW sig (raw)':>12}")
    print("-" * 95)
    results = []
    for m in candidates:
        if m not in df.columns:
            continue
        a = df.filter(pl.col("condition") == "no TRAK")[m].to_numpy()
        b = df.filter(pl.col("condition") == "TRAK2")[m].to_numpy()
        a = a[~np.isnan(a)]
        b = b[~np.isnan(b)]
        if a.size < 3 or b.size < 3:
            continue
        d = cohens_d(a, b)
        mw = stats.mannwhitneyu(a, b, alternative="two-sided").pvalue
        w = stats.ttest_ind(a, b, equal_var=False).pvalue
        sig = "***" if mw < 0.001 else "**" if mw < 0.01 else "*" if mw < 0.05 else "ns"
        results.append((m, d, mw, w, sig))
    # Sort by MW p
    results.sort(key=lambda r: r[2])
    for m, d, mw, w, sig in results:
        print(f"{m:<50}  {d:+6.2f}  {mw:8.4f}  {w:8.4f}  {sig:>12}")

    # -------- 3. Per-plate Cohen's d for the top-ranking metric --------
    print("\n\n" + "=" * 78)
    print("Per-plate Cohen's d · no TRAK vs TRAK2 · top candidates")
    print("=" * 78)
    top_metrics = [r[0] for r in results[:6]]
    for m in top_metrics:
        print(f"\n{m}:")
        by_plate_d = []
        for plate in sorted(df["plate"].unique().to_list()):
            a = df.filter((pl.col("condition") == "no TRAK") & (pl.col("plate") == plate))[m].to_numpy()
            b = df.filter((pl.col("condition") == "TRAK2") & (pl.col("plate") == plate))[m].to_numpy()
            a = a[~np.isnan(a)]
            b = b[~np.isnan(b)]
            if a.size == 0 or b.size == 0:
                print(f"  {plate:<40s}  n_no={a.size:3d}  n_T2={b.size:3d}  (skip)")
                continue
            d = cohens_d(a, b) if a.size >= 2 and b.size >= 2 else np.nan
            print(f"  {plate:<40s}  n_no={a.size:3d}  n_T2={b.size:3d}  "
                  f"mean_no={a.mean():7.3f}  mean_T2={b.mean():7.3f}  d={d:+6.2f}")
            by_plate_d.append(d)
        if by_plate_d:
            arr = np.array(by_plate_d)
            print(f"  plate d's: mean={arr.mean():+5.2f}, sign consistency: "
                  f"{(arr > 0).sum()}/{arr.size} positive")


if __name__ == "__main__":
    main()
