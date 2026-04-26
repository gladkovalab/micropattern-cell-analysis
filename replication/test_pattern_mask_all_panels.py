"""Test pattern-mask + Mark-style single-zone metrics on all reviewer-flagged panels.

Panels to test:
- Fig 4B: TRAK isoform (mito), no TRAK vs TRAK2 (already done), plus all 3 pairs
- Fig 4C / S11 D: TRAK1 helix muts, T1 wt vs mDRH, T1 mDRH vs mDRH/dSp
- Fig 4D / S11 E: TRAK2 helix muts, TRAK2 vs mDRH, mDRH vs mSpindly
- Fig 4E / S11 F: MAPK9 siRNA, multiple pairs

Compare Mark-crop, Mark-pattern, plus a few single-zone metrics under pattern mask.
"""
from __future__ import annotations

import pathlib
import sys
import warnings

import numpy as np
import polars as pl

REPO = pathlib.Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / "replication"))

from test_single_zone_metrics import add_radial_derived  # noqa: E402
from replicate_stats import ConditionData, nested_oneway_anova, sidak_pairwise  # noqa: E402
from evaluate_metrics import cohens_d  # noqa: E402

warnings.filterwarnings("ignore")

COMBINED = REPO / "replication" / "overnight_out" / "combined.csv"

SHEET_PAIRS = {
    "TRAK isoform (mito)": [("no TRAK", "TRAK1"), ("no TRAK", "TRAK2"), ("TRAK1", "TRAK2")],
    "TRAK1 helix muts": [("T1 wt", "T1 mDRH"), ("T1 mDRH", "T1 mDRH / dSp")],
    "TRAK2 helix muts": [("TRAK2", "TRAK2 mDRH"), ("TRAK2 mDRH", "TRAK2 mDRH mSpindly")],
    "MAPK9 siRNA": [
        ("ctrl ctrl", "MAPK9 ctrl"), ("ctrl ctrl", "ctrl Ars"), ("ctrl ctrl", "MAPK9 Ars"),
        ("MAPK9 ctrl", "MAPK9 Ars"), ("ctrl Ars", "MAPK9 Ars"),
    ],
}


def _collect(df: pl.DataFrame, metric: str) -> list[ConditionData]:
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


def test_pair(df: pl.DataFrame, metric: str, pair: tuple[str, str], family: int) -> dict | None:
    conds = _collect(df, metric)
    name_to_idx = {c.name: i for i, c in enumerate(conds)}
    if pair[0] not in name_to_idx or pair[1] not in name_to_idx:
        return None
    i, j = name_to_idx[pair[0]], name_to_idx[pair[1]]
    a = nested_oneway_anova(conds)
    r = sidak_pairwise(conds, a, pairs=[(i, j)])[0]
    p_raw = r["p_sidak"]
    p_sidak = 1 - (1 - p_raw) ** family if np.isfinite(p_raw) else np.nan
    d = cohens_d(conds[i].all_cells, conds[j].all_cells)
    return {"d": d, "p_sidak": p_sidak}


METRIC_SET_MITO = [
    # Crop (existing baseline)
    "zsum_crop_perinuclear_5um_pct",
    "zsum_crop_peripheral_5um_pct",
    # Pattern-masked Mark-style
    "zsum_pattern_perinuclear_5um_pct",
    "zsum_pattern_peripheral_5um_pct",
    # Pattern-masked single-zone
    "zsum_pattern_mean_dist_to_nucleus_um",
    "zsum_pattern_median_dist_to_nucleus_um",
    "zsum_pattern_apical_fraction_pct",
    "zsum_pattern_com_vs_pattern_offset_um",
    # MaxIP pattern (sensitivity)
    "maxip_pattern_perinuclear_5um_pct",
    "maxip_pattern_mean_dist_to_nucleus_um",
]


def main():
    full = pl.read_csv(COMBINED)
    for proj in ("zsum", "maxip"):
        for mask in ("crop", "pattern"):
            full = add_radial_derived(full, proj, mask)

    for sheet, pairs in SHEET_PAIRS.items():
        family = len(pairs)
        sub = full.filter(pl.col("sheet") == sheet)
        print(f"\n\n{'=' * 78}\n{sheet}  (n={sub.height}, family m={family})\n{'=' * 78}")

        for pair in pairs:
            print(f"\n--- {pair[0]} vs {pair[1]} ---")
            print(f"{'Metric':<52}  {'d':>6}  {'p (Šídák)':>10}  sig")
            print("-" * 80)
            results = []
            for m in METRIC_SET_MITO:
                if m not in sub.columns:
                    continue
                r = test_pair(sub, m, pair, family)
                if r is not None:
                    results.append((m, r["d"], r["p_sidak"]))
            results.sort(key=lambda r: r[2] if np.isfinite(r[2]) else 1.0)
            for m, d, p in results:
                sig = "***" if p < 0.001 else "**" if p < 0.01 else "*" if p < 0.05 else "ns"
                print(f"{m:<52}  {d:+6.2f}  {p:10.4f}  {sig}")


if __name__ == "__main__":
    main()
