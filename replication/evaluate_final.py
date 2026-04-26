"""Multi-sheet evaluator for the final pipeline.

Loads `overnight_final_out/combined_raw.csv`, joins plate/well/sheet/condition
metadata from the canonical `overnight_out/combined.csv`, runs nested ANOVA +
Šídák for the panel-specific pair families Mark uses in his Prism analyses,
and writes a wide ranked CSV plus a focused report on the keeper metrics.

Pair families per sheet match Mark's Prism selections.
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

NEW_CSV = REPO / "replication" / "overnight_final_out" / "combined_raw.csv"
OLD_CSV = REPO / "replication" / "overnight_out" / "combined.csv"
PEROX_CSV = REPO / "replication" / "overnight_final_out" / "peroxisome_metadata.csv"
OUT_DIR = REPO / "replication" / "overnight_final_out"


def load_metadata() -> pl.DataFrame:
    """Concatenate the canonical metadata CSV with the peroxisome addendum."""
    base = pl.read_csv(OLD_CSV).select(
        ["path", "plate", "well", "sheet", "condition"])
    parts = [base]
    if PEROX_CSV.exists():
        parts.append(pl.read_csv(PEROX_CSV).select(
            ["path", "plate", "well", "sheet", "condition"]))
    sixtymer = REPO / "replication" / "overnight_final_out" / "sixtymer_metadata.csv"
    if sixtymer.exists():
        parts.append(pl.read_csv(sixtymer).select(
            ["path", "plate", "well", "sheet", "condition"]))
    return pl.concat(parts, how="vertical")


SHEET_PAIRS = {
    "TRAK isoform (mito)": [
        ("no TRAK", "TRAK1"), ("no TRAK", "TRAK2"), ("TRAK1", "TRAK2"),
    ],
    "TRAK isoform (peroxisome)": [
        ("no TRAK", "TRAK1"), ("no TRAK", "TRAK2"), ("TRAK1", "TRAK2"),
    ],
    "TRAK isoform (60mer)": [
        ("no TRAK", "TRAK1"), ("no TRAK", "TRAK2"), ("TRAK1", "TRAK2"),
    ],
    "TRAK1 helix muts": [
        ("T1 wt", "T1 mDRH"), ("T1 mDRH", "T1 mDRH / dSp"),
    ],
    "TRAK2 helix muts": [
        ("TRAK2", "TRAK2 mDRH"), ("TRAK2 mDRH", "TRAK2 mDRH mSpindly"),
    ],
    "MAPK9 siRNA": [
        ("ctrl ctrl", "MAPK9 ctrl"),
        ("ctrl ctrl", "ctrl Ars"),
        ("ctrl ctrl", "MAPK9 Ars"),
        ("MAPK9 ctrl", "MAPK9 Ars"),
        ("ctrl Ars", "MAPK9 Ars"),
    ],
}

KEEPER_SCALARS = [
    # Mark baselines
    "{p}_perinuclear_5um_pct", "{p}_peripheral_5um_pct",
    "{p}_mean_dist_to_nucleus_um",
    # Y-axis projection
    "{p}_y_gini", "{p}_y_entropy", "{p}_y_sd_um", "{p}_y_skew", "{p}_y_mean_um",
    # Wedge-r polar
    "{p}_wedge_r_gini", "{p}_wedge_r_entropy", "{p}_wedge_r_ks_vs_uniform",
    "{p}_wedge_r_sd_um", "{p}_wedge_r_skew", "{p}_wedge_r_mean_um",
    "{p}_wedge_r_q25_um", "{p}_wedge_r_q50_um", "{p}_wedge_r_q75_um",
    "{p}_wedge_r_20_35um_frac_pct", "{p}_wedge_r_35_55um_frac_pct",
    "{p}_wedge_mt_apex_elongation",
    "{p}_wedge_mt_apex_lam_min_um2",
    "{p}_wedge_frac_pct",
]


def collect(df: pl.DataFrame, metric: str) -> list[ConditionData]:
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


def test_pair(df, metric, pair, family):
    conds = collect(df, metric)
    n2i = {c.name: i for i, c in enumerate(conds)}
    if pair[0] not in n2i or pair[1] not in n2i:
        return None
    i, j = n2i[pair[0]], n2i[pair[1]]
    # Try nested-ANOVA-based pairwise test (preferred when ≥2 plates)
    p = np.nan
    try:
        a = nested_oneway_anova(conds)
        r = sidak_pairwise(conds, a, pairs=[(i, j)])[0]
        if np.isfinite(r["p_sidak"]):
            p = 1 - (1 - r["p_sidak"]) ** family
    except Exception:
        pass
    # Fallback for single-plate datasets: Welch t-test on pooled cells +
    # Šídák correction over the requested family size.
    if not np.isfinite(p):
        try:
            from scipy import stats
            a_vals = conds[i].all_cells
            b_vals = conds[j].all_cells
            if len(a_vals) > 1 and len(b_vals) > 1:
                t = stats.ttest_ind(a_vals, b_vals, equal_var=False)
                p_raw = float(t.pvalue)
                if np.isfinite(p_raw):
                    p = 1 - (1 - p_raw) ** family
        except Exception:
            pass
    return {"d": cohens_d(conds[i].all_cells, conds[j].all_cells), "p": p}


def main():
    if not NEW_CSV.exists():
        print(f"Not found: {NEW_CSV}"); return 1
    new = pl.read_csv(NEW_CSV)
    df = new.join(load_metadata(), on="path", how="left").filter(
        pl.col("condition").is_not_null())
    print(f"Loaded {df.height} cells across "
          f"{sorted(df['sheet'].unique().to_list())}")

    OUT_DIR.mkdir(parents=True, exist_ok=True)

    # Build the keeper scalar list (zsum + maxip variants)
    keeper_metrics = []
    for tmpl in KEEPER_SCALARS:
        for proj in ("zsum", "maxip"):
            keeper_metrics.append(tmpl.format(p=proj))
    keeper_metrics = [m for m in keeper_metrics if m in df.columns]

    rows = []
    for sheet, pairs in SHEET_PAIRS.items():
        sheet_df = df.filter(pl.col("sheet") == sheet)
        if sheet_df.height == 0:
            continue
        family = len(pairs)
        for pair in pairs:
            for m in keeper_metrics:
                r = test_pair(sheet_df, m, pair, family)
                if r and np.isfinite(r["p"]):
                    rows.append({"sheet": sheet,
                                 "pair": f"{pair[0]} vs {pair[1]}",
                                 "metric": m, **r})
    summary = pl.from_dicts(rows)
    summary.write_csv(OUT_DIR / "evaluation_summary.csv")
    print(f"Wrote {len(rows)} rows to {OUT_DIR / 'evaluation_summary.csv'}")

    # Pretty per-sheet, per-pair table
    for sheet, pairs in SHEET_PAIRS.items():
        sheet_rows = summary.filter(pl.col("sheet") == sheet)
        if sheet_rows.height == 0:
            continue
        print(f"\n\n{'=' * 90}\n{sheet}  (family m={len(pairs)})\n{'=' * 90}")
        for pair in pairs:
            label = f"{pair[0]} vs {pair[1]}"
            sub = sheet_rows.filter(pl.col("pair") == label).sort("p")
            print(f"\n--- {label} (top 12) ---")
            print(f"{'Metric':<48}  {'d':>6}  {'p':>10}  sig")
            print("-" * 75)
            for r in sub.head(12).iter_rows(named=True):
                p = r["p"]
                sig = "***" if p < 0.001 else "**" if p < 0.01 else "*" if p < 0.05 else "ns"
                print(f"{r['metric']:<48}  {r['d']:+6.2f}  {p:10.5f}  {sig}")


if __name__ == "__main__":
    main()
