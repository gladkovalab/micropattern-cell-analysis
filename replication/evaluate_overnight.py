"""Evaluate the overnight_run metrics the same way derived_metrics.py
evaluates CSV-derived metrics: nested one-way ANOVA + Šídák pairwise over
the sheet-specific pair family Mark uses in his Prism files.

Reads replication/overnight_out/combined.csv and the comparisons table,
scores every metric × pair, writes ranked CSVs, prints a headline summary
comparing MaxIP variants vs z-sum variants vs Mark's baseline.
"""
from __future__ import annotations

import pathlib
import sys
import warnings
from itertools import combinations

import fastexcel
import numpy as np
import polars as pl

REPO = pathlib.Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / "replication"))

from replicate_stats import (  # noqa: E402
    ConditionData,
    mixedlm_condition_test,
    nested_oneway_anova,
    sidak_pairwise,
)
from evaluate_metrics import cohens_d  # noqa: E402

warnings.filterwarnings("ignore")

OVERNIGHT_CSV = REPO / "replication" / "overnight_out" / "combined.csv"
OUT_DIR = REPO / "replication" / "overnight_eval_out"


# Pair families mirror Mark's Prism selections (same as derived_metrics.py)
SHEET_PAIRS: dict[str, list[tuple[str, str]]] = {
    "TRAK isoform (mito)": [("no TRAK", "TRAK1"), ("no TRAK", "TRAK2"), ("TRAK1", "TRAK2")],
    "TRAK1 helix muts": [("T1 wt", "T1 mDRH"), ("T1 mDRH", "T1 mDRH / dSp")],
    "TRAK2 helix muts": [("TRAK2", "TRAK2 mDRH"), ("TRAK2 mDRH", "TRAK2 mDRH mSpindly")],
    "MAPK9 siRNA": [
        # Fig 4E peripheral family (A-B, A-C, A-D in Mark's Prism column order)
        ("ctrl ctrl", "MAPK9 ctrl"), ("ctrl ctrl", "ctrl Ars"), ("ctrl ctrl", "MAPK9 Ars"),
        # Fig S11 F perinuclear family (A-C, B-D, C-D)
        ("MAPK9 ctrl", "MAPK9 Ars"), ("ctrl Ars", "MAPK9 Ars"),
    ],
}


def _collect(long: pl.DataFrame, sheet: str, metric: str) -> list[ConditionData]:
    sub = long.filter(
        (pl.col("sheet") == sheet) &
        pl.col(metric).is_not_null() &
        pl.col(metric).is_not_nan()
    )
    conds: list[ConditionData] = []
    for cond_name in sorted(sub["condition"].unique().to_list()):
        g = sub.filter(pl.col("condition") == cond_name)
        if g.height == 0:
            continue
        plate_cells: dict[str, np.ndarray] = {}
        for plate, grp in g.group_by("plate"):
            plate_key = plate[0] if isinstance(plate, tuple) else plate
            plate_cells[plate_key] = grp[metric].to_numpy().astype(float)
        conds.append(ConditionData(name=cond_name, plate_cells=plate_cells))
    return conds


def evaluate_pair(long: pl.DataFrame, sheet: str, metric: str,
                  pair: tuple[str, str], family_size: int) -> dict | None:
    all_conds = _collect(long, sheet, metric)
    if len(all_conds) < 2:
        return None
    name_to_idx = {c.name: i for i, c in enumerate(all_conds)}
    if pair[0] not in name_to_idx or pair[1] not in name_to_idx:
        return None
    i, j = name_to_idx[pair[0]], name_to_idx[pair[1]]
    ci, cj = all_conds[i], all_conds[j]
    if not ci.plate_cells or not cj.plate_cells:
        return None

    a = nested_oneway_anova(all_conds)
    # compute raw pairwise with m=1 (we apply Šídák family correction below)
    clr = sidak_pairwise(all_conds, a, pairs=[(i, j)])[0]
    p_raw = clr["p_sidak"]  # raw, because m=1
    p_sidak = 1 - (1 - p_raw) ** family_size if p_raw is not None and np.isfinite(p_raw) else np.nan
    return {
        "sheet": sheet,
        "metric": metric,
        "pair": f"{pair[0]} vs {pair[1]}",
        "family": family_size,
        "n_i": ci.all_cells.size,
        "n_j": cj.all_cells.size,
        "mean_i": float(ci.all_cells.mean()),
        "mean_j": float(cj.all_cells.mean()),
        "diff": float(ci.all_cells.mean() - cj.all_cells.mean()),
        "cohens_d": cohens_d(ci.all_cells, cj.all_cells),
        "f_classical": a["F"],
        "p_classical_anova": a["p"],
        "t_classical": clr["t"],
        "p_raw": p_raw,
        "p_sidak": p_sidak,
    }


def expand_metrics(df: pl.DataFrame) -> pl.DataFrame:
    """Add derived composite metrics per (projection, mask): diff, ratio,
    peri-share. These match the polarization story from the CSV-derived
    memo but computed from every (projection, mask) combination."""
    exprs = []
    for proj in ("zsum", "maxip"):
        for mask in ("crop", "pattern"):
            peri = f"{proj}_{mask}_peripheral_5um_pct"
            nuc = f"{proj}_{mask}_perinuclear_5um_pct"
            if peri in df.columns and nuc in df.columns:
                eps = 1e-9
                exprs.append((pl.col(peri) - pl.col(nuc)).alias(f"{proj}_{mask}_peri_minus_nuc"))
                exprs.append((pl.col(peri) / (pl.col(nuc) + eps)).alias(f"{proj}_{mask}_peri_over_nuc"))
                exprs.append((100 * pl.col(peri) / (pl.col(peri) + pl.col(nuc) + eps)).alias(f"{proj}_{mask}_peripheral_share_pct"))
    if exprs:
        df = df.with_columns(exprs)
    return df


def main():
    if not OVERNIGHT_CSV.exists():
        print(f"Not found: {OVERNIGHT_CSV}")
        return 1
    df = pl.read_csv(OVERNIGHT_CSV)
    df = expand_metrics(df)
    meta_cols = {"path", "template_matching_score", "lateral_pixel_pitch_um",
                 "sheet", "condition", "plate", "well",
                 "zsum_bg_threshold", "maxip_bg_threshold"}
    metric_cols = [c for c in df.columns if c not in meta_cols]
    print(f"Loaded {df.height} cells × {len(metric_cols)} metrics from {OVERNIGHT_CSV.relative_to(REPO)}")

    rows = []
    for sheet, pairs in SHEET_PAIRS.items():
        m = len(pairs)
        for pair in pairs:
            for metric in metric_cols:
                try:
                    r = evaluate_pair(df, sheet, metric, pair, family_size=m)
                    if r is not None:
                        rows.append(r)
                except Exception as e:
                    pass
    summary = pl.from_dicts(rows)
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    summary.write_csv(OUT_DIR / "summary.csv")
    print(f"Wrote {len(rows)} rows to {OUT_DIR / 'summary.csv'}")

    # Headline report per panel: compare Mark, proposed diff, ratio, plus best MaxIP variant
    MARK_S11 = "zsum_crop_perinuclear_5um_pct"   # Mark Fig S11 baseline (z-sum)
    MARK_4 = "zsum_crop_peripheral_5um_pct"      # Mark Fig 4 baseline (z-sum; denoised not re-run here)
    CORE = [
        ("Mark S11 (zsum perinuclear 5µm)",      MARK_S11),
        ("Mark 4   (zsum peripheral 5µm)",       MARK_4),
        ("Diff (zsum crop)",                     "zsum_crop_peri_minus_nuc"),
        ("Ratio (zsum crop)",                    "zsum_crop_peri_over_nuc"),
        ("Diff (maxip crop)",                    "maxip_crop_peri_minus_nuc"),
        ("Ratio (maxip crop)",                   "maxip_crop_peri_over_nuc"),
        ("MaxIP peripheral 5µm alone",           "maxip_crop_peripheral_5um_pct"),
        ("MaxIP perinuclear 5µm alone",          "maxip_crop_perinuclear_5um_pct"),
    ]
    print()
    for (sheet, pairs) in SHEET_PAIRS.items():
        for pair in pairs:
            pair_label = f"{pair[0]} vs {pair[1]}"
            print(f"\n--- {sheet} · {pair_label} (m={len(pairs)}) ---")
            for nice, metric in CORE:
                r = summary.filter((pl.col("sheet") == sheet) &
                                   (pl.col("pair") == pair_label) &
                                   (pl.col("metric") == metric))
                if r.height == 0:
                    continue
                d = r["cohens_d"].item()
                p = r["p_sidak"].item()
                sig = "***" if p < 0.001 else "**" if p < 0.01 else "*" if p < 0.05 else "ns"
                print(f"  {nice:<38}  d={d:+6.2f}  p={p:.4f}  {sig}")

    # Also rank the best single metric per (sheet, pair) across EVERYTHING
    print("\n=== Best-single-metric per (sheet, pair) across ALL computed metrics ===")
    for (sheet, pairs) in SHEET_PAIRS.items():
        for pair in pairs:
            pair_label = f"{pair[0]} vs {pair[1]}"
            sub = summary.filter((pl.col("sheet") == sheet) &
                                 (pl.col("pair") == pair_label) &
                                 pl.col("p_sidak").is_not_null() &
                                 pl.col("p_sidak").is_not_nan())
            if sub.height == 0:
                continue
            best = sub.sort("p_sidak").head(5)
            print(f"\n  {sheet} · {pair_label}:")
            for r in best.iter_rows(named=True):
                sig = "***" if r["p_sidak"] < 0.001 else "**" if r["p_sidak"] < 0.01 else "*" if r["p_sidak"] < 0.05 else "ns"
                print(f"    d={r['cohens_d']:+6.2f}  p={r['p_sidak']:.4f}  {sig:<3}  {r['metric']}")


if __name__ == "__main__":
    sys.exit(main())
