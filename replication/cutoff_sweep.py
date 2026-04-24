"""Cutoff-combination sweep for the peripheral÷perinuclear ratio and
peripheral−perinuclear difference.

Mark's pipeline outputs cumulative signal within 1, 2, 3, 4, 5 µm of each zone
boundary. So for each cell we can form:
    diff(X, Y)  = peripheral_Xum − perinuclear_Yum
    ratio(X, Y) = peripheral_Xum / perinuclear_Yum
with (X, Y) both ranging over {1..5} — 25 combinations. Mark's existing metric
corresponds to (X=5, Y=5). We score every (X, Y) with the same nested ANOVA +
Šídák framework as `derived_metrics.py`, rank by Šídák p on the key per-panel
claim, and print the top hits.

This runs on Mark's pre-existing z-sum CSVs only (no ND2s needed). The MaxIP
variant would require re-running the pipeline with each cutoff — not cheap.
The z-sum sweep is a free win.
"""
from __future__ import annotations

import pathlib
import sys
import warnings
from itertools import combinations, product

import fastexcel
import numpy as np
import polars as pl

REPO = pathlib.Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / "replication"))

from derived_metrics import _find_well_dir, _read_csv_values, DISTANCES  # noqa: E402
from replicate_stats import ConditionData, nested_oneway_anova, sidak_pairwise  # noqa: E402
from evaluate_metrics import cohens_d  # noqa: E402

warnings.filterwarnings("ignore")

OUT_DIR = REPO / "replication" / "cutoff_sweep_out"
OUT_DIR.mkdir(parents=True, exist_ok=True)

COMPARISONS_XLSX = REPO / "config" / "Comparisons_table_v3.xlsx"

# Same sheet/pair families as the main pitch.
SHEET_PAIRS: dict[str, list[tuple[str, str]]] = {
    "TRAK isoform (mito)": [("no TRAK", "TRAK1"), ("no TRAK", "TRAK2"), ("TRAK1", "TRAK2")],
    "TRAK1 helix muts": [("T1 wt", "T1 mDRH"), ("T1 mDRH", "T1 mDRH / dSp")],
    "TRAK2 helix muts": [("TRAK2", "TRAK2 mDRH"), ("TRAK2 mDRH", "TRAK2 mDRH mSpindly")],
    "MAPK9 siRNA": [
        ("ctrl ctrl", "MAPK9 ctrl"), ("ctrl ctrl", "ctrl Ars"), ("ctrl ctrl", "MAPK9 Ars"),
        ("MAPK9 ctrl", "MAPK9 Ars"), ("ctrl Ars", "MAPK9 Ars"),
    ],
}


def gather_cells() -> pl.DataFrame:
    """One row per (sheet, condition, plate, well, cell_idx) with cumulative
    peripheral_Xum and perinuclear_Yum percent-of-total for X, Y in 1..5."""
    reader = fastexcel.read_excel(COMPARISONS_XLSX)
    rows = []
    for sheet_name in reader.sheet_names:
        if sheet_name not in SHEET_PAIRS:
            continue
        df = pl.read_excel(COMPARISONS_XLSX, sheet_name=sheet_name)
        plate_col = df.columns[0]
        for cond in df.columns[1:]:
            for rec in df.iter_rows(named=True):
                plate = rec[plate_col]
                well = rec[cond]
                if not well:
                    continue
                well_dir = _find_well_dir(plate, well)
                if well_dir is None:
                    continue
                vals = _read_csv_values(well_dir / "template_matching.csv")
                if vals is None:
                    continue
                n_cells = next(iter(vals.values())).size
                for i in range(n_cells):
                    row = {"sheet": sheet_name, "condition": cond, "plate": plate,
                           "well": well, "cell_idx": i}
                    for d in DISTANCES:
                        row[f"peripheral_{d}um"] = float(vals[f"peripheral_{d}um"][i]) if f"peripheral_{d}um" in vals else None
                        row[f"perinuclear_{d}um"] = float(vals[f"perinuclear_{d}um"][i]) if f"perinuclear_{d}um" in vals else None
                    rows.append(row)
    return pl.from_dicts(rows)


def add_cutoff_combos(df: pl.DataFrame) -> pl.DataFrame:
    eps = 1e-9
    exprs = []
    for x, y in product(DISTANCES, DISTANCES):
        peri = f"peripheral_{x}um"
        nuc = f"perinuclear_{y}um"
        if peri in df.columns and nuc in df.columns:
            exprs.append((pl.col(peri) - pl.col(nuc)).alias(f"diff_p{x}_n{y}"))
            exprs.append((pl.col(peri) / (pl.col(nuc) + eps)).alias(f"ratio_p{x}_n{y}"))
    return df.with_columns(exprs) if exprs else df


def collect_conds(long: pl.DataFrame, sheet: str, metric: str) -> list[ConditionData]:
    sub = long.filter((pl.col("sheet") == sheet) & pl.col(metric).is_not_null() &
                      pl.col(metric).is_not_nan())
    conds = []
    for cond_name in sorted(sub["condition"].unique().to_list()):
        g = sub.filter(pl.col("condition") == cond_name)
        if g.height == 0:
            continue
        plate_cells = {}
        for plate, grp in g.group_by("plate"):
            plate_key = plate[0] if isinstance(plate, tuple) else plate
            plate_cells[plate_key] = grp[metric].to_numpy().astype(float)
        conds.append(ConditionData(name=cond_name, plate_cells=plate_cells))
    return conds


def evaluate(long: pl.DataFrame, sheet: str, metric: str,
             pair: tuple[str, str], m: int) -> dict | None:
    conds = collect_conds(long, sheet, metric)
    if len(conds) < 2:
        return None
    name_to_idx = {c.name: i for i, c in enumerate(conds)}
    if pair[0] not in name_to_idx or pair[1] not in name_to_idx:
        return None
    i, j = name_to_idx[pair[0]], name_to_idx[pair[1]]
    ci, cj = conds[i], conds[j]
    if not ci.plate_cells or not cj.plate_cells:
        return None
    a = nested_oneway_anova(conds)
    clr = sidak_pairwise(conds, a, pairs=[(i, j)])[0]
    p_raw = clr["p_sidak"]  # raw because m=1
    if p_raw is None or not np.isfinite(p_raw):
        return None
    p_sidak = 1 - (1 - p_raw) ** m
    return {
        "sheet": sheet, "pair": f"{pair[0]} vs {pair[1]}", "metric": metric,
        "cohens_d": cohens_d(ci.all_cells, cj.all_cells),
        "mean_i": float(ci.all_cells.mean()), "mean_j": float(cj.all_cells.mean()),
        "n_i": ci.all_cells.size, "n_j": cj.all_cells.size,
        "p_raw": p_raw, "p_sidak": p_sidak,
    }


def main():
    print("Gathering per-cell cumulative values from Mark's CSVs…")
    long = gather_cells()
    print(f"  {long.height} cells across {long['sheet'].n_unique()} sheets")
    long = add_cutoff_combos(long)
    metrics = [c for c in long.columns if c.startswith(("diff_p", "ratio_p"))]
    print(f"  {len(metrics)} cutoff-combination metrics added")

    results = []
    for sheet, pairs in SHEET_PAIRS.items():
        m = len(pairs)
        for pair in pairs:
            for metric in metrics:
                r = evaluate(long, sheet, metric, pair, m=m)
                if r is not None:
                    results.append(r)
    df = pl.from_dicts(results)
    df.write_csv(OUT_DIR / "cutoff_sweep.csv")
    print(f"Wrote {OUT_DIR / 'cutoff_sweep.csv'} ({df.height} rows)")

    # Headline per (sheet, pair): baseline (5,5) + best (X,Y) for diff and ratio
    print("\n=== Best (X µm peripheral, Y µm perinuclear) per (sheet, pair) ===")
    pl.Config.set_tbl_rows(80); pl.Config.set_tbl_width_chars(200); pl.Config.set_fmt_str_lengths(40)
    for (sheet, pair_str), g in df.group_by(["sheet", "pair"]):
        g = g.sort("p_sidak")
        baseline_diff = g.filter(pl.col("metric") == "diff_p5_n5")
        baseline_ratio = g.filter(pl.col("metric") == "ratio_p5_n5")
        best_diff = g.filter(pl.col("metric").str.starts_with("diff_")).head(1)
        best_ratio = g.filter(pl.col("metric").str.starts_with("ratio_")).head(1)
        print(f"\n{sheet}  ·  {pair_str}")
        for label, sub in (("5,5 diff ", baseline_diff), ("best diff", best_diff),
                           ("5,5 ratio", baseline_ratio), ("best ratio", best_ratio)):
            if sub.height == 0:
                continue
            r = sub.row(0, named=True)
            sig = "***" if r["p_sidak"] < 0.001 else "**" if r["p_sidak"] < 0.01 else "*" if r["p_sidak"] < 0.05 else "ns"
            print(f"  {label}: {r['metric']:<20}  d={r['cohens_d']:+6.2f}  p={r['p_sidak']:.4f}  {sig}")


if __name__ == "__main__":
    sys.exit(main())
