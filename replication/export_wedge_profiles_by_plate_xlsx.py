"""Export the wedge-r profile source data with one worksheet per
comparison sheet, where every (condition, plate_date) combination is its
own column. Cell values are per-bin means across cells of that
(condition, plate) group; a parallel `_sem` column is included so the
user can plot mean ± SEM bands.

Layout per worksheet:
  bin_lo_um | bin_hi_um | bin_center_um |
    no TRAK_250612_mean | no TRAK_250612_sem |
    no TRAK_250710_mean | no TRAK_250710_sem | ...

A small metadata block above the table records cell counts per group.
"""
from __future__ import annotations

import argparse
import pathlib
import re
import sys

import numpy as np
import polars as pl
import xlsxwriter

REPO = pathlib.Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / "replication"))
from plot_metrics import (  # noqa: E402
    SHEET_CONFIG, load_template_matching, join_with_metadata,
)

SHEET_PREFIX = {
    "TRAK isoform (mito)":       "mito",
    "TRAK isoform (peroxisome)": "perox",
    "TRAK isoform (60mer)":      "60mer",
    "TRAK1 helix muts":          "T1helix",
    "TRAK2 helix muts":          "T2helix",
    "MAPK9 siRNA":               "MAPK9",
}


def plate_date(plate: str) -> str:
    """Pull the leading YYMMDD out of e.g. '250612_patterned_plate_3'."""
    m = re.match(r"(\d{6})", plate)
    return m.group(1) if m else plate


def safe_ws_name(prefix: str) -> str:
    return prefix[:31]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out",
                    default="replication/figures_wedge_r_ks/wedge_r_profiles_by_plate.xlsx")
    args = ap.parse_args()

    df = load_template_matching(pathlib.Path(
        "replication/wedge_r_ks_out_all_denoised/by_well"))
    df = join_with_metadata(df, REPO / "config/Comparisons_table_v3.xlsx")

    bin_cols = [f"wedge_r_{i:02d}_{i+1:02d}um_pct" for i in range(60)]
    bin_cols = [c for c in bin_cols if c in df.columns]
    n_bins = len(bin_cols)

    df = df.with_columns(
        pl.col("plate").map_elements(plate_date, return_dtype=pl.String)
        .alias("plate_date")
    )

    out_path = pathlib.Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    wb = xlsxwriter.Workbook(str(out_path))

    bold = wb.add_format({"bold": True})
    bold_grey = wb.add_format({"bold": True, "bg_color": "#EEEEEE"})
    grey = wb.add_format({"bg_color": "#EEEEEE"})
    bin_fmt = wb.add_format({"num_format": "0.0"})
    val_fmt = wb.add_format({"num_format": "0.000000"})

    summary_rows = []
    for sheet_label, prefix in SHEET_PREFIX.items():
        cfg = SHEET_CONFIG[sheet_label]
        cond_order = cfg["conditions"]
        sub = df.filter(pl.col("sheet") == sheet_label)
        if sub.height == 0:
            continue

        plates = sorted(sub["plate_date"].unique().to_list())
        ws = wb.add_worksheet(safe_ws_name(prefix))

        # --- Top metadata block: group cell counts ---
        ws.write(0, 0, f"{sheet_label}", bold)
        ws.write(1, 0, "n_cells per (condition, plate_date):", bold_grey)
        ws.write(2, 0, "condition", bold_grey)
        for j, plate in enumerate(plates):
            ws.write(2, 1 + j, plate, bold_grey)
        for i, cond in enumerate(cond_order):
            ws.write(3 + i, 0, cond, grey)
            for j, plate in enumerate(plates):
                n = sub.filter((pl.col("condition") == cond) &
                               (pl.col("plate_date") == plate)).height
                ws.write(3 + i, 1 + j, n, grey)
        header_row = 3 + len(cond_order) + 2  # blank gap before the data table

        # --- Data table ---
        # Column 0: bin_lo, 1: bin_hi, 2: bin_center, then mean/sem per (cond, plate).
        ws.write(header_row, 0, "bin_lo_um", bold)
        ws.write(header_row, 1, "bin_hi_um", bold)
        ws.write(header_row, 2, "bin_center_um", bold)

        col_map: list[tuple[str, str, str]] = []  # (condition, plate, kind)
        col = 3
        for cond in cond_order:
            for plate in plates:
                ws.write(header_row, col, f"{cond}_{plate}_mean", bold)
                ws.write(header_row, col + 1, f"{cond}_{plate}_sem", bold)
                col_map.append((cond, plate, "mean")); col += 1
                col_map.append((cond, plate, "sem"));  col += 1

        # Pre-compute per-(cond, plate) means/sems for every bin.
        n_groups = len([1 for cond in cond_order for _ in plates])
        means = np.full((n_bins, n_groups), np.nan)
        sems  = np.full((n_bins, n_groups), np.nan)
        ns    = np.zeros(n_groups, dtype=int)
        idx = 0
        for cond in cond_order:
            for plate in plates:
                grp = sub.filter((pl.col("condition") == cond) &
                                 (pl.col("plate_date") == plate))
                ns[idx] = grp.height
                if grp.height > 0:
                    arr = grp.select(bin_cols).to_numpy().astype(float)
                    means[:, idx] = np.nanmean(arr, axis=0)
                    if grp.height > 1:
                        sems[:, idx] = (np.nanstd(arr, axis=0, ddof=1) /
                                        np.sqrt(grp.height))
                idx += 1

        for b in range(n_bins):
            row = header_row + 1 + b
            ws.write(row, 0, b, bin_fmt)
            ws.write(row, 1, b + 1, bin_fmt)
            ws.write(row, 2, b + 0.5, bin_fmt)
            for k, (_, _, kind) in enumerate(col_map):
                grp_idx = k // 2
                v = means[b, grp_idx] if kind == "mean" else sems[b, grp_idx]
                if np.isfinite(v):
                    ws.write(row, 3 + k, float(v), val_fmt)
                else:
                    ws.write_blank(row, 3 + k, None, val_fmt)

        ws.set_column(0, 2, 12, bin_fmt)
        ws.set_column(3, 3 + len(col_map) - 1, 18, val_fmt)
        ws.freeze_panes(header_row + 1, 3)

        summary_rows.append({
            "sheet": sheet_label,
            "worksheet": safe_ws_name(prefix),
            "conditions": ", ".join(cond_order),
            "plate_dates": ", ".join(plates),
            "n_groups": len(col_map) // 2,
            "total_cells": sub.height,
        })
        print(f"  {sheet_label}: {len(cond_order)} conds × {len(plates)} plates "
              f"= {len(col_map)//2} groups, {sub.height} cells")

    # Index sheet at the front.
    idx_ws = wb.add_worksheet("index")
    idx_ws.write_row(0, 0, ["sheet", "worksheet", "conditions",
                            "plate_dates", "n_groups", "total_cells"], bold)
    for i, r in enumerate(summary_rows):
        idx_ws.write_row(1 + i, 0,
                         [r["sheet"], r["worksheet"], r["conditions"],
                          r["plate_dates"], r["n_groups"], r["total_cells"]])
    idx_ws.set_column(0, 5, 30)
    wb.worksheets_objs.insert(0, wb.worksheets_objs.pop())  # move to front

    wb.close()
    print(f"\nwrote {out_path}")


if __name__ == "__main__":
    main()
