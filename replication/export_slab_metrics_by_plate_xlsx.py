"""Export the two slab metrics — centrosomal [18, 33) µm and peripheral
[41, 56) µm — broken out by (condition, plate_date) for each comparison
sheet.

One worksheet per comparison. Each worksheet has two side-by-side
blocks: the peripheral metric on the left and the centrosomal metric on
the right. Within each block, every (condition, plate_date) tuple gets
its own column. Cells are filled with the per-cell metric values down
each column; columns have variable length because each (condition, plate)
has a different cell count.

Plate ordering: plates that appear across the most conditions first
(common plates), then plates unique to a single condition last (these
are the "outlier" replicates contributing only to one isoform). Within
each tier, ordered by plate date.
"""
from __future__ import annotations

import argparse
import collections
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

CENTROSOMAL = (18, 33)
PERIPHERAL  = (41, 56)

SHEET_PREFIX = {
    "TRAK isoform (mito)":       "mito",
    "TRAK isoform (peroxisome)": "perox",
    "TRAK isoform (60mer)":      "60mer",
    "TRAK1 helix muts":          "T1helix",
    "TRAK2 helix muts":          "T2helix",
    "MAPK9 siRNA":               "MAPK9",
}


def plate_date(plate: str) -> str:
    m = re.match(r"(\d{6})", plate)
    return m.group(1) if m else plate


def add_slab_columns(df: pl.DataFrame) -> pl.DataFrame:
    centro_cols = [f"wedge_r_{i:02d}_{i+1:02d}um_pct"
                   for i in range(*CENTROSOMAL)
                   if f"wedge_r_{i:02d}_{i+1:02d}um_pct" in df.columns]
    periph_cols = [f"wedge_r_{i:02d}_{i+1:02d}um_pct"
                   for i in range(*PERIPHERAL)
                   if f"wedge_r_{i:02d}_{i+1:02d}um_pct" in df.columns]
    return df.with_columns([
        pl.sum_horizontal([pl.col(c) for c in centro_cols]).alias("centrosomal_pct"),
        pl.sum_horizontal([pl.col(c) for c in periph_cols]).alias("peripheral_pct"),
    ])


def order_plates_globally(sheet_df: pl.DataFrame) -> list[str]:
    """Return plate dates ordered by (n_conditions desc, plate_date asc),
    so common plates come first and condition-unique plates last."""
    counts: dict[str, set[str]] = collections.defaultdict(set)
    for plate, cond in sheet_df.select(["plate_date", "condition"]).iter_rows():
        counts[plate].add(cond)
    return sorted(counts.keys(), key=lambda p: (-len(counts[p]), p))


def write_block(ws, df_groups, col_offset, header_row, label,
                bold, banner_fmt, val_fmt, max_n):
    """Write one metric block. df_groups = list of (cond, plate, values)."""
    n_cols = len(df_groups)
    ws.merge_range(header_row, col_offset, header_row, col_offset + n_cols - 1,
                   label, banner_fmt)
    for k, (cond, plate, _) in enumerate(df_groups):
        ws.write(header_row + 1, col_offset + k, f"{cond} {plate}", bold)
    for k, (_, _, vals) in enumerate(df_groups):
        for r, v in enumerate(vals):
            if np.isfinite(v):
                ws.write(header_row + 2 + r, col_offset + k, float(v), val_fmt)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out",
                    default="replication/figures_wedge_r_ks/slab_metrics_by_plate.xlsx")
    args = ap.parse_args()

    df = load_template_matching(pathlib.Path(
        "replication/wedge_r_ks_out_all_denoised/by_well"))
    df = join_with_metadata(df, REPO / "config/Comparisons_table_v3.xlsx")
    df = add_slab_columns(df)
    df = df.with_columns(
        pl.col("plate").map_elements(plate_date, return_dtype=pl.String)
        .alias("plate_date")
    )

    out_path = pathlib.Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    wb = xlsxwriter.Workbook(str(out_path))

    bold = wb.add_format({"bold": True})
    banner_fmt = wb.add_format({"bold": True, "bg_color": "#D9E1F2",
                                "align": "center"})
    val_fmt = wb.add_format({"num_format": "0.000"})

    summary_rows = []
    for sheet_label, prefix in SHEET_PREFIX.items():
        cfg = SHEET_CONFIG[sheet_label]
        cond_order = cfg["conditions"]
        sub = df.filter(pl.col("sheet") == sheet_label)
        if sub.height == 0:
            continue

        plate_order = order_plates_globally(sub)

        # Build the (condition, plate) groups in the user's preferred order:
        # group by condition first, plates in the global order skipping
        # those without cells for that condition.
        peripheral_groups: list[tuple[str, str, np.ndarray]] = []
        centrosomal_groups: list[tuple[str, str, np.ndarray]] = []
        for cond in cond_order:
            for plate in plate_order:
                grp = sub.filter((pl.col("condition") == cond) &
                                 (pl.col("plate_date") == plate))
                if grp.height == 0:
                    continue
                peripheral_groups.append(
                    (cond, plate, grp["peripheral_pct"].to_numpy().astype(float)))
                centrosomal_groups.append(
                    (cond, plate, grp["centrosomal_pct"].to_numpy().astype(float)))

        max_n = max((len(g[2]) for g in peripheral_groups), default=0)

        ws = wb.add_worksheet(prefix[:31])
        ws.write(0, 0, sheet_label, bold)

        # Side-by-side blocks: peripheral first, gap, centrosomal.
        n_periph = len(peripheral_groups)
        n_centro = len(centrosomal_groups)
        gap = 1  # single empty column between blocks

        write_block(ws, peripheral_groups, col_offset=0, header_row=2,
                    label=f"peripheral (41–55 µm) — % intensity",
                    bold=bold, banner_fmt=banner_fmt, val_fmt=val_fmt,
                    max_n=max_n)
        write_block(ws, centrosomal_groups, col_offset=n_periph + gap,
                    header_row=2,
                    label=f"centrosomal (18–32 µm) — % intensity",
                    bold=bold, banner_fmt=banner_fmt, val_fmt=val_fmt,
                    max_n=max_n)

        ws.set_column(0, n_periph - 1, 16, val_fmt)
        ws.set_column(n_periph, n_periph + gap - 1, 2)
        ws.set_column(n_periph + gap, n_periph + gap + n_centro - 1, 16, val_fmt)
        ws.freeze_panes(4, 0)

        summary_rows.append({
            "sheet": sheet_label,
            "worksheet": prefix,
            "n_groups_per_block": n_periph,
            "max_cells_per_group": max_n,
            "total_cells": sub.height,
        })
        print(f"  {sheet_label}: {n_periph} (cond, plate) groups · "
              f"max {max_n} cells/group · {sub.height} total cells")

    # Index sheet at the front.
    idx_ws = wb.add_worksheet("index")
    idx_ws.write_row(0, 0, ["sheet", "worksheet", "n_groups_per_block",
                            "max_cells_per_group", "total_cells"], bold)
    for i, r in enumerate(summary_rows):
        idx_ws.write_row(1 + i, 0,
                         [r["sheet"], r["worksheet"], r["n_groups_per_block"],
                          r["max_cells_per_group"], r["total_cells"]])
    idx_ws.set_column(0, 4, 28)
    wb.worksheets_objs.insert(0, wb.worksheets_objs.pop())

    wb.close()
    print(f"\nwrote {out_path}")


if __name__ == "__main__":
    main()
