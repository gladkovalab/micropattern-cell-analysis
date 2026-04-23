"""Rebuild replication/overnight_out/combined.csv from the per-well CSVs,
dropping the sheet/condition/plate/well tags (which were wrong when a well
is referenced by multiple sheets) and replacing them by joining against the
comparisons table at load time.

Every ND2 file has metrics that are a property of the image alone; the
(sheet, condition) label is a property of the comparisons mapping. We
generate one row per (path, sheet, condition) combination so a well shared
across two sheets contributes to both.
"""
from __future__ import annotations

import pathlib
import re
import sys

import fastexcel
import polars as pl

REPO = pathlib.Path(__file__).resolve().parent.parent
BY_WELL = REPO / "replication" / "overnight_out" / "by_well"
OUT_CSV = REPO / "replication" / "overnight_out" / "combined.csv"
COMPARISONS_XLSX = REPO / "config" / "Comparisons_table_v3.xlsx"

TARGET_SHEETS = {
    "TRAK isoform (mito)", "TRAK1 helix muts", "TRAK2 helix muts", "MAPK9 siRNA",
}


def build_comparison_index() -> list[dict]:
    """Flatten the comparisons workbook into one row per (sheet, condition, plate, well)."""
    reader = fastexcel.read_excel(COMPARISONS_XLSX)
    out = []
    for sheet_name in reader.sheet_names:
        if sheet_name not in TARGET_SHEETS:
            continue
        df = pl.read_excel(COMPARISONS_XLSX, sheet_name=sheet_name)
        plate_col = df.columns[0]
        for cond in df.columns[1:]:
            for row in df.iter_rows(named=True):
                plate = row[plate_col]
                well = row[cond]
                if not plate or not well:
                    continue
                out.append({"sheet": sheet_name, "condition": cond,
                            "plate": plate, "well": well})
    return out


def plate_and_well_from_path(path_str: str) -> tuple[str, str] | None:
    """Extract (plate, well_prefix) from an ND2 path like
    .../patterned_data/<plate>/<well_dir>/cellN.nd2. well_prefix is the
    short well ID (e.g. 'B06') from the first component of the well dir name."""
    parts = pathlib.Path(path_str).parts
    # walk back until we find /patterned_data/
    for i, p in enumerate(parts):
        if p == "patterned_data":
            if i + 2 < len(parts):
                plate = parts[i + 1]
                well_dir = parts[i + 2]
                m = re.match(r"([A-Za-z]\d+)", well_dir)
                if m:
                    return plate, m.group(1)
    return None


def main():
    # Load all per-well CSVs into one frame (stripping stale tags)
    frames = []
    for csv in sorted(BY_WELL.rglob("metrics.csv")):
        df = pl.read_csv(csv)
        drop = [c for c in ("sheet", "condition", "plate", "well") if c in df.columns]
        if drop:
            df = df.drop(drop)
        frames.append(df)
    if not frames:
        print("no per-well CSVs found")
        return 1
    cells = pl.concat(frames, how="diagonal_relaxed")
    print(f"Loaded {cells.height} cells from {len(frames)} well CSVs")

    # Add (plate, well) derived from path
    plates, wells = [], []
    for p in cells["path"].to_list():
        pw = plate_and_well_from_path(p)
        plates.append(pw[0] if pw else None)
        wells.append(pw[1] if pw else None)
    cells = cells.with_columns([
        pl.Series("plate", plates),
        pl.Series("well", wells),
    ])

    # Cross-join with comparisons index so every (cell, sheet/condition) assignment
    # becomes its own row. A cell referenced by 2 sheets now contributes 2 rows.
    index = pl.from_dicts(build_comparison_index())
    print(f"Comparison index has {index.height} (sheet, cond, plate, well) rows")
    combined = cells.join(index, on=["plate", "well"], how="inner")
    print(f"After join: {combined.height} rows")

    combined.write_csv(OUT_CSV)
    print(f"Wrote {OUT_CSV}")

    # Summary
    print("\nCells per (sheet, condition, plate):")
    print(combined.group_by(["sheet", "condition"]).agg(
        pl.col("plate").n_unique().alias("n_plates"),
        pl.col("path").len().alias("n_rows")
    ).sort(["sheet", "condition"]))


if __name__ == "__main__":
    sys.exit(main())
