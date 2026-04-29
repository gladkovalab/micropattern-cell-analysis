"""Export the wedge-r profile source data as an XLSX with one worksheet
per (sheet, condition).

Reads `replication/figures_wedge_r_ks/wedge_r_profiles_source.csv`
(produced by export_wedge_profiles_csv.py) and pivots it so each comparison
sheet/condition lives in its own tab. Replaces the SEM column with
explicit upper/lower band columns (mean ± SEM) for both 488 mito and the
405 nuclear / perinuclear halo curves.
"""
from __future__ import annotations

import argparse
import pathlib
import re

import polars as pl
import xlsxwriter

# Short worksheet name per (sheet, condition). Excel limits sheet names
# to 31 chars and forbids /\?*[]. Keys must be exactly the labels used in
# SHEET_CONFIG.
SHEET_PREFIX = {
    "TRAK isoform (mito)":       "mito",
    "TRAK isoform (peroxisome)": "perox",
    "TRAK isoform (60mer)":      "60mer",
    "TRAK1 helix muts":          "T1helix",
    "TRAK2 helix muts":          "T2helix",
    "MAPK9 siRNA":               "MAPK9",
}


def safe_sheet_name(prefix: str, condition: str) -> str:
    cond = re.sub(r"[\\/?*\[\]]", "_", condition).strip()
    cond = re.sub(r"\s+", "_", cond)
    name = f"{prefix}_{cond}"
    return name[:31]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--in-csv",
                    default="replication/figures_wedge_r_ks/wedge_r_profiles_source.csv")
    ap.add_argument("--out",
                    default="replication/figures_wedge_r_ks/wedge_r_profiles_source.xlsx")
    args = ap.parse_args()

    df = pl.read_csv(args.in_csv)

    # Add the upper/lower band columns the user wants in place of *_sem.
    df = df.with_columns([
        (pl.col("mito_mean_pct") + pl.col("mito_sem_pct")).alias("mito_upper_pct"),
        (pl.col("mito_mean_pct") - pl.col("mito_sem_pct")).alias("mito_lower_pct"),
        (pl.col("nuc_mask_mean_pct") + pl.col("nuc_mask_sem_pct")).alias("nuc_mask_upper_pct"),
        (pl.col("nuc_mask_mean_pct") - pl.col("nuc_mask_sem_pct")).alias("nuc_mask_lower_pct"),
        (pl.col("perinuc_halo_mean_pct") + pl.col("perinuc_halo_sem_pct")).alias("perinuc_halo_upper_pct"),
        (pl.col("perinuc_halo_mean_pct") - pl.col("perinuc_halo_sem_pct")).alias("perinuc_halo_lower_pct"),
    ])

    # Column order each worksheet ends up with.
    worksheet_cols = [
        "bin_lo_um", "bin_hi_um", "bin_center_um",
        "n_cells_mito", "n_cells_nuc",
        "mito_mean_pct", "mito_upper_pct", "mito_lower_pct",
        "nuc_mask_mean_pct", "nuc_mask_upper_pct", "nuc_mask_lower_pct",
        "perinuc_halo_mean_pct", "perinuc_halo_upper_pct", "perinuc_halo_lower_pct",
    ]

    out_path = pathlib.Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    wb = xlsxwriter.Workbook(str(out_path))
    written = []
    seen_names: set[str] = set()
    for sheet_label in SHEET_PREFIX:
        prefix = SHEET_PREFIX[sheet_label]
        sub_sheet = df.filter(pl.col("sheet") == sheet_label)
        if sub_sheet.height == 0:
            continue
        for cond in sub_sheet["condition"].unique(maintain_order=True).to_list():
            cdf = sub_sheet.filter(pl.col("condition") == cond).select(
                ["sheet", "condition"] + worksheet_cols
            ).sort("bin_lo_um")
            ws_name = safe_sheet_name(prefix, cond)
            base = ws_name
            i = 1
            while ws_name in seen_names:
                suf = f"_{i}"
                ws_name = base[: 31 - len(suf)] + suf
                i += 1
            seen_names.add(ws_name)
            cdf.write_excel(
                workbook=wb, worksheet=ws_name,
                table_style="Table Style Medium 9",
                autofit=True,
            )
            written.append((sheet_label, cond, ws_name, cdf.height))

    # Index sheet first (placed at the front by writing it last and
    # reordering), so reviewers see a guide.
    idx_df = pl.DataFrame({
        "sheet": [w[0] for w in written],
        "condition": [w[1] for w in written],
        "worksheet": [w[2] for w in written],
        "n_bins": [w[3] for w in written],
    })
    idx_df.write_excel(workbook=wb, worksheet="index",
                       table_style="Table Style Medium 4", autofit=True)
    wb.worksheets_objs.insert(0, wb.worksheets_objs.pop())  # move index to front

    wb.close()
    print(f"wrote {out_path}  ({len(written)} condition worksheets + 1 index)")
    for s, c, w, n in written:
        print(f"  [{w:31s}] {s} :: {c}  ({n} bins)")


if __name__ == "__main__":
    main()
