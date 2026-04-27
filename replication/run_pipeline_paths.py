"""Run the patched template_matching_bulk pipeline on a curated set of cell
paths (e.g. only `denoised/` ND2s for a specific sheet).

Mark's `template_matching_bulk.main` walks a directory tree, which mixes
raw + denoised cells and processes the entire tree below the given root.
For testing the wedge-r KS patch we want to control which cells get
processed — this thin driver does that:

  1. Reads `replication/overnight_out/combined.csv` for sheet/condition/
     plate/well metadata.
  2. Selects target wells matching `--sheet`.
  3. Walks each well's `denoised/` subdir (or `--variant raw` for the
     non-denoised ND2s) and runs `score_template_match` on every cell.
  4. Writes one `template_matching.csv` per well, mirroring Mark's output
     layout under the directory passed as `--out-root`.

Resumable: cells whose JSON checkpoint already exists under
`{out_root}/by_well/.../cells/` are skipped. SMB drops only lose the
in-progress cell.

Example:
    pixi run python replication/run_pipeline_paths.py \\
        --sheet "TRAK isoform (mito)" --variant denoised \\
        --out-root replication/wedge_r_ks_out
"""
from __future__ import annotations

import argparse
import json
import os
import pathlib
import sys
import time
import traceback

import fastexcel
import numpy as np
import polars as pl

REPO = pathlib.Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))


def load_comparisons_table(xlsx: pathlib.Path) -> pl.DataFrame:
    """Long-form (plate, well, sheet, condition) view of Mark's
    `config/Comparisons_table_v3.xlsx`. Each xlsx sheet maps one
    comparison family — rows are plates, columns are conditions, cell
    values are well codes (e.g. "B02"). This unpivots into per-well rows."""
    fe = fastexcel.read_excel(str(xlsx))
    rows = []
    for sheet_name in fe.sheet_names:
        df = fe.load_sheet_by_name(sheet_name).to_polars()
        plate_col = df.columns[0]
        for record in df.iter_rows(named=True):
            plate = record[plate_col]
            if not plate:
                continue
            for cond in df.columns[1:]:
                well = record[cond]
                if well:
                    rows.append({"plate": plate, "well": well,
                                 "sheet": sheet_name, "condition": cond})
    return pl.from_dicts(rows)


def discover_cells(comparisons_xlsx: pathlib.Path, sheet: str, variant: str,
                   data_root: pathlib.Path) -> list[pathlib.Path]:
    """Return the list of ND2 paths to process. `variant` is 'denoised' or
    'raw'."""
    df = load_comparisons_table(comparisons_xlsx).filter(pl.col("sheet") == sheet)
    target_wells = sorted({(plate, well) for plate, well
                           in df.select(["plate", "well"]).unique().iter_rows()})
    cells: list[pathlib.Path] = []
    for plate, well_short in target_wells:
        plate_dir = data_root / plate
        if not plate_dir.exists():
            print(f"  MISSING plate dir: {plate_dir}", flush=True)
            continue
        # well metadata stores e.g. "B06"; on disk this is e.g. "B06_250528_..."
        well_prefix = f"{well_short}_"
        for well_dir in plate_dir.iterdir():
            if not (well_dir.is_dir() and well_dir.name.startswith(well_prefix)):
                continue
            scan_dir = well_dir / "denoised" if variant == "denoised" else well_dir
            if not scan_dir.exists():
                continue
            for fn in sorted(scan_dir.iterdir()):
                if not fn.name.endswith(".nd2"):
                    continue
                if not fn.name.lower().startswith("cell"):
                    continue
                if variant == "raw" and fn.parent.name == "denoised":
                    continue
                cells.append(fn)
    return cells


def run(cells: list[pathlib.Path], out_root: pathlib.Path,
        data_root: pathlib.Path):
    import template_matching_bulk as tmb  # imports patched pipeline

    template_hat = tmb.get_template_hat(1326)
    template = tmb.get_padded_template_at_width(1326)

    records_by_well: dict[pathlib.Path, list[dict]] = {}
    for i, img_path in enumerate(cells, 1):
        well_dir = img_path.parent if img_path.parent.name != "denoised" else img_path.parent.parent
        rel = well_dir.resolve().relative_to(data_root.resolve())
        chk_dir = out_root / "by_well" / rel / "cells"
        chk_dir.mkdir(parents=True, exist_ok=True)
        chk_path = chk_dir / f"{img_path.stem}.json"

        if chk_path.exists():
            try:
                m = json.loads(chk_path.read_text())
                records_by_well.setdefault(well_dir, []).append(m)
                print(f"  [{i}/{len(cells)}] CACHED {img_path.relative_to(data_root)}",
                      flush=True)
                continue
            except Exception:
                pass

        try:
            t0 = time.time()
            out = tmb.score_template_match(img_path, template_hat=template_hat,
                                           template=template)
            out["path"] = str(img_path)
            out["template_matching_score"] = out.pop("score")
            # Coerce xarray scalars -> floats for CSV / JSON safety.
            for k, v in list(out.items()):
                if hasattr(v, "values") and not isinstance(v, (str, bytes)):
                    try:
                        out[k] = float(v.values.item())
                    except Exception:
                        pass
                elif isinstance(v, np.generic):
                    out[k] = float(v) if v.dtype.kind in "fc" else int(v)
            chk_path.write_text(json.dumps(out))
            records_by_well.setdefault(well_dir, []).append(out)
            print(f"  [{i}/{len(cells)}] OK ({time.time()-t0:.1f}s) "
                  f"{img_path.relative_to(data_root)}", flush=True)
        except Exception as e:
            print(f"  [{i}/{len(cells)}] ERR {img_path}: {e}", flush=True)
            traceback.print_exc()

    # Write per-well CSV mirroring Mark's layout, plus aggregated CSV.
    for well_dir, recs in records_by_well.items():
        rel = well_dir.resolve().relative_to(data_root.resolve())
        # Add Mark's percent_total derivations (his main() does this in-CSV).
        df = pl.from_dicts(recs)
        for d in [1, 2, 3, 4, 5]:
            cols = []
            if f"peripheral_{d}um_sum" in df.columns and "mitochondria_sum" in df.columns:
                cols.append((pl.col(f"peripheral_{d}um_sum") /
                             pl.col("mitochondria_sum") * 100)
                            .alias(f"peripheral_{d}um_percent_total"))
            if f"perinuclear_{d}um_sum" in df.columns and "mitochondria_sum" in df.columns:
                cols.append((pl.col(f"perinuclear_{d}um_sum") /
                             pl.col("mitochondria_sum") * 100)
                            .alias(f"perinuclear_{d}um_percent_total"))
            if f"peripheral_{d}um_simple_sum" in df.columns and "mitochondria_sum" in df.columns:
                cols.append((pl.col(f"peripheral_{d}um_simple_sum") /
                             pl.col("mitochondria_sum") * 100)
                            .alias(f"peripheral_{d}um_simple_percent_total"))
            if cols:
                df = df.with_columns(cols)
        out_csv = out_root / "by_well" / rel / "template_matching.csv"
        out_csv.parent.mkdir(parents=True, exist_ok=True)
        df.write_csv(out_csv)
        print(f"[run_pipeline_paths] wrote {out_csv} ({df.height} cells)", flush=True)

    all_recs = [r for recs in records_by_well.values() for r in recs]
    if all_recs:
        df_all = pl.from_dicts(all_recs)
        for d in [1, 2, 3, 4, 5]:
            cols = []
            if f"peripheral_{d}um_sum" in df_all.columns and "mitochondria_sum" in df_all.columns:
                cols.append((pl.col(f"peripheral_{d}um_sum") /
                             pl.col("mitochondria_sum") * 100)
                            .alias(f"peripheral_{d}um_percent_total"))
            if f"perinuclear_{d}um_sum" in df_all.columns and "mitochondria_sum" in df_all.columns:
                cols.append((pl.col(f"perinuclear_{d}um_sum") /
                             pl.col("mitochondria_sum") * 100)
                            .alias(f"perinuclear_{d}um_percent_total"))
            if cols:
                df_all = df_all.with_columns(cols)
        agg_csv = out_root / "combined.csv"
        df_all.write_csv(agg_csv)
        print(f"[run_pipeline_paths] wrote {agg_csv} ({df_all.height} cells total)",
              flush=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sheet", required=True,
                    help='Sheet to process, e.g. "TRAK isoform (mito)".')
    ap.add_argument("--variant", choices=["raw", "denoised"], default="denoised",
                    help="Which ND2 set to process per well.")
    ap.add_argument("--comparisons-xlsx",
                    default=str(REPO / "config" / "Comparisons_table_v3.xlsx"),
                    help="Authoritative sheet/condition/plate/well map.")
    ap.add_argument("--data-root",
                    default=str(REPO / "mark_data" / "patterned_data"))
    ap.add_argument("--out-root",
                    default=str(REPO / "replication" / "wedge_r_ks_out"))
    args = ap.parse_args()

    data_root = pathlib.Path(args.data_root).resolve()
    os.environ["MICROPATTERN_DATA_ROOT"] = str(data_root)
    out_root = pathlib.Path(args.out_root).resolve()
    out_root.mkdir(parents=True, exist_ok=True)

    cells = discover_cells(pathlib.Path(args.comparisons_xlsx).resolve(),
                           args.sheet, args.variant, data_root)
    print(f"[run_pipeline_paths] sheet={args.sheet!r} variant={args.variant} "
          f"-> {len(cells)} cells", flush=True)
    if not cells:
        sys.exit("no cells matched; exiting")

    run(cells, out_root, data_root)


if __name__ == "__main__":
    main()
