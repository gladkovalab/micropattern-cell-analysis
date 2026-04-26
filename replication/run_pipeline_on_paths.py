"""Run final_pipeline.process_cell on an arbitrary list of cell paths.

The standard pipeline driver (`final_pipeline.run`) walks well directories
and skips `denoised/` subdirectories, so it can't process denoised ND2s
even when they're listed explicitly in metadata. This wrapper bypasses
that walk and processes every path in --paths-csv directly.

Same per-cell JSON checkpointing, per-well metrics.csv, and combined_raw.csv
as the canonical run.
"""
from __future__ import annotations
import argparse, json, pathlib, sys, traceback
import polars as pl

REPO = pathlib.Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / "replication"))

from final_pipeline import process_cell  # noqa: E402
import template_matching_bulk as tmb  # noqa: E402


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--paths-csv", required=True,
                    help="CSV with at least a 'path' column listing ND2 files")
    ap.add_argument("--out-root", required=True)
    args = ap.parse_args()

    out_root = pathlib.Path(args.out_root)
    out_root.mkdir(parents=True, exist_ok=True)
    df = pl.read_csv(args.paths_csv)
    paths = [pathlib.Path(p) for p in df["path"].to_list()]
    print(f"[run_pipeline_on_paths] {len(paths)} cells to process", flush=True)

    template_hat = tmb.get_template_hat(1326)
    template = tmb.get_padded_template_at_width(1326)
    data_root = pathlib.Path(tmb.DATA_ROOT).resolve()

    records_by_well: dict[pathlib.Path, list[dict]] = {}
    for i, img_path in enumerate(paths, 1):
        well_dir = img_path.parent              # may be the `denoised/` subdir
        try:
            rel = well_dir.resolve().relative_to(data_root)
        except ValueError:
            rel = pathlib.Path("misc") / well_dir.name
        cell_chk = out_root / "by_well" / rel / "cells" / f"{img_path.stem}.json"
        cell_chk.parent.mkdir(parents=True, exist_ok=True)
        if cell_chk.exists():
            try:
                m = json.loads(cell_chk.read_text())
                records_by_well.setdefault(well_dir, []).append(m)
                print(f"  [{i}/{len(paths)}] CACHED {img_path}", flush=True)
                continue
            except Exception:
                pass
        try:
            m = process_cell(img_path, template_hat=template_hat,
                             template=template)
            records_by_well.setdefault(well_dir, []).append(m)
            cell_chk.write_text(json.dumps(m))
            print(f"  [{i}/{len(paths)}] OK {img_path}", flush=True)
        except Exception as e:
            print(f"  [{i}/{len(paths)}] ERR {img_path}: {e}", flush=True)
            traceback.print_exc()

    for well_dir, recs in records_by_well.items():
        try:
            rel = well_dir.resolve().relative_to(data_root)
        except ValueError:
            rel = pathlib.Path("misc") / well_dir.name
        out_csv = out_root / "by_well" / rel / "metrics.csv"
        out_csv.parent.mkdir(parents=True, exist_ok=True)
        pl.from_dicts(recs).write_csv(out_csv)
        print(f"[wrote] {out_csv} ({len(recs)} cells)", flush=True)

    all_recs = [r for recs in records_by_well.values() for r in recs]
    if all_recs:
        combined = out_root / "combined_raw.csv"
        pl.from_dicts(all_recs).write_csv(combined)
        print(f"[wrote] {combined} ({len(all_recs)} cells total)", flush=True)


if __name__ == "__main__":
    main()
