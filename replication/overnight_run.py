"""Overnight metric pipeline: runs metric_pipeline.process_cell on every cell
in the Fig 4 / Fig S11 target sheets, writing per-well CSVs as they complete
and a progress log with timestamps so we can check status in the morning.

Design notes:
- Sequential, one cell at a time. Initial attempts at parallel SMB reads
  saturated the network; serial is not slower in practice.
- Checkpoints per well: the per-well CSV is written immediately after the
  last cell in that well finishes. If the run is interrupted, re-running
  will skip any well whose CSV already exists (idempotent resume).
- Per-cell errors are logged and the run continues (e.g., template-match
  failures on a single cell do not abort the run).
- Reads from mark_data/ (read-only) only; writes to replication/overnight_out/.
"""
from __future__ import annotations

import datetime
import pathlib
import sys
import time
import traceback
from dataclasses import dataclass

import fastexcel
import nd2  # noqa: F401  (import here so failures surface immediately)
import polars as pl

REPO = pathlib.Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / "replication"))

import template_matching_bulk as tmb
from metric_pipeline import process_cell, iter_target_nd2s


OUT_ROOT = REPO / "replication" / "overnight_out"
LOG_PATH = OUT_ROOT / "progress.log"

# Panels whose wells we want to cover. One entry per (sheet_name, conditions).
# Conditions listed here match the column headers in Comparisons_table_v3.xlsx.
# Only conditions listed will be processed (so we skip sheet columns that
# don't belong to our target panels).
TARGET_SHEETS: dict[str, list[str]] = {
    "TRAK isoform (mito)": ["no TRAK", "TRAK1", "TRAK2"],
    "TRAK1 helix muts": ["T1 wt", "T1 mDRH", "T1 mDRH / dSp"],
    "TRAK2 helix muts": ["TRAK2", "TRAK2 mDRH", "TRAK2 mDRH mSpindly"],
    "MAPK9 siRNA": ["ctrl ctrl", "ctrl Ars", "MAPK9 ctrl", "MAPK9 Ars"],
}

COMPARISONS_XLSX = REPO / "config" / "Comparisons_table_v3.xlsx"


@dataclass
class WellTarget:
    sheet: str
    condition: str
    plate: str
    well: str
    raw_dir: pathlib.Path


def _log(msg: str) -> None:
    ts = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    line = f"[{ts}] {msg}"
    print(line, flush=True)
    LOG_PATH.parent.mkdir(parents=True, exist_ok=True)
    with LOG_PATH.open("a") as f:
        f.write(line + "\n")


def collect_wells() -> list[WellTarget]:
    """Read Comparisons_table_v3.xlsx and resolve each (sheet, plate, well) to
    an ND2 directory under the MICROPATTERN_DATA_ROOT tree."""
    data_root = pathlib.Path(tmb.DATA_ROOT).resolve()
    out: list[WellTarget] = []
    reader = fastexcel.read_excel(COMPARISONS_XLSX)
    for sheet_name in reader.sheet_names:
        if sheet_name not in TARGET_SHEETS:
            continue
        df = pl.read_excel(COMPARISONS_XLSX, sheet_name=sheet_name)
        plate_col = df.columns[0]
        wanted_conds = TARGET_SHEETS[sheet_name]
        for cond in df.columns[1:]:
            if cond not in wanted_conds:
                continue
            for row in df.iter_rows(named=True):
                plate = row[plate_col]
                well = row[cond]
                if not well or not plate:
                    continue
                plate_dir = data_root / plate
                if not plate_dir.is_dir():
                    _log(f"  [miss] plate dir not found: {plate_dir}")
                    continue
                candidate = None
                for sub in plate_dir.iterdir():
                    if sub.is_dir() and sub.name.startswith(well + "_"):
                        candidate = sub
                        break
                if candidate is None:
                    _log(f"  [miss] no well dir starting with {well}_ in {plate_dir}")
                    continue
                out.append(WellTarget(sheet=sheet_name, condition=cond, plate=plate,
                                       well=well, raw_dir=candidate))
    return out


def process_well(target: WellTarget, template_hat, template) -> tuple[int, int]:
    """Process every ND2 under target.raw_dir. Writes per-cell CSVs as they
    complete (so a mid-well crash doesn't lose already-processed cells) and
    returns (n_ok, n_err). The well is 'done' (skippable on resume) when its
    `done.marker` exists, written only after all cells attempted."""
    rel = target.raw_dir.resolve().relative_to(pathlib.Path(tmb.DATA_ROOT).resolve())
    well_out_dir = OUT_ROOT / "by_well" / rel
    well_out_dir.mkdir(parents=True, exist_ok=True)
    done_marker = well_out_dir / "done.marker"
    well_csv = well_out_dir / "metrics.csv"
    cells_dir = well_out_dir / "cells"
    cells_dir.mkdir(parents=True, exist_ok=True)

    if done_marker.exists():
        _log(f"  SKIP (already done): {rel}")
        return 0, 0
    # Back-compat: pre-patch runs wrote well_csv without a done.marker and
    # without a cells/ subdir. Honour those as completed and fast-forward.
    has_old_csv = well_csv.exists() and not any(cells_dir.iterdir())
    if has_old_csv:
        done_marker.touch()
        _log(f"  SKIP (pre-patch well_csv complete): {rel}")
        return 0, 0

    cells = list(iter_target_nd2s(target.raw_dir))
    if not cells:
        _log(f"  (no cells under {rel})")
        done_marker.touch()
        return 0, 0

    n_ok = n_err = 0
    for i, img_path in enumerate(cells, 1):
        per_cell_csv = cells_dir / f"{img_path.stem}.csv"
        if per_cell_csv.exists():
            n_ok += 1
            _log(f"  SKIP {i}/{len(cells)} {img_path.name} (cell already done)")
            continue
        t0 = time.time()
        try:
            m = process_cell(img_path, template_hat=template_hat, template=template,
                             out_root=OUT_ROOT, save_projections=False)
            m["sheet"] = target.sheet
            m["condition"] = target.condition
            m["plate"] = target.plate
            m["well"] = target.well
            # Atomic write: tmp → rename
            tmp = per_cell_csv.with_suffix(".csv.tmp")
            pl.from_dicts([m]).write_csv(tmp)
            tmp.rename(per_cell_csv)
            n_ok += 1
            elapsed = time.time() - t0
            _log(f"  OK  {i}/{len(cells)}  {img_path.name}  ({elapsed:.1f}s)")
        except Exception as e:
            n_err += 1
            _log(f"  ERR {i}/{len(cells)}  {img_path.name}: {e}")
            traceback.print_exc()

    # Aggregate per-cell CSVs into the well CSV
    per_cell_files = sorted(cells_dir.glob("*.csv"))
    if per_cell_files:
        frames = [pl.read_csv(p) for p in per_cell_files]
        pl.concat(frames, how="diagonal_relaxed").write_csv(well_csv)
        _log(f"  wrote {well_csv.relative_to(OUT_ROOT)} ({len(per_cell_files)} cells)")
    done_marker.touch()
    return n_ok, n_err


def aggregate_combined() -> None:
    """Combine every well CSV written so far into a single combined CSV."""
    rows = []
    for csv in (OUT_ROOT / "by_well").rglob("metrics.csv"):
        try:
            df = pl.read_csv(csv)
            rows.extend(df.to_dicts())
        except Exception as e:
            _log(f"  aggregate skip {csv}: {e}")
    if rows:
        combined = OUT_ROOT / "combined.csv"
        pl.from_dicts(rows).write_csv(combined)
        _log(f"  aggregated → {combined} ({len(rows)} cells)")


def main() -> int:
    OUT_ROOT.mkdir(parents=True, exist_ok=True)
    _log("=" * 60)
    _log("overnight_run starting")
    _log(f"  DATA_ROOT = {tmb.DATA_ROOT}")
    _log(f"  OUT_ROOT  = {OUT_ROOT}")

    targets = collect_wells()
    _log(f"  {len(targets)} target wells across {len(TARGET_SHEETS)} sheets")
    for t in targets[:3]:
        _log(f"    sample: {t.sheet} / {t.plate} / {t.well} → {t.raw_dir}")

    _log("initializing template…")
    template_hat = tmb.get_template_hat(1326)
    template = tmb.get_padded_template_at_width(1326)

    t_start = time.time()
    total_ok = total_err = 0
    for i, target in enumerate(targets, 1):
        _log("-" * 60)
        _log(f"[{i}/{len(targets)}] {target.sheet} · {target.plate} · {target.well}")
        ok, err = process_well(target, template_hat, template)
        total_ok += ok
        total_err += err
        aggregate_combined()  # refresh combined.csv after every well
        elapsed_h = (time.time() - t_start) / 3600
        _log(f"  running totals: ok={total_ok} err={total_err} elapsed={elapsed_h:.2f}h")

    _log("=" * 60)
    _log(f"DONE. total ok={total_ok}, err={total_err}, wall={(time.time()-t_start)/3600:.2f}h")
    return 0


if __name__ == "__main__":
    sys.exit(main())
