"""Export the wedge-r 1D profile source data as a long-format CSV.

For every (sheet, condition, radial bin), writes the per-cell mean and
SEM of:
  - 488 mitochondria intensity profile (% per bin)        ['mito']
  - 405 nuclear mask radial distribution (% per bin)      ['nuc_mask']
  - 5 µm perinuclear halo radial distribution (% per bin) ['perinuc_halo']

These are exactly the curves drawn in the per-sheet split figures
(profiles_*_with_nuclear.png) and the bands-overlay figure.
"""
from __future__ import annotations

import argparse
import pathlib
import sys

import numpy as np
import polars as pl

REPO = pathlib.Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / "replication"))
from plot_metrics import (  # noqa: E402
    SHEET_CONFIG, load_template_matching, join_with_metadata,
)
import plot_all_with_nuclear as pawn  # noqa: E402


def mito_cols() -> list[str]:
    return [f"wedge_r_{i:02d}_{i+1:02d}um_pct" for i in range(60)]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out",
                    default="replication/figures_wedge_r_ks/wedge_r_profiles_source.csv")
    args = ap.parse_args()

    df = load_template_matching(pathlib.Path(
        "replication/wedge_r_ks_out_all_denoised/by_well"))
    df = join_with_metadata(df, REPO / "config/Comparisons_table_v3.xlsx")

    mcols = [c for c in mito_cols() if c in df.columns]
    n_bins = len(mcols)

    cache: dict[str, tuple[np.ndarray, np.ndarray]] = {}

    rows: list[dict] = []
    for sheet, cfg in SHEET_CONFIG.items():
        sheet_df = df.filter(pl.col("sheet") == sheet)
        if sheet_df.height == 0:
            print(f"  skipping {sheet}: no rows")
            continue
        print(f"  {sheet}: {sheet_df.height} cells")
        for cond in cfg["conditions"]:
            sub = sheet_df.filter(pl.col("condition") == cond)
            if sub.height == 0:
                continue

            # --- Mito profile (per-bin %), already per-cell in the CSV.
            mito = sub.select(mcols).to_numpy().astype(float)
            n_mito = mito.shape[0]

            # --- Nuclear mask & perinuclear halo profiles, computed from
            # the saved 405 Z-sum projections (cached across sheets).
            nuc_rows, halo_rows = [], []
            for path in sub["path"].to_list():
                if path not in cache:
                    cache[path] = pawn.nuclear_mask_profiles(path)
                np_, hp_ = cache[path]
                nuc_rows.append(np_); halo_rows.append(hp_)
            nuc = np.vstack(nuc_rows) if nuc_rows else np.zeros((0, n_bins))
            halo = np.vstack(halo_rows) if halo_rows else np.zeros((0, n_bins))
            valid = (np.isfinite(nuc).all(axis=1) if nuc.size
                     else np.array([], bool))
            nuc_v = nuc[valid] if nuc.size else nuc
            halo_v = halo[valid] if halo.size else halo
            n_nuc = nuc_v.shape[0]

            for b in range(n_bins):
                m_mean = float(np.nanmean(mito[:, b])) if n_mito else float("nan")
                m_sem = (float(np.nanstd(mito[:, b], ddof=1) /
                               max(np.sqrt(n_mito), 1.0))
                         if n_mito > 1 else float("nan"))
                if n_nuc:
                    n_mean = float(np.nanmean(nuc_v[:, b]))
                    h_mean = float(np.nanmean(halo_v[:, b]))
                    n_sem = (float(np.nanstd(nuc_v[:, b], ddof=1) /
                                   max(np.sqrt(n_nuc), 1.0))
                             if n_nuc > 1 else float("nan"))
                    h_sem = (float(np.nanstd(halo_v[:, b], ddof=1) /
                                   max(np.sqrt(n_nuc), 1.0))
                             if n_nuc > 1 else float("nan"))
                else:
                    n_mean = h_mean = n_sem = h_sem = float("nan")

                rows.append({
                    "sheet": sheet,
                    "condition": cond,
                    "n_cells_mito": n_mito,
                    "n_cells_nuc": n_nuc,
                    "bin_lo_um": b,
                    "bin_hi_um": b + 1,
                    "bin_center_um": b + 0.5,
                    "mito_mean_pct": m_mean,
                    "mito_sem_pct": m_sem,
                    "nuc_mask_mean_pct": n_mean,
                    "nuc_mask_sem_pct": n_sem,
                    "perinuc_halo_mean_pct": h_mean,
                    "perinuc_halo_sem_pct": h_sem,
                })

    out = pl.from_dicts(rows)
    out_path = pathlib.Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out.write_csv(out_path, float_precision=6)
    print(f"\nwrote {out_path}  ({out.height} rows · "
          f"{out['sheet'].n_unique()} sheets · "
          f"{out['condition'].n_unique()} unique conditions)")


if __name__ == "__main__":
    main()
