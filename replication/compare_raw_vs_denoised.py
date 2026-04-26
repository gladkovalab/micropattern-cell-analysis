"""Compare raw vs denoised inputs on the same 5 Fig 4B cells.

Runs `final_pipeline.process_cell` twice per cell — once on the raw ND2,
once on the denoised version (`<well>/denoised/<cellname> - Denoised.nd2`).
Reports a side-by-side table of the keeper metrics (Y-Gini, wedge-r-Gini,
plus a few cross-checks) for both projections.

Picks 5 cells across plates/conditions deterministically. Output is a
markdown table printed to stdout plus a CSV at
`replication/overnight_final_out/raw_vs_denoised.csv`.
"""
from __future__ import annotations

import pathlib
import sys

import polars as pl

REPO = pathlib.Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))

from replication.final_pipeline import process_cell  # noqa: E402
import template_matching_bulk as tmb  # noqa: E402

OUT = REPO / "replication" / "overnight_final_out" / "raw_vs_denoised.csv"

# 5 cells stratified across plates and conditions
CELLS = [
    ("250612_patterned_plate_3", "B02_NoV_250616",        "no TRAK", "cell5"),
    ("250612_patterned_plate_3", "B04_TRAK2_250616",      "TRAK2",   "cell1"),
    ("250612_patterned_plate_3", "B03_TRAK1_250616",      "TRAK1",   "cell7"),
    ("250710_patterned_plate_9_good", "C04_250718_TRAK2", "TRAK2",   "Cell5"),
    ("250731_patterned_plate_11_good", "F05_250808_TRAK2_wt", "TRAK2", "Cell9"),
]

KEEPER = [
    "zsum_y_gini", "maxip_y_gini",
    "zsum_wedge_r_gini", "maxip_wedge_r_gini",
    "zsum_wedge_r_entropy", "maxip_wedge_r_entropy",
    "zsum_wedge_r_sd_um", "maxip_wedge_r_sd_um",
    "zsum_perinuclear_5um_pct", "maxip_perinuclear_5um_pct",
    "zsum_peripheral_5um_pct", "maxip_peripheral_5um_pct",
    "wedge_opening_deg",
    "nuc_area_um2", "nuc_solidity",
]


def main():
    template_hat = tmb.get_template_hat(1326)
    template = tmb.get_padded_template_at_width(1326)
    root = pathlib.Path("/Volumes/valelab/_for_Mark/patterned_data")

    rows = []
    for plate, well_dir_name, cond, cell_stem in CELLS:
        raw_path = root / plate / well_dir_name / f"{cell_stem}.nd2"
        den_path = root / plate / well_dir_name / "denoised" / f"{cell_stem} - Denoised.nd2"
        if not raw_path.exists() or not den_path.exists():
            print(f"MISSING: {raw_path} or {den_path}")
            continue
        print(f"\n--- {plate.split('_patterned_plate_')[-1]} {well_dir_name.split('_')[0]} "
              f"{cond} {cell_stem} ---")
        try:
            m_raw = process_cell(raw_path, template_hat=template_hat, template=template)
            m_den = process_cell(den_path, template_hat=template_hat, template=template)
        except Exception as e:
            print(f"  ERR: {e}")
            continue
        row = {
            "plate": plate.split("_patterned_plate_")[-1],
            "well": well_dir_name.split("_")[0],
            "condition": cond,
            "cell": cell_stem,
        }
        for k in KEEPER:
            row[f"raw_{k}"] = m_raw.get(k)
            row[f"den_{k}"] = m_den.get(k)
            if m_raw.get(k) is not None and m_den.get(k) is not None:
                row[f"delta_{k}"] = m_den[k] - m_raw[k]
        rows.append(row)
        # Pretty print delta on the keeper metrics
        for k in KEEPER:
            r, d = m_raw.get(k), m_den.get(k)
            if r is None or d is None:
                continue
            try:
                pct = (d - r) / r * 100 if abs(r) > 1e-9 else float("nan")
            except Exception:
                pct = float("nan")
            print(f"  {k:<35s} raw={r:>10.4f}  den={d:>10.4f}  Δ={d-r:+8.4f} ({pct:+6.1f}%)")

    if rows:
        OUT.parent.mkdir(parents=True, exist_ok=True)
        pl.from_dicts(rows).write_csv(OUT)
        print(f"\nWrote {OUT}")

        # Print summary table for the headline metrics
        df = pl.from_dicts(rows)
        print("\n" + "=" * 100)
        print("Summary: per-cell delta on the headline metrics (denoised − raw)")
        print("=" * 100)
        for k in ["maxip_wedge_r_gini", "zsum_wedge_r_gini",
                  "maxip_y_gini", "zsum_y_gini",
                  "maxip_perinuclear_5um_pct", "zsum_perinuclear_5um_pct"]:
            raw_col = f"raw_{k}"
            den_col = f"den_{k}"
            if raw_col not in df.columns or den_col not in df.columns:
                continue
            print(f"\n{k}:")
            for r in df.iter_rows(named=True):
                rv = r.get(raw_col); dv = r.get(den_col)
                if rv is None or dv is None:
                    continue
                pct = (dv - rv) / rv * 100 if abs(rv) > 1e-9 else float("nan")
                lbl = f"{r['plate']:>7} {r['well']:<4} {r['condition']:<8} {r['cell']:<8}"
                print(f"  {lbl}  raw={rv:>9.4f}  den={dv:>9.4f}  Δ={dv-rv:+8.4f} ({pct:+6.1f}%)")


if __name__ == "__main__":
    main()
