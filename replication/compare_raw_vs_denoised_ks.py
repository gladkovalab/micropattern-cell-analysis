"""Denoising sensitivity check on 10 diagnostic cells, focused on the two
new headline metrics: wedge_r_gini and wedge_r_ks_vs_uniform.

Re-runs `final_pipeline.process_cell` on raw and denoised ND2s for the same
10 cells and reports per-cell deltas. The 10 cells are stratified across
sheets and conditions so the result is representative of the full pitch.

Output:
  - replication/overnight_final_out/raw_vs_denoised_ks.csv
  - markdown table to stdout
"""
from __future__ import annotations
import pathlib, sys
import numpy as np
import polars as pl

REPO = pathlib.Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))
from replication.final_pipeline import process_cell  # noqa: E402
import template_matching_bulk as tmb  # noqa: E402

OUT = REPO / "replication" / "overnight_final_out" / "raw_vs_denoised_ks.csv"
ROOT = pathlib.Path("/Volumes/valelab/_for_Mark/patterned_data")

# 10 diagnostic cells, stratified across sheets and conditions.
# Format: (plate, well_dir, condition, cell_stem, denoised_suffix)
# denoised_suffix is " - Denoised.nd2" except plate_12 which uses
# " - Denoised2.nd2".
CELLS = [
    # TRAK isoform (mito) — 3 conditions
    ("250612_patterned_plate_3", "B02_NoV_250616",         "TRAK iso mito · no TRAK", "cell5",  " - Denoised.nd2"),
    ("250612_patterned_plate_3", "B03_TRAK1_250616",       "TRAK iso mito · TRAK1",   "cell7",  " - Denoised.nd2"),
    ("250612_patterned_plate_3", "B04_TRAK2_250616",       "TRAK iso mito · TRAK2",   "cell1",  " - Denoised.nd2"),
    # TRAK1 helix muts — wt + mDRH
    ("250612_patterned_plate_3", "B06_250617_TRAK1_mDRH_dSp", "TRAK1 muts · mDRH/dSp", "Cell1", " - Denoised.nd2"),
    ("250731_patterned_plate_11_good", "E05_250808_TRAK1_wt", "TRAK1 muts · wt",      "Cell2",  " - Denoised.nd2"),
    # TRAK2 helix muts — wt + Spindly rescue
    ("250731_patterned_plate_11_good", "F05_250808_TRAK2_wt", "TRAK2 muts · wt",      "Cell1",  " - Denoised.nd2"),
    ("250731_patterned_plate_11_good", "E07_250811_TRAK2_mDRH_Sp", "TRAK2 muts · mDRH+mSp", "Cell3", " - Denoised.nd2"),
    # MAPK9 siRNA — ctrl-ctrl, ctrl-Ars, MAPK9-Ars (plate 12 uses 'Denoised.nd2' suffix on these wells)
    ("250807_patterned_plate_12", "F04_250814_ctrl_siRNA_ctrl",  "MAPK9 · ctrl ctrl", "cell11", " - Denoised.nd2"),
    ("250807_patterned_plate_12", "B04_250813_ctrl_siRNA_Ars",   "MAPK9 · ctrl Ars",  "Cell10", " - Denoised.nd2"),
    ("250807_patterned_plate_12", "C04_250814_MAPK9_siRNA_Ars",  "MAPK9 · MAPK9 Ars", "cell1",  " - Denoised.nd2"),
]

# Headline metrics — both projections of the two new clustering metrics
HEADLINE = [
    "zsum_wedge_r_gini",          "maxip_wedge_r_gini",
    "zsum_wedge_r_ks_vs_uniform", "maxip_wedge_r_ks_vs_uniform",
]


def find_existing(parent, well_prefix):
    """Find the actual well dir matching the prefix on disk."""
    for d in parent.iterdir():
        if d.is_dir() and d.name.startswith(well_prefix):
            return d
    return None


def main():
    template_hat = tmb.get_template_hat(1326)
    template = tmb.get_padded_template_at_width(1326)

    rows = []
    for plate, well_prefix, label, stem, suffix in CELLS:
        plate_dir = ROOT / plate
        well_dir = find_existing(plate_dir, well_prefix.split("_")[0] + "_")
        if well_dir is None:
            print(f"  MISSING well: {plate}/{well_prefix}")
            continue
        raw_path = well_dir / f"{stem}.nd2"
        den_path = well_dir / "denoised" / f"{stem}{suffix}"
        if not raw_path.exists():
            # Try common stem variants — different plates use cell vs Cell
            for alt in (stem.replace("cell", "Cell"), stem.replace("Cell", "cell")):
                if (well_dir / f"{alt}.nd2").exists():
                    raw_path = well_dir / f"{alt}.nd2"
                    den_path = well_dir / "denoised" / f"{alt}{suffix}"
                    break
        if not raw_path.exists() or not den_path.exists():
            print(f"  MISSING: {raw_path.name} or {den_path.name} in {well_dir.name}")
            continue
        print(f"\n--- {label} ({plate.split('_plate_')[-1]} {well_dir.name.split('_')[0]} {raw_path.stem}) ---")
        try:
            m_raw = process_cell(raw_path, template_hat=template_hat, template=template)
            m_den = process_cell(den_path, template_hat=template_hat, template=template)
        except Exception as e:
            print(f"  ERR: {e}")
            continue
        row = {"plate": plate.split("_plate_")[-1],
               "well": well_dir.name.split("_")[0],
               "label": label, "cell": raw_path.stem}
        for k in HEADLINE:
            r, d = m_raw.get(k), m_den.get(k)
            row[f"raw_{k}"] = r
            row[f"den_{k}"] = d
            if r is not None and d is not None:
                row[f"delta_{k}"] = d - r
                row[f"pct_{k}"] = (d - r) / r * 100 if abs(r) > 1e-9 else None
            print(f"  {k:<32s} raw={r:.4f}  den={d:.4f}  Δ={d-r:+.4f} "
                  f"({(d-r)/r*100:+.2f}%)" if r else f"  {k}: missing")
        rows.append(row)

    if not rows:
        print("No cells processed!")
        return

    OUT.parent.mkdir(parents=True, exist_ok=True)
    pl.from_dicts(rows).write_csv(OUT)
    print(f"\nWrote {OUT}")

    # Summary
    print("\n" + "=" * 78)
    print("Summary across {} cells: |Δ| per-cell on the four headline metrics".format(len(rows)))
    print("=" * 78)
    print(f"{'Metric':<32} {'mean Δ':>9} {'sd Δ':>8} {'max |Δ|':>9}  {'mean |%|':>9}")
    print("-" * 78)
    for k in HEADLINE:
        deltas = np.array([r.get(f"delta_{k}") for r in rows
                           if r.get(f"delta_{k}") is not None], dtype=float)
        pcts = np.array([abs(r.get(f"pct_{k}")) for r in rows
                         if r.get(f"pct_{k}") is not None], dtype=float)
        if deltas.size == 0:
            continue
        print(f"{k:<32} {deltas.mean():>+9.4f} {deltas.std():>8.4f} "
              f"{np.abs(deltas).max():>9.4f}  {pcts.mean():>9.2f}")


if __name__ == "__main__":
    main()
