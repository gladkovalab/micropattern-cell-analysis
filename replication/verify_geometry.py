"""Verify geometry on 20 random Fig 4B cells using `final_pipeline.py`'s helpers.

Picks 20 cells stratified across the three Fig 4B conditions and 5 plates,
opens each ND2, runs the SAME template_matching → cropping → nuclear seg →
pattern extreme → wedge construction logic that `final_pipeline.process_cell`
will run, and renders one PNG per cell with:
  * 488 (mito) signal as backdrop
  * Orange pattern outline (find_contours of shifted_template, in crop frame)
  * Cyan nucleus segmentation contour
  * Four pattern extreme points (red bottom, white top, yellow left, magenta right)
  * Two wedge-boundary rays from apex through left and right
  * Light-red overlay on pixels OUTSIDE the wedge

Output: replication/overnight_final_out/figures/verify_geometry_4B/<idx>_<...>.png

Run with the SMB mount available:
  MICROPATTERN_DATA_ROOT=/Volumes/valelab/_for_Mark/patterned_data \\
    pixi run python replication/verify_geometry.py
"""
from __future__ import annotations

import pathlib
import re
import sys

import nd2
import numpy as np
import polars as pl
import skimage
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.ndimage import center_of_mass

REPO = pathlib.Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))

import template_matching_bulk as tmb  # noqa: E402
from replication.final_pipeline import (  # noqa: E402
    _pattern_extremes, _build_wedge_geometry,
)

OLD_CSV = REPO / "replication" / "overnight_out" / "combined.csv"
OUT_DIR = REPO / "replication" / "overnight_final_out" / "figures" / "verify_geometry_4B"
SHEET = "TRAK isoform (mito)"
N_CELLS = 20
SEED = 4242  # so the picks are reproducible


def slugify(s: str) -> str:
    return re.sub(r"[^A-Za-z0-9]+", "_", s).strip("_")


def pick_stratified(df: pl.DataFrame, n: int) -> list[dict]:
    """Stratified sample across (plate, condition) cells. Falls back to any
    available stratum if some have too few cells."""
    rng = np.random.default_rng(SEED)
    sub = df.filter(pl.col("sheet") == SHEET)
    by_stratum: dict[tuple, list[dict]] = {}
    for r in sub.iter_rows(named=True):
        key = (r["plate"], r["condition"])
        by_stratum.setdefault(key, []).append(
            {"plate": r["plate"], "well": r["well"],
             "condition": r["condition"], "path": r["path"]})
    # round-robin one cell at a time
    out: list[dict] = []
    keys = list(by_stratum.keys())
    rng.shuffle(keys)
    while len(out) < n:
        progress = False
        for k in keys:
            if not by_stratum[k]:
                continue
            i = int(rng.integers(0, len(by_stratum[k])))
            out.append(by_stratum[k].pop(i))
            progress = True
            if len(out) >= n:
                break
        if not progress:
            break
    return out


def render_one(ax, img_path, template_hat, template):
    key = tmb.cluster_key(img_path)
    offset = tmb.offset_overrides.get(key, [128, 128])
    roi = tmb.roi_overrides.get(key, None)

    img = nd2.imread(img_path, xarray=True)
    zsum = img.sum(axis=0)

    max_coords = tmb.get_template_center(img, img_path, template_hat=template_hat,
                                         offset=offset, roi=roi)
    shifted_template = np.roll(template, (max_coords[0] - 1024, max_coords[1] - 1024),
                               axis=(0, 1))

    y_start = max_coords[0] - 512 + offset[0]
    y_end = max_coords[0] + 512 + offset[0]
    x_start = max_coords[1] - 512 + offset[1]
    x_end = max_coords[1] + 512 + offset[1]
    cropped_zsum = zsum.isel(Y=slice(y_start, y_end), X=slice(x_start, x_end))
    nuc = cropped_zsum.sel(C="405").to_numpy()
    mito = cropped_zsum.sel(C="488").to_numpy()

    # Orange contour (authoritative pattern outline in crop frame)
    contour = skimage.measure.find_contours(shifted_template)[0].copy()
    contour[:, 0] -= max_coords[0] - 512
    contour[:, 1] -= max_coords[1] - 512

    # Pattern extremes — same call as in final_pipeline.process_cell
    ext = _pattern_extremes(shifted_template, max_coords)

    # Nucleus segmentation
    nuc_mask = nuc > skimage.filters.threshold_otsu(nuc)
    nl = skimage.measure.label(nuc_mask)
    nprops = skimage.measure.regionprops(nl)
    nmax = int(np.argmax([p.area for p in nprops])) + 1
    nuc_mask = nl == nmax
    nuc_com = center_of_mass(nuc_mask)

    pitch = img.metadata["metadata"].channels[0].volume.axesCalibration[0]

    # Wedge geometry — same call as in final_pipeline.process_cell
    wedge_mask, _, opening_rad, _, _ = _build_wedge_geometry(nuc.shape, ext, pitch)

    vmax = np.percentile(mito, 99.5)
    ax.imshow(mito, cmap="Greens_r", vmax=vmax)

    # Wedge: dim outside in red
    outside = ~wedge_mask
    red_overlay = np.zeros((*wedge_mask.shape, 4))
    red_overlay[outside] = [1, 0, 0, 0.22]
    ax.imshow(red_overlay)

    # Pattern outline
    ax.plot(contour[:, 1], contour[:, 0], color="orange", linewidth=1.0,
            alpha=0.95)
    # Nucleus
    for c in skimage.measure.find_contours(nuc_mask.astype(float), 0.5):
        ax.plot(c[:, 1], c[:, 0], color="cyan", linewidth=1.0)

    # Pattern extreme points
    for name, color, marker in [("left", "yellow", "<"), ("right", "magenta", ">"),
                                ("top", "white", "^"), ("bottom", "red", "v")]:
        pt = ext[name]
        ax.scatter([pt[1]], [pt[0]], color=color, s=55, marker=marker,
                   edgecolor="black", linewidth=0.6, zorder=5)

    # Wedge boundary rays from apex through left + right
    apex = ext["bottom"]
    H, W = mito.shape
    for name, color in [("left", "yellow"), ("right", "magenta")]:
        pt = ext[name]
        dy = pt[0] - apex[0]
        dx = pt[1] - apex[1]
        ey = apex[0] + dy * 6.0
        ex = apex[1] + dx * 6.0
        ax.plot([apex[1], ex], [apex[0], ey], color=color, linewidth=1.0,
                linestyle="--", alpha=0.85, zorder=4)

    ax.set_xlim(0, W - 1)
    ax.set_ylim(H - 1, 0)
    ax.axis("off")
    return float(np.degrees(opening_rad))


def main():
    df = pl.read_csv(OLD_CSV)
    sample = pick_stratified(df, N_CELLS)
    print(f"Picked {len(sample)} cells:")
    for s in sample:
        print(f"  {s['plate']}/{s['well']}/{pathlib.Path(s['path']).name} ({s['condition']})")

    template_hat = tmb.get_template_hat(1326)
    template = tmb.get_padded_template_at_width(1326)

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    for i, cell in enumerate(sample, 1):
        fig, ax = plt.subplots(1, 1, figsize=(8, 8))
        try:
            opening = render_one(ax, pathlib.Path(cell["path"]),
                                 template_hat, template)
            label = (f"Fig 4B verify {i:02d}/20 · "
                     f"{cell['plate'].split('_patterned_plate_')[-1]} "
                     f"{cell['well']} · {cell['condition']} · "
                     f"{pathlib.Path(cell['path']).stem} · "
                     f"wedge {opening:.0f}°\n"
                     "488 signal · orange = pattern · cyan = nucleus · "
                     "red v / white ^ / yellow < / magenta > = extremes · "
                     "red shading = excluded from wedge")
            ax.set_title(label, fontsize=8)
        except Exception as e:
            ax.text(0.5, 0.5, f"ERR\n{e}", transform=ax.transAxes,
                    ha="center", va="center", fontsize=10)
            ax.axis("off")

        plate_short = cell["plate"].split("_patterned_plate_")[-1]
        out_name = (f"{i:02d}_p{plate_short}_{cell['well']}_"
                    f"{slugify(cell['condition'])}_"
                    f"{pathlib.Path(cell['path']).stem}.png")
        out_path = OUT_DIR / out_name
        plt.tight_layout()
        fig.savefig(out_path, dpi=130, bbox_inches="tight")
        plt.close(fig)
        print(f"  [{i:02d}/{len(sample)}] {out_path}")

    print(f"\n{len(sample)} images written to {OUT_DIR}/")


if __name__ == "__main__":
    main()
