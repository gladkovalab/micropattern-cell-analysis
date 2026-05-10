"""Offline twin of plot_wedge_illustration.py — renders the wedge geometry
on a real cell using the saved cropped projections (no ND2 / SMB needed).

The pipeline saves cropped 1024×1024 multi-channel Z-sum projections under
`projections/...<cell>.nc` and the bg-subtracted 488 MaxIP under
`projections/...<cell>_488_bg_subtracted.nc`. Both are already template-
aligned, so the wedge constants from template_matching_bulk apply
directly without re-running template matching.

Default cell: a TRAK2 mito plate-3 cell with strong perinuclear pile-up,
chosen so the wedge ROI is visually obvious.

CLI:
    pixi run python analysis/plot_wedge_illustration_offline.py
"""
from __future__ import annotations

import argparse
import pathlib
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import skimage
import xarray as xr

REPO = pathlib.Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))
import template_matching_bulk as tmb  # noqa: E402

PITCH_UM = 0.065  # constant across the dataset


def render(zsum_nc: pathlib.Path, mito_nc: pathlib.Path, out_path: pathlib.Path):
    zsum = xr.open_dataset(zsum_nc)
    nuc = zsum.sel(C="405").to_array().squeeze().values.reshape((1024, 1024))

    if mito_nc.exists():
        mito_da = xr.open_dataset(mito_nc).to_array().squeeze()
        mito = mito_da.values
        if mito.ndim != 2:
            mito = mito.reshape(mito.shape[-2:])
        mito_label = "488 MaxIP, bg-subtracted (metric input)"
    else:
        mito = zsum.sel(C="488").to_array().squeeze().values.reshape((1024, 1024))
        mito_label = "488 Z-sum (bg-subtracted MaxIP not on disk)"

    # Pattern contour: regenerate the template directly so we don't need
    # the original ND2 or template-matching state. The pipeline crop puts
    # the template centre at (512, 512) by construction, so the template
    # contour translates to the cropped frame by subtracting (1024-512).
    template = tmb.get_padded_template_at_width(1326)
    contour = skimage.measure.find_contours(template)[0].copy()
    contour[:, 0] -= 1024 - 512
    contour[:, 1] -= 1024 - 512

    # Nucleus segmentation (same recipe as the pipeline).
    nuc_mask = nuc > skimage.filters.threshold_otsu(nuc)
    label = skimage.measure.label(nuc_mask)
    props = skimage.measure.regionprops(label)
    if props:
        largest = int(np.argmax([p.area for p in props])) + 1
        nuc_mask = label == largest

    wedge_mask, _, _, _, _ = tmb._get_wedge_geometry(mito.shape, PITCH_UM)
    apex = tmb.WEDGE_APEX
    left = tmb.WEDGE_LEFT
    right = tmb.WEDGE_RIGHT
    a_left = float(np.arctan2(left[0] - apex[0], left[1] - apex[1]))
    a_right = float(np.arctan2(right[0] - apex[0], right[1] - apex[1]))
    opening_deg = abs(np.degrees(a_left - a_right))
    wedge_frac = 100 * wedge_mask.sum() / wedge_mask.size

    fig, ax = plt.subplots(figsize=(8.5, 8.5))
    vmax = np.percentile(mito, 99.5)
    ax.imshow(mito, cmap="Greens_r", vmax=vmax if vmax > 0 else None)
    # Translucent red outside the wedge.
    outside = ~wedge_mask
    overlay = np.zeros((*wedge_mask.shape, 4))
    overlay[outside] = [1, 0, 0, 0.32]
    ax.imshow(overlay)
    # Pattern outline.
    ax.plot(contour[:, 1], contour[:, 0], color="orange", linewidth=1.6,
            alpha=0.95, label="pattern outline")
    # Nucleus contour.
    for c in skimage.measure.find_contours(nuc_mask.astype(float), 0.5):
        ax.plot(c[:, 1], c[:, 0], color="cyan", linewidth=1.4,
                label="_nuc" if c is not None else None)
    ax.plot([], [], color="cyan", linewidth=1.4, label="nucleus (Otsu)")
    # Apex / L / R anchors.
    ax.scatter([apex[1]], [apex[0]], color="red", s=160, marker="v",
               edgecolor="black", linewidth=1, zorder=5,
               label=f"apex {apex}")
    ax.scatter([left[1]], [left[0]], color="yellow", s=130, marker="<",
               edgecolor="black", linewidth=1, zorder=5,
               label=f"L tangent {left}")
    ax.scatter([right[1]], [right[0]], color="magenta", s=130, marker=">",
               edgecolor="black", linewidth=1, zorder=5,
               label=f"R tangent {right}")
    # Bounding rays from apex through L and R, extended to image bounds.
    H, W = mito.shape
    for pt, color in ((left, "yellow"), (right, "magenta")):
        dy = pt[0] - apex[0]
        dx = pt[1] - apex[1]
        ey, ex = apex[0] + dy * 6.0, apex[1] + dx * 6.0
        ax.plot([apex[1], ex], [apex[0], ey], color=color, linewidth=1.8,
                linestyle="--", alpha=0.95, zorder=4)

    ax.set_xlim(0, W - 1)
    ax.set_ylim(H - 1, 0)
    ax.set_title(
        f"{zsum_nc.parent.name}/{zsum_nc.stem}  ·  pitch={PITCH_UM:.3f} µm/px\n"
        f"wedge: apex→L/R cone, opening ≈ {opening_deg:.1f}°  ·  "
        f"covers {wedge_frac:.1f}% of crop  ·  "
        f"display = {mito_label}",
        fontsize=10)
    ax.legend(loc="upper right", fontsize=9, framealpha=0.9)
    ax.axis("off")

    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=140, bbox_inches="tight")
    plt.close(fig)
    print(f"Wrote {out_path}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "--zsum-nc",
        default=str(REPO / "projections" / "250612_patterned_plate_3" /
                    "B04_TRAK2_250616" / "denoised" /
                    "cell1 - Denoised.nc"),
        help="Path to the multi-channel cropped Z-sum projection .nc.")
    ap.add_argument(
        "--mito-nc", default=None,
        help="Optional bg-subtracted 488 MaxIP .nc. Inferred from --zsum-nc "
             "if omitted (suffix `_488_bg_subtracted.nc`).")
    ap.add_argument(
        "--out",
        default=str(REPO / "analysis" / "figures_wedge_r_ks" /
                    "wedge_illustration_offline.png"))
    args = ap.parse_args()

    zsum_nc = pathlib.Path(args.zsum_nc).resolve()
    mito_nc = (pathlib.Path(args.mito_nc).resolve() if args.mito_nc else
               zsum_nc.with_name(zsum_nc.stem + "_488_bg_subtracted.nc"))
    render(zsum_nc, mito_nc, pathlib.Path(args.out).resolve())


if __name__ == "__main__":
    main()
