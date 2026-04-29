"""Render the canonical wedge geometry on a real cell, for the manuscript
reviewer document.

Loads one ND2, runs the same template-matching + crop logic the pipeline
uses, then overlays:

  * the orange pattern contour,
  * the segmented nucleus,
  * the apex / L / R anchor points (yellow / magenta / red),
  * the two wedge bounding rays from apex through L and R,
  * a translucent red overlay on the cropped image *outside* the wedge,

so the reviewer can see at a glance which pixels contribute to the
wedge-r KS metric. Anchor coordinates and the wedge mask come from the
constants and helpers in `template_matching_bulk` itself, so the
illustration is guaranteed consistent with the metric the pipeline
emits.

CLI:
    pixi run python replication/plot_wedge_illustration.py \\
        --cell mark_data/patterned_data/250612_patterned_plate_3/B04_TRAK2_250616/cell1.nd2 \\
        --out  replication/figures_wedge_r_ks/wedge_illustration.png
"""
from __future__ import annotations

import argparse
import pathlib
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import nd2
import numpy as np
import skimage
from scipy.ndimage import center_of_mass

REPO = pathlib.Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))

import template_matching_bulk as tmb  # noqa: E402


def render(cell_path: pathlib.Path, out_path: pathlib.Path):
    template_hat = tmb.get_template_hat(1326)
    template = tmb.get_padded_template_at_width(1326)
    key = tmb.cluster_key(cell_path)
    offset = tmb.offset_overrides.get(key, [128, 128])
    roi = tmb.roi_overrides.get(key, None)

    img = nd2.imread(cell_path, xarray=True)
    zsum = img.sum(axis=0)

    max_coords = tmb.get_template_center(img, cell_path,
                                         template_hat=template_hat,
                                         offset=offset, roi=roi)
    shifted_template = np.roll(template,
                               (max_coords[0] - 1024, max_coords[1] - 1024),
                               axis=(0, 1))

    y_start = max_coords[0] - 512 + offset[0]
    y_end = max_coords[0] + 512 + offset[0]
    x_start = max_coords[1] - 512 + offset[1]
    x_end = max_coords[1] + 512 + offset[1]
    cropped_zsum = zsum.isel(Y=slice(y_start, y_end), X=slice(x_start, x_end))
    nuc = cropped_zsum.sel(C="405").to_numpy()
    # Display channel: MaxIP of mitochondria (matches what the metric is
    # actually computed on).
    cropped_max = img.isel(Y=slice(y_start, y_end), X=slice(x_start, x_end)).max(axis=0)
    mito = cropped_max.sel(C="488").to_numpy().astype(np.float64)

    # Pattern outline in crop coords.
    contour = skimage.measure.find_contours(shifted_template)[0].copy()
    contour[:, 0] -= max_coords[0] - 512
    contour[:, 1] -= max_coords[1] - 512

    # Nucleus segmentation (same recipe as the pipeline).
    nuc_mask = nuc > skimage.filters.threshold_otsu(nuc)
    nl = skimage.measure.label(nuc_mask)
    nprops = skimage.measure.regionprops(nl)
    nmax = int(np.argmax([p.area for p in nprops])) + 1
    nuc_mask = nl == nmax

    pitch = img.metadata["metadata"].channels[0].volume.axesCalibration[0]
    wedge_mask, _, _, _, vol_arc = tmb._get_wedge_geometry(mito.shape, pitch)
    apex = tmb.WEDGE_APEX
    left = tmb.WEDGE_LEFT
    right = tmb.WEDGE_RIGHT
    a_left = float(np.arctan2(left[0] - apex[0], left[1] - apex[1]))
    a_right = float(np.arctan2(right[0] - apex[0], right[1] - apex[1]))
    opening_deg = abs(np.degrees(a_left - a_right))
    wedge_frac = 100 * wedge_mask.sum() / wedge_mask.size

    fig, ax = plt.subplots(figsize=(7.5, 7.5))
    ax.imshow(mito, cmap="Greens_r", vmax=np.percentile(mito, 99.5))
    # Translucent red outside the wedge.
    outside = ~wedge_mask
    overlay = np.zeros((*wedge_mask.shape, 4))
    overlay[outside] = [1, 0, 0, 0.32]
    ax.imshow(overlay)
    # Pattern outline.
    ax.plot(contour[:, 1], contour[:, 0], color="orange", linewidth=1.4,
            alpha=0.95, label="pattern outline")
    # Nucleus.
    for c in skimage.measure.find_contours(nuc_mask.astype(float), 0.5):
        ax.plot(c[:, 1], c[:, 0], color="cyan", linewidth=1.4)
    # Anchors.
    ax.scatter([apex[1]], [apex[0]], color="red", s=140, marker="v",
               edgecolor="black", linewidth=1, zorder=5,
               label=f"apex {apex}")
    ax.scatter([left[1]], [left[0]], color="yellow", s=110, marker="<",
               edgecolor="black", linewidth=1, zorder=5,
               label=f"L {left}")
    ax.scatter([right[1]], [right[0]], color="magenta", s=110, marker=">",
               edgecolor="black", linewidth=1, zorder=5,
               label=f"R {right}")
    # Bounding rays from apex through L and R.
    H, W = mito.shape
    for name, pt, color in (("left", left, "yellow"),
                            ("right", right, "magenta")):
        dy = pt[0] - apex[0]
        dx = pt[1] - apex[1]
        ey, ex = apex[0] + dy * 6.0, apex[1] + dx * 6.0
        ax.plot([apex[1], ex], [apex[0], ey], color=color, linewidth=1.8,
                linestyle="--", alpha=0.95, zorder=4)
    ax.set_xlim(0, W - 1)
    ax.set_ylim(H - 1, 0)
    ax.set_title(
        f"{cell_path.name}  ·  pitch={pitch:.3f} µm/px\n"
        f"opening ≈ {opening_deg:.1f}°  ·  wedge covers "
        f"{wedge_frac:.1f}% of crop",
        fontsize=10)
    ax.legend(loc="upper right", fontsize=8, framealpha=0.9)
    ax.axis("off")

    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=140, bbox_inches="tight")
    plt.close(fig)
    print(f"Wrote {out_path}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "--cell",
        default=str(REPO / "mark_data" / "patterned_data" /
                    "250612_patterned_plate_3" / "B04_TRAK2_250616" /
                    "cell1.nd2"),
        help="ND2 file to render the wedge over.")
    ap.add_argument(
        "--out",
        default=str(REPO / "replication" / "wedge_illustration.png"),
        help="Output PNG path.")
    args = ap.parse_args()
    render(pathlib.Path(args.cell).resolve(), pathlib.Path(args.out).resolve())


if __name__ == "__main__":
    main()
