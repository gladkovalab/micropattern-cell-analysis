"""Visualize the wedge-restricted radial metric's angular sector.

v3 — corrections:
  * Orange contour is the authoritative pattern outline (unchanged, as in
    Mark's original pipeline — find_contours(shifted_template) then shifted by
    -(max_coords - 512)). Pattern extreme POINTS are re-extracted from a
    correctly-aligned pattern mask so they sit on the orange contour.
  * Three origin variants in one figure:
      A. Nucleus CoM
      B. Pattern CoM
      C. Pattern bottom-extremum (red marker)
  * All wedges sweep UPWARD toward the arch, bounded by rays through the
    leftmost and rightmost pattern extreme points.

Cell used: plate 3 B04 TRAK2 cell1.
"""
from __future__ import annotations

import pathlib
import sys

import nd2
import numpy as np
import skimage
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.ndimage import center_of_mass

REPO = pathlib.Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))

import template_matching_bulk as tmb  # noqa: E402

CELL = pathlib.Path("/Volumes/valelab/_for_Mark/patterned_data/"
                    "250612_patterned_plate_3/B04_TRAK2_250616/cell1.nd2")
OUT = REPO / "replication" / "overnight_fig4b_out" / "figures" / "wedge_illustration_v3.png"


def pattern_extremes(pattern_mask):
    ys, xs = np.where(pattern_mask)
    return {
        "bottom": (int(ys[np.argmax(ys)]), int(xs[np.argmax(ys)])),
        "top":    (int(ys[np.argmin(ys)]), int(xs[np.argmin(ys)])),
        "left":   (int(ys[np.argmin(xs)]), int(xs[np.argmin(xs)])),
        "right":  (int(ys[np.argmax(xs)]), int(xs[np.argmax(xs)])),
    }


def upper_wedge_mask(shape, origin_yx, left_pt, right_pt):
    H, W = shape
    Y, X = np.mgrid[:H, :W]
    ang = np.arctan2(Y - origin_yx[0], X - origin_yx[1])

    def pt_angle(pt):
        return float(np.arctan2(pt[0] - origin_yx[0], pt[1] - origin_yx[1]))

    a_left = pt_angle(left_pt)
    a_right = pt_angle(right_pt)
    a_up = -np.pi / 2
    lo = min(a_left, a_right)
    hi = max(a_left, a_right)
    if lo <= a_up <= hi:
        mask = (ang >= lo) & (ang <= hi)
    else:
        mask = (ang <= lo) | (ang >= hi)
    return mask, (a_left, a_right)


def render_wedge(ax, mito, contour, nuc_mask, nuc_com, pattern_com, ext,
                 origin_yx, origin_label, title):
    vmax = np.percentile(mito, 99.5)
    ax.imshow(mito, cmap="Greens_r", vmax=vmax)

    # Build and overlay wedge
    wedge_mask_img, angs = upper_wedge_mask(mito.shape, origin_yx,
                                            ext["left"], ext["right"])
    outside = ~wedge_mask_img
    red_overlay = np.zeros((*wedge_mask_img.shape, 4))
    red_overlay[outside] = [1, 0, 0, 0.35]
    ax.imshow(red_overlay)

    # Pattern outline (authoritative)
    ax.plot(contour[:, 1], contour[:, 0], color="orange", linewidth=1.4, alpha=0.95)
    # Nucleus contour
    for c in skimage.measure.find_contours(nuc_mask.astype(float), 0.5):
        ax.plot(c[:, 1], c[:, 0], color="cyan", linewidth=1.4)
    # Nucleus CoM (blue +) and pattern CoM (white x)
    ax.scatter([nuc_com[1]], [nuc_com[0]], color="blue", s=80, marker="+",
               linewidth=2.5, zorder=6)
    ax.scatter([pattern_com[1]], [pattern_com[0]], color="white", s=60,
               marker="x", linewidth=2, zorder=6)
    # Pattern extreme points — on the orange contour
    for name, color, marker in [("left", "yellow", "<"), ("right", "magenta", ">"),
                                ("top", "white", "^"), ("bottom", "red", "v")]:
        pt = ext[name]
        ax.scatter([pt[1]], [pt[0]], color=color, s=110, marker=marker,
                   edgecolor="black", linewidth=1, zorder=5)

    # Highlight active origin
    ax.scatter([origin_yx[1]], [origin_yx[0]], s=260, facecolors="none",
               edgecolors="lime", linewidth=2.5, zorder=7)

    # Wedge boundary rays from active origin through left + right
    H, W = mito.shape
    for name, color in [("left", "yellow"), ("right", "magenta")]:
        pt = ext[name]
        dy = pt[0] - origin_yx[0]
        dx = pt[1] - origin_yx[1]
        ey, ex = origin_yx[0] + dy * 6.0, origin_yx[1] + dx * 6.0
        ax.plot([origin_yx[1], ex], [origin_yx[0], ey], color=color, linewidth=1.8,
                linestyle="--", alpha=0.95, zorder=4)

    opening = np.degrees(abs(angs[0] - angs[1]))
    frac = 100 * wedge_mask_img.sum() / wedge_mask_img.size
    ax.set_xlim(0, W - 1)
    ax.set_ylim(H - 1, 0)
    ax.set_title(f"{title}\norigin: {origin_label} · opening ≈ {opening:.0f}° · "
                 f"{frac:.1f}% of crop",
                 fontsize=9)
    ax.axis("off")


def main():
    template_hat = tmb.get_template_hat(1326)
    template = tmb.get_padded_template_at_width(1326)
    key = tmb.cluster_key(CELL)
    offset = tmb.offset_overrides.get(key, [128, 128])
    roi = tmb.roi_overrides.get(key, None)

    img = nd2.imread(CELL, xarray=True)
    zsum = img.sum(axis=0)

    max_coords = tmb.get_template_center(img, CELL, template_hat=template_hat,
                                         offset=offset, roi=roi)
    shifted_template = np.roll(template, (max_coords[0] - 1024, max_coords[1] - 1024),
                               axis=(0, 1))

    # Crop mito/nuclear image in the ORIGINAL frame
    y_start = max_coords[0] - 512 + offset[0]
    y_end = max_coords[0] + 512 + offset[0]
    x_start = max_coords[1] - 512 + offset[1]
    x_end = max_coords[1] + 512 + offset[1]
    cropped_zsum = zsum.isel(Y=slice(y_start, y_end), X=slice(x_start, x_end))
    nuc = cropped_zsum.sel(C="405").to_numpy()
    mito = cropped_zsum.sel(C="488").to_numpy()

    # ORANGE CONTOUR — authoritative pattern outline in the crop's frame
    # (unchanged from Mark's original logic)
    contour = skimage.measure.find_contours(shifted_template)[0].copy()
    contour[:, 0] -= max_coords[0] - 512
    contour[:, 1] -= max_coords[1] - 512

    # CORRECTLY-ALIGNED pattern mask — slice shifted_template in its own frame
    # (pattern center is at max_coords in shifted_template).
    pattern_mask = shifted_template[
        max_coords[0] - 512:max_coords[0] + 512,
        max_coords[1] - 512:max_coords[1] + 512
    ] > 0
    pattern_com = center_of_mass(pattern_mask)
    ext = pattern_extremes(pattern_mask)

    # Nucleus segmentation in the crop frame (as usual)
    nuc_mask = nuc > skimage.filters.threshold_otsu(nuc)
    nl = skimage.measure.label(nuc_mask)
    nprops = skimage.measure.regionprops(nl)
    nmax = int(np.argmax([p.area for p in nprops])) + 1
    nuc_mask = nl == nmax
    nuc_com = center_of_mass(nuc_mask)

    pitch = img.metadata["metadata"].channels[0].volume.axesCalibration[0]

    # Render three panels: nucleus CoM / pattern CoM / pattern bottom anchor
    fig, axes = plt.subplots(1, 3, figsize=(17, 6))
    render_wedge(axes[0], mito, contour, nuc_mask, nuc_com, pattern_com, ext,
                 nuc_com, "nucleus CoM",
                 "A. Nucleus CoM anchor")
    render_wedge(axes[1], mito, contour, nuc_mask, nuc_com, pattern_com, ext,
                 pattern_com, "pattern CoM",
                 "B. Pattern CoM anchor")
    render_wedge(axes[2], mito, contour, nuc_mask, nuc_com, pattern_com, ext,
                 ext["bottom"], "pattern bottom",
                 "C. Pattern bottom-extremum anchor")

    fig.suptitle(f"Wedge v3 · {CELL.name} (plate 3 TRAK2) · pitch={pitch:.3f} µm/px · "
                 f"red = excluded",
                 fontsize=10)
    plt.tight_layout(rect=[0, 0.02, 1, 0.96])
    OUT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT, dpi=130, bbox_inches="tight")
    plt.close(fig)
    print(f"Wrote {OUT}")
    for name in ("bottom", "left", "right", "top"):
        pt = ext[name]
        dy_um = (pt[0] - nuc_com[0]) * pitch
        dx_um = (pt[1] - nuc_com[1]) * pitch
        print(f"  pattern {name}: dy={dy_um:+.1f} µm, dx={dx_um:+.1f} µm from nucleus CoM")


if __name__ == "__main__":
    main()
