"""Visualize three candidate polar-coordinate schemes for the wedge-restricted
radial metric. All three share the SAME wedge (Panel C from wedge v3: anchored
at the pattern bottom-extremum, rays to left/right pattern extremes). They
differ in where the polar ORIGIN sits — the center of the concentric arcs
used to compute r for each pixel.

  Scheme 1: polar origin = wedge apex (red arrow, pattern bottom). Most
            coordinate-consistent option: wedge and radius share the same
            origin, so r and θ live in the same polar system.

  Scheme 2: polar origin = arch-top (white ^). Places r=0 at the periphery,
            r increasing toward the nucleus/stalk.

  Scheme 3: polar origin = nucleus CoM (blue +). Matches the existing Mark-
            style EDT radial, but restricted to the wedge.

Each panel shows the wedge boundary rays plus concentric circles at 5/10/15/
20/25/30/35/40 µm from that scheme's origin, clipped visually to the wedge.
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
from matplotlib.patches import Circle
from scipy.ndimage import center_of_mass

REPO = pathlib.Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))
import template_matching_bulk as tmb  # noqa: E402

CELL = pathlib.Path("/Volumes/valelab/_for_Mark/patterned_data/"
                    "250612_patterned_plate_3/B04_TRAK2_250616/cell1.nd2")
OUT = REPO / "replication" / "overnight_fig4b_out" / "figures" / "polar_schemes.png"

RING_RADII_UM = [5, 10, 15, 20, 25, 30, 35, 40]


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
    return mask


def render_scheme(ax, mito, contour, nuc_mask, ext, wedge_mask_img,
                  wedge_apex, polar_origin, pitch_um, title, origin_label):
    vmax = np.percentile(mito, 99.5)
    ax.imshow(mito, cmap="Greens_r", vmax=vmax)

    # Wedge overlay (dim red outside the wedge)
    outside = ~wedge_mask_img
    red_overlay = np.zeros((*wedge_mask_img.shape, 4))
    red_overlay[outside] = [1, 0, 0, 0.22]
    ax.imshow(red_overlay)

    # Pattern outline
    ax.plot(contour[:, 1], contour[:, 0], color="orange", linewidth=1.4, alpha=0.95)
    # Nucleus outline
    for c in skimage.measure.find_contours(nuc_mask.astype(float), 0.5):
        ax.plot(c[:, 1], c[:, 0], color="cyan", linewidth=1.4)

    # Concentric circles around polar origin (the sketch of the polar scheme)
    for r_um in RING_RADII_UM:
        r_px = r_um / pitch_um
        circ = Circle((polar_origin[1], polar_origin[0]), r_px,
                      fill=False, edgecolor="white", linewidth=0.9,
                      linestyle="-", alpha=0.85)
        ax.add_patch(circ)
        # Label the ring near the top of the image where it intersects
        lx = polar_origin[1]
        ly = polar_origin[0] - r_px
        if 0 <= ly < mito.shape[0] and 0 <= lx < mito.shape[1]:
            ax.text(lx + 3, ly - 3, f"{r_um}",
                    color="white", fontsize=6, ha="left", va="bottom",
                    alpha=0.9)

    # Wedge boundary rays (from wedge apex through left and right pattern extremes)
    H, W = mito.shape
    for name, color in [("left", "yellow"), ("right", "magenta")]:
        pt = ext[name]
        dy = pt[0] - wedge_apex[0]
        dx = pt[1] - wedge_apex[1]
        ey = wedge_apex[0] + dy * 8.0
        ex = wedge_apex[1] + dx * 8.0
        ax.plot([wedge_apex[1], ex], [wedge_apex[0], ey], color=color,
                linewidth=1.8, linestyle="--", alpha=0.95, zorder=4)

    # Pattern extreme points
    for name, color, marker in [("left", "yellow", "<"), ("right", "magenta", ">"),
                                ("top", "white", "^"), ("bottom", "red", "v")]:
        pt = ext[name]
        ax.scatter([pt[1]], [pt[0]], color=color, s=100, marker=marker,
                   edgecolor="black", linewidth=1, zorder=5)

    # Nucleus CoM and pattern CoM reference markers
    ax.scatter([center_of_mass(nuc_mask)[1]], [center_of_mass(nuc_mask)[0]],
               color="blue", s=70, marker="+", linewidth=2.2, zorder=6)

    # Highlight the polar origin
    ax.scatter([polar_origin[1]], [polar_origin[0]], s=260,
               facecolors="none", edgecolors="lime", linewidth=2.8, zorder=7)

    ax.set_xlim(0, W - 1)
    ax.set_ylim(H - 1, 0)
    ax.set_title(f"{title}\npolar origin: {origin_label}\n"
                 f"r = distance (µm) from origin; rings at 5, 10, …, 40 µm",
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
    y_start = max_coords[0] - 512 + offset[0]
    y_end = max_coords[0] + 512 + offset[0]
    x_start = max_coords[1] - 512 + offset[1]
    x_end = max_coords[1] + 512 + offset[1]
    cropped_zsum = zsum.isel(Y=slice(y_start, y_end), X=slice(x_start, x_end))
    nuc = cropped_zsum.sel(C="405").to_numpy()
    mito = cropped_zsum.sel(C="488").to_numpy()

    contour = skimage.measure.find_contours(shifted_template)[0].copy()
    contour[:, 0] -= max_coords[0] - 512
    contour[:, 1] -= max_coords[1] - 512

    pattern_mask = shifted_template[
        max_coords[0] - 512:max_coords[0] + 512,
        max_coords[1] - 512:max_coords[1] + 512
    ] > 0
    ext = pattern_extremes(pattern_mask)

    nuc_mask = nuc > skimage.filters.threshold_otsu(nuc)
    nl = skimage.measure.label(nuc_mask)
    nprops = skimage.measure.regionprops(nl)
    nmax = int(np.argmax([p.area for p in nprops])) + 1
    nuc_mask = nl == nmax
    nuc_com = center_of_mass(nuc_mask)

    pitch = img.metadata["metadata"].channels[0].volume.axesCalibration[0]

    # Wedge apex = pattern bottom extremum (Panel C)
    wedge_apex = ext["bottom"]
    wedge_mask_img = upper_wedge_mask(mito.shape, wedge_apex, ext["left"], ext["right"])

    fig, axes = plt.subplots(1, 3, figsize=(17, 6))
    render_scheme(axes[0], mito, contour, nuc_mask, ext, wedge_mask_img,
                  wedge_apex, wedge_apex, pitch,
                  "1. Polar origin = wedge apex",
                  "pattern bottom (red v)")
    render_scheme(axes[1], mito, contour, nuc_mask, ext, wedge_mask_img,
                  wedge_apex, ext["top"], pitch,
                  "2. Polar origin = arch-top",
                  "top of arch (white ^)")
    render_scheme(axes[2], mito, contour, nuc_mask, ext, wedge_mask_img,
                  wedge_apex, nuc_com, pitch,
                  "3. Polar origin = nucleus CoM",
                  "nucleus CoM (blue +)")

    fig.suptitle(f"Polar-coordinate schemes for wedge-restricted radial metric · "
                 f"{CELL.name} · pitch={pitch:.3f} µm/px",
                 fontsize=10)
    plt.tight_layout(rect=[0, 0.02, 1, 0.95])
    OUT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT, dpi=130, bbox_inches="tight")
    plt.close(fig)
    print(f"Wrote {OUT}")

    # Print where the wedge apex, arch-top, nucleus CoM sit in µm from each other
    def d_um(a, b):
        return np.hypot((a[0]-b[0])*pitch, (a[1]-b[1])*pitch)
    print(f"  wedge apex (red) → arch top (white): {d_um(wedge_apex, ext['top']):.1f} µm")
    print(f"  wedge apex (red) → nucleus CoM: {d_um(wedge_apex, nuc_com):.1f} µm")
    print(f"  arch top → nucleus CoM: {d_um(ext['top'], nuc_com):.1f} µm")


if __name__ == "__main__":
    main()
