"""Geometry verification: where does the wedge apex sit relative to the
cropped image midpoint?

The pipeline crops every cell into a 1024×1024 slice with the template-matched
pattern center at slice position (512, 512) — i.e. the image midpoint. The
wedge is then anchored at the *pattern bottom extremum* (the stalk tip),
which sits at a fixed offset from the pattern center determined by the rigid
template shape.

This script:
  1. Computes the pattern-bottom offset analytically from the canonical
     template (deterministic — same for every cell).
  2. Renders one representative cell with both anchors overlaid:
       - red: current wedge apex (pattern bottom) at slice (896, 512)
       - yellow: image-midpoint apex at slice (512, 512)
  3. Shows the two wedge cones side by side so the geometric difference is
     immediately visible.
"""
from __future__ import annotations
import os, sys, pathlib, numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import nd2, skimage, polars as pl

REPO = pathlib.Path("/Users/gladkoc/Dev/micropattern_cell_analysis")
sys.path.insert(0, str(REPO))
os.chdir(REPO)
os.environ.setdefault("MICROPATTERN_DATA_ROOT",
                     "/Volumes/valelab/_for_Mark/patterned_data")
import template_matching_bulk as tmb  # noqa: E402

OUT_DIR = REPO / "replication" / "overnight_final_out" / "geometry_verification"
OUT_DIR.mkdir(parents=True, exist_ok=True)


def template_extremes():
    template = tmb.get_padded_template_at_width(1326)
    ys, xs = np.where(template > 0)
    # Replicate _pattern_extremes in the 1024-slice frame (slice spans
    # template[1024-512:1024+512]) so we get the actual pixel chosen as
    # 'left/right/top/bottom' — these are single-pixel argmin/argmax picks,
    # not centroids.
    slice_ys = ys - (1024 - 512)
    slice_xs = xs - (1024 - 512)
    mask = (slice_ys >= 0) & (slice_ys < 1024) & (slice_xs >= 0) & (slice_xs < 1024)
    sy = slice_ys[mask]; sx = slice_xs[mask]
    return {
        "bottom_y_offset": int(ys.max() - 1024),
        "top_y_offset":    int(ys.min() - 1024),
        "left_x_offset":   int(xs.min() - 1024),
        "right_x_offset":  int(xs.max() - 1024),
        "x_center_offset": int((xs.max() + xs.min()) // 2 - 1024),
        # Slice-frame extremes (the actual pixels the pipeline picks)
        "ext_bottom": (int(sy[np.argmax(sy)]), int(sx[np.argmax(sy)])),
        "ext_top":    (int(sy[np.argmin(sy)]), int(sx[np.argmin(sy)])),
        "ext_left":   (int(sy[np.argmin(sx)]), int(sx[np.argmin(sx)])),
        "ext_right":  (int(sy[np.argmax(sx)]), int(sx[np.argmax(sx)])),
        "template": template,
    }


def build_wedge_mask(apex, left, right, shape=(1024, 1024)):
    """Replicate `_build_wedge_geometry`'s wedge mask given any apex."""
    H, W = shape
    Y, X = np.mgrid[:H, :W]
    ang = np.arctan2(Y - apex[0], X - apex[1])
    a_left = float(np.arctan2(left[0] - apex[0], left[1] - apex[1]))
    a_right = float(np.arctan2(right[0] - apex[0], right[1] - apex[1]))
    a_up = -np.pi / 2
    lo, hi = min(a_left, a_right), max(a_left, a_right)
    if lo <= a_up <= hi:
        wedge = (ang >= lo) & (ang <= hi)
        opening = hi - lo
    else:
        wedge = (ang <= lo) | (ang >= hi)
        opening = (2 * np.pi) - (hi - lo)
    return wedge, np.degrees(opening)


def find_arch_tangent_points(template, apex, arch_y_max):
    """Tangent points from `apex` to the pattern's outer boundary in the
    portion above the apex (pattern pixels with y ≤ arch_y_max). The
    tangent is the most-extreme-angle pixel within that subset, so a
    wedge with rays through these points encloses the entire upper portion
    of the pattern.

    To exclude the narrow stalk near the apex (which would force
    near-horizontal rays), arch_y_max is set well above the apex.
    Using arch_y_max = apex_y - stalk_buffer (default ~80 px ≈ 5 µm above
    apex) keeps the search to the wide arch+lower-arm region and excludes
    the stalk.

    Returns (left_tangent, right_tangent) as (y, x) tuples in slice coords.
    """
    ys, xs = np.where(template > 0)
    sy = ys - (1024 - 512)
    sx = xs - (1024 - 512)
    in_slice = (sy >= 0) & (sy < 1024) & (sx >= 0) & (sx < 1024)
    sy, sx = sy[in_slice], sx[in_slice]
    arch = sy <= arch_y_max
    sy_a, sx_a = sy[arch], sx[arch]
    ang = np.arctan2(sy_a - apex[0], sx_a - apex[1])
    i_left = int(np.argmin(ang))
    i_right = int(np.argmax(ang))
    return ((int(sy_a[i_left]), int(sx_a[i_left])),
            (int(sy_a[i_right]), int(sx_a[i_right])))


def find_stalk_top_y(template):
    """Find the y row at the bottom of the donut/arms (top of the stalk).
    The pattern's per-row x-extent jumps from ~115 px (stalk) to ~470 px
    (donut) between y=400 and y=327. We pick the first y (walking up from
    apex) at which extent exceeds 200 px — that's a strong enough threshold
    to skip the stalk entirely and land at the lowest donut-leg row."""
    ys, xs = np.where(template > 0)
    sy = ys - (1024 - 512)
    sx = xs - (1024 - 512)
    in_slice = (sy >= 0) & (sy < 1024) & (sx >= 0) & (sx < 1024)
    sy, sx = sy[in_slice], sx[in_slice]
    by_row = {}
    for y, x in zip(sy, sx):
        rng = by_row.setdefault(int(y), [10**9, -10**9])
        rng[0] = min(rng[0], int(x))
        rng[1] = max(rng[1], int(x))
    for y in sorted(by_row.keys(), reverse=True):
        lo, hi = by_row[y]
        if hi - lo > 200:
            return y
    return min(by_row.keys())


def render_overlay(cell_path: str, label: str, te: dict, pitch_um: float):
    """Render a figure with the current wedge and a tangent-based
    alternative overlaid using semi-transparent fills."""
    template = te["template"]
    template_hat = tmb.get_template_hat(1326)
    img = nd2.imread(cell_path, xarray=True)
    key = tmb.cluster_key(cell_path)
    offset = tmb.offset_overrides.get(key, [128, 128])
    roi = tmb.roi_overrides.get(key, None)
    max_coords = tmb.get_template_center(img, cell_path,
                                         template_hat=template_hat,
                                         offset=offset, roi=roi)
    zsum = img.sum(axis=0).sel(C="488").to_numpy().astype(np.float64)
    y0 = max_coords[0] - 512 + offset[0]
    x0 = max_coords[1] - 512 + offset[1]
    crop = zsum[y0:y0 + 1024, x0:x0 + 1024]

    # CURRENT wedge: apex at pattern_bottom pixel, rays through leftmost/rightmost-x
    apex_cur = te["ext_bottom"]                            # (896, 499)
    L_cur = te["ext_left"]                                 # leftmost-x pixel
    R_cur = te["ext_right"]                                # rightmost-x pixel

    # ALTERNATIVE wedge: apex at lateral midpoint, 25 µm below vertical midpoint;
    # tangent points cover the full pattern outline including the lower arms,
    # excluding only the narrow stalk near the apex (which would force
    # near-horizontal rays). The stalk-cutoff y is determined automatically
    # from where the pattern's x-extent first exceeds ~5 µm.
    apex_alt = (round(511.5 + 25.0 / pitch_um), 512)       # (896, 512)
    stalk_top_y = find_stalk_top_y(template)
    L_tan, R_tan = find_arch_tangent_points(template, apex_alt, stalk_top_y)

    mask_cur, open_cur = build_wedge_mask(apex_cur, L_cur, R_cur)
    mask_alt, open_alt = build_wedge_mask(apex_alt, L_tan, R_tan)

    SC = 1.4
    fig, ax = plt.subplots(figsize=(11, 11))
    v0, v1 = np.percentile(crop, (1, 99))
    ax.imshow(crop, cmap="gray", vmin=v0, vmax=v1)

    # Pattern outline
    shifted = np.roll(template, (max_coords[0] - 1024, max_coords[1] - 1024),
                      axis=(0, 1))
    pattern_mask = shifted[max_coords[0] - 512:max_coords[0] + 512,
                           max_coords[1] - 512:max_coords[1] + 512] > 0
    contours = skimage.measure.find_contours(pattern_mask.astype(float), 0.5)
    for c in contours:
        ax.plot(c[:, 1], c[:, 0], "w-", linewidth=0.8, alpha=0.65)

    # Wedge fills (semi-transparent overlays)
    rgba_cur = np.zeros((1024, 1024, 4))
    rgba_cur[mask_cur] = [1.0, 0.20, 0.20, 0.30]   # red, alpha 0.30
    ax.imshow(rgba_cur)
    rgba_alt = np.zeros((1024, 1024, 4))
    rgba_alt[mask_alt] = [0.10, 0.55, 1.0, 0.30]   # blue, alpha 0.30
    ax.imshow(rgba_alt)

    # Apex markers (different for each wedge)
    ax.plot(apex_cur[1], apex_cur[0], marker="*", color="red",
            markersize=22, mew=1.5, markeredgecolor="black",
            label=f"current apex {apex_cur}  (pattern_bottom pixel)")
    ax.plot(apex_alt[1], apex_alt[0], marker="*", color="dodgerblue",
            markersize=22, mew=1.5, markeredgecolor="black",
            label=f"alternative apex {apex_alt}  (lateral midpoint, 25 µm below)")
    # Boundary points
    ax.plot(L_cur[1], L_cur[0], "o", color="red", markersize=9,
            markeredgecolor="black",
            label=f"current L,R: {L_cur}, {R_cur}  (leftmost/rightmost-x pixels)")
    ax.plot(R_cur[1], R_cur[0], "o", color="red", markersize=9,
            markeredgecolor="black")
    ax.plot(L_tan[1], L_tan[0], "s", color="dodgerblue", markersize=9,
            markeredgecolor="black",
            label=f"tangent L,R: {L_tan}, {R_tan}  (upper-arch tangents from alt apex)")
    ax.plot(R_tan[1], R_tan[0], "s", color="dodgerblue", markersize=9,
            markeredgecolor="black")
    ax.plot(512, 512, "+", color="yellow", markersize=18, mew=2.0,
            label="image midpoint (512, 512)")

    ax.set_xlim(0, 1023); ax.set_ylim(1023, 0)
    ax.set_xticks([0, 256, 512, 768, 1023])
    ax.set_yticks([0, 256, 512, 768, 1023])
    ax.grid(alpha=0.25)
    ax.legend(loc="upper right", fontsize=9 * SC)

    extra = mask_alt & ~mask_cur                # area only in alt wedge
    lost = mask_cur & ~mask_alt                  # area only in current wedge
    pct_gain = 100 * extra.sum() / mask_cur.sum()
    pct_lost = 100 * lost.sum() / mask_cur.sum()
    ax.set_title(
        f"Wedge overlay — {label}\n"
        f"Red fill: current wedge (apex at pattern_bottom; rays through leftmost/rightmost-x pixels)\n"
        f"Blue fill: alt wedge (apex at midpoint+25µm-down; rays tangent to FULL ARCH+ARMS, stalk excluded)\n"
        f"opening: current={open_cur:.2f}°  alt={open_alt:.2f}°  (Δ={open_alt-open_cur:+.2f}°)   "
        f"area gain: +{pct_gain:.1f}%   area dropped: −{pct_lost:.1f}%",
        fontsize=10 * SC)
    out = OUT_DIR / f"{label}_wedge_overlay.png"
    plt.tight_layout()
    fig.savefig(out, dpi=130)
    plt.close(fig)
    print(f"  wrote {out.name}")
    return {"apex_cur": apex_cur, "apex_alt": apex_alt,
            "L_cur": L_cur, "R_cur": R_cur,
            "L_tan": L_tan, "R_tan": R_tan,
            "open_cur": open_cur, "open_alt": open_alt,
            "area_gain_pct": float(pct_gain),
            "area_lost_pct": float(pct_lost)}


def render_cell(cell_path: str, label: str, te: dict, pitch_um: float):
    template = te["template"]
    template_hat = tmb.get_template_hat(1326)
    img = nd2.imread(cell_path, xarray=True)

    key = tmb.cluster_key(cell_path)
    offset = tmb.offset_overrides.get(key, [128, 128])
    roi = tmb.roi_overrides.get(key, None)
    max_coords = tmb.get_template_center(img, cell_path,
                                         template_hat=template_hat,
                                         offset=offset, roi=roi)

    # Crop the same way the pipeline does
    zsum = img.sum(axis=0).sel(C="488").to_numpy().astype(np.float64)
    y0 = max_coords[0] - 512 + offset[0]
    x0 = max_coords[1] - 512 + offset[1]
    crop = zsum[y0:y0 + 1024, x0:x0 + 1024]

    # Pattern center in slice frame is always (512, 512); pattern bottom
    # (current wedge apex) is at (512 + bottom_y_offset, 512 + x_center_offset).
    pat_cy, pat_cx = 512, 512
    apex_cur = (pat_cy + te["bottom_y_offset"], pat_cx + te["x_center_offset"])
    apex_alt = (512, 512)  # exact image midpoint

    # Build a contour of the pattern outline for overlay
    shifted = np.roll(template, (max_coords[0] - 1024, max_coords[1] - 1024),
                      axis=(0, 1))
    pattern_mask = shifted[max_coords[0] - 512:max_coords[0] + 512,
                           max_coords[1] - 512:max_coords[1] + 512] > 0
    contours = skimage.measure.find_contours(pattern_mask.astype(float), 0.5)

    # Wedge geometry (rays through pattern_left and pattern_right extremes)
    L = (pat_cy + 0, pat_cx + te["left_x_offset"])   # leftmost pattern point
    R = (pat_cy + 0, pat_cx + te["right_x_offset"])  # rightmost
    # The actual pipeline uses the topmost/leftmost/rightmost pattern pixels
    # at their actual y-row; we'll use horizontal extremes at y=512 for
    # illustration. Wedge rays go from apex through L and R.

    def wedge_lines(apex):
        """Compute end points where the rays from apex through L and R hit
        the slice boundaries (for drawing)."""
        ay, ax = apex
        out = []
        for tgt in (L, R):
            ty, tx = tgt
            dy, dx = ty - ay, tx - ax
            if dy == 0 and dx == 0:
                continue
            # Extend ray until it leaves the [0, 1023] box
            t_vals = []
            if dy != 0:
                t_vals += [(0 - ay) / dy, (1023 - ay) / dy]
            if dx != 0:
                t_vals += [(0 - ax) / dx, (1023 - ax) / dx]
            t_vals = [t for t in t_vals if t > 0]
            if not t_vals:
                continue
            t = min(t_vals)
            out.append(((ay, ax), (ay + t * dy, ax + t * dx)))
        return out

    SC = 1.4
    fig, axes = plt.subplots(1, 2, figsize=(15, 7.5))
    for ax, apex, color, name in [
        (axes[0], apex_cur, "red", "current: wedge apex at pattern-bottom"),
        (axes[1], apex_alt, "yellow", "alternative: wedge apex at image midpoint"),
    ]:
        v0, v1 = np.percentile(crop, (1, 99))
        ax.imshow(crop, cmap="gray", vmin=v0, vmax=v1)
        # Pattern outline
        for c in contours:
            ax.plot(c[:, 1], c[:, 0], "w-", linewidth=0.6, alpha=0.7)
        # Image midpoint cross
        ax.plot(512, 512, "+", color="yellow", markersize=18, mew=2,
                label="image midpoint (512, 512)")
        # Wedge apex
        ax.plot(apex[1], apex[0], "*", color=color, markersize=22, mew=1.5,
                markeredgecolor="black", label=f"wedge apex {apex}")
        # Wedge rays
        for (a, b) in wedge_lines(apex):
            ax.plot([a[1], b[1]], [a[0], b[0]], color=color, linewidth=1.6,
                    linestyle="--", alpha=0.85)
        ax.set_title(name, fontsize=11 * SC)
        ax.legend(loc="upper right", fontsize=8 * SC)
        ax.set_xlim(0, 1023); ax.set_ylim(1023, 0)
        ax.set_xticks([0, 256, 512, 768, 1023])
        ax.set_yticks([0, 256, 512, 768, 1023])
        ax.grid(alpha=0.25)

    dy_um = te["bottom_y_offset"] * pitch_um
    dx_um = te["x_center_offset"] * pitch_um
    fig.suptitle(
        f"Geometry verification — {label}\n"
        f"Apex displacement (current vs midpoint): "
        f"Δy = {te['bottom_y_offset']:+d} px ({dy_um:+.2f} µm),  "
        f"Δx = {te['x_center_offset']:+d} px ({dx_um:+.2f} µm),  "
        f"|Δ| = {np.hypot(te['bottom_y_offset'], te['x_center_offset'])*pitch_um:.2f} µm",
        fontsize=11 * SC)
    plt.tight_layout(rect=[0, 0.02, 1, 0.92])
    out = OUT_DIR / f"{label}_geometry_verification.png"
    fig.savefig(out, dpi=130)
    plt.close(fig)
    print(f"  wrote {out.name}")


def main():
    te = template_extremes()
    print("Pattern extremes in template frame (offsets from template center):")
    print(f"  y_min (top):    {te['top_y_offset']:+d} px")
    print(f"  y_max (bottom): {te['bottom_y_offset']:+d} px")
    print(f"  x_min (left):   {te['left_x_offset']:+d} px")
    print(f"  x_max (right):  {te['right_x_offset']:+d} px")
    print(f"  x_center bias:  {te['x_center_offset']:+d} px")

    df = pl.read_csv("replication/overnight_final_out/combined_raw.csv")
    pitch = float(df["lateral_pixel_pitch_um"].mean())
    print(f"\nMean pixel pitch: {pitch:.5f} µm/px")
    print(f"Wedge-apex displacement from image midpoint:")
    print(f"  Δy = {te['bottom_y_offset']:+d} px = {te['bottom_y_offset']*pitch:+.2f} µm  (DOWN)")
    print(f"  Δx = {te['x_center_offset']:+d} px = {te['x_center_offset']*pitch:+.2f} µm  (essentially 0)")
    print(f"  |Δ| = {np.hypot(te['bottom_y_offset'], te['x_center_offset'])*pitch:.2f} µm")
    print()
    print("This displacement is DETERMINISTIC: it depends only on the rigid template,")
    print("so it is identical for every one of the 494 cells in the dataset.")

    # Render one representative cell each from no TRAK / TRAK1 / TRAK2 (mito)
    targets = [
        ("noTRAK_plate1_F03",
         "/Volumes/valelab/_for_Mark/patterned_data/250521_patterned_plate_1/F03_no_TRAK_combined/Cell1.nd2"),
    ]
    # Pick one that we know exists
    samples = df.sort("path").head(1)
    cell_path = samples["path"][0]
    label = pathlib.Path(cell_path).parent.parent.name + "__" + \
            pathlib.Path(cell_path).parent.name + "__" + \
            pathlib.Path(cell_path).stem
    print(f"\nRendering one representative cell: {label}")
    render_cell(cell_path, label, te, pitch)

    # New overlay figure: midpoint-apex wedge with upper-arch tangent
    print(f"\nRendering wedge-overlay figure (upper-arch tangent alternative)...")
    info = render_overlay(cell_path, label, te, pitch)
    print(f"  current   apex / L / R:   {info['apex_cur']} / {info['L_cur']} / {info['R_cur']}")
    print(f"  alt       apex / L / R:   {info['apex_alt']} / {info['L_tan']} / {info['R_tan']}")
    print(f"  opening:                  cur={info['open_cur']:.3f}°, "
          f"alt={info['open_alt']:.3f}°  "
          f"(Δ={info['open_alt']-info['open_cur']:+.3f}°)")
    print(f"  alt wedge gains:          +{info['area_gain_pct']:.2f}% area  "
          f"(area only in current: −{info['area_lost_pct']:.2f}%)")


if __name__ == "__main__":
    main()
