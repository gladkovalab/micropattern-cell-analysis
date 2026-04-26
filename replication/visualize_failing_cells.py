"""Diagnose the 5 cells that failed in the whole-dataset run.

For each failing cell, render a 2-panel figure:
  - left:  full-frame DAPI z-sum with the template-match center marked,
           pattern outline overlaid, and override coordinate (if any) shown
  - right: cropped DAPI in the matched frame, with attempted nucleus seg

The ND2 reads use the same code path the pipeline uses, so the visual
mirrors what the pipeline 'sees' before the failure.
"""
from __future__ import annotations
import os, pathlib, sys, traceback
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import nd2, skimage

REPO = pathlib.Path("/Users/gladkoc/Dev/micropattern_cell_analysis")
sys.path.insert(0, str(REPO))
os.chdir(REPO)
os.environ.setdefault("MICROPATTERN_DATA_ROOT", "/Volumes/valelab/_for_Mark/patterned_data")

import template_matching_bulk as tmb  # noqa: E402

OUT_DIR = REPO / "replication" / "overnight_final_out" / "failing_cells_diagnostic"
OUT_DIR.mkdir(parents=True, exist_ok=True)

FAILING = [
    "/Volumes/valelab/_for_Mark/patterned_data/250731_patterned_plate_11_good/B03_250806_ctrl_siRNA_Ars/Cell1.nd2",
    "/Volumes/valelab/_for_Mark/patterned_data/250731_patterned_plate_11_good/E03_250806_ctrl_siRNA_ctrl/Cell12.nd2",
    "/Volumes/valelab/_for_Mark/patterned_data/250731_patterned_plate_11_good/E03_250806_ctrl_siRNA_ctrl/Cell2.nd2",
    "/Volumes/valelab/_for_Mark/patterned_data/250731_patterned_plate_11_good/F05_250808_TRAK2_wt/Cell12.nd2",
    "/Volumes/valelab/_for_Mark/patterned_data/250807_patterned_plate_12/G04_250814_MAPK9_siRNA_ctrl/Cell6.nd2",
]


def has_override(img_path: str) -> tuple[bool, tuple | None]:
    """Try several key conventions to locate an override for this raw ND2."""
    p = pathlib.Path(img_path)
    candidates = [
        tmb.cluster_key(img_path),  # raw — the buggy lookup
        tmb.cluster_key(str(p.parent / "denoised" / f"{p.stem} - Denoised.nd2")),
    ]
    for k in candidates:
        if k in tmb.coordinate_overrides_dict:
            return True, tmb.coordinate_overrides_dict[k]
    return False, None


def visualize(img_path: str):
    p = pathlib.Path(img_path)
    plate = p.parent.parent.name
    well = p.parent.name
    label = f"{plate}__{well}__{p.stem}"
    out_png = OUT_DIR / f"{label}.png"

    template_hat = tmb.get_template_hat(1326)
    template = tmb.get_padded_template_at_width(1326)

    img = nd2.imread(img_path, xarray=True)
    zsum = img.sum(axis=0).sel(C="405").to_numpy().astype(np.float64)

    # Run template matching to get the auto-detected pattern center
    key = tmb.cluster_key(img_path)
    offset = tmb.offset_overrides.get(key, [128, 128])
    roi = tmb.roi_overrides.get(key, None)
    try:
        max_coords = tmb.get_template_center(img, img_path,
                                             template_hat=template_hat,
                                             offset=offset, roi=roi)
        match_ok = True
    except Exception:
        max_coords = None
        match_ok = False

    has_ov, ov_coord = has_override(img_path)

    fig, axes = plt.subplots(1, 2, figsize=(14, 7))

    # Left: full-frame DAPI z-sum with annotations
    ax = axes[0]
    vmin, vmax = np.percentile(zsum, (1, 99.5))
    ax.imshow(zsum, cmap="gray", vmin=vmin, vmax=vmax)
    title_bits = [f"{plate}/{well}/{p.stem}", f"DAPI z-sum, full frame {zsum.shape}"]
    if match_ok:
        ax.plot(max_coords[1], max_coords[0], "r+", markersize=20, mew=2,
                label=f"auto match @ ({max_coords[0]}, {max_coords[1]})")
        try:
            shifted_template = np.roll(template,
                                       (max_coords[0] - 1024, max_coords[1] - 1024),
                                       axis=(0, 1))
            contours = skimage.measure.find_contours(shifted_template)
            if contours:
                contour = contours[0]
                ax.plot(contour[:, 1], contour[:, 0], "r-", linewidth=1.0,
                        alpha=0.7, label="auto template outline")
        except Exception:
            pass
    else:
        title_bits.append("(template match raised)")
    if has_ov:
        ov_x, ov_y = ov_coord
        ax.plot(ov_x, ov_y, "y*", markersize=22, mew=2,
                markerfacecolor="yellow", markeredgecolor="black",
                label=f"override @ ({ov_y}, {ov_x})")
        try:
            shifted_template_ov = np.roll(template,
                                          (ov_y - 1024, ov_x - 1024),
                                          axis=(0, 1))
            contours = skimage.measure.find_contours(shifted_template_ov)
            if contours:
                contour = contours[0]
                ax.plot(contour[:, 1], contour[:, 0], "y--", linewidth=1.2,
                        alpha=0.85, label="override template outline")
        except Exception:
            pass
    ax.legend(loc="upper right", fontsize=8)
    ax.set_title("\n".join(title_bits), fontsize=10)
    ax.set_axis_off()

    # Right: crop in the (auto or override) frame
    ax2 = axes[1]
    use_coord = ov_coord and (ov_coord[1], ov_coord[0]) or (max_coords if match_ok else None)
    if use_coord:
        cy, cx = (ov_coord[1], ov_coord[0]) if has_ov else (max_coords[0], max_coords[1])
        y0, y1 = cy - 512 + offset[0], cy + 512 + offset[0]
        x0, x1 = cx - 512 + offset[1], cx + 512 + offset[1]
        y0, y1 = max(0, y0), min(zsum.shape[0], y1)
        x0, x1 = max(0, x0), min(zsum.shape[1], x1)
        crop = zsum[y0:y1, x0:x1]
        if crop.size > 0:
            vmin2, vmax2 = np.percentile(crop, (1, 99.5))
            ax2.imshow(crop, cmap="gray", vmin=vmin2, vmax=vmax2)
            ax2.set_title(
                f"crop @ {'OVERRIDE' if has_ov else 'auto'} center  "
                f"shape={crop.shape}", fontsize=10)
        else:
            ax2.text(0.5, 0.5, "crop is empty\n(center near image edge)",
                     ha="center", va="center", transform=ax2.transAxes)
    else:
        ax2.text(0.5, 0.5, "no usable center\n(auto match failed, no override)",
                 ha="center", va="center", transform=ax2.transAxes)
    ax2.set_axis_off()

    suptitle = f"FAIL DIAGNOSTIC — {label}"
    if has_ov:
        suptitle += "  ★ OVERRIDE EXISTS (under denoised key, not picked up by raw run)"
    else:
        suptitle += "  ✗ no override found"
    fig.suptitle(suptitle, fontsize=11)
    plt.tight_layout(rect=[0, 0.02, 1, 0.95])
    fig.savefig(out_png, dpi=130)
    plt.close(fig)
    print(f"  wrote {out_png.name}  match_ok={match_ok}  override={has_ov}")
    return {"label": label, "match_ok": match_ok, "has_override": has_ov,
            "override_coord_xy": ov_coord, "auto_max_coords_yx": max_coords}


def main():
    print(f"Visualizing {len(FAILING)} failing cells → {OUT_DIR}")
    results = []
    for p in FAILING:
        try:
            results.append(visualize(p))
        except Exception:
            print(f"  ERR on {p}:")
            traceback.print_exc()
    print()
    print("Summary:")
    print(f"  {'label':<70} match_ok  override")
    for r in results:
        print(f"  {r['label']:<70} {str(r['match_ok']):<8} {r['has_override']}")


if __name__ == "__main__":
    main()
