"""For the 2 cells without overrides, compare template matching on the raw
ND2 vs the denoised ND2. The denoised version may give a usable center that
we can promote to a manual override on the raw run.
"""
from __future__ import annotations
import os, pathlib, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import nd2, skimage

REPO = pathlib.Path("/Users/gladkoc/Dev/micropattern_cell_analysis")
sys.path.insert(0, str(REPO))
os.chdir(REPO)
os.environ.setdefault("MICROPATTERN_DATA_ROOT",
                     "/Volumes/valelab/_for_Mark/patterned_data")

import template_matching_bulk as tmb  # noqa: E402

OUT_DIR = REPO / "replication" / "overnight_final_out" / "failing_cells_diagnostic"
OUT_DIR.mkdir(parents=True, exist_ok=True)

PAIRS = [
    {
        "label": "plate_11_F05_Cell12",
        "raw":      "/Volumes/valelab/_for_Mark/patterned_data/250731_patterned_plate_11_good/F05_250808_TRAK2_wt/Cell12.nd2",
        "denoised": "/Volumes/valelab/_for_Mark/patterned_data/250731_patterned_plate_11_good/F05_250808_TRAK2_wt/denoised/Cell12 - Denoised.nd2",
    },
    {
        "label": "plate_12_G04_Cell6",
        "raw":      "/Volumes/valelab/_for_Mark/patterned_data/250807_patterned_plate_12/G04_250814_MAPK9_siRNA_ctrl/Cell6.nd2",
        "denoised": "/Volumes/valelab/_for_Mark/patterned_data/250807_patterned_plate_12/G04_250814_MAPK9_siRNA_ctrl/denoised/Cell6 - Denoised.nd2",
    },
]


def match_with_score(img, img_path, template_hat, template):
    """Run template matching using the same FFT method as the pipeline.
    Returns (max_coords_yx, contrast=peak/abs_median, zsum_dapi)."""
    key = tmb.cluster_key(img_path)
    offset = tmb.offset_overrides.get(key, [128, 128])
    roi = tmb.roi_overrides.get(key, None)
    matching = tmb.match_template(img, template_hat=template_hat, offset=offset)
    if roi is not None:
        sub = matching[roi[0], roi[1]]
    else:
        sub = matching
    yx_local = np.unravel_index(np.argmax(sub), sub.shape)
    if roi is not None:
        oy = roi[0].start or 0
        ox = roi[1].start or 0
        yx = (yx_local[0] + oy, yx_local[1] + ox)
    else:
        yx = yx_local
    peak = float(sub.max())
    contrast = peak / max(float(np.abs(np.median(sub))), 1e-9)
    zsum = img.sum(axis=0).sel(C="405").to_numpy().astype(np.float64)
    return yx, contrast, zsum


def render_pair(pair, template_hat, template):
    fig, axes = plt.subplots(2, 2, figsize=(15, 14))
    rows = [("RAW", pair["raw"]), ("DENOISED", pair["denoised"])]
    summary = {"label": pair["label"]}
    for row_idx, (tag, path) in enumerate(rows):
        ax_full = axes[row_idx, 0]
        ax_crop = axes[row_idx, 1]
        try:
            img = nd2.imread(path, xarray=True)
            yx, peak, zsum = match_with_score(img, path, template_hat, template)
            summary[f"{tag}_yx"] = yx
            summary[f"{tag}_peak"] = peak
        except Exception as e:
            ax_full.text(0.5, 0.5, f"ERROR: {e}", ha="center", va="center",
                         transform=ax_full.transAxes)
            ax_crop.set_axis_off()
            continue

        # Full-frame DAPI z-sum with match center + template outline
        vmin, vmax = np.percentile(zsum, (1, 99.5))
        ax_full.imshow(zsum, cmap="gray", vmin=vmin, vmax=vmax)
        ax_full.plot(yx[1], yx[0], "r+", markersize=22, mew=2,
                     label=f"match @ ({yx[0]}, {yx[1]})  peak NCC={peak:.3f}")
        try:
            shifted = np.roll(template, (yx[0] - 1024, yx[1] - 1024), axis=(0, 1))
            cs = skimage.measure.find_contours(shifted)
            if cs:
                ax_full.plot(cs[0][:, 1], cs[0][:, 0], "r-", linewidth=1.0,
                             alpha=0.75, label="template outline")
        except Exception:
            pass
        ax_full.legend(loc="upper right", fontsize=10)
        ax_full.set_title(f"{tag}  full frame {zsum.shape}", fontsize=11)
        ax_full.set_axis_off()

        # Cropped panel at the matched center
        offset = [128, 128]
        y0, y1 = yx[0] - 512 + offset[0], yx[0] + 512 + offset[0]
        x0, x1 = yx[1] - 512 + offset[1], yx[1] + 512 + offset[1]
        y0, y1 = max(0, y0), min(zsum.shape[0], y1)
        x0, x1 = max(0, x0), min(zsum.shape[1], x1)
        crop = zsum[y0:y1, x0:x1]
        if crop.size > 0:
            v0, v1 = np.percentile(crop, (1, 99.5))
            ax_crop.imshow(crop, cmap="gray", vmin=v0, vmax=v1)
            ax_crop.set_title(
                f"{tag}  crop @ matched center  shape={crop.shape}  "
                f"{(y1-y0)}×{(x1-x0)}", fontsize=11)
        else:
            ax_crop.text(0.5, 0.5, "empty crop", ha="center", va="center",
                         transform=ax_crop.transAxes)
        ax_crop.set_axis_off()

    out = OUT_DIR / f"{pair['label']}_raw_vs_denoised.png"
    fig.suptitle(pair["label"], fontsize=13)
    plt.tight_layout(rect=[0, 0.02, 1, 0.96])
    fig.savefig(out, dpi=130)
    plt.close(fig)
    return summary, out


def main():
    template_hat = tmb.get_template_hat(1326)
    template = tmb.get_padded_template_at_width(1326)
    print("Cell                          | RAW peak  RAW (y,x)        | DENOISED peak  DENOISED (y,x)")
    print("-" * 100)
    for pair in PAIRS:
        s, out = render_pair(pair, template_hat, template)
        raw_yx = s.get("RAW_yx"); raw_p = s.get("RAW_peak")
        d_yx = s.get("DENOISED_yx"); d_p = s.get("DENOISED_peak")
        raw_str = f"{raw_p:.3f}  {raw_yx}" if raw_p is not None else "FAILED"
        d_str = f"{d_p:.3f}  {d_yx}" if d_p is not None else "FAILED"
        print(f"{s['label']:<28} | RAW: {raw_str:<28} | DENOISED: {d_str}")
        print(f"  → wrote {out.name}")


if __name__ == "__main__":
    main()
