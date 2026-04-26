"""Validate Mark's nucleus segmentation on a stratified sample of Fig 4B cells.

For each sampled cell:
  - Load ND2, z-sum the 405 and 488 channels
  - Run Mark's segmentation (Otsu on 405 zsum + largest connected component)
  - Compute quantitative seg quality metrics
  - Compute alternative segmentations (triangle, Li) for cross-method comparison
  - Save a PNG overlay (405 with nucleus contour; 488 with nucleus + perinuclear
    zone + pattern arch) and a per-cell JSON with seg metrics

Stratified sample: 2 cells per (plate × condition) for the TRAK isoform (mito)
sheet — 11 wells × 2 cells = up to 22 cells.

Idempotent per cell: writes <out>/<plate>/<well>/<cell>.{png,json}; skips if both
exist. SMB-drop-tolerant.
"""
from __future__ import annotations

import json
import pathlib
import sys
import traceback

import numpy as np
import polars as pl
import skimage
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrow  # noqa: F401
from scipy.ndimage import distance_transform_edt

import nd2

REPO = pathlib.Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))
import template_matching_bulk as tmb  # noqa: E402

COMBINED = REPO / "replication" / "overnight_out" / "combined.csv"
OUT_DIR = REPO / "replication" / "nucleus_seg_validation"
SHEET = "TRAK isoform (mito)"
CELLS_PER_STRATUM = 2


def pick_sample(df: pl.DataFrame) -> list[dict]:
    """Return up to 2 cells per (plate, condition) for the target sheet."""
    sub = df.filter(pl.col("sheet") == SHEET)
    sample = []
    for (plate, cond), grp in sub.group_by(["plate", "condition"]):
        # deterministic: take first 2 (alphabetically by path)
        picks = grp.sort("path").head(CELLS_PER_STRATUM)
        for r in picks.iter_rows(named=True):
            sample.append({
                "plate": r["plate"],
                "well": r["well"],
                "condition": r["condition"],
                "path": r["path"],
            })
    return sample


def segment_nucleus_otsu(nuc_zsum: np.ndarray) -> tuple[np.ndarray, float, dict]:
    """Mark's pipeline: Otsu + label + largest component. Also return diagnostics."""
    thr = skimage.filters.threshold_otsu(nuc_zsum)
    binary = nuc_zsum > thr
    label = skimage.measure.label(binary)
    props = skimage.measure.regionprops(label)
    if len(props) == 0:
        return np.zeros_like(binary), thr, {"n_components": 0}
    areas = [p.area for p in props]
    largest_idx = int(np.argmax(areas)) + 1
    mask = label == largest_idx
    # Diagnostic: total components and area distribution
    diag = {
        "n_components": len(props),
        "largest_frac_of_total_signal": float(areas[largest_idx - 1] / sum(areas)),
        "n_components_above_0_1_frac": int(sum(1 for a in areas if a / max(areas) > 0.1)),
    }
    return mask, float(thr), diag


def segment_nucleus_alt(nuc_zsum: np.ndarray, method: str) -> np.ndarray:
    if method == "triangle":
        thr = skimage.filters.threshold_triangle(nuc_zsum)
    elif method == "li":
        thr = skimage.filters.threshold_li(nuc_zsum)
    elif method == "mean":
        thr = nuc_zsum.mean()
    else:
        raise ValueError(method)
    binary = nuc_zsum > thr
    label = skimage.measure.label(binary)
    props = skimage.measure.regionprops(label)
    if len(props) == 0:
        return np.zeros_like(binary)
    largest_idx = int(np.argmax([p.area for p in props])) + 1
    return label == largest_idx


def seg_quality(mask: np.ndarray, pitch_um: float) -> dict:
    if not mask.any():
        return {"area_um2": 0.0, "valid": False}
    props = skimage.measure.regionprops(mask.astype(int))
    p = props[0]
    return {
        "area_um2": float(p.area * pitch_um ** 2),
        "area_px": int(p.area),
        "solidity": float(p.solidity),
        "euler_number": int(p.euler_number),
        "eccentricity": float(p.eccentricity),
        "bbox_aspect": float(max(p.bbox[2] - p.bbox[0], p.bbox[3] - p.bbox[1])
                             / max(1, min(p.bbox[2] - p.bbox[0], p.bbox[3] - p.bbox[1]))),
        "major_axis_um": float(p.major_axis_length * pitch_um),
        "minor_axis_um": float(p.minor_axis_length * pitch_um),
        "valid": True,
    }


def edge_vs_interior(nuc_zsum: np.ndarray, mask: np.ndarray) -> dict:
    """Is the Otsu cut at a sensible intensity level? Compare intensities in a
    shell just inside vs just outside the mask."""
    if not mask.any():
        return {"edge_ratio": float("nan")}
    # 2-pixel shell just outside
    edt_out = distance_transform_edt(~mask)
    outer_shell = (edt_out > 0) & (edt_out <= 2)
    # 2-pixel shell just inside
    edt_in = distance_transform_edt(mask)
    inner_shell = (edt_in > 0) & (edt_in <= 2)
    if not outer_shell.any() or not inner_shell.any():
        return {"edge_ratio": float("nan")}
    inner_mean = float(nuc_zsum[inner_shell].mean())
    outer_mean = float(nuc_zsum[outer_shell].mean())
    return {
        "inner_shell_mean": inner_mean,
        "outer_shell_mean": outer_mean,
        "edge_ratio": inner_mean / outer_mean if outer_mean > 0 else float("nan"),
    }


def render_overlay(img_path, nuc_zsum, mito_zsum, mask_otsu, mask_tri, mask_li,
                   contour, pitch_um, out_png):
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))

    # Panel 1: 405 zsum with Otsu segmentation contour
    ax = axes[0]
    vmax = np.percentile(nuc_zsum, 99.5)
    ax.imshow(nuc_zsum, cmap="Blues_r", vmax=vmax)
    # Otsu contour (cyan), triangle (yellow), Li (magenta)
    for msk, color, lbl in [(mask_otsu, "cyan", "Otsu (Mark)"),
                            (mask_tri, "yellow", "triangle"),
                            (mask_li, "magenta", "Li")]:
        if msk.any():
            contours = skimage.measure.find_contours(msk.astype(float), 0.5)
            for c in contours:
                ax.plot(c[:, 1], c[:, 0], color=color, linewidth=1.5, label=lbl)
                lbl = None
    ax.legend(loc="lower right", fontsize=8)
    ax.set_title(f"405 z-sum · nucleus seg\n{img_path.name}")
    ax.axis("off")

    # Panel 2: 488 zsum with Otsu nucleus contour + 5um perinuclear zone
    ax = axes[1]
    vmax = np.percentile(mito_zsum, 99.5)
    ax.imshow(mito_zsum, cmap="Greens_r", vmax=vmax)
    if mask_otsu.any():
        # Nucleus contour (cyan)
        contours = skimage.measure.find_contours(mask_otsu.astype(float), 0.5)
        for c in contours:
            ax.plot(c[:, 1], c[:, 0], color="cyan", linewidth=1.5)
        # 5 um perinuclear ring (red)
        edt = distance_transform_edt(~mask_otsu)
        ring_5um = edt * pitch_um
        ring_contour = skimage.measure.find_contours((ring_5um < 5).astype(float), 0.5)
        for c in ring_contour:
            ax.plot(c[:, 1], c[:, 0], color="red", linewidth=1.0, linestyle="--")
    # pattern arch contour (orange, dotted)
    ax.plot(contour[:, 1], contour[:, 0], color="orange", linewidth=0.8, linestyle=":", alpha=0.7)
    ax.set_title("488 z-sum · mito with seg + 5µm zone")
    ax.axis("off")

    # Panel 3: Seg disagreement heatmap — where do methods differ?
    ax = axes[2]
    disagreement = mask_otsu.astype(int) + mask_tri.astype(int) + mask_li.astype(int)
    # 0 = no method includes, 3 = all methods include
    ax.imshow(disagreement, cmap="RdYlGn", vmin=0, vmax=3)
    ax.set_title("method agreement (0=none, 3=all)")
    ax.axis("off")

    plt.tight_layout()
    fig.savefig(out_png, dpi=110, bbox_inches="tight")
    plt.close(fig)


def process_one(cell_info: dict, template_hat, template) -> dict | None:
    img_path = pathlib.Path(cell_info["path"])
    rel_out = f"{cell_info['plate']}_{cell_info['well']}_{img_path.stem}"
    out_png = OUT_DIR / f"{rel_out}.png"
    out_json = OUT_DIR / f"{rel_out}.json"
    if out_png.exists() and out_json.exists():
        return json.loads(out_json.read_text())

    img = nd2.imread(img_path, xarray=True)
    zsum = img.sum(axis=0)  # (C, Y, X)

    # Match Mark's crop framing
    key = tmb.cluster_key(img_path)
    offset = tmb.offset_overrides.get(key, [128, 128])
    roi = tmb.roi_overrides.get(key, None)
    max_coords = tmb.get_template_center(img, img_path, template_hat=template_hat,
                                         offset=offset, roi=roi)
    shifted_template = np.roll(template, (max_coords[0] - 1024, max_coords[1] - 1024),
                               axis=(0, 1))
    y_start, y_end = max_coords[0] - 512 + offset[0], max_coords[0] + 512 + offset[0]
    x_start, x_end = max_coords[1] - 512 + offset[1], max_coords[1] + 512 + offset[1]

    cropped_zsum = zsum.isel(Y=slice(y_start, y_end), X=slice(x_start, x_end))
    nuc = cropped_zsum.sel(C="405").to_numpy()
    mito = cropped_zsum.sel(C="488").to_numpy()

    # Pattern arch contour (for overlay)
    contour = skimage.measure.find_contours(shifted_template)[0].copy()
    contour[:, 0] -= max_coords[0] - 512
    contour[:, 1] -= max_coords[1] - 512

    # Segmentations
    mask_otsu, thr_otsu, otsu_diag = segment_nucleus_otsu(nuc)
    mask_tri = segment_nucleus_alt(nuc, "triangle")
    mask_li = segment_nucleus_alt(nuc, "li")

    pitch_um = img.metadata["metadata"].channels[0].volume.axesCalibration[0]

    q_otsu = seg_quality(mask_otsu, pitch_um)
    q_tri = seg_quality(mask_tri, pitch_um)
    q_li = seg_quality(mask_li, pitch_um)
    edge = edge_vs_interior(nuc, mask_otsu)

    agreement = {}
    if mask_otsu.any() and mask_tri.any():
        inter = (mask_otsu & mask_tri).sum()
        union = (mask_otsu | mask_tri).sum()
        agreement["iou_otsu_triangle"] = float(inter / union)
    if mask_otsu.any() and mask_li.any():
        inter = (mask_otsu & mask_li).sum()
        union = (mask_otsu | mask_li).sum()
        agreement["iou_otsu_li"] = float(inter / union)

    record = {
        **cell_info,
        "pitch_um": pitch_um,
        "otsu_threshold": thr_otsu,
        "otsu_diagnostics": otsu_diag,
        "otsu_quality": q_otsu,
        "triangle_quality": q_tri,
        "li_quality": q_li,
        "edge_vs_interior": edge,
        "agreement": agreement,
    }

    render_overlay(img_path, nuc, mito, mask_otsu, mask_tri, mask_li,
                   contour, pitch_um, out_png)
    out_json.write_text(json.dumps(record, indent=2))
    return record


def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    df = pl.read_csv(COMBINED)
    sample = pick_sample(df)
    print(f"[validate_nucleus_seg] {len(sample)} cells to process", flush=True)

    template_hat = tmb.get_template_hat(1326)
    template = tmb.get_padded_template_at_width(1326)

    records = []
    for i, cell in enumerate(sample, 1):
        try:
            r = process_one(cell, template_hat, template)
            records.append(r)
            if r is not None:
                q = r["otsu_quality"]
                print(f"  [{i}/{len(sample)}] OK  {cell['plate']}/{cell['well']}/"
                      f"{pathlib.Path(cell['path']).name}  "
                      f"area={q['area_um2']:.1f}µm²  "
                      f"solidity={q.get('solidity', 0):.3f}  "
                      f"n_components={r['otsu_diagnostics']['n_components']}",
                      flush=True)
        except Exception as e:
            print(f"  [{i}/{len(sample)}] ERR  {cell['path']}: {e}", flush=True)
            traceback.print_exc()

    # Write summary CSV
    flat_rows = []
    for r in records:
        if r is None:
            continue
        flat = {
            "plate": r["plate"], "well": r["well"], "condition": r["condition"],
            "cell": pathlib.Path(r["path"]).name,
            "area_um2": r["otsu_quality"]["area_um2"],
            "solidity": r["otsu_quality"].get("solidity", float("nan")),
            "euler_number": r["otsu_quality"].get("euler_number", -999),
            "eccentricity": r["otsu_quality"].get("eccentricity", float("nan")),
            "major_axis_um": r["otsu_quality"].get("major_axis_um", float("nan")),
            "minor_axis_um": r["otsu_quality"].get("minor_axis_um", float("nan")),
            "n_components": r["otsu_diagnostics"]["n_components"],
            "largest_frac": r["otsu_diagnostics"]["largest_frac_of_total_signal"],
            "n_components_significant": r["otsu_diagnostics"]["n_components_above_0_1_frac"],
            "edge_ratio": r["edge_vs_interior"].get("edge_ratio", float("nan")),
            "iou_otsu_triangle": r["agreement"].get("iou_otsu_triangle", float("nan")),
            "iou_otsu_li": r["agreement"].get("iou_otsu_li", float("nan")),
            "area_triangle_um2": r["triangle_quality"]["area_um2"],
            "area_li_um2": r["li_quality"]["area_um2"],
        }
        flat_rows.append(flat)
    if flat_rows:
        pl.from_dicts(flat_rows).write_csv(OUT_DIR / "summary.csv")
        print(f"\nSummary CSV: {OUT_DIR / 'summary.csv'}")
        print(f"Overlays: {OUT_DIR}/*.png")


if __name__ == "__main__":
    main()
