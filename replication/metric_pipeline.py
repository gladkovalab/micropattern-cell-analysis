"""Companion pipeline that augments Mark's per-cell analysis with additional
candidate distribution metrics — including MaxIP-based variants — for
exploring whether a better-motivated metric recovers stronger effects than
Mark's perinuclear_5um_percent_total.

Reuses template matching, cropping, nuclear segmentation, and the pattern
top-arch contour from `template_matching_bulk` so the spatial framework is
identical to Mark's. Only the projection choice (z-sum vs MaxIP) and the
metric extraction change.

Outputs to replication/metrics_out/{plate}/{well}/metrics.csv — never writes
into mark_data/.
"""
from __future__ import annotations

import argparse
import pathlib
import sys
import traceback

import nd2
import numpy as np
import polars as pl
import skimage
import xarray as xr
from scipy.ndimage import distance_transform_edt, center_of_mass

REPO = pathlib.Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))

import template_matching_bulk as tmb  # noqa: E402

DEFAULT_OUT = REPO / "replication" / "metrics_out"


# ------------- projection helpers -------------

def _sum_and_maxip(img_xr) -> tuple[xr.DataArray, xr.DataArray]:
    """Return (z-sum, MaxIP) xarrays over the Z dimension for all channels."""
    zsum = img_xr.sum(axis=0)  # axis 0 is Z in the (Z, C, Y, X) layout
    zmax = img_xr.max(axis=0)
    return zsum, zmax


def _bg_subtracted(mito_proj: np.ndarray, strip_px: int = 128) -> tuple[np.ndarray, float]:
    """Mark's background subtraction: 99.99th percentile of the left and right
    `strip_px` columns, with a robustness check for asymmetry.
    Returns (bg-subtracted image clipped to [0, ∞), threshold)."""
    img = mito_proj.astype(np.float64)
    left = np.percentile(img[:, :strip_px], 99.99)
    right = np.percentile(img[:, -strip_px:], 99.99)
    # Stretch01-compatible asymmetry check: if after normalizing to [0,1] by
    # max, the percentile difference exceeds 0.1, use min instead of the pooled
    # percentile. Approximate by translating the check into raw units.
    peak = img.max() - img.min()
    if peak > 0 and abs(left - right) / peak > 0.1:
        bg = min(left, right)
    else:
        both = np.concatenate((img[:, :strip_px], img[:, -strip_px:]), axis=1)
        bg = np.percentile(both, 99.99)
    out = np.clip(img - bg, 0, None)
    return out, float(bg)


# ------------- metric helpers -------------

def _gini(values: np.ndarray) -> float:
    """Gini coefficient over non-negative values (e.g., pixel intensities).
    Uses the standard sorted formula; ignores zeros only if all-zero."""
    v = np.asarray(values, dtype=np.float64).ravel()
    v = v[v >= 0]
    n = v.size
    if n == 0 or v.sum() == 0:
        return float("nan")
    v = np.sort(v)
    idx = np.arange(1, n + 1)
    return float((2 * (idx * v).sum() - (n + 1) * v.sum()) / (n * v.sum()))


def _radial_bin_fractions(intensity: np.ndarray, nuc_edt_um: np.ndarray, edges_um):
    """Fraction of total intensity in each annular bin (perinuclear distance)."""
    total = intensity.sum()
    if total <= 0:
        return [float("nan")] * (len(edges_um) - 1)
    out = []
    for lo, hi in zip(edges_um[:-1], edges_um[1:]):
        if hi is None:  # last open-ended bin
            mask = nuc_edt_um >= lo
        else:
            mask = (nuc_edt_um >= lo) & (nuc_edt_um < hi)
        out.append(float((mask * intensity).sum() / total * 100))
    return out


def _com_offset_um(intensity: np.ndarray, origin_yx: tuple[float, float],
                   pitch_um: float) -> tuple[float, float, float]:
    """Center-of-mass of the intensity image minus a reference origin.
    Returns (dy_um, dx_um, |d|_um). The reference origin is (y, x) in pixels."""
    if intensity.sum() <= 0:
        return float("nan"), float("nan"), float("nan")
    cy, cx = center_of_mass(intensity)
    dy = (cy - origin_yx[0]) * pitch_um
    dx = (cx - origin_yx[1]) * pitch_um
    return float(dy), float(dx), float(np.hypot(dy, dx))


def _apical_fraction(intensity: np.ndarray, pattern_com_y: float) -> float:
    """Fraction of intensity above the pattern CoM Y (closer to arch)."""
    total = intensity.sum()
    if total <= 0:
        return float("nan")
    ys = np.arange(intensity.shape[0])[:, None]
    above = intensity[ys[:, 0] < pattern_com_y, :].sum()
    return float(above / total * 100)


def _weighted_mean(values: np.ndarray, weights: np.ndarray) -> float:
    w = weights.sum()
    if w <= 0:
        return float("nan")
    return float((values * weights).sum() / w)


def _weighted_quantile(values: np.ndarray, weights: np.ndarray, q: float) -> float:
    w = weights.sum()
    if w <= 0:
        return float("nan")
    idx = np.argsort(values)
    vs = values[idx]
    ws = weights[idx]
    cum = np.cumsum(ws) / w
    return float(vs[np.searchsorted(cum, q)])


# ------------- main per-cell routine -------------

def process_cell(
    img_path: pathlib.Path,
    *,
    template_hat,
    template,
    out_root: pathlib.Path,
    save_projections: bool = False,
) -> dict:
    key = tmb.cluster_key(img_path)
    offset = tmb.offset_overrides.get(key, [128, 128])
    roi = tmb.roi_overrides.get(key, None)

    img = nd2.imread(img_path, xarray=True)

    zsum, zmax = _sum_and_maxip(img)

    # ---------- template matching on 640 z-sum (match Mark exactly) ----------
    # Reconstruct Mark's `score_template_match` template-matching prelude.
    sumproj = zsum.sel(C="640").to_numpy()[offset[0]:2048 + offset[0], offset[1]:2048 + offset[1]]
    sumproj_threshold = skimage.filters.threshold_otsu(sumproj)
    sumproj_thresholded = sumproj > sumproj_threshold
    max_coords = tmb.get_template_center(img, img_path, template_hat=template_hat,
                                         offset=offset, roi=roi)
    shifted_template = np.roll(template, (max_coords[0] - 1024, max_coords[1] - 1024),
                               axis=(0, 1))
    score = np.sum(sumproj_thresholded & shifted_template) / np.sum(shifted_template > 0)

    # ---------- crops ----------
    y_start, y_end = max_coords[0] - 512 + offset[0], max_coords[0] + 512 + offset[0]
    x_start, x_end = max_coords[1] - 512 + offset[1], max_coords[1] + 512 + offset[1]

    cropped_zsum = zsum.isel(Y=slice(y_start, y_end), X=slice(x_start, x_end))
    cropped_zmax = zmax.isel(Y=slice(y_start, y_end), X=slice(x_start, x_end))

    # template contour in crop coords
    contour = skimage.measure.find_contours(shifted_template)[0].copy()
    contour[:, 0] -= max_coords[0] - 512
    contour[:, 1] -= max_coords[1] - 512

    # ---------- nuclear segmentation (from z-sum 405) ----------
    nuc = cropped_zsum.sel(C="405").to_numpy()
    nuc_thresh = skimage.filters.threshold_otsu(nuc)
    nuc_mask = nuc > nuc_thresh
    nuc_label = skimage.measure.label(nuc_mask)
    nuc_props = skimage.measure.regionprops(nuc_label)
    nuc_max = int(np.argmax([p.area for p in nuc_props])) + 1
    nuc_mask = nuc_label == nuc_max
    nuc_edt_px = distance_transform_edt(np.invert(nuc_mask))
    nuc_com = center_of_mass(nuc_mask)  # (y, x) pixels

    # top-arch edt (Mark's arch index range 1083:1951)
    arch_px = np.zeros_like(nuc_mask)
    arch_px[
        np.round(contour[1083:1951, 0]).astype(int),
        np.round(contour[1083:1951, 1]).astype(int),
    ] = True
    arch_edt_px = distance_transform_edt(np.invert(arch_px))

    # pattern mask (inside template)
    pattern_mask = np.zeros_like(nuc_mask)
    rr = np.round(contour[:, 0]).astype(int)
    cc = np.round(contour[:, 1]).astype(int)
    # Create a filled pattern mask via the template: template is 2048×2048
    # and shifted so that pattern center is at `max_coords`. Crop to 1024.
    pattern_mask_big = shifted_template[y_start:y_end, x_start:x_end] > 0
    pattern_com = center_of_mass(pattern_mask_big)

    pitch_um = img.metadata["metadata"].channels[0].volume.axesCalibration[0]
    d_um = lambda px: px * pitch_um  # noqa: E731

    metrics: dict = {
        "path": str(img_path),
        "template_matching_score": float(score),
        "lateral_pixel_pitch_um": pitch_um,
    }

    # ---------- compute metrics for each projection variant ----------
    for proj_name, proj_xr in [("zsum", cropped_zsum), ("maxip", cropped_zmax)]:
        mito_raw = proj_xr.sel(C="488").to_numpy().astype(np.float64)
        mito_bg, bg_thr = _bg_subtracted(mito_raw)
        metrics[f"{proj_name}_bg_threshold"] = bg_thr
        # Weighted by intensity (bg-subtracted). Also do within-pattern variant.
        for mask_name, mask in [("crop", np.ones_like(mito_bg, dtype=bool)),
                                ("pattern", pattern_mask_big)]:
            I = mito_bg * mask
            total = I.sum()
            metrics[f"{proj_name}_{mask_name}_total_signal"] = float(total)
            if total <= 0:
                continue

            # Mark-style zones (to verify replication when zsum+crop)
            peri_mask = (nuc_edt_px < 5 / pitch_um) & mask
            per_mask = (arch_edt_px <= 5 / pitch_um) & mask
            metrics[f"{proj_name}_{mask_name}_perinuclear_5um_pct"] = float(
                (peri_mask * I).sum() / total * 100)
            metrics[f"{proj_name}_{mask_name}_peripheral_5um_pct"] = float(
                (per_mask * I).sum() / total * 100)

            # Radial bin profile (µm from nucleus)
            edges = [0, 2, 5, 10, 15, None]
            bins = _radial_bin_fractions(I, nuc_edt_px * pitch_um, edges)
            metrics[f"{proj_name}_{mask_name}_radial_0_2um_pct"] = bins[0]
            metrics[f"{proj_name}_{mask_name}_radial_2_5um_pct"] = bins[1]
            metrics[f"{proj_name}_{mask_name}_radial_5_10um_pct"] = bins[2]
            metrics[f"{proj_name}_{mask_name}_radial_10_15um_pct"] = bins[3]
            metrics[f"{proj_name}_{mask_name}_radial_ge15um_pct"] = bins[4]

            # CoM offset from nucleus CoM (y is apical→basal on pattern; positive
            # dy means mito below nucleus = away from arch)
            dy, dx, mag = _com_offset_um(I, nuc_com, pitch_um)
            metrics[f"{proj_name}_{mask_name}_com_dy_um"] = dy
            metrics[f"{proj_name}_{mask_name}_com_dx_um"] = dx
            metrics[f"{proj_name}_{mask_name}_com_offset_um"] = mag

            # Offset from pattern CoM (normalizes out cell-size variation)
            dy_p, dx_p, mag_p = _com_offset_um(I, pattern_com, pitch_um)
            metrics[f"{proj_name}_{mask_name}_com_vs_pattern_dy_um"] = dy_p
            metrics[f"{proj_name}_{mask_name}_com_vs_pattern_offset_um"] = mag_p

            # Apical fraction (above pattern CoM)
            metrics[f"{proj_name}_{mask_name}_apical_fraction_pct"] = _apical_fraction(I, pattern_com[0])

            # Weighted distance statistics (intensity-weighted distance to nucleus & arch)
            nuc_d_um = nuc_edt_px * pitch_um
            arch_d_um = arch_edt_px * pitch_um
            metrics[f"{proj_name}_{mask_name}_mean_dist_to_nucleus_um"] = _weighted_mean(nuc_d_um.ravel(), I.ravel())
            metrics[f"{proj_name}_{mask_name}_median_dist_to_nucleus_um"] = _weighted_quantile(nuc_d_um.ravel(), I.ravel(), 0.5)
            metrics[f"{proj_name}_{mask_name}_q90_dist_to_nucleus_um"] = _weighted_quantile(nuc_d_um.ravel(), I.ravel(), 0.9)
            metrics[f"{proj_name}_{mask_name}_mean_dist_to_arch_um"] = _weighted_mean(arch_d_um.ravel(), I.ravel())
            metrics[f"{proj_name}_{mask_name}_median_dist_to_arch_um"] = _weighted_quantile(arch_d_um.ravel(), I.ravel(), 0.5)

            # Gini over within-mask intensities
            metrics[f"{proj_name}_{mask_name}_gini"] = _gini(I[mask])

    # ---------- optional save of projections ----------
    if save_projections:
        rel = pathlib.Path(img_path).resolve().relative_to(pathlib.Path(tmb.DATA_ROOT).resolve())
        dest = out_root / "projections" / rel.with_suffix("")
        dest.parent.mkdir(parents=True, exist_ok=True)
        cropped_zsum.to_netcdf(dest.with_suffix(".zsum.nc"))
        cropped_zmax.to_netcdf(dest.with_suffix(".maxip.nc"))

    return metrics


# ------------- orchestration -------------

def iter_target_nd2s(root: pathlib.Path):
    """Walk `root`, yield ND2 files whose name starts with 'Cell' (case-insensitive).
    Prunes descent into MaxIP/MaxIPs/Excluded_cells/denoised subtrees so nested
    directories like `Excluded_cells/B06_250528/` are also skipped (the original
    `continue` only skipped direct children)."""
    skip = {"MaxIP", "MaxIPs", "Excluded_cells", "denoised"}
    for dirpath, dirnames, filenames in root.walk():
        dirnames[:] = [d for d in dirnames if d not in skip]
        for fn in filenames:
            if fn.endswith(".nd2") and fn.lower().startswith("cell"):
                yield dirpath / fn


def run(root: pathlib.Path, out_root: pathlib.Path, *, save_projections: bool = False) -> int:
    out_root.mkdir(parents=True, exist_ok=True)
    template_hat = tmb.get_template_hat(1326)
    template = tmb.get_padded_template_at_width(1326)

    records_by_well: dict[pathlib.Path, list[dict]] = {}
    cells = list(iter_target_nd2s(root))
    print(f"[metric_pipeline] {len(cells)} cells to process under {root}", flush=True)
    for i, img_path in enumerate(cells, 1):
        well_dir = img_path.parent
        try:
            m = process_cell(img_path, template_hat=template_hat, template=template,
                             out_root=out_root, save_projections=save_projections)
            records_by_well.setdefault(well_dir, []).append(m)
            print(f"  [{i}/{len(cells)}] OK {img_path.relative_to(root)}", flush=True)
        except Exception as e:
            print(f"  [{i}/{len(cells)}] ERR {img_path}: {e}", flush=True)
            traceback.print_exc()

    # write one CSV per well
    for well_dir, recs in records_by_well.items():
        rel = well_dir.resolve().relative_to(pathlib.Path(tmb.DATA_ROOT).resolve())
        out_csv = out_root / "by_well" / rel / "metrics.csv"
        out_csv.parent.mkdir(parents=True, exist_ok=True)
        pl.from_dicts(recs).write_csv(out_csv)
        print(f"[metric_pipeline] wrote {out_csv} ({len(recs)} cells)")

    # write combined CSV across all wells processed this run
    all_recs = [r for recs in records_by_well.values() for r in recs]
    if all_recs:
        combined = out_root / "combined" / f"{root.name}.csv"
        combined.parent.mkdir(parents=True, exist_ok=True)
        pl.from_dicts(all_recs).write_csv(combined)
        print(f"[metric_pipeline] wrote {combined} ({len(all_recs)} cells total)")
    return 0


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("root", help="Directory of ND2 files (walks recursively, skips denoised/MaxIPs/Excluded_cells)")
    ap.add_argument("--out-root", default=str(DEFAULT_OUT))
    ap.add_argument("--save-projections", action="store_true")
    args = ap.parse_args()
    return run(pathlib.Path(args.root).resolve(), pathlib.Path(args.out_root).resolve(),
               save_projections=args.save_projections)


if __name__ == "__main__":
    sys.exit(main())
