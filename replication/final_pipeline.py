"""Streamlined pipeline for the Fig 4 / S11 alt-metric pitch.

Computes only the metrics that earned their place in the v2 evaluation,
plus Mark's published baselines for back-compat. Drops every experimental
column that didn't perform (X-projection, fine nucleus-anchored radial,
angular sector Gini, multiple full-crop moment tensors, legacy
wrong-frame pattern_mask).

Per cell, per projection (zsum + maxip):

  Mark baselines (regression check)
    perinuclear_5um_pct, peripheral_5um_pct
    mean_dist_to_nucleus_um

  Y-axis projection (image Y vs pattern CoM Y)
    y_gini, y_entropy, y_sd_um, y_skew, y_mean_um
    60-bin y_profile_*um_pct  (image-Y axis, ±30 µm from pattern CoM)

  Wedge-r polar (Scheme 1 — pattern bottom apex, upper wedge)
    wedge_r_gini, wedge_r_entropy, wedge_r_ks_vs_uniform
    wedge_r_mean_um, wedge_r_sd_um, wedge_r_skew
    wedge_r_q25_um, wedge_r_q50_um, wedge_r_q75_um
    wedge_r_20_35um_frac_pct (perinuclear band)
    wedge_r_35_55um_frac_pct (peripheral/arch band)
    wedge_mt_apex_elongation, wedge_mt_apex_lam_max_um2, wedge_mt_apex_lam_min_um2
    60-bin wedge_r_*_um_pct (1 µm bins, 0..60)

  Diagnostics (cell-level, not per-projection)
    nuc_area_um2, nuc_solidity, nuc_eccentricity, nuc_n_components,
      nuc_largest_area_frac, nuc_euler_number
    pattern_{bottom,top,left,right}_{dy,dx}_um_from_nuc
    wedge_opening_deg, wedge_px_fraction
    lateral_pixel_pitch_um

Performance notes:
  * Y, X meshgrids and the pattern-extreme extraction are done ONCE per cell
    and shared across all metric helpers.
  * The wedge mask + the r-distance map are computed ONCE per cell (both
    z-sum and MaxIP reuse them — they're geometry, not channel-dependent).
  * Histograms use np.bincount on integer-binned indices (faster than
    a Python loop of mask+sum, used elsewhere in the v2 pipeline).
  * The legacy wrong-frame pattern_mask_big is not built at all.

Entry points:
  * `process_cell(img_path, ..., target_wells=...)` — single cell
  * `run(root, out_root, target_wells)` — driver with per-cell checkpoint
  * CLI: --sheets 'Sheet name' [...] reads `combined.csv` to pick wells
"""
from __future__ import annotations

import argparse
import json
import pathlib
import sys
import traceback

import nd2
import numpy as np
import polars as pl
import skimage
from scipy.ndimage import distance_transform_edt, center_of_mass

REPO = pathlib.Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))

import template_matching_bulk as tmb  # noqa: E402

DEFAULT_OUT = REPO / "replication" / "overnight_final_out"

# Y-projection: 1 µm bins, ±30 µm around pattern CoM Y
Y_HALF_RANGE_UM = 30.0
Y_STEP_UM = 1.0
# Wedge-r: 1 µm bins, 0..60 µm
WEDGE_R_MAX_UM = 60.0
WEDGE_R_STEP_UM = 1.0


# -------------------------------------------------------------------- helpers

def _bg_subtracted(mito_proj: np.ndarray, strip_px: int = 128) -> tuple[np.ndarray, float]:
    """Mark's 99.99-percentile background subtraction (preserved from his pipeline)."""
    img = mito_proj.astype(np.float64)
    left = np.percentile(img[:, :strip_px], 99.99)
    right = np.percentile(img[:, -strip_px:], 99.99)
    peak = img.max() - img.min()
    if peak > 0 and abs(left - right) / peak > 0.1:
        bg = min(left, right)
    else:
        both = np.concatenate((img[:, :strip_px], img[:, -strip_px:]), axis=1)
        bg = np.percentile(both, 99.99)
    return np.clip(img - bg, 0, None), float(bg)


def _gini_of_bins(pct_bins: np.ndarray) -> float:
    p = np.asarray(pct_bins, dtype=np.float64)
    p = p[~np.isnan(p)]
    p = p[p >= 0]
    if p.size == 0 or p.sum() <= 0:
        return float("nan")
    p = np.sort(p)
    n = p.size
    idx = np.arange(1, n + 1)
    return float((2 * (idx * p).sum() - (n + 1) * p.sum()) / (n * p.sum()))


def _entropy_of_bins(pct_bins: np.ndarray) -> float:
    p = np.asarray(pct_bins, dtype=np.float64)
    p = p[~np.isnan(p)]
    s = p.sum()
    if s <= 0:
        return float("nan")
    p = p / s
    p = p[p > 0]
    return float(-(p * np.log(p)).sum())


def _ks_vs_uniform(pct_bins: np.ndarray, vol_bins: np.ndarray) -> float:
    p = np.asarray(pct_bins, dtype=np.float64)
    v = np.asarray(vol_bins, dtype=np.float64)
    if np.nansum(p) <= 0 or v.sum() <= 0:
        return float("nan")
    p = np.where(np.isnan(p), 0.0, p)
    cdf_obs = np.cumsum(p) / p.sum()
    cdf_uni = np.cumsum(v) / v.sum()
    return float(np.max(np.abs(cdf_obs - cdf_uni)))


def _bincount_intensity(intensity_flat: np.ndarray, bin_idx_flat: np.ndarray,
                        n_bins: int, mask_flat: np.ndarray | None = None) -> np.ndarray:
    """Sum intensity per bin index using np.bincount. ~5× faster than the
    mask+sum loop used in earlier pipelines."""
    if mask_flat is not None:
        I = intensity_flat * mask_flat
    else:
        I = intensity_flat
    return np.bincount(bin_idx_flat, weights=I, minlength=n_bins)[:n_bins]


# ---------------------------------------------------------- pattern + wedge

def _pattern_extremes(shifted_template: np.ndarray, max_coords) -> dict:
    """Correctly-aligned: slice shifted_template in its own frame so the
    extracted points sit on the authoritative orange contour."""
    pattern_mask = shifted_template[
        max_coords[0] - 512:max_coords[0] + 512,
        max_coords[1] - 512:max_coords[1] + 512
    ] > 0
    ys, xs = np.where(pattern_mask)
    if ys.size == 0:
        return None
    return {
        "bottom": (int(ys[np.argmax(ys)]), int(xs[np.argmax(ys)])),
        "top":    (int(ys[np.argmin(ys)]), int(xs[np.argmin(ys)])),
        "left":   (int(ys[np.argmin(xs)]), int(xs[np.argmin(xs)])),
        "right":  (int(ys[np.argmax(xs)]), int(xs[np.argmax(xs)])),
    }


def _build_wedge_geometry(shape: tuple, ext: dict, pitch_um: float):
    """Compute the wedge mask + per-pixel r in µm from the wedge apex.
    Returns (wedge_mask, r_um, opening_rad, dy_grid_um, dx_grid_um) with all
    arrays the same shape as `shape`."""
    H, W = shape
    Y_idx, X_idx = np.mgrid[:H, :W]
    apex = ext["bottom"]
    dy_um = (Y_idx - apex[0]) * pitch_um
    dx_um = (X_idx - apex[1]) * pitch_um
    r_um = np.hypot(dy_um, dx_um)

    # Wedge angle bounds (rays from apex through left/right pattern extremes,
    # arc that contains "up" = dy<0)
    ang = np.arctan2(Y_idx - apex[0], X_idx - apex[1])
    a_left = float(np.arctan2(ext["left"][0] - apex[0], ext["left"][1] - apex[1]))
    a_right = float(np.arctan2(ext["right"][0] - apex[0], ext["right"][1] - apex[1]))
    a_up = -np.pi / 2
    lo = min(a_left, a_right)
    hi = max(a_left, a_right)
    if lo <= a_up <= hi:
        wedge = (ang >= lo) & (ang <= hi)
        opening = hi - lo
    else:
        wedge = (ang <= lo) | (ang >= hi)
        opening = (2 * np.pi) - (hi - lo)

    return wedge, r_um, float(opening), dy_um, dx_um


# --------------------------------------------------------- per-projection

def _projection_metrics(I: np.ndarray, *, proj_name: str,
                        nuc_edt_um: np.ndarray, arch_edt_um: np.ndarray,
                        pattern_com_y: float, pitch_um: float,
                        wedge_mask: np.ndarray, r_um: np.ndarray,
                        wedge_apex_dy: np.ndarray, wedge_apex_dx: np.ndarray
                        ) -> dict:
    """Compute the keeper metrics for one projection (zsum or maxip)."""
    out: dict = {}
    total = float(I.sum())
    out[f"{proj_name}_total_signal"] = total
    if total <= 0:
        return out

    # ---- Mark baselines ----
    out[f"{proj_name}_perinuclear_5um_pct"] = float(((nuc_edt_um < 5) * I).sum() / total * 100)
    out[f"{proj_name}_peripheral_5um_pct"] = float(((arch_edt_um <= 5) * I).sum() / total * 100)
    # Intensity-weighted mean distance to nucleus (most-cited single-zone)
    out[f"{proj_name}_mean_dist_to_nucleus_um"] = float(
        (I * nuc_edt_um).sum() / total)

    # ---- Y-axis projection ----
    # Collapse X to get I(y); coordinate is dy from pattern CoM
    I_y = I.sum(axis=1)  # shape (H,)
    H = I_y.size
    y_px = np.arange(H)
    dy_um_1d = (y_px - pattern_com_y) * pitch_um
    edges = np.arange(-Y_HALF_RANGE_UM, Y_HALF_RANGE_UM + Y_STEP_UM, Y_STEP_UM)
    n_y = len(edges) - 1
    # Map each y-row to a bin index; rows outside the range get bin = -1
    bin_idx = np.floor((dy_um_1d - edges[0]) / Y_STEP_UM).astype(int)
    in_range = (bin_idx >= 0) & (bin_idx < n_y)
    y_profile = np.zeros(n_y)
    y_profile[bin_idx[in_range]] += I_y[in_range]
    # Note: multiple rows can map to the same bin (when binning is finer than
    # row spacing they map 1:1; coarser, multiple rows accumulate). We use
    # add.at to handle the ambiguity correctly.
    y_profile = np.zeros(n_y)
    np.add.at(y_profile, bin_idx[in_range], I_y[in_range])
    y_total = float(y_profile.sum())
    if y_total > 0:
        y_pct = y_profile / y_total * 100
        for i, pct in enumerate(y_pct):
            out[f"{proj_name}_y_profile_{int(edges[i]):+04d}um_pct"] = float(pct)
        # Y-axis moments (computed on the I_y signal directly, not the binned
        # profile — same answer if pixel pitch is fine enough)
        I_y_in = I_y[in_range].astype(np.float64)
        dy_in = dy_um_1d[in_range]
        w = float(I_y_in.sum())
        mean_y = float((I_y_in * dy_in).sum() / w)
        var_y = float((I_y_in * (dy_in - mean_y) ** 2).sum() / w)
        sd_y = float(np.sqrt(max(0.0, var_y)))
        out[f"{proj_name}_y_mean_um"] = mean_y
        out[f"{proj_name}_y_sd_um"] = sd_y
        out[f"{proj_name}_y_skew"] = float(
            (I_y_in * ((dy_in - mean_y) / sd_y) ** 3).sum() / w) if sd_y > 0 else float("nan")
        out[f"{proj_name}_y_gini"] = _gini_of_bins(y_pct)
        out[f"{proj_name}_y_entropy"] = _entropy_of_bins(y_pct)

    # ---- Wedge-r polar (Scheme 1) ----
    # Bin r within the wedge using bincount (faster than per-bin masks)
    r_max = WEDGE_R_MAX_UM
    bin_idx_r = np.floor(r_um / WEDGE_R_STEP_UM).astype(int)
    n_r = int(r_max / WEDGE_R_STEP_UM)
    in_r = (bin_idx_r >= 0) & (bin_idx_r < n_r) & wedge_mask
    I_w = I * wedge_mask
    total_in = float(I_w.sum())
    out[f"{proj_name}_wedge_frac_pct"] = float(total_in / total * 100)
    if total_in <= 0:
        return out

    r_profile = np.bincount(
        bin_idx_r[in_r].ravel(),
        weights=I_w[in_r].ravel(),
        minlength=n_r,
    )[:n_r]
    r_pct = r_profile / total_in * 100
    vol_arc = np.bincount(bin_idx_r[in_r].ravel(),
                          weights=np.ones_like(I_w[in_r].ravel()),
                          minlength=n_r)[:n_r]
    for i in range(n_r):
        out[f"{proj_name}_wedge_r_{i:02d}_{i+1:02d}um_pct"] = float(r_pct[i])

    # Moments and quantiles on raw (pixel-level) r within the wedge
    I_flat = I_w.ravel()
    r_flat = r_um.ravel()
    valid = I_flat > 0
    if valid.any():
        I_v = I_flat[valid]
        r_v = r_flat[valid]
        w = float(I_v.sum())
        mean_r = float((I_v * r_v).sum() / w)
        var_r = float((I_v * (r_v - mean_r) ** 2).sum() / w)
        sd_r = float(np.sqrt(max(0.0, var_r)))
        out[f"{proj_name}_wedge_r_mean_um"] = mean_r
        out[f"{proj_name}_wedge_r_sd_um"] = sd_r
        if sd_r > 0:
            out[f"{proj_name}_wedge_r_skew"] = float(
                (I_v * ((r_v - mean_r) / sd_r) ** 3).sum() / w)
        else:
            out[f"{proj_name}_wedge_r_skew"] = float("nan")
        # Quantiles
        idx = np.argsort(r_v)
        cum = np.cumsum(I_v[idx]) / w
        for q in (0.25, 0.50, 0.75):
            i = int(np.searchsorted(cum, q))
            i = min(i, idx.size - 1)
            out[f"{proj_name}_wedge_r_q{int(q*100):02d}_um"] = float(r_v[idx[i]])

    # Distribution-shape scalars
    out[f"{proj_name}_wedge_r_gini"] = _gini_of_bins(r_pct)
    out[f"{proj_name}_wedge_r_entropy"] = _entropy_of_bins(r_pct)
    out[f"{proj_name}_wedge_r_ks_vs_uniform"] = _ks_vs_uniform(r_pct, vol_arc)

    # Band fractions (perinuclear / peripheral)
    bin_lo_idx = lambda lo: int(lo / WEDGE_R_STEP_UM)
    out[f"{proj_name}_wedge_r_20_35um_frac_pct"] = float(r_pct[bin_lo_idx(20):bin_lo_idx(35)].sum())
    out[f"{proj_name}_wedge_r_35_55um_frac_pct"] = float(r_pct[bin_lo_idx(35):bin_lo_idx(55)].sum())

    # Wedge-restricted moment tensor (anchored at wedge apex)
    if w > 0:
        # Use the precomputed dy/dx grids from the wedge geometry
        dy_w = wedge_apex_dy * wedge_mask
        dx_w = wedge_apex_dx * wedge_mask
        syy = float((I_w * dy_w * dy_w).sum() / w)
        sxx = float((I_w * dx_w * dx_w).sum() / w)
        sxy = float((I_w * dy_w * dx_w).sum() / w)
        cov = np.array([[syy, sxy], [sxy, sxx]])
        eigvals, _ = np.linalg.eigh(cov)
        lam_min, lam_max = float(eigvals[0]), float(eigvals[-1])
        out[f"{proj_name}_wedge_mt_apex_lam_max_um2"] = lam_max
        out[f"{proj_name}_wedge_mt_apex_lam_min_um2"] = lam_min
        out[f"{proj_name}_wedge_mt_apex_elongation"] = float(np.sqrt(
            max(lam_max, 1e-12) / max(lam_min, 1e-12)))

    return out


# ---------------------------------------------------------- per-cell entry

def process_cell(img_path: pathlib.Path, *, template_hat, template) -> dict:
    key = tmb.cluster_key(img_path)
    offset = tmb.offset_overrides.get(key, [128, 128])
    roi = tmb.roi_overrides.get(key, None)

    img = nd2.imread(img_path, xarray=True)
    zsum = img.sum(axis=0)
    zmax = img.max(axis=0)

    max_coords = tmb.get_template_center(img, img_path, template_hat=template_hat,
                                         offset=offset, roi=roi)
    shifted_template = np.roll(template, (max_coords[0] - 1024, max_coords[1] - 1024),
                               axis=(0, 1))

    y_start, y_end = max_coords[0] - 512 + offset[0], max_coords[0] + 512 + offset[0]
    x_start, x_end = max_coords[1] - 512 + offset[1], max_coords[1] + 512 + offset[1]
    cropped_zsum = zsum.isel(Y=slice(y_start, y_end), X=slice(x_start, x_end))
    cropped_zmax = zmax.isel(Y=slice(y_start, y_end), X=slice(x_start, x_end))

    contour = skimage.measure.find_contours(shifted_template)[0].copy()
    contour[:, 0] -= max_coords[0] - 512
    contour[:, 1] -= max_coords[1] - 512

    nuc = cropped_zsum.sel(C="405").to_numpy()
    nuc_mask = nuc > skimage.filters.threshold_otsu(nuc)
    nuc_label = skimage.measure.label(nuc_mask)
    nuc_props = skimage.measure.regionprops(nuc_label)
    if len(nuc_props) == 0:
        raise ValueError("nuclear segmentation produced 0 components")
    nuc_max = int(np.argmax([p.area for p in nuc_props])) + 1
    nuc_mask = nuc_label == nuc_max
    nuc_edt_px = distance_transform_edt(np.invert(nuc_mask))
    nuc_com = center_of_mass(nuc_mask)

    arch_px = np.zeros_like(nuc_mask)
    arch_px[np.round(contour[1083:1951, 0]).astype(int),
            np.round(contour[1083:1951, 1]).astype(int)] = True
    arch_edt_px = distance_transform_edt(np.invert(arch_px))

    pitch_um = img.metadata["metadata"].channels[0].volume.axesCalibration[0]
    nuc_edt_um = nuc_edt_px * pitch_um
    arch_edt_um = arch_edt_px * pitch_um

    # Pattern geometry (correctly-aligned). Pattern CoM Y in the crop frame is
    # 512 (because shifted_template's pattern is at max_coords, and the crop
    # starts at max_coords - 512 + offset; the pattern's actual position in the
    # crop is 512 - offset on the Y axis. Use that for the Y-projection origin
    # so it lines up with the orange contour.)
    ext = _pattern_extremes(shifted_template, max_coords)
    if ext is None:
        raise ValueError("pattern mask produced 0 nonzero pixels")
    # In the correctly-aligned mask, the pattern center is at (512, 512).
    pattern_com_y = 512.0  # by construction

    # Wedge geometry shared across both projections
    wedge_mask, r_um, opening_rad, dy_grid_um, dx_grid_um = _build_wedge_geometry(
        nuc.shape, ext, pitch_um)

    # Cell-level diagnostics
    nuc_region = skimage.measure.regionprops(nuc_mask.astype(int))[0]
    metrics: dict = {
        "path": str(img_path),
        "lateral_pixel_pitch_um": pitch_um,
        "nuc_area_um2": float(nuc_region.area * pitch_um ** 2),
        "nuc_solidity": float(nuc_region.solidity),
        "nuc_eccentricity": float(nuc_region.eccentricity),
        "nuc_euler_number": int(nuc_region.euler_number),
        "nuc_n_components": int(len(nuc_props)),
        "nuc_largest_area_frac": float(
            nuc_props[nuc_max - 1].area / sum(p.area for p in nuc_props)),
        "wedge_opening_deg": float(np.degrees(opening_rad)),
        "wedge_px_fraction": float(wedge_mask.sum() / wedge_mask.size),
    }
    for name in ("bottom", "top", "left", "right"):
        pt = ext[name]
        metrics[f"pattern_{name}_dy_um_from_nuc"] = float((pt[0] - nuc_com[0]) * pitch_um)
        metrics[f"pattern_{name}_dx_um_from_nuc"] = float((pt[1] - nuc_com[1]) * pitch_um)

    for proj_name, proj_xr in [("zsum", cropped_zsum), ("maxip", cropped_zmax)]:
        mito_raw = proj_xr.sel(C="488").to_numpy().astype(np.float64)
        I, bg_thr = _bg_subtracted(mito_raw)
        metrics[f"{proj_name}_bg_threshold"] = bg_thr
        m = _projection_metrics(
            I,
            proj_name=proj_name,
            nuc_edt_um=nuc_edt_um,
            arch_edt_um=arch_edt_um,
            pattern_com_y=pattern_com_y,
            pitch_um=pitch_um,
            wedge_mask=wedge_mask,
            r_um=r_um,
            wedge_apex_dy=dy_grid_um,
            wedge_apex_dx=dx_grid_um,
        )
        metrics.update(m)

    return metrics


# ---------------------------------------------------------- orchestration

def discover_wells(combined_csv: pathlib.Path,
                   sheets: list[str]) -> list[tuple[str, str]]:
    df = pl.read_csv(combined_csv).filter(pl.col("sheet").is_in(sheets))
    pairs = set()
    for plate, well in df.select(["plate", "well"]).unique().iter_rows():
        pairs.add((plate, f"{well}_"))
    return sorted(pairs)


def iter_cells(root: pathlib.Path, target_wells: list[tuple[str, str]]):
    skip = {"MaxIP", "MaxIPs", "Excluded_cells", "denoised"}
    for plate, well_prefix in target_wells:
        plate_dir = root / plate
        if not plate_dir.exists():
            print(f"  MISSING plate dir: {plate_dir}", flush=True)
            continue
        for well_dir in plate_dir.iterdir():
            if not well_dir.is_dir() or not well_dir.name.startswith(well_prefix):
                continue
            for dirpath, dirnames, filenames in well_dir.walk():
                dirnames[:] = [d for d in dirnames if d not in skip]
                for fn in filenames:
                    if fn.endswith(".nd2") and fn.lower().startswith("cell"):
                        yield dirpath / fn


def run(root: pathlib.Path, out_root: pathlib.Path,
        target_wells: list[tuple[str, str]]) -> int:
    out_root.mkdir(parents=True, exist_ok=True)
    template_hat = tmb.get_template_hat(1326)
    template = tmb.get_padded_template_at_width(1326)

    cells = list(iter_cells(root, target_wells))
    print(f"[final_pipeline] {len(cells)} cells to process under {root}", flush=True)

    records_by_well: dict[pathlib.Path, list[dict]] = {}
    for i, img_path in enumerate(cells, 1):
        well_dir = img_path.parent
        rel = well_dir.resolve().relative_to(pathlib.Path(tmb.DATA_ROOT).resolve())
        cell_chk = out_root / "by_well" / rel / "cells" / f"{img_path.stem}.json"
        cell_chk.parent.mkdir(parents=True, exist_ok=True)

        if cell_chk.exists():
            try:
                m = json.loads(cell_chk.read_text())
                records_by_well.setdefault(well_dir, []).append(m)
                print(f"  [{i}/{len(cells)}] CACHED {img_path.relative_to(root)}", flush=True)
                continue
            except Exception:
                pass

        try:
            m = process_cell(img_path, template_hat=template_hat, template=template)
            records_by_well.setdefault(well_dir, []).append(m)
            cell_chk.write_text(json.dumps(m))
            print(f"  [{i}/{len(cells)}] OK {img_path.relative_to(root)}", flush=True)
        except Exception as e:
            print(f"  [{i}/{len(cells)}] ERR {img_path}: {e}", flush=True)
            traceback.print_exc()

    for well_dir, recs in records_by_well.items():
        rel = well_dir.resolve().relative_to(pathlib.Path(tmb.DATA_ROOT).resolve())
        out_csv = out_root / "by_well" / rel / "metrics.csv"
        out_csv.parent.mkdir(parents=True, exist_ok=True)
        pl.from_dicts(recs).write_csv(out_csv)
        (out_root / "by_well" / rel / "done.marker").touch()
        print(f"[final_pipeline] wrote {out_csv} ({len(recs)} cells)", flush=True)

    all_recs = [r for recs in records_by_well.values() for r in recs]
    if all_recs:
        combined = out_root / "combined_raw.csv"
        pl.from_dicts(all_recs).write_csv(combined)
        print(f"[final_pipeline] wrote {combined} ({len(all_recs)} cells total)",
              flush=True)
    return 0


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", default=str(pathlib.Path(tmb.DATA_ROOT)))
    ap.add_argument("--out-root", default=str(DEFAULT_OUT))
    ap.add_argument("--combined-csv",
                    default=str(REPO / "replication" / "overnight_out" / "combined.csv"),
                    help="reference CSV used to look up plate/well/sheet/condition")
    ap.add_argument("--sheets", nargs="+", required=True,
                    help="sheets to process, e.g. 'TRAK isoform (mito)' "
                         "'TRAK1 helix muts' 'TRAK2 helix muts' 'MAPK9 siRNA' "
                         "'TRAK isoform (peroxisomes)'")
    args = ap.parse_args()
    target_wells = discover_wells(pathlib.Path(args.combined_csv), args.sheets)
    print(f"[final_pipeline] sheets: {args.sheets}")
    print(f"[final_pipeline] {len(target_wells)} wells across "
          f"{len({p for p, _ in target_wells})} plates", flush=True)
    return run(pathlib.Path(args.root).resolve(),
               pathlib.Path(args.out_root).resolve(),
               target_wells)


if __name__ == "__main__":
    sys.exit(main())
