"""Focused Fig 4B rerun with the extended metric set the user requested:

- Fine 1-µm radial profile (nucleus-anchored, 0-25 µm)
- Radial moments: σ_r, skew_r, kurt_r
- Radial quantiles: Q10, Q25, Q50, Q75, Q90
- Radial Gini (fine), radial entropy, radial KS vs uniform-by-area
- 2D moment-tensor eigenvalues around nucleus CoM + pattern CoM (elongation, λ_max, λ_min)
- Angular-sector Gini (n=8) around nucleus CoM + pattern CoM
- Y-axis projection: Y-CoM, σ_y, skew_y, kurt_y, apical/basal ratio, Y-Gini, Y-entropy,
  + 1-µm Y-profile across pattern span
- X-axis projection: same scalars (symmetry sanity check)

All metrics computed on both z-sum and MaxIP projections, crop mask only (pattern
mask was shown to be confounded — it discards legitimate fan-shape signal).
Existing Mark-style metrics (peri_5/nuc_5, 5-bin radial profile, CoM, Q50/Q90,
apical_fraction, pixel-intensity Gini) preserved for direct comparison.

Processes only the 11 wells involved in the TRAK isoform (mito) sheet (Fig 4B).
Output: replication/overnight_fig4b_out/
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
import xarray as xr
from scipy.ndimage import distance_transform_edt, center_of_mass

REPO = pathlib.Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))

import template_matching_bulk as tmb  # noqa: E402
from replication.metric_pipeline import (  # noqa: E402
    _bg_subtracted, _gini, _com_offset_um, _apical_fraction,
    _weighted_mean, _weighted_quantile, _radial_bin_fractions,
)

DEFAULT_OUT = REPO / "replication" / "overnight_fig4b_out"
SHEET = "TRAK isoform (mito)"

# 11 wells × 3 conditions (no TRAK, TRAK1, TRAK2) for Fig 4B
TARGET_WELLS = [
    # (plate, well-dir-prefix) — resolved against patterned_data root
    ("250521_patterned_plate_1", "B06_"),
    ("250606_patterned_plate_2", "C02_"),
    ("250612_patterned_plate_3", "B02_"),
    ("250612_patterned_plate_3", "B03_"),
    ("250612_patterned_plate_3", "B04_"),
    ("250710_patterned_plate_9_good", "C02_"),
    ("250710_patterned_plate_9_good", "C03_"),
    ("250710_patterned_plate_9_good", "C04_"),
    ("250731_patterned_plate_11_good", "D06_"),
    ("250731_patterned_plate_11_good", "E05_"),
    ("250731_patterned_plate_11_good", "F05_"),
]

FINE_R_MAX_UM = 25.0
FINE_R_STEP_UM = 1.0
Y_HALF_RANGE_UM = 30.0
Y_STEP_UM = 1.0
N_SECTORS = 8


# -------- New metric helpers --------

def _fine_radial_profile(intensity, nuc_edt_um, max_r_um=FINE_R_MAX_UM, step_um=FINE_R_STEP_UM):
    edges = np.arange(0, max_r_um + step_um, step_um)
    total = intensity.sum()
    n_bins = len(edges) - 1
    out_pct = np.full(n_bins, np.nan)
    out_vol_px = np.zeros(n_bins)
    if total <= 0:
        return out_pct, out_vol_px, edges
    for i, (lo, hi) in enumerate(zip(edges[:-1], edges[1:])):
        mask = (nuc_edt_um >= lo) & (nuc_edt_um < hi)
        out_pct[i] = (mask * intensity).sum() / total * 100
        out_vol_px[i] = mask.sum()
    return out_pct, out_vol_px, edges


def _radial_moments(intensity, nuc_edt_um):
    I = intensity.ravel()
    r = nuc_edt_um.ravel()
    w = I.sum()
    nan = dict(mean_r=np.nan, sd_r=np.nan, skew_r=np.nan, kurt_r=np.nan)
    if w <= 0:
        return nan
    mean_r = float((I * r).sum() / w)
    var_r = float((I * (r - mean_r) ** 2).sum() / w)
    sd_r = float(np.sqrt(max(0.0, var_r)))
    if sd_r <= 0:
        return dict(mean_r=mean_r, sd_r=sd_r, skew_r=np.nan, kurt_r=np.nan)
    skew_r = float((I * ((r - mean_r) / sd_r) ** 3).sum() / w)
    kurt_r = float((I * ((r - mean_r) / sd_r) ** 4).sum() / w - 3)
    return dict(mean_r=mean_r, sd_r=sd_r, skew_r=skew_r, kurt_r=kurt_r)


def _radial_quantiles(intensity, nuc_edt_um, qs=(0.10, 0.25, 0.50, 0.75, 0.90)):
    I = intensity.ravel()
    r = nuc_edt_um.ravel()
    w = I.sum()
    if w <= 0:
        return {f"q{int(q*100):02d}": np.nan for q in qs}
    idx = np.argsort(r)
    rs = r[idx]
    ws = I[idx]
    cum = np.cumsum(ws) / w
    out = {}
    for q in qs:
        out[f"q{int(q*100):02d}"] = float(rs[np.searchsorted(cum, q)])
    return out


def _gini_of_bins(pct_bins):
    p = np.asarray(pct_bins, dtype=np.float64)
    p = p[~np.isnan(p)]
    p = p[p >= 0]
    if p.size == 0 or p.sum() <= 0:
        return np.nan
    p = np.sort(p)
    n = p.size
    idx = np.arange(1, n + 1)
    return float((2 * (idx * p).sum() - (n + 1) * p.sum()) / (n * p.sum()))


def _entropy_of_bins(pct_bins):
    p = np.asarray(pct_bins, dtype=np.float64)
    p = p[~np.isnan(p)]
    s = p.sum()
    if s <= 0:
        return np.nan
    p = p / s
    p = p[p > 0]
    return float(-(p * np.log(p)).sum())


def _ks_vs_uniform(pct_bins, vol_px):
    p = np.asarray(pct_bins, dtype=np.float64)
    v = np.asarray(vol_px, dtype=np.float64)
    if np.nansum(p) <= 0 or v.sum() <= 0:
        return np.nan
    p = np.where(np.isnan(p), 0, p)
    cdf_obs = np.cumsum(p) / p.sum()
    cdf_uni = np.cumsum(v) / v.sum()
    return float(np.max(np.abs(cdf_obs - cdf_uni)))


def _moment_tensor(intensity, origin_yx, pitch_um):
    I = intensity
    w = I.sum()
    nan = dict(lam_max_um2=np.nan, lam_min_um2=np.nan, elongation=np.nan, orient_rad=np.nan)
    if w <= 0:
        return nan
    Y, X = np.mgrid[:I.shape[0], :I.shape[1]]
    dy = (Y - origin_yx[0]) * pitch_um
    dx = (X - origin_yx[1]) * pitch_um
    syy = (I * dy * dy).sum() / w
    sxx = (I * dx * dx).sum() / w
    sxy = (I * dy * dx).sum() / w
    cov = np.array([[syy, sxy], [sxy, sxx]])
    eigvals, eigvecs = np.linalg.eigh(cov)
    lam_min, lam_max = float(eigvals[0]), float(eigvals[-1])
    elongation = float(np.sqrt(max(lam_max, 1e-12) / max(lam_min, 1e-12)))
    v = eigvecs[:, -1]
    orient = float(np.arctan2(v[0], v[1]))
    return dict(lam_max_um2=lam_max, lam_min_um2=lam_min,
                elongation=elongation, orient_rad=orient)


def _angular_sector_gini(intensity, origin_yx, n_sectors=N_SECTORS):
    I = intensity
    total = I.sum()
    if total <= 0:
        return np.nan
    Y, X = np.mgrid[:I.shape[0], :I.shape[1]]
    theta = np.arctan2(Y - origin_yx[0], X - origin_yx[1])  # [-pi, pi]
    theta_n = (theta + np.pi) / (2 * np.pi)  # [0, 1)
    sector = np.clip((theta_n * n_sectors).astype(int), 0, n_sectors - 1)
    sums = np.zeros(n_sectors)
    for s in range(n_sectors):
        sums[s] = (I * (sector == s)).sum()
    return _gini(sums)


def _axis_projection_stats(intensity, origin_1d, pitch_um, axis, half_range_um=Y_HALF_RANGE_UM,
                           step_um=Y_STEP_UM):
    """axis=0: project along X (return Y-distribution). axis=1: project along Y
    (return X-distribution). origin_1d is the reference position on the retained
    axis, in pixels."""
    I_proj = intensity.sum(axis=axis)
    total = I_proj.sum()
    edges = np.arange(-half_range_um, half_range_um + step_um, step_um)
    n_bins = len(edges) - 1
    profile = np.full(n_bins, np.nan)
    nan_scalars = dict(mean_du=np.nan, sd_u=np.nan, skew_u=np.nan, kurt_u=np.nan,
                       apical_pct=np.nan, basal_pct=np.nan, apical_basal_ratio=np.nan,
                       gini=np.nan, entropy=np.nan)
    if total <= 0:
        return profile, edges, nan_scalars
    u_px = np.arange(len(I_proj))
    du_um = (u_px - origin_1d) * pitch_um  # signed distance from origin
    mean_u = float((I_proj * du_um).sum() / total)
    var_u = float((I_proj * (du_um - mean_u) ** 2).sum() / total)
    sd_u = float(np.sqrt(max(0.0, var_u)))
    if sd_u > 0:
        skew_u = float((I_proj * ((du_um - mean_u) / sd_u) ** 3).sum() / total)
        kurt_u = float((I_proj * ((du_um - mean_u) / sd_u) ** 4).sum() / total - 3)
    else:
        skew_u = kurt_u = np.nan
    # "apical" = signal on the negative side of origin (toward arch for Y-axis);
    # "basal" = on positive side. For X-axis, these are "left" vs "right" halves.
    apical = float(I_proj[du_um < 0].sum() / total * 100)
    basal = float(I_proj[du_um >= 0].sum() / total * 100)
    for i, (lo, hi) in enumerate(zip(edges[:-1], edges[1:])):
        mask = (du_um >= lo) & (du_um < hi)
        profile[i] = I_proj[mask].sum() / total * 100
    gini = _gini_of_bins(profile)
    entropy = _entropy_of_bins(profile)
    return profile, edges, dict(mean_du=mean_u, sd_u=sd_u, skew_u=skew_u, kurt_u=kurt_u,
                                apical_pct=apical, basal_pct=basal,
                                apical_basal_ratio=apical / max(basal, 1e-9),
                                gini=gini, entropy=entropy)


# -------- Per-cell processing --------

def process_cell(img_path, *, template_hat, template, save_projections=False) -> dict:
    key = tmb.cluster_key(img_path)
    offset = tmb.offset_overrides.get(key, [128, 128])
    roi = tmb.roi_overrides.get(key, None)

    img = nd2.imread(img_path, xarray=True)
    zsum = img.sum(axis=0)
    zmax = img.max(axis=0)

    sumproj = zsum.sel(C="640").to_numpy()[offset[0]:2048 + offset[0], offset[1]:2048 + offset[1]]
    sumproj_thresh = sumproj > skimage.filters.threshold_otsu(sumproj)
    max_coords = tmb.get_template_center(img, img_path, template_hat=template_hat,
                                         offset=offset, roi=roi)
    shifted_template = np.roll(template, (max_coords[0] - 1024, max_coords[1] - 1024), axis=(0, 1))
    score = np.sum(sumproj_thresh & shifted_template) / np.sum(shifted_template > 0)

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
    nuc_max = int(np.argmax([p.area for p in nuc_props])) + 1
    nuc_mask = nuc_label == nuc_max
    nuc_edt_px = distance_transform_edt(np.invert(nuc_mask))
    nuc_com = center_of_mass(nuc_mask)

    arch_px = np.zeros_like(nuc_mask)
    arch_px[np.round(contour[1083:1951, 0]).astype(int),
            np.round(contour[1083:1951, 1]).astype(int)] = True
    arch_edt_px = distance_transform_edt(np.invert(arch_px))

    pattern_mask_big = shifted_template[y_start:y_end, x_start:x_end] > 0
    pattern_com = center_of_mass(pattern_mask_big)

    pitch_um = img.metadata["metadata"].channels[0].volume.axesCalibration[0]
    nuc_d_um = nuc_edt_px * pitch_um
    arch_d_um = arch_edt_px * pitch_um

    # Nucleus seg quality (saved for sanity-check; cheap)
    nuc_region = skimage.measure.regionprops(nuc_mask.astype(int))[0]
    seg_quality = dict(
        nuc_area_um2=float(nuc_region.area * pitch_um ** 2),
        nuc_solidity=float(nuc_region.solidity),
        nuc_eccentricity=float(nuc_region.eccentricity),
        nuc_euler_number=int(nuc_region.euler_number),
        nuc_n_components=int(len(nuc_props)),
        nuc_largest_area_frac=float(nuc_props[nuc_max - 1].area / sum(p.area for p in nuc_props)),
    )

    metrics: dict = {
        "path": str(img_path),
        "template_matching_score": float(score),
        "lateral_pixel_pitch_um": pitch_um,
        **seg_quality,
    }

    for proj_name, proj_xr in [("zsum", cropped_zsum), ("maxip", cropped_zmax)]:
        mito_raw = proj_xr.sel(C="488").to_numpy().astype(np.float64)
        mito_bg, bg_thr = _bg_subtracted(mito_raw)
        metrics[f"{proj_name}_bg_threshold"] = bg_thr
        I = mito_bg  # crop mask only — drop pattern mask (confounded per discussion)
        total = float(I.sum())
        metrics[f"{proj_name}_total_signal"] = total
        if total <= 0:
            continue

        # ---- Mark-style existing metrics (for direct comparison) ----
        peri = (nuc_edt_px < 5 / pitch_um)
        per = (arch_edt_px <= 5 / pitch_um)
        metrics[f"{proj_name}_perinuclear_5um_pct"] = float((peri * I).sum() / total * 100)
        metrics[f"{proj_name}_peripheral_5um_pct"] = float((per * I).sum() / total * 100)
        # 5-bin coarse radial profile (Mark's existing)
        for label, lo, hi in [("0_2", 0, 2), ("2_5", 2, 5), ("5_10", 5, 10),
                              ("10_15", 10, 15), ("ge15", 15, None)]:
            if hi is None:
                m = nuc_d_um >= lo
            else:
                m = (nuc_d_um >= lo) & (nuc_d_um < hi)
            metrics[f"{proj_name}_radial_{label}um_pct"] = float((m * I).sum() / total * 100)
        # CoM offsets
        dy, dx, mag = _com_offset_um(I, nuc_com, pitch_um)
        metrics[f"{proj_name}_com_dy_um"] = dy
        metrics[f"{proj_name}_com_dx_um"] = dx
        metrics[f"{proj_name}_com_offset_um"] = mag
        dy_p, dx_p, mag_p = _com_offset_um(I, pattern_com, pitch_um)
        metrics[f"{proj_name}_com_vs_pattern_dy_um"] = dy_p
        metrics[f"{proj_name}_com_vs_pattern_offset_um"] = mag_p
        metrics[f"{proj_name}_apical_fraction_pct"] = _apical_fraction(I, pattern_com[0])
        metrics[f"{proj_name}_mean_dist_to_nucleus_um"] = _weighted_mean(nuc_d_um.ravel(), I.ravel())
        metrics[f"{proj_name}_mean_dist_to_arch_um"] = _weighted_mean(arch_d_um.ravel(), I.ravel())
        metrics[f"{proj_name}_pixel_gini"] = _gini(I)

        # ---- NEW: fine radial profile ----
        fine_pct, fine_vol, edges = _fine_radial_profile(I, nuc_d_um)
        for i, pct in enumerate(fine_pct):
            lo, hi = int(edges[i]), int(edges[i + 1])
            metrics[f"{proj_name}_fine_r_{lo:02d}_{hi:02d}um_pct"] = float(pct)

        # ---- NEW: radial moments ----
        rm = _radial_moments(I, nuc_d_um)
        for k, v in rm.items():
            metrics[f"{proj_name}_{k}_um" if k in ("mean_r", "sd_r") else f"{proj_name}_{k}"] = v

        # ---- NEW: radial quantiles ----
        rq = _radial_quantiles(I, nuc_d_um, qs=(0.10, 0.25, 0.50, 0.75, 0.90))
        for q, v in rq.items():
            metrics[f"{proj_name}_dist_to_nuc_{q}_um"] = v

        # ---- NEW: radial Gini / entropy / KS-vs-uniform (over the fine profile) ----
        metrics[f"{proj_name}_radial_gini_fine"] = _gini_of_bins(fine_pct)
        metrics[f"{proj_name}_radial_entropy_fine"] = _entropy_of_bins(fine_pct)
        metrics[f"{proj_name}_radial_ks_vs_uniform"] = _ks_vs_uniform(fine_pct, fine_vol)

        # ---- NEW: moment tensor around nucleus CoM + pattern CoM ----
        mt_nuc = _moment_tensor(I, nuc_com, pitch_um)
        for k, v in mt_nuc.items():
            metrics[f"{proj_name}_mt_nuc_{k}"] = v
        mt_pat = _moment_tensor(I, pattern_com, pitch_um)
        for k, v in mt_pat.items():
            metrics[f"{proj_name}_mt_pat_{k}"] = v

        # ---- NEW: angular-sector Gini ----
        metrics[f"{proj_name}_angular_gini_nuc"] = _angular_sector_gini(I, nuc_com)
        metrics[f"{proj_name}_angular_gini_pat"] = _angular_sector_gini(I, pattern_com)

        # ---- NEW: Y-axis projection (axis=1 collapses X, yielding Y distribution) ----
        y_prof, y_edges, y_scalars = _axis_projection_stats(I, pattern_com[0], pitch_um, axis=1)
        for k, v in y_scalars.items():
            metrics[f"{proj_name}_y_{k}"] = v
        for i, pct in enumerate(y_prof):
            lo = int(y_edges[i])
            metrics[f"{proj_name}_y_profile_{lo:+04d}um_pct"] = float(pct)

        # ---- NEW: X-axis projection (axis=0 collapses Y, yielding X distribution) ----
        x_prof, x_edges, x_scalars = _axis_projection_stats(I, pattern_com[1], pitch_um, axis=0)
        for k, v in x_scalars.items():
            metrics[f"{proj_name}_x_{k}"] = v
        for i, pct in enumerate(x_prof):
            lo = int(x_edges[i])
            metrics[f"{proj_name}_x_profile_{lo:+04d}um_pct"] = float(pct)

    return metrics


# -------- Orchestration --------

def iter_cells(root: pathlib.Path):
    skip = {"MaxIP", "MaxIPs", "Excluded_cells", "denoised"}
    for plate, well_prefix in TARGET_WELLS:
        plate_dir = root / plate
        if not plate_dir.exists():
            print(f"  MISSING plate dir: {plate_dir}", flush=True)
            continue
        for well_dir in plate_dir.iterdir():
            if not well_dir.is_dir() or not well_dir.name.startswith(well_prefix):
                continue
            # walk, skipping excluded subdirs
            for dirpath, dirnames, filenames in well_dir.walk():
                dirnames[:] = [d for d in dirnames if d not in skip]
                for fn in filenames:
                    if fn.endswith(".nd2") and fn.lower().startswith("cell"):
                        yield dirpath / fn


def run(root: pathlib.Path, out_root: pathlib.Path) -> int:
    out_root.mkdir(parents=True, exist_ok=True)
    template_hat = tmb.get_template_hat(1326)
    template = tmb.get_padded_template_at_width(1326)

    cells = list(iter_cells(root))
    print(f"[fig4b_pipeline] {len(cells)} cells to process under {root}", flush=True)

    records_by_well: dict[pathlib.Path, list[dict]] = {}
    for i, img_path in enumerate(cells, 1):
        well_dir = img_path.parent
        # per-cell checkpoint (SMB-drop tolerant)
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
        print(f"[fig4b_pipeline] wrote {out_csv} ({len(recs)} cells)", flush=True)

    all_recs = [r for recs in records_by_well.values() for r in recs]
    if all_recs:
        combined = out_root / "combined_raw.csv"
        pl.from_dicts(all_recs).write_csv(combined)
        print(f"[fig4b_pipeline] wrote {combined} ({len(all_recs)} cells total)", flush=True)
    return 0


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", default=str(pathlib.Path(tmb.DATA_ROOT)),
                    help="patterned_data root (defaults to template_matching_bulk.DATA_ROOT)")
    ap.add_argument("--out-root", default=str(DEFAULT_OUT))
    args = ap.parse_args()
    return run(pathlib.Path(args.root).resolve(), pathlib.Path(args.out_root).resolve())


if __name__ == "__main__":
    sys.exit(main())
