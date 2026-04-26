"""Extended-metric rerun for Fig S11 D / E / F (and also Fig 4E).

Processes the sheets Mark annotates in his Prism files for the three
reviewer-flagged panels beyond 4B:
  - TRAK1 helix muts  (Fig 4C / S11 D)
  - TRAK2 helix muts  (Fig 4D / S11 E)
  - MAPK9 siRNA       (Fig 4E / S11 F)

Inherits every metric from fig4b_pipeline.process_cell, plus new
wedge-restricted radial metrics suggested by the user: instead of an
angle-isotropic radial profile, compute radial statistics only within the
angular wedge subtended by the cell's actual extent (defined by the three
outermost points of the pattern mask — bottom, left, right).

The wedge is expressed as an angular mask around the NUCLEUS CoM. Pixels
outside the wedge don't contribute to the wedge-restricted radial metrics.

Output: replication/overnight_s11_out/ with the same layout as
overnight_fig4b_out (per-well JSON checkpoints, per-well metrics.csv, and
a combined_raw.csv at the top level).
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
from replication.metric_pipeline import (  # noqa: E402
    _bg_subtracted, _gini, _com_offset_um, _apical_fraction,
    _weighted_mean,
)
from replication.fig4b_pipeline import (  # noqa: E402
    _fine_radial_profile, _radial_moments, _radial_quantiles,
    _gini_of_bins, _entropy_of_bins, _ks_vs_uniform,
    _moment_tensor, _angular_sector_gini, _axis_projection_stats,
    FINE_R_MAX_UM, FINE_R_STEP_UM, Y_HALF_RANGE_UM, Y_STEP_UM, N_SECTORS,
)

DEFAULT_OUT = REPO / "replication" / "overnight_s11_out"

# Wells for S11 D / E / F. Discovered by inspecting combined.csv with:
#   awk -F, '$sheet!="TRAK isoform (mito)" { print $plate "/" $well "\t" $sheet "\t" $condition }'
# Populated programmatically at runtime from the existing combined.csv.


def discover_wells(old_combined: pathlib.Path, sheets: list[str]) -> list[tuple[str, str]]:
    """Return unique (plate, well-prefix) pairs for the given sheets."""
    df = pl.read_csv(old_combined).filter(pl.col("sheet").is_in(sheets))
    pairs = set()
    for plate, well in df.select(["plate", "well"]).unique().iter_rows():
        pairs.add((plate, f"{well}_"))
    return sorted(pairs)


# -------- New: wedge-restricted radial metrics --------

def _pattern_extremes(pattern_mask: np.ndarray) -> dict:
    """Return the outermost bottom/left/right pixel coords of the pattern mask
    (non-zero pixels). Coords are (y, x) in pixels."""
    ys, xs = np.where(pattern_mask)
    if ys.size == 0:
        return {"bottom": None, "left": None, "right": None, "top": None}
    # "bottom" = max y (image coords; y increases downward)
    ib = int(np.argmax(ys))
    il = int(np.argmin(xs))
    ir = int(np.argmax(xs))
    it = int(np.argmin(ys))
    return {
        "bottom": (int(ys[ib]), int(xs[ib])),
        "left":   (int(ys[il]), int(xs[il])),
        "right":  (int(ys[ir]), int(xs[ir])),
        "top":    (int(ys[it]), int(xs[it])),
    }


def _wedge_mask(shape: tuple, origin_yx: tuple, outer_points: dict) -> np.ndarray:
    """Build an angular wedge mask. The wedge is the minor-arc span from
    `left` through `bottom` to `right` around origin_yx. (Pixels whose angle
    from origin falls in that arc are included.) Pixels above the arc are
    excluded.

    Wedge is defined so that 'below' and 'lateral' regions — where the cell
    fan lives — are included, and the region directly above the nucleus (toward
    the arch) is excluded only if the 'top' of the pattern is used as the
    upper boundary. Here we use left/bottom/right to form a ~270° wedge
    excluding the 90° directly above the nucleus.
    """
    H, W = shape
    Y, X = np.mgrid[:H, :W]
    dy = Y - origin_yx[0]
    dx = X - origin_yx[1]
    # Angle measured counterclockwise from +x axis; y increases downward in image coords
    # so angle convention: atan2(dy, dx). Below origin = positive dy = positive angle.
    ang = np.arctan2(dy, dx)  # in [-pi, pi]

    def pt_angle(pt):
        return float(np.arctan2(pt[0] - origin_yx[0], pt[1] - origin_yx[1]))

    a_left = pt_angle(outer_points["left"])
    a_bot = pt_angle(outer_points["bottom"])
    a_right = pt_angle(outer_points["right"])

    # We want the wedge going counterclockwise from a_right → a_bot → a_left
    # (i.e. the wedge that INCLUDES bottom). Since y is flipped, "bottom"
    # has the largest positive angle among the three.
    # The wedge is: angles in the arc from a_right (smallest) → a_left (largest),
    # passing through a_bot. Equivalently, include pixels whose angle is between
    # min(a_right, a_left) and max(a_right, a_left) along the arc containing a_bot.
    # Simplest: build a mask using "angle from origin to any of the three points
    # in a range that contains the bottom point's angle".
    lo = min(a_right, a_left)
    hi = max(a_right, a_left)
    # If the bottom angle is within [lo, hi], include [lo, hi]. Otherwise include
    # the complement (wrap-around case).
    if lo <= a_bot <= hi:
        mask = (ang >= lo) & (ang <= hi)
    else:
        mask = (ang <= lo) | (ang >= hi)
    return mask


def _wedge_restricted_radial(intensity, nuc_edt_um, wedge_mask,
                             max_r_um=FINE_R_MAX_UM, step_um=FINE_R_STEP_UM):
    """Radial metrics restricted to the wedge."""
    I = intensity * wedge_mask
    total = I.sum()
    if total <= 0:
        return None
    # Fine radial profile within wedge
    edges = np.arange(0, max_r_um + step_um, step_um)
    n_bins = len(edges) - 1
    pct = np.full(n_bins, np.nan)
    vol = np.zeros(n_bins)
    for i, (lo, hi) in enumerate(zip(edges[:-1], edges[1:])):
        m = (nuc_edt_um >= lo) & (nuc_edt_um < hi) & wedge_mask
        pct[i] = (m * intensity).sum() / total * 100
        vol[i] = m.sum()
    # Moments
    r = nuc_edt_um.ravel()
    w = I.ravel()
    total_in = float(w.sum())
    mean_r = float((w * r).sum() / total_in) if total_in > 0 else np.nan
    var_r = float((w * (r - mean_r) ** 2).sum() / total_in) if total_in > 0 else np.nan
    sd_r = float(np.sqrt(max(0.0, var_r))) if np.isfinite(var_r) else np.nan
    if sd_r > 0:
        skew_r = float((w * ((r - mean_r) / sd_r) ** 3).sum() / total_in)
        kurt_r = float((w * ((r - mean_r) / sd_r) ** 4).sum() / total_in - 3)
    else:
        skew_r = kurt_r = np.nan
    # Gini / entropy / KS
    gini = _gini_of_bins(pct)
    entropy = _entropy_of_bins(pct)
    ks = _ks_vs_uniform(pct, vol)
    # Fraction of total cell signal captured by the wedge
    wedge_frac = float(total / max(intensity.sum(), 1e-9) * 100)

    out = {
        "wedge_frac_pct": wedge_frac,
        "wedge_mean_r_um": mean_r, "wedge_sd_r_um": sd_r,
        "wedge_skew_r": skew_r, "wedge_kurt_r": kurt_r,
        "wedge_radial_gini": gini, "wedge_radial_entropy": entropy,
        "wedge_radial_ks_vs_uniform": ks,
    }
    for i, p in enumerate(pct):
        lo, hi = int(edges[i]), int(edges[i + 1])
        out[f"wedge_fine_r_{lo:02d}_{hi:02d}um_pct"] = float(p)
    return out


# -------- Per-cell routine (extends fig4b_pipeline.process_cell) --------

def process_cell(img_path, *, template_hat, template) -> dict:
    key = tmb.cluster_key(img_path)
    offset = tmb.offset_overrides.get(key, [128, 128])
    roi = tmb.roi_overrides.get(key, None)

    img = nd2.imread(img_path, xarray=True)
    zsum = img.sum(axis=0)
    zmax = img.max(axis=0)

    sumproj = zsum.sel(C="640").to_numpy()[offset[0]:2048 + offset[0], offset[1]:2048 + offset[1]]
    _ = sumproj > skimage.filters.threshold_otsu(sumproj)
    max_coords = tmb.get_template_center(img, img_path, template_hat=template_hat,
                                         offset=offset, roi=roi)
    shifted_template = np.roll(template, (max_coords[0] - 1024, max_coords[1] - 1024), axis=(0, 1))

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

    nuc_region = skimage.measure.regionprops(nuc_mask.astype(int))[0]
    seg_quality = dict(
        nuc_area_um2=float(nuc_region.area * pitch_um ** 2),
        nuc_solidity=float(nuc_region.solidity),
        nuc_eccentricity=float(nuc_region.eccentricity),
        nuc_euler_number=int(nuc_region.euler_number),
        nuc_n_components=int(len(nuc_props)),
        nuc_largest_area_frac=float(nuc_props[nuc_max - 1].area / sum(p.area for p in nuc_props)),
    )

    # Pattern extremes and wedge mask (for wedge-restricted metrics)
    ext = _pattern_extremes(pattern_mask_big)
    wedge = _wedge_mask(pattern_mask_big.shape, nuc_com, ext)
    # Save pattern extreme coords in µm from nuc CoM for QC
    def _rel_um(pt):
        if pt is None:
            return (np.nan, np.nan)
        return (float((pt[0] - nuc_com[0]) * pitch_um),
                float((pt[1] - nuc_com[1]) * pitch_um))

    metrics: dict = {
        "path": str(img_path),
        "lateral_pixel_pitch_um": pitch_um,
        **seg_quality,
        "wedge_px_fraction": float(wedge.sum() / wedge.size),
    }
    for name in ("bottom", "left", "right", "top"):
        dy_um, dx_um = _rel_um(ext[name])
        metrics[f"pattern_{name}_dy_um_from_nuc"] = dy_um
        metrics[f"pattern_{name}_dx_um_from_nuc"] = dx_um

    for proj_name, proj_xr in [("zsum", cropped_zsum), ("maxip", cropped_zmax)]:
        mito_raw = proj_xr.sel(C="488").to_numpy().astype(np.float64)
        mito_bg, bg_thr = _bg_subtracted(mito_raw)
        metrics[f"{proj_name}_bg_threshold"] = bg_thr
        I = mito_bg
        total = float(I.sum())
        metrics[f"{proj_name}_total_signal"] = total
        if total <= 0:
            continue

        # Mark-style existing metrics
        peri = (nuc_edt_px < 5 / pitch_um)
        per = (arch_edt_px <= 5 / pitch_um)
        metrics[f"{proj_name}_perinuclear_5um_pct"] = float((peri * I).sum() / total * 100)
        metrics[f"{proj_name}_peripheral_5um_pct"] = float((per * I).sum() / total * 100)
        for label, lo, hi in [("0_2", 0, 2), ("2_5", 2, 5), ("5_10", 5, 10),
                              ("10_15", 10, 15), ("ge15", 15, None)]:
            if hi is None:
                m = nuc_d_um >= lo
            else:
                m = (nuc_d_um >= lo) & (nuc_d_um < hi)
            metrics[f"{proj_name}_radial_{label}um_pct"] = float((m * I).sum() / total * 100)
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

        # Extended metrics (Fig 4B set)
        fine_pct, fine_vol, edges = _fine_radial_profile(I, nuc_d_um)
        for i, pct in enumerate(fine_pct):
            lo, hi = int(edges[i]), int(edges[i + 1])
            metrics[f"{proj_name}_fine_r_{lo:02d}_{hi:02d}um_pct"] = float(pct)
        rm = _radial_moments(I, nuc_d_um)
        for k, v in rm.items():
            metrics[f"{proj_name}_{k}_um" if k in ("mean_r", "sd_r") else f"{proj_name}_{k}"] = v
        rq = _radial_quantiles(I, nuc_d_um, qs=(0.10, 0.25, 0.50, 0.75, 0.90))
        for q, v in rq.items():
            metrics[f"{proj_name}_dist_to_nuc_{q}_um"] = v
        metrics[f"{proj_name}_radial_gini_fine"] = _gini_of_bins(fine_pct)
        metrics[f"{proj_name}_radial_entropy_fine"] = _entropy_of_bins(fine_pct)
        metrics[f"{proj_name}_radial_ks_vs_uniform"] = _ks_vs_uniform(fine_pct, fine_vol)
        mt_nuc = _moment_tensor(I, nuc_com, pitch_um)
        for k, v in mt_nuc.items():
            metrics[f"{proj_name}_mt_nuc_{k}"] = v
        mt_pat = _moment_tensor(I, pattern_com, pitch_um)
        for k, v in mt_pat.items():
            metrics[f"{proj_name}_mt_pat_{k}"] = v
        metrics[f"{proj_name}_angular_gini_nuc"] = _angular_sector_gini(I, nuc_com)
        metrics[f"{proj_name}_angular_gini_pat"] = _angular_sector_gini(I, pattern_com)
        y_prof, y_edges, y_scalars = _axis_projection_stats(I, pattern_com[0], pitch_um, axis=1)
        for k, v in y_scalars.items():
            metrics[f"{proj_name}_y_{k}"] = v
        for i, pct in enumerate(y_prof):
            lo = int(y_edges[i])
            metrics[f"{proj_name}_y_profile_{lo:+04d}um_pct"] = float(pct)
        x_prof, x_edges, x_scalars = _axis_projection_stats(I, pattern_com[1], pitch_um, axis=0)
        for k, v in x_scalars.items():
            metrics[f"{proj_name}_x_{k}"] = v
        for i, pct in enumerate(x_prof):
            lo = int(x_edges[i])
            metrics[f"{proj_name}_x_profile_{lo:+04d}um_pct"] = float(pct)

        # NEW: wedge-restricted radial metrics
        wres = _wedge_restricted_radial(I, nuc_d_um, wedge)
        if wres is not None:
            for k, v in wres.items():
                metrics[f"{proj_name}_{k}"] = v

    return metrics


# -------- Orchestration --------

def iter_cells(root: pathlib.Path, target_wells: list[tuple[str, str]]):
    skip = {"MaxIP", "MaxIPs", "Excluded_cells", "denoised"}
    for plate, well_prefix in target_wells:
        plate_dir = root / plate
        if not plate_dir.exists():
            continue
        for well_dir in plate_dir.iterdir():
            if not well_dir.is_dir() or not well_dir.name.startswith(well_prefix):
                continue
            for dirpath, dirnames, filenames in well_dir.walk():
                dirnames[:] = [d for d in dirnames if d not in skip]
                for fn in filenames:
                    if fn.endswith(".nd2") and fn.lower().startswith("cell"):
                        yield dirpath / fn


def run(root: pathlib.Path, out_root: pathlib.Path, sheets: list[str]) -> int:
    out_root.mkdir(parents=True, exist_ok=True)
    template_hat = tmb.get_template_hat(1326)
    template = tmb.get_padded_template_at_width(1326)

    old_combined = REPO / "replication" / "overnight_out" / "combined.csv"
    target_wells = discover_wells(old_combined, sheets)
    print(f"[s11_pipeline] target sheets: {sheets}")
    print(f"[s11_pipeline] {len(target_wells)} wells across {len({p for p, _ in target_wells})} plates",
          flush=True)

    cells = list(iter_cells(root, target_wells))
    print(f"[s11_pipeline] {len(cells)} cells to process under {root}", flush=True)

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
        print(f"[s11_pipeline] wrote {out_csv} ({len(recs)} cells)", flush=True)

    all_recs = [r for recs in records_by_well.values() for r in recs]
    if all_recs:
        combined = out_root / "combined_raw.csv"
        pl.from_dicts(all_recs).write_csv(combined)
        print(f"[s11_pipeline] wrote {combined} ({len(all_recs)} cells total)", flush=True)
    return 0


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", default=str(pathlib.Path(tmb.DATA_ROOT)))
    ap.add_argument("--out-root", default=str(DEFAULT_OUT))
    ap.add_argument("--sheets", nargs="+",
                    default=["TRAK1 helix muts", "TRAK2 helix muts", "MAPK9 siRNA"])
    args = ap.parse_args()
    return run(pathlib.Path(args.root).resolve(), pathlib.Path(args.out_root).resolve(),
               args.sheets)


if __name__ == "__main__":
    sys.exit(main())
