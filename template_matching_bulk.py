import marimo as mo
import cairosvg
import skimage
import numpy as np
import pymupdf
import io
import nd2
import matplotlib.pyplot as plt
import os
import pathlib
import polars as pl
import sys
import traceback
import netCDF4
import argparse
from scipy.ndimage import distance_transform_edt
from matplotlib.backends.backend_pdf import PdfPages

# Root under which ND2 data lives. Mark's original cluster path is the default
# so his override-dict keys (which are absolute cluster paths) continue to
# match; override with MICROPATTERN_DATA_ROOT when running off-cluster.
CLUSTER_DATA_ROOT = "/groups/vale/valelab/_for_Mark/patterned_data"
DATA_ROOT = os.environ.get("MICROPATTERN_DATA_ROOT", CLUSTER_DATA_ROOT)

def cluster_key(img_path):
    """Translate an on-disk path to the cluster-absolute path Mark used as
    override-dict keys, so lookups work regardless of where the data is mounted.
    """
    p = pathlib.Path(img_path).resolve()
    try:
        rel = p.relative_to(pathlib.Path(DATA_ROOT).resolve())
        return str(pathlib.Path(CLUSTER_DATA_ROOT) / rel)
    except ValueError:
        return str(img_path)

def get_coordinate_overrides_dict():
    schema = {
        "path": pl.String,
        "x": pl.Int16,
        "y": pl.Int16
    }
    coordinate_overrides_df = pl.read_csv(
        "coordinate_overrides.csv",
        has_header=False,
        new_columns = ["path", "x", "y"],
        schema=schema
    )

    coordinate_overrides_dict = {}
    for row in coordinate_overrides_df.iter_rows():
        coordinate_overrides_dict[row[0]] = (row[1],row[2])

    return coordinate_overrides_dict

coordinate_overrides_dict = get_coordinate_overrides_dict()

def top_coordinate_overrides_to_template_center(path, *, offset = None):
    # Distance in pixels from the top of the template to the center of the template
    top_to_center = 385 # 1024 - np.argmax(template[:,1024])
    key = cluster_key(path)
    if offset is None:
        offset = offset_overrides.get(key, [128, 128])
    top_x, top_y = coordinate_overrides_dict[key]
    # return in the same order as max_match_template
    return top_y + top_to_center - offset[0], top_x - offset[1]

def find_override_key(path):
    """Return the override-dict key matching `path`, falling back to the
    denoised-counterpart key if the raw path isn't directly registered.
    Raw runs read `Cell N.nd2`, but overrides are often recorded against
    `denoised/Cell N - Denoised.nd2` (or `… Denoised2.nd2` on plate 12)."""
    raw_key = cluster_key(path)
    if raw_key in coordinate_overrides_dict:
        return raw_key
    p = pathlib.Path(path)
    for suffix in (" - Denoised.nd2", " - Denoised2.nd2"):
        candidate = p.parent / "denoised" / f"{p.stem}{suffix}"
        k = cluster_key(str(candidate))
        if k in coordinate_overrides_dict:
            return k
    return None


def get_template_center(img, path, *, template_hat = None, offset=None, roi=None):
    key = find_override_key(path)
    if key is not None:
        return top_coordinate_overrides_to_template_center(key, offset=offset)
    return max_match_template(img, template_hat = template_hat, offset=offset, roi=roi)

def get_template_at_width(width):
    file = pymupdf.open("single_pattern.ai")
    png_bytes = cairosvg.svg2png(file[0].get_svg_image(), output_width=width, output_height=width)
    with io.BytesIO() as buf:
        buf.write(png_bytes)
        template_img = skimage.io.imread(buf)
    return template_img


def get_padded_template_at_width(template_width, *, base_template=None):
    if base_template is None:
        base_template = get_template_at_width(template_width)[:,:,0]
    pad = (2048-template_width)//2
    template = np.pad(base_template,(pad, pad))
    return template

def get_template_hat(template_width):
    """
    get frequency space template at size 2048
    """
    template = get_padded_template_at_width(template_width)
    #pad = (2048-template_width)//2
    #template = np.pad(get_template_at_width(template_width)[:,:,0],(pad, pad))
    template_flipped = np.flip(template, axis=(0,1))
    template_hat = np.fft.fft2(template_flipped)
    return template_hat

offset_overrides = {
    "/groups/vale/valelab/_for_Mark/patterned_data/250731_patterned_plate_11_good/E05_250808_TRAK1_wt/Cell2.nd2": [256,128],
    "/groups/vale/valelab/_for_Mark/patterned_data/250731_patterned_plate_11_good/E05_250808_TRAK1_wt/Cell2 - Denoised.nd2": [256,128],
    "/groups/vale/valelab/_for_Mark/patterned_data/250731_patterned_plate_11_good/E05_250808_TRAK1_wt/Cell5.nd2": [204,128],
    "/groups/vale/valelab/_for_Mark/patterned_data/250731_patterned_plate_11_good/E05_250808_TRAK1_wt/Cell5 - Denoised.nd2": [204,128],

    "/groups/vale/valelab/_for_Mark/patterned_data/250731_patterned_plate_11_good/E05_250808_TRAK1_wt/Cell8.nd2": [256,128],
    "/groups/vale/valelab/_for_Mark/patterned_data/250731_patterned_plate_11_good/E05_250808_TRAK1_wt/Cell8 - Denoised.nd2": [256,128],
    "/groups/vale/valelab/_for_Mark/patterned_data/250612_patterned_plate_3/B06_250617_TRAK1_mDRH_dSp/Cell12.nd2": [64,128]
}

def get_img_template_sum_projection(img):
    channel_640 = img.sel(C="640")
    if channel_640.ndim > 3:
        # If there are two channel 640s, then use the first one
        # TODO: Make this more specific to the channels
        channel_640 = channel_640.isel(C=0)
    img_template_sum_projection = np.sum(img.sel(C="640"), axis=0)
    assert img_template_sum_projection.ndim == 2
    return img_template_sum_projection


def get_image_hat(img, offset=[128, 128]):
    img_template_sum_projection = get_img_template_sum_projection(img)
    img_template_sum_projection_norm = img_template_sum_projection / np.max(img_template_sum_projection)
    #img_template_sum_projection_hat = np.abs(np.fft.fft2(img_template_sum_projection))
    img_template_sum_projection_norm_2048 = img_template_sum_projection_norm[offset[0]:2048+offset[0],offset[1]:2048+offset[1]]
    img_template_sum_projection_norm_2048_hat = np.fft.fft2(img_template_sum_projection_norm_2048)
    return img_template_sum_projection_norm_2048_hat

"""
def get_image_hat(img):
    channel_640 = img.sel(C="640")
    if channel_640.ndim > 3:
        # If there are two channel 640s, then use the first one
        # TODO: Make this more specific to the channels
        channel_640 = channel_640.isel(C=0)
    img_template_sum_projection = np.sum(channel_640, axis=0)
    assert img_template_sum_projection.ndim == 2
    img_template_sum_projection_norm = img_template_sum_projection / np.max(img_template_sum_projection)
    #img_template_sum_projection_hat = np.abs(np.fft.fft2(img_template_sum_projection))
    img_template_sum_projection_norm_2048 = img_template_sum_projection_norm[:2048,:2048]
    img_template_sum_projection_norm_2048_hat = np.fft.fft2(img_template_sum_projection_norm_2048)
    return img_template_sum_projection_norm_2048_hat
"""

def match_template(img, *, template_hat = None, offset=None):
    if isinstance(img, str) or isinstance(img, pathlib.Path):
        img_path = img
        offset = offset_overrides.get(cluster_key(img_path), [128, 128])
        img = nd2.imread(img_path, xarray=True)
    if template_hat is None:
        template_hat = get_template_hat(1326)
    img_template_hat = get_image_hat(img, offset=offset)
    template_matching = np.fft.fftshift(np.real(np.fft.ifft2(template_hat * img_template_hat)))
    return template_matching

roi_overrides = {
    "/groups/vale/valelab/_for_Mark/patterned_data/250731_patterned_plate_11_good/F06_250811_TRAK1_mDRH_dSp/Cell3.nd2": [slice(None), slice(0,1200)]
}

def max_match_template(img, *, template_hat = None, offset=None, roi=None):
    template_matching = match_template(img, template_hat = template_hat, offset=offset)

    if roi is not None:
       template_matching = template_matching[*roi]

    max_idx = np.argmax(template_matching)
    out = np.unravel_index(max_idx, template_matching.shape)

    if roi is not None:
        # reshift based on start
        if roi[0].start is not None:
            out = (out[0] + roi[0].start, out[1])
        if roi[1].start is not None:
            out = (out[0], out[1] + roi[1].start)

    return out

def stretch01(img, *, min_percentile=0.1, max_percentile=99.9):
    _min = np.percentile(img, min_percentile)
    _max = np.percentile(img, max_percentile)
    return np.clip((img - _min)/(_max - _min), 0, 1)

def make_rgb(R, G, B):
    RGB = np.zeros([3, *R.shape[-2:]], dtype="float32")
    if R is not None:
        RGB[0,:,:] = R
    if G is not None:
        RGB[1,:,:] = G
    if B is not None:
        RGB[2,:,:] = B
    return np.permute_dims(RGB,(1,2,0))

def draw_scale_bar(pixel_length):
    plt.plot([800, 800+pixel_length], [900, 900], color="white")
    plt.text(790, 950, "5 μm", color="white")


# ----------------------------------------------------------------------
# Wedge-r polar metric — added 2026-04-26
# Coordinates below are in the cropped-image frame (1024x1024), where
# the pattern center sits at (512, 512) by construction. The wedge is a
# fixed cone from the pattern's bottom apex through the upper-left and
# upper-right tangent points — pegged to the rigid micropattern, so the
# geometry is identical for every cell.
# ----------------------------------------------------------------------
WEDGE_APEX = (896, 512)    # (y, x) — pattern bottom extreme
WEDGE_LEFT = (373, 281)    # (y, x) — upper-left tangent
WEDGE_RIGHT = (374, 742)   # (y, x) — upper-right tangent
WEDGE_R_MAX_UM = 60.0
WEDGE_R_STEP_UM = 1.0
WEDGE_N_BINS = int(WEDGE_R_MAX_UM / WEDGE_R_STEP_UM)

# Two 15-µm radial slabs straddling the TRAK1/TRAK2 isobestic point at
# r ≈ 36.8 µm on the 60mer comparison sheet, with a symmetric 4 µm
# exclusion gap on either side. The inner slab samples the
# apex/centrosomal side of the perinuclear region and the outer slab
# samples the rim — together they form a directional companion to the
# unsigned wedge-r KS metric. See analysis/HANDOFF_v4.md §3 for the
# isobestic-point derivation.
WEDGE_CENTROSOMAL_BINS = (18, 33)   # [lo, hi)  — % wedge intensity in this slab
WEDGE_PERIPHERAL_BINS  = (41, 56)   # [lo, hi)  — % wedge intensity in this slab

# Empirical reference CDF for "passive cytoplasmic fill": per-cell
# wedge-r CDFs averaged across the 60mer no-TRAK condition (n=13) from
# the v3 whole-dataset run (analysis/overnight_final_out).
_REF_CDF_60MER_NOTRAK = np.array([
    0.000103829845, 0.000469579903, 0.001069454900, 0.002137864065, 0.003420138193, 0.006024210367,
    0.009229198131, 0.012847223574, 0.017205047483, 0.022979989072, 0.029759357091, 0.037543832800,
    0.045123871418, 0.052820943994, 0.063116039466, 0.073638884012, 0.082496243332, 0.093136748467,
    0.106977810043, 0.124315424593, 0.146213142644, 0.166457629086, 0.185143607069, 0.201790127430,
    0.219382599816, 0.236713316688, 0.262022446018, 0.290991208900, 0.314743648317, 0.337914281272,
    0.362441917900, 0.389370616054, 0.414549696500, 0.443135468372, 0.470286422415, 0.498785208113,
    0.529507112973, 0.562214707185, 0.596735497686, 0.633601782715, 0.670488483934, 0.706088578479,
    0.740015134297, 0.774175247822, 0.806066807562, 0.837525495249, 0.864098811278, 0.888682845265,
    0.909734033801, 0.928625176076, 0.944508945154, 0.957637877330, 0.967540833814, 0.976631144914,
    0.984397267326, 0.991116542234, 0.995923098873, 0.998877404052, 0.999967215325, 1.000000000000,
])

_wedge_geometry_cache: dict = {}

def _get_wedge_geometry(shape, pitch_um):
    """Build (wedge_mask, r_um, bin_idx_r, in_r, vol_arc) once per
    (shape, pitch). The wedge is a polygonal cone from WEDGE_APEX through
    WEDGE_LEFT/RIGHT, sweeping the upper hemicircle (dy < 0). Cached
    across cells because the pattern is rigid."""
    key = (tuple(shape), round(float(pitch_um), 6))
    if key in _wedge_geometry_cache:
        return _wedge_geometry_cache[key]
    H, W = shape
    Y_idx, X_idx = np.mgrid[:H, :W]
    apex = WEDGE_APEX
    dy_um = (Y_idx - apex[0]) * pitch_um
    dx_um = (X_idx - apex[1]) * pitch_um
    r_um = np.hypot(dy_um, dx_um)
    ang = np.arctan2(Y_idx - apex[0], X_idx - apex[1])
    a_left = float(np.arctan2(WEDGE_LEFT[0] - apex[0], WEDGE_LEFT[1] - apex[1]))
    a_right = float(np.arctan2(WEDGE_RIGHT[0] - apex[0], WEDGE_RIGHT[1] - apex[1]))
    # Pick the arc containing the L/R midpoint as seen from the apex, so
    # the branch-cut handling is correct for any apex/tangent geometry —
    # not just the upward-opening config we ship with.
    mid_y = 0.5 * (WEDGE_LEFT[0] + WEDGE_RIGHT[0])
    mid_x = 0.5 * (WEDGE_LEFT[1] + WEDGE_RIGHT[1])
    a_mid = float(np.arctan2(mid_y - apex[0], mid_x - apex[1]))
    lo, hi = min(a_left, a_right), max(a_left, a_right)
    if lo <= a_mid <= hi:
        wedge_mask = (ang >= lo) & (ang <= hi)
    else:
        wedge_mask = (ang <= lo) | (ang >= hi)
    bin_idx_r = np.floor(r_um / WEDGE_R_STEP_UM).astype(int)
    in_r = (bin_idx_r >= 0) & (bin_idx_r < WEDGE_N_BINS) & wedge_mask
    vol_arc = np.bincount(bin_idx_r[in_r].ravel(),
                          minlength=WEDGE_N_BINS)[:WEDGE_N_BINS].astype(np.float64)
    out = (wedge_mask, r_um, bin_idx_r, in_r, vol_arc)
    _wedge_geometry_cache[key] = out
    return out


def wedge_r_profile(intensity, *, shape=None, pitch_um=None, wedge_geom=None):
    """Per-bin intensity within the wedge as a (WEDGE_N_BINS,) ndarray of
    percentages summing to 100, or all-NaN if no signal lies in the wedge
    radial window. Normalization uses only pixels with bin in [0, N_BINS),
    so pixels at r >= WEDGE_R_MAX_UM (which would otherwise leak into the
    denominator without contributing to any bin) are excluded."""
    if wedge_geom is None:
        wedge_geom = _get_wedge_geometry(shape, pitch_um)
    _, _, bin_idx_r, in_r, _ = wedge_geom
    I = np.asarray(intensity, dtype=np.float64)
    profile = np.bincount(bin_idx_r[in_r].ravel(),
                          weights=I[in_r].ravel(),
                          minlength=WEDGE_N_BINS)[:WEDGE_N_BINS]
    total = float(profile.sum())
    if total <= 0:
        return np.full(WEDGE_N_BINS, np.nan, dtype=float)
    return profile / total * 100.0


def wedge_r_cdf(profile_pct):
    """Convert a per-bin % profile into a normalized CDF (cumsum / sum)."""
    p = np.asarray(profile_pct, dtype=np.float64)
    p = np.where(np.isnan(p), 0.0, p)
    s = p.sum()
    if s <= 0:
        return np.full_like(p, np.nan)
    return np.cumsum(p) / s


def ks_vs_uniform(profile_pct, vol_arc):
    """KS distance between the per-cell wedge-r CDF and the analytical
    area-uniform CDF derived from the wedge's per-bin pixel volume.
    Returns a scalar in [0, 1]; NaN if the wedge has no signal."""
    p = np.asarray(profile_pct, dtype=np.float64)
    v = np.asarray(vol_arc, dtype=np.float64)
    if np.nansum(p) <= 0 or v.sum() <= 0:
        return float("nan")
    p = np.where(np.isnan(p), 0.0, p)
    cdf_obs = np.cumsum(p) / p.sum()
    cdf_uni = np.cumsum(v) / v.sum()
    return float(np.max(np.abs(cdf_obs - cdf_uni)))


def ks_vs_60mer_noTRAK(profile_pct):
    """KS distance between the per-cell wedge-r CDF and the empirical
    60mer no-TRAK reference CDF (passive cytoplasmic fill baseline)."""
    cdf = wedge_r_cdf(profile_pct)
    if not np.isfinite(cdf).all():
        return float("nan")
    return float(np.max(np.abs(cdf - _REF_CDF_60MER_NOTRAK)))


def score_template_match(img_path, *, template_hat = None, template = None):
    # Load image
    img = nd2.imread(img_path, xarray=True)

    # Get offset and ROI overrides
    key = cluster_key(img_path)
    offset = offset_overrides.get(key, [128, 128])
    roi = roi_overrides.get(key, None)
    print(f"{img_path = }")
    #print(f"{offset = }")
    #print(f"{roi = }")

    # Get channel 640, resolving duplicates
    sumproj = get_img_template_sum_projection(img)[offset[0]:2048+offset[0], offset[1]:2048+offset[1]]
    sumproj_threshold = skimage.filters.threshold_otsu(sumproj.to_numpy())
    sumproj_thresholded = sumproj > sumproj_threshold
 
    # Match
    max_coords = get_template_center(img, img_path, template_hat=template_hat, offset=offset, roi=roi)
    #max_coords = max_match_template(img, template_hat=template_hat, offset=offset, roi=roi)
    shifted_template = np.roll(template, (max_coords[0] - 1024, max_coords[1] - 1024), axis=(0,1))
    shifted_template_contour = skimage.measure.find_contours(shifted_template)
    score = np.sum(sumproj_thresholded & shifted_template)/(np.sum(shifted_template > 0))
    score = score.values.item()

    relative_path = pathlib.Path(img_path).resolve().relative_to(pathlib.Path(DATA_ROOT).resolve())
    proj_path = pathlib.Path("projections",*relative_path.parts).with_suffix(".nc")
    proj_path.parent.mkdir(parents=True, exist_ok=True)

    y_start, y_end = max_coords[0]-512+offset[0], max_coords[0]+512+offset[0]
    x_start, x_end = max_coords[1]-512+offset[1], max_coords[1]+512+offset[1]
    
    cropped_proj_img = np.sum(img[:,:,y_start:y_end, x_start:x_end], axis=0)
    if cropped_proj_img.size == 0:
         print(f"Warning: Empty cropped_proj_img for {img_path}. {img.shape=}, {y_start=}:{y_end}, {x_start=}:{x_end}")
         raise ValueError(f"Empty cropped projection for {img_path}")

    cropped_proj_img.to_netcdf(proj_path)
    cropped_template_contour = shifted_template_contour[0].copy()
    cropped_template_contour[:,0] -= (max_coords[0]-512)
    cropped_template_contour[:,1] -= (max_coords[1]-512)

    cropped_nuc_proj = cropped_proj_img.sel(C="405")
    if cropped_nuc_proj.size == 0:
        raise ValueError(f"Empty cropped nuclear projection for {img_path}")
    nuc_proj_threshold = skimage.filters.threshold_otsu(cropped_nuc_proj.to_numpy())
    cropped_nuc_mask = cropped_nuc_proj > nuc_proj_threshold
    cropped_nuc_label = skimage.measure.label(cropped_nuc_mask)
    cropped_nuc_props = skimage.measure.regionprops(cropped_nuc_label)
    cropped_nuc_max_area = np.argmax([p.area for p in cropped_nuc_props])
    cropped_nuc_mask = (cropped_nuc_label == cropped_nuc_max_area+1)
    cropped_nuc_edt = distance_transform_edt(np.invert(cropped_nuc_mask))

    top_arch_mask = np.zeros_like(cropped_nuc_mask)
    top_arch_mask[
        np.round(cropped_template_contour[1083:1951,0]).astype("int"),
        np.round(cropped_template_contour[1083:1951,1]).astype("int")
    ] = True
    cropped_arch_edt = distance_transform_edt(np.invert(top_arch_mask))

    acute_arch_mask = np.zeros_like(cropped_nuc_mask)
    acute_arch_mask[
        np.round(cropped_template_contour[1300:1734,0]).astype("int"),
        np.round(cropped_template_contour[1300:1734,1]).astype("int")
    ] = True
    cropped_acute_arch_edt = distance_transform_edt(np.invert(acute_arch_mask))

    metadata = img.metadata["metadata"]
    lateral_pixel_pitch = metadata.channels[0].volume.axesCalibration[0]
    distances_um = [1, 2, 3, 4, 5]

    # MaxIP across Z for the mitochondria channel (replaces Mark's z-sum,
    # which is denoiser-sensitive at ~9% per cell on the wedge-r metrics;
    # MaxIP keeps drift below 1.1% — see HANDOFF_v3.md §1).
    cropped_proj_mitochondria = img[:, :, y_start:y_end, x_start:x_end].sel(C="488").max(axis=0)
    cropped_proj_mitochondria_stretched = stretch01(cropped_proj_mitochondria)

    # Assume the left and right edges consist of background
    cropped_background = np.concatenate((cropped_proj_mitochondria_stretched[:,:128], cropped_proj_mitochondria_stretched[:,-128:]), axis=1)
    cropped_background_unstretched = np.concatenate((cropped_proj_mitochondria[:,:128], cropped_proj_mitochondria[:,-128:]), axis=1)
    left_percentile = np.percentile(cropped_proj_mitochondria_stretched[:,:128], 99.99)
    right_percentile = np.percentile(cropped_proj_mitochondria_stretched[:,-128:], 99.99)
    left_percentile_unstretched = np.percentile(cropped_proj_mitochondria[:,:128], 99.99)
    right_percentile_unstretched = np.percentile(cropped_proj_mitochondria[:,-128:], 99.99)
    if abs(left_percentile - right_percentile) > 0.1:
        cropped_background_threshold = min(left_percentile, right_percentile)
        cropped_background_threshold_unstretched = min(left_percentile_unstretched, right_percentile_unstretched)
    else:
        cropped_background_threshold = np.percentile(cropped_background, 99.99)
        cropped_background_threshold_unstretched = np.percentile(cropped_background_unstretched, 99.99)
    cropped_proj_mitochondria_bg_subtracted = cropped_proj_mitochondria_stretched - cropped_background_threshold
    cropped_proj_mitochondria_bg_subtracted_unstretched = cropped_proj_mitochondria - cropped_background_threshold_unstretched
    # Set negative values to 0
    cropped_proj_mitochondria_bg_subtracted = np.clip(cropped_proj_mitochondria_bg_subtracted, 0, None)
    cropped_proj_mitochondria_bg_subtracted_unstretched = np.clip(cropped_proj_mitochondria_bg_subtracted_unstretched, 0, None)

    proj_mitochondria_bg_subtracted_path = pathlib.Path("projections",*relative_path.parts).with_suffix("")
    proj_mitochondria_bg_subtracted_path = proj_mitochondria_bg_subtracted_path.with_stem(
        f"{proj_mitochondria_bg_subtracted_path.stem}_488_bg_subtracted"
    )
    proj_mitochondria_bg_subtracted_path = proj_mitochondria_bg_subtracted_path.with_suffix(".nc")
    cropped_proj_mitochondria_bg_subtracted_unstretched.to_netcdf(proj_mitochondria_bg_subtracted_path)

    mitochondria_sum = np.sum(cropped_proj_mitochondria_bg_subtracted)

    dist_results = {}
    for d in distances_um:
        d_pixels = d / lateral_pixel_pitch
        d_perinuclear_mask = cropped_nuc_edt < d_pixels
        d_peripheral_mask = (cropped_arch_edt <= cropped_nuc_edt) & np.invert(d_perinuclear_mask)
        d_peripheral_d_um_mask = (cropped_arch_edt <= d_pixels) & d_peripheral_mask
        d_peripheral_simple_mask = cropped_arch_edt <= d_pixels
        d_acute_peripheral_mask = (cropped_acute_arch_edt <= cropped_nuc_edt) & np.invert(d_perinuclear_mask) & (cropped_arch_edt <= d_pixels)

        d_perinuclear_sum = np.sum(d_perinuclear_mask * cropped_proj_mitochondria_bg_subtracted)
        d_peripheral_d_um_sum = np.sum(d_peripheral_d_um_mask * cropped_proj_mitochondria_bg_subtracted)
        d_peripheral_simple_sum = np.sum(d_peripheral_simple_mask * cropped_proj_mitochondria_bg_subtracted)
        d_acute_peripheral_sum = np.sum(d_acute_peripheral_mask * cropped_proj_mitochondria_bg_subtracted)

        dist_results[f"perinuclear_{d}um_sum"] = d_perinuclear_sum
        dist_results[f"peripheral_{d}um_sum"] = d_peripheral_d_um_sum
        dist_results[f"peripheral_{d}um_simple_sum"] = d_peripheral_simple_sum
        dist_results[f"acute_peripheral_{d}um_sum"] = d_acute_peripheral_sum
        
        # Percentages relative to (peripheral + perinuclear)
        dist_results[f"peripheral_{d}um_percent"] = d_peripheral_d_um_sum / (d_peripheral_d_um_sum + d_perinuclear_sum + 1e-10) * 100
        dist_results[f"peripheral_{d}um_simple_percent"] = d_peripheral_simple_sum / (d_peripheral_simple_sum + d_perinuclear_sum + 1e-10) * 100
        dist_results[f"acute_peripheral_{d}um_percent"] = d_acute_peripheral_sum / (d_acute_peripheral_sum + d_perinuclear_sum + 1e-10) * 100

    # Specifically retain 5um versions for plotting and existing metrics
    perinuclear_space_distance_pixels = 5 / lateral_pixel_pitch
    perinuclear_mask = cropped_nuc_edt < perinuclear_space_distance_pixels
    peripheral_mask = (cropped_arch_edt <= cropped_nuc_edt) & np.invert(perinuclear_mask)
    peripheral_5um_mask  = (cropped_arch_edt <= perinuclear_space_distance_pixels) & peripheral_mask
    peripheral_mask_extended = peripheral_mask & (cropped_arch_edt <= perinuclear_space_distance_pixels * 1.75)
    # Simple does peripheral mask does not depend on the nucleus position
    peripheral_5um_simple_mask = cropped_arch_edt <= perinuclear_space_distance_pixels

    acute_peripheral_mask = (cropped_acute_arch_edt <= cropped_nuc_edt) & np.invert(perinuclear_mask) & (cropped_arch_edt <= perinuclear_space_distance_pixels)

    perinuclear_mitochondria = perinuclear_mask * cropped_proj_mitochondria_bg_subtracted
    peripheral_5um_mitochondria = peripheral_5um_mask * cropped_proj_mitochondria_bg_subtracted

    peripheral_contour = skimage.measure.find_contours(peripheral_mask_extended)[0]
    peripheral_5um_contour = skimage.measure.find_contours(peripheral_5um_mask)[0]
    acute_peripheral_contour = skimage.measure.find_contours(acute_peripheral_mask)[0]
    perinuclear_contour = skimage.measure.find_contours(perinuclear_mask)[0]

    perinuclear_sum = dist_results["perinuclear_5um_sum"]
    peripheral_sum = np.sum(peripheral_mask_extended * cropped_proj_mitochondria_bg_subtracted)
    peripheral_5um_sum = dist_results["peripheral_5um_sum"]
    peripheral_5um_simple_sum = dist_results["peripheral_5um_simple_sum"]
    acute_peripheral_sum = dist_results["acute_peripheral_5um_sum"]

    pp_sum = perinuclear_sum + peripheral_sum
    peripheral_percent = peripheral_sum / (pp_sum + 1e-10) * 100
    peripheral_5um_percent = dist_results["peripheral_5um_percent"]
    peripheral_5um_simple_percent = dist_results["peripheral_5um_simple_percent"]
    acute_peripheral_percent = dist_results["acute_peripheral_5um_percent"]


    # draw contour around 0.5 um of nucleus
    cropped_nuclear_contour = get_nuclear_contour(cropped_nuc_edt < 0.5/lateral_pixel_pitch)

    # Plot figure to PDF
    # moved up
    # relative_path = pathlib.Path(img_path).relative_to("/groups/vale/valelab/_for_Mark/patterned_data")
    pdf_path = pathlib.Path("template_matching",*relative_path.parts).with_suffix(".pdf")
    pdf_path.parent.mkdir(parents=True, exist_ok=True)
    with PdfPages(pdf_path) as pdf:
        # Template matching figure
        fig = plt.figure()
        #plt.imshow(shifted_template)
        #plt.imshow(sumproj_thresholded) # , alpha=0.5)
        plt.imshow(stretch01(sumproj))
        plt.plot(shifted_template_contour[0][:,1], shifted_template_contour[0][:,0], color="black")
        plt.scatter(max_coords[1], max_coords[0])
        plt.annotate(text="{:.3f}%".format(score*100), xy=(max_coords[1] + 100, max_coords[0]), color="yellow")
        plt.title(img_path, loc="right")
        pdf.savefig()
        plt.close()

        fig = plt.figure()
        # 488 viz uses the same MaxIP that the metrics consume so the PDF
        # image matches the perinuclear/peripheral/wedge numbers below.
        # 405 keeps Z-sum (Mark's nuclear-segmentation projection).
        cropped_rgb = make_rgb(
           stretch01(cropped_proj_mitochondria),
           None,
           stretch01(cropped_proj_img.sel(C="405"))
        )
        plt.imshow(cropped_rgb)
        plt.plot(cropped_template_contour[:,1], cropped_template_contour[:,0], color="white")
        # Top Arch
        plt.plot(cropped_template_contour[1083:1951,1], cropped_template_contour[1083:1951,0], color="magenta")
        # Top Left point
        plt.scatter(cropped_template_contour[1083,1], cropped_template_contour[1083,0]+60, color="white")
        # Top Right point
        plt.scatter(cropped_template_contour[1951,1], cropped_template_contour[1951,0]+60, color="white")
        # Bottom Middle point
        plt.scatter(cropped_template_contour[12,1], cropped_template_contour[12,0], color="white")
        # Bottom Left point
        plt.scatter(cropped_template_contour[1083,1], cropped_template_contour[12,0], color="white")
        # Bottom Right point
        plt.scatter(cropped_template_contour[1951,1], cropped_template_contour[12,0], color="white")
        draw_scale_bar(perinuclear_space_distance_pixels)
        pdf.savefig()
        plt.close()

        fig = plt.figure()
        plt.imshow(make_rgb(
            stretch01(-cropped_arch_edt),
            stretch01(-cropped_nuc_edt),
            stretch01(-cropped_arch_edt)
        ))
        plt.plot(cropped_template_contour[1083:1951,1],cropped_template_contour[1083:1951,0], color="black")
        draw_scale_bar(perinuclear_space_distance_pixels)
        pdf.savefig()
        plt.close()

        fig = plt.figure()
        plt.imshow(cropped_arch_edt <= cropped_nuc_edt)
        plt.plot(cropped_template_contour[:,1], cropped_template_contour[:,0], color="white")
        plt.plot(cropped_template_contour[1083:1951,1],cropped_template_contour[1083:1951,0], color="magenta", alpha=0.5)
        plt.plot(cropped_nuclear_contour[:,1], cropped_nuclear_contour[:,0], color="blue", alpha=0.5)
        draw_scale_bar(perinuclear_space_distance_pixels)
        pdf.savefig()
        plt.close()

        fig = plt.figure()
        plt.imshow(cropped_rgb)
        plt.plot(peripheral_contour[:,1],  peripheral_contour[:,0], color="yellow", linestyle="dotted")
        plt.plot(peripheral_5um_contour[:,1],  peripheral_5um_contour[:,0], color="yellow")
        plt.plot(acute_peripheral_contour[:,1],  acute_peripheral_contour[:,0], color="yellow", linestyle="dotted")
        plt.plot(cropped_template_contour[1083:1951,1],cropped_template_contour[1083:1951,0], color="magenta", alpha=0.5)
        plt.plot(cropped_template_contour[1300:1734,1],cropped_template_contour[1300:1734,0], color="cyan", alpha=0.5)
        plt.plot(perinuclear_contour[:,1], perinuclear_contour[:,0], color="blue")
        pdf.savefig()
        plt.close()

        fig = plt.figure()
        # Draw peripheral mitochondria as yellow
        plt.imshow(make_rgb(
            peripheral_5um_mitochondria,
            peripheral_5um_mitochondria,
            perinuclear_mitochondria
        ))
        plt.plot(cropped_template_contour[1083:1951,1],cropped_template_contour[1083:1951,0], color="white", alpha=0.5)
        plt.plot(cropped_template_contour[1300:1734,1],cropped_template_contour[1300:1734,0], color="cyan", alpha=0.5)
        plt.plot(cropped_nuclear_contour[:,1], cropped_nuclear_contour[:,0], color="white", alpha=0.5)
        print(f"Debug: {peripheral_5um_percent= }, {peripheral_5um_sum= }, {mitochondria_sum= }, {perinuclear_sum= }")
        plt.title("P5um/(P5um+N): {:.1f}%, P5um/Crop: {:.1f}%, N/Crop: {:.1f}%".format(peripheral_5um_percent, peripheral_5um_sum/mitochondria_sum*100, perinuclear_sum/mitochondria_sum*100))
        draw_scale_bar(perinuclear_space_distance_pixels)
        pdf.savefig()
        plt.close()

    # --- Wedge-r KS metric (added 2026-04-26 — see HANDOFF_v3.md §1).
    # The wedge is fixed by the rigid micropattern (same for every cell);
    # we evaluate the intensity-weighted CDF along its radial axis, then
    # report KS distance to (i) an analytical area-uniform sector and
    # (ii) the empirical 60mer no-TRAK reference CDF.
    wedge_geom = _get_wedge_geometry(
        cropped_proj_mitochondria_bg_subtracted.shape, lateral_pixel_pitch)
    wedge_profile = wedge_r_profile(
        cropped_proj_mitochondria_bg_subtracted, wedge_geom=wedge_geom)
    wedge_results = {
        "wedge_r_ks_vs_uniform": ks_vs_uniform(wedge_profile, wedge_geom[4]),
        "wedge_r_ks_vs_60merNoTRAK": ks_vs_60mer_noTRAK(wedge_profile),
        "wedge_r_centrosomal_18_33um_pct":
            float(wedge_profile[slice(*WEDGE_CENTROSOMAL_BINS)].sum()),
        "wedge_r_peripheral_41_56um_pct":
            float(wedge_profile[slice(*WEDGE_PERIPHERAL_BINS)].sum()),
    }
    for i in range(WEDGE_N_BINS):
        v = wedge_profile[i]
        wedge_results[f"wedge_r_{i:02d}_{i+1:02d}um_pct"] = (
            float(v) if np.isfinite(v) else float("nan"))

    output = {
            "score": score,
            "mitochondria_sum": mitochondria_sum,
            "cropped_background_threshold": cropped_background_threshold,
            **dist_results,
            **wedge_results,
    }

    return output

def get_nuclear_contour(nuclear_mask):
    nuclear_contours = skimage.measure.find_contours(nuclear_mask)
    nuclear_contour_index = np.argmax([len(contour) for contour in nuclear_contours])
    return nuclear_contours[nuclear_contour_index]

def main(root_path, keep_sums=False, only_simple=True, remove_acute=True, only_total=True):
    pl.Config.set_tbl_cell_alignment("RIGHT")

    #Set up template
    template_hat = get_template_hat(1326)
    template = get_padded_template_at_width(1326)

    print(f"Scanning {root_path}")

    # img_path = "/groups/vale/valelab/_for_Mark/patterned_data/250521_patterned_plate_1/B06_250528_TRAK1-wt/Cell8.nd2"
    for (dirpath, dirnames, filenames) in pathlib.Path(root_path).walk():
        if dirpath.parts[-1] in ["MaxIP", "MaxIPs", "Excluded_cells"]:
            continue

        records = []
        relative_path = dirpath.resolve().relative_to(pathlib.Path(DATA_ROOT).resolve())
        csv_path = pathlib.Path("template_matching", *relative_path.parts, "template_matching.csv")
        xlsx_path = pathlib.Path("template_matching", *relative_path.parts, "template_matching.xlsx")
        print(csv_path)
        print(xlsx_path)

        for filename in filenames:
            if filename.endswith(".nd2") and (filename.startswith("Cell") or filename.startswith("cell")):
                img_path = dirpath / filename
                try:
                    print(img_path)
                    output = score_template_match(img_path, template_hat=template_hat, template=template)
                    output["path"] = str(img_path)
                    output["template_matching_score"] = output.pop("score")
                    records.append(output)
                except Exception as e:
                    print(f"An error occurred with {img_path}: {e}")
                    traceback.print_exc()
                    records.append({"path": str(img_path)})

        if not records:
            continue

        print("Writing data frame")
        score_df = pl.from_dicts(records)

        # Add percent total columns for each distance
        for d in [1, 2, 3, 4, 5]:
            cols_to_add = []
            if f"peripheral_{d}um_sum" in score_df.columns:
                cols_to_add.append((pl.col(f"peripheral_{d}um_sum") / pl.col("mitochondria_sum") * 100).alias(f"peripheral_{d}um_percent_total"))
            if f"peripheral_{d}um_simple_sum" in score_df.columns:
                cols_to_add.append((pl.col(f"peripheral_{d}um_simple_sum") / pl.col("mitochondria_sum") * 100).alias(f"peripheral_{d}um_simple_percent_total"))
            if f"perinuclear_{d}um_sum" in score_df.columns:
                cols_to_add.append((pl.col(f"perinuclear_{d}um_sum") / pl.col("mitochondria_sum") * 100).alias(f"perinuclear_{d}um_percent_total"))
            
            if cols_to_add:
                score_df = score_df.with_columns(cols_to_add)

        if not keep_sums:
            sum_cols = [col for col in score_df.columns if col.endswith("_sum")]
            score_df = score_df.drop(sum_cols)

        if remove_acute:
            acute_cols = [col for col in score_df.columns if "acute_" in col]
            score_df = score_df.drop(acute_cols)

        if only_simple:
            # Drop non-simple peripheral columns, but keep metadata and perinuclear reference
            cols_to_drop = []
            for col in score_df.columns:
                if "peripheral_" in col and "_simple" not in col:
                    cols_to_drop.append(col)
            score_df = score_df.drop(cols_to_drop)

        if only_total:
            # Keep only columns with "_total" or essential metadata
            metadata_cols = ["path", "template_matching_score", "cropped_background_threshold"]
            total_cols = [col for col in score_df.columns if "_total" in col]
            # Ensure we only try to select columns that exist
            cols_to_select = [col for col in metadata_cols if col in score_df.columns] + total_cols
            score_df = score_df.select(cols_to_select)

        csv_path.parent.mkdir(parents=True, exist_ok=True)
        score_df.write_csv(csv_path)
        score_df.write_excel(xlsx_path)

        def right_abbreviate(text: str, max_len: int = 20) -> str:
            if len(text) > max_len:
                return "..." + text[-(max_len - 3):] # Keep the end of the string
            return text

        score_df_abbreviated = score_df.with_columns(
            pl.col("path").map_elements(right_abbreviate).alias("path")
        )
        print(score_df_abbreviated)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Analyze micropatterned cell images.")
    parser.add_argument("root_path", help="Path to the root directory containing .nd2 files.")
    parser.add_argument("--keep-sums", action="store_true", help="Keep the _sum columns in the output (default: False).")
    parser.add_argument("--include-complex", action="store_false", dest="only_simple", help="Include non-simple peripheral measurements (default: simple only).")
    parser.add_argument("--include-acute", action="store_false", dest="remove_acute", help="Include acute peripheral measurements (default: removed).")
    parser.add_argument("--include-all-percents", action="store_false", dest="only_total", help="Include all percentage measurements (default: only _total).")
    parser.set_defaults(only_simple=True, remove_acute=True, only_total=True)
    args = parser.parse_args()
    main(args.root_path, keep_sums=args.keep_sums, only_simple=args.only_simple, remove_acute=args.remove_acute, only_total=args.only_total)
