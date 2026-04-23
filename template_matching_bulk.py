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

def get_template_center(img, path, *, template_hat = None, offset=None, roi=None):
    if cluster_key(path) in coordinate_overrides_dict:
        return top_coordinate_overrides_to_template_center(str(path), offset=offset)
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

    cropped_proj_mitochondria = cropped_proj_img.sel(C="488")
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
        cropped_rgb = make_rgb(
           stretch01(cropped_proj_img.sel(C="488")),
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

    output = {
            "score": score,
            "mitochondria_sum": mitochondria_sum,
            "cropped_background_threshold": cropped_background_threshold,
            **dist_results
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
