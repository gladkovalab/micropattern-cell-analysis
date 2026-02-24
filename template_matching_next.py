import marimo

__generated_with = "0.17.8"
app = marimo.App(width="medium")


@app.cell
def _():
    import marimo as mo
    import cairosvg
    import skimage
    import numpy as np
    import pymupdf
    import io
    return cairosvg, io, np, pymupdf, skimage


@app.cell
def _():
    import nd2
    return (nd2,)


@app.cell
def _():
    import matplotlib.pyplot as plt
    return (plt,)


@app.cell
def _():
    import pathlib
    return (pathlib,)


@app.cell
def _():
    import polars as pl
    return (pl,)


@app.cell
def _():
    import xarray
    import netCDF4
    return (xarray,)


@app.cell
def _():
    from scipy.ndimage import distance_transform_edt
    return (distance_transform_edt,)


@app.cell
def _(cairosvg, io, pymupdf, skimage):
    def get_template_at_width(width):
        file = pymupdf.open("single_pattern.ai")
        png_bytes = cairosvg.svg2png(file[0].get_svg_image(), output_width=width, output_height=width)
        with io.BytesIO() as buf:
            buf.write(png_bytes)
            template_img = skimage.io.imread(buf)
        return template_img
    return (get_template_at_width,)


@app.cell
def _(get_template_at_width, np, skimage):
    def get_padded_template_at_width(template_width, *, base_template=None):
        if base_template is None:
            base_template = get_template_at_width(template_width)[:,:,0]
        pad = (2048-template_width)//2
        template = np.pad(base_template,(pad, pad)).astype("int8")
        dilated_template = skimage.morphology.isotropic_dilation(template, 50)
        template[(dilated_template > 0) & np.invert(template > 0)] = 0
        return template
    return (get_padded_template_at_width,)


@app.cell
def _(get_padded_template_at_width, plt):
    plt.imshow(get_padded_template_at_width(1326))
    return


@app.cell
def _(get_padded_template_at_width, plt):
    plt.plot(get_padded_template_at_width(1326)[:,1000])
    return


@app.cell
def _():
    100 & 1
    return


@app.cell
def _(get_padded_template_at_width, np):
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
    return (get_template_hat,)


@app.cell
def _(get_template_hat):
    template_hat = get_template_hat(1326)
    return


@app.cell
def _(img_path):
    offset_overrides = {
        "/groups/vale/valelab/_for_Mark/patterned_data/250731_patterned_plate_11_good/E05_250808_TRAK1_wt/Cell2.nd2": [256,128],
        "/groups/vale/valelab/_for_Mark/patterned_data/250731_patterned_plate_11_good/E05_250808_TRAK1_wt/Cell2 - Denoised.nd2": [256,128],
        "/groups/vale/valelab/_for_Mark/patterned_data/250731_patterned_plate_11_good/E05_250808_TRAK1_wt/Cell5.nd2": [204,128],
        "/groups/vale/valelab/_for_Mark/patterned_data/250731_patterned_plate_11_good/E05_250808_TRAK1_wt/Cell5 - Denoised.nd2": [204,128],

        "/groups/vale/valelab/_for_Mark/patterned_data/250731_patterned_plate_11_good/E05_250808_TRAK1_wt/Cell8.nd2": [256,128],
        "/groups/vale/valelab/_for_Mark/patterned_data/250731_patterned_plate_11_good/E05_250808_TRAK1_wt/Cell8 - Denoised.nd2": [256,128],
        "/groups/vale/valelab/_for_Mark/patterned_data/250612_patterned_plate_3/B06_250617_TRAK1_mDRH_dSp/Cell12.nd2": [64,128]
    }
    offset = offset_overrides.get(img_path, [128, 128])
    offset
    return (offset,)


@app.cell
def _():
    roi_overrides = {
        "/groups/vale/valelab/_for_Mark/patterned_data/250731_patterned_plate_11_good/F06_250811_TRAK1_mDRH_dSp/Cell3.nd2": [slice(None), slice(0,1200)]
    }
    return (roi_overrides,)


@app.cell
def _(nd2):
    #img_path = "/groups/vale/valelab/_for_Mark/patterned_data/250521_patterned_plate_1/B06_250528_TRAK1-wt/Cell3.nd2"
    #img_path = "/groups/vale/valelab/_for_Mark/patterned_data/250521_patterned_plate_1/B06_250528_TRAK1-wt/Cell8.nd2"
    #img_path = "/groups/vale/valelab/_for_Mark/patterned_data/250710_patterned_plate_9_good/C02_250718_NoV/cell8.nd2"
    #img_path = "/groups/vale/valelab/_for_Mark/patterned_data/250710_patterned_plate_9_good/C03_250718_TRAK1/Cell5.nd2"
    #img_path = "/groups/vale/valelab/_for_Mark/patterned_data/250710_patterned_plate_9_good/C03_250718_TRAK1/Cell2.nd2"
    # Low signal, need custom crop
    #img_path = "/groups/vale/valelab/_for_Mark/patterned_data/250731_patterned_plate_11_good/E05_250808_TRAK1_wt/Cell2.nd2"
    #img_path = "/groups/vale/valelab/_for_Mark/patterned_data/250731_patterned_plate_11_good/E05_250808_TRAK1_wt/Cell2 - Denoised.nd2"
    #img_path = "/groups/vale/valelab/_for_Mark/patterned_data/250731_patterned_plate_11_good/E05_250808_TRAK1_wt/Cell5.nd2"
    #img_path = "/groups/vale/valelab/_for_Mark/patterned_data/250731_patterned_plate_11_good/E05_250808_TRAK1_wt/Cell8.nd2"
    #img_path = "/groups/vale/valelab/_for_Mark/patterned_data/250612_patterned_plate_3/B06_250617_TRAK1_mDRH_dSp/Cell12.nd2"
    #img_path = "/groups/vale/valelab/_for_Mark/patterned_data/250731_patterned_plate_11_good/F06_250811_TRAK1_mDRH_dSp/Cell3.nd2"
    #img_path = "/groups/vale/valelab/_for_Mark/patterned_data/250710_patterned_plate_9_good/G09_250718_MAPK9_siRNA_Ars/denoised/Cell9 - Denoised.nd2"
    img_path = "/groups/vale/valelab/_for_Mark/patterned_data/250626_patterned_plate_7/F04_250703_TRAK2_wt_peroxisome/denoised/Cell6 - Denoised.nd2"
    img = nd2.imread(img_path, xarray=True)

    # check this one
    # /groups/vale/valelab/_for_Mark/patterned_data/250521_patterned_plate_1/B06_TRAK1_wt_combined/Cell12_1.nd2
    return img, img_path


@app.cell
def _(img):
    img
    return


@app.cell
def _(np, offset):
    def get_image_hat(img, offset=offset):
        img_template_sum_projection = np.sum(img.sel(C="640"), axis=0)
        img_template_sum_projection_norm = img_template_sum_projection / np.max(img_template_sum_projection)
        #img_template_sum_projection_hat = np.abs(np.fft.fft2(img_template_sum_projection))
        img_template_sum_projection_norm_2048 = img_template_sum_projection_norm[offset[0]:2048+offset[0],offset[1]:2048+offset[1]]
        img_template_sum_projection_norm_2048_hat = np.fft.fft2(img_template_sum_projection_norm_2048)
        return img_template_sum_projection_norm_2048_hat

    #img_template_hat = get_image_hat(img)
    return (get_image_hat,)


@app.cell
def _(get_image_hat, get_template_hat, nd2, np):
    def match_template(img, *, template_hat = None):
        if isinstance(img, str):
            img_path = img
            img = nd2.imread(img, xarray=True)
        if template_hat is None:
            template_hat = get_template_hat(1326)
        img_template_hat = get_image_hat(img)
        template_matching = np.fft.fftshift(np.real(np.fft.ifft2(template_hat * img_template_hat)))
        return template_matching
    return (match_template,)


@app.cell
def _(match_template, np):
    def max_match_template(img, *, template_hat = None, roi: tuple[slice,slice] = None):
        template_matching = match_template(img, template_hat = template_hat)
        if roi is not None:
            template_matching = template_matching[*roi]
        max_idx = np.argmax(template_matching)
        out = np.unravel_index(max_idx, template_matching.shape)
        if roi is not None:
            # reshift based on 
            if roi[0].start is not None:
                out = (out[0] + roi[0].start, out[1])
            if roi[1].start is not None:
                out = (out[0], out[1] + roi[1].start)
        return out
    return (max_match_template,)


@app.cell
def _(img_path, match_template):
    #template_matching = np.fft.fftshift(np.real(np.fft.ifft2(template_hat * img_template_hat)))
    template_matching = match_template(img_path)
    return (template_matching,)


@app.cell
def _(img, offset, plt, stretch01):
    plt.imshow(stretch01(img[:,:,offset[0]:2048+offset[0],offset[1]:2048+offset[1]].sum(axis=0).sel(C="640")))
    return


@app.cell
def _(img):
    img.shape
    return


@app.cell
def _(plt, template_matching):
    plt.imshow(template_matching)
    return


@app.cell
def _(img_path, max_match_template, roi_overrides):
    max_coords = max_match_template(img_path, roi = roi_overrides.get(img_path,(slice(None),slice(None))))
    max_coords
    return (max_coords,)


@app.cell
def _(get_padded_template_at_width, max_coords, np, plt):
    shifted_template = get_padded_template_at_width(1326)
    shifted_template = np.roll(shifted_template, (max_coords[0] - 1024,max_coords[1] - 1024), axis=(0,1))
    plt.imshow(shifted_template)
    plt.scatter(max_coords[1], max_coords[0])
    return (shifted_template,)


@app.cell
def _(get_padded_template_at_width, np):
    template_coords = np.where(get_padded_template_at_width(1326))
    template_coords_bounds = (
        min(template_coords[0]),
        max(template_coords[0]),
        min(template_coords[1]),
        max(template_coords[1])
    )
    return (template_coords_bounds,)


@app.cell
def _(template_coords_bounds):
    template_coords_bounds[1] - template_coords_bounds[0], template_coords_bounds[3] - template_coords_bounds[2]
    return


@app.cell
def _(sumproj):
    sumproj
    return


@app.cell
def _(max_coords, plt, sumproj, sumproj_threshold):
    plt.imshow(sumproj[max_coords[0]-512:max_coords[0]+512, max_coords[1]-512:max_coords[1]+512] > sumproj_threshold)
    return


@app.cell
def _(img, max_coords, np):
    cropped_proj_img = np.sum(img[:,:,max_coords[0]-512+128:max_coords[0]+512+128, max_coords[1]-512+128:max_coords[1]+512+128], axis=0)
    #cropped_proj_img.to_netcdf("test.nc")
    return (cropped_proj_img,)


@app.cell
def _():
    return


@app.cell
def _(xarray):
    reloaded = xarray.open_dataarray("test.nc")
    reloaded
    return (reloaded,)


@app.cell
def _(plt, reloaded):
    plt.imshow(reloaded.sel(C="561"), vmax=10000)
    return


@app.cell
def _(np):
    def stretch01(img, *, min_percentile=0.1, max_percentile=99.9):
        _min = np.percentile(img, min_percentile)
        _max = np.percentile(img, max_percentile)
        return np.clip((img - _min)/(_max - _min), 0, 1)
    return (stretch01,)


@app.cell
def _(cropped_proj_img):
    cropped_proj_img.sel(C="561").min()
    return


@app.cell
def _(np):
    def make_rgb(R, G, B):
        RGB = np.zeros([3, *R.shape[-2:]], dtype="float32")
        if R is not None:
            RGB[0,:,:] = R
        if G is not None:
            RGB[1,:,:] = G
        if B is not None:
            RGB[2,:,:] = B
        return np.permute_dims(RGB,(1,2,0))
    return (make_rgb,)


@app.cell
def _(max_coords, shifted_template_contour):
    cropped_template_contour = shifted_template_contour[0].copy()
    cropped_template_contour[:,0] -= (max_coords[0]-512)
    cropped_template_contour[:,1] -= (max_coords[1]-512)
    return (cropped_template_contour,)


@app.cell
def _(cropped_proj_img, cropped_template_contour, make_rgb, plt, stretch01):
    #cropped_RGB = np.zeros([3, *cropped_proj_img.shape[1:]], dtype="float32")
    #cropped_RGB[0,:,:] = stretch01(cropped_proj_img.sel(C="488"), max_percentile=100)
    #cropped_RGB[1,:,:] = 0 # stretch01(cropped_proj_img.sel(C="561"), max_percentile=99)
    #cropped_RGB[2,:,:] = stretch01(cropped_proj_img.sel(C="405"), min_percentile=10, max_percentile=99.9)
    #plt.imshow(np.permute_dims(cropped_RGB,(1,2,0)))
    plt.imshow(make_rgb(
        stretch01(cropped_proj_img.sel(C="488")),
        stretch01(cropped_proj_img.sel(C="640")),
        stretch01(cropped_proj_img.sel(C="405"))
    ))
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
    return


@app.cell
def _(cropped_proj_img, cropped_template_contour, np, plt, skimage):
    cropped_nuc_proj_threshold = skimage.filters.threshold_otsu(cropped_proj_img.sel(C="405").to_numpy())
    cropped_nuc_mask = cropped_proj_img.sel(C="405") > cropped_nuc_proj_threshold
    cropped_nuc_label = skimage.measure.label(cropped_nuc_mask)
    cropped_nuc_props = skimage.measure.regionprops(cropped_nuc_label)
    cropped_nuc_max_area = np.argmax([p.area for p in cropped_nuc_props])
    cropped_nuc_mask = (cropped_nuc_label == cropped_nuc_max_area+1)
    plt.imshow(cropped_nuc_mask)
    plt.plot(cropped_template_contour[1083:1951,1],cropped_template_contour[1083:1951,0], color="white")
    return (cropped_nuc_mask,)


@app.cell
def _(cropped_nuc_mask, distance_transform_edt, np, plt):
    cropped_nuc_edt = distance_transform_edt(np.invert(cropped_nuc_mask))
    plt.imshow(cropped_nuc_edt)
    return (cropped_nuc_edt,)


@app.cell
def _(
    cropped_nuc_mask,
    cropped_template_contour,
    distance_transform_edt,
    np,
    plt,
):
    top_arch_mask = np.zeros_like(cropped_nuc_mask)
    top_arch_mask[
        np.round(cropped_template_contour[1083:1951,0]).astype("int"),
        np.round(cropped_template_contour[1083:1951,1]).astype("int")
    ] = True
    cropped_arch_edt = distance_transform_edt(np.invert(top_arch_mask))
    plt.imshow(cropped_arch_edt)
    return (cropped_arch_edt,)


@app.cell
def _(cropped_nuc_edt, np, skimage):
    def get_nuclear_contour(nuclear_mask):
        nuclear_contours = skimage.measure.find_contours(nuclear_mask)
        nuclear_contour_index = np.argmax([len(contour) for contour in nuclear_contours])
        return nuclear_contours[nuclear_contour_index]

    cropped_nuclear_contour = get_nuclear_contour(cropped_nuc_edt < 7.69)
    return (cropped_nuclear_contour,)


@app.cell
def _(
    cropped_arch_edt,
    cropped_nuc_edt,
    cropped_template_contour,
    make_rgb,
    plt,
    stretch01,
):
    plt.imshow(make_rgb(
        stretch01(-cropped_arch_edt),
        stretch01(-cropped_nuc_edt),
        stretch01(-cropped_arch_edt)
    ))
    plt.plot(cropped_template_contour[1083:1951,1],cropped_template_contour[1083:1951,0], color="black")
    return


@app.cell
def _(
    cropped_arch_edt,
    cropped_nuc_edt,
    cropped_template_contour,
    perinuclear_space_distance_pixels,
    plt,
):
    plt.imshow(((cropped_arch_edt <= cropped_nuc_edt) & (cropped_arch_edt < 70)) < perinuclear_space_distance_pixels)
    plt.plot(cropped_template_contour[:,1], cropped_template_contour[:,0], color="white")
    return


@app.cell
def _(
    cropped_arch_edt,
    cropped_nuc_edt,
    cropped_nuclear_contour,
    cropped_proj_img,
    cropped_template_contour,
    make_rgb,
    np,
    perinuclear_space_distance_pixels,
    plt,
    stretch01,
):
    cropped_proj_mitochondria = cropped_proj_img.sel(C="488")
    cropped_proj_mitochondria_streched = stretch01(cropped_proj_img.sel(C="488"))
    #perinuclear_mitochondria = (cropped_arch_edt > cropped_nuc_edt) * cropped_proj_mitochondria_streched
    perinuclear_mask = cropped_nuc_edt < perinuclear_space_distance_pixels
    peripheral_mask = cropped_arch_edt <= perinuclear_space_distance_pixels
    perinuclear_mitochondria = perinuclear_mask * cropped_proj_mitochondria_streched
    peripheral_mitochondria = ((cropped_arch_edt <= cropped_nuc_edt) & peripheral_mask & np.invert(perinuclear_mask)) * cropped_proj_mitochondria_streched
    plt.imshow(make_rgb(
        peripheral_mitochondria,
        peripheral_mitochondria,
        perinuclear_mitochondria
    ))
    plt.plot(cropped_template_contour[1083:1951,1],cropped_template_contour[1083:1951,0], color="white", alpha=0.5)
    plt.plot(cropped_nuclear_contour[:,1], cropped_nuclear_contour[:,0], color="white", alpha=0.5)
    plt.plot([800, 800+perinuclear_space_distance_pixels],[900, 900], color="white")
    plt.text(790, 950, "5 μm", color="white")
    return (
        cropped_proj_mitochondria,
        cropped_proj_mitochondria_streched,
        perinuclear_mask,
        perinuclear_mitochondria,
        peripheral_mask,
        peripheral_mitochondria,
    )


@app.cell
def _(cropped_template_contour, peripheral_mask, plt):
    plt.imshow(peripheral_mask)
    plt.plot(cropped_template_contour[:,1], cropped_template_contour[:,0], color="white")
    return


@app.cell
def _(cropped_arch_edt, cropped_nuc_edt, cropped_template_contour, plt):
    plt.imshow((cropped_nuc_edt < 76.923) & (cropped_arch_edt > cropped_nuc_edt))
    plt.plot(cropped_template_contour[:,1], cropped_template_contour[:,0], color="white")
    return


@app.cell
def _(img):
    metadata = img.metadata["metadata"]
    return (metadata,)


@app.cell
def _(metadata):
    lateral_pixel_pitch = metadata.channels[0].volume.axesCalibration[0]
    perinuclear_space_distance_um = 5 # micrometers
    perinuclear_space_distance_pixels = perinuclear_space_distance_um/lateral_pixel_pitch
    return (perinuclear_space_distance_pixels,)


@app.cell
def _(metadata):
    metadata
    return


@app.cell
def _(cropped_proj_mitochondria_streched, plt, skimage):
    cropped_proj_mitochondria_streched_threshold = skimage.filters.threshold_otsu(cropped_proj_mitochondria_streched.to_numpy())
    cropped_proj_mitochondria_streched_thresholded = cropped_proj_mitochondria_streched > cropped_proj_mitochondria_streched_threshold
    plt.imshow(cropped_proj_mitochondria_streched_thresholded)
    return (cropped_proj_mitochondria_streched_thresholded,)


@app.cell
def _(cropped_proj_mitochondria, plt):
    plt.imshow(cropped_proj_mitochondria)
    return


@app.cell
def _(
    cropped_proj_mitochondria_streched_thresholded,
    np,
    perinuclear_mitochondria,
    peripheral_mitochondria,
):
    np.sum(cropped_proj_mitochondria_streched_thresholded * perinuclear_mitochondria), np.sum(cropped_proj_mitochondria_streched_thresholded * peripheral_mitochondria)
    return


@app.cell
def _(cropped_proj_mitochondria_streched_thresholded, plt):
    plt.imshow(cropped_proj_mitochondria_streched_thresholded)
    return


@app.cell
def _(cropped_proj_mitochondria_streched, plt):
    plt.imshow(cropped_proj_mitochondria_streched)
    return


@app.cell
def _(cropped_template_contour):
    cropped_template_contour[1083,:]
    return


@app.cell
def _(cropped_proj_mitochondria_streched, np):
    np.histogram(cropped_proj_mitochondria_streched[:,:275])
    return


@app.cell
def _(cropped_proj_mitochondria_streched, np):
    np.histogram(cropped_proj_mitochondria_streched)
    return


@app.cell
def _(cropped_proj_mitochondria_streched, plt):
    plt.imshow(cropped_proj_mitochondria_streched > 0.08)
    return


@app.cell
def _(cropped_proj_mitochondria_streched, np, plt):
    cropped_background = np.concatenate((cropped_proj_mitochondria_streched[:,:128], cropped_proj_mitochondria_streched[:,-128:]), axis=1)
    left_percentile = np.percentile(cropped_proj_mitochondria_streched[:,:128], 99.99)
    right_percentile = np.percentile(cropped_proj_mitochondria_streched[:,-128:], 99.99)
    if abs(left_percentile - right_percentile) > 0.1:
        cropped_background_threshold = min(left_percentile, right_percentile)
    else:
        cropped_background_threshold = np.percentile(cropped_background, 99.99)
    cropped_proj_mitochondria_streched_background_subtracted = cropped_proj_mitochondria_streched - cropped_background_threshold
    cropped_proj_mitochondria_streched_background_subtracted = np.clip(cropped_proj_mitochondria_streched_background_subtracted, 0, None)
    plt.imshow(cropped_proj_mitochondria_streched_background_subtracted)
    return (
        cropped_background,
        cropped_background_threshold,
        cropped_proj_mitochondria_streched_background_subtracted,
    )


@app.cell
def _(cropped_proj_mitochondria_streched, np):
    np.percentile(cropped_proj_mitochondria_streched[:,:128], 99.99)
    return


@app.cell
def _(cropped_proj_mitochondria_streched, np):
    np.percentile(cropped_proj_mitochondria_streched[:,-128:], 99.99)
    return


@app.cell
def _(cropped_background, cropped_proj_mitochondria_streched, np):
    np.percentile(cropped_background[cropped_background < np.max(cropped_proj_mitochondria_streched[:,-128:]).item()], 99.99)
    return


@app.cell
def _(cropped_proj_mitochondria_streched, np):
    np.max(cropped_proj_mitochondria_streched[:,-128:])
    return


@app.cell
def _(cropped_background, np):
    np.std(cropped_background)
    return


@app.cell
def _(cropped_proj_mitochondria_streched, np):
    np.percentile(cropped_proj_mitochondria_streched[:,:128],50)
    return


@app.cell
def _(cropped_background, np):
    outlier_limit = np.std(cropped_background)*3 + np.mean(cropped_background)
    return (outlier_limit,)


@app.cell
def _(cropped_background, np):
    np.histogram(cropped_background)
    return


@app.cell
def _(cropped_background, np):
    np.std(cropped_background)
    return


@app.cell
def _(outlier_limit):
    outlier_limit
    return


@app.cell
def _(cropped_background, np, outlier_limit):
    np.percentile(cropped_background[cropped_background < outlier_limit.item()], 99.99)
    return


@app.cell
def _(outlier_limit):
    outlier_limit.item()
    return


@app.cell
def _(cropped_background_threshold):
    cropped_background_threshold
    return


@app.cell
def _(
    cropped_proj_mitochondria_streched_background_subtracted,
    np,
    peripheral_mask,
):
    bg_sub_peripheral_per = np.sum(cropped_proj_mitochondria_streched_background_subtracted * peripheral_mask) / np.sum(cropped_proj_mitochondria_streched_background_subtracted)
    bg_sub_peripheral_per
    return (bg_sub_peripheral_per,)


@app.cell
def _(
    cropped_proj_mitochondria_streched_background_subtracted,
    np,
    perinuclear_mask,
):
    bg_sub_perinuclear_per = np.sum(cropped_proj_mitochondria_streched_background_subtracted * perinuclear_mask) / np.sum(cropped_proj_mitochondria_streched_background_subtracted)
    bg_sub_perinuclear_per
    return (bg_sub_perinuclear_per,)


@app.cell
def _(bg_sub_perinuclear_per, bg_sub_peripheral_per):
    bg_sub_peripheral_per / bg_sub_perinuclear_per
    return


@app.cell
def _(
    cropped_proj_mitochondria_streched_background_subtracted,
    np,
    perinuclear_mask,
    peripheral_mask,
):
    np.sum(cropped_proj_mitochondria_streched_background_subtracted * np.invert(peripheral_mask | perinuclear_mask)) / np.sum(cropped_proj_mitochondria_streched_background_subtracted)
    return


@app.cell
def _():
    return


@app.cell
def _(cropped_proj_mitochondria_streched_background_subtracted, np):
    np.sum(cropped_proj_mitochondria_streched_background_subtracted[:,:128]), np.sum(cropped_proj_mitochondria_streched_background_subtracted[:,-128:]), np.sum(cropped_proj_mitochondria_streched_background_subtracted[:,128:-128])
    return


@app.cell
def _(cropped_proj_mitochondria_streched, np):
    np.sum(np.clip(cropped_proj_mitochondria_streched[:,:275] - 0.1, 0, None))
    return


@app.cell
def _(cropped_proj_mitochondria_streched):
    cropped_proj_mitochondria_streched[:,:275]
    return


@app.cell
def _(cropped_template_contour):
    cropped_template_contour[12,:]
    return


@app.cell
def _(xarr):
    xarr
    return


@app.cell
def _(img):
    img
    return


@app.cell
def _(img, np, plt, skimage):
    sumproj = np.sum(img[:,1,128:2048+128,128:2048+128], axis=0)
    sumproj_threshold = skimage.filters.threshold_otsu(sumproj.to_numpy())
    sumproj_thresholded = sumproj > sumproj_threshold
    plt.imshow(sumproj_thresholded)
    return sumproj, sumproj_threshold, sumproj_thresholded


@app.cell
def _():
    128+ 2048
    return


@app.cell
def _(plt, shifted_template, sumproj_thresholded):
    plt.imshow(shifted_template)
    plt.imshow(sumproj_thresholded, alpha=0.5)
    return


@app.cell
def _(shifted_template, skimage):
    shifted_template_contour = skimage.measure.find_contours(shifted_template)
    shifted_template_contour
    return (shifted_template_contour,)


@app.cell
def _(plt, shifted_template_contour, sumproj_thresholded):
    plt.imshow(sumproj_thresholded)
    plt.plot(shifted_template_contour[0][:,1], shifted_template_contour[0][:,0], color="black")
    return


@app.cell
def _(plt, shifted_template_contour, stretch01, sumproj):
    plt.imshow(stretch01(sumproj))
    plt.plot(shifted_template_contour[0][:,1], shifted_template_contour[0][:,0], color="black")
    return


@app.cell
def _(shifted_template_contour):
    shifted_template_contour[0]
    return


@app.cell
def _(shifted_template_contour):
    shifted_template_contour[0]
    return


@app.cell
def _(np, plt, shifted_template, sumproj):
    fig = plt.gcf()

    axes_coords = [0.1, 0.1, 0.8, 0.8]

    ax_polar = fig.add_axes([axes_coords[0] + 0.05, axes_coords[1] - 0.175, *axes_coords[2:]], projection = 'polar')
    ax_polar.patch.set_alpha(0)
    #ax_polar.plot(theta, r)
    ax_polar.set_ylim(30, 41)
    ax_polar.set_yticks(np.arange(30, 41, 2))
    ax_polar.set_yticklabels([])
    ax_polar.set_rlabel_position(-22.5)  # get radial labels away from plotted line
    ax_polar.grid(True)
    ax_polar.set_title("Polar", va = 'bottom')

    ax_image = fig.add_axes(axes_coords)
    ax_image.imshow(sumproj, alpha = .5)
    ax_image.imshow(shifted_template, alpha=0.5)
    ax_image.axis('off')  # don't show the axes ticks/lines/etc. associated with the image


    #fig.add_axes(axes_coords)

    #plt.imshow(shifted_template)
    #plt.imshow(sumproj, alpha=0.5)
    fig
    return


@app.cell
def _(img, np):
    np.sum(img[:,1,:2048,:2048], axis=0)
    return


@app.cell
def _(np, shifted_template, sumproj_thresholded):
    score = (np.sum(sumproj_thresholded & shifted_template)/(np.sum(shifted_template > 0))).values.item()
    score
    return (score,)


@app.cell
def _(np, shifted_template, sumproj_thresholded):
    confusion_matrix = np.reshape([
        np.sum(sumproj_thresholded  & shifted_template).values.item(), np.sum( sumproj_thresholded & ~shifted_template).values.item(),
        np.sum(~sumproj_thresholded & shifted_template).values.item(), np.sum(~sumproj_thresholded & ~shifted_template).values.item(),
    ], (2,2))
    confusion_matrix
    return (confusion_matrix,)


@app.cell
def _(confusion_matrix, np):
    positive_predictive_value = confusion_matrix[0,0]/np.sum(confusion_matrix[:,0])
    positive_predictive_value
    return


@app.cell
def _(confusion_matrix, np):
    accuracy = (confusion_matrix[0,0] + confusion_matrix[1,1])/np.sum(confusion_matrix)
    accuracy
    return


@app.cell
def _(confusion_matrix, np):
    false_omission_rate = confusion_matrix[0,1] / np.sum(confusion_matrix[:,1])
    false_omission_rate
    return


@app.cell
def _(confusion_matrix, np):
    negative_predictive_value = confusion_matrix[1,1] / np.sum(confusion_matrix[:,1])
    negative_predictive_value
    return


@app.cell
def _(plt, shifted_template, sumproj_thresholded):
    plt.imshow(sumproj_thresholded & shifted_template)
    return


@app.cell
def _():
    2048*2048
    return


@app.cell
def _(img, np, plt):
    mip = np.sum(img.sel(C="405"), axis=0)
    print(np.max(mip))
    plt.imshow(mip > 5200)
    return (mip,)


@app.cell
def _(img, plt, skimage):
    nuc_proj = img.sel(C="405").sum(axis=0)[128:2048+128, 128:2048+128]
    nuc_proj_threshold = skimage.filters.threshold_otsu(nuc_proj.to_numpy())
    plt.imshow(nuc_proj > nuc_proj_threshold)
    return nuc_proj, nuc_proj_threshold


@app.cell
def _(img, plt, skimage):
    mitochondria_proj = img.sel(C="488").sum(axis=0)[128:2048+128, 128:2048+128]
    mitochondria_proj_threshold = skimage.filters.threshold_otsu(mitochondria_proj.to_numpy())
    plt.imshow(mitochondria_proj > mitochondria_proj_threshold)
    return mitochondria_proj, mitochondria_proj_threshold


@app.cell
def _(img, plt, skimage):
    trak_proj = img.sel(C="561").sum(axis=0)[128:2048+128, 128:2048+128]
    trak_proj_threshold = skimage.filters.threshold_otsu(trak_proj.to_numpy())
    plt.imshow(trak_proj > trak_proj_threshold)
    return trak_proj, trak_proj_threshold


@app.cell
def _(np, trak_proj, trak_proj_threshold):
    np.median(trak_proj), trak_proj_threshold
    return


@app.cell
def _(
    mitochondria_proj,
    mitochondria_proj_threshold,
    np,
    nuc_proj,
    nuc_proj_threshold,
    plt,
    shifted_template_contour,
    trak_proj,
    trak_proj_threshold,
):
    RGB = np.zeros([3, *trak_proj.shape], dtype="float32")
    RGB[0] = mitochondria_proj > mitochondria_proj_threshold
    RGB[1] = 0 if np.median(trak_proj) > trak_proj_threshold else trak_proj > trak_proj_threshold
    RGB[2] = nuc_proj > nuc_proj_threshold
    plt.imshow(np.permute_dims(RGB, (1,2,0)))
    plt.plot(shifted_template_contour[0][:,1], shifted_template_contour[0][:,0], color="white")
    return (RGB,)


@app.cell
def _(shifted_template_contour):
    shifted_template_contour[0][:,1]
    return


@app.cell
def _(np):
    def get_weighted_center(proj, mask):
        weights = np.ravel(proj)[np.ravel(mask)]
        sum = np.sum(weights)
        coords = np.where(mask)
        weighted_center = np.sum(coords[0] * weights) / sum, np.sum(coords[1] * weights) / sum
        return weighted_center
    return (get_weighted_center,)


@app.cell
def _(np):
    def get_weighted_center_and_variance(proj, mask):
        weights = np.ravel(proj)[np.ravel(mask)]
        sum = np.sum(weights)
        coords = np.where(mask)
        weighted_center = np.sum(coords[0] * weights) / sum, np.sum(coords[1] * weights) / sum
        weighted_variance = np.sum((coords[0] - weighted_center[0])**2 * weights)  / sum, np.sum((coords[1] - weighted_center[1])**2 * weights) / sum
        return weighted_center, weighted_variance
    return (get_weighted_center_and_variance,)


@app.cell
def _(mitochondria_mask, mitochondria_proj, np):
    weights = np.ravel(mitochondria_proj)[np.ravel(mitochondria_mask)]
    weights
    return


@app.cell
def _(mitochondria_mask, np):
    coords = np.where(mitochondria_mask)
    return (coords,)


@app.cell
def _(mitochondria_proj):
    mitochondria_proj
    return


@app.cell
def _(mitochondria_proj):
    mitochondria_proj.as_numpy
    return


@app.cell
def _(mitochondria_mask, mitochondria_proj, pl, skimage):
    mitochondria_labels = skimage.morphology.label(mitochondria_mask)
    mitochondria_properties = skimage.measure.regionprops(mitochondria_labels, intensity_image=mitochondria_proj.to_numpy())
    mitochondira_proptable = skimage.measure.regionprops_table(mitochondria_labels, intensity_image=mitochondria_proj.to_numpy(), properties=["area", "mean_intensity"])
    pl.dataframe.DataFrame(mitochondira_proptable)
    return mitochondria_labels, mitochondria_properties


@app.cell
def _(mitochondria_labels, mitochondria_properties, np, plt):
    mitochondria_labels_filtered = np.zeros_like(mitochondria_labels)
    for prop in mitochondria_properties:
        if prop.area > 1:
            mitochondria_labels_filtered[mitochondria_labels == prop.label] = prop.label
    plt.imshow(mitochondria_labels_filtered)
    plt.scatter(np.where(mitochondria_labels_filtered)[1], np.where(mitochondria_labels_filtered)[0])
    return


@app.cell
def _(template_matching):
    template_matching
    return


@app.cell
def _(shifted_template_contour):
    (
        min(shifted_template_contour[0][:,0]),
        max(shifted_template_contour[0][:,0]),
        min(shifted_template_contour[0][:,1]),
        max(shifted_template_contour[0][:,1])
    )
    return


@app.cell
def _(coords, mitochondria_mask, plt):
    plt.imshow(mitochondria_mask)
    plt.scatter(coords[1], coords[0])
    return


@app.cell
def _():
    return


@app.cell
def _(max_coords):
    bottom_coords = max_coords[0] + 350, max_coords[1]
    return (bottom_coords,)


@app.cell
def _(max_coords):
    max_coords
    return


@app.cell
def _(get_weighted_center_and_variance, np, nuc_mask):
    np.sqrt(get_weighted_center_and_variance(nuc_mask, nuc_mask)[1])
    return


@app.cell
def _(
    get_weighted_center_and_variance,
    mitochondria_mask,
    mitochondria_proj,
    np,
):
    np.sqrt(get_weighted_center_and_variance(mitochondria_proj, mitochondria_mask)[1])
    return


@app.cell
def _(bottom_coords, get_weighted_center, np, nuc_proj, nuc_proj_threshold):
    nuc_mask = nuc_proj > nuc_proj_threshold
    nuc_vec = get_weighted_center(nuc_proj, nuc_mask)
    nuc_vec2 = np.array(nuc_vec) - np.array(bottom_coords)
    return nuc_mask, nuc_vec, nuc_vec2


@app.cell
def _(bottom_coords, np, nuc_vec):
    np.array(nuc_vec) - np.array(bottom_coords)
    return


@app.cell
def _(bottom_coords, nuc_vec):
    nuc_vec[1] - bottom_coords[1], nuc_vec[0] - bottom_coords[0]
    return


@app.cell
def _(
    bottom_coords,
    get_weighted_center,
    mitochondria_proj,
    mitochondria_proj_threshold,
    np,
):
    mitochondria_mask = mitochondria_proj > mitochondria_proj_threshold
    mitochondira_vec = get_weighted_center(mitochondria_proj, mitochondria_mask)
    mitochondira_vec2 = np.array(mitochondira_vec) - np.array(bottom_coords)
    return mitochondira_vec, mitochondira_vec2, mitochondria_mask


@app.cell
def _(mitochondria_proj, plt):
    plt.imshow(mitochondria_proj, vmax=9000)
    return


@app.cell
def _(nuc_proj, plt):
    plt.imshow(nuc_proj, vmax=7500)
    return


@app.cell
def _(mitochondria_proj):
    mitochondria_proj
    return


@app.cell
def _(mitochondria_mask, mitochondria_proj, np, plt):
    mitochondria_proj_histogram = np.histogram(np.ravel(mitochondria_proj)[np.ravel(mitochondria_mask)])
    plt.stairs(mitochondria_proj_histogram[0], mitochondria_proj_histogram[1])
    return


@app.cell
def _(
    RGB,
    bottom_coords,
    mitochondira_vec,
    np,
    nuc_vec,
    plt,
    shifted_template_contour,
):
    plt.imshow(np.permute_dims(RGB, (1,2,0)))
    plt.plot(shifted_template_contour[0][:,1], shifted_template_contour[0][:,0], color="white")
    plt.scatter(bottom_coords[1], bottom_coords[0])
    plt.scatter(nuc_vec[1], nuc_vec[0])
    plt.scatter(mitochondira_vec[1], mitochondira_vec[0])
    plt.arrow(
        bottom_coords[1], bottom_coords[0],
        nuc_vec[1] - bottom_coords[1], nuc_vec[0] - bottom_coords[0],
        color="red",
        width=10
    )
    plt.arrow(
        bottom_coords[1], bottom_coords[0],
        mitochondira_vec[1] - bottom_coords[1], mitochondira_vec[0] - bottom_coords[0],
        color="red",
        width=10,
        linestyle="--"
    )
    return


@app.cell
def _(np, shifted_template_contour):
    np.argmin(shifted_template_contour[0][:,1]), np.argmax(shifted_template_contour[0][:,1])
    return


@app.cell
def _(shifted_template_contour):
    shifted_template_contour[0]
    return


@app.cell
def _(RGB, np, plt, shifted_template_contour):
    plt.imshow(np.permute_dims(RGB, (1,2,0)))
    plt.plot(shifted_template_contour[0][1083:1951,1], shifted_template_contour[0][1083:1951,0], color="white")
    plt.plot(shifted_template_contour[0][1083:1951,1], shifted_template_contour[0][1083:1951,0]+200, color="white")
    return


@app.cell
def _(bottom_coords, plt, shifted_template, shifted_template_contour):
    plt.imshow(shifted_template)
    plt.plot(shifted_template_contour[0][1083:1951,1], shifted_template_contour[0][1083:1951,0], color="white")
    plt.scatter(shifted_template_contour[0][1083,1], shifted_template_contour[0][1083,0]+60, color="white")
    plt.scatter(shifted_template_contour[0][1951,1], shifted_template_contour[0][1951,0]+60, color="white")
    plt.scatter(bottom_coords[1], bottom_coords[0]+20, color="white")
    return


@app.cell
def _(shifted_template_contour):
    shifted_template_contour[0][1083,1]
    return


@app.cell
def _(shifted_template_contour):
    shifted_template_contour[0][1951,1]
    return


@app.cell
def _(plt, shifted_template):
    plt.imshow(shifted_template[:,:718])
    return


app._unparsable_cell(
    r"""
     shifted_template_contour[0][1083:1951,:]
    """,
    name="_"
)


@app.cell
def _():
    return


@app.cell
def _(shifted_template_contour):
    len(shifted_template_contour[0][:,1])
    return


@app.cell
def _(np, nuc_vec2):
    nuc_vec_length = np.hypot(nuc_vec2[1], nuc_vec2[0])
    nuc_vec_length
    return


@app.cell
def _(mitochondira_vec2, np):
    mitochondira_vec_length = np.hypot(mitochondira_vec2[1], mitochondira_vec2[0])
    mitochondira_vec_length
    return


@app.cell
def _(mitochondira_vec2, np, nuc_vec2):
    diff_length = np.hypot(mitochondira_vec2[1]-nuc_vec2[1], mitochondira_vec2[0]-nuc_vec2[0])
    diff_length
    return


@app.cell
def _(trak_proj_threshold):
    trak_proj_threshold
    return


@app.cell
def _():
    2290 - 2048
    return


@app.cell
def _(mip, np):
    mip_hist = np.histogram(mip[:], range(4000, 7000, 100))
    mip_hist
    return (mip_hist,)


@app.cell
def _(mip_hist, plt):
    plt.stairs(mip_hist[0], mip_hist[1])
    return


@app.cell
def _(img):
    img.C
    return


@app.cell
def _(img):
    img.sel(C="640").ndim
    return


@app.cell
def _(img, np, plt):
    plt.imshow(np.max(img[:,3,:,:], axis=0))
    return


@app.cell
def _(img):
    img.Z
    return


@app.cell
def _(img):
    img.C
    return


@app.cell
def _(img_path, pathlib):
    pdf_path = pathlib.Path("template_matching",*pathlib.Path(img_path).parts[-3:]).with_suffix(".pdf")
    pdf_path
    return (pdf_path,)


@app.cell
def _(pdf_path):
    pdf_path.parent
    return


@app.cell
def _(img_path):
    img_path
    return


@app.cell
def _(score):
    str(score)
    return


@app.cell
def _(img_path, pathlib):
    for (dirpath, dirnames, filenames) in pathlib.Path(img_path).parent.parent.parent.walk():
        for filename in filenames:
            if filename.endswith(".nd2") and filename.startswith("Cell"):
                print(dirpath / filename)
    return


@app.cell
def _():
    float('nan')
    return


@app.cell
def _(img):
    len(img.sel(C="640").dims) > 3
    return


@app.cell
def _(nd2):
    img2 = nd2.imread("/groups/vale/valelab/_for_Mark/patterned_data/250521_patterned_plate_1/B06_250529_TRAK1-wt/Cell10.nd2", xarray=True, dask=True)
    return


@app.cell
def _():
    return


@app.cell
def _(nd2, pathlib):
    def shit_data_scan():
        root_path = "/groups/vale/valelab/_for_Mark/patterned_data"
        for (dirpath, dirnames, filenames) in pathlib.Path(root_path).walk():
            relative_path = dirpath.relative_to(root_path)
            relative_depth = len(relative_path.parts)
            if relative_depth != 2:
                #print(f"Ignoring {dirpath=}")
                continue
            for filename in filenames:
                 if filename.endswith(".nd2") and (filename.startswith("Cell") or filename.startswith("cell")):
                    img = nd2.imread(
                        pathlib.Path(dirpath) / filename,
                        xarray=True,
                        dask=True
                    )
                    try:
                        if len(img.sel(C="640").dims) > 3:
                            print(pathlib.Path(dirpath) / filename)
                    except Exception as e:
                        print(f"Error: {e}")
                        print(pathlib.Path(dirpath) / filename)
        print("Scan complete!")

    shit_data_scan()
    return


@app.cell
def _(nd2):
    img3 = nd2.imread("/groups/vale/valelab/_for_Mark/patterned_data/250626_patterned_plate_7/C09_250630_ctrl_siRNA/Cell6_DIC.nd2", xarray=True, dask=True)
    img3
    return


@app.cell
def _(img, plt):
    img488 = img.sel(C="488").sum(axis=0)
    plt.imshow(img488)
    return (img488,)


@app.cell
def _(max_coords):
    max_coords
    return


@app.cell
def _(img488):
    img488
    return


@app.cell
def _():
    # https://hhmionline-my.sharepoint.com/:x:/r/personal/gladkovac_hhmi_org/_layouts/15/Doc.aspx?sourcedoc=%7BBA3E9993-FDDF-46A6-8F7F-FA6F7EA4E62D%7D&file=Comparisons_table.xlsx&action=default&mobileredirect=true&DefaultItemOpen=1&web=1
    # TRAK1 Helix Muts
    # Overlay sum projection over many cells
    # minimum goal is the overlay
    return


if __name__ == "__main__":
    app.run()
