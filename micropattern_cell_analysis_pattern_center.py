import marimo

__generated_with = "0.14.15"
app = marimo.App(width="medium", app_title="Micropattern Cell Analysis")


@app.cell
def _():
    import marimo as mo
    import os
    import pathlib
    import nd2
    import numpy as np
    import matplotlib.pyplot as plt
    from pathlib import Path
    import skimage
    import xarray as xr
    return Path, mo, nd2, np, plt, skimage


@app.cell
def _(mo):
    mo.md(r"""https://resisted-curiosity-682.notion.site/Micropatterned-cell-analysis-1fc79054849480e887f6d45ba3aeecfb""")
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
    1. File structure:
        - parent folder corresponds to a 96 well plate that was plated and fixed on the same day
            - \\[prfs.hhmi.org](http://prfs.hhmi.org/)\valelab\\Gaby\Vale\imaging\2025\250521_round_E_patterned_1
        - subfolders correspond to individual wells with different conditions (in this case expressing different variants of TRAK); in the name they contain information about the date imaged and the condition
        - each subfolder contains .nd2 stacks corresponding to a cell that was acquired

    2. Data for each cell:
        - 4 colour z-stacks through single patterned cells as .nd2 files
        405 - nuclear stain (Hoechst dye)
        488 - organelle of interest = mitochondria or peroxisomes
        561 - expressed TRAK protein that is expected to affect distribution
        640 - micro pattern visualised by Fibronectin-647
        - we might consider processing by denoising using NIS Elements; this could be very effective to boost our signal:noise ratio
    """
    )
    return


@app.cell
def _():
    data_path = "valelab/Gaby/Vale/imaging/2025/250521_patterned_plate_1"
    return (data_path,)


@app.cell
def _(Path, data_path):
    datasets = [str(d) for d in Path(data_path).iterdir() if d.is_dir()]
    return (datasets,)


@app.cell
def _(datasets, mo):
    dataset_dropdown = mo.ui.dropdown(options=datasets, label="Select Dataset", value=datasets[0])
    dataset_dropdown
    return (dataset_dropdown,)


@app.cell
def _(Path, dataset_dropdown):
    nd2_images = [str(image) for image in Path(dataset_dropdown.selected_key).iterdir() if image.suffix == ".nd2" and image.name.startswith("Cell")]
    return (nd2_images,)


@app.cell
def _(mo, nd2_images):
    images_dropdown = mo.ui.dropdown(options=nd2_images, label = "Select Image", value=nd2_images[0])
    images_dropdown
    return (images_dropdown,)


@app.cell
def _(Path, images_dropdown, nd2):
    image_path = Path(images_dropdown.selected_key)
    image = nd2.imread(image_path, xarray=True, dask=True)
    return (image,)


@app.cell
def _(Path, data_path, nd2):
    cell_1 = nd2.imread(Path(data_path)/"B06_250528_TRAK1-wt/Cell1.nd2", xarray=True)
    return


@app.cell
def _(np):
    def scale(arr):
        min = np.min(arr)
        max = np.max(arr)
        return (arr - min)/(max-min)
    return (scale,)


@app.cell
def _(image, mo):
    channel_dropdown = mo.ui.dropdown(
        options=[str(c) for c in image.C.values],
        value=image.C.values[0],
        label="Channel"
    )
    channel_dropdown
    return (channel_dropdown,)


@app.cell
def _(image):
    image
    return


@app.cell
def _(image, mo):
    z_slider = mo.ui.slider(steps=image.Z.values, full_width=True, label="Z")
    return (z_slider,)


@app.cell
def _(mo):
    image_scale_slider = mo.ui.range_slider(
        orientation="vertical",
        start=0.0,
        stop=1.0,
        step=0.01,
        full_width=True,
        show_value=True)
    return (image_scale_slider,)


@app.cell
def _(channel_dropdown, image, z_slider):
    image_CZ = image.sel(C=channel_dropdown.selected_key, Z=z_slider.value)
    return (image_CZ,)


@app.cell
def _(image_CZ, image_scale_slider, plt, scale):
    def imshow_cz():
        plt.imshow(
            scale(image_CZ),
            vmin=image_scale_slider.value[0],
            vmax=image_scale_slider.value[1]
        )
        #plt.scatter(centroid[0], centroid[1], color='red', marker='x')
        return plt.gca()

    return (imshow_cz,)


@app.cell
def _(
    channel_dropdown,
    dataset_dropdown,
    image_scale_slider,
    images_dropdown,
    imshow_cz,
    mo,
    z_slider,
):
    mo.vstack([
        dataset_dropdown,
        images_dropdown,
        mo.hstack([
            imshow_cz(),
            image_scale_slider
        ]),
        mo.hstack([channel_dropdown,z_slider])
    ])
    return


@app.cell
def _(mean_values, polar_bar):
    polar_bar(mean_values)
    return


@app.cell
def _(mean_values, np, plt):
    plt.plot(np.linspace(0, 360, 37)[0:36], mean_values)
    return


@app.cell
def _(pattern_mip_scaled, plt, rp):
    plt.imshow(pattern_mip_scaled, vmax=0.2)
    plt.scatter(rp[0].centroid[0], rp[0].centroid[1], color='red', marker='x')
    return


@app.cell
def _(pattern_mip_scaled, skimage):
    pattern_mip_scaled_threshold = skimage.filters.threshold_otsu(pattern_mip_scaled.values)
    return (pattern_mip_scaled_threshold,)


@app.cell
def _(pattern_mip_scaled, pattern_mip_scaled_threshold, plt):
    plt.imshow(pattern_mip_scaled > pattern_mip_scaled_threshold)
    return


@app.cell
def _(pattern_mip_scaled, pattern_mip_scaled_threshold, plt, skimage):
    pattern_mip_scaled_dilated = skimage.morphology.dilation(pattern_mip_scaled, footprint=skimage.morphology.disk(radius=20))
    pattern_mip_scaled_dilated_binary = pattern_mip_scaled_dilated > pattern_mip_scaled_threshold
    plt.imshow(pattern_mip_scaled_dilated_binary)
    return (pattern_mip_scaled_dilated_binary,)


@app.cell
def _(pattern_mip_scaled_dilated_binary, plt, skimage):
    pattern_mip_scaled_dilated_label = skimage.measure.label(pattern_mip_scaled_dilated_binary)
    plt.imshow(pattern_mip_scaled_dilated_label)
    return (pattern_mip_scaled_dilated_label,)


@app.cell
def _(pattern_mip_scaled_dilated_label, skimage):
    pattern_rp = skimage.measure.regionprops(pattern_mip_scaled_dilated_label)
    return (pattern_rp,)


@app.cell
def _(np, pattern_mip_scaled_dilated_label, pattern_rp):
    pattern_areas = [region.area for region in pattern_rp]
    max_pattern_area = max(pattern_areas)
    arg_max_pattern_area = np.argmax(np.array(pattern_areas))
    pattern_binary = pattern_mip_scaled_dilated_label == pattern_rp[arg_max_pattern_area].label
    return arg_max_pattern_area, pattern_binary


@app.cell
def _(arg_max_pattern_area, pattern_binary, pattern_rp, plt):
    pattern_centroid = pattern_rp[arg_max_pattern_area].centroid
    plt.imshow(pattern_binary)
    plt.scatter(pattern_centroid[1], pattern_centroid[0], color='magenta', marker='x')
    return


@app.cell
def _(pattern_binary):
    pattern_binary
    return


if __name__ == "__main__":
    app.run()
