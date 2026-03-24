import marimo

__generated_with = "0.14.15"
app = marimo.App(
    width="medium",
    app_title="Micropattern Cell Analysis",
    layout_file="layouts/micropattern_cell_analysis.slides.json",
)


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
    return Path, mo, nd2, np, plt, skimage, xr


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
def _(centroid, image_CZ, image_scale_slider, plt, scale):
    def imshow_cz():
        plt.imshow(
            scale(image_CZ),
            vmin=image_scale_slider.value[0],
            vmax=image_scale_slider.value[1]
        )
        plt.scatter(centroid[0], centroid[1], color='red', marker='x')
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
    return (pattern_centroid,)


@app.cell
def _(pattern_binary):
    pattern_binary
    return


@app.cell
def _(image, np, plt, scale, skimage):
    nucleus_mip = image.sel(C="405").max(axis=0)
    nucleus_mip_scaled = scale(nucleus_mip)
    plt.imshow(nucleus_mip_scaled, vmax=0.01)
    threshold_nucleus_mip_scaled = 0.01
    nucleus_mip_binary = nucleus_mip_scaled > threshold_nucleus_mip_scaled
    rp = skimage.measure.regionprops(nucleus_mip_binary.compute().values.astype(np.uint8))
    centroid = rp[0].centroid
    plt.scatter(centroid[0], centroid[1], color='red', marker='x')
    return centroid, rp


@app.cell
def _(image, plt, scale):
    pattern_mip = image.sel(C="640").max(axis=0)
    pattern_mip_scaled = scale(pattern_mip)
    plt.imshow(pattern_mip)
    return pattern_mip, pattern_mip_scaled


@app.cell
def _(pattern_mip, plt, radius_mask):
    plt.imshow(pattern_mip * radius_mask)
    return


@app.cell
def _(image_CZ, plt, radius_mask, scale):
    plt.imshow(scale(image_CZ * radius_mask), vmax=0.2)
    return


@app.cell
def _(channel_dropdown, image, plt, scale):
    image_C_sum = image.sel(C=channel_dropdown.selected_key).sum(axis=0).compute()
    plt.imshow(scale(image_C_sum), vmax=0.2)
    return (image_C_sum,)


@app.cell
def _(mo):
    radius_range_slider = mo.ui.range_slider(start=0, stop=26, show_value=True, full_width=True)
    radius_range_slider
    return (radius_range_slider,)


@app.cell
def _(pattern_mip):
    pattern_mip.compute()
    return


@app.cell
def _(image):
    image.X.values
    return


@app.cell
def _(image):
    image.coords
    return


@app.cell
def _(image, xr):
    XY = xr.broadcast(image.Y,image.X)
    return (XY,)


@app.cell
def _(XY):
    XY[0].values,XY[1].values
    return


@app.cell
def _(XY, image, np, pattern_centroid, plt):
    theta = np.atan2(
            XY[0].values-image.Y.values[round(pattern_centroid[0])],
            XY[1].values-image.X.values[round(pattern_centroid[1])]
        )
    plt.imshow(
        theta,
        cmap="hsv",
        vmin=-np.pi,
        vmax=np.pi
    )
    return (theta,)


@app.cell
def _(radius_range_slider):
    radius_range_slider.value
    return


@app.cell
def _(XY, image, np, pattern_centroid, radius_range_slider):
    radius = np.hypot(
        XY[0].values-image.Y.values[round(pattern_centroid[0])],
        XY[1].values-image.X.values[round(pattern_centroid[1])]
    )
    radii_bins = np.arange(radius_range_slider.value[0], radius_range_slider.value[1], 2)
    radius_mask = (radius >= radius_range_slider.value[0]) & (radius < radius_range_slider.value[1])
    return radii_bins, radius, radius_mask


@app.cell
def _(radii_bins):
    radii_bins
    return


@app.cell
def _(theta):
    theta.min()
    return


@app.cell
def _(XY):
    XY[1].values
    return


@app.cell
def _(image):
    image.X
    return


@app.cell
def _(centroid, image):
    image.X.values[round(centroid[0])]
    return


@app.cell
def _(theta):
    theta
    return


@app.cell
def _(np):
    theta_bins = np.linspace(-np.pi, np.pi, 37)
    return (theta_bins,)


@app.cell
def _(theta_bins):
    theta_bins
    return


@app.cell
def _(radius_mask, theta, theta_bins):
    theta_groups = [None]*(len(theta_bins)-1)
    for i in range(len(theta_bins)-1):
        theta_groups[i] = (theta >= theta_bins[i]) & (theta < theta_bins[i+1]) & radius_mask
    return (theta_groups,)


@app.cell
def _(np, plt, theta_groups):
    theta_group_map = np.zeros(theta_groups[0].shape, dtype="uint8")
    for j in range(len(theta_groups)):
        theta_group_map[theta_groups[j]] = (j+1)
    plt.imshow(theta_group_map)
    return


@app.cell
def _(radii_bins, radius):
    def get_radii_groups():
        radii_groups = [None]*(len(radii_bins)-1)
        for i in range(len(radii_bins)-1):
            radii_groups[i] = (radius >= radii_bins[i]) & (radius < radii_bins[i+1])
        return radii_groups

    radii_groups = get_radii_groups()
    return (radii_groups,)


@app.cell
def _(np, plt, radii_groups):
    def show_radii_group_map():
        radii_group_map = np.zeros(radii_groups[0].shape, dtype="uint8")
        for j in range(len(radii_groups)):
            radii_group_map[radii_groups[j]] = (j+1)
        return plt.imshow(radii_group_map)

    show_radii_group_map()
    return


@app.cell
def _(np, plt, radii_groups, theta_groups):
    def show_theta_radii_map():
        counter = 1
        theta_radii_group_map = np.zeros(radii_groups[0].shape, dtype="uint8")
        for i in range(len(theta_groups)):
            for j in range(len(radii_groups)):
                theta_radii_group_map[theta_groups[i] * radii_groups[j]] = np.random.randint(1,100)
                counter += 1
        return plt.imshow(theta_radii_group_map)

    show_theta_radii_map()
    return


@app.cell
def _(image_C_sum, np, plt, radii_groups, theta_groups):
    def get_theta_radii_means(image):
        theta_radii_means = np.zeros((len(theta_groups), len(radii_groups)))
        for i in range(len(theta_groups)):
            for j in range(len(radii_groups)):
                theta_radii_means[i,j] = image.values[theta_groups[i] * radii_groups[j]].mean()
        return theta_radii_means

    theta_radii_means = get_theta_radii_means(image_C_sum)
    plt.imshow(theta_radii_means)
    return (theta_radii_means,)


@app.cell
def _(theta_radii_means):
    theta_radii_means_sig_count = (theta_radii_means > 5500).sum(axis=0).reshape(1,6)
    theta_radii_means_sig_count
    return


@app.cell
def _(plt, theta_radii_means):
    plt.imshow(theta_radii_means.max(axis=0).reshape(1,6))
    return


@app.cell
def _(theta_radii_means):
    theta_radii_means
    return


@app.cell
def _(np, plt, radii_groups, scale, theta_groups, theta_radii_means):
    def get_theta_radii_mean_map():
        theta_radii_mean_map = np.ones(radii_groups[0].shape)* theta_radii_means.min()
        for i in range(len(theta_groups)):
            for j in range(len(radii_groups)):
                theta_radii_mean_map[theta_groups[i] * radii_groups[j]] = theta_radii_means[i,j]
        return theta_radii_mean_map

    theta_radii_mean_map = get_theta_radii_mean_map()
    plt.imshow(scale(theta_radii_mean_map))
    return (theta_radii_mean_map,)


@app.cell
def _(theta_radii_mean_map):
    theta_radii_mean_map
    return


@app.cell
def _(image_C_sum, plt, scale):
    plt.imshow(scale(image_C_sum), vmax=0.2)
    return


@app.cell
def _(image_CZ, plt, theta_groups):
    plt.imshow(image_CZ.values * theta_groups[3])
    return


@app.cell
def _(theta_groups):
    theta_groups[0].shape
    return


@app.cell
def _(image_CZ):
    image_CZ.shape
    return


@app.cell
def _(image_CZ, np, theta_groups):
    mean_values = np.zeros(len(theta_groups))
    for i2 in range(len(theta_groups)):
        mean_values[i2] = (image_CZ.values[theta_groups[i2]]).mean()
    return (mean_values,)


@app.cell
def _(mean_values):
    mean_values
    return


@app.cell
def _(np, plt):
    def polar_bar(values, min_to_zero=False):
        if min_to_zero:
            values = values - values.min()
        fig = plt.figure(figsize=[5,5])
        ax = fig.add_axes([0.1,0.1,0.8,0.8], polar=True)
        ax.bar(np.linspace(0, 2*np.pi, 37)[0:36], values, width=np.pi/36*2, bottom=0)
        ax.set_theta_offset(np.pi)
        ax.set_theta_direction("clockwise")
        plt.ylim(0, values.max())
        return fig
    return (polar_bar,)


@app.cell
def _():
    x = 5
    return (x,)


@app.cell
def _(x):
    x + 3
    return


@app.cell
def _(x):
    2*x
    return


@app.cell
def _(
    channel_dropdown,
    nd2,
    np,
    radius_range_slider,
    scale,
    skimage,
    theta_bins,
    xr,
):
    def analyze_cell(image_path):
        image = nd2.imread(image_path, xarray=True, dask=True)
        image_C_sum = image.sel(C=channel_dropdown.selected_key).sum(axis=0).compute()

        pattern_mip = image.sel(C="640").max(axis=0)
        pattern_mip_scaled = scale(pattern_mip)

       #return pattern_mip_scaled

        pattern_mip_scaled_dilated = skimage.morphology.dilation(pattern_mip_scaled, footprint=skimage.morphology.disk(radius=20))
        pattern_mip_scaled_threshold = skimage.filters.threshold_otsu(pattern_mip_scaled.values)
        pattern_mip_scaled_dilated_binary = pattern_mip_scaled_dilated > pattern_mip_scaled_threshold
        pattern_mip_scaled_dilated_label = skimage.measure.label(pattern_mip_scaled_dilated_binary)

        pattern_rp = skimage.measure.regionprops(pattern_mip_scaled_dilated_label)

        pattern_areas = [region.area for region in pattern_rp]
        max_pattern_area = max(pattern_areas)
        arg_max_pattern_area = np.argmax(np.array(pattern_areas))
        pattern_binary = pattern_mip_scaled_dilated_label == pattern_rp[arg_max_pattern_area].label

        pattern_centroid = pattern_rp[arg_max_pattern_area].centroid

        XY = xr.broadcast(image.Y,image.X)

        radius = np.hypot(
            XY[0].values-image.Y.values[round(pattern_centroid[0])],
            XY[1].values-image.X.values[round(pattern_centroid[1])]
        )
        radii_bins = np.arange(radius_range_slider.value[0], radius_range_slider.value[1], 2)
        radius_mask = (radius >= radius_range_slider.value[0]) & (radius < radius_range_slider.value[1])

        def get_radii_groups():
            radii_groups = [None]*(len(radii_bins)-1)
            for i in range(len(radii_bins)-1):
                radii_groups[i] = (radius >= radii_bins[i]) & (radius < radii_bins[i+1])
            return radii_groups

        radii_groups = get_radii_groups()

        theta = np.atan2(
            XY[0].values-image.Y.values[round(pattern_centroid[0])],
            XY[1].values-image.X.values[round(pattern_centroid[1])]
        )

        theta_groups = [None]*(len(theta_bins)-1)
        for i in range(len(theta_bins)-1):
            theta_groups[i] = (theta >= theta_bins[i]) & (theta < theta_bins[i+1]) & radius_mask

        def get_theta_radii_means(image):
            theta_radii_means = np.zeros((len(theta_groups), len(radii_groups)))
            for i in range(len(theta_groups)):
                for j in range(len(radii_groups)):
                    theta_radii_means[i,j] = image.values[theta_groups[i] * radii_groups[j]].mean()
            return theta_radii_means

        theta_radii_means = get_theta_radii_means(image_C_sum)

        return {
            "theta_radii_means": theta_radii_means,
            "pattern_mip_scaled_dilated_binary": pattern_mip_scaled_dilated_binary,
            "pattern_centroid": pattern_centroid,
            "image_C_sum": image_C_sum
        }

    return (analyze_cell,)


@app.cell
def _(nd2_images):
    nd2_images
    return


@app.cell
def _(nd2_images):
    nd2_images
    return


@app.cell
def _(plt, results_1, results_slider):
    plt.imshow(results_1[results_slider.value]["theta_radii_means"])
    return


@app.cell
def _(mo, results_1):
    results_slider = mo.ui.slider(start=0, stop=len(results_1)-1)
    results_slider
    return (results_slider,)


@app.cell
def _(plt, results_1, results_slider):
    plt.plot(results_1[results_slider.value]["theta_radii_means"])
    return


@app.cell
def _(plt, results_1, results_slider):
    plt.plot(results_1[results_slider.value]["theta_radii_means"].T)
    return


@app.cell
def _(datasets):
    datasets
    return


@app.cell
def _(Path, analyze_cell):
    def get_results(dataset):
        nd2_images = [
            str(image) for image in Path(dataset).iterdir()
            if image.suffix == ".nd2" and image.name.startswith("Cell")
        ]
        results = [None] * len(nd2_images)
        for n, nd2_image in enumerate(nd2_images):
            print(n, nd2_images)
            try:
                results[n] = analyze_cell(nd2_image)
            except:
                print("Problem processing image ", n, nd2_image)
        return results
    return (get_results,)


@app.cell
def _(datasets, get_results):
    results_0 = get_results(datasets[0])
    return (results_0,)


@app.cell
def _(datasets, get_results):
    results_2 = get_results(datasets[2])
    return (results_2,)


@app.cell
def _(plt):
    def plot_results(results):
        fig = plt.figure(figsize=(12,3))
        ax = fig.subplots(1, len(results))
        for i in range(len(results)):
            print(i)
            ax[i].imshow(results[i]["theta_radii_means"])
        return fig
    return (plot_results,)


@app.cell
def _(plot_results, results_0):
    plot_results(results_0)
    return


@app.cell
def _(plot_results, results_2):
    plot_results(results_2)
    return


@app.cell
def _(data_path, nd2):
    test_image = nd2.imread(f"{data_path}/B06_250529_TRAK1-wt/Cell1-a.nd2", xarray=True)
    return (test_image,)


@app.cell
def _(plt, test_image):
    plt.imshow(test_image[:,0,:,:].max(axis=0))
    return


@app.cell
def _(plt, test_image):
    plt.imshow(test_image[:,1,:,:].max(axis=0))
    return


@app.cell
def _(plt, test_image):
    plt.imshow(test_image[:,2,:,:].max(axis=0))
    return


@app.cell
def _(plt, scale, test_image):
    plt.imshow(scale(test_image[:,3,:,:].max(axis=0)), vmax=0.001)
    return


@app.cell
def _(test_image):
    test_image
    return


@app.cell
def _(datasets):
    datasets
    return


@app.cell
def _(Path, datasets, nd2):
    def check_images():
        for dataset in datasets:
            nd2_images = [
                str(image) for image in Path(dataset).iterdir()
                if image.suffix == ".nd2"
            ]
            for nd2_image in nd2_images:
                image = nd2.imread(nd2_image, dask=True, xarray=True)
                print(nd2_image)
                print(image.C.coords)

    check_images()
    return


@app.cell
def _(image):
    v = image.C.coords
    return (v,)


@app.cell
def _(v):
    type(v)
    return


@app.cell
def _(datasets):
    datasets
    return


@app.cell
def _():
    return


if __name__ == "__main__":
    app.run()
