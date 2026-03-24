import marimo

__generated_with = "0.17.8"
app = marimo.App(
    width="medium",
    app_title="Micropattern Cell Analysis",
    layout_file="layouts/micropattern_cell_analysis_viewer.slides.json",
)


@app.cell
def _():
    # Configuration
    # Initial starting path for the file browser (CHANGE ME)
    initial_path_str = "/groups/vale/valelab/_for_Mark/patterned_data/250521_patterned_plate_1/B06_TRAK1_wt_combined/"
    return (initial_path_str,)


@app.cell(hide_code=True)
def _(Path, initial_path_str, mo):
    file_browser = mo.ui.file_browser(
        initial_path=Path(initial_path_str),
        selection_mode = "file",
        multiple=False
    )
    file_browser
    return (file_browser,)


@app.cell(hide_code=True)
def _(
    channel_dropdown,
    image_scale_slider,
    imshow_cz,
    mo,
    selection,
    z_slider,
):
    def show_viewer():
        if len(selection) == 0:
            print("Please select a file")
            return None
        else:
            return mo.vstack([
                mo.hstack([
                    imshow_cz(),
                    image_scale_slider
                ]),
                mo.hstack([channel_dropdown,z_slider])
            ])

    show_viewer()
    return


@app.cell(hide_code=True)
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
    return Path, mo, nd2, np, plt


@app.cell(hide_code=True)
def _():
    # do not worry about clustered and dispersed
    # focus on patterned_data
    # 250710_patterned_plate_9_good, C2, C3, C4
    # Take any nd2 file, do not look at excluded cell folder,
    # cell_filtering.txt is informational only, no need to parse it
    # ignore tile.nd2, those tiles are big
    # only look for nd2 files that start with cell
    # NoV = No virus
    # _for_Mark is a duplicate of raw data
    # D05 was imaged on multiple days
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    https://resisted-curiosity-682.notion.site/Final-folder-structure-26879054849480ac8473e5427c072825
    """)
    return


@app.cell(hide_code=True)
def _():
    folder_structure = """
    - no TRAK / TRAK1 / TRAK2 - mito
        - 250612_patterned_plate_3 (B02 / B03 / B04)
        - 250710_patterned_plate_9 (C02 (n=8) / C03 (n=8) / C04 (n=11))
        - 250731_patterned_plate_11 (D06 / E05 / F05)
    - no TRAK / TRAK1 / TRAK2 - peroxisome
        - 250606_patterned_plate_2 (G10 / F06 / G06)
        - 250612_patterned_plate_3 (E02 / E03 / E04)
        - 250626_patterned_plate_7 (F02 / G03 / F04)
        - 250710_patterned_plate_9 (D02 / D03 / D04)
    - no TRAK / TRAK1 / TRAK2 - 60mers
        - 250620_patterned_plate_5 (E02 / E04 / E03) - cGP80s
        - 250624_patterned_plate_6 (E07 / E08 / E09) - cGP200s (better)
    - TRAK2: wt / DRH / DRH+Spindly
        - 250606_patterned_plate_2 (C02 / E04 / D05)
        - 250612_patterned_plate_3 (B04 / B07 / B08)
        - 250710_patterned_plate_9 (C04 / C07 / C08)
        - 250731_patterned_plate_11 (F05 / D07 / E07)
    - TRAK1: wt / DRH / DRH+Spindly
        - 250521_patterned_plate_1 (B06 / E06 / D07) - pilot / may not have enough cells
        - 250612_patterned_plate_3 (B03 / B05 / B06)
        - 250710_patterned_plate_9 (C03 / C05 / C06)
        - 250731_patterned_plate_11 (E05 / E06 / F06)
    - TRAK2: wt / S84A / S84E
        - 250606_patterned_plate_2 (C02 / B05 / C05)
        - 250612_patterned_plate_3 (B04 / B09 / B10)
        - 250710_patterned_plate_9 (C04 / C09 / C10)
        - 250731_patterned_plate_11 (F05 / F07 / D08)
    - TRAK2: wt -/+ Ars / S84A -/+ Ars # Christina needs to double check this
        - 250626_patterned_plate_7 (B02 / B03 / B09 / B10) -  pilot / may not have enough
        - 250710_patterned_plate_9 (F02 / F03 / F08 / F09)
        - 250731_patterned_plate_11 (E04 / F04 / B04 / C04)
        - 250807_patterned_plate_12 (F05 / G05 / B05 / C05))
    - wt cells: ctrl siRNA -/+ Ars / MAPK9 siRNA -/+ Ars
        - 250710_patterned_plate_9 (F05 / G03 / F11 / G09)
        - 250724_patterned_plate_10 (D05 / E05 / E02 / B03)
        - 250731_patterned_plate_11 ( E03 / F03 / B03 / C03)
        - 250807_patterned_plate_12 (F04 / G04 / B04 / C04)
    - wt cells: +/- Ars
        - 250618_patterned_plate_4 (B02 / B08)
        - 250624_patterned_plate_6 (B06 / D03)
        - 250701_patterned_plate_8 (B06 / D06)
    """
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    * Metric between clustered and dispersed mitochondira
    * Get probability density distribution between the two
    * Sum across the intensities
    * 60mers do not have a triplicate
    * Just extract the density map
    """)
    return


@app.cell(hide_code=True)
def _(file_browser):
    selection = file_browser.value
    selection
    return (selection,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    https://resisted-curiosity-682.notion.site/Micropatterned-cell-analysis-1fc79054849480e887f6d45ba3aeecfb
    """)
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
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
    """)
    return


@app.cell(hide_code=True)
def _(nd2, selection):
    if len(selection) > 0:
        image_path = selection[0].path
        image = nd2.imread(image_path, xarray=True, dask=True)
    return image, image_path


@app.cell(hide_code=True)
def _(np):
    def scale(arr):
        min = np.min(arr)
        max = np.max(arr)
        return (arr - min)/(max-min)
    return (scale,)


@app.cell(hide_code=True)
def _(image, mo, selection):
    if len(selection) > 0:
        channel_dropdown = mo.ui.dropdown(
            options=[str(c) for c in image.C.values],
            value=image.C.values[0],
            label="Channel"
        )
        channel_dropdown
    return (channel_dropdown,)


@app.cell(hide_code=True)
def _(image, mo, selection):
    if len(selection) > 0:
        z_slider = mo.ui.slider(steps=image.Z.values, full_width=True, label="Z")
    return (z_slider,)


@app.cell(hide_code=True)
def _(mo, selection):
    if len(selection) > 0:
        image_scale_slider = mo.ui.range_slider(
            orientation="vertical",
            start=0.0,
            stop=1.0,
            step=0.01,
            full_width=True,
            show_value=True)
    return (image_scale_slider,)


@app.cell(hide_code=True)
def _(channel_dropdown, image, selection, z_slider):
    if len(selection) > 0:
        image_CZ = image.sel(C=channel_dropdown.selected_key, Z=z_slider.value)
    return (image_CZ,)


@app.cell(hide_code=True)
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
def _(image_path):
    image_path
    return


@app.cell
def _(image):
    image.Z
    return


@app.cell
def _():
    return


if __name__ == "__main__":
    app.run()
