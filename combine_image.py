import marimo

__generated_with = "0.17.8"
app = marimo.App(width="medium")


@app.cell
def _():
    import xarray
    import netCDF4
    import pathlib
    import matplotlib.pyplot as plt
    import numpy as np
    import marimo as mo
    return mo, np, pathlib, plt, xarray


@app.cell
def _():
    #data = xarray.load_dataarray("projections/250710_patterned_plate_9_good/C03_250718_TRAK1/Cell1.nc")
    return


@app.cell
def _():
    return


@app.cell
def _(pathlib, root_path, stretch01, xarray):
    data = None
    for (dirpath, dirnames, filenames) in pathlib.Path(root_path).walk():
        for file in filenames:
            if "Denoised" not in file:
                print(file)
                if data is None:
                    data = stretch01(xarray.load_dataarray(dirpath / file).drop_indexes(["Y", "X"]))
                else:
                    data += stretch01(xarray.load_dataarray(dirpath / file).drop_indexes(["Y", "X"]))
    return data, dirpath


@app.cell
def _(data):
    data
    return


@app.cell
def _(data):
    data
    return


@app.cell
def _(data):
    data
    return


@app.cell
def _(np):
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
    return make_rgb, stretch01


@app.cell
def _(data, make_rgb, stretch01):
    rgb_data = make_rgb(
        stretch01(data.sel(C="488")),
        stretch01(data.sel(C="488")),
        stretch01(data.sel(C="405"))
    )
    return (rgb_data,)


@app.cell
def _():
    import pymupdf
    import io
    import cairosvg
    import skimage


    def get_template_at_width(width):
        file = pymupdf.open("single_pattern.ai")
        png_bytes = cairosvg.svg2png(file[0].get_svg_image(), output_width=width, output_height=width)
        with io.BytesIO() as buf:
            buf.write(png_bytes)
            template_img = skimage.io.imread(buf)
        return template_img
    return get_template_at_width, skimage


@app.cell
def _(get_template_at_width):
    template = get_template_at_width(1326)[:,:,0]
    return (template,)


@app.cell
def _(skimage, template):
    template_contour = skimage.measure.find_contours(template)[0]
    template_contour -= 1326//2
    template_contour += 512
    return (template_contour,)


@app.cell
def _():
    #root_path = pathlib.Path("projections/250710_patterned_plate_9_good/C03_250718_TRAK1/")
    #root_path = pathlib.Path("projections/250710_patterned_plate_9_good/C05_250718_TRAK1_mDRH/")
    #root_path = pathlib.Path("projections/250710_patterned_plate_9_good/C06_250718_TRAK1_mDRH_dSp/")
    return


@app.cell(hide_code=True)
def _(mo, pathlib):
    file_browser = mo.ui.file_browser(
        initial_path=pathlib.Path("projections"),
        selection_mode = "directory",
        multiple=False
    )
    file_browser
    return (file_browser,)


@app.cell
def _(data, plt, stretch01, template_contour):
    plt.imshow(stretch01(data.sel(C="488")))
    plt.plot(template_contour[:,1], template_contour[:,0], color="white")
    return


@app.cell
def _(plt, rgb_data, template_contour):
    plt.imshow(rgb_data)
    plt.plot(template_contour[:,1], template_contour[:,0], color="white")
    return


@app.cell
def _(data, plt, stretch01, template_contour):
    plt.imshow(stretch01(data.sel(C="405")))
    plt.plot(template_contour[:,1], template_contour[:,0], color="white")
    return


@app.cell
def _(file_browser):
    root_path = file_browser.value[0].path
    return (root_path,)


@app.cell
def _(data):
    data
    return


@app.cell
def _(data, plt):
    plt.imshow(data.isel(C="488"))
    return


@app.cell
def _(data, plt):
    plt.imshow(data.isel(C=1))
    return


@app.cell
def _(np, plt, template_contour):
    fig = plt.figure(figsize=(10.24, 10.24), dpi=100)
    ax = plt.Axes(fig, [0., 0., 1., 1.])
    ax.set_axis_off()
    fig.add_axes(ax)
    ax.imshow(np.zeros((1024, 1024)), cmap="gray")
    ax.plot(template_contour[:, 1], template_contour[:, 0], color="white")
    fig.savefig("template_1024x1024.tiff", dpi=100)
    return


@app.cell
def _(dirpath, plt, xarray):
    plt.imshow(xarray.load_dataarray(dirpath / "Cell9 - Denoised.nc").drop_indexes(["Y", "X"]).sel(C="405"))
    return


@app.cell
def _():
    return


if __name__ == "__main__":
    app.run()
