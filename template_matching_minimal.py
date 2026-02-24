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
    return


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
def _(get_template_at_width, np):
    def get_padded_template_at_width(template_width, *, base_template=None):
        if base_template is None:
            base_template = get_template_at_width(template_width)[:,:,0]
        pad = (2048-template_width)//2
        template = np.pad(base_template,(pad, pad))
        return template
    return (get_padded_template_at_width,)


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
def _(nd2):
    #img_path = "/groups/vale/valelab/_for_Mark/patterned_data/250521_patterned_plate_1/B06_250528_TRAK1-wt/Cell3.nd2"
    #img_path = "/groups/vale/valelab/_for_Mark/patterned_data/250521_patterned_plate_1/B06_250528_TRAK1-wt/Cell8.nd2"
    #img_path = "/groups/vale/valelab/_for_Mark/patterned_data/250710_patterned_plate_9_good/C02_250718_NoV/cell8.nd2"
    #img_path = "/groups/vale/valelab/_for_Mark/patterned_data/250612_patterned_plate_3/B06_250617_TRAK1_mDRH_dSp/denoised/Cell12 - Denoised.nd2"
    img_path = "/groups/vale/valelab/_for_Mark/patterned_data/250710_patterned_plate_9_good/G09_250718_MAPK9_siRNA_Ars/denoised/Cell9 - Denoised.nd2"
    img = nd2.imread(img_path, xarray=True)
    return img, img_path


@app.cell
def _(img):
    img
    return


@app.cell
def _(img, np, plt):
    plt.imshow(np.sum(img.sel(C="640"), axis=0), vmin=0, vmax=10000)
    # 1082, 407
    plt.scatter(1082, 407)
    return


@app.cell
def _(get_padded_template_at_width, plt):
    shifted_template2 = get_padded_template_at_width(1326)
    #shifted_template2 = np.roll(shifted_template, (0 - 1024,0 - 1024), axis=(0,1))
    plt.imshow(shifted_template2)
    plt.scatter(1024,1024)
    return (shifted_template2,)


@app.cell
def _(shifted_template2):
    shifted_template2.shape
    return


@app.cell
def _(np):
    def get_image_hat(img):
        img_template_sum_projection = np.sum(img.sel(C="640"), axis=0)
        img_template_sum_projection_norm = img_template_sum_projection / np.max(img_template_sum_projection)
        #img_template_sum_projection_hat = np.abs(np.fft.fft2(img_template_sum_projection))
        img_template_sum_projection_norm_2048 = img_template_sum_projection_norm[:2048,:2048]
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
    def max_match_template(img, *, template_hat = None):
        template_matching = match_template(img, template_hat = template_hat)
        max_idx = np.argmax(template_matching)
        return np.unravel_index(max_idx, template_matching.shape)
    return (max_match_template,)


@app.cell
def _(img_path, match_template):
    #template_matching = np.fft.fftshift(np.real(np.fft.ifft2(template_hat * img_template_hat)))
    template_matching = match_template(img_path)
    return (template_matching,)


@app.cell
def _(plt, template_matching):
    plt.imshow(template_matching)
    return


@app.cell
def _(img_path, max_match_template):
    max_coords = max_match_template(img_path)
    return (max_coords,)


@app.cell
def _(get_padded_template_at_width, max_coords, np, plt):
    shifted_template = get_padded_template_at_width(1326)
    shifted_template = np.roll(shifted_template, (max_coords[0] - 1024,max_coords[1] - 1024), axis=(0,1))
    plt.imshow(shifted_template)
    plt.scatter(max_coords[1], max_coords[0])
    return (shifted_template,)


@app.cell
def _(img, np, plt):
    sumproj = np.sum(img[:,1,:2048,:2048], axis=0) > 6100
    plt.imshow(sumproj)
    return (sumproj,)


@app.cell
def _(plt, shifted_template, sumproj):
    plt.imshow(shifted_template)
    plt.imshow(sumproj, alpha=0.5)
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
def _(np, shifted_template, sumproj):
    score = (np.sum(sumproj & shifted_template)/(np.sum(shifted_template > 0))).values.item()
    score
    return (score,)


@app.cell
def _(plt, shifted_template, sumproj):
    plt.imshow(sumproj & shifted_template)
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
def _(img):
    img.sel(C="405")
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
    img4 = nd2.imread("/groups/vale/valelab/_for_Mark/patterned_data/250731_patterned_plate_11_good/E05_250808_TRAK1_wt/denoised/Cell2 - Denoised.nd2", xarray=True)
    return (img4,)


@app.cell
def _(img4, np):
    np.sum(img4, axis=0)
    return


@app.cell
def _(img, plt):
    plt.imshow(img.sel(C="488").sum(axis=0), vmax=8000)
    return


@app.cell
def _():
    return


if __name__ == "__main__":
    app.run()
