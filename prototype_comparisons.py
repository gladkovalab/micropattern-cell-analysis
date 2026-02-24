import marimo

__generated_with = "0.17.8"
app = marimo.App()


@app.cell
def _():
    import marimo as mo
    import matplotlib.pyplot as plt
    import numpy as np
    import xarray as xr
    import skimage.io
    from pathlib import Path
    from comparison_loader import (
        load_comparisons,
        find_well_directory,
        list_cells,
        load_cell,
        gather_datasets,
        aggregate_datasets,
        sum_channel,
        convert_to_uint16
    )
    return (
        Path,
        aggregate_datasets,
        convert_to_uint16,
        find_well_directory,
        gather_datasets,
        list_cells,
        load_cell,
        load_comparisons,
        mo,
        plt,
        skimage,
        sum_channel,
    )


@app.cell
def _(load_comparisons):
    dfs = load_comparisons()
    return (dfs,)


@app.cell
def _(dfs, mo):
    sheet_selector = mo.ui.dropdown(options=list(dfs.keys()), label="Select Sheet", value=list(dfs.keys())[0] if dfs else None)
    sheet_selector
    return (sheet_selector,)


@app.cell
def _(dfs, sheet_selector):
    selected_df = dfs[sheet_selector.value] if sheet_selector.value else None

    if selected_df is not None:
        # The first column is the plate name (often unnamed)
        plate_col = selected_df.columns[0]
        plates = selected_df[plate_col].to_list()
        # Other columns are conditions
        conditions = selected_df.columns[1:]
    else:
        plates = []
        conditions = []
        plate_col = None
    return conditions, plate_col, selected_df


@app.cell
def _(conditions, mo):
    condition_selector = mo.ui.dropdown(options=conditions, label="Select Condition")
    denoised_toggle = mo.ui.checkbox(label="Denoised", value=False)
    run_button = mo.ui.run_button(label="Aggregate Images")

    mo.hstack([condition_selector, denoised_toggle, run_button])
    return condition_selector, denoised_toggle, run_button


@app.cell
def _(
    aggregate_datasets,
    condition_selector,
    denoised_toggle,
    gather_datasets,
    mo,
    plate_col,
    run_button,
    selected_df,
):
    mo.stop(not run_button.value or not condition_selector.value)

    all_datasets = gather_datasets(selected_df, plate_col, condition_selector.value, denoised_toggle.value)

    aggregated_ds = aggregate_datasets(all_datasets)

    if aggregated_ds is None:
        mo.output.replace(mo.md("No cells found for this condition."))
    else:
        mo.output.replace(mo.md(f"**Aggregated {len(all_datasets)} cells.**"))
    return (all_datasets,)


@app.cell
def _(all_datasets):
    all_datasets
    return


@app.cell
def _(find_well_directory, list_cells):
    def list_datasets(selected_df, plate_col, condition, denoised):
        all_datasets = []
        # Iterate over all rows in the selected dataframe
        for row in selected_df.iter_rows(named=True):
            plate_name = row[plate_col]
            well_id = row[condition]

            if not well_id:
                continue

            well_dir = find_well_directory(plate_name, well_id)
            if well_dir:
                cells = list_cells(well_dir, denoised=denoised)
                for cell_path in cells:
                    try:
                        all_datasets.append(cell_path)
                    except Exception as e:
                        print(f"Error loading {cell_path}: {e}")
        return all_datasets
    return (list_datasets,)


@app.cell
def _(
    condition_selector,
    denoised_toggle,
    list_datasets,
    plate_col,
    selected_df,
):
    _list_datasets = list_datasets(selected_df, plate_col, condition_selector.value, denoised_toggle.value)
    _list_datasets
    return


@app.cell
def _(all_datasets, sum_channel):
    ch488 = sum_channel(all_datasets, "488")
    return (ch488,)


@app.cell
def _(ch488, plt):
    plt.imshow(ch488)
    return


@app.cell
def _(ch488):
    ch488.max() - ch488.min()
    return


@app.cell
def _(ch488, convert_to_uint16, plt):
    ch488_uint16 = convert_to_uint16(ch488)
    plt.imshow(ch488_uint16)
    return


@app.cell
def _(cell_selector, load_cell, mo):
    mo.stop(not cell_selector or not cell_selector.value)

    ds = load_cell(cell_selector.value)

    # Coordinates for channels
    channels = ds.coords['C'].values.tolist()
    channel_selector = mo.ui.dropdown(options=channels, label="Select Channel", value=channels[0] if channels else None)

    channel_selector
    return channel_selector, ds


@app.cell
def _(channel_selector, ds, mo, plt):
    mo.stop(not channel_selector.value)

    # The variable name in these NetCDFs is often '__xarray_dataarray_variable__'
    var_name = list(ds.data_vars)[0]
    data = ds[var_name].sel(C=channel_selector.value)

    fig, ax = plt.subplots(figsize=(8, 8))
    im = ax.imshow(data.values, cmap='gray')
    ax.set_title(f"{channel_selector.value}")
    plt.colorbar(im, ax=ax)

    mo.as_html(fig)
    return


@app.cell
def _(mo):
    batch_run_button = mo.ui.run_button(label="Generate All Comparison Projections")
    mo.vstack([
        mo.md("## Batch Processing"),
        mo.md("Click the button below to generate projections for all sheets and conditions."),
        batch_run_button
    ])
    return (batch_run_button,)


@app.cell
def _(
    Path,
    batch_run_button,
    convert_to_uint16,
    dfs,
    gather_datasets,
    mo,
    skimage,
    sum_channel,
):
    mo.stop(not batch_run_button.value)

    _output_dir = Path("comparison_projections")
    _output_dir.mkdir(exist_ok=True)

    _log = mo.output
    _log.append("")
    _log.append(f"Starting batch processing... Output: {_output_dir}")

    for _sheet_name, _df in dfs.items():
        _sheet_dir = _output_dir / _sheet_name
        _sheet_dir.mkdir(exist_ok=True)
        _log.append(f"Processing sheet: {_sheet_name}")

        _plate_col = _df.columns[0]
        _conditions = _df.columns[1:]

        for _condition in _conditions:
            # Clean condition name for filesystem
            _cond_safe = _condition.replace("/", "_").replace("\\", "_")
            _cond_dir = _sheet_dir / _cond_safe
            _cond_dir.mkdir(exist_ok=True)

            _datasets = gather_datasets(_df, _plate_col, _condition, denoised=False)

            if not _datasets:
                continue

            # 488
            _sum_488 = sum_channel(_datasets, "488")
            if _sum_488 is not None:
                _uint16_488 = convert_to_uint16(_sum_488)
                skimage.io.imsave(_cond_dir / "sum_488.tif", _uint16_488, check_contrast=False)

            # 405
            _sum_405 = sum_channel(_datasets, "405")
            if _sum_405 is not None:
                _uint16_405 = convert_to_uint16(_sum_405)
                skimage.io.imsave(_cond_dir / "sum_405.tif", _uint16_405, check_contrast=False)

    _log.append("Batch processing complete.")
    return


@app.cell
def _():
    return


if __name__ == "__main__":
    app.run()
