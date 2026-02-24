import polars as pl
import fastexcel
import xarray as xr
import numpy as np
from pathlib import Path

PROJECTIONS_DIR = Path("projections")

def load_comparisons(file_path: str | Path = "/groups/vale/valelab/_for_Mark/analysis/Comparisons_table_v3.xlsx") -> dict[str, pl.DataFrame]:
    """Reads all sheets from the comparisons Excel table using fastexcel and polars."""
    excel_reader = fastexcel.read_excel(file_path)
    sheets = {}
    for sheet_name in excel_reader.sheet_names:
        sheets[sheet_name] = pl.read_excel(file_path, sheet_name=sheet_name)
    return sheets

def find_well_directory(plate_name: str, well_id: str) -> Path | None:
    """Finds the directory matching the plate and well ID."""
    plate_dir = PROJECTIONS_DIR / plate_name
    if not plate_dir.exists():
        return None
    
    # Well ID is usually the start of the directory name (e.g., B06_...)
    for subdir in plate_dir.iterdir():
        if subdir.is_dir() and subdir.name.startswith(well_id):
            return subdir
    return None

def list_cells(well_dir: Path, denoised: bool = False, pattern: str = "*.nc") -> list[Path]:
    """Lists .nc files in the well directory or its denoised subdirectory matching a pattern."""
    target_dir = well_dir / "denoised" if denoised else well_dir
    if not target_dir.exists():
        return []
    return sorted(list(target_dir.glob(pattern)))

def load_cell(file_path: Path) -> xr.Dataset:
    """Loads a cell projection NetCDF file."""
    return xr.open_dataset(file_path)

def gather_datasets(selected_df, plate_col, condition, denoised, pattern="*.nc"):
    all_datasets = []
    # Iterate over all rows in the selected dataframe
    for row in selected_df.iter_rows(named=True):
        plate_name = row[plate_col]
        well_id = row[condition]

        if not well_id:
            continue

        well_dir = find_well_directory(plate_name, well_id)
        if well_dir:
            cells = list_cells(well_dir, denoised=denoised, pattern=pattern)
            for cell_path in cells:
                try:
                    loaded_ds = load_cell(cell_path)
                    all_datasets.append(loaded_ds)
                except Exception as e:
                    print(f"Error loading {cell_path}: {e}")
    return all_datasets

def aggregate_datasets(all_datasets):
    aggregated = None
    first_da = None

    for dataset in all_datasets:
        da = dataset.to_dataarray().sum(axis=0)
        if aggregated is None:
            aggregated = da.values.astype(np.float64)
            first_da = da
        else:
            aggregated += da.values

    if first_da is not None:
        result = first_da.copy()
        result.values = aggregated
        return result
    return None

def sum_channel(all_datasets, channel="488"):
    if not all_datasets:
        return None

    arrays = []
    max_h, max_w = 0, 0

    for dataset in all_datasets:
        try:
            # Replicate the logic: sum axis 0 (likely Z or variables), then select channel
            # If it's a Dataset with one variable, to_dataarray() adds a 'variable' dimension.
            da = dataset.to_dataarray().sum(axis=0)
            if "C" in da.coords and da.coords["C"].values.size == 1:
                 # Scalar coordinate
                 if str(da.coords["C"].values) != channel:
                      continue
                 arr = da.to_numpy()
            elif "C" in da.dims:
                 arr = da.sel(C=channel).to_numpy()
            else:
                 # No C dimension/coord, assume it matches if we only have one?
                 # Actually better to be safe.
                 continue

            if np.issubdtype(arr.dtype, np.floating):
                 pass
                 #arr = arr * 65535.0

            arrays.append(arr)
            
            h, w = arr.shape
            max_h = max(max_h, h)
            max_w = max(max_w, w)
        except Exception:
            # Channel might not exist or other error
            continue

    if not arrays:
        return None

    ch = np.zeros((max_h, max_w), dtype=np.float64)

    for arr in arrays:
        h, w = arr.shape
        # Center alignment
        y_off = (max_h - h) // 2
        x_off = (max_w - w) // 2
        ch[y_off:y_off+h, x_off:x_off+w] += arr

    return ch

def convert_to_uint16(arr, stretch=False):
    if stretch:
        _min = arr.min()
        _max = arr.max()
        if _max == _min:
            return np.zeros_like(arr, dtype=np.uint16)
        return ((arr - _min) * (np.iinfo(np.uint16).max / (_max - _min))).astype(np.uint16)
    else:
        return arr.astype(np.uint16)

if __name__ == "__main__":
    dfs = load_comparisons()
    # Test finding a directory
    plate = "250521_patterned_plate_1"
    well = "B06"
    wdir = find_well_directory(plate, well)
    print(f"Well dir: {wdir}")
    if wdir:
        cells = list_cells(wdir)
        print(f"First cell: {cells[0] if cells else 'None'}")
