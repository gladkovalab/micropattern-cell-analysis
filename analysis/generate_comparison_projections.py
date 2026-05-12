from pathlib import Path
import skimage.io
import numpy as np
from comparison_loader import (
    load_comparisons,
    find_well_directory,
    list_cells,
    load_cell,
    sum_channel,
    convert_to_uint16
)

def save_mean_image(datasets, channel, output_path):
    n = len(datasets)
    if n == 0:
        return
    mean_img = sum_channel(datasets, channel)
    if mean_img is not None:
        mean_img /= n
        uint16_img = convert_to_uint16(mean_img, stretch=False)
        skimage.io.imsave(output_path, uint16_img, check_contrast=False)

def process_datasets(well_to_datasets, output_dir, prefix=""):
    # Condition-wide
    all_datasets = []
    for datasets in well_to_datasets.values():
        all_datasets.extend(datasets)
    
    if not all_datasets:
        return

    # Channels to process
    channels = ["488", "405"]
    for ch in channels:
        suffix = f"_{prefix}" if prefix else ""
        save_mean_image(all_datasets, ch, output_dir / f"mean_{ch}{suffix}.tif")

    # Per-well
    wells_dir = output_dir / "individual_wells"
    wells_dir.mkdir(exist_ok=True)
    for (plate, well), datasets in well_to_datasets.items():
        well_prefix = f"{plate}_{well}"
        for ch in channels:
            suffix = f"_{prefix}" if prefix else ""
            save_mean_image(datasets, ch, wells_dir / f"{well_prefix}_mean_{ch}{suffix}.tif")

def gather_well_data(df, plate_col, condition, denoised, pattern="*.nc", exclude_bg_subtracted=False):
    well_to_datasets = {}
    for row in df.iter_rows(named=True):
        plate_name = row[plate_col]
        well_id = row[condition]
        if not well_id:
            continue
        
        well_dir = find_well_directory(plate_name, well_id)
        if well_dir:
            cells = list_cells(well_dir, denoised=denoised, pattern=pattern)
            if exclude_bg_subtracted:
                cells = [c for c in cells if "_bg_subtracted" not in c.name]
            
            datasets = []
            for cell_path in cells:
                try:
                    ds = load_cell(cell_path)
                    datasets.append(ds)
                except Exception as e:
                    print(f"Error loading {cell_path}: {e}")
            
            if datasets:
                well_to_datasets[(plate_name, well_id)] = datasets
    return well_to_datasets

def main():
    dfs = load_comparisons()
    output_base = Path("comparison_projections")
    output_base.mkdir(exist_ok=True)
    
    print(f"Starting batch processing... Output: {output_base}")
    
    for sheet_name, df in dfs.items():
        sheet_dir = output_base / sheet_name
        sheet_dir.mkdir(exist_ok=True)
        print(f"Processing sheet: {sheet_name}")
        
        plate_col = df.columns[0]
        conditions = df.columns[1:]
        
        for condition in conditions:
            cond_safe = condition.replace("/", "_").replace("\\", "_")
            cond_dir = sheet_dir / cond_safe
            cond_dir.mkdir(exist_ok=True)
            print(f"  Condition: {condition}")

            # 1. Standard projections
            well_data = gather_well_data(df, plate_col, condition, denoised=False, pattern="*.nc", exclude_bg_subtracted=True)
            process_datasets(well_data, cond_dir)

            # 2. Denoised standard projections
            denoised_dir = cond_dir / "denoised"
            denoised_dir.mkdir(exist_ok=True)
            well_data_denoised = gather_well_data(df, plate_col, condition, denoised=True, pattern="*.nc", exclude_bg_subtracted=True)
            process_datasets(well_data_denoised, denoised_dir)

            # 3. BG Subtracted (Newly created files)
            # These are typically in the main well dir, not denoised (based on template_matching_bulk.py)
            well_data_bg = gather_well_data(df, plate_col, condition, denoised=False, pattern="*_bg_subtracted.nc")
            process_datasets(well_data_bg, cond_dir, prefix="bg_subtracted")

            # 4. Denoised BG Subtracted
            well_data_bg_denoised = gather_well_data(df, plate_col, condition, denoised=True, pattern="*_bg_subtracted.nc")
            process_datasets(well_data_bg_denoised, denoised_dir, prefix="bg_subtracted")

    print("Batch processing complete.")

if __name__ == "__main__":
    main()