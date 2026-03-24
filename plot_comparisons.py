import polars as pl
import matplotlib.pyplot as plt
import numpy as np
import pathlib
from pathlib import Path
from comparison_loader import load_comparisons
import matplotlib.colors as mcolors

TEMPLATE_MATCHING_DIR = Path("template_matching")

def find_template_matching_csv(plate_name: str, well_id: str) -> Path | None:
    """Finds the template_matching.csv matching the plate and well ID."""
    plate_dir = TEMPLATE_MATCHING_DIR / plate_name
    if not plate_dir.exists():
        return None
    
    # Well ID is usually the start of the directory name (e.g., B06_...)
    for subdir in plate_dir.iterdir():
        if subdir.is_dir() and subdir.name.startswith(well_id):
            csv_path = subdir / "template_matching.csv"
            if csv_path.exists():
                return csv_path
    return None

def get_shades(base_color, n):
    rgb = mcolors.to_rgb(base_color)
    hsv = mcolors.rgb_to_hsv(rgb)
    shades = []
    if n == 1:
        return [rgb]
    for i in range(n):
        # Vary value/brightness and saturation a bit
        s = hsv[1] * (0.6 + 0.4 * (i / (n - 1)))
        v = hsv[2] * (0.6 + 0.4 * (i / (n - 1)))
        shades.append(mcolors.hsv_to_rgb((hsv[0], s, v)))
    return shades

def plot_metric(dfs, metric_name, output_dir):
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    base_colors = list(mcolors.TABLEAU_COLORS.values())

    for sheet_name, df in dfs.items():
        print(f"Plotting sheet: {sheet_name}")
        plate_col = df.columns[0]
        conditions = df.columns[1:]
        
        fig, ax = plt.subplots(figsize=(12, 6))
        
        all_x_labels = []
        
        for cond_idx, condition in enumerate(conditions):
            base_color = base_colors[cond_idx % len(base_colors)]
            
            # Gather all data for this condition
            condition_data = [] # List of (well_mean, well_cells, shade)
            
            wells_found = []
            for row in df.iter_rows(named=True):
                plate_name = row[plate_col]
                well_id = row[condition]
                if not well_id:
                    continue
                
                csv_path = find_template_matching_csv(plate_name, well_id)
                if csv_path:
                    try:
                        # We need to read the CSV. 
                        # Note: template_matching_bulk.py might have multiple rows per CSV if it walked subdirs, 
                        # but usually it's one CSV per well directory.
                        well_df = pl.read_csv(csv_path)
                        if metric_name in well_df.columns:
                            values = well_df[metric_name].drop_nans().to_numpy()
                            if len(values) > 0:
                                wells_found.append(values)
                    except Exception as e:
                        print(f"  Error reading {csv_path}: {e}")

            if not wells_found:
                continue

            n_wells = len(wells_found)
            shades = get_shades(base_color, n_wells)
            
            x_base = cond_idx
            
            all_well_means = []
            all_condition_cells = []
            
            for i, well_cells in enumerate(wells_found):
                shade = shades[i]
                well_mean = np.mean(well_cells)
                all_well_means.append(well_mean)
                all_condition_cells.extend(well_cells)
                
                # Jitter for individual cells
                jitter = np.random.normal(0, 0.05, size=len(well_cells))
                ax.scatter(np.full_like(well_cells, x_base) + jitter, well_cells, 
                           color=shade, alpha=0.2, s=20, edgecolors='none')
                
                # Well mean
                ax.scatter([x_base], [well_mean], color=shade, alpha=1.0, s=100, marker='o', edgecolors='white', zorder=3)

            # Overall condition stats
            cond_mean = np.mean(all_well_means)
            cond_sem = np.std(all_well_means) / np.sqrt(len(all_well_means)) if len(all_well_means) > 1 else 0
            
            ax.errorbar(x_base, cond_mean, yerr=cond_sem, fmt='none', ecolor='black', capsize=10, elinewidth=2, zorder=4)
            ax.hlines(cond_mean, x_base - 0.2, x_base + 0.2, colors='black', linewidth=3, zorder=4)
            
            all_x_labels.append(condition)

        ax.set_xticks(range(len(all_x_labels)))
        ax.set_xticklabels(all_x_labels, rotation=45, ha='right')
        ax.set_title(f"{sheet_name} - {metric_name}")
        ax.set_ylabel(metric_name)
        
        plt.tight_layout()
        plt.savefig(output_dir / f"{sheet_name}_{metric_name}.pdf")
        plt.savefig(output_dir / f"{sheet_name}_{metric_name}.png")
        plt.close()

def main():
    try:
        dfs = load_comparisons()
    except Exception as e:
        print(f"Could not load comparisons: {e}")
        return

    metrics = [
        "peripheral_1um_simple_percent_total",
        "peripheral_2um_simple_percent_total",
        "peripheral_3um_simple_percent_total",
        "peripheral_4um_simple_percent_total",
        "peripheral_5um_simple_percent_total",
        "perinuclear_1um_percent_total",
        "perinuclear_2um_percent_total",
        "perinuclear_3um_percent_total",
        "perinuclear_4um_percent_total",
        "perinuclear_5um_percent_total"
    ]
    
    output_dir = Path("comparison_plots")
    
    for metric in metrics:
        print(f"Processing metric: {metric}")
        plot_metric(dfs, metric, output_dir)

if __name__ == "__main__":
    main()
