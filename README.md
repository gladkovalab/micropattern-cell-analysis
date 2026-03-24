# Micropattern Cell Analysis

Python pipeline for quantifying mitochondrial distribution in fluorescence microscopy images of micropatterned cells. The pipeline localizes fibronectin micropatterns via template matching, segments nuclei, and computes perinuclear vs. peripheral mitochondrial signal ratios.

## Overview

Cells are cultured on fibronectin micropatterns that constrain cell shape. This pipeline:

1. Localizes the micropattern center in each ND2 image via FFT-based template matching
2. Extracts a 1024×1024 px crop centered on the pattern
3. Segments the nucleus (405 nm / DAPI channel)
4. Computes Euclidean distance transforms from the nuclear boundary and pattern arch
5. Subtracts background from the mitochondrial channel (488 nm)
6. Quantifies mitochondrial signal in perinuclear (<5 µm from nucleus) and peripheral (<5 µm from pattern arch) zones
7. Outputs per-cell metrics to CSV/Excel and diagnostic figures to PDF

An optional second stage aggregates per-cell projections into mean images per experimental condition for visual comparison.

## Installation

Dependencies are managed with [Pixi](https://pixi.sh). With Pixi installed:

```bash
pixi install
```

This installs Python 3.12 and all required packages into an isolated environment on Linux or macOS.

## Usage

### Batch analysis

Process all ND2 files under a root directory:

```bash
pixi run python template_matching_bulk.py /path/to/patterned_data/plate_name
```

**Options:**

| Flag | Description |
|------|-------------|
| `--keep-sums` | Retain raw signal sum columns (default: dropped) |
| `--include-complex` | Include non-simple peripheral measurements |
| `--include-acute` | Include acute peripheral zone measurements |
| `--include-all-percents` | Include all percentage columns (default: `_total` only) |

**Outputs** (written relative to the working directory):

- `projections/{plate}/{well}/*.nc` — per-cell cropped z-sum projections (NetCDF4)
- `projections/{plate}/{well}/*_488_bg_subtracted.nc` — background-subtracted mitochondrial projections
- `template_matching/{plate}/{well}/template_matching.csv` — per-cell metrics
- `template_matching/{plate}/{well}/template_matching.xlsx` — same data as Excel
- `template_matching/{plate}/{well}/*.pdf` — diagnostic figures per cell

### Cluster submission

For large datasets, `bsub_analysis.sh` submits one LSF job per directory using the paths listed in `config/20251229_paths_for_analysis.txt`:

```bash
bash bsub_analysis.sh
```

To target a different set of directories, replace the contents of `config/20251229_paths_for_analysis.txt` with the desired paths (one per line), or edit `bsub_analysis.sh` directly. Each job requests 8 cores via `bsub -n 8 -P vale`.

### Comparison projections

After running the batch analysis, generate mean projection images per experimental condition:

```bash
pixi run python generate_comparison_projections.py
```

This reads condition/well assignments from `config/Comparisons_table_v3.xlsx` and writes mean TIFF images to `comparison_projections/`.

### Interactive notebooks

The Marimo notebooks can be launched for interactive exploration:

```bash
pixi run marimo edit micropattern_cell_analysis_viewer.py     # browse raw ND2 files
pixi run marimo edit prototype_comparisons.py                 # explore condition comparisons
pixi run marimo edit micropattern_cell_analysis.py            # single-cell analysis
```

Or launch the viewer via the configured Pixi task:

```bash
pixi run notebook
```

## Script Reference

| Script | Description |
|--------|-------------|
| `template_matching_bulk.py` | **Main pipeline.** Batch processes ND2 files: template matching, cropping, nuclear segmentation, background subtraction, zone quantification, output to NetCDF4/CSV/Excel/PDF. |
| `comparison_loader.py` | Utility module. Provides `load_comparisons()`, `find_well_directory()`, `list_cells()`, `load_cell()`, and channel aggregation helpers used by other scripts. |
| `generate_comparison_projections.py` | Aggregates per-cell projections into mean images per experimental condition. Handles standard, denoised, and background-subtracted variants. |
| `combine_analysis.py` | Merges per-well `template_matching.csv` files into a single summary Excel workbook. |
| `plot_comparisons.py` | Plots template matching scores and per-condition metrics from CSV output. |
| `micropattern_cell_analysis.py` | Marimo notebook for interactive single-cell analysis and visualization. |
| `micropattern_cell_analysis_batch.py` | Marimo notebook version of the batch pipeline with enhanced visualization. |
| `micropattern_cell_analysis_viewer.py` | Marimo notebook for browsing ND2 files with channel and z-slice selection. |
| `micropattern_cell_analysis_pattern_center.py` | Marimo notebook for validating and manually setting pattern center coordinates. |
| `prototype_comparisons.py` | Marimo notebook for comparing projections across experimental conditions. |

## Configuration

### Coordinate overrides

When template matching fails or requires manual correction, add entries to `coordinate_overrides.csv` (no header):

```
/path/to/Cell1.nd2,x_pixels,y_pixels
```

The `x` and `y` values specify the top of the pattern (not the center); the pipeline converts these to pattern center coordinates automatically.

### Offset overrides

For images where the pattern is not near the default image offset, add entries to the `offset_overrides` dict in `template_matching_bulk.py` (lines 80–89). The default offset is `[128, 128]` pixels.

### ROI overrides

To restrict template matching to a sub-region of an image, add entries to the `roi_overrides` dict in `template_matching_bulk.py` (lines 137–139).

## Data Organization

**Input** (ND2 files):
```
/path/to/patterned_data/
└── {plate_name}/
    └── {well_id}_{condition}/
        ├── Cell1.nd2
        ├── Cell2.nd2
        └── ...
```

**Output:**
```
projections/
└── {plate_name}/{well_id}_{condition}/
    ├── Cell1.nc                        # all-channel cropped projection
    ├── Cell1_488_bg_subtracted.nc      # background-subtracted mito channel
    └── denoised/                       # denoised variants (if acquired)

template_matching/
└── {plate_name}/{well_id}_{condition}/
    ├── Cell1.pdf                       # diagnostic figures
    ├── template_matching.csv
    └── template_matching.xlsx

comparison_projections/
└── {condition}/
    ├── mean_488.tif                    # mean mitochondrial projection
    ├── mean_405.tif                    # mean nuclear projection
    └── individual_wells/
        └── {plate}_{well}_mean_*.tif
```

## Output Metrics

Each row in `template_matching.csv` corresponds to one cell and includes:

| Column | Description |
|--------|-------------|
| `path` | Path to the source ND2 file |
| `template_matching_score` | Fractional overlap of Otsu-thresholded image with the pattern template (quality control) |
| `cropped_background_threshold` | Background level estimated from image edges |
| `peripheral_{d}um_simple_percent_total` | Peripheral mitochondrial signal within *d* µm of arch / total crop signal (%) |
| `perinuclear_{d}um_percent_total` | Perinuclear mitochondrial signal within *d* µm of nucleus / total crop signal (%) |

Distances *d* = 1, 2, 3, 4, 5 µm are computed for each zone.

## Methods

Detailed methods suitable for a Methods section are provided in:

- [`methods/2026_03_06_methods.md`](methods/2026_03_06_methods.md) — template matching, segmentation, and quantification
- [`methods/2026_03_06_methods_comparison_projections.md`](methods/2026_03_06_methods_comparison_projections.md) — comparison projection generation
- [`methods/2026_03_06_references.md`](methods/2026_03_06_references.md) — package citations

## Dependencies

Key packages (see `pixi.toml` for full list):

| Package | Purpose |
|---------|---------|
| `nd2` | Reading Nikon ND2 microscopy files |
| `numpy` / `scipy` | FFT-based template matching, distance transforms |
| `scikit-image` | Otsu thresholding, contour finding, region properties |
| `xarray` / `netCDF4` | Labeled image arrays, projection file I/O |
| `polars` | DataFrame operations and Excel output |
| `matplotlib` | Diagnostic figure generation |
| `cairosvg` / `pymupdf` | Rasterizing the SVG micropattern template |
| `marimo` | Reactive Python notebooks for interactive exploration |

## Authors

Mark Kittisopikul — [kittisopikulm@janelia.hhmi.org](mailto:kittisopikulm@janelia.hhmi.org)

## License

Copyright © 2025 Howard Hughes Medical Institute

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

- Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
- Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
- Neither the name of HHMI nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

See [`LICENSE.md`](LICENSE.md) for the full license text.
