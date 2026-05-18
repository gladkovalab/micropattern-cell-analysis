# Micropattern Cell Analysis

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.20269804.svg)](https://doi.org/10.5281/zenodo.20269804)

Python pipeline for quantifying mitochondrial distribution in fluorescence microscopy images of micropatterned cells. The pipeline localizes fibronectin micropatterns via template matching, segments nuclei, and computes per-cell radial slab metrics of mitochondrial intensity used in the accompanying paper.

Sister repos under [`gladkovalab`](https://github.com/gladkovalab):

- [`synthetic-cargo-accumulation-pipeline`](https://github.com/gladkovalab/synthetic-cargo-accumulation-pipeline) — image-analysis pipeline quantifying Miro1 synthetic-cargo distribution in fixed microscopy (nuclear segmentation, perinuclear Gini, edge-spot detection)
- [`synthetic-cargo-particle-tracking`](https://github.com/gladkovalab/synthetic-cargo-particle-tracking) — single-particle tracking of the same Miro1 cargo across TRAK isoform conditions (TrackMate output)

## Overview

Cells are cultured on fibronectin micropatterns that constrain cell shape into a reproducible geometry. For each ND2 image, the pipeline:

1. Localizes the micropattern centre via FFT-based template matching against an SVG of the pattern.
2. Extracts a 1024×1024 px crop centred on the pattern.
3. Segments the nucleus (405 nm) with Otsu thresholding.
4. Computes Euclidean distance transforms from the nuclear boundary and the pattern arch.
5. Subtracts background from the mitochondrial channel (488 nm).
6. For each cell, computes the wedge-r profile (% mitochondrial intensity in 1 µm radial bins from the pattern arch) and two radial slab metrics:
   - **Centrosomal slab**: % wedge intensity in [18, 33) µm
   - **Peripheral slab**: % wedge intensity in [41, 56) µm
7. Writes per-cell metrics to CSV/Excel plus diagnostic PDFs.

The per-plate CSVs are then aggregated into the slab-metrics Excel workbook used as the paper's data source.

## Installation

Dependencies are managed with [Pixi](https://pixi.sh). With Pixi installed:

```bash
pixi install
```

This installs Python 3.12 and all required packages into an isolated environment on Linux or macOS.

### Data-root setup

The driver and pipeline expect raw ND2 files under `mark_data/patterned_data/{plate}/{well}_{condition}/*.nd2`. `mark_data/` is gitignored — make it a symlink (or directory) pointing at wherever the ND2 data actually lives:

```bash
ln -s /path/to/patterned_data mark_data/patterned_data
```

Or override per invocation with `--data-root /path/to/patterned_data` on `analysis/run_pipeline_paths.py`.

## Usage

### Paper-regen workflow (canonical)

End-to-end, from an empty checkout to the paper data tables and figures, on a machine with the ND2 data mounted under `mark_data/patterned_data/`:

```bash
# 1. Run the patched pipeline for every comparison sheet. Writes per-well CSVs
#    under analysis/wedge_r_ks_out_all_denoised/by_well/ and per-cell NetCDF
#    projections under projections/. Allow several hours per sheet.
nohup caffeinate -dimsu bash -c '
  for sheet in "TRAK isoform (peroxisome)" "TRAK isoform (60mer)" \
               "TRAK1 helix muts" "TRAK2 helix muts" \
               "MAPK9 siRNA" "TRAK isoform (mito)"; do
    pixi run python analysis/run_pipeline_paths.py \
      --sheet "$sheet" --variant denoised \
      --out-root analysis/wedge_r_ks_out_all_denoised
  done
' > analysis/wedge_r_ks_out_all_denoised/stdout_overnight.log 2>&1 &

# 2. Paper figures (per-sheet 6-panel with profile + CDF + scalar strips).
for label in "trak_isoform_mito:TRAK isoform (mito)" \
             "trak_isoform_peroxisome:TRAK isoform (peroxisome)" \
             "trak_isoform_60mer:TRAK isoform (60mer)" \
             "trak1_helix_muts:TRAK1 helix muts" \
             "trak2_helix_muts:TRAK2 helix muts" \
             "mapk9_sirna:MAPK9 siRNA"; do
  slug=${label%%:*}; sheet=${label#*:}
  pixi run python analysis/plot_metrics.py --sheet "$sheet" \
    --out analysis/figures_wedge_r_ks/${slug}.png
done

# 3. Source data exports (XLSX/CSV).
pixi run python analysis/export_slab_metrics_by_plate_xlsx.py
pixi run python analysis/export_wedge_profiles_csv.py
pixi run python analysis/export_wedge_profiles_xlsx.py
pixi run python analysis/export_wedge_profiles_by_plate_xlsx.py

# 4. Supporting figures (profile-with-bands, nuclear overlays, wedge illustration).
pixi run python analysis/plot_profiles_with_bands.py
pixi run python analysis/plot_all_with_nuclear.py
pixi run python analysis/plot_60mer_with_nuclear.py
pixi run python analysis/plot_wedge_illustration_offline.py

# 5. Per-condition mean-projection TIFFs.
# Writes to comparison_projections/{sheet}/{condition}/.  The paper-relevant
# output is mean_488_bg_subtracted.tif (mean over the bg-subtracted 488
# MaxIPs — same MaxIPs the wedge-r / slab metrics quantify).  Also emits
# mean_405.tif from the z-sum (nuclear overlay; not bg-subtracted) and a
# pair of legacy z-sum-based 488 TIFFs which are *not* what the paper uses.
pixi run python analysis/generate_comparison_projections.py

# 6. (Optional) Audit: per-sheet nested-ANOVA + Šídák stats as text.
pixi run python analysis/stats_summary.py
```

Each script has `--help` and writes its output under `analysis/figures_wedge_r_ks/` by default. See `analysis/HANDOFF_v4.md` §6 for the methodology trail and the rationale for the loop above.

### Single-plate diagnostic (debugging only — not the paper-regen path)

Process all ND2 files under a single plate directory, bypassing the driver:

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
- `template_matching/{plate}/{well}/template_matching.csv` — per-cell metrics (one row per cell, includes `wedge_r_*_pct` columns)
- `template_matching/{plate}/{well}/template_matching.xlsx` — same data as Excel
- `template_matching/{plate}/{well}/*.pdf` — diagnostic figures per cell

#### Cluster submission

For large datasets, `bsub_analysis.sh` submits one LSF job per directory using the paths listed in `config/20251229_paths_for_analysis.txt`:

```bash
bash bsub_analysis.sh
```

To target a different set of directories, edit `config/20251229_paths_for_analysis.txt` (one path per line) or `bsub_analysis.sh` directly. Each job requests 8 cores via `bsub -n 8 -P vale`.

## Script Reference

| Script | Description |
|--------|-------------|
| `template_matching_bulk.py` | **Main pipeline.** Batch processes ND2 files: template matching, cropping, nuclear segmentation, background subtraction, wedge-r profile and slab metric computation, output to NetCDF4/CSV/Excel/PDF. |
| `bsub_analysis.sh` | LSF cluster submission wrapper around `template_matching_bulk.py`. |
| `analysis/run_pipeline_paths.py` | Driver — walks the comparisons table and invokes the pipeline per cell, writing the combined per-well CSV consumed by everything else under `analysis/`. |
| `analysis/plot_metrics.py` | Shared utilities (`SHEET_CONFIG`, `load_template_matching`, `join_with_metadata`) plus the per-sheet 6-panel figure generator (`make_figure`). |
| `analysis/export_slab_metrics_by_plate_xlsx.py` | Per-comparison Excel of the two slab metrics, one column per (condition, plate-date). |
| `analysis/export_wedge_profiles_csv.py` | Long-format CSV of per-bin profile means / SEMs across (sheet, condition, radial bin). |
| `analysis/export_wedge_profiles_xlsx.py` | XLSX of the same source data, one worksheet per (sheet, condition). |
| `analysis/export_wedge_profiles_by_plate_xlsx.py` | XLSX broken out per plate-date, with `_mean` and `_sem` columns for plotting bands. |
| `analysis/plot_60mer_with_nuclear.py` | Three-panel figure (no TRAK / TRAK1 / TRAK2) of the 60mer wedge-r profile with the per-condition nuclear (405) profile overlaid as a dashed line. |
| `analysis/plot_all_with_nuclear.py` | Per-sheet split views: one panel per condition, with 488 mitochondria profile, 405 nuclear mask radial distribution, and the two slab bands shaded. |
| `analysis/plot_profiles_with_bands.py` | 2×3 grid of wedge-r profiles per sheet with the two slabs as grey shadings. |
| `analysis/plot_wedge_illustration_offline.py` | Canonical wedge-on-cell illustration; reads cached projections, no ND2 needed. |
| `analysis/generate_comparison_projections.py` | Per-condition mean projection TIFFs. Walks the comparisons table and accumulates four variants per (sheet, condition): z-sum, denoised z-sum, bg-subtracted MaxIP, and denoised bg-subtracted MaxIP. The paper figures use `mean_488_bg_subtracted.tif` (the bg-subtracted MaxIP mean — same source as the wedge-r / slab quantification); the z-sum 488 variants are kept for backward compatibility but should not be used in figures. |
| `analysis/comparison_loader.py` | Utility module backing the projection script: comparisons-table loading, well-dir resolution, per-cell NetCDF loading, per-channel center-padded summing, uint16 conversion. |
| `analysis/stats_summary.py` | Audit dump of nested-ANOVA + Šídák pairwise stats per sheet (the numbers shown in figure brackets). |
| `analysis/HANDOFF_v4.md` | Methodology trail: how the slab edges were chosen, the isobestic-point derivation, validation against the prior perinuclear metric. |
| `analysis/WEDGE_R_KS.md` | Wedge-r KS metric methodology and geometry. |

## Configuration

### Coordinate overrides

When template matching fails or requires manual correction, add entries to `coordinate_overrides.csv` (no header):

```
/path/to/Cell1.nd2,x_pixels,y_pixels
```

The `x` and `y` values specify the top of the pattern (not the centre); the pipeline converts these to pattern-centre coordinates automatically.

### Offset overrides

For images where the pattern is not near the default image offset, add entries to the `offset_overrides` dict in `template_matching_bulk.py` (search the source for `offset_overrides = {`). The default offset is `[128, 128]` pixels.

### ROI overrides

To restrict template matching to a sub-region of an image, add entries to the `roi_overrides` dict in `template_matching_bulk.py` (search the source for `roi_overrides = {`).

### Slab band definitions

The two slab bands are defined as module constants in `template_matching_bulk.py`:

```python
WEDGE_CENTROSOMAL_BINS = (18, 33)   # [lo, hi)  — % wedge intensity in this slab
WEDGE_PERIPHERAL_BINS  = (41, 56)   # [lo, hi)  — % wedge intensity in this slab
```

The slab values are computed per-cell from the `wedge_r_NN_NN+1um_pct` columns in the pipeline output.

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
    ├── Cell1.nc                       # all-channel cropped projection
    ├── Cell1_488_bg_subtracted.nc     # background-subtracted mito channel
    └── denoised/                      # denoised variants (if acquired)

template_matching/
└── {plate_name}/{well_id}_{condition}/
    ├── Cell1.pdf                      # diagnostic figures
    ├── template_matching.csv
    └── template_matching.xlsx
```

## Output Metrics

Each row in `template_matching.csv` corresponds to one cell. The columns most relevant to the paper:

| Column | Description |
|--------|-------------|
| `path` | Path to the source ND2 file |
| `template_matching_score` | Fractional overlap of Otsu-thresholded image with the pattern template (quality control) |
| `wedge_r_NN_NN+1um_pct` | % mitochondrial wedge intensity in radial bin `[NN, NN+1)` µm (one column per 1 µm bin, e.g. `wedge_r_00_01um_pct` … `wedge_r_55_56um_pct`) |
| `wedge_r_centrosomal_18_33um_pct` | Centrosomal slab metric: sum of wedge bins in [18, 33) µm |
| `wedge_r_peripheral_41_56um_pct` | Peripheral slab metric: sum of wedge bins in [41, 56) µm |

The legacy `peripheral_{d}um_*` and `perinuclear_{d}um_*` columns (d = 1..5 µm) are also produced by the pipeline but are not used in the paper's headline analysis.

## Methods

Detailed methods suitable for a Methods section are provided in:

- [`methods/2026_05_12_methods.md`](methods/2026_05_12_methods.md) — template matching, segmentation, quantification, and mean-projection imaging
- [`methods/2026_03_06_references.md`](methods/2026_03_06_references.md) — package citations

## Dependencies

Key packages (see `pixi.toml` for the full list):

| Package | Purpose |
|---------|---------|
| `nd2` | Reading Nikon ND2 microscopy files |
| `numpy` / `scipy` | FFT-based template matching, distance transforms |
| `scikit-image` | Otsu thresholding, contour finding, region properties |
| `xarray` / `netCDF4` | Labelled image arrays, projection file I/O |
| `polars` / `fastexcel` / `xlsxwriter` | DataFrame operations and Excel output |
| `matplotlib` | Diagnostic figure generation |
| `cairosvg` / `pymupdf` | Rasterizing the SVG micropattern template |

## Authors

- Mark Kittisopikul — [kittisopikulm@janelia.hhmi.org](mailto:kittisopikulm@janelia.hhmi.org) — pipeline
- William Grant — [hello@wpg.io](mailto:hello@wpg.io) — slab metrics and paper-output curation

## License

Copyright © 2025 Howard Hughes Medical Institute

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

- Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
- Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
- Neither the name of HHMI nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

See [`LICENSE.md`](LICENSE.md) for the full license text.
