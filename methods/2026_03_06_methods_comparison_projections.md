# Methods – Comparison Projection Generation
*2026-03-06*

**Experimental Condition Layout**
- Experimental conditions and well assignments across plates were specified in a structured Excel spreadsheet (`Comparisons_table_v3.xlsx`), read using fastexcel (v0.18) and Polars (v1.36.1).

**Loading Aligned Projections**
- Per-cell cropped z-sum projection images (produced by the template matching pipeline) were loaded from NetCDF4 files using xarray (v2025.7.1). Both standard and denoised image variants were handled separately.

**Mean Projection Computation**
- For each experimental condition, per-cell 2D projection arrays from the 488 nm (mitochondria) and 405 nm (nucleus) channels were spatially aligned by center-padding to a common bounding box and summed, then divided by cell count to produce a mean projection image per condition. This was performed both across all wells and per-well.

**Background-Subtracted Projections**
- Background-subtracted 488 nm projections (saved separately during template matching) were aggregated identically to produce mean background-subtracted images per condition.

**Output**
- Mean projections were saved as 16-bit TIFF files using scikit-image (v0.25.2), organized by condition and further subdivided by individual well, enabling both condition-level and well-level comparisons.
