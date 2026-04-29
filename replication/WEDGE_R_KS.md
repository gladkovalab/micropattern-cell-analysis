# Wedge-r KS metric

Reviewer-facing description of the metric added to
`template_matching_bulk.score_template_match` on this branch.

## What it measures

For each cell, after MaxIP across Z and Mark's existing background
subtraction on the 488 (mitochondria) channel, we compute an
intensity-weighted CDF along a fixed radial axis and report two
Kolmogorov–Smirnov scalar distances:

  * `wedge_r_ks_vs_uniform` — KS distance between the per-cell CDF and an
    analytical area-uniform sector reference. Larger values mean the
    intensity is more concentrated than a uniform fill of the same wedge
    (perinuclear pile-up, peripheral pile-up, or any other deviation
    from `r²/R²` accumulation).
  * `wedge_r_ks_vs_60merNoTRAK` — KS distance between the per-cell CDF and
    an empirical reference CDF averaged across 13 60mer no-TRAK control
    cells. The 60mer no-TRAK condition is the cleanest available proxy
    for passive cytoplasmic fill on this micropattern, so this is the
    biological null.

Both scalars live in `[0, 1]` and are unsigned. They are reported per
cell alongside `perinuclear_5um_percent_total` and
`peripheral_5um_percent_total` (now MaxIP-based, see "Projection change"
below).

## Wedge geometry

The wedge is a polygonal cone fixed by the rigid micropattern and
identical for every cell. In the cropped-image frame (1024×1024, with
the pattern centered at `(512, 512)`):

| Anchor | (y, x) | role                          |
|---     |---     |---                            |
| apex   | (896, 512) | pattern bottom extreme        |
| L      | (373, 281) | upper-left tangent            |
| R      | (374, 742) | upper-right tangent           |

The wedge sweeps the upper hemicircle (the pattern's arch) from the
apex through the L and R tangent points, with an opening angle of
~47.6°. Per-pixel `r` is the µm distance from the apex; per-bin volume
`vol_arc` is the geometric pixel count per 1 µm shell within the wedge,
used as the area-uniform reference.

The geometry is hardcoded as module-level constants
(`WEDGE_APEX`, `WEDGE_LEFT`, `WEDGE_RIGHT` in `template_matching_bulk.py`)
and built once per `(shape, pixel_pitch)` combination via
`_get_wedge_geometry`. The committed `replication/wedge_illustration.png`
shows the canonical wedge overlaid on a real cell (plate 3 / B04 /
cell1, TRAK2 condition); regenerate it with
`pixi run python replication/plot_wedge_illustration.py`.

## Projection change

Mark's original pipeline used a Z-sum projection for the 488 channel.
The denoising sensitivity check across 10 cells (zsum_wedge_r_gini and
zsum_wedge_r_ks_vs_uniform shifted ~9% per cell on raw → denoised
input; the same scalars on MaxIP shifted < 1.1%). On that basis the
patched pipeline switches the 488-channel projection to MaxIP for all
downstream metrics. The 405 nuclear-segmentation channel keeps Z-sum so
the existing Otsu segmentation is undisturbed.

This means `perinuclear_5um_percent_total` and
`peripheral_5um_percent_total` now report MaxIP-based values; they will
not match Mark's previously published Z-sum numbers. Within-condition
ranking is preserved and effect sizes are typically larger.

## Empirical reference CDF (60mer no-TRAK)

`_REF_CDF_60MER_NOTRAK` is a 60-element constant in
`template_matching_bulk.py`, baked at the values produced by the v3
whole-dataset run (`replication/overnight_final_out/combined_raw.csv`,
filtered to `sheet == "TRAK isoform (60mer)" && condition == "no TRAK"`,
n = 13 cells, per-cell wedge-r CDFs averaged element-wise).

If the reference ever needs to be regenerated, the script
`replication/ks_vs_60mer_reference.py` on the `wpg/alt-metrics` branch
shows how. The constant should be regenerated whenever the wedge
geometry, projection, or background subtraction changes.

## Validation

A single TRAK1 mito cell (plate 3 / B03 / cell1) processed by this
patched pipeline reports `wedge_r_ks_vs_uniform = 0.3968`; the v3
reference value for the same cell is `0.3950`. The 0.5% offset is
consistent with the fact that this branch uses Mark's stretch-then-
percentile background subtraction whereas v3 used a raw-percentile
variant. KS rankings within a sheet are preserved.

## Output columns added

Per cell, the patched `score_template_match` now emits these new
columns in addition to Mark's existing output:

  * `wedge_r_ks_vs_uniform`            — scalar in `[0, 1]`
  * `wedge_r_ks_vs_60merNoTRAK`        — scalar in `[0, 1]`
  * `wedge_r_NN_NN+1um_pct` for `NN ∈ {00, 01, …, 59}`
                                       — per-bin intensity % within the
                                         radial window [0, 60) µm; sums
                                         to 100 (or all-NaN when the
                                         wedge has no signal in window)

The 60 per-bin profile columns let downstream scripts (e.g.
`replication/plot_metrics.py`) reconstruct the wedge-r profile and CDF
without re-running the pipeline.
