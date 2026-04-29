# Wedge-r KS reviewer rebuttal — handoff

This document captures everything done in the 2026-04-28/29 session on
top of the wedge-r KS work from `0532e5f`. The driving goal was to
respond to a reviewer who flagged the original unsigned KS metric as
direction-agnostic, and to find a directional companion that survives
the dataset's statistics. The result is two new pipeline-emitted
scalars (centrosomal and peripheral slabs) plus the analysis trail that
justifies their definition.

---

## 1. What changed in the pipeline

### `template_matching_bulk.py`

* Two new module-level constants document the slab definition:
    ```python
    WEDGE_CENTROSOMAL_BINS = (18, 33)   # [lo, hi)
    WEDGE_PERIPHERAL_BINS  = (41, 56)   # [lo, hi)
    ```
* `score_template_match` now emits two extra columns alongside the
  existing wedge-r KS metrics:
    ```
    wedge_r_centrosomal_18_33um_pct  = sum of wedge_r_18_19..32_33 (% of wedge total)
    wedge_r_peripheral_41_56um_pct   = sum of wedge_r_41_42..55_56 (% of wedge total)
    ```
  Both inherit the all-NaN-when-empty convention from the per-bin
  profile, so downstream stats can rely on the same NaN handling as the
  KS scalars.

### `replication/run_pipeline_paths.py`

* Added a `_backfill_slabs` helper that re-derives the slab columns
  from the per-bin `wedge_r_NN_NN+1um_pct` columns when assembling the
  per-well CSV and the combined CSV. This means cached JSON
  checkpoints from earlier runs (which predate the new pipeline output)
  are upgraded transparently — no need to re-process cells.

### `replication/plot_metrics.py`

* Restructured the per-sheet figure into a 3-row grid (`figsize=(15,
  14)`):
    * row 0 — wedge-r profile + wedge-r CDF (unchanged)
    * rows 1+2 — six strip-plot panels in a 3×2 layout: peripheral
      5 µm, perinuclear 5 µm, **centrosomal slab (18–33 µm, %
      wedge)**, **peripheral slab (41–56 µm, % wedge)**, KS vs
      area-uniform, KS vs 60mer no-TRAK
* Same nested-ANOVA + Šídák brackets are applied to the new panels
  automatically — no per-metric special-casing.

---

## 2. How the slab bin edges were chosen

We did this end-to-end in the session — the short version:

1. We started from a single 10 µm slab centred on the distal arch
   point at r ≈ 50 µm (`[45, 55)`). It worked on most sheets but only
   reached *trend* (`p = 0.067`) for TRAK1 vs TRAK2 on the mito sheet.
2. We swept all 10 µm windows (`replication/sweep_window_mito.py`) and
   found two distinct zones of significance: an inner mid-zone window
   around `[19, 29) – [22, 32)` and an outer rim-adjacent window around
   `[35, 45) – [40, 50)`. The two zones flip directions: TRAK2
   accumulates inner, TRAK1 accumulates outer. They are direct
   readouts of the underlying biology and exactly what the unsigned
   KS metric was averaging together.
3. To pick non-arbitrary bin edges, we computed the **TRAK1/TRAK2
   isobestic point** on the 60mer profile (the reference comparison
   sheet for this metric):
    * isobestic crossings at r ≈ 16.0 µm and r ≈ 36.8 µm
    * the 36.8 µm crossing sits in the natural separation between the
      inner and outer significance zones, so we used it as the centre
      of the gap between the two slabs
4. Final iso-centred slab definition:
    * gap: 8 µm wide, centred on r = 36.83 µm → bin edges 33 and 41
    * inner slab `[18, 33)` µm — 15 µm wide, "centrosomal"
    * outer slab `[41, 56)` µm — 15 µm wide, "peripheral"

The figure that visualises this (60mer wedge-r profile with the two
isobestic crossings marked) is regenerable via
`replication/plot_60mer_with_nuclear.py`.

---

## 3. Why the slabs beat the original perinuclear-5µm metric

This was the core of the rebuttal. The original `perinuclear_5um_pct`
metric is a 2D 5 µm shell around the nuclear boundary in image space.
We reconstructed where that shell sits in wedge-r space using the
saved 405 Z-sum projections (`replication/plot_all_with_nuclear.py`):

* the binary nuclear mask peaks at r ≈ 28 µm
* the 5 µm dilation halo peaks at r ≈ 32 µm and spans r ≈ 16–40 µm
* **51%** of the halo signal sits in the centrosomal slab
* **27%** of the halo signal sits in the iso-point gap `[33, 41)`
* only **6%** sits in the peripheral slab

The 27% leakage into the gap is the source of the resolution loss:
that's the radial transition zone where TRAK1 and TRAK2 cross over,
and including it averages the two opposite biological behaviours into
one number. The slab metric stops just before the cross-over, so the
TRAK2 perinuclear signal isn't diluted by rim-side TRAK1 leakage.

Concrete numbers on the 60mer comparison:

| metric                       | TRAK1   | TRAK2   | p (TRAK1 vs TRAK2) |
| ---                          | ---     | ---     | ---                |
| `perinuclear_5um_pct_total`  | 22.1%   | 52.9%   | *** 1.4e-04        |
| **`wedge_r_centrosomal...`** | 14.9%   | 56.2%   | **** 1.4e-05**     |

Same direction, ~10× tighter p-value.

---

## 4. Validation

Spot-check on a representative TRAK2 mito plate-3 cell (used as the
canonical wedge-illustration cell): the new pipeline computes
`centrosomal = 65.6%`, `peripheral = 1.9%` from the per-bin profile —
matches summing `wedge_r_18_19um_pct .. wedge_r_32_33um_pct` from the
existing checkpoint exactly.

Across 50 sampled cells, all centrosomal/peripheral values are in
[0, 100] and the bin slices `WEDGE_CENTROSOMAL_BINS` /
`WEDGE_PERIPHERAL_BINS` cover bins 18..32 and 41..55 respectively
(as designed).

The full mito-sheet pairwise stats reproduce the alt-metrics evaluator
within rounding (e.g. `wedge_r_ks_vs_uniform` TRAK1 vs TRAK2: this
branch p=0.0013, alt-metrics p=0.0013, matches the original commit
message).

---

## 5. Figures and source data shipped to the collaborator

All under `replication/figures_wedge_r_ks/`:

| file                                       | what                                                                                   |
| ---                                        | ---                                                                                    |
| `trak_isoform_mito.png`                    | mito 6-panel figure (now with centrosomal/peripheral slab strip plots populated)       |
| `trak_isoform_peroxisome.png`              | same layout, peroxisome sheet                                                          |
| `trak_isoform_60mer.png`                   | same layout, 60mer sheet                                                               |
| `trak1_helix_muts.png`                     | same layout, TRAK1 helix muts                                                          |
| `trak2_helix_muts.png`                     | same layout, TRAK2 helix muts                                                          |
| `mapk9_sirna.png`                          | same layout, MAPK9 siRNA                                                               |
| `profiles_with_bands_iso.png`              | 6-panel grid of wedge-r profiles per condition with the two slabs as grey slabs       |
| `profiles_*_with_nuclear.png`              | 6 per-sheet split views (one panel per condition) with the binary nuclear mask + 5 µm halo overlaid |
| `wedge_illustration_offline.png`           | canonical wedge-on-cell illustration (template outline, apex, L/R, wedge mask)        |
| `band_18_33um_all.png` / `band_41_56um_all.png` | 6-sheet grids of just the slab strip plots, large-font labels for emailing             |
| `band_40_50um_all.png` / `band_40_55um_all.png` | window-position sensitivity figures                                                   |
| `sweep_window_mito.png`                    | -log10(p) vs window position on the mito sheet                                        |
| `wedge_r_profiles_source.csv` / `.xlsx`    | long-format source data for the 1D plots (per-bin mean/SEM per condition per sheet)   |
| `wedge_r_profiles_by_plate.xlsx`           | same data per (condition, plate_date), one worksheet per comparison                   |
| `slab_metrics_by_plate.xlsx`               | per-cell raw values for the two slab metrics, broken down by (condition, plate)       |

---

## 6. Reproducing the figures from a fresh checkout

The pipeline output `replication/wedge_r_ks_out_all_denoised/` is
gitignored. To rebuild it:

```sh
# Run the patched pipeline against the SMB-mounted data root.
nohup caffeinate -dimsu bash -c '
  for sheet in "TRAK isoform (peroxisome)" "TRAK isoform (60mer)" \
               "TRAK1 helix muts" "TRAK2 helix muts" \
               "MAPK9 siRNA" "TRAK isoform (mito)"; do
    pixi run python replication/run_pipeline_paths.py \
      --sheet "$sheet" --variant denoised \
      --out-root replication/wedge_r_ks_out_all_denoised
  done
' > replication/wedge_r_ks_out_all_denoised/stdout_overnight.log 2>&1 &
```

Then the figures:

```sh
for label in "trak_isoform_mito:TRAK isoform (mito)" \
             "trak_isoform_peroxisome:TRAK isoform (peroxisome)" \
             "trak_isoform_60mer:TRAK isoform (60mer)" \
             "trak1_helix_muts:TRAK1 helix muts" \
             "trak2_helix_muts:TRAK2 helix muts" \
             "mapk9_sirna:MAPK9 siRNA"; do
  slug=${label%%:*}; sheet=${label#*:}
  pixi run python replication/plot_metrics.py \
    --template-matching replication/wedge_r_ks_out_all_denoised/by_well \
    --sheet "$sheet" \
    --out replication/figures_wedge_r_ks/${slug}.png
done
```

Each post-hoc analysis script under `replication/` is independent and
reads only from `replication/wedge_r_ks_out_all_denoised/by_well/**`
plus `projections/**` for the nuclear-mask figures. None of them
re-runs cells or modifies the pipeline.

---

## 7. Post-hoc scripts inventory

All under `replication/`. None of these are used by the batch
pipeline — they're exploration/figure-building tools:

| script                                    | role                                                                      |
| ---                                       | ---                                                                       |
| `signed_ks_analysis.py`                   | signed KS sign distribution per condition (showed sign is always + vs uniform, mostly + vs 60mer) |
| `stats_summary.py`                        | text dump of nested-ANOVA + Šídák stats per sheet                         |
| `plot_bin_20_30um.py`                     | 6-sheet grid for an arbitrary inner radial slab                           |
| `plot_arch_band_polar.py`                 | parametric `(BIN_LO, BIN_HI)` 6-sheet grid; used for [40,50), [40,55), [45,55), [18,33), [41,56) windows |
| `sweep_window_mito.py`                    | -log10(p) vs 10 µm window position on the mito sheet                      |
| `plot_profiles_with_bands.py`             | wedge-r profile per condition with both slab grey shadings                |
| `plot_60mer_with_nuclear.py`              | single-sheet 60mer split + nuclear mask overlay (early version; superseded by `plot_all_with_nuclear.py`) |
| `plot_all_with_nuclear.py`                | per-sheet split views with binary nuclear mask + 5 µm halo overlays       |
| `plot_wedge_illustration_offline.py`      | offline twin of the original `plot_wedge_illustration.py`; reads saved projections, no SMB needed |
| `export_wedge_profiles_csv.py`            | long-format CSV of per-bin profile means/SEMs                             |
| `export_wedge_profiles_xlsx.py`           | XLSX, one worksheet per (sheet, condition)                                |
| `export_wedge_profiles_by_plate_xlsx.py`  | XLSX per comparison, columns broken out per (condition, plate_date)       |
| `export_slab_metrics_by_plate_xlsx.py`    | XLSX per comparison, two slab blocks side-by-side, per-cell raw values per (condition, plate) |

---

## 8. Open follow-ups

* The 60mer empirical reference CDF baked into
  `template_matching_bulk._REF_CDF_60MER_NOTRAK` is still the v3-era
  constant; if anyone regenerates it from this run's data, update both
  the constant and a short note here.
* The wedge-r profile is *not* volume-normalized in the displayed
  curves (intensity per shell is divided by total wedge intensity, not
  by per-bin pixel volume). `ks_vs_uniform` correctly handles the
  geometric weighting via `vol_arc`. If a reviewer asks for a
  density-form profile, a 3-line tweak in `wedge_r_profile` would do
  it (see session notes).
* The 5 conditions × 6 sheets pipeline took ~24 wall-clock hours over
  SMB. If we re-run for any reason, set the output root to a fresh
  directory so the existing checkpoints aren't overwritten.

---

Last updated: 2026-04-29 — William Grant
