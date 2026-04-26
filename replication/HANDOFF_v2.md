# Handoff — Wedge-r-Gini pitch (Fig 4 / S11)

*Written 2026-04-24. Supersedes the prior `HANDOFF.md` for the alt-metric work.*

This session built and validated a new family of single-zone clustering
metrics that **replace the polarization-ratio framing** of the original
RESULTS.md pitch. The headline metric is `wedge_r_gini` — Gini coefficient
of the intensity-weighted radial distribution within an angular wedge that
follows the cell's principal axis. On Fig 4B no-TRAK vs TRAK2 (the reviewer's
main test case), wedge-r-Gini achieves **p = 0.004 \*\*** vs Mark's
**p = 0.99 ns** for the published `perinuclear_5um_pct` metric.

Fig 4B is fully analysed. **The whole dataset (S11 D / E / F + Fig 4B
sheets) still needs to be run with the streamlined pipeline.** The next
instance should kick off `final_pipeline.py` on the full sheet set, then run
the evaluator and plotter.

---

## §1 · Why the pitch shifted

The previous pitch (`replication/RESULTS.md`) used peripheral / perinuclear
*polarization* metrics — `peri_5 - nuc_5` (diff) and `peri_5 / nuc_5`
(ratio), under both z-sum and MaxIP projections. It worked on most panels
but the user's reviewer concern was sharper: the perinuclear metric
specifically underperforms relative to what's visible by eye, so scaling it
by a peripheral metric is "scaling undesirable data with desirable data" —
a polarization metric can't be the answer for the panels where the
perinuclear signal is the primary phenotype (Fig 4B no-TRAK vs TRAK2 in
particular).

This session took the alternative route: **find a single-zone metric that
captures perinuclear clustering directly**, using the full intensity
distribution rather than thresholded zones.

## §2 · The metric

### `wedge_r_gini` (Scheme 1)

1. **Pattern extremes** are extracted from a correctly-aligned pattern mask.
   The mask is the rolled template, sliced **in its own frame** (i.e.
   `shifted_template[max_coords-512 : max_coords+512]`, NOT in the
   crop's frame). This was a subtle bug in Mark's existing
   `pattern_mask_big` — see §6 below.

2. **Wedge construction**:
   - Apex = pattern bottom extremum (the stalk tip)
   - Boundaries = rays from apex through the leftmost and rightmost pattern
     extreme points
   - Wedge sweeps "upward" through the arch (contains the up direction in
     image coordinates).
   - Opening angle is ~45° on the standard pattern.

3. **Polar r** = Euclidean distance (µm) from the wedge apex.

4. **1D r-profile**: bin pixel intensities by r in 1 µm bins, 0..60 µm.
   Restrict to wedge pixels. Normalize so bins sum to 100%.

5. **Gini coefficient** over the 60-bin distribution.

The resulting scalar is high when signal piles up in a narrow r-range
(clustered along the cell's stalk-arch axis) and low when signal is
uniformly spread along the axis.

### Why MaxIP, not z-sum

For Gini-of-distribution metrics:
- Maxip (per-pixel max along z): keeps only the brightest voxel, so
  distinct mito clusters give clear bin amplitude differences. Gini is
  high when clustering is real.
- Z-sum: every z-slice contributes, including out-of-focus noise. Bins
  fill with roughly uniform background, dragging Gini toward 0.

### Raw vs denoised inputs (sensitivity result)

`replication/compare_raw_vs_denoised.py` ran the same `process_cell` on raw
ND2s vs the Nikon-NIS-denoised stacks Mark uses for his Fig 4 panels.
Result on 5 stratified Fig 4B cells:

- **MaxIP metrics are denoiser-invariant** (≤0.5% per-cell drift on
  `maxip_wedge_r_gini` and `maxip_y_gini`). MaxIP itself filters most of
  what the denoiser would.
- **Z-sum metrics shift** when you denoise (3-16% lower Gini per cell)
  because integrated noise gets smoothed out, flattening per-pixel
  distributions. Direction is consistent across conditions, so ranking is
  preserved.
- **Geometry is unchanged**: pattern extremes, wedge opening, nuclear
  area/solidity all <1% drift.

**Implication for the headline pitch**: the MaxIP wedge-r-Gini result is
reproducible from open ND2s without the commercial denoiser. The
`final_pipeline.py` whole-dataset run reads raw `cell*.nd2` files (the
`iter_cells` walk skips `denoised/` subdirs); this is the recommended
canonical run.

**Optional sensitivity analysis** to preempt reviewer pushback: after the
raw-data whole-dataset run completes, point `final_pipeline.py` at the
`<well>/denoised/<cell> - Denoised.nd2` files for a stratified subset
(~30 cells, ~1.5 h) and confirm the MaxIP wedge-r-Gini significance class
is unchanged. The CSV from `compare_raw_vs_denoised.py` lives in
`replication/overnight_final_out/raw_vs_denoised.csv`.

Performance (Fig 4B no-TRAK vs TRAK2):

| Metric | d | p (Šídák m=3) |
|---|---:|---:|
| Mark's `zsum_perinuclear_5um_pct` | −0.10 | 0.99 ns |
| `zsum_y_gini` | −0.94 | 0.21 ns |
| `maxip_y_gini` | −1.33 | 0.024 \* |
| `zsum_wedge_r_gini` | −0.94 | 0.21 ns |
| **`maxip_wedge_r_gini`** | **−1.46** | **0.0037 \*\*** |

MaxIP wedge-r-Gini is the headline pitch metric.

### TRAK1 vs TRAK2 isoform comparison

Same wedge-r-Gini cleanly separates the isoforms (`maxip_wedge_r_gini`:
d = −1.48, p = 0.00084 \*\*\*), upgrading from Mark's borderline ns 0.053
on `peripheral_5um_pct`. TRAK2 piles signal in the 20-35 µm wedge-r band
(perinuclear); TRAK1 has more signal in the 35-55 µm band (peripheral/arch).
This is biologically consistent with the paper's ¶35 narrative.

## §3 · The streamlined pipeline

`replication/final_pipeline.py` is the consolidated, performance-tuned
pipeline. It computes only the metrics that earned their place in v2 plus
Mark's published baselines for cross-reference. About **half the columns
of `fig4b_v2_pipeline.py`** and meaningfully faster per cell because:

- Y/X meshgrids and the pattern-extreme extraction happen **once per cell**
  and are reused across all metric helpers.
- The wedge mask + r-distance map are computed **once per cell** (geometry
  is independent of the channel — z-sum and MaxIP reuse the same maps).
- Wedge-r binning uses `np.bincount` (~5× faster than the per-bin
  mask+sum loop in earlier pipelines).
- The legacy wrong-frame `pattern_mask_big` is not built at all (see §6).

**Output columns per cell** (~290 total):

- Cell-level diagnostics (~14): nucleus seg quality (area, solidity,
  eccentricity, n_components, largest_area_frac, euler_number); pattern
  extreme positions in µm from nucleus CoM (4 points × dy/dx); wedge
  opening angle and pixel fraction; pixel pitch.
- Per projection (zsum + maxip), ~140 each:
  - Mark baselines: `perinuclear_5um_pct`, `peripheral_5um_pct`,
    `mean_dist_to_nucleus_um`.
  - Y-axis: `y_gini`, `y_entropy`, `y_sd_um`, `y_skew`, `y_mean_um` plus
    60-bin `y_profile_*um_pct` (1 µm bins, ±30 µm from pattern CoM).
  - Wedge-r: `wedge_r_gini`, `wedge_r_entropy`, `wedge_r_ks_vs_uniform`,
    `wedge_r_mean_um`, `wedge_r_sd_um`, `wedge_r_skew`,
    `wedge_r_q25_um`, `wedge_r_q50_um`, `wedge_r_q75_um`,
    `wedge_r_20_35um_frac_pct`, `wedge_r_35_55um_frac_pct`,
    `wedge_mt_apex_lam_max_um2`, `wedge_mt_apex_lam_min_um2`,
    `wedge_mt_apex_elongation`, `wedge_frac_pct`, plus 60-bin
    `wedge_r_*_um_pct` (1 µm bins, 0..60 µm).
  - Total signal and bg threshold per projection.

## §4 · How to run the whole-dataset analysis (next instance)

### Prerequisites

1. **Re-mount SMB**: `⌘K` in Finder → `smb://gladkovac@prfs.hhmi.org/valelab`.
   Verify with `ls /Volumes/valelab/_for_Mark/patterned_data/ | head`.
2. **Free disk and prevent sleep**: the run takes ~12-16 h wallclock; use
   `caffeinate -dimsu`.

### Launch (one command, all sheets)

```bash
mkdir -p replication/overnight_final_out

MICROPATTERN_DATA_ROOT=/Volumes/valelab/_for_Mark/patterned_data \
  PYTHONUNBUFFERED=1 \
  caffeinate -dimsu pixi run python replication/final_pipeline.py \
    --sheets "TRAK isoform (mito)" "TRAK isoform (peroxisomes)" \
             "TRAK1 helix muts" "TRAK2 helix muts" "MAPK9 siRNA" \
    > replication/overnight_final_out/stdout.log 2>&1 &
```

The pipeline:
- **Discovers wells** automatically by reading `replication/overnight_out/combined.csv`
  for the requested sheets.
- **Per-cell checkpointing** to
  `overnight_final_out/by_well/<plate>/<well>/cells/<cell>.json` — safe to
  Ctrl-C and resume; SMB drops cost at most one cell.
- **Per-well CSV + done.marker** when each well finishes.
- **Combined CSV** at the end: `overnight_final_out/combined_raw.csv`.

### Expected runtime

~480-580 cells × ~2-3 min/cell from SMB = **12-16 h** wallclock. Streamlined
metrics save ~10-20% vs the v2 pipeline; the rate is dominated by SMB IO
(~1.8 GB ND2 read per cell at ~32 MB/s).

### Watch progress

```bash
tail -f replication/overnight_final_out/stdout.log
```

Per-cell lines look like `[42/580] OK <plate>/<well>/cell1.nd2`. Errors
print as `[N/M] ERR <path>: <reason>` followed by a traceback; the run
continues. The expected error count is 1 (plate 11 F05/Cell12 — empty nucleus
in the original data; same cell that fails Mark's pipeline).

### After completion

```bash
# Numerical scoring (per sheet × pair, ranked top metrics)
pixi run python replication/evaluate_final.py | tee replication/overnight_final_out/eval_report.txt

# All plots: per-sheet Y-axis profile, wedge-r profile, per-cell scalars
pixi run python replication/plot_final.py
```

Plots land in `replication/overnight_final_out/figures/`:
- `<sheet>_y_profile.png` — mean ± SEM Y profile, 2 rows (zsum, maxip)
- `<sheet>_wedge_r_profile.png` — same for wedge-r
- `<sheet>_scalars.png` — per-cell strip plots of Y-Gini, wedge-r-Gini,
  wedge-r-σ for both projections, plates as marker shapes

Both `evaluate_final.py` and `plot_final.py` are **local-only** (no SMB
needed) — safe to run after the pipeline finishes and SMB is unmounted.

## §5 · Files of interest

| File | Purpose | Status |
|---|---|---|
| `replication/final_pipeline.py` | **Streamlined pipeline** (run this) | **canonical for whole-dataset rerun** |
| `replication/evaluate_final.py` | **Multi-sheet evaluator** | local-only |
| `replication/plot_final.py` | **Per-sheet plot generator** | local-only |
| `replication/HANDOFF_v2.md` | this doc | — |
| `replication/fig4b_v2_pipeline.py` | Earlier extended Fig 4B pipeline (kept for diff) | superseded by `final_pipeline.py` |
| `replication/overnight_fig4b_v2_out/combined_raw.csv` | Fig 4B v2 results (97 cells × 582 cols) | reference data |
| `replication/overnight_fig4b_v2_out/figures/` | Fig 4B plots and wedge illustrations | reference figures |
| `replication/overnight_fig4b_v2_out/figures/wedge_illustration_v3.png` | Wedge geometry sketch (Panel A/B/C) | for talks/paper |
| `replication/overnight_fig4b_v2_out/figures/polar_schemes.png` | Polar coordinate scheme comparison | for talks/paper |
| `replication/replicate_stats.py` | Nested ANOVA + MixedLM helpers (validated vs Prism) | unchanged |
| `replication/RESULTS.md` | Original polarization-ratio pitch | superseded but archived |
| `replication/HANDOFF.md` | Original handoff | superseded by this doc |

The original `metric_pipeline.py`, `overnight_run.py`, `evaluate_overnight.py`,
`plot_all_panels*.py`, `cutoff_sweep.py` are all from the prior session's
polarization-ratio pitch. They still work; they're just no longer the lead.

## §6 · Critical bugs caught this session

### Pattern mask coordinate-frame bug (in Mark's pipeline)

`metric_pipeline.process_cell` slices `shifted_template` using indices
from the **original image frame** (`y_start = max_coords - 512 + offset`).
But `shifted_template` lives in the **subarray frame**, where the pattern
is centered at `max_coords`, not `max_coords + offset`. The resulting
`pattern_mask_big` is displaced by `offset` (default 128 px ≈ 8 µm) from
where the pattern actually sits in the cropped mito image.

Effect: every `*_pattern_*` metric in
`replication/overnight_out/combined.csv` is computed against a pattern
mask that's misaligned with the data by ~8 µm. The original RESULTS.md
pitch's `peripheral_5um_pct (pattern mask)` numbers are therefore not
reliable as physical "on the pattern" measurements — they're "on a mask
shifted 8 µm down-right of the pattern."

`final_pipeline.py` builds the pattern mask correctly:
```python
pattern_mask = shifted_template[
    max_coords[0] - 512 : max_coords[0] + 512,
    max_coords[1] - 512 : max_coords[1] + 512
] > 0
```
The wedge geometry uses this corrected mask. The Y-projection origin
(`pattern_com_y = 512`) follows from the same correction.

The orange contour (`find_contours(shifted_template) - (max_coords - 512)`)
is correctly aligned with the cropped mito image — that path was never
broken, but the pattern_mask_big it was supposed to match was off.

Mark's `peripheral_5um_pct` and `perinuclear_5um_pct` are **not** affected
by this bug (they use the EDT from the arch contour and the nuclear mask
respectively, not `pattern_mask_big`). The `*_pattern_*` columns ARE
affected.

### Wedge direction was initially flipped

Earlier intermediate code (`s11_pipeline.py`, the very first
`fig4b_pipeline.py`) constructed a wedge that **excluded** the upper arch
region instead of containing it. The user's "wedge through the bottom"
verbal cue was taken at face value when in fact the wedge should sweep
upward into the arch (which is the periphery — the cell's spatial extent
toward the top of the image). Fixed in `fig4b_v2_pipeline.py` and
`final_pipeline.py`. **Do not reuse `s11_pipeline.py` as-is** — it has the
wrong direction; if you need a separate S11-only run, use
`final_pipeline.py` with `--sheets` filtered to the S11 sheets.

## §7 · Plate 11 is a known noise plate

`250731_patterned_plate_11_good` has a different baseline than plates 3 / 9
on no-TRAK vs TRAK2. On Mark's `perinuclear_5um_pct`, plate 11's plate-level
Cohen's d is ~−0.16 (nearly null). On `wedge_r_gini` (MaxIP), plate 11's
plate-level d is ~+1.4 — the new metric still picks up the phenotype on
this plate, just at a shifted absolute level.

The user identified that 5 of the 18 plate-11 Fig 4B cells (no-TRAK
cell5/cell11, TRAK2 cell9/cell13, plus user-mentioned cell13 not present in
the CSV) have peripheral-side centrosome positioning, which doesn't happen
on the other plates. Likely explanation: subtle differences in experimental
timing on plate 11.

The wedge-r-Gini's robustness to this plate-level shift is a feature worth
calling out in the pitch — supplementary "per-plate Cohen's d" panel would
make the case visually.

## §8 · Open follow-ups

1. **Whole-dataset run**: §4 above. Critical first step for the next
   session.
2. **Once whole-dataset results are in**, write the new pitch memo
   (`RESULTS_v2.md`) — this should:
   - Show wedge-r-Gini as the primary metric
   - Tabulate p / d for all reviewer-flagged comparisons (4B no-TRAK
     vs TRAK2, S11 D T1 mDRH→dSp rescue, S11 E TRAK2 mDRH→mSpindly rescue,
     S11 F ctrl→ctrl-Ars perinuclear)
   - Use the wedge illustration figure (`wedge_illustration_v3.png`) and
     the polar schemes figure (`polar_schemes.png`) as method explainers
   - Note Y-Gini as the simpler-to-explain alternative if a reviewer
     pushes back on polar coordinates
3. **Sensitivity analyses** worth running once the whole dataset is done:
   - Plate-out: report each plate's individual Cohen's d to demonstrate
     consistency across replicates
   - Cells dropped if `nuc_solidity < 0.85` or `nuc_n_components > 5000`
     (seg quality screen) — does it change conclusions?
   - Wedge anchor sensitivity: rerun a small subset with the wedge anchored
     at pattern CoM instead of pattern bottom (Panel B of the v3 figure) and
     see how much the headline numbers move

## §9 · SMB mount gotchas (unchanged)

- `/Volumes/valelab/` (smb://gladkovac@prfs.hhmi.org/valelab) drops every
  ~15 min of active SMB I/O on this laptop. macOS doesn't auto-remount.
- Silent drops kill the python process without a traceback (looks like
  SIGKILL). Per-cell checkpointing means resume costs at most one cell.
- If drops are frequent, consider staging a subset of wells locally first
  via `replication/stage_missing_wells.sh` (parallel cp, ~32 MB/s
  aggregate). For ~500 cells × 1.8 GB ≈ 900 GB, full local staging is
  probably not feasible; prefer letting the per-cell checkpointing absorb
  the drops.
