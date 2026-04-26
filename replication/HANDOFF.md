# Handoff — Fig 4 / S11 alt-metric pitch

*Written 2026-04-23. Updated to reflect complete MaxIP pipeline run.*

Start here when you come back. This is a quick reorientation, **not** a
repeat of `RESULTS.md` (which is the actual pitch memo for the coauthors).

---

## Where things stand

**Branch**: `wpg/alt-metrics` pushed to
`git@github.com:gladkovalab/micropattern_cell_analysis.git`. Upstream is
`JaneliaSciComp/micropattern_cell_analysis`.

**Main deliverable**: `replication/RESULTS.md` — pitch memo for Mark and
the senior author. §0 is the paper-claims audit (confirms every "recovered
significance" corresponds to a main-text claim the paper makes). §2 has
per-panel tables with 5 metric options each. §3 is the diff-vs-ratio
trade-off and recommended pitch.

**Pipeline run**: complete. 429 cells × 89 metrics across all 49 wells
(TRAK1 helix, TRAK2 helix, TRAK isoform mito, MAPK9 siRNA). Per-cell
data in `replication/overnight_out/combined.csv`. Scored results in
`replication/overnight_eval_out/summary.csv`. Zero processing errors in
the final run.

**Headline findings** (all in `RESULTS.md` §0):

1. **MaxIP projection + peri-nuc polarization metrics** strictly improve
   over Mark's z-sum perinuclear metric on every Fig S11 pair, and crucially
   recover significance on three comparisons the paper text claims as
   significant but Mark's current metric misses:
   - Fig 4D TRAK2 vs mDRH peripheral: ns 0.39 → **\* 0.01** (MaxIP diff)
   - Fig S11 E TRAK2 mDRH vs mSpindly perinuclear: ns 0.052 → **\*\*\* 0.0003** (MaxIP ratio)
   - Fig S11 F ctrl vs ctrl-Ars perinuclear: ns 0.083 → **\*\*\* 0.0007** (MaxIP ratio)

2. **Cutoff-combination sweep** (5 × 5 grid of zone radii): (5 µm, 5 µm)
   is the best-or-tied-best combination on virtually every pair. Defends
   Mark's methodological choice. No cherry-picking.

3. **Trade-off between diff and ratio**: diff catches TRAK2 wt→mDRH
   (both Fig 4D and S11 E primary claim), ratio is much stronger on MAPK9
   + Fig 4B. Neither strictly dominates. Recommended primary: MaxIP diff
   (preserves the TRAK2 claim which the paper relies on) with an honest
   "ratio is stronger on MAPK9" note.

---

## Open follow-up work

### High priority

**Zone-area normalization** (queued task #19). Mark's metric is
zone_signal / total_cell_signal. A "density" version would be
zone_signal / zone_area, making the polarization metric area-independent
(expected value 0 for diff, 1.0 for ratio under uniform distribution).
Whether this matters empirically depends on whether nucleus/zone size
systematically varies by condition.

To test: patch `metric_pipeline.py` to save
`perinuclear_5um_area_px`, `peripheral_5um_area_px`, and
`nucleus_area_px`, re-run on all 49 wells (~15 hours), rebuild_combined
and evaluate. Expected outcome: second-order effect that doesn't
overturn the headline MaxIP findings, but good defensive paragraph for
the pitch.

### Lower priority

- **Rerun with plate-weighted nested ANOVA** for exact Prism p-value
  parity. Current framework uses cell-count-weighted; 1.1–1.5× offset on
  raw p; no significance class flips. A 20-line change in
  `replicate_stats.py` if reviewers demand exact parity.
- **MaxIP on peroxisome data (Fig S11 C)**. Not processed overnight.
  Paper's explicit null claim — if MaxIP happens to overturn it, that's
  a real finding rather than a pitch.

---

## How the code is organized

All new work is under `replication/`. Only edit to Mark's pipeline:
`template_matching_bulk.py` (added `MICROPATTERN_DATA_ROOT` env var so
the pipeline runs off-cluster; all other logic unchanged).

### Data flow

```
raw ND2s (mark_data/ or local_staged/)
        │
        ├── template_matching_bulk.py  ─→  mark_data/analysis/260224/**/template_matching.csv
        │       (Mark's original pipeline, lightly patched)                │
        │                                                                  ▼
        │                                           replication/derived_metrics.py
        │                                           → per-cell diff/ratio from z-sum CSVs
        │                                           → Mark baseline vs our alternatives
        │                                           → replication/derived_metrics_out/*.csv
        │
        └── replication/metric_pipeline.py (process_cell: z-sum + MaxIP)
                │                                                          │
                ▼                                                          ▼
        replication/overnight_run.py                    replication/plot_all_panels.py
        (per-well driver with checkpointing)            → replication/figures/*_alt_metrics.{png,pdf}
                │                                        (three-row: Mark + z-sum diff + z-sum ratio)
                ▼
        replication/overnight_out/by_well/**/metrics.csv
                │
                ▼
        replication/rebuild_combined.py  (fix cross-sheet tags)
                │
                ▼
        replication/overnight_out/combined.csv
                │
                ├── replication/evaluate_overnight.py
                │   → replication/overnight_eval_out/summary.csv
                │
                ├── replication/cutoff_sweep.py
                │   → replication/cutoff_sweep_out/cutoff_sweep.csv
                │
                └── replication/plot_all_panels_maxip.py
                    → replication/figures/*_alt_metrics_maxip.{png,pdf}
                    (five-row: Mark + z-sum diff/ratio + MaxIP diff/ratio)
```

### Files that matter (mental model)

| File | Purpose |
|---|---|
| `replication/RESULTS.md` | the pitch, hand to coauthors |
| `replication/HANDOFF.md` | this doc |
| `replication/replicate_stats.py` | nested ANOVA + MixedLM helpers; validated vs Prism |
| `replication/derived_metrics.py` | metrics from Mark's CSVs (no re-imaging) |
| `replication/metric_pipeline.py` | metrics from raw ND2s (z-sum + MaxIP) |
| `replication/overnight_run.py` | driver that runs metric_pipeline per well |
| `replication/stage_missing_wells.sh` | parallel cp SMB → local (for when SMB is flaky) |
| `replication/rebuild_combined.py` | fix cross-sheet tags in combined.csv |
| `replication/evaluate_overnight.py` | scoring framework over overnight data |
| `replication/cutoff_sweep.py` | (X, Y) zone-radius grid search |
| `replication/plot_all_panels.py` | three-row plots (Mark + z-sum diff/ratio) |
| `replication/plot_all_panels_maxip.py` | five-row plots adding MaxIP rows |

---

## Reproduce the full story from scratch

If the pipeline needs to be re-run (e.g. for zone-area normalization):

```bash
# 1. Mount valelab (⌘K in Finder, smb://gladkovac@prfs.hhmi.org/valelab)
# 2. Verify
ls /Volumes/valelab/_for_Mark/patterned_data/ | head

# 3. Run the pipeline against SMB (or via local_staged for robustness)
MICROPATTERN_DATA_ROOT="$(pwd)/mark_data/patterned_data" \
  PYTHONUNBUFFERED=1 \
  nohup pixi run python replication/overnight_run.py \
  >> replication/overnight_out/stdout.log 2>&1 &

# caffeinate to prevent sleep
caffeinate -dimsu -w <pid>

# 4. Consolidate and score
pixi run python replication/rebuild_combined.py
pixi run python replication/evaluate_overnight.py
pixi run python replication/cutoff_sweep.py
pixi run python replication/plot_all_panels.py
pixi run python replication/plot_all_panels_maxip.py
```

The `overnight_run.py` is **idempotent per cell** — if the mount drops
mid-cell, only that one cell is lost. Already-processed cells persist as
`by_well/{plate}/{well}/cells/*.csv` and get aggregated into a
per-well CSV when the well finishes (with a `done.marker`).

---

## Context / provenance gotchas

- **Pair families differ per panel**. Fig S11 D/E/4C/4D use m=2 adjacent
  pairs. Fig 4B / S11 C use m=3 all-pairs. Fig 4E peripheral uses m=3
  against a reference (A-B, A-C, A-D in Prism column order). Fig S11 F
  perinuclear uses m=3 of (A-C, B-D, C-D). See per-panel comments in
  `plot_all_panels.py`.
- **MAPK9 Prism column order differs from the Comparisons sheet.** Prism
  A/B/C/D = (ctrl ctrl, MAPK9 ctrl, ctrl Ars, MAPK9 Ars); sheet columns
  are (ctrl ctrl, ctrl Ars, MAPK9 ctrl, MAPK9 Ars). Lost an hour to this.
- **Plate 9 has older CSV schema** (`perinuclear_percent_total` not
  `perinuclear_5um_percent_total`). `derived_metrics.py` canonicalizes
  old → new internally.
- **Cross-sheet well aliasing**: many wells are referenced by multiple
  sheets under different condition labels (e.g. plate 3 B03 is both
  `TRAK isoform (mito) · TRAK1` and `TRAK1 helix muts · T1 wt`). My
  first overnight_run dropped the duplicate tag; `rebuild_combined.py`
  fixes this.
- **Paper Fig S11 F uses m=3** Šídák (published p-values 0.0827, 0.7397,
  0.7738), while the on-disk Prism snapshot uses m=2 (p-values 0.056,
  0.63). The published figure is authoritative — `plot_all_panels*.py`
  use m=3 for this panel.

---

## SMB mount gotchas (only matters if re-running against SMB)

- `/Volumes/valelab/` drops silently every ~15 min of active SMB I/O on
  this laptop. macOS SMB sessions don't auto-remount after disconnection.
- Silent drops kill the python process without a traceback (looks like
  SIGKILL).
- `overnight_run.py` has per-cell checkpointing; safe to Ctrl-C and
  resume. Just re-run and it picks up where it left off (done.marker
  files + per-cell CSVs).
- `stage_missing_wells.sh` provides a workaround for persistent flakiness:
  parallel cp (3 streams) into `replication/local_staged/`, ~32 MB/s
  aggregate. Processing from local_staged is fully network-independent.

---

## Reference paths

- Manuscript: `manuscript/aeh1475_CombinedPDF_v3.pdf` (Fig 4 page 20,
  Fig S11 page 45). Main text claims for each pair in ¶35-37 of
  `aeh1475_ArticleContent_v2.docx`.
- Mark's Prism files: `mark_data/analysis/260224/prism_plots/{...}/`.
- Comparisons table: `config/Comparisons_table_v3.xlsx`.
