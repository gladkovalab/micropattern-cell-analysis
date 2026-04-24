# Handoff — Fig 4 / S11 alt-metric pitch

*Written 2026-04-23 after pushing `wpg/alt-metrics` to the `gladkovalab` fork.*

Start here when you come back. This is a quick reorientation, **not** a
repeat of `RESULTS.md` (which is the actual pitch memo for the coauthors).

---

## Where we left off

**Branch**: `wpg/alt-metrics` pushed to
`git@github.com:gladkovalab/micropattern_cell_analysis.git`. Upstream
remote is `JaneliaSciComp`.

**Main deliverable**: `replication/RESULTS.md` — the pitch memo for Mark
and the senior author. Contents in order:
- §0 Paper-claims audit (confirms every "recovered significance" maps to
  a main-text claim the paper makes that Mark's metric misses)
- §1 Candidate metrics defined
- §2 Per-panel Šídák-corrected p-values (Mark vs diff vs ratio, on all
  seven Fig 4 / Fig S11 panels)
- §3 Diff vs ratio trade-off
- §4 Methodology footnote on classical ANOVA vs MixedLM

**Headline result**: For **TRAK1 wt → mDRH** (the panel closest to
reviewers' complaint about weak significance), **MaxIP-based polarization
metrics push Cohen's d from Mark's −0.75 to +1.51** (p < 0.001), solely
by changing the projection method (per-pixel z-max instead of z-sum).
This beats even Mark's denoised-z-sum peripheral metric — without the
Nikon commercial denoiser dependency.

---

## The unfinished bit — TRAK2 MaxIP data

**Status**: 6 of 12 TRAK2 helix-muts wells are missing MaxIP data because
the SMB mount to `/Volumes/valelab/` drops silently every ~15 min of
active I/O on this laptop. The `RESULTS.md` TRAK2 numbers are currently
from Mark's z-sum CSVs only; the MaxIP numbers there are partial and
unreliable for TRAK2.

**Missing wells** (5 plates × 2 conditions = 6 wells, ~106 GB raw ND2):

| Plate | Well | Condition |
|---|---|---|
| 250612_patterned_plate_3 | B07 | TRAK2 mDRH |
| 250612_patterned_plate_3 | B08 | TRAK2 mDRH / mSpindly |
| 250710_patterned_plate_9_good | C07 | TRAK2 mDRH |
| 250710_patterned_plate_9_good | C08 | TRAK2 mDRH / mSpindly |
| 250731_patterned_plate_11_good | D07 | TRAK2 mDRH |
| 250731_patterned_plate_11_good | E07 | TRAK2 mDRH / mSpindly |
| 250606_patterned_plate_2 | D05 | TRAK2 mDRH / mSpindly |

(That's 7 wells — one of them on plate 2 wasn't even processed over SMB.)

**To finish it**:
```
# 1. Mount valelab (⌘K in Finder, smb://gladkovac@prfs.hhmi.org/valelab)
# 2. Stage all 7 wells to local disk (~1 hour at ~32 MB/s with 3-way parallelism)
bash replication/stage_missing_wells.sh
tail -f replication/local_staged/stage.log   # watch progress

# 3. Run the pipeline against local data (~1 hour, no network)
MICROPATTERN_DATA_ROOT="$(pwd)/replication/local_staged/patterned_data" \
  PYTHONUNBUFFERED=1 \
  nohup pixi run python replication/overnight_run.py \
  >> replication/overnight_out/stdout.log 2>&1 &

# 4. Rebuild the combined CSV (fixes cross-sheet tagging) and re-score
pixi run python replication/rebuild_combined.py
pixi run python replication/evaluate_overnight.py

# 5. Update plots + memo with the completed TRAK2 MaxIP numbers
pixi run python replication/plot_all_panels.py
# manually refresh the TRAK2 rows in RESULTS.md

# 6. Delete the staged copy when done (saves ~110 GB)
rm -rf replication/local_staged/
```

The `overnight_run.py` is **idempotent per cell** — if the mount drops
mid-cell, only that one cell is lost; all previously-completed cells
persist as `replication/overnight_out/by_well/{plate}/{well}/cells/*.csv`
and get aggregated into a per-well CSV when the well finishes.

---

## How the code is organized

All new work is under `replication/`. Modified repo files: only
`template_matching_bulk.py` (added `MICROPATTERN_DATA_ROOT` env var so the
pipeline runs off-cluster; all other logic unchanged).

### Data flow

```
raw ND2s in mark_data/ (or local_staged/)
        │
        ▼
 template_matching_bulk.py  ← Mark's original (lightly patched)
        │                     per-cell perinuclear_5um / peripheral_5um
        ▼
 mark_data/analysis/260224/**/template_matching.csv
        │
        ▼
 replication/derived_metrics.py
        │     per-cell derived metrics (diff, ratio, share, bins, …)
        │     nested ANOVA + MixedLM + Šídák across pair family
        ▼
 replication/derived_metrics_out/{per_cell,per_metric_summary}.csv
        │
        ▼
 replication/plot_all_panels.py  →  replication/figures/*.png,pdf
        │
        ▼
 replication/RESULTS.md  (the pitch)
```

The overnight branch is parallel to this — it re-processes raw ND2s to
get MaxIP metrics alongside z-sum:

```
raw ND2s → replication/metric_pipeline.py (process_cell: z-sum + MaxIP)
         → replication/overnight_out/by_well/**/metrics.csv
         → replication/rebuild_combined.py (joins across sheets correctly)
         → replication/overnight_out/combined.csv
         → replication/evaluate_overnight.py (same stats as above)
         → replication/overnight_eval_out/summary.csv
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
| `replication/stage_missing_wells.sh` | parallel cp SMB → local |
| `replication/rebuild_combined.py` | fix cross-sheet tags in combined.csv |
| `replication/evaluate_overnight.py` | scoring framework over overnight data |
| `replication/plot_all_panels.py` | Fig 4 + S11 side-by-side plots |

---

## Context / provenance gotchas

- **Pairs differ per panel**. The Šídák family Mark uses in each Prism file is different. Fig S11 D/E/4C/4D use m=2 adjacent pairs. Fig 4B / S11 C / TRAK isoform-mito use m=3 all-pairs. Fig 4E peripheral uses m=3 of (A-B, A-C, A-D) against a reference. Fig S11 F perinuclear uses m=3 of (A-C, B-D, C-D). See per-panel comments in `plot_all_panels.py`.
- **MAPK9 Prism column order** is NOT the Comparisons sheet order. Prism A/B/C/D = (ctrl ctrl, MAPK9 ctrl, ctrl Ars, MAPK9 Ars); sheet columns are (ctrl ctrl, ctrl Ars, MAPK9 ctrl, MAPK9 Ars). Lost an hour to this.
- **Plate 9 schema** is the older one (`perinuclear_percent_total` not `perinuclear_5um_percent_total`). `derived_metrics.py` canonicalises old → new internally.
- **Cross-sheet well aliasing**: many wells are referenced by multiple sheets under different condition labels (e.g. plate 3 B03 is both `TRAK isoform (mito) · TRAK1` and `TRAK1 helix muts · T1 wt`). The first overnight_run dropped the duplicate tag; `rebuild_combined.py` fixes this by joining per-cell metrics against the comparisons table.
- **Paper Fig S11 F uses m=3** Šídák (published p-values 0.0827, 0.7397, 0.7738), while the on-disk Prism snapshot uses m=2 (p-values 0.056, 0.63). The published figure is authoritative.

---

## One thing to verify next session

The classical nested-ANOVA p-values in my framework are ~1.1-1.5× Mark's
Prism values on the same pair (cell-count-weighted vs Prism's plate-level
weighting in the SS partition). Significance class agrees on every
comparison I checked. If a reviewer specifically pushes back on p-value
precision, the fix is to switch `nested_oneway_anova` in
`replicate_stats.py` from cell-weighted SS to plate-weighted SS. I didn't
do it because the calls of significance all match — but it's a ~20-line
change if needed.

---

## Optional extensions (not required for the pitch)

1. **Add MaxIP rows to `plot_all_panels.py` per-panel plots.** Currently
   the committed plots show Mark / z-sum diff / z-sum ratio. Adding
   MaxIP diff and MaxIP ratio rows would surface the TRAK1 MaxIP result
   visually. Needs data from `overnight_out/combined.csv` merged into the
   plotting dataframe.
2. **Rerun with plate-weighted nested ANOVA** for exact Prism parity.
3. **MaxIP on peroxisome data (Fig S11 C)** — not processed overnight. If
   MaxIP sharpens the peroxisome distribution as much as it did for
   mitochondria, the paper's null claim for peroxisomes might flip. That
   would be a real finding, not a pitch issue.
4. **Commit the figures as SVG instead of PNG+PDF** for smaller diffs if
   you regenerate often.
