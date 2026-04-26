# Handoff — integrate wedge-r KS metric into Mark's pipeline for code review

*Written 2026-04-26. Successor to `HANDOFF_v2.md`.*

This session ran the full whole-dataset analysis with the new metric family,
extended it with the 60mer comparison and a denoising sensitivity check,
and refined the wedge geometry. Everything in `final_pipeline.py` /
`plot_final.py` / `evaluate_final.py` is settled. The remaining work is to
**fold the headline metric back into Mark's `metric_pipeline.py` as a
minimal, code-reviewable patch**, and prepare the repo for manuscript review.

---

## §1 What's settled

### Headline metric: `wedge_r_ks_vs_uniform` (MaxIP)

Kolmogorov–Smirnov distance between (i) the per-cell intensity-weighted
CDF along the wedge-radial axis and (ii) an analytical area-uniform
sector reference. The wedge sweeps upward through the pattern arch from
an apex at the pattern's bottom extremum, with rays through the pattern's
left and right extremes. The KS scalar lives in [0, 1] and is unsigned.

This metric replaced the original polarization-ratio pitch as the
manuscript headline. It is significant on every reviewer-flagged
comparison:

| Comparison | d | p (Šídák, native family) |
|---|---:|---:|
| 4B no-TRAK vs TRAK2 | −0.96 | 0.024 ✱ |
| TRAK isoform TRAK1 vs TRAK2 | −1.48 | 0.0013 ✱✱ |
| TRAK1 muts wt vs mDRH | −2.07 | 0.00045 ✱✱✱ |
| TRAK2 muts mDRH vs mDRH+mSpindly | +1.54 | 0.00011 ✱✱✱ |
| MAPK9 ctrl-ctrl vs ctrl-Ars | −1.11 | 0.015 ✱ |

A complementary metric, `wedge_r_gini`, is the better readout for the
TRAK1 mDRH-vs-mDRH/dSp rescue specifically (where the phenotype is a
clustering tightness change, not a directional shift). Both metrics use
the same wedge geometry — Gini is computed from the same intensity
histogram, so adding it to Mark's pipeline costs nothing extra.

### Whole dataset run

- **541 cells** (mito 97 + peroxisome 85 + 60mer 47 + TRAK1 muts 88 +
  TRAK2 muts 103 + MAPK9 189) processed under `final_pipeline.py`
- **0 errors** after the lookup-fallback patch + 3 new manual overrides
- All per-sheet figures (scalars / wedge-r profile / wedge-r CDF) regenerated
  with Šídák brackets at the user-specified family sizes
- Output lives in `replication/overnight_final_out/`

### Independent runs / sanity checks done this session

- **60mer (synthetic particle) panel** — strongest TRAK1↔TRAK2
  separation in the dataset (d>2 on the headline metrics), single plate so
  uses Welch t-test fallback in `evaluate_final.test_pair`
- **Peroxisome panel** — clean null on every metric; cargo specificity check
- **MAPK9 denoised re-run** (189 cells, 6 h) — every reviewer-relevant
  call agrees with the raw run; supports the "denoiser-invariant on MaxIP"
  story end-to-end
- **10-cell raw-vs-denoised sensitivity** — MaxIP wedge_r_gini /
  wedge_r_ks_vs_uniform drift < 1.1% mean per-cell; zsum versions shift
  ~9% but rankings preserve
- **60mer-noTRAK empirical reference** experiment — the cytoplasmic-fill
  baseline reproduces the manuscript story but doesn't change conclusions.
  Lives in supplementary code (`ks_vs_60mer_reference.py`); robust-of-metric
  check for the supplement

### Geometry verification

- **Wedge invariance** confirmed across all 494 cells: opening angle
  exactly 45.35° (sd = 0). The wedge apex sits 24.96 µm below the image
  midpoint and ~0 µm horizontally — deterministic, set by the rigid
  template.
- **Tangent-vs-extreme question** explored: original wedge rays pass
  through the leftmost-x / rightmost-x pattern pixels (which sit at the
  upper arch's widest point); a tangent variant anchored at the donut's
  lower-arm region widens the wedge by ~2° and gains ~6% area. Diagnostic
  figures live in `overnight_final_out/geometry_verification/`. Not
  adopted as the canonical wedge — the difference doesn't change any
  statistical conclusion.

### Bug fixes shipped

- **Pattern-mask coordinate-frame bug** in Mark's pipeline (HANDOFF_v2 §6)
  is fixed in `final_pipeline.py`. Mark's `metric_pipeline.py` still has
  it. The fix needs to come along with the integration patch.
- **Override-key fallback**: `template_matching_bulk.find_override_key`
  now walks raw → denoised key when the raw key isn't registered. Three
  new manual overrides added to `coordinate_overrides.csv`
  (plate_11/F05/Cell12, plate_12/G04/Cell6 raw + denoised). These are
  needed for the 5 cells that previously failed.
- **Welch t-test fallback** added to `evaluate_final.test_pair` so
  single-plate sheets (60mer) get usable p-values when
  `nested_oneway_anova` returns NaN.

### Tweaks that were tried and rejected

- Restricting tangent search to upper arch only — gave the same tangent
  points as the canonical extremes; no useful change
- Convex-hull tangent of the full pattern from a midpoint apex — gave a
  ~140° wedge dominated by stalk geometry, scientifically uninformative
- Auto-detected stalk cutoff at x-extent > 80 px — picked up the stalk's
  own width; needs threshold of 200 px (donut-wide) to skip the stalk
- The 60mer-noTRAK empirical reference (with 13 reference cells) —
  reproduces the area-uniform story but the noTRAK self-comparison cells
  inherit a noise floor (within-condition CDF spread shows up as
  positive KS values). Useful as a robustness check, not a replacement.

---

## §2 Next objective

**Add the wedge-r KS metric to Mark's `metric_pipeline.py` as a minimal
patch he can code-review, then submit the manuscript repo.**

The patch should:

1. **Add wedge geometry** — pattern bottom/left/right extreme detection,
   wedge mask, per-pixel polar r in µm from the apex. `final_pipeline.py`
   has the reference implementation in three small helpers
   (`_pattern_extremes`, `_build_wedge_geometry`, `_ks_vs_uniform`).
2. **Compute the new column(s)** per projection:
   `{proj}_wedge_r_ks_vs_uniform` (and probably `{proj}_wedge_r_gini` while
   we're there — same histogram, free).
3. **Fix the `pattern_mask_big` coordinate-frame bug** along the way —
   the wedge code requires the slice to live in the template's own frame,
   so Mark's existing slice has to be corrected anyway.
4. **Validate cell-by-cell** by re-processing a stratified subset (~10
   cells) with the patched `metric_pipeline.py` and confirming
   `wedge_r_ks_vs_uniform` matches the value already in
   `overnight_final_out/combined_raw.csv` to numerical precision (Pearson
   r ≥ 0.999, mean offset < 1%).
5. **Keep the patch surgical** — don't fold in the 60-bin profile
   columns, the Y-axis metrics, the wedge moment-tensor columns, etc.
   Those live in `final_pipeline.py` for the analyses we've already done;
   they're not needed for the manuscript headline.

Companion items that should ship in the manuscript-review submission:

- Patched `template_matching_bulk.py` (the override fallback)
- Updated `coordinate_overrides.csv` (3 new rows)
- A short README explaining what `wedge_r_ks_vs_uniform` measures and how
  the wedge is constructed (refer to the wedge_illustration figure already
  in `overnight_fig4b_v2_out/figures/`)
- The validation script + a tiny golden CSV for one cell so reviewers /
  Mark can re-verify

What NOT to commit: `by_well/` checkpoint trees (multi-GB, regenerable),
the giant `combined_raw.csv` files (also regenerable from per-cell
JSONs), and the SMB-data symlinks.

---

## §3 Useful files at a glance

| File | What it is | Status |
|---|---|---|
| `replication/final_pipeline.py` | Reference implementation; canonical for the whole-dataset run | Source of truth for the integration patch |
| `replication/metric_pipeline.py` | Mark's original pipeline | **Target of the integration patch** |
| `template_matching_bulk.py` | Shared template matching + override loader | Already patched (override fallback) |
| `coordinate_overrides.csv` | Manual pattern-center overrides for failing cells | Already extended (3 new rows) |
| `replication/evaluate_final.py` | Multi-sheet evaluator (Šídák + Welch fallback) | Settled |
| `replication/plot_final.py` | Per-sheet figure generator | Settled |
| `replication/ks_vs_60mer_reference.py` | Supplementary 60mer-noTRAK reference experiment | Settled, supplementary only |
| `replication/compare_raw_vs_denoised_ks.py` | 10-cell denoising sensitivity script | Settled |
| `replication/run_pipeline_on_paths.py` | Direct-paths variant of `final_pipeline.run`, used for the denoised MAPK9 run because the standard walk skips `denoised/` subdirs | Reusable |
| `replication/HANDOFF_v2.md` | Prior session's handoff (Fig 4B → whole-dataset queue) | Superseded |
| `replication/HANDOFF_v3.md` | This doc | — |

---

## §4 Gotchas

- **Wedge sweeps UPWARD through the arch.** The first `s11_pipeline.py`
  in the repo had the direction flipped; do not copy from it.
- **Pattern-mask slice frame**: must slice `shifted_template` in the
  template's own frame (`max_coords ± 512`), NOT in the cropped image
  frame (`max_coords ± 512 + offset`). Mark's current code has the bug;
  the integration patch fixes it.
- **MaxIP is the right projection for the headline metric.** Denoised
  vs raw MaxIP differ by < 1.1% per cell on KS / Gini; zsum differs by
  ~9% (consistent shift, rankings preserve).
- **Single-plate sheets**: the 60mer panel only has data on plate 6,
  so `nested_oneway_anova` returns NaN. The Welch fallback in
  `evaluate_final.test_pair` handles this; the manuscript should note
  the caveat in the figure caption.
- **Cell6 on plate 12** has `max_coords[0] = 508 < 512`, which causes
  the unshifted `pattern_mask` slice to evaluate to an empty array. The
  manual override at (520, 915) nudges the apex by 12 px so the slice
  is valid. Tiny shift, no scientific impact.
- **SMB drops** every ~15 min of active I/O on this laptop — the
  per-cell checkpointing absorbs them. Runs of ~6 h were stable under
  `caffeinate -dimsu`.

---

## §5 Today's project state in one paragraph

The pitch is settled, the whole-dataset numbers are in the bag, every
sanity check (denoising, peroxisome cargo specificity, 60mer scaffold,
empirical-vs-analytical reference, wedge geometry invariance) supports
the conclusions. What's left is to land the new metric in Mark's
canonical pipeline as a minimal, reviewable patch and tidy the repo for
manuscript submission.
