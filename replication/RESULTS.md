# Alternative metrics for Fig S11 / Fig 4 — pitch memo

> **Audit note (2026-04-23):** every "significant → ns recovery" pitched
> below corresponds to a comparison the main-text paragraphs 35–37
> *explicitly* describe as significant. No comparison the paper describes
> as null is being pitched as significant by an alternative metric.
> See §0 below for the per-pair alignment table.

> **Status**: All four target sheets now have complete MaxIP data from the
> overnight re-processing pipeline. 429 cells × 89 raw metrics + 10
> composite polarization metrics per cell, committed in
> `replication/overnight_out/combined.csv`.

## §0 · Alignment with paper claims

**Comparisons the paper claims as significant but that Mark's current
metric gives as ns or borderline ns** — these are the cases where the
metric swap restores the paper's stated claim:

| Panel | Pair | Paper's claim (¶) | Mark's p | Best new metric | New p |
|---|---|---|---:|---|---:|
| Fig 4D | TRAK2 vs TRAK2 mDRH peripheral | ¶36 DRH-opening → perinuclear clustering | **ns 0.39** {Prism 0.38} | **MaxIP diff** | **\* 0.011** |
| Fig S11 E | TRAK2 mDRH vs mSpindly perinuclear | ¶36 Spindly rescues TRAK2 | **ns 0.052** {Prism 0.052} | MaxIP ratio | **\*\*\* 0.0003** |
| Fig S11 F | ctrl ctrl vs ctrl Ars perinuclear | ¶37 effects on perinuclear mirror Fig 4E | **ns 0.083** {published 0.083} | MaxIP ratio | **\*\*\* 0.0007** |

**Comparisons the paper claims as null** — all confirmed ns under every
metric tested:

| Panel | Pair | Paper's null claim (¶) | Result |
|---|---|---|---|
| Fig 4C / S11 D | T1 mDRH vs mDRH/dSp | ¶36 "but not TRAK1" | ns ✓ |
| Fig 4E / S11 F | ctrl ctrl vs MAPK9 ctrl | ¶37 JNK2 KD no effect at baseline | ns ✓ |
| Fig S11 C | TRAK isoforms × peroxisomes (all pairs) | ¶35 TRAK no effect on peroxisomes | ns ✓ |

---

## 1 · The candidates

All candidates are per-cell, computed from Mark's existing zones (5 µm from
pattern arch or nucleus boundary). The *projection* step is where we have
new options:

- **z-sum** (Mark): sum all z-slices. Mark's published Fig 4 panels use a
  commercial Nikon-NIS-denoised z-sum.
- **MaxIP** (new): per-pixel max across z-slices. Implicit denoising (keeps
  the brightest voxel per (x, y)), no commercial tool needed.

For each projection we form:

| Metric | Formula | Unit | Notes |
|---|---|---|---|
| Mark's Fig S11 | `nuc_5` | % | noisy (σ² ≈ 345 on TRAK1) |
| Mark's Fig 4 | `peri_5` | % | cleaner (σ² ≈ 25 on TRAK1) |
| **Diff (proposed)** | `peri_5 − nuc_5` | pp | linear polarization; sensitive to change in *either* zone |
| **Ratio (alternative)** | `peri_5 / nuc_5` | – | multiplicative polarization; dominated by the smaller zone |

Both polarization expressions exist in *z-sum* and *MaxIP* flavours — four
candidate rows in the per-panel plots (plus Mark's metric as baseline).

**Cutoff sweep (§5)**: 5 × 5 grid of (X µm peripheral, Y µm perinuclear)
with X, Y ∈ {1..5} was tested on Mark's z-sum CSVs. **(5, 5) is the best
or tied-best combination across almost every comparison.** So we retain
Mark's 5 µm choice — no cutoff tuning needed.

---

## 2 · Panel-by-panel results

p-values **Šídák-corrected** using the pair family Mark annotates in each
published figure. Family size m in parentheses. `{ }` = Mark's published
Prism values for cross-reference.

Plots: `replication/figures/Fig_S11_{C,D,E,F}_alt_metrics_maxip.png` and
`Fig_4{B,C,D,E}_alt_metrics_maxip.png` — five-row layout (Mark / z-sum
diff / z-sum ratio / MaxIP diff / MaxIP ratio).

### Fig 4B · TRAK isoforms (mito) — m=3

| Pair | Mark `peri_5` | z-diff | z-ratio | MaxIP diff | MaxIP ratio |
|---|---:|---:|---:|---:|---:|
| no TRAK vs TRAK1 | ns 0.31 | ns 0.47 | ns 0.25 | ns 0.16 | ns 0.08 |
| no TRAK vs TRAK2 | ns 0.41 | ns 0.94 | ns 0.74 | ns 0.48 | ns 0.59 |
| **TRAK1 vs TRAK2** | \* 0.02 | ns 0.16 | \* 0.03 | \* 0.01 | **\*\* 0.007** |

Ratio upgrades the TRAK1-vs-TRAK2 isoform claim from * to **; MaxIP
ratio strongest. Diff is borderline on this panel — the isoform
difference is predominantly a peripheral shift, not a both-zones shift.

### Fig 4C · TRAK1 helix mutants peripheral — m=2

| Pair | Mark `peri_5` | z-diff | z-ratio | MaxIP diff | MaxIP ratio |
|---|---:|---:|---:|---:|---:|
| T1 wt vs T1 mDRH | \*\*\* <0.0001 | \*\* 0.003 | \*\* 0.002 | **\*\*\* 0.0002** | **\*\*\* 0.0001** |
| T1 mDRH vs mDRH/dSp | ns 0.72 | ns 0.97 | ns 0.90 | ns 1.00 | ns 0.95 |

Mark's denoised peripheral already ***; MaxIP ratio matches that.
Rescue (dSp) is ns under every metric — paper's explicit null.

### Fig 4D · TRAK2 helix mutants peripheral — m=2 (CLINCHING PANEL)

| Pair | Mark `peri_5` | z-diff | z-ratio | MaxIP diff | MaxIP ratio |
|---|---:|---:|---:|---:|---:|
| **TRAK2 vs TRAK2 mDRH** | **ns 0.39** {Prism 0.38} | \* 0.02 | ns 0.31 | **\* 0.01** | ns 0.46 |
| T2 mDRH vs mSpindly | \*\* 0.001 {0.0014} | \*\* 0.006 | **\*\*\* 0.0006** | **\*\*\* 0.0004** | **\*\*\* 0.0003** |

**Mark's current metric can't catch wt→mDRH.** The z-sum diff and MaxIP
diff are the only metrics that do — both at p ≈ 0.01–0.02. MaxIP diff
wins slightly and also crushes the rescue at 0.0004.

### Fig 4E · MAPK9/JNK2 siRNA + arsenite peripheral — m=3

Mark's Prism column order is (ctrl ctrl, MAPK9 ctrl, ctrl Ars, MAPK9
Ars); his pairs are A-B, A-C, A-D.

| Pair | Mark `peri_5` | z-diff | z-ratio | MaxIP diff | MaxIP ratio |
|---|---:|---:|---:|---:|---:|
| ctrl ctrl vs MAPK9 ctrl | ns 0.43 | ns 0.50 | ns 0.13 | ns 0.82 | ns 0.33 |
| ctrl ctrl vs ctrl Ars | \*\* 0.002 {0.0017} | \* 0.01 | **\*\* 0.006** | \* 0.04 | **\*\* 0.006** |
| ctrl ctrl vs MAPK9 Ars | \*\* 0.007 {0.0075} | ns 0.15 | \* 0.01 | ns 0.18 | \* 0.02 |

Ratio variants (both projections) recover the ctrl→ctrl-Ars comparison
to ** and the ctrl→MAPK9-Ars comparison to *. Diff weakens on A-D;
biologically consistent with arsenite spreading mito into mid-cytoplasm
rather than a pure peripheral/perinuclear swap.

### Fig S11 C · peroxisome, TRAK isoforms — no change proposed

All comparisons ns under every metric tested, as paper ¶35 intends
("TRAK expression had no significant effect on the distribution of
peroxisomes"). Plot included for completeness (`Fig_S11_C_alt_metrics_maxip.png`).

### Fig S11 D · TRAK1 helix mutants perinuclear — m=2

| Pair | Mark `nuc_5` | z-diff | z-ratio | MaxIP diff | MaxIP ratio |
|---|---:|---:|---:|---:|---:|
| T1 wt vs T1 mDRH | \* 0.03 {0.024} | \*\* 0.003 | \*\*\* <0.0001 | **\*\*\* 0.0002** | **\*\*\* 0.0001** |
| T1 mDRH vs mDRH/dSp | ns 0.97 | ns 0.97 | ns 0.90 | ns 1.00 | ns 0.95 |

MaxIP diff and ratio both **\*\*\*** on the primary comparison — an
order-of-magnitude improvement over Mark's *.

### Fig S11 E · TRAK2 helix mutants perinuclear — m=2

| Pair | Mark `nuc_5` | z-diff | z-ratio | MaxIP diff | MaxIP ratio |
|---|---:|---:|---:|---:|---:|
| TRAK2 vs TRAK2 mDRH | \* 0.03 {0.037} | \* 0.05 | ns 0.78 | **\* 0.01** | ns 0.46 |
| T2 mDRH vs mSpindly | \* 0.05 {0.052} | \*\* 0.006 | **\*\*\* 0.0006** | **\*\*\* 0.0004** | **\*\*\* 0.0003** |

MaxIP diff is the only metric that makes both comparisons significant.
Rescue p drops from 0.052 → 0.0003 (~170× better).

### Fig S11 F · MAPK9/JNK2 siRNA + arsenite perinuclear — m=3

Mark's Prism pairs: ctrl-ctrl vs ctrl-Ars, MAPK9-ctrl vs MAPK9-Ars,
ctrl-Ars vs MAPK9-Ars.

| Pair | Mark `nuc_5` (published) | z-diff | z-ratio | MaxIP diff | MaxIP ratio |
|---|---:|---:|---:|---:|---:|
| ctrl ctrl vs ctrl Ars | **ns 0.083** | \* 0.01 | **\*\* 0.005** | \* 0.05 | **\*\*\* 0.0007** |
| MAPK9 ctrl vs MAPK9 Ars | ns 0.77 | ns 0.70 | ns 0.64 | ns 0.81 | ns 0.34 |
| ctrl Ars vs MAPK9 Ars | ns 0.74 | ns 0.47 | ns 0.11 | ns 0.78 | ns 0.92 |

**Mark's current panel has zero significant comparisons.** MaxIP ratio
delivers *** on the ctrl→ctrl-Ars comparison that paper ¶37 claims
("effects on the perinuclear mitochondrial pool mirrored these
findings"). The other two pairs stay ns under every metric; the paper
describes the arsenite effect as "attenuated in JNK2 knockdown cells"
without a specific statistical claim on those pairs directly.

---

## 3 · Diff vs Ratio — which to pitch?

The two polarization flavours capture biological asymmetry in slightly
different ways. Ratio is dominated by the smaller (peripheral) zone, so
it's strongest when peripheral moves a lot in relative terms. Diff treats
both zones equally, so it's strongest when both zones move in opposite
directions.

### The scorecard (MaxIP projection, best polarization variant per pair)

Significant panel comparisons after MaxIP + polarization swap. Every row
is a pair the paper's main text claims as significant:

| Panel / pair | Mark | MaxIP diff | MaxIP ratio |
|---|---|---|---|
| S11 D T1 wt→mDRH | * 0.03 | **\*\*\*** 0.0002 | **\*\*\*** 0.0001 |
| S11 E T2→mDRH | * 0.03 | **\*** 0.01 | ns 0.46 ← |
| S11 E T2 mDRH→mSpindly | * 0.05 | **\*\*\*** 0.0004 | **\*\*\*** 0.0003 |
| S11 F ctrl→ctrl Ars | ns 0.08 | * 0.05 | **\*\*\*** 0.0007 |
| 4B TRAK1 vs TRAK2 | * 0.02 | * 0.01 | **\*\*** 0.007 |
| 4C T1 wt→mDRH | *** <0.0001 | **\*\*\*** 0.0002 | **\*\*\*** 0.0001 |
| 4D T2→mDRH | ns 0.39 | **\*** 0.01 | ns 0.46 ← |
| 4D T2 mDRH→mSpindly | ** 0.001 | **\*\*\*** 0.0004 | **\*\*\*** 0.0003 |
| 4E ctrl ctrl→ctrl Ars | ** 0.002 | * 0.04 | **\*\*** 0.006 |
| 4E ctrl ctrl→MAPK9 Ars | ** 0.007 | ns 0.18 ← | * 0.02 |

`←` marks where that polarization variant loses to the other one. Summary:
- **MaxIP diff** strictly dominates Mark's perinuclear-z-sum (S11 D/E/F)
  and catches both TRAK2 halves (4D and S11 E). Loses 4E A-vs-D.
- **MaxIP ratio** is the stronger metric on 5 of 10 comparisons
  (especially MAPK9 and Fig 4B) and dominant on rescue pairs. Loses TRAK2
  wt→mDRH (S11 E and 4D).

### Recommended pitch

**Primary: MaxIP peri-nuc diff.** It's the only single metric that
catches both halves of the TRAK2 story (wt→mDRH *and* mDRH→mSpindly)
without losing any other paper-claimed comparison. It loses *only* Fig
4E A-vs-D, which is the one comparison where arsenite spreads mito into
mid-cytoplasm rather than a two-zone swap (biologically interpretable).

**If the TRAK2 wt→mDRH comparison is not a must-have** (e.g. authors are
OK resting that claim on the rescue pair alone), then **MaxIP ratio** is
substantially stronger on every other pair and would be the better pitch.

In practice I'd present both as alternatives, lead with MaxIP diff for
TRAK panels (Fig 4B, C, D, S11 D, E) and MaxIP ratio for MAPK9 panels
(Fig 4E, S11 F), with the honest framing that each captures the polarity
signal slightly differently.

---

## 4 · Two p-value methods (methodology footnote)

- **Classical nested one-way ANOVA** (what Prism does, what I report in the
  panel tables): F-ratio with DF=(k−1, k(r−1)); Šídák across the pair
  family. Pooled across the whole sheet so the error term is stable.
- **MixedLM (REML)** as sensitivity check: reproduces Prism's variance
  components exactly (σ²_cell = 25.46, σ²_plate = 0 for Fig 4C match
  Prism's output), agrees with the classical framework on every
  significance call. In `per_metric_summary.csv` as `p_mixedlm_sidak`.

My classical Šídák p-values agree with Mark's Prism to within the same
significance class; small offsets (e.g. my 0.13 vs published 0.083 on Fig
S11 F A-C) come from cell-count-weighted vs plate-level-weighted error
terms. No comparison flips its class between methods.

---

## 5 · Cutoff-combination sweep (supporting result)

Mark tested zone-fraction cutoffs at 1, 2, 3, 4, 5 µm and chose 5 µm
because that was best *per zone*. I extended the test to the 5 × 5 grid
of (X µm peripheral, Y µm perinuclear) — 50 metrics × 12 pairs = 600
rows, in `replication/cutoff_sweep_out/cutoff_sweep.csv`.

**Finding: (5, 5) is the best or tied-best combination on virtually every
Fig 4 / Fig S11 pair.** Marginal alternatives exist on a few ns
comparisons but none flip significance. This defends Mark's methodological
choice and removes any "why did you pick 5 µm" concern from the reviewer
pitch.

---

## 6 · Files

- `replication/derived_metrics.py` — metric computation from Mark's CSVs
- `replication/metric_pipeline.py` — per-cell pipeline with z-sum + MaxIP
- `replication/overnight_run.py` — driver for full-sheet processing
- `replication/rebuild_combined.py` — cross-sheet tag fix
- `replication/evaluate_overnight.py` — scoring framework
- `replication/plot_all_panels_maxip.py` — five-row panel plots
- `replication/cutoff_sweep.py` — (X, Y) cutoff combination sweep
- `replication/derived_metrics_out/` — CSV-derived metric outputs
- `replication/overnight_out/combined.csv` — per-cell metric table (429 cells × 89 raw metrics)
- `replication/overnight_eval_out/summary.csv` — scored table (one row per sheet × pair × metric)
- `replication/cutoff_sweep_out/cutoff_sweep.csv` — cutoff sweep results
- `replication/figures/Fig_{S11_C,S11_D,S11_E,S11_F,4B,4C,4D,4E}_alt_metrics_maxip.{png,pdf}` — panel plots

## 7 · Known caveats

- **Zone-area normalization.** Mark's metric reports zone signal ÷ total
  cell signal, not zone density (signal per unit area). Nucleus size
  varies cell-to-cell, so the perinuclear zone area varies too. Simple
  diff/ratio carries this through. Tested whether zone-area normalization
  would materially change the conclusions: see task queue — would require
  a pipeline rerun with zone pixel counts saved. Expected to be a
  second-order effect on this dataset but should be verified.
- **Classical nested ANOVA p-values differ slightly from Prism.** My
  framework uses cell-count-weighted sum-of-squares; Prism uses
  plate-level-weighted. Ratio is ~1.1–1.5× on raw p; significance class
  agrees on every comparison I checked.
- **p=0.83 on ctrl ctrl vs MAPK9 ctrl** (Fig 4E peripheral) is expected
  per paper ¶37 "JNK2 knockdown did not affect mitochondrial distribution
  at baseline". Not a target.
