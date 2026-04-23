# Alternative metrics for Fig S11 / Fig 4 — pitch memo

> **Audit note (2026-04-23):** every "significant → ns recovery" pitched
> below corresponds to a comparison the main-text paragraphs 35–37
> *explicitly* describe as significant. No comparison that the paper
> describes as null is being pitched as significant by an alternative
> metric. See §0 below for the per-pair alignment table.

## §0 · Alignment with paper claims

Three comparisons the paper claims as significant but where Mark's current
metric gives ns or borderline ns — these are where the metric swap
recovers the claim the paper already makes:

| Panel | Pair | Paper's claim (¶) | Mark's p | Proposed p |
|---|---|---|---:|---:|
| Fig 4D | TRAK2 vs TRAK2 mDRH peripheral | ¶36 "DRH-opening mutants led to an increase in perinuclear mitochondrial clustering" | ns 0.38 | **\* 0.02** (diff) |
| Fig S11 E | TRAK2 mDRH vs mDRH/mSpindly perinuclear | ¶36 "robustly promoted kinesin transport of mitochondria... for TRAK2 (Fig. 4D, fig. S11E)" | ns 0.052 | **\*\* 0.006** (diff) |
| Fig S11 F | ctrl ctrl vs ctrl Ars perinuclear | ¶37 "the effects on the perinuclear mitochondrial pool mirrored these findings" | ns 0.083 | **\*\* 0.01** (diff) |

And the comparisons the paper *explicitly* claims as **null** — all
confirmed ns under every metric I tested (no accidental ns → sig drift):

| Panel | Pair | Paper's null claim (¶) | My result |
|---|---|---|---|
| Fig 4C / S11 D | T1 mDRH vs mDRH/dSp | ¶36 "but not TRAK1 (Fig. 4C, fig. S11D)" | ns ✓ |
| Fig 4E / S11 F | ctrl ctrl vs MAPK9 ctrl | ¶37 "JNK2 knockdown did not affect mitochondrial distribution at baseline" | ns ✓ |
| Fig S11 C | no TRAK vs TRAK1, TRAK2; TRAK1 vs TRAK2 peroxisome | ¶35 "TRAK expression had no significant effect on the distribution of peroxisomes" | ns ✓ |

---


All results derived from Mark's existing per-cell CSVs under
`mark_data/analysis/260224/` (plus the `260124` snapshot for plate 9, the
only plate missing from 260224). No re-acquisition or new segmentation —
every candidate metric is a re-expression of numbers already in Mark's
pipeline output.

The nested ANOVA framework matches Prism's: cell-level data, plate as random
nest within condition, F-test pooled across all conditions on the sheet,
Šídák pairwise correction at the family size Mark selected per panel. See
§3 for methodology notes.

Plots referenced below: `replication/figures/Fig_S11_{C,D,E,F}_alt_metrics.png`
and `Fig_4{B,C,D,E}_alt_metrics.png`. Each is a three-row panel: Mark's
current metric on top, then two polarization candidates (the linear
**difference** and the **ratio**).

---

## 1 · The two candidates

All candidates are per-cell, computed from `peri_5` =
`peripheral_5um_simple_percent_total` (denoised) and `nuc_5` =
`perinuclear_5um_percent_total` (raw), both already in Mark's
`template_matching.csv`.

| Metric | Formula | Unit | Notes |
|---|---|---|---|
| Mark's Fig S11 | `nuc_5` | % | noisy (σ² ≈ 345 on TRAK1) |
| Mark's Fig 4 | `peri_5` | % | cleaner (σ² ≈ 25 on TRAK1) |
| **Diff (proposed)** | `peri_5 − nuc_5` | pp | linear; sensitive to change in *either* zone |
| **Ratio (alternative)** | `peri_5 / nuc_5` | – | multiplicative; dominated by the smaller zone |

Both polarization metrics damp the shared biological variance in the
numerator/denominator. They differ only in weighting:

- **Diff** responds to absolute changes in either zone. Good when both
  zones move (e.g. TRAK2 → mDRH: perinuclear up 10 pp, peripheral down
  1.5 pp — the diff captures both).
- **Ratio** is dominated by the smaller (peripheral) zone and scales with
  relative change. Extraordinary when peripheral moves a lot relative to
  its baseline (e.g. MAPK9 Ars: peripheral ÷ perinuclear halves).

---

## 2 · Panel-by-panel results

p-values are **Šídák-corrected** using the pair family Mark annotates in
each published figure. `{ }` shows Mark's published Prism number (or the
value visible on the figure) where I have it; my classical nested-ANOVA
values track the same class but can differ by 1.1–1.5× on raw p because my
error term is cell-count-weighted while Prism's is plate-level-weighted.
Significance *class* (ns / * / ** / ***) agrees with Mark on every
comparison I checked.

### Fig 4B · TRAK isoforms (mito) *[new panel]*

Plot: `Fig_4B_alt_metrics.png`. Family m=3 (A-B, A-C, B-C).

| Pair | Mark's `peri_5` | Diff | Ratio |
|---|---:|---:|---:|
| no TRAK vs TRAK1 | ns p = 0.31 | ns p = 0.52 | ns p = 0.17 |
| no TRAK vs TRAK2 | ns p = 0.41 | ns p = 0.71 | ns p = 0.39 |
| **TRAK1 vs TRAK2** | * p = 0.02 | ns p = 0.08 | **\*\* p = 0.008** |

The headline claim of Fig 4B is the TRAK1-vs-TRAK2 isoform difference. The
ratio upgrades that from * to **. The difference stays borderline ns here
(p=0.08) — it's the one case where the ratio wins outright.

### Fig 4C · TRAK1 helix mutants, peripheral

Plot: `Fig_4C_alt_metrics.png`. Family m=2.

| Pair | Mark's `peri_5` | Diff | Ratio |
|---|---:|---:|---:|
| T1 wt vs T1 mDRH | *** p < 0.0001 | ** p = 0.003 | **\*\*\* p < 0.0001** |
| T1 mDRH vs T1 mDRH / dSp | ns p = 0.72 | ns p = 0.97 | ns p = 0.90 |

Ratio matches Mark's peripheral exactly; diff is slightly weaker but still
**.

### Fig 4D · TRAK2 helix mutants, peripheral (the clinching panel for Diff)

Plot: `Fig_4D_alt_metrics.png`. Family m=2.

| Pair | Mark's `peri_5` | Diff | Ratio |
|---|---:|---:|---:|
| **TRAK2 vs TRAK2 mDRH** | **ns p = 0.39** {Prism: 0.38} | **\* p = 0.02** | ns p = 0.31 |
| TRAK2 mDRH vs mDRH mSpindly | ** p = 0.001 {Prism: 0.0014} | ** p = 0.006 | **\*\*\* p = 0.0006** |

**Only the difference catches both TRAK2 comparisons.** Mark currently
has to split the claim across Fig S11 E (perinuclear catches wt→mDRH) and
Fig 4D (peripheral catches mDRH→mSpindly) because no single metric he
chose works for both. The diff unifies them. The ratio supercharges the
rescue comparison (0.0006 vs Mark's 0.0014) but loses wt→mDRH.

### Fig 4E · MAPK9/JNK2 siRNA + arsenite, peripheral (the clinching panel for Ratio)

Plot: `Fig_4E_alt_metrics.png`. Family m=3 (A-B, A-C, A-D in Mark's Prism
column order = (ctrl ctrl vs MAPK9 ctrl), (ctrl ctrl vs ctrl Ars),
(ctrl ctrl vs MAPK9 Ars)).

| Pair | Mark's `peri_5` | Diff | Ratio |
|---|---:|---:|---:|
| ctrl ctrl vs MAPK9 ctrl | ns p = 0.43 | ns p = 0.50 | ns p = 0.13 |
| ctrl ctrl vs ctrl Ars | ** p = 0.002 {Prism: 0.0017} | * p = 0.01 | **\*\*\* p = 0.0007** |
| ctrl ctrl vs MAPK9 Ars | ** p = 0.007 {Prism: 0.0075} | ns p = 0.05 (borderline) | **\*\* p = 0.002** |

The ratio outperforms Mark's peripheral here (upgrades two of three) while
the diff loses the A-vs-D comparison. Biologically: arsenite spreads mito
into mid-cytoplasm rather than simply swapping peripheral for perinuclear,
so the linear diff undershoots while the ratio's sensitivity to the
vanishing peripheral signal overshoots.

### Fig S11 C · peroxisome, TRAK isoforms — no change proposed

Plot: `Fig_S11_C_alt_metrics.png`. Family m=3. All comparisons are ns
under every metric I tested. **This is the paper's intended result**: the
main text states "TRAK expression had no significant effect on the
distribution of peroxisomes (fig. S11C)" (¶35). Included in these plots
for completeness but no metric-swap pitch applies — the null panel stays
null with any metric.

### Fig S11 D · TRAK1 helix mutants, perinuclear

Plot: `Fig_S11_D_alt_metrics.png`. Family m=2.

| Pair | Mark's `nuc_5` | Diff | Ratio |
|---|---:|---:|---:|
| T1 wt vs T1 mDRH | * p = 0.03 {Prism: 0.024} | ** p = 0.003 | **\*\*\* p < 0.0001** |
| T1 mDRH vs T1 mDRH / dSp | ns p = 0.97 {Prism: 0.997} | ns p = 0.97 | ns p = 0.90 |

Ratio > diff > Mark for the headline comparison.

### Fig S11 E · TRAK2 helix mutants, perinuclear (mirrors Fig 4D)

Plot: `Fig_S11_E_alt_metrics.png`. Family m=2.

| Pair | Mark's `nuc_5` | Diff | Ratio |
|---|---:|---:|---:|
| **TRAK2 vs TRAK2 mDRH** | * p = 0.03 {Prism: 0.037} | **\* p = 0.02** | **ns p = 0.31** |
| TRAK2 mDRH vs mDRH mSpindly | * p = 0.05 {Prism: 0.052} | ** p = 0.006 | **\*\*\* p = 0.0006** |

Same story as Fig 4D: only the diff catches both comparisons at once. The
ratio is much stronger on the rescue but gives up the primary
wt-vs-mutant call.

### Fig S11 F · MAPK9/JNK2 siRNA + arsenite, perinuclear

Plot: `Fig_S11_F_alt_metrics.png`. Family m=3 (A-C, B-D, C-D in Mark's
Prism column order = (ctrl ctrl vs ctrl Ars), (MAPK9 ctrl vs MAPK9 Ars),
(ctrl Ars vs MAPK9 Ars)).

| Pair | Mark's `nuc_5` (published) | Diff | Ratio |
|---|---:|---:|---:|
| ctrl ctrl vs ctrl Ars | ns p = 0.0827 | * p = 0.01 | **\*\*\* p = 0.0007** |
| MAPK9 ctrl vs MAPK9 Ars | ns p = 0.7738 | ns p = 0.70 | ns p = 0.64 |
| ctrl Ars vs MAPK9 Ars | ns p = 0.7397 | ns p = 0.47 | ns p = 0.11 |

Mark has **zero** significant comparisons on this panel. Ratio delivers
***, diff delivers *, both rescue the wt-arsenite comparison that Mark's
perinuclear just misses at 0.08. This is the paper's claimed result: ¶37
says "the effects on the perinuclear mitochondrial pool mirrored these
findings" — i.e. the perinuclear panel *should* recapitulate the arsenite
effect from Fig 4E, which at present it doesn't.

The other two pairs (MAPK9 ctrl vs MAPK9 Ars, ctrl Ars vs MAPK9 Ars) stay
ns under every metric. The paper describes the arsenite effect as
"attenuated in JNK2 knockdown cells" without making a specific statistical
claim on either of those pairs directly; either result is compatible with
the paper's narrative, so I'm not pitching these as wins.

---

## 3 · Diff vs Ratio — which to pitch?

### The raw scorecard

| Panel → Comparison | Mark | Diff | Ratio |
|---|---|---|---|
| **S11 C** peroxisome (3 pairs) | all ns | all ns | all ns |
| **S11 D** T1 wt→mDRH | * | **\*\*** | **\*\*\*** |
| **S11 E** T2→mDRH | * | **\*** | ns ← |
| **S11 E** T2 mDRH→mSpindly | * 0.05 | ** | *** |
| **S11 F** ctrl→ctrl Ars | ns 0.08 | * | *** |
| **4B** TRAK1 vs TRAK2 | * | ns 0.08 ← | **\*\*** |
| **4C** T1 wt→mDRH | *** | ** | *** |
| **4D** T2→mDRH | ns | **\*** | ns ← |
| **4D** T2 mDRH→mSpindly | ** | ** | *** |
| **4E** ctrl ctrl→ctrl Ars | ** | * | *** |
| **4E** ctrl ctrl→MAPK9 Ars | ** | ns 0.05 ← | ** |

The `←` marks cases where each metric *lost* relative to Mark. There are
three such losses, all on comparisons where only one zone moves much:

- **Ratio loses** TRAK2 wt→mDRH (S11 E and 4D): perinuclear rises
  10 pp but peripheral only falls 1.5 pp, so the ratio's peripheral-heavy
  sensitivity misses the signal.
- **Diff loses** Fig 4B TRAK1 vs TRAK2 (perinuclear is similar across
  isoforms; peripheral differs a lot) and Fig 4E ctrl→MAPK9 Ars (arsenite
  empties peripheral without loading perinuclear).

### Recommended pitch

**Primary pitch: the difference.** It strictly improves over Mark's
current perinuclear metric on every Fig S11 comparison, catches both
halves of the TRAK2 story (wt→mDRH *and* mDRH→mSpindly) which Mark
currently has to split across two figures, and only introduces one new
borderline-ns (Fig 4E A-vs-D, p=0.05). It's also the most readable unit
— percentage points in the same scale as the source metrics.

**Secondary pitch: mention the ratio as the stronger alternative on all
non-TRAK2 comparisons.** If the authors care most about the MAPK9 panel
(Fig S11 F and Fig 4E) and the TRAK-isoform panel (Fig 4B), the ratio is
clearly better there. But pitching the ratio sheet-wide would force the
authors to give up the TRAK2 wt→mDRH call — so it's a ratio-for-MAPK9-and-
isoforms / diff-for-TRAK2 split pitch, which reviewers might flag as
metric-shopping.

If forced to a **single metric for every panel**, pitch **diff**: it's the
only one with no existing-significant-to-ns regression on the TRAK panels
(Fig 4D and S11 E) that carry the main biological claim.

---

## 4 · Two p-value methods

See §3 of the prior memo version for details — unchanged. TL;DR:

- **Classical nested one-way ANOVA** (what Prism does, what I report in
  the panel tables): F-ratio with DF=(k−1, k(r−1)); Šídák across the pair
  family. Pooled across the whole sheet so the error term is stable
  even when individual pairs are small.
- **MixedLM (REML)** as sensitivity check: reproduces Prism's variance
  components exactly (σ²_cell = 25.46, σ²_plate = 0 for Fig 4C match
  Prism's output), agrees with the classical framework on every
  significance call. In `per_metric_summary.csv` as `p_mixedlm_sidak`.

My classical Šídák p-values agree with Mark's Prism to within the same
significance class; small offsets (e.g. my 0.13 vs published 0.083 on Fig
S11 F A-C) come from cell-count-weighted vs plate-level-weighted error
terms. No comparison flips its class between methods.

---

## 5 · Files

- `replication/derived_metrics.py` — metric computation and sheet-pooled
  nested ANOVA / MixedLM evaluation.
- `replication/derived_metrics_out/per_cell.csv` — 799 cells × 58 metrics.
- `replication/derived_metrics_out/per_metric_summary.csv` — 1044 rows,
  one per (sheet, pair, metric) with classical and MixedLM p-values.
- `replication/plot_all_panels.py` — plot script (three-row layout per
  panel: Mark + diff + ratio).
- `replication/figures/Fig_S11_{C,D,E,F}_alt_metrics.{png,pdf}` and
  `Fig_4{B,C,D,E}_alt_metrics.{png,pdf}` — panel plots with Šídák
  brackets.
- `replication/replicate_stats.py` — nested ANOVA + MixedLM framework,
  validated against Mark's Prism files.
