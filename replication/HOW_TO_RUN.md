# How to run the wedge-r KS analysis

End-to-end reproduction for the analysis archived to Zenodo and cited by
the manuscript. A fresh operator with a clean checkout of this branch
should be able to follow this top-to-bottom and produce the same per-cell
CSVs, per-sheet figures, and source-data XLSXs.

If anything fails, fix the underlying issue rather than the script —
the pipeline output is the artifact, so silent fallbacks would
contaminate it.

---

## 1. What this analysis is

For each cell on a fixed-geometry micropattern, the patched pipeline:

1. Reads the ND2, applies Mark's per-channel projection (MaxIP on the
   488 mitochondria channel; Z-sum on the 405 nuclear-segmentation
   channel — see `WEDGE_R_KS.md` for the rationale).
2. Runs Mark's existing template-matching + background subtraction.
3. Computes an intensity-weighted radial profile within a fixed
   polygonal wedge over the micropattern arch (apex at the bottom
   extreme of the pattern, two tangent points at the top corners).
4. Emits per-cell scalars: `wedge_r_ks_vs_uniform`,
   `wedge_r_ks_vs_60merNoTRAK`, the two slab metrics
   `wedge_r_centrosomal_18_33um_pct` and `wedge_r_peripheral_41_56um_pct`,
   plus the 60-bin profile `wedge_r_NN_NN+1um_pct` for `NN ∈ {00..59}`.
5. Preserves Mark's existing output columns (`peripheral_5um_*`,
   `perinuclear_5um_*`, etc.), now derived from the MaxIP projection.

`WEDGE_R_KS.md` is the reviewer-facing description of the metric.
`HANDOFF_v4.md` is the session log for the slab-metric extension and
the rebuttal analysis — not required reading for reproduction.

---

## 2. Prerequisites

1. **Repo on this branch.** `wpg/wedge-r-ks`. The patched pipeline is in
   `template_matching_bulk.py`; everything else added for this analysis
   lives under `replication/`.
2. **Pixi environment.** All commands below are run via `pixi run`,
   which transparently activates the locked environment. No manual
   `conda activate` step.
3. **Git LFS.** Used for the committed `wedge_illustration.png`. If you
   cloned without LFS, run `pixi run git lfs pull` once.
4. **Data root.** The raw ND2s for the six comparison sheets live under
   `mark_data/patterned_data/` (in this repo, `mark_data` is a symlink
   into the lab fileshare). If your data lives elsewhere, point the
   driver at it via `--data-root /abs/path/to/patterned_data`. **Do not
   modify anything under that root** — the pipeline is read-only on the
   data side.
5. **Authoritative metadata.** `config/Comparisons_table_v3.xlsx`
   (committed to upstream/main) maps each (sheet, condition) pair to
   `(plate, well)` pairs. The driver reads it directly; do not maintain
   a parallel metadata CSV.
6. **macOS only:** if running locally, prefix long commands with
   `caffeinate -dimsu` so the machine doesn't sleep mid-run (see
   §4). On Linux this is unnecessary.

Quick environment check:

```sh
pixi run python -c "import template_matching_bulk as tmb; \
print('apex:', tmb.WEDGE_APEX, 'L:', tmb.WEDGE_LEFT, 'R:', tmb.WEDGE_RIGHT); \
print('centrosomal bins:', tmb.WEDGE_CENTROSOMAL_BINS, \
      'peripheral bins:', tmb.WEDGE_PERIPHERAL_BINS)"
```

Expected output:

```
apex: (896, 512) L: (373, 281) R: (374, 742)
centrosomal bins: (18, 33) peripheral bins: (41, 56)
```

If any constant differs, you are not on the right branch.

---

## 3. The six comparison sheets

The full analysis processes these six sheets from
`config/Comparisons_table_v3.xlsx` (`--variant denoised` in all cases —
denoised ND2s live under each well's `denoised/` subdirectory):

| `--sheet` argument              | description                          |
| ---                             | ---                                  |
| `TRAK isoform (mito)`           | TRAK isoforms vs the mito construct  |
| `TRAK isoform (peroxisome)`     | same isoforms, peroxisome target     |
| `TRAK isoform (60mer)`          | 60mer scaffold series                 |
| `TRAK1 helix muts`              | TRAK1 helix mutation panel           |
| `TRAK2 helix muts`              | TRAK2 helix mutation panel           |
| `MAPK9 siRNA`                   | MAPK9 knockdown                       |

If you want to smoke-test before the overnight run, start with
`TRAK isoform (mito)` — smallest sheet, ~70 cells, finishes in well
under an hour off SMB.

---

## 4. Run the pipeline

The output root must be a fresh, empty directory. Existing per-cell JSON
checkpoints under `{out_root}/by_well/.../cells/` are honored — that's
how the run is resumable across SMB drops — but it also means a stale
cache will silently skip cells. **For an archive-quality run, start
from an empty `--out-root`.**

```sh
OUTDIR=replication/wedge_r_ks_out_all_denoised
mkdir -p "$OUTDIR"

# If $OUTDIR already exists from a prior run, wipe the cache first:
# rm -rf "$OUTDIR/by_well"

nohup caffeinate -dimsu bash -c '
  set -e
  OUTDIR=replication/wedge_r_ks_out_all_denoised
  for sheet in "TRAK isoform (mito)" \
               "TRAK isoform (peroxisome)" \
               "TRAK isoform (60mer)" \
               "TRAK1 helix muts" \
               "TRAK2 helix muts" \
               "MAPK9 siRNA"; do
    echo "=== $(date) === starting: $sheet ==="
    pixi run python replication/run_pipeline_paths.py \
      --sheet "$sheet" --variant denoised --out-root "$OUTDIR"
    echo "=== $(date) === completed: $sheet ==="
  done
  echo "=== $(date) === ALL SHEETS DONE ==="
' > "$OUTDIR/stdout_overnight.log" 2>&1 &

echo "pid=$!"
```

Notes:

- Wall-clock: ~24 h over SMB for all six sheets on a 2024 MacBook Pro.
  Locally-mounted data is closer to 2–4 h. Most of the wall time is ND2
  decode.
- The driver discovers cells from
  `config/Comparisons_table_v3.xlsx`, falls back to the well's denoised
  ND2s on disk, and skips wells that aren't actually present.
- Per-well CSVs land at
  `$OUTDIR/by_well/{plate}/{well_dir}/template_matching.csv` and a
  flat `$OUTDIR/combined.csv` is rewritten at the end of each sheet.
- Tail the log:
  `tail -f replication/wedge_r_ks_out_all_denoised/stdout_overnight.log`.

If the pipeline aborts mid-sheet, just re-run the same command — cached
cells are skipped via their `cells/*.json` checkpoints.

---

## 5. Generate the figures

Once the pipeline run is complete (look for `ALL SHEETS DONE` in the
log), generate one figure per sheet:

```sh
mkdir -p replication/figures_wedge_r_ks

for label in "trak_isoform_mito:TRAK isoform (mito)" \
             "trak_isoform_peroxisome:TRAK isoform (peroxisome)" \
             "trak_isoform_60mer:TRAK isoform (60mer)" \
             "trak1_helix_muts:TRAK1 helix muts" \
             "trak2_helix_muts:TRAK2 helix muts" \
             "mapk9_sirna:MAPK9 siRNA"; do
  slug=${label%%:*}
  sheet=${label#*:}
  pixi run python replication/plot_metrics.py \
    --template-matching replication/wedge_r_ks_out_all_denoised/by_well \
    --sheet "$sheet" \
    --out "replication/figures_wedge_r_ks/${slug}.png"
done
```

Each PNG is a 3-row figure: wedge-r profile + CDF on top, six scalar
strip plots below (peripheral 5 µm, perinuclear 5 µm, centrosomal slab,
peripheral slab, KS vs area-uniform, KS vs 60mer no-TRAK). Per-plate
shading and nested-ANOVA + Šídák brackets are applied automatically.

To regenerate the canonical wedge-on-cell illustration:

```sh
pixi run python replication/plot_wedge_illustration.py
# writes replication/wedge_illustration.png
```

(Requires read access to the source ND2; use
`plot_wedge_illustration_offline.py` if running from saved projections
only.)

---

## 6. Validation / spot checks

A single-cell sanity check before trusting the full output:

```sh
pixi run python - <<'PY'
import polars as pl
df = pl.read_csv("replication/wedge_r_ks_out_all_denoised/combined.csv")
hit = df.filter(
    pl.col("path").str.contains("/plate_3_") &
    pl.col("path").str.contains("/B03_") &
    pl.col("path").str.contains("cell1")
)
print(hit.select([
    "wedge_r_ks_vs_uniform",
    "wedge_r_centrosomal_18_33um_pct",
    "wedge_r_peripheral_41_56um_pct",
]))
PY
```

Expected (matches the values reported in `WEDGE_R_KS.md` and
`HANDOFF_v4.md` within rounding):

- `wedge_r_ks_vs_uniform`            ≈ 0.397 (this branch reports 0.3968)
- `wedge_r_centrosomal_18_33um_pct`  ≈ 65.6  (TRAK2 mito plate-3 cell)
- `wedge_r_peripheral_41_56um_pct`   ≈ 1.9

For the headline TRAK1 vs TRAK2 statistic on the mito sheet, the
expected `wedge_r_ks_vs_uniform` nested-ANOVA + Šídák pairwise p-value
is **0.0013**. The same p-value is printed on the `trak_isoform_mito.png`
strip plot.

Any deviation > a few percent on these spot-check values means the
pipeline is reading a different input set or running with different
constants — investigate before treating the run as final.

---

## 7. What gets archived

Inputs needed to rerun:

- `template_matching_bulk.py` (the patched pipeline)
- `replication/run_pipeline_paths.py` (driver)
- `replication/plot_metrics.py` (figure generator)
- `replication/plot_wedge_illustration.py` (geometry illustration)
- `replication/WEDGE_R_KS.md` (metric definition)
- `replication/HOW_TO_RUN.md` (this file)
- `config/Comparisons_table_v3.xlsx`
- `coordinate_overrides.csv`

Outputs to archive:

- `replication/wedge_r_ks_out_all_denoised/combined.csv` — per-cell
  results, all six sheets concatenated.
- `replication/wedge_r_ks_out_all_denoised/by_well/**/template_matching.csv`
  — per-well CSVs mirroring Mark's existing layout.
- `replication/figures_wedge_r_ks/{slug}.png` — one figure per sheet.
- `replication/wedge_illustration.png` — geometry illustration.

The per-cell JSON checkpoints under `cells/` are intermediate cache —
keep them for resumability during the run, but they don't need to land
in Zenodo (`combined.csv` is the authoritative per-cell artifact).

---

## 8. Common pitfalls

- **Empty `--out-root` is required for archive runs.** If you pass an
  out-root that already has a `by_well/` tree, the driver will reuse
  cached JSONs — fine for resuming an interrupted run, but it means
  code changes to the pipeline are silently skipped. Always wipe
  `{out_root}/by_well` between code revisions.
- **`mark_data` is a symlink.** A `git clone` will leave a dangling
  symlink. Either recreate it (`ln -s /path/to/lab/patterned_data
  mark_data`) or pass `--data-root` to `run_pipeline_paths.py`.
- **SMB stalls.** The driver writes one JSON per cell, so an SMB drop
  loses at most the in-progress cell. Just rerun the same command.
- **`caffeinate` is macOS-only.** On Linux, drop the `caffeinate -dimsu`
  prefix; everything else is the same.
- **Don't write into `mark_data/`.** All outputs go under `replication/`.
- **`config/Comparisons_table_v3.xlsx` is authoritative.** Do not
  introduce a parallel `metadata.csv` — the driver reads the xlsx
  directly via `fastexcel`.
