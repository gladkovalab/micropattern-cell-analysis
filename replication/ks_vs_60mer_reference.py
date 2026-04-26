"""Re-compute wedge-r KS using the 60mer no-TRAK condition as the empirical
"true uniform" reference, instead of the analytical area-uniform sector.

The 60mer (synthetic particle) under no-TRAK is the cleanest available
proxy for a passive cytoplasmic fill on this micropattern, so its mean
wedge-r CDF is a defensible empirical reference for "no biased transport".

Per cell, the new metric is:
    KS_vs_60merNoTRAK = max_i |cdf_cell[i] - mean_cdf_60merNoTRAK[i]|

Computed for both projections (zsum, MaxIP) on every cell in the merged
494+47-cell dataset. Output: per-sheet scalar figures with Šídák brackets
matching the canonical plot family.
"""
from __future__ import annotations
import pathlib, re, sys
import numpy as np
import polars as pl
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

REPO = pathlib.Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / "replication"))

from plot_final import (  # noqa: E402
    SHEET_CONFIG, EVAL_FAMILY_M, CONDITION_COLORS,
    extract_profile, slug, _format_p, load,
)
from evaluate_final import test_pair  # noqa: E402

OUT_DIR = REPO / "replication" / "overnight_final_out" / "figures_ks_vs_60mer"
EVAL_OUT = REPO / "replication" / "overnight_final_out" / "evaluation_summary_ks_vs_60mer.csv"


def per_cell_cdf(profile: np.ndarray) -> np.ndarray:
    """Convert a (n_cells, 60) per-bin % matrix into per-cell CDFs."""
    row_sums = np.nansum(profile, axis=1, keepdims=True)
    out = np.full_like(profile, np.nan, dtype=float)
    ok = (row_sums.flatten() > 0)
    out[ok] = np.nancumsum(profile[ok] / row_sums[ok], axis=1)
    return out


def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    df = load()
    print(f"Loaded {df.height} cells across "
          f"{sorted(df['sheet'].unique().to_list())}")

    # --- 1. Build the empirical reference CDFs (zsum, maxip) from 60mer no-TRAK
    ref_60 = df.filter((pl.col("sheet") == "TRAK isoform (60mer)") &
                       (pl.col("condition") == "no TRAK"))
    print(f"60mer no-TRAK reference: n = {ref_60.height} cells")

    ref_cdf = {}
    for proj in ("zsum", "maxip"):
        centers, prof = extract_profile(ref_60, proj, "wedge_r")
        cdf = per_cell_cdf(prof)
        ref_cdf[proj] = np.nanmean(cdf, axis=0)
    print(f"Reference CDF anchors (zsum):  bin5={ref_cdf['zsum'][5]:.3f}  "
          f"bin30={ref_cdf['zsum'][30]:.3f}  bin55={ref_cdf['zsum'][55]:.3f}")
    print(f"Reference CDF anchors (maxip): bin5={ref_cdf['maxip'][5]:.3f}  "
          f"bin30={ref_cdf['maxip'][30]:.3f}  bin55={ref_cdf['maxip'][55]:.3f}")

    # Compare to analytical sector reference at same anchors
    n_bins = ref_cdf["zsum"].size
    cdf_sector = (np.arange(1, n_bins + 1) ** 2) / (n_bins ** 2)
    print(f"Analytical (area-uniform):     bin5={cdf_sector[5]:.3f}  "
          f"bin30={cdf_sector[30]:.3f}  bin55={cdf_sector[55]:.3f}")

    # --- 2. Compute the new metric per cell, per projection
    new_cols = []
    for proj in ("zsum", "maxip"):
        centers, prof = extract_profile(df, proj, "wedge_r")
        cdf = per_cell_cdf(prof)
        ks = np.nanmax(np.abs(cdf - ref_cdf[proj]), axis=1)
        col_name = f"{proj}_wedge_r_ks_vs_60merNoTRAK"
        new_cols.append(pl.Series(col_name, ks))
    df = df.with_columns(new_cols)

    # --- 3. Evaluate Šídák stats per sheet for BOTH the new and old metrics
    rows = []
    metric_set = ("zsum_wedge_r_ks_vs_60merNoTRAK",
                  "maxip_wedge_r_ks_vs_60merNoTRAK",
                  "zsum_wedge_r_ks_vs_uniform",
                  "maxip_wedge_r_ks_vs_uniform")
    for sheet, cfg in SHEET_CONFIG.items():
        sheet_df = df.filter(pl.col("sheet") == sheet)
        if sheet_df.height == 0:
            continue
        for pair in cfg["pairs"]:
            for m in metric_set:
                r = test_pair(sheet_df, m, pair, cfg["family_m"])
                if r is None:
                    continue
                rows.append({"sheet": sheet,
                             "pair": f"{pair[0]} vs {pair[1]}",
                             "metric": m, **r})
    eval_df = pl.from_dicts(rows)
    eval_df.write_csv(EVAL_OUT)
    print(f"\nWrote {EVAL_OUT} ({eval_df.height} rows)")

    # --- 4. Plot per-sheet 1×2 scalar figures with Šídák brackets
    SC = 1.5  # text scale (matches canonical scalar figures)
    rng = np.random.default_rng(0)
    for sheet, cfg in SHEET_CONFIG.items():
        sheet_df = df.filter(pl.col("sheet") == sheet)
        if sheet_df.height == 0:
            continue
        conditions = cfg["conditions"]
        pairs = cfg["pairs"]
        family_m = cfg["family_m"]
        cond_idx = {c: i for i, c in enumerate(conditions)}
        plates = sorted(sheet_df["plate"].unique().to_list())
        plate_markers = dict(zip(plates,
                                 ["o","s","D","^","v","P","X","<",">","*"]))
        color_map = {c: CONDITION_COLORS[i % len(CONDITION_COLORS)]
                     for i, c in enumerate(conditions)}

        # 2 rows × 2 cols: top = new metric (vs 60mer-noTRAK),
        # bottom = old metric (vs analytical area-uniform)
        fig, axes = plt.subplots(2, 2, figsize=(11, 13), sharex=True)
        row_specs = [
            (axes[0], "_ks_vs_60merNoTRAK",  "wedge-r KS vs 60mer-noTRAK"),
            (axes[1], "_ks_vs_uniform",      "wedge-r KS vs area-uniform (original)"),
        ]
        for row, suffix, ylab_root in row_specs:
            for ax, proj in zip(row, ("zsum", "maxip")):
                metric = f"{proj}_wedge_r{suffix}"
                for j, cond in enumerate(conditions):
                    sub = sheet_df.filter(pl.col("condition") == cond)
                    for plate in plates:
                        psub = sub.filter(pl.col("plate") == plate)
                        if psub.height == 0:
                            continue
                        vals = psub[metric].to_numpy()
                        xs = j + rng.uniform(-0.18, 0.18, size=len(vals))
                        ax.scatter(xs, vals, marker=plate_markers[plate],
                                   s=42, color=color_map[cond],
                                   edgecolor="black", linewidth=0.4,
                                   alpha=0.85)
                    mn = sub[metric].mean()
                    if mn is not None:
                        ax.hlines(mn, j - 0.3, j + 0.3, color="black",
                                  linewidth=2.5, zorder=5)
                ax.set_xticks(range(len(conditions)))
                ax.set_xticklabels(conditions, rotation=15, ha="right",
                                   fontsize=8 * SC)
                ax.set_ylabel(f"{ylab_root} ({proj.upper()})",
                              fontsize=10 * SC)
                ax.tick_params(axis="y", labelsize=9 * SC)
                ax.grid(axis="y", alpha=0.3)

                # Fix y-axis to [0, 1] and stack brackets within that range
                ax.set_ylim(0.0, 1.0)
                data_max = float(
                    sheet_df[metric].drop_nulls().drop_nans().max() or 0.0)
                n_pairs_eff = max(len(pairs), 1)
                bracket_base = max(data_max + 0.03, 0.78)
                bracket_step = min(0.05, (0.97 - bracket_base) / n_pairs_eff)
                tick = bracket_step * 0.18
                for k, (a, b) in enumerate(pairs):
                    if a not in cond_idx or b not in cond_idx:
                        continue
                    x1, x2 = sorted((cond_idx[a], cond_idx[b]))
                    row_e = eval_df.filter(
                        (pl.col("sheet") == sheet) &
                        (pl.col("metric") == metric) &
                        (pl.col("pair") == f"{a} vs {b}"))
                    if row_e.height == 0:
                        continue
                    p = row_e["p"][0]
                    if p is None or not np.isfinite(p):
                        continue
                    txt, _ = _format_p(float(p))
                    y_bar = bracket_base + bracket_step * (k + 0.5)
                    ax.plot([x1, x1, x2, x2],
                            [y_bar - tick, y_bar, y_bar, y_bar - tick],
                            lw=0.9, color="black")
                    ax.text((x1 + x2) / 2, y_bar + bracket_step * 0.05, txt,
                            ha="center", va="bottom", fontsize=7.5 * SC)
                ax.set_ylim(0.0, 1.0)  # reaffirm in case anything nudged it

        fig.suptitle(
            f"{sheet} (n={sheet_df.height})  ·  "
            f"wedge-r KS — new (top) vs original (bottom)  ·  "
            f"Šídák m={family_m}",
            fontsize=11 * SC)
        plt.tight_layout(rect=[0, 0.02, 1, 0.96])
        out = OUT_DIR / f"{slug(sheet)}_ks_vs_60mer.png"
        fig.savefig(out, dpi=130)
        plt.close(fig)
        print(f"  wrote {out.name}")


if __name__ == "__main__":
    main()
