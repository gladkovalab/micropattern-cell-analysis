"""Exploratory metric (polar-frame analog of the 20–30 µm slab):

  arch_band_pct = sum(wedge_r_45_46um_pct .. wedge_r_54_55um_pct)

The distal-most arch point sits at r ≈ 50.01 µm from the wedge apex
(769.43 pix × 0.065 µm/pix; pitch is constant across the dataset). The
±5 µm window in radial bins becomes bins 45..54 inclusive — a 10-µm slab
centered on the arch.

Same 2x3 layout / nested-ANOVA + Šídák statistics as plot_bin_20_30um.py,
just a different bin window.
"""
from __future__ import annotations

import argparse
import pathlib
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import polars as pl

REPO = pathlib.Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / "replication"))
from plot_metrics import (  # noqa: E402
    SHEET_CONFIG, CONDITION_COLORS, _plate_shades, _test_pair, _format_p,
    load_template_matching, join_with_metadata,
)

# Window: 5 µm either side of the distal-most arch point.
# r_arch_max derived from template_contour[1083:1951] → 769.43 px from apex
# (896, 512) in the cropped frame. At pitch 0.065 µm/pix that is 50.01 µm.
ARCH_R_UM = 50.0
BAND_HALFWIDTH_UM = 5.0
BIN_LO = int(round(ARCH_R_UM - BAND_HALFWIDTH_UM))   # 45
BIN_HI = int(round(ARCH_R_UM + BAND_HALFWIDTH_UM))   # 55  (exclusive)
METRIC_NAME = f"wedge_r_{BIN_LO:02d}_{BIN_HI:02d}um_pct"
METRIC_LABEL = (f"% mito intensity\nin arch ±{int(BAND_HALFWIDTH_UM)} µm "
                f"({BIN_LO}–{BIN_HI} µm)")


def add_metric(df: pl.DataFrame) -> pl.DataFrame:
    cols = [f"wedge_r_{i:02d}_{i+1:02d}um_pct" for i in range(BIN_LO, BIN_HI)]
    cols = [c for c in cols if c in df.columns]
    if not cols:
        raise RuntimeError(
            f"no wedge_r_NN_NN+1um_pct columns in [{BIN_LO}..{BIN_HI}) found")
    return df.with_columns(
        pl.sum_horizontal([pl.col(c) for c in cols]).alias(METRIC_NAME)
    )


def panel(ax, sheet_df: pl.DataFrame, sheet: str, rng):
    cfg = SHEET_CONFIG[sheet]
    conditions = cfg["conditions"]
    pairs = cfg["pairs"]
    family_m = cfg["family_m"]

    color_map = {c: CONDITION_COLORS[i % len(CONDITION_COLORS)]
                 for i, c in enumerate(conditions)}
    plates = sorted(sheet_df["plate"].unique().to_list())
    plate_idx = {p: i for i, p in enumerate(plates)}
    plate_palette = {cond: _plate_shades(color_map[cond], len(plates))
                     for cond in conditions}

    all_vals: list[float] = []
    for j, cond in enumerate(conditions):
        cond_sub = sheet_df.filter(pl.col("condition") == cond)
        plate_means = []
        for plate in plates:
            psub = cond_sub.filter(pl.col("plate") == plate)
            vals = psub[METRIC_NAME].to_numpy().astype(float)
            vals = vals[np.isfinite(vals)]
            if vals.size == 0:
                continue
            shade = plate_palette[cond][plate_idx[plate]]
            xs = j + rng.uniform(-0.18, 0.18, size=vals.size)
            ax.scatter(xs, vals, s=22, color=shade,
                       edgecolor="black", linewidth=0.3, alpha=0.85)
            plate_means.append((j + rng.uniform(-0.08, 0.08),
                                float(np.mean(vals)), shade))
            all_vals.extend(vals.tolist())
        for mx, mv, sh in plate_means:
            ax.scatter(mx, mv, s=85, color=sh, marker="D",
                       edgecolor="black", linewidth=0.8, zorder=6)
        cond_vals = cond_sub[METRIC_NAME].to_numpy().astype(float)
        cond_vals = cond_vals[np.isfinite(cond_vals)]
        if cond_vals.size > 0:
            ax.hlines(float(np.mean(cond_vals)), j - 0.32, j + 0.32,
                      color="black", linewidth=2.0, zorder=7)
    ax.set_xticks(range(len(conditions)))
    ax.set_xticklabels(conditions, rotation=20, ha="right", fontsize=16)
    n_total = sheet_df.height
    ax.set_title(f"{sheet}  (n={n_total})", fontsize=10)
    ax.grid(axis="y", alpha=0.3)

    if not all_vals:
        return
    data_max = max(all_vals)
    data_min = min(all_vals)
    span = max(data_max - data_min, 1e-6)
    bracket_base = data_max + span * 0.05
    bracket_step = span * 0.10
    ax.set_ylim(data_min - span * 0.05,
                bracket_base + bracket_step * (len(pairs) + 0.5))

    for k, (a, b) in enumerate(pairs):
        p = _test_pair(sheet_df, METRIC_NAME, (a, b), family_m)
        if not np.isfinite(p):
            continue
        xa = conditions.index(a)
        xb = conditions.index(b)
        y = bracket_base + bracket_step * (k + 0.3)
        tick = bracket_step * 0.18
        ax.plot([xa, xa, xb, xb], [y - tick, y, y, y - tick],
                color="black", linewidth=0.8)
        ax.text((xa + xb) / 2, y + bracket_step * 0.05, _format_p(p),
                ha="center", va="bottom", fontsize=14)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--template-matching",
                    default="replication/wedge_r_ks_out_all_denoised/by_well")
    ap.add_argument("--comparisons-xlsx",
                    default=str(REPO / "config/Comparisons_table_v3.xlsx"))
    ap.add_argument("--out",
                    default="replication/figures_wedge_r_ks/arch_band_polar_all.png")
    ap.add_argument("--title", default=None,
                    help="Override the figure suptitle.")
    args = ap.parse_args()

    df = load_template_matching(pathlib.Path(args.template_matching))
    df = join_with_metadata(df, pathlib.Path(args.comparisons_xlsx))
    df = add_metric(df)

    sheet_order = [
        "TRAK isoform (mito)",
        "TRAK isoform (peroxisome)",
        "TRAK isoform (60mer)",
        "TRAK1 helix muts",
        "TRAK2 helix muts",
        "MAPK9 siRNA",
    ]

    print(f"=== per-sheet stats: {METRIC_LABEL.replace(chr(10), ' ')} ===")
    for sheet in sheet_order:
        cfg = SHEET_CONFIG.get(sheet)
        if cfg is None:
            continue
        sub = df.filter(pl.col("sheet") == sheet)
        if sub.height == 0:
            continue
        print(f"\n  {sheet}  (n={sub.height})")
        for cond in cfg["conditions"]:
            v = sub.filter(pl.col("condition") == cond)[METRIC_NAME].to_numpy().astype(float)
            v = v[np.isfinite(v)]
            if v.size == 0:
                continue
            sem = v.std(ddof=1) / max(np.sqrt(v.size), 1) if v.size > 1 else float("nan")
            print(f"    {cond:25s}: n={v.size:3d}  mean={v.mean():7.3f}%  sem={sem:6.3f}")
        for pair in cfg["pairs"]:
            p = _test_pair(sub, METRIC_NAME, pair, cfg["family_m"])
            print(f"    {pair[0]:>18s}  vs  {pair[1]:<18s}  →  {_format_p(p)}")

    rng = np.random.default_rng(0)
    fig, axes = plt.subplots(2, 3, figsize=(15, 9))
    for ax, sheet in zip(axes.flat, sheet_order):
        sub = df.filter(pl.col("sheet") == sheet)
        if sub.height == 0:
            ax.text(0.5, 0.5, f"no rows: {sheet}", ha="center", va="center",
                    transform=ax.transAxes)
            ax.set_xticks([])
            continue
        panel(ax, sub, sheet, rng)
        ax.set_ylabel(METRIC_LABEL, fontsize=9)

    suptitle = args.title or (
        f"Wedge-r {BIN_LO}–{BIN_HI} µm slab (±{int(BAND_HALFWIDTH_UM)} µm of distal arch "
        f"@ r={ARCH_R_UM:.0f} µm) — fraction of wedge intensity\n"
        "nested ANOVA + Šídák pairwise (Welch fallback for single-plate)"
    )
    fig.suptitle(suptitle, fontsize=14)
    plt.tight_layout(rect=[0, 0, 1, 0.95])
    out = pathlib.Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=140)
    plt.close(fig)
    print(f"\nwrote {out}")


if __name__ == "__main__":
    main()
