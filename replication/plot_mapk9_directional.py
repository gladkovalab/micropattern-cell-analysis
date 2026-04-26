"""MAPK9 scalar figure with directional metrics.

Replaces the {Gini, σ} scalars on the standard `MAPK9_siRNA_scalars.png` with
metrics that pick up the Ars phenotype (peripheral depletion → perinuclear
pile-up): peripheral_5um_pct, perinuclear_5um_pct, wedge_r_ks_vs_uniform.
Both zsum and MaxIP projections, mirroring the layout of plot_final.plot_scalars.
"""
from __future__ import annotations

import pathlib

import numpy as np
import polars as pl
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from plot_final import load, OUT_DIR, CONDITION_COLORS

SHEET = "MAPK9 siRNA"
CONDITIONS = ["ctrl ctrl", "ctrl Ars", "MAPK9 ctrl", "MAPK9 Ars"]

METRICS = [
    ("zsum_peripheral_5um_pct",       "peripheral 5 µm % (zsum)"),
    ("maxip_peripheral_5um_pct",      "peripheral 5 µm % (MaxIP)"),
    ("zsum_perinuclear_5um_pct",      "perinuclear 5 µm % (zsum)"),
    ("maxip_perinuclear_5um_pct",     "perinuclear 5 µm % (MaxIP)"),
    ("zsum_wedge_r_ks_vs_uniform",    "wedge-r KS vs uniform (zsum)"),
    ("maxip_wedge_r_ks_vs_uniform",   "wedge-r KS vs uniform (MaxIP)"),
]


def main():
    df = load().filter(pl.col("sheet") == SHEET)
    plates = sorted(df["plate"].unique().to_list())
    plate_markers = dict(zip(plates,
                             ["o", "s", "D", "^", "v", "P", "X", "<", ">", "*"]))
    color_map = {c: CONDITION_COLORS[i % len(CONDITION_COLORS)]
                 for i, c in enumerate(CONDITIONS)}

    fig, axes = plt.subplots(3, 2, figsize=(11, 12))
    rng = np.random.default_rng(0)
    for ax, (m, label) in zip(axes.flat, METRICS):
        if m not in df.columns:
            ax.text(0.5, 0.5, f"no column {m}",
                    transform=ax.transAxes, ha="center")
            continue
        for j, cond in enumerate(CONDITIONS):
            sub = df.filter(pl.col("condition") == cond)
            for plate in plates:
                psub = sub.filter(pl.col("plate") == plate)
                if psub.height == 0:
                    continue
                vals = psub[m].to_numpy()
                xs = j + rng.uniform(-0.18, 0.18, size=len(vals))
                ax.scatter(xs, vals, marker=plate_markers[plate], s=42,
                           color=color_map[cond], edgecolor="black",
                           linewidth=0.4, alpha=0.85,
                           label=plate if ax is axes.flat[0] else None)
            mn = sub[m].mean()
            if mn is not None:
                ax.hlines(mn, j - 0.3, j + 0.3, color="black",
                          linewidth=2.5, zorder=5)
        ax.set_xticks(range(len(CONDITIONS)))
        ax.set_xticklabels(CONDITIONS, rotation=15, ha="right", fontsize=8)
        ax.set_ylabel(label)
        ax.grid(axis="y", alpha=0.3)
    axes.flat[0].legend(loc="best", fontsize=6, title="plate",
                        title_fontsize=7, ncol=2)
    fig.suptitle(f"{SHEET} · directional scalars (peripheral / perinuclear / wedge-r KS)",
                 fontsize=11)
    plt.tight_layout(rect=[0, 0.02, 1, 0.97])
    out = OUT_DIR / "MAPK9_siRNA_scalars_directional.png"
    fig.savefig(out, dpi=130)
    plt.close(fig)
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
