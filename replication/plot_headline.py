"""Plot per-cell distributions for the headline Fig S11 comparison (TRAK1 helix
muts) under Mark's perinuclear metric vs the best-performing derived metric,
so the improvement is visible to the reader.
"""
from __future__ import annotations

import pathlib
import sys

import matplotlib.pyplot as plt
import numpy as np
import polars as pl

REPO = pathlib.Path(__file__).resolve().parent.parent
OUT = REPO / "replication" / "derived_metrics_out"
FIG = REPO / "replication" / "figures"
FIG.mkdir(parents=True, exist_ok=True)


def panel(ax, long: pl.DataFrame, sheet: str, metric: str, conditions: list[str],
          *, y_label: str, title: str):
    sub = long.filter((pl.col("sheet") == sheet) & pl.col(metric).is_not_null())
    colors = plt.get_cmap("tab10").colors
    for i, cond in enumerate(conditions):
        g = sub.filter(pl.col("condition") == cond)
        if g.height == 0:
            continue
        # per-plate means as circles; per-cell as jittered dots
        plates = sorted(g["plate"].unique().to_list())
        x = np.full(g.height, i, dtype=float) + np.random.uniform(-0.18, 0.18, size=g.height)
        ax.scatter(x, g[metric].to_numpy(), alpha=0.25, s=18, color=colors[i],
                   edgecolors="none", label=None)
        for p_idx, plate in enumerate(plates):
            gp = g.filter(pl.col("plate") == plate)
            m = gp[metric].mean()
            ax.scatter([i + (p_idx - (len(plates) - 1) / 2) * 0.08], [m],
                       s=60, color=colors[i], edgecolor="black", linewidths=1.0, zorder=3)
        # condition mean
        mu = g[metric].mean()
        sem = (g.group_by("plate").agg(pl.col(metric).mean())[metric].std() /
               np.sqrt(max(1, g["plate"].n_unique())))
        ax.errorbar([i], [mu], yerr=[sem], fmt="_", color="black", capsize=8,
                    elinewidth=2, markersize=28, zorder=4)

    ax.set_xticks(range(len(conditions)))
    ax.set_xticklabels(conditions, rotation=25, ha="right")
    ax.set_ylabel(y_label)
    ax.set_title(title)


def main():
    long = pl.read_csv(OUT / "per_cell.csv")
    summary = pl.read_csv(OUT / "per_metric_summary.csv")

    def p_for(sheet, pair, metric):
        r = summary.filter((pl.col("sheet") == sheet) & (pl.col("pair") == pair) &
                           (pl.col("metric") == metric))
        return r["p_classical_sidak"].item() if r.height else float("nan")

    sheet = "TRAK1 helix muts"
    conditions = ["T1 wt", "T1 mDRH", "T1 mDRH / dSp"]
    pair = "T1 wt vs T1 mDRH"

    metrics_to_show = [
        ("raw_perinuclear_5um", "% mito within 5 µm of nucleus (raw z-sum)",
         "Mark's Fig S11 D baseline"),
        ("den_peripheral_5um", "% mito within 5 µm of arch (denoised)",
         "Mark's Fig 4C baseline"),
        ("den_peri_over_nuc", "(peripheral 5 µm) ÷ (perinuclear 5 µm)",
         "Proposed polarization ratio"),
        ("den_peri_minus_nuc", "(peripheral − perinuclear) (pp)",
         "Proposed polarization diff."),
    ]
    fig, axes = plt.subplots(1, len(metrics_to_show), figsize=(5 * len(metrics_to_show), 5),
                             sharex=True)
    for ax, (m, ylabel, subtitle) in zip(axes, metrics_to_show):
        p = p_for(sheet, pair, m)
        title = f"{subtitle}\n{m}\n(Šídák p[{pair}] = {p:.4f})"
        panel(ax, long, sheet, m, conditions, y_label=ylabel, title=title)

    fig.suptitle("TRAK1 helix muts — per-cell distributions across 4 plates", y=1.02)
    fig.tight_layout()
    out_path = FIG / "trak1_helix_muts_metrics.png"
    fig.savefig(out_path, dpi=140, bbox_inches="tight")
    fig.savefig(FIG / "trak1_helix_muts_metrics.pdf", bbox_inches="tight")
    print(f"Wrote {out_path}")


if __name__ == "__main__":
    sys.exit(main())
