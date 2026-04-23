"""Generate side-by-side panel plots for Fig S11 and Fig 4, comparing Mark's
current metric against the proposed `peri − nuc` polarization difference,
with Šídák-corrected p-value brackets for the pairs Mark reports in his
Prism files.

Each figure corresponds to one Fig panel (S11 C/D/E/F or 4 C/D/E). Rows
within a figure are metrics (Mark's vs proposed).
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

# Panels: (fig_label, sheet_name, condition order, pairs to annotate, mark_metric, mark_ylabel, mark_title)
# Pair choices match what Mark actually tested in the corresponding Prism file:
#   - helix-muts and peroxisome panels: adjacent pairs on the sheet
#   - MAPK9 perinuclear (Fig S11 F): A-C (knockdown effect ctrl), B-D (knockdown effect Ars)
#   - MAPK9 peripheral  (Fig 4E):    A-B (Ars in wt), A-C (knockdown), A-D (combined)
PANELS = [
    # --- Fig S11 (perinuclear baseline) ---
    ("Fig S11 C  peroxisome  TRAK isoforms",
     "TRAK isoform (peroxisome)",
     ["no TRAK", "TRAK1", "TRAK2"],
     [("no TRAK", "TRAK1"), ("no TRAK", "TRAK2"), ("TRAK1", "TRAK2")],
     "raw_perinuclear_5um",
     "% mito within 5 µm of nucleus",
     "Mark's current metric — perinuclear 5 µm (raw z-sum)"),
    ("Fig S11 D  TRAK1 helix mutants",
     "TRAK1 helix muts",
     ["T1 wt", "T1 mDRH", "T1 mDRH / dSp"],
     [("T1 wt", "T1 mDRH"), ("T1 mDRH", "T1 mDRH / dSp")],
     "raw_perinuclear_5um",
     "% mito within 5 µm of nucleus",
     "Mark's current metric — perinuclear 5 µm (raw z-sum)"),
    ("Fig S11 E  TRAK2 helix mutants",
     "TRAK2 helix muts",
     ["TRAK2", "TRAK2 mDRH", "TRAK2 mDRH mSpindly"],
     [("TRAK2", "TRAK2 mDRH"), ("TRAK2 mDRH", "TRAK2 mDRH mSpindly")],
     "raw_perinuclear_5um",
     "% mito within 5 µm of nucleus",
     "Mark's current metric — perinuclear 5 µm (raw z-sum)"),
    ("Fig S11 F  MAPK9/JNK2 siRNA + arsenite",
     "MAPK9 siRNA",
     ["ctrl ctrl", "ctrl Ars", "MAPK9 ctrl", "MAPK9 Ars"],
     # Published Fig S11 F pairs (A-C, B-D, C-D in Mark's Prism column order):
     #   ctrl ctrl vs ctrl Ars, MAPK9 ctrl vs MAPK9 Ars, ctrl Ars vs MAPK9 Ars
     # m=3 Šídák family
     [("ctrl ctrl", "ctrl Ars"), ("MAPK9 ctrl", "MAPK9 Ars"),
      ("ctrl Ars", "MAPK9 Ars")],
     "raw_perinuclear_5um",
     "% mito within 5 µm of nucleus",
     "Mark's current metric — perinuclear 5 µm (raw z-sum)"),

    # --- Fig 4 (peripheral-denoised baseline) ---
    ("Fig 4B  TRAK isoforms (mito)",
     "TRAK isoform (mito)",
     ["no TRAK", "TRAK1", "TRAK2"],
     # Prism m=3: A-B, A-C, B-C (no TRAK vs TRAK1, no TRAK vs TRAK2, TRAK1 vs TRAK2)
     [("no TRAK", "TRAK1"), ("no TRAK", "TRAK2"), ("TRAK1", "TRAK2")],
     "den_peripheral_5um",
     "% mito within 5 µm of arch",
     "Mark's current metric — peripheral 5 µm (denoised)"),
    ("Fig 4C  TRAK1 helix mutants",
     "TRAK1 helix muts",
     ["T1 wt", "T1 mDRH", "T1 mDRH / dSp"],
     [("T1 wt", "T1 mDRH"), ("T1 mDRH", "T1 mDRH / dSp")],
     "den_peripheral_5um",
     "% mito within 5 µm of arch",
     "Mark's current metric — peripheral 5 µm (denoised)"),
    ("Fig 4D  TRAK2 helix mutants",
     "TRAK2 helix muts",
     ["TRAK2", "TRAK2 mDRH", "TRAK2 mDRH mSpindly"],
     [("TRAK2", "TRAK2 mDRH"), ("TRAK2 mDRH", "TRAK2 mDRH mSpindly")],
     "den_peripheral_5um",
     "% mito within 5 µm of arch",
     "Mark's current metric — peripheral 5 µm (denoised)"),
    ("Fig 4E  MAPK9/JNK2 siRNA + arsenite",
     "MAPK9 siRNA",
     ["ctrl ctrl", "ctrl Ars", "MAPK9 ctrl", "MAPK9 Ars"],
     # Mark's Prism (A-B, A-C, A-D in his order ctrl/MAPK9/ctrl-Ars/MAPK9-Ars)
     # = each perturbation vs the double control
     [("ctrl ctrl", "MAPK9 ctrl"), ("ctrl ctrl", "ctrl Ars"), ("ctrl ctrl", "MAPK9 Ars")],
     "den_peripheral_5um",
     "% mito within 5 µm of arch",
     "Mark's current metric — peripheral 5 µm (denoised)"),
]

PROPOSED_ROWS = [
    ("den_peri_minus_nuc",
     "peripheral − perinuclear  (pp)",
     "Proposed — peripheral 5 µm − perinuclear 5 µm (denoised)"),
    ("den_peri_over_nuc",
     "peripheral ÷ perinuclear",
     "Alternative — peripheral 5 µm ÷ perinuclear 5 µm (denoised)"),
]


def p_fmt(p):
    if p < 0.0001: return "< 0.0001"
    if p < 0.001: return f"{p:.4f}"
    if p < 0.01: return f"{p:.3f}"
    return f"{p:.2f}"


def sig(p):
    if p < 0.001: return "***"
    if p < 0.01: return "**"
    if p < 0.05: return "*"
    return "ns"


def add_bracket(ax, x0, x1, y, label, *, pad=0.02):
    h = 0.03 * (ax.get_ylim()[1] - ax.get_ylim()[0])
    ax.plot([x0, x0, x1, x1], [y, y + h, y + h, y], color="black", lw=1.0)
    ax.text((x0 + x1) / 2, y + h + pad * (ax.get_ylim()[1] - ax.get_ylim()[0]),
            label, ha="center", va="bottom", fontsize=9)


def sidak(p_raw: float, m: int) -> float:
    return 1 - (1 - p_raw) ** m


def draw_subplot(ax, long: pl.DataFrame, summary: pl.DataFrame, sheet: str,
                 metric: str, conditions: list[str], pairs: list[tuple[str, str]],
                 ylabel: str, title: str):
    sub = long.filter((pl.col("sheet") == sheet) & pl.col(metric).is_not_null() &
                      pl.col(metric).is_not_nan())
    colors = plt.get_cmap("tab10").colors
    for i, cond in enumerate(conditions):
        g = sub.filter(pl.col("condition") == cond)
        if g.height == 0:
            continue
        plates = sorted(g["plate"].unique().to_list())
        x = np.full(g.height, i, dtype=float) + np.random.uniform(-0.15, 0.15, size=g.height)
        ax.scatter(x, g[metric].to_numpy(), alpha=0.22, s=14,
                   color=colors[i % len(colors)], edgecolors="none")
        for p_idx, plate in enumerate(plates):
            gp = g.filter(pl.col("plate") == plate)
            m = gp[metric].mean()
            ax.scatter([i + (p_idx - (len(plates) - 1) / 2) * 0.08], [m],
                       s=60, color=colors[i % len(colors)], edgecolor="black",
                       linewidths=1.0, zorder=3)
        mu = g[metric].mean()
        sem = (g.group_by("plate").agg(pl.col(metric).mean())[metric].std() /
               np.sqrt(max(1, g["plate"].n_unique())))
        ax.errorbar([i], [mu], yerr=[sem], fmt="_", color="black", capsize=8,
                    elinewidth=2, markersize=28, zorder=4)

    ys_all = sub[metric].to_numpy()
    y_top = np.nanmax(ys_all) if ys_all.size else 1.0
    y_bot = np.nanmin(ys_all) if ys_all.size else 0.0
    rng = y_top - y_bot
    ax.set_ylim(y_bot - 0.05 * rng, y_top + 0.35 * rng * (1 + 0.3 * len(pairs)))
    y_base = y_top + 0.05 * rng

    m = len(pairs)
    for k, (a, b) in enumerate(pairs):
        if a not in conditions or b not in conditions:
            continue
        ia, ib = conditions.index(a), conditions.index(b)
        pair_label = f"{a} vs {b}"
        r = summary.filter((pl.col("sheet") == sheet) & (pl.col("pair") == pair_label) &
                           (pl.col("metric") == metric))
        if r.height == 0:
            continue
        p_raw = r["p_classical_sidak"].item()  # single-pair, so actually raw
        if p_raw is None or p_raw != p_raw:
            continue
        p = sidak(p_raw, m)
        label = f"{sig(p)}   p = {p_fmt(p)}"
        bracket_y = y_base + k * 0.10 * rng
        add_bracket(ax, ia, ib, bracket_y, label)

    ax.set_xticks(range(len(conditions)))
    ax.set_xticklabels(conditions, rotation=25, ha="right")
    ax.set_ylabel(ylabel)
    ax.set_title(title, fontsize=11, loc="left")


def plot_panel(long, summary, fig_label, sheet, conditions, pairs,
               mark_metric, mark_ylabel, mark_title, out_path):
    n_rows = 1 + len(PROPOSED_ROWS)
    fig, axes = plt.subplots(n_rows, 1,
                             figsize=(2.0 + 1.4 * len(conditions), 3.6 * n_rows),
                             squeeze=False)
    draw_subplot(axes[0, 0], long, summary, sheet, mark_metric,
                 conditions, pairs, mark_ylabel, mark_title)
    for idx, (metric, ylabel, title) in enumerate(PROPOSED_ROWS, start=1):
        draw_subplot(axes[idx, 0], long, summary, sheet, metric,
                     conditions, pairs, ylabel, title)
    fig.suptitle(fig_label, fontsize=13, y=0.995)
    fig.tight_layout(rect=[0, 0, 1, 0.99])
    fig.savefig(out_path, dpi=140, bbox_inches="tight")
    fig.savefig(out_path.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)
    print(f"Wrote {out_path}")


def main():
    long = pl.read_csv(OUT / "per_cell.csv")
    summary = pl.read_csv(OUT / "per_metric_summary.csv")
    for fig_label, sheet, conds, pairs, mark_metric, mark_ylabel, mark_title in PANELS:
        slug = fig_label.split("  ")[0].replace(" ", "_")
        out = FIG / f"{slug}_alt_metrics.png"
        plot_panel(long, summary, fig_label, sheet, conds, pairs,
                   mark_metric, mark_ylabel, mark_title, out)


if __name__ == "__main__":
    sys.exit(main())
