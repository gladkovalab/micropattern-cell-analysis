"""Fig 4 / Fig S11 panel plots with MaxIP rows added.

Extension of `plot_all_panels.py` that reads BOTH data sources:
  - replication/derived_metrics_out/per_cell.csv (CSV-derived z-sum metrics)
  - replication/overnight_out/combined.csv (raw-ND2-derived z-sum + MaxIP metrics)

Produces a 5-row per-panel plot:
  1. Mark's current metric (z-sum perinuclear for S11, z-sum peripheral for 4)
  2. Proposed  — z-sum peripheral − perinuclear (diff)
  3. Alternative — z-sum peripheral ÷ perinuclear (ratio)
  4. MaxIP diff
  5. MaxIP ratio

Outputs: replication/figures/*_alt_metrics_maxip.{png,pdf}

Same Šídák family sizes and pair choices as `plot_all_panels.py`.
"""
from __future__ import annotations

import pathlib
import sys

import fastexcel
import matplotlib.pyplot as plt
import numpy as np
import polars as pl

REPO = pathlib.Path(__file__).resolve().parent.parent
OUT_DERIVED = REPO / "replication" / "derived_metrics_out" / "per_cell.csv"
OUT_OVERNIGHT = REPO / "replication" / "overnight_out" / "combined.csv"
FIG = REPO / "replication" / "figures"
FIG.mkdir(parents=True, exist_ok=True)

PANELS = [
    ("Fig S11 C  peroxisome  TRAK isoforms", "TRAK isoform (peroxisome)",
     ["no TRAK", "TRAK1", "TRAK2"],
     [("no TRAK", "TRAK1"), ("no TRAK", "TRAK2"), ("TRAK1", "TRAK2")],
     "fig_s11"),
    ("Fig S11 D  TRAK1 helix mutants", "TRAK1 helix muts",
     ["T1 wt", "T1 mDRH", "T1 mDRH / dSp"],
     [("T1 wt", "T1 mDRH"), ("T1 mDRH", "T1 mDRH / dSp")],
     "fig_s11"),
    ("Fig S11 E  TRAK2 helix mutants", "TRAK2 helix muts",
     ["TRAK2", "TRAK2 mDRH", "TRAK2 mDRH mSpindly"],
     [("TRAK2", "TRAK2 mDRH"), ("TRAK2 mDRH", "TRAK2 mDRH mSpindly")],
     "fig_s11"),
    ("Fig S11 F  MAPK9/JNK2 siRNA + arsenite", "MAPK9 siRNA",
     ["ctrl ctrl", "ctrl Ars", "MAPK9 ctrl", "MAPK9 Ars"],
     [("ctrl ctrl", "ctrl Ars"), ("MAPK9 ctrl", "MAPK9 Ars"), ("ctrl Ars", "MAPK9 Ars")],
     "fig_s11"),
    ("Fig 4B  TRAK isoforms (mito)", "TRAK isoform (mito)",
     ["no TRAK", "TRAK1", "TRAK2"],
     [("no TRAK", "TRAK1"), ("no TRAK", "TRAK2"), ("TRAK1", "TRAK2")],
     "fig_4"),
    ("Fig 4C  TRAK1 helix mutants", "TRAK1 helix muts",
     ["T1 wt", "T1 mDRH", "T1 mDRH / dSp"],
     [("T1 wt", "T1 mDRH"), ("T1 mDRH", "T1 mDRH / dSp")],
     "fig_4"),
    ("Fig 4D  TRAK2 helix mutants", "TRAK2 helix muts",
     ["TRAK2", "TRAK2 mDRH", "TRAK2 mDRH mSpindly"],
     [("TRAK2", "TRAK2 mDRH"), ("TRAK2 mDRH", "TRAK2 mDRH mSpindly")],
     "fig_4"),
    ("Fig 4E  MAPK9/JNK2 siRNA + arsenite", "MAPK9 siRNA",
     ["ctrl ctrl", "ctrl Ars", "MAPK9 ctrl", "MAPK9 Ars"],
     [("ctrl ctrl", "MAPK9 ctrl"), ("ctrl ctrl", "ctrl Ars"), ("ctrl ctrl", "MAPK9 Ars")],
     "fig_4"),
]

# For S11 (perinuclear) panels Mark's baseline is raw z-sum perinuclear.
# For Fig 4 (peripheral) panels Mark's baseline is denoised z-sum peripheral.
MARK_METRIC = {
    "fig_s11": ("raw_perinuclear_5um", "% mito within 5 µm of nucleus",
                "Mark's current metric — perinuclear 5 µm (raw z-sum)"),
    "fig_4": ("den_peripheral_5um", "% mito within 5 µm of arch",
              "Mark's current metric — peripheral 5 µm (denoised z-sum)"),
}


def p_fmt(p):
    if p is None or p != p:
        return "—"
    if p < 0.0001: return "< 0.0001"
    if p < 0.001: return f"{p:.4f}"
    if p < 0.01: return f"{p:.3f}"
    return f"{p:.2f}"


def sig(p):
    if p is None or p != p: return "n/a"
    if p < 0.001: return "***"
    if p < 0.01: return "**"
    if p < 0.05: return "*"
    return "ns"


def sidak(p_raw: float, m: int) -> float:
    return 1 - (1 - p_raw) ** m


def nested_t(conds_cells: dict[str, np.ndarray]) -> float | None:
    """Return raw p for all-sheet nested one-way ANOVA + pairwise t between
    first two populated conditions. NB: this is ONLY used to recompute p-values
    lazily when we have combined-csv data without a pre-computed summary.
    For the main memo we rely on summary.csv / per_metric_summary.csv."""
    # Simplified: if we have ≥ 2 conds with plate info, approximate Šídák.
    # We're not actually recomputing here — plots pull from summary files.
    return None


def load_metrics() -> pl.DataFrame:
    """Build a unified per-cell table combining derived + overnight.
    Adds composite metrics (diff, ratio) for both z-sum and MaxIP."""
    derived = pl.read_csv(OUT_DERIVED)
    overnight = pl.read_csv(OUT_OVERNIGHT) if OUT_OVERNIGHT.exists() else None

    if overnight is not None and overnight.height:
        # overnight has one row per (sheet, condition, plate, well, path) already
        # (from rebuild_combined) — add composites
        eps = 1e-9
        extra = []
        for proj in ("zsum", "maxip"):
            for mask in ("crop",):
                peri = f"{proj}_{mask}_peripheral_5um_pct"
                nuc = f"{proj}_{mask}_perinuclear_5um_pct"
                if peri in overnight.columns and nuc in overnight.columns:
                    extra.append((pl.col(peri) - pl.col(nuc)).alias(f"{proj}_{mask}_peri_minus_nuc"))
                    extra.append((pl.col(peri) / (pl.col(nuc) + eps)).alias(f"{proj}_{mask}_peri_over_nuc"))
        if extra:
            overnight = overnight.with_columns(extra)

    return derived, overnight


def get_classical_sidak_p(summary: pl.DataFrame, sheet: str, pair: tuple, metric: str,
                          m: int) -> float | None:
    """Look up the raw p in summary and apply Šídák with family size m."""
    pair_label = f"{pair[0]} vs {pair[1]}"
    r = summary.filter(
        (pl.col("sheet") == sheet) &
        (pl.col("pair") == pair_label) &
        (pl.col("metric") == metric)
    )
    if r.height == 0:
        return None
    p_raw = r["p_classical_sidak"].item() if "p_classical_sidak" in r.columns else r["p_sidak"].item()
    if p_raw is None or p_raw != p_raw:
        return None
    return sidak(p_raw, m)


def add_bracket(ax, x0, x1, y, label, pad=0.02):
    h = 0.03 * (ax.get_ylim()[1] - ax.get_ylim()[0])
    ax.plot([x0, x0, x1, x1], [y, y + h, y + h, y], color="black", lw=1.0)
    ax.text((x0 + x1) / 2, y + h + pad * (ax.get_ylim()[1] - ax.get_ylim()[0]),
            label, ha="center", va="bottom", fontsize=9)


def draw_row(ax, df: pl.DataFrame, sheet: str, metric: str, conditions: list[str],
             pairs: list[tuple[str, str]], ylabel: str, title: str,
             summary_df: pl.DataFrame):
    sub = df.filter((pl.col("sheet") == sheet) & pl.col(metric).is_not_null() &
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
        std = g.group_by("plate").agg(pl.col(metric).mean())[metric].std()
        n_p = max(1, g["plate"].n_unique())
        sem = (std / np.sqrt(n_p)) if (std is not None and np.isfinite(std)) else 0.0
        ax.errorbar([i], [mu], yerr=[sem], fmt="_", color="black", capsize=8,
                    elinewidth=2, markersize=28, zorder=4)

    if sub.height == 0:
        ax.text(0.5, 0.5, "(no data for this metric)", ha="center", va="center",
                transform=ax.transAxes, fontsize=11, color="gray")
        ax.set_xticks(range(len(conditions)))
        ax.set_xticklabels(conditions, rotation=25, ha="right")
        ax.set_ylabel(ylabel)
        ax.set_title(title, fontsize=11, loc="left")
        return

    ys_all = sub[metric].to_numpy()
    y_top = np.nanmax(ys_all)
    y_bot = np.nanmin(ys_all)
    rng = max(y_top - y_bot, 1e-6)
    ax.set_ylim(y_bot - 0.05 * rng, y_top + 0.35 * rng * (1 + 0.3 * len(pairs)))
    y_base = y_top + 0.05 * rng

    m = len(pairs)
    for k, (a, b) in enumerate(pairs):
        if a not in conditions or b not in conditions:
            continue
        ia, ib = conditions.index(a), conditions.index(b)
        p = get_classical_sidak_p(summary_df, sheet, (a, b), metric, m)
        if p is None:
            continue
        label = f"{sig(p)}   p = {p_fmt(p)}"
        bracket_y = y_base + k * 0.10 * rng
        add_bracket(ax, ia, ib, bracket_y, label)

    ax.set_xticks(range(len(conditions)))
    ax.set_xticklabels(conditions, rotation=25, ha="right")
    ax.set_ylabel(ylabel)
    ax.set_title(title, fontsize=11, loc="left")


def build_summary() -> pl.DataFrame:
    """Concatenate derived + overnight summaries so every metric has a lookup."""
    dfs = []
    derived_sum = REPO / "replication" / "derived_metrics_out" / "per_metric_summary.csv"
    if derived_sum.exists():
        d = pl.read_csv(derived_sum).select([
            pl.col("sheet"), pl.col("pair"), pl.col("metric"),
            pl.col("p_classical_sidak")
        ])
        dfs.append(d)
    overnight_sum = REPO / "replication" / "overnight_eval_out" / "summary.csv"
    if overnight_sum.exists():
        o = pl.read_csv(overnight_sum).select([
            pl.col("sheet"), pl.col("pair"), pl.col("metric"),
            pl.col("p_raw").alias("p_classical_sidak")
        ])
        dfs.append(o)
    return pl.concat(dfs, how="diagonal_relaxed") if dfs else pl.DataFrame()


def main():
    derived, overnight = load_metrics()
    summary = build_summary()

    for fig_label, sheet, conditions, pairs, fig_group in PANELS:
        mark_metric, mark_ylabel, mark_title = MARK_METRIC[fig_group]

        rows = [(derived, mark_metric, mark_ylabel, mark_title)]
        rows.append((derived, "den_peri_minus_nuc", "peripheral − perinuclear (pp)",
                     "z-sum diff (denoised) — peripheral 5 µm − perinuclear 5 µm"))
        rows.append((derived, "den_peri_over_nuc", "peripheral ÷ perinuclear",
                     "z-sum ratio (denoised)"))
        if overnight is not None:
            rows.append((overnight, "maxip_crop_peri_minus_nuc",
                         "peripheral − perinuclear (pp)",
                         "MaxIP diff — peripheral 5 µm − perinuclear 5 µm"))
            rows.append((overnight, "maxip_crop_peri_over_nuc",
                         "peripheral ÷ perinuclear",
                         "MaxIP ratio"))

        n_rows = len(rows)
        fig, axes = plt.subplots(n_rows, 1,
                                 figsize=(2.0 + 1.4 * len(conditions), 3.6 * n_rows),
                                 squeeze=False)
        for r_idx, (src_df, metric, ylabel, title) in enumerate(rows):
            draw_row(axes[r_idx, 0], src_df, sheet, metric, conditions, pairs,
                     ylabel, title, summary)
        fig.suptitle(fig_label, fontsize=13, y=0.995)
        fig.tight_layout(rect=[0, 0, 1, 0.99])
        slug = fig_label.split("  ")[0].replace(" ", "_")
        out = FIG / f"{slug}_alt_metrics_maxip.png"
        fig.savefig(out, dpi=140, bbox_inches="tight")
        fig.savefig(out.with_suffix(".pdf"), bbox_inches="tight")
        plt.close(fig)
        print(f"Wrote {out}")


if __name__ == "__main__":
    sys.exit(main())
