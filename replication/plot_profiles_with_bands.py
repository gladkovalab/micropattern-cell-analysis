"""Six-panel grid of wedge-r intensity profiles (one per sheet) with the
two manuscript-candidate bands highlighted as grey slabs:
  - [17, 32) µm  perinuclear/mid-zone slab
  - [40, 55) µm  rim-zone slab

Profile panel mirrors the upper-left panel of plot_metrics.make_figure
(mean ± SEM per condition). Layout: 2 rows × 3 cols.
"""
from __future__ import annotations

import argparse
import pathlib
import re
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import polars as pl

REPO = pathlib.Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / "replication"))
from plot_metrics import (  # noqa: E402
    SHEET_CONFIG, CONDITION_COLORS, load_template_matching, join_with_metadata,
)

INNER_BAND = (17, 32)   # µm: perinuclear/mid-zone slab
OUTER_BAND = (40, 55)   # µm: rim-zone slab


def wedge_cols(df: pl.DataFrame) -> list[str]:
    pat = re.compile(r"^wedge_r_(\d{2})_(\d{2})um_pct$")
    matches = []
    for c in df.columns:
        m = pat.match(c)
        if m:
            matches.append((int(m.group(1)), c))
    matches.sort()
    return [c for _, c in matches]


def panel(ax, sheet_df: pl.DataFrame, sheet: str, cols: list[str]):
    cfg = SHEET_CONFIG[sheet]
    conditions = cfg["conditions"]
    color_map = {c: CONDITION_COLORS[i % len(CONDITION_COLORS)]
                 for i, c in enumerate(conditions)}
    n_bins = len(cols)
    centers = np.array([i + 0.5 for i in range(n_bins)])

    # Grey slabs first so they sit underneath the lines.
    for lo, hi in (INNER_BAND, OUTER_BAND):
        ax.axvspan(lo, hi, color="0.85", zorder=0, linewidth=0)

    for cond in conditions:
        sub = sheet_df.filter(pl.col("condition") == cond)
        if sub.height == 0:
            continue
        prof = sub.select(cols).to_numpy()
        mean = np.nanmean(prof, axis=0)
        sem = np.nanstd(prof, axis=0, ddof=1) / np.sqrt(np.maximum(prof.shape[0], 1))
        col = color_map[cond]
        ax.plot(centers, mean, color=col, lw=1.6,
                label=f"{cond} (n={sub.height})")
        ax.fill_between(centers, mean - sem, mean + sem,
                        color=col, alpha=0.18, linewidth=0)

    # Annotate each band centre once (only on the first panel).
    ax.set_xlabel("wedge-r (µm from apex)", fontsize=9)
    ax.set_ylabel("mean intensity per bin (% of wedge total)", fontsize=9)
    ax.set_title(f"{sheet}  (n={sheet_df.height})", fontsize=10)
    ax.legend(fontsize=7, loc="upper right")
    ax.grid(alpha=0.3)
    ax.set_xlim(0, n_bins)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--template-matching",
                    default="replication/wedge_r_ks_out_all_denoised/by_well")
    ap.add_argument("--comparisons-xlsx",
                    default=str(REPO / "config/Comparisons_table_v3.xlsx"))
    ap.add_argument("--out",
                    default="replication/figures_wedge_r_ks/profiles_with_bands.png")
    args = ap.parse_args()

    df = load_template_matching(pathlib.Path(args.template_matching))
    df = join_with_metadata(df, pathlib.Path(args.comparisons_xlsx))

    sheet_order = [
        "TRAK isoform (mito)",
        "TRAK isoform (peroxisome)",
        "TRAK isoform (60mer)",
        "TRAK1 helix muts",
        "TRAK2 helix muts",
        "MAPK9 siRNA",
    ]
    cols = wedge_cols(df)

    fig, axes = plt.subplots(2, 3, figsize=(16, 9), sharex=True)
    for ax, sheet in zip(axes.flat, sheet_order):
        sub = df.filter(pl.col("sheet") == sheet)
        if sub.height == 0:
            ax.text(0.5, 0.5, f"no rows: {sheet}", ha="center", va="center",
                    transform=ax.transAxes)
            continue
        panel(ax, sub, sheet, cols)

    # Single shared legend entry for the slabs (annotated as a sub-title).
    fig.suptitle(
        f"Wedge-r 1D intensity profile per sheet  ·  shaded slabs: "
        f"[{INNER_BAND[0]}, {INNER_BAND[1]}) µm (perinuclear/mid-zone) "
        f"and [{OUTER_BAND[0]}, {OUTER_BAND[1]}) µm (rim-zone)",
        fontsize=11)
    plt.tight_layout(rect=[0, 0, 1, 0.96])
    out = pathlib.Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=140)
    plt.close(fig)
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
