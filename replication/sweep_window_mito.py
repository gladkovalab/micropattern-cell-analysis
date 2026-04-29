"""Sweep all 10-µm windows on the wedge-r profile for the mito sheet.

For each window [lo, lo+10) µm with lo in 0..50, compute the % of wedge
intensity in that window per cell, run the nested-ANOVA + Šídák test for
each TRAK pair, and tabulate the p-values. Also produce a small
significance-vs-window plot so we can see where the TRAK1 vs TRAK2 split
becomes significant and how that compares to the other two pairs.
"""
from __future__ import annotations

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
    SHEET_CONFIG, load_template_matching, join_with_metadata, _test_pair,
)

WINDOW_UM = 10
SHEET = "TRAK isoform (mito)"


def main():
    df = load_template_matching(pathlib.Path(
        "replication/wedge_r_ks_out_all_denoised/by_well"))
    df = join_with_metadata(df, REPO / "config/Comparisons_table_v3.xlsx")
    sub = df.filter(pl.col("sheet") == SHEET)
    cfg = SHEET_CONFIG[SHEET]

    rows = []
    for lo in range(0, 60 - WINDOW_UM + 1):
        hi = lo + WINDOW_UM
        cols = [f"wedge_r_{i:02d}_{i+1:02d}um_pct" for i in range(lo, hi)]
        cols = [c for c in cols if c in sub.columns]
        if len(cols) < WINDOW_UM:
            continue
        metric_name = f"_band_{lo:02d}_{hi:02d}"
        cur = sub.with_columns(
            pl.sum_horizontal([pl.col(c) for c in cols]).alias(metric_name))
        ps = {}
        means = {}
        for cond in cfg["conditions"]:
            v = cur.filter(pl.col("condition") == cond)[metric_name].to_numpy().astype(float)
            v = v[np.isfinite(v)]
            means[cond] = float(v.mean()) if v.size else float("nan")
        for pair in cfg["pairs"]:
            ps[f"{pair[0]} vs {pair[1]}"] = _test_pair(
                cur, metric_name, pair, cfg["family_m"])
        rows.append({"lo": lo, "hi": hi, **{f"mean_{c}": means[c] for c in cfg["conditions"]},
                     **ps})

    out = pl.from_dicts(rows)
    print("=== window sweep: % mito intensity in 10-µm slab vs window position ===")
    print(f"sheet={SHEET!r}, family m={cfg['family_m']}\n")
    pair_cols = [f"{a} vs {b}" for a, b in cfg["pairs"]]

    # text table: window, means, three p-values
    print(f"{'window (µm)':<14} {'mean noTRAK':>12} {'mean TRAK1':>12} "
          f"{'mean TRAK2':>12}  ", end="")
    for pc in pair_cols:
        print(f"{pc:>22}  ", end="")
    print()
    for r in out.iter_rows(named=True):
        print(f"  [{r['lo']:2d},{r['hi']:2d})    "
              f"{r['mean_no TRAK']:>11.2f}% {r['mean_TRAK1']:>11.2f}% "
              f"{r['mean_TRAK2']:>11.2f}%  ", end="")
        for pc in pair_cols:
            p = r[pc]
            mark = (("***" if p < 0.001 else "** " if p < 0.01
                     else "*  " if p < 0.05 else "ns ")
                    if np.isfinite(p) else "?  ")
            print(f"{mark} {p:>10.4g}      ", end="")
        print()

    # Figure: -log10(p) vs window center, one line per pair
    fig, ax = plt.subplots(figsize=(11, 5))
    centers = (out["lo"].to_numpy() + out["hi"].to_numpy()) / 2.0
    colors = {"no TRAK vs TRAK1": "#4c78a8",
              "no TRAK vs TRAK2": "#e15759",
              "TRAK1 vs TRAK2":   "#59a14f"}
    for pc in pair_cols:
        ps = out[pc].to_numpy()
        ps = np.where(ps > 0, ps, np.nan)  # avoid log(0)
        ax.plot(centers, -np.log10(ps), label=pc, lw=2, color=colors.get(pc))
    for thr, label in [(0.05, "p=0.05"), (0.01, "p=0.01"), (0.001, "p=0.001")]:
        ax.axhline(-np.log10(thr), color="gray", linestyle="--", alpha=0.5,
                   linewidth=0.8)
        ax.text(60, -np.log10(thr), f" {label}", va="center", fontsize=8,
                color="gray")
    # mark the windows we already named
    for win_lo, win_label in [(20, "20–30 µm"), (45, "arch ±5 µm")]:
        ax.axvspan(win_lo, win_lo + WINDOW_UM, alpha=0.08, color="orange")
        ax.text(win_lo + WINDOW_UM / 2, ax.get_ylim()[1] * 0.95,
                win_label, ha="center", va="top", fontsize=8, color="darkorange")
    ax.set_xlabel("window center (µm from wedge apex)")
    ax.set_ylabel("-log10(p)  [nested ANOVA + Šídák m=3]")
    ax.set_title(f"{SHEET}: significance vs 10-µm window position")
    ax.set_xlim(0, 60)
    ax.legend(loc="upper right", fontsize=9)
    ax.grid(alpha=0.3)
    plt.tight_layout()
    out_png = pathlib.Path("replication/figures_wedge_r_ks/sweep_window_mito.png")
    fig.savefig(out_png, dpi=140)
    plt.close(fig)
    print(f"\nwrote {out_png}")


if __name__ == "__main__":
    main()
