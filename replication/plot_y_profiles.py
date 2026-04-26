"""Visualise the Y-axis projection story for Fig 4B.

Panel A: per-condition mean Y-profile (intensity-fraction vs dY-from-pattern-CoM)
         with shaded ±SEM band. Separate rows for z-sum and MaxIP.
Panel B: per-cell Y-Gini scatter, colored by plate, grouped by condition.
Panel C: per-cell Y-σ, Y-entropy, Y-skewness scatter — companion summary metrics.
Panel D: per-plate mean Y-profiles for no-TRAK vs TRAK2 (confirm effect is not
         driven by one plate alone).
"""
from __future__ import annotations

import pathlib
import re
import sys

import numpy as np
import polars as pl
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

REPO = pathlib.Path(__file__).resolve().parent.parent
NEW_CSV = REPO / "replication" / "overnight_fig4b_out" / "combined_raw.csv"
OLD_CSV = REPO / "replication" / "overnight_out" / "combined.csv"
OUT_DIR = REPO / "replication" / "overnight_fig4b_out" / "figures"

SHEET = "TRAK isoform (mito)"
COND_ORDER = ["no TRAK", "TRAK1", "TRAK2"]
COND_COLOR = {"no TRAK": "#4c78a8", "TRAK1": "#59a14f", "TRAK2": "#e15759"}


def load() -> pl.DataFrame:
    new = pl.read_csv(NEW_CSV)
    old = pl.read_csv(OLD_CSV).filter(pl.col("sheet") == SHEET).select(
        ["path", "plate", "well", "sheet", "condition"])
    df = new.join(old, on="path", how="left").filter(pl.col("condition").is_not_null())
    return df


def extract_profile(df: pl.DataFrame, proj: str, axis: str) -> tuple[np.ndarray, np.ndarray]:
    """axis = 'y' or 'x'. Returns (bin_centers_um, profile_matrix (n_cells, n_bins))."""
    pat = re.compile(rf"^{proj}_{axis}_profile_([+-]\d{{3,4}})um_pct$")
    matches = []
    for c in df.columns:
        m = pat.match(c)
        if m:
            matches.append((int(m.group(1)), c))
    matches.sort()
    bin_edges = np.array([e for e, _ in matches])  # bin LOWER edges in µm
    centers = bin_edges + 0.5  # 1 µm bins
    cols = [c for _, c in matches]
    profile = df.select(cols).to_numpy()  # (n_cells, n_bins)
    return centers, profile


def plot_profile_rows(df: pl.DataFrame, out_png: pathlib.Path):
    fig, axes = plt.subplots(2, 1, figsize=(8, 7), sharex=True)
    for ax, proj in zip(axes, ["zsum", "maxip"]):
        centers, profile = extract_profile(df, proj, "y")
        for cond in COND_ORDER:
            mask = df["condition"].to_numpy() == cond
            P = profile[mask]
            n = mask.sum()
            mean = np.nanmean(P, axis=0)
            sem = np.nanstd(P, axis=0) / np.sqrt(n)
            ax.plot(centers, mean, color=COND_COLOR[cond], linewidth=2,
                    label=f"{cond} (n={n})")
            ax.fill_between(centers, mean - sem, mean + sem,
                            color=COND_COLOR[cond], alpha=0.2)
        ax.axvline(0, color="black", linewidth=0.6, linestyle=":", alpha=0.5)
        ax.set_ylabel("% of total mito signal / µm")
        ax.set_title(f"{proj.upper()} Y-profile (mean ± SEM)")
        ax.legend(loc="upper right", fontsize=9)
        ax.grid(alpha=0.3)
    axes[1].set_xlabel("Y position (µm from pattern CoM; positive = basal)")
    axes[0].annotate("← apical (arch)", xy=(-28, axes[0].get_ylim()[1] * 0.92),
                     fontsize=8, color="gray")
    axes[0].annotate("basal (nucleus) →", xy=(18, axes[0].get_ylim()[1] * 0.92),
                     fontsize=8, color="gray")
    plt.tight_layout()
    fig.savefig(out_png, dpi=130)
    plt.close(fig)


def plot_per_cell_scalars(df: pl.DataFrame, out_png: pathlib.Path):
    """Y-Gini, Y-entropy, Y-σ, Y-skew per cell, colored by plate."""
    metrics = [
        ("zsum_y_gini", "Y-Gini (zsum)"),
        ("maxip_y_gini", "Y-Gini (MaxIP)"),
        ("zsum_y_entropy", "Y-entropy (zsum)"),
        ("maxip_y_entropy", "Y-entropy (MaxIP)"),
        ("zsum_y_sd_u", "Y-σ (zsum, µm)"),
        ("maxip_y_sd_u", "Y-σ (MaxIP, µm)"),
    ]
    plates = sorted(df["plate"].unique().to_list())
    plate_markers = dict(zip(plates, ["o", "s", "D", "^", "v", "P", "X"]))
    fig, axes = plt.subplots(3, 2, figsize=(10, 11))
    for ax, (metric, label) in zip(axes.flat, metrics):
        # jittered scatter by condition, colored by plate
        cond_x = {c: i for i, c in enumerate(COND_ORDER)}
        for cond in COND_ORDER:
            sub = df.filter(pl.col("condition") == cond)
            for plate in plates:
                psub = sub.filter(pl.col("plate") == plate)
                if psub.height == 0:
                    continue
                vals = psub[metric].to_numpy()
                xs = cond_x[cond] + np.random.uniform(-0.15, 0.15, size=len(vals))
                ax.scatter(xs, vals, marker=plate_markers[plate], s=40,
                           color=COND_COLOR[cond], edgecolor="black", linewidth=0.4,
                           alpha=0.8, label=plate if ax is axes.flat[0] else None)
            # mean bar
            mean = sub[metric].mean()
            if mean is not None:
                ax.hlines(mean, cond_x[cond] - 0.25, cond_x[cond] + 0.25,
                          color="black", linewidth=2.5, zorder=5)
        ax.set_xticks(list(cond_x.values()))
        ax.set_xticklabels(COND_ORDER)
        ax.set_ylabel(label)
        ax.grid(axis="y", alpha=0.3)
    axes.flat[0].legend(loc="best", fontsize=7, title="plate", title_fontsize=8)
    plt.tight_layout()
    fig.savefig(out_png, dpi=130)
    plt.close(fig)


def plot_per_plate_profiles(df: pl.DataFrame, out_png: pathlib.Path):
    """Per-plate mean Y-profiles to check plate-consistency of the effect."""
    plates = sorted(df["plate"].unique().to_list())
    fig, axes = plt.subplots(len(plates), 2, figsize=(10, 2.5 * len(plates)),
                             sharex=True)
    if len(plates) == 1:
        axes = np.array([axes])
    for row, plate in enumerate(plates):
        for col, proj in enumerate(["zsum", "maxip"]):
            ax = axes[row, col]
            centers, profile = extract_profile(df, proj, "y")
            plate_mask = df["plate"].to_numpy() == plate
            n_by_cond = {}
            for cond in ("no TRAK", "TRAK2"):  # the pair of interest
                m = plate_mask & (df["condition"].to_numpy() == cond)
                if not m.any():
                    continue
                P = profile[m]
                n_by_cond[cond] = m.sum()
                mean = np.nanmean(P, axis=0)
                sem = np.nanstd(P, axis=0) / np.sqrt(m.sum())
                ax.plot(centers, mean, color=COND_COLOR[cond], linewidth=1.8,
                        label=f"{cond} (n={m.sum()})")
                ax.fill_between(centers, mean - sem, mean + sem,
                                color=COND_COLOR[cond], alpha=0.2)
            ax.axvline(0, color="black", linewidth=0.5, linestyle=":", alpha=0.4)
            ax.set_title(f"{plate.replace('250', '').replace('_patterned_plate', ' plate ')[:40]} — {proj.upper()}",
                         fontsize=9)
            ax.grid(alpha=0.3)
            if col == 0:
                ax.set_ylabel("% signal / µm")
            if row == len(plates) - 1:
                ax.set_xlabel("Y (µm from pattern CoM)")
            ax.legend(loc="upper right", fontsize=7)
    plt.tight_layout()
    fig.savefig(out_png, dpi=130)
    plt.close(fig)


def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    df = load()
    print(f"Loaded {df.height} cells across {df['condition'].unique().to_list()}")
    print(f"Plates: {df['plate'].unique().to_list()}")

    plot_profile_rows(df, OUT_DIR / "fig4B_y_profile_mean_sem.png")
    plot_per_cell_scalars(df, OUT_DIR / "fig4B_y_scalars_per_cell.png")
    plot_per_plate_profiles(df, OUT_DIR / "fig4B_y_profile_per_plate.png")

    print(f"\nPlots written to {OUT_DIR}/")
    for p in sorted(OUT_DIR.glob("fig4B_y_*.png")):
        print(f"  {p.name}")


if __name__ == "__main__":
    main()
