"""Plot the 1D wedge-r profiles with mean ± SEM, mirror of plot_y_profiles.py."""
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
NEW_CSV = REPO / "replication" / "overnight_fig4b_v2_out" / "combined_raw.csv"
OLD_CSV = REPO / "replication" / "overnight_out" / "combined.csv"
OUT_DIR = REPO / "replication" / "overnight_fig4b_v2_out" / "figures"

SHEET = "TRAK isoform (mito)"
COND_ORDER = ["no TRAK", "TRAK1", "TRAK2"]
COND_COLOR = {"no TRAK": "#4c78a8", "TRAK1": "#59a14f", "TRAK2": "#e15759"}


def load() -> pl.DataFrame:
    new = pl.read_csv(NEW_CSV)
    old = pl.read_csv(OLD_CSV).filter(pl.col("sheet") == SHEET).select(
        ["path", "plate", "well", "sheet", "condition"])
    return new.join(old, on="path", how="left").filter(pl.col("condition").is_not_null())


def extract_wedge_r_profile(df: pl.DataFrame, proj: str):
    pat = re.compile(rf"^{proj}_wedge_r_(\d{{2}})_(\d{{2}})um_pct$")
    matches = []
    for c in df.columns:
        m = pat.match(c)
        if m:
            matches.append((int(m.group(1)), c))
    matches.sort()
    centers = np.array([lo + 0.5 for lo, _ in matches])
    cols = [c for _, c in matches]
    profile = df.select(cols).to_numpy()
    return centers, profile


def plot_profile_rows(df: pl.DataFrame, out_png: pathlib.Path):
    fig, axes = plt.subplots(2, 1, figsize=(8.5, 7), sharex=True)
    for ax, proj in zip(axes, ["zsum", "maxip"]):
        centers, profile = extract_wedge_r_profile(df, proj)
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
        # Annotate the band boundaries used for the fraction metrics
        for x, label in [(20, "perinuclear band starts"), (35, "peripheral band starts")]:
            ax.axvline(x, color="gray", linewidth=0.6, linestyle="--", alpha=0.5)
        ax.set_ylabel("% of wedge mito signal / µm")
        ax.set_title(f"{proj.upper()} wedge-r profile (mean ± SEM)")
        ax.legend(loc="upper right", fontsize=9)
        ax.grid(alpha=0.3)
    axes[1].set_xlabel("r (µm from wedge apex = pattern bottom; arch ≈ 45-50 µm)")
    axes[0].annotate("← stalk", xy=(2, axes[0].get_ylim()[1] * 0.92),
                     fontsize=8, color="gray")
    axes[0].annotate("nucleus zone", xy=(24, axes[0].get_ylim()[1] * 0.92),
                     fontsize=8, color="gray")
    axes[0].annotate("arch →", xy=(46, axes[0].get_ylim()[1] * 0.92),
                     fontsize=8, color="gray")
    plt.tight_layout()
    fig.savefig(out_png, dpi=130)
    plt.close(fig)


def plot_per_cell_scalars(df: pl.DataFrame, out_png: pathlib.Path):
    metrics = [
        ("zsum_wedge_r_gini", "wedge-r-Gini (zsum)"),
        ("maxip_wedge_r_gini", "wedge-r-Gini (MaxIP)"),
        ("zsum_wedge_r_entropy", "wedge-r-entropy (zsum)"),
        ("maxip_wedge_r_entropy", "wedge-r-entropy (MaxIP)"),
        ("zsum_wedge_r_sd_um", "wedge-r-σ (zsum, µm)"),
        ("maxip_wedge_r_sd_um", "wedge-r-σ (MaxIP, µm)"),
    ]
    plates = sorted(df["plate"].unique().to_list())
    plate_markers = dict(zip(plates, ["o", "s", "D", "^", "v", "P", "X"]))
    fig, axes = plt.subplots(3, 2, figsize=(10, 11))
    cond_x = {c: i for i, c in enumerate(COND_ORDER)}
    for ax, (metric, label) in zip(axes.flat, metrics):
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
    plates = sorted(df["plate"].unique().to_list())
    fig, axes = plt.subplots(len(plates), 2, figsize=(10, 2.5 * len(plates)),
                             sharex=True)
    if len(plates) == 1:
        axes = np.array([axes])
    for row, plate in enumerate(plates):
        for col, proj in enumerate(["zsum", "maxip"]):
            ax = axes[row, col]
            centers, profile = extract_wedge_r_profile(df, proj)
            plate_mask = df["plate"].to_numpy() == plate
            for cond in ("no TRAK", "TRAK2"):
                m = plate_mask & (df["condition"].to_numpy() == cond)
                if not m.any():
                    continue
                P = profile[m]
                mean = np.nanmean(P, axis=0)
                sem = np.nanstd(P, axis=0) / np.sqrt(m.sum())
                ax.plot(centers, mean, color=COND_COLOR[cond], linewidth=1.8,
                        label=f"{cond} (n={m.sum()})")
                ax.fill_between(centers, mean - sem, mean + sem,
                                color=COND_COLOR[cond], alpha=0.2)
            ax.set_title(f"{plate.replace('250', '').replace('_patterned_plate', ' plate ')[:40]} — {proj.upper()}",
                         fontsize=9)
            ax.grid(alpha=0.3)
            if col == 0:
                ax.set_ylabel("% / µm")
            if row == len(plates) - 1:
                ax.set_xlabel("wedge-r (µm)")
            ax.legend(loc="upper right", fontsize=7)
    plt.tight_layout()
    fig.savefig(out_png, dpi=130)
    plt.close(fig)


def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    df = load()
    print(f"Loaded {df.height} cells")

    plot_profile_rows(df, OUT_DIR / "fig4B_wedge_r_profile_mean_sem.png")
    plot_per_cell_scalars(df, OUT_DIR / "fig4B_wedge_r_scalars_per_cell.png")
    plot_per_plate_profiles(df, OUT_DIR / "fig4B_wedge_r_profile_per_plate.png")

    print(f"\nPlots written to {OUT_DIR}/")
    for p in sorted(OUT_DIR.glob("fig4B_wedge_r_*.png")):
        print(f"  {p.name}")


if __name__ == "__main__":
    main()
