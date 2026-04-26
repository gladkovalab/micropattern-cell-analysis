"""Per-sheet, per-pair plot generator for the final pipeline output.

For every sheet × pair, produces:
  * `{sheet}_y_profile.png`     — Y-axis profile mean ± SEM (zsum + maxip rows),
                                   one curve per condition in the sheet
  * `{sheet}_wedge_r_profile.png` — wedge-r profile mean ± SEM (zsum + maxip rows)
  * `{sheet}_scalars.png`       — per-cell strip plots of Y-Gini, wedge-r-Gini,
                                   wedge-r-σ for each projection, conditions side
                                   by side, plates as marker shapes

Reads `overnight_final_out/combined_raw.csv` joined with the canonical
`overnight_out/combined.csv` for sheet/condition/plate/well metadata.
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
NEW_CSV = REPO / "replication" / "overnight_final_out" / "combined_raw.csv"
OLD_CSV = REPO / "replication" / "overnight_out" / "combined.csv"
PEROX_CSV = REPO / "replication" / "overnight_final_out" / "peroxisome_metadata.csv"
EVAL_CSV = REPO / "replication" / "overnight_final_out" / "evaluation_summary.csv"
OUT_DIR = REPO / "replication" / "overnight_final_out" / "figures"

CONDITION_COLORS = [
    "#4c78a8", "#59a14f", "#e15759", "#f28e2b", "#b07aa1", "#76b7b2",
    "#edc949", "#9c755f",
]

# Per-sheet condition order, pairs to annotate, and Šídák family size.
SHEET_CONFIG = {
    "TRAK isoform (mito)": {
        "conditions": ["no TRAK", "TRAK1", "TRAK2"],
        "pairs": [("no TRAK", "TRAK1"), ("no TRAK", "TRAK2"), ("TRAK1", "TRAK2")],
        "family_m": 3,
    },
    "TRAK isoform (peroxisome)": {
        "conditions": ["no TRAK", "TRAK1", "TRAK2"],
        "pairs": [("no TRAK", "TRAK1"), ("no TRAK", "TRAK2"), ("TRAK1", "TRAK2")],
        "family_m": 3,
    },
    "TRAK isoform (60mer)": {
        "conditions": ["no TRAK", "TRAK1", "TRAK2"],
        "pairs": [("no TRAK", "TRAK1"), ("no TRAK", "TRAK2"), ("TRAK1", "TRAK2")],
        "family_m": 3,
    },
    "TRAK1 helix muts": {
        "conditions": ["T1 wt", "T1 mDRH", "T1 mDRH / dSp"],
        "pairs": [("T1 wt", "T1 mDRH"), ("T1 mDRH", "T1 mDRH / dSp")],
        "family_m": 2,
    },
    "TRAK2 helix muts": {
        "conditions": ["TRAK2", "TRAK2 mDRH", "TRAK2 mDRH mSpindly"],
        "pairs": [("TRAK2", "TRAK2 mDRH"),
                  ("TRAK2 mDRH", "TRAK2 mDRH mSpindly")],
        "family_m": 2,
    },
    "MAPK9 siRNA": {
        "conditions": ["ctrl ctrl", "ctrl Ars", "MAPK9 ctrl", "MAPK9 Ars"],
        "pairs": [("ctrl ctrl", "ctrl Ars"),
                  ("ctrl ctrl", "MAPK9 ctrl"),
                  ("MAPK9 ctrl", "MAPK9 Ars")],
        "family_m": 3,
    },
}

# Family size used by evaluate_final.py when it wrote evaluation_summary.csv.
EVAL_FAMILY_M = {
    "TRAK isoform (mito)": 3,
    "TRAK isoform (peroxisome)": 3,
    "TRAK isoform (60mer)": 3,
    "TRAK1 helix muts": 2,
    "TRAK2 helix muts": 2,
    "MAPK9 siRNA": 5,
}


def slug(s: str) -> str:
    s = re.sub(r"[^A-Za-z0-9_]+", "_", s)
    return s.strip("_")


def load() -> pl.DataFrame:
    new = pl.read_csv(NEW_CSV)
    cols = ["path", "plate", "well", "sheet", "condition"]
    base = pl.read_csv(OLD_CSV).select(cols)
    parts = [base]
    if PEROX_CSV.exists():
        parts.append(pl.read_csv(PEROX_CSV).select(cols))
    sixtymer = REPO / "replication" / "overnight_final_out" / "sixtymer_metadata.csv"
    if sixtymer.exists():
        parts.append(pl.read_csv(sixtymer).select(cols))
    meta = pl.concat(parts, how="vertical") if len(parts) > 1 else base
    return new.join(meta, on="path", how="left").filter(
        pl.col("condition").is_not_null())


def extract_profile(df: pl.DataFrame, proj: str, kind: str):
    """kind = 'y' (image-Y axis) or 'wedge_r' (polar r from wedge apex)."""
    if kind == "y":
        pat = re.compile(rf"^{proj}_y_profile_([+-]\d{{3,4}})um_pct$")
        matches = []
        for c in df.columns:
            m = pat.match(c)
            if m:
                matches.append((int(m.group(1)), c))
        matches.sort()
        centers = np.array([lo + 0.5 for lo, _ in matches])
    elif kind == "wedge_r":
        pat = re.compile(rf"^{proj}_wedge_r_(\d{{2}})_(\d{{2}})um_pct$")
        matches = []
        for c in df.columns:
            m = pat.match(c)
            if m:
                matches.append((int(m.group(1)), c))
        matches.sort()
        centers = np.array([lo + 0.5 for lo, _ in matches])
    else:
        raise ValueError(kind)
    cols = [c for _, c in matches]
    profile = df.select(cols).to_numpy()
    return centers, profile


def plot_profile(df_sheet: pl.DataFrame, conditions: list[str], kind: str,
                 title_prefix: str, out_png: pathlib.Path):
    fig, axes = plt.subplots(2, 1, figsize=(9, 7), sharex=True)
    color_map = {c: CONDITION_COLORS[i % len(CONDITION_COLORS)]
                 for i, c in enumerate(conditions)}
    for ax, proj in zip(axes, ["zsum", "maxip"]):
        centers, profile = extract_profile(df_sheet, proj, kind)
        if centers.size == 0:
            ax.text(0.5, 0.5, f"no {proj} {kind} columns",
                    transform=ax.transAxes, ha="center")
            continue
        for cond in conditions:
            mask = df_sheet["condition"].to_numpy() == cond
            n = mask.sum()
            if n == 0:
                continue
            P = profile[mask]
            mean = np.nanmean(P, axis=0)
            sem = np.nanstd(P, axis=0) / np.sqrt(n)
            ax.plot(centers, mean, color=color_map[cond], linewidth=2,
                    label=f"{cond} (n={n})")
            ax.fill_between(centers, mean - sem, mean + sem,
                            color=color_map[cond], alpha=0.18)
        ax.set_ylabel("% mito signal / µm")
        ax.set_title(f"{title_prefix} · {proj.upper()}")
        ax.legend(loc="upper right", fontsize=8)
        ax.grid(alpha=0.3)
    if kind == "y":
        axes[1].set_xlabel("Y position (µm from pattern CoM; positive = stalk side)")
    else:
        axes[1].set_xlabel("r (µm from wedge apex; arch ≈ 45-50 µm)")
        for ax in axes:
            for x in (20, 35):
                ax.axvline(x, color="gray", linewidth=0.6, linestyle="--", alpha=0.4)
    plt.tight_layout()
    fig.savefig(out_png, dpi=130)
    plt.close(fig)


def plot_wedge_r_cdf(df_sheet: pl.DataFrame, conditions: list[str],
                     title_prefix: str, out_png: pathlib.Path):
    """CDF view of the wedge-r distribution that the KS metric is built on.

    Per-cell CDFs are the cumsum of the 60 wedge_r_*_um_pct bins (already
    normalized to sum=100% per cell). Mean ± SEM is computed per-bin across
    cells in each condition. The dashed black reference is the analytical
    area-uniform CDF for a 45° sector, ((i+1)/60)² — validated to match
    stored per-cell KS to within ±0.016 across all 404 cells (r=0.9999).

    Right column shows (CDF − reference); KS = max|Δ| is the height of the
    most extreme excursion of each curve from zero.
    """
    fig, axes = plt.subplots(2, 2, figsize=(13, 8), sharex=True)
    color_map = {c: CONDITION_COLORS[i % len(CONDITION_COLORS)]
                 for i, c in enumerate(conditions)}
    for row, proj in enumerate(["zsum", "maxip"]):
        centers, profile = extract_profile(df_sheet, proj, "wedge_r")
        if centers.size == 0:
            for col in (0, 1):
                axes[row, col].text(0.5, 0.5, f"no {proj} columns",
                                    transform=axes[row, col].transAxes,
                                    ha="center")
            continue
        n_bins = profile.shape[1]
        row_sums = np.nansum(profile, axis=1, keepdims=True)
        ok = (row_sums.flatten() > 0)
        cdf_obs_all = np.full_like(profile, np.nan, dtype=float)
        cdf_obs_all[ok] = np.nancumsum(profile[ok] / row_sums[ok], axis=1)
        cdf_uni = (np.arange(1, n_bins + 1) ** 2) / (n_bins ** 2)

        ax_cdf = axes[row, 0]
        ax_dev = axes[row, 1]
        for cond in conditions:
            mask = df_sheet["condition"].to_numpy() == cond
            n = int(mask.sum())
            if n == 0:
                continue
            C = cdf_obs_all[mask]
            mean = np.nanmean(C, axis=0)
            sem = np.nanstd(C, axis=0) / np.sqrt(n)
            ax_cdf.plot(centers, mean, color=color_map[cond], linewidth=2,
                        label=f"{cond} (n={n})")
            ax_cdf.fill_between(centers, mean - sem, mean + sem,
                                color=color_map[cond], alpha=0.18)
            dev = mean - cdf_uni
            ax_dev.plot(centers, dev, color=color_map[cond], linewidth=2,
                        label=f"{cond} (n={n})")
            ax_dev.fill_between(centers, dev - sem, dev + sem,
                                color=color_map[cond], alpha=0.18)

        ax_cdf.plot(centers, cdf_uni, color="black", linestyle="--",
                    linewidth=1.5, label="area-uniform reference")
        ax_cdf.set_ylabel(f"CDF · {proj.upper()}")
        ax_cdf.legend(loc="lower right", fontsize=7)
        ax_cdf.grid(alpha=0.3)
        ax_cdf.set_ylim(0, 1.02)

        ax_dev.axhline(0, color="black", linestyle="--", linewidth=1)
        ax_dev.set_ylabel(f"CDF − reference · {proj.upper()}")
        ax_dev.grid(alpha=0.3)

        for ax in (ax_cdf, ax_dev):
            for x in (20, 35):
                ax.axvline(x, color="gray", linewidth=0.7, linestyle=":")

    axes[0, 0].set_title("CDF of wedge-r intensity")
    axes[0, 1].set_title("CDF − area-uniform reference  (KS = max|Δ|)")
    axes[1, 0].set_xlabel("r (µm from wedge apex; arch ≈ 45-50 µm)")
    axes[1, 1].set_xlabel("r (µm from wedge apex; arch ≈ 45-50 µm)")
    fig.suptitle(title_prefix, fontsize=11)
    plt.tight_layout(rect=[0, 0.02, 1, 0.96])
    fig.savefig(out_png, dpi=130)
    plt.close(fig)


def make_p_lookup(eval_df: pl.DataFrame):
    def lookup(sheet, metric, cond_a, cond_b, family_m):
        eval_m = EVAL_FAMILY_M.get(sheet)
        if eval_m is None:
            return None
        sub = eval_df.filter((pl.col("sheet") == sheet) &
                             (pl.col("metric") == metric))
        for label in (f"{cond_a} vs {cond_b}", f"{cond_b} vs {cond_a}"):
            r = sub.filter(pl.col("pair") == label)
            if r.height > 0:
                p_eval = r["p"][0]
                if not np.isfinite(p_eval):
                    return None
                if family_m == eval_m:
                    return float(p_eval)
                return float(1 - (1 - p_eval) ** (family_m / eval_m))
        return None
    return lookup


def _format_p(p: float) -> tuple[str, str]:
    sig = "***" if p < 0.001 else "**" if p < 0.01 else "*" if p < 0.05 else "ns"
    if p < 0.001:
        ptxt = f"p<0.001"
    else:
        ptxt = f"p={p:.3f}"
    return f"{ptxt} {sig}", sig


def plot_scalars(df_sheet: pl.DataFrame, sheet: str, conditions: list[str],
                 pairs: list[tuple[str, str]], family_m: int,
                 p_lookup, title_prefix: str, out_png: pathlib.Path):
    metrics = [
        ("zsum_peripheral_5um_pct",    "peripheral 5 µm % (zsum)"),
        ("maxip_peripheral_5um_pct",   "peripheral 5 µm % (MaxIP)"),
        ("zsum_perinuclear_5um_pct",   "perinuclear 5 µm % (zsum)"),
        ("maxip_perinuclear_5um_pct",  "perinuclear 5 µm % (MaxIP)"),
        ("zsum_wedge_r_gini",          "wedge-r Gini (zsum)"),
        ("maxip_wedge_r_gini",         "wedge-r Gini (MaxIP)"),
        ("zsum_wedge_r_ks_vs_uniform", "wedge-r KS vs uniform (zsum)"),
        ("maxip_wedge_r_ks_vs_uniform","wedge-r KS vs uniform (MaxIP)"),
    ]
    plates = sorted(df_sheet["plate"].unique().to_list())
    plate_markers = dict(zip(plates,
                             ["o", "s", "D", "^", "v", "P", "X", "<", ">", "*"]))
    color_map = {c: CONDITION_COLORS[i % len(CONDITION_COLORS)]
                 for i, c in enumerate(conditions)}
    cond_idx = {c: i for i, c in enumerate(conditions)}

    # Text scale: 1.5× the previous sizes
    SC = 1.5
    fig, axes = plt.subplots(4, 2, figsize=(13, 22))
    rng = np.random.default_rng(0)
    for ax, (m, label) in zip(axes.flat, metrics):
        if m not in df_sheet.columns:
            ax.text(0.5, 0.5, f"no column {m}",
                    transform=ax.transAxes, ha="center", fontsize=10 * SC)
            continue
        for j, cond in enumerate(conditions):
            sub = df_sheet.filter(pl.col("condition") == cond)
            for plate in plates:
                psub = sub.filter(pl.col("plate") == plate)
                if psub.height == 0:
                    continue
                vals = psub[m].to_numpy()
                xs = j + rng.uniform(-0.18, 0.18, size=len(vals))
                ax.scatter(xs, vals, marker=plate_markers[plate], s=42,
                           color=color_map[cond], edgecolor="black",
                           linewidth=0.4, alpha=0.85)
            mn = sub[m].mean()
            if mn is not None:
                ax.hlines(mn, j - 0.3, j + 0.3, color="black",
                          linewidth=2.5, zorder=5)
        ax.set_xticks(range(len(conditions)))
        ax.set_xticklabels(conditions, rotation=15, ha="right",
                           fontsize=8 * SC)
        ax.set_ylabel(label, fontsize=10 * SC)
        ax.tick_params(axis="y", labelsize=9 * SC)
        ax.grid(axis="y", alpha=0.3)

        # Significance brackets stacked above the data
        y_lo, y_hi = ax.get_ylim()
        y_range = max(y_hi - y_lo, 1e-9)
        bracket_step = 0.11 * y_range
        for k, (a, b) in enumerate(pairs):
            if a not in cond_idx or b not in cond_idx:
                continue
            x1, x2 = sorted((cond_idx[a], cond_idx[b]))
            p = p_lookup(sheet, m, a, b, family_m)
            if p is None:
                continue
            txt, _ = _format_p(p)
            y_bar = y_hi + bracket_step * (k + 0.45)
            tick = bracket_step * 0.18
            ax.plot([x1, x1, x2, x2],
                    [y_bar - tick, y_bar, y_bar, y_bar - tick],
                    lw=0.9, color="black")
            ax.text((x1 + x2) / 2, y_bar + bracket_step * 0.05, txt,
                    ha="center", va="bottom", fontsize=7.5 * SC)
        ax.set_ylim(y_lo, y_hi + bracket_step * (len(pairs) + 0.6))

    fig.suptitle(f"{title_prefix}  ·  Šídák m={family_m}",
                 fontsize=11 * SC)
    plt.tight_layout(rect=[0, 0.01, 1, 0.97])
    fig.savefig(out_png, dpi=130)
    plt.close(fig)


def main():
    if not NEW_CSV.exists():
        print(f"Not found: {NEW_CSV}"); return 1
    df = load()
    print(f"Loaded {df.height} cells across {sorted(df['sheet'].unique().to_list())}")

    OUT_DIR.mkdir(parents=True, exist_ok=True)

    eval_df = pl.read_csv(EVAL_CSV) if EVAL_CSV.exists() else None
    p_lookup = make_p_lookup(eval_df) if eval_df is not None else lambda *a, **k: None

    for sheet in sorted(df["sheet"].unique().to_list()):
        sheet_df = df.filter(pl.col("sheet") == sheet)
        cfg = SHEET_CONFIG.get(sheet)
        if cfg is None:
            conditions = sorted(sheet_df["condition"].unique().to_list())
            pairs = []
            family_m = 1
        else:
            conditions = cfg["conditions"]
            pairs = cfg["pairs"]
            family_m = cfg["family_m"]
        slug_name = slug(sheet)
        title = f"{sheet} (n={sheet_df.height})"
        print(f"\nPlotting {sheet} (n={sheet_df.height}, conditions={conditions})")

        plot_profile(sheet_df, conditions, "y", f"{title} · Y-axis profile",
                     OUT_DIR / f"{slug_name}_y_profile.png")
        plot_profile(sheet_df, conditions, "wedge_r",
                     f"{title} · wedge-r profile",
                     OUT_DIR / f"{slug_name}_wedge_r_profile.png")
        plot_wedge_r_cdf(sheet_df, conditions, f"{title} · wedge-r CDF",
                         OUT_DIR / f"{slug_name}_wedge_r_cdf.png")
        plot_scalars(sheet_df, sheet, conditions, pairs, family_m, p_lookup,
                     f"{title} · per-cell scalars",
                     OUT_DIR / f"{slug_name}_scalars.png")

    print(f"\nFigures written to {OUT_DIR}/")
    for p in sorted(OUT_DIR.glob("*.png")):
        print(f"  {p.name}")


if __name__ == "__main__":
    main()
