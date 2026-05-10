"""Per-sheet split view of the wedge-r mito profile + nuclear-mask overlay.

For each of the 6 comparison sheets, generates a 1×N panel figure (one
panel per condition) showing:

  - solid coloured curve: 488 mitochondria mean ± SEM (% of wedge total)
  - dashed black curve:   binary nuclear mask radial distribution
                          (Otsu segmentation of 405 Z-sum, same as Mark's
                          pipeline; mask % per 1 µm shell)
  - shaded grey slabs:    [18, 33) µm centrosomal and [41, 56) µm peripheral
                          (the iso-centred bands)

Writes one PNG per sheet to analysis/figures_wedge_r_ks/.
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
import skimage
from scipy.ndimage import distance_transform_edt
import xarray as xr

REPO = pathlib.Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / "analysis"))
from plot_metrics import (  # noqa: E402
    SHEET_CONFIG, CONDITION_COLORS, load_template_matching, join_with_metadata,
)
import template_matching_bulk as tmb  # noqa: E402

INNER_BAND = (18, 33)
OUTER_BAND = (41, 56)
PITCH_UM = 0.065

SHEET_OUT = {
    "TRAK isoform (mito)":       "profiles_mito_with_nuclear.png",
    "TRAK isoform (peroxisome)": "profiles_peroxisome_with_nuclear.png",
    "TRAK isoform (60mer)":      "profiles_60mer_with_nuclear.png",
    "TRAK1 helix muts":          "profiles_trak1_helix_with_nuclear.png",
    "TRAK2 helix muts":          "profiles_trak2_helix_with_nuclear.png",
    "MAPK9 siRNA":               "profiles_mapk9_with_nuclear.png",
}


def proj_path(nd2_path: str) -> pathlib.Path:
    parts = pathlib.Path(nd2_path).parts
    idx = parts.index("patterned_data")
    return REPO / "projections" / pathlib.Path(*parts[idx + 1:]).with_suffix(".nc")


def nuclear_mask_profiles(nd2_path: str) -> tuple[np.ndarray, np.ndarray]:
    """Wedge-r profiles of:
      1. the binary nuclear mask (Otsu + largest CC of 405 Z-sum)
      2. the 5 µm dilation halo around it — exactly the extra pixels the
         perinuclear_5um metric includes beyond the bare nucleus.
    Returns (nuc_profile, halo_profile), both NaN if segmentation fails.
    """
    pp = proj_path(nd2_path)
    nan = np.full(60, np.nan)
    if not pp.exists():
        return nan, nan
    da = xr.open_dataset(pp)
    nuc = da.sel(C="405").to_array().squeeze().values.reshape((1024, 1024))
    thresh = skimage.filters.threshold_otsu(nuc)
    mask = nuc > thresh
    label = skimage.measure.label(mask)
    props = skimage.measure.regionprops(label)
    if not props:
        return nan, nan
    largest = int(np.argmax([p.area for p in props])) + 1
    nuc_mask = (label == largest)

    # perinuclear-5µm region (matches Mark's metric): pixels with EDT-to-
    # nucleus < 5 µm. Inside-the-nucleus EDT is 0, so this is nucleus +
    # 5 µm halo. Halo alone = perinuclear region minus nucleus.
    edt = distance_transform_edt(np.invert(nuc_mask))
    perinuc = edt < (5.0 / PITCH_UM)
    halo = perinuc & np.invert(nuc_mask)

    geom = tmb._get_wedge_geometry(nuc_mask.shape, PITCH_UM)
    nuc_prof = tmb.wedge_r_profile(nuc_mask.astype(np.float64), wedge_geom=geom)
    halo_prof = tmb.wedge_r_profile(halo.astype(np.float64), wedge_geom=geom)
    return nuc_prof, halo_prof


def mito_cols(df: pl.DataFrame) -> list[str]:
    cols = [f"wedge_r_{i:02d}_{i+1:02d}um_pct" for i in range(60)]
    return [c for c in cols if c in df.columns]


def per_condition(sheet_df: pl.DataFrame, mcols: list[str], cache: dict):
    out = {}
    for cond in sheet_df["condition"].unique().to_list():
        sub = sheet_df.filter(pl.col("condition") == cond)
        mito = sub.select(mcols).to_numpy().astype(float)
        nuc_rows, halo_rows = [], []
        for path in sub["path"].to_list():
            if path not in cache:
                cache[path] = nuclear_mask_profiles(path)
            np_, hp_ = cache[path]
            nuc_rows.append(np_); halo_rows.append(hp_)
        nuc = np.vstack(nuc_rows) if nuc_rows else np.zeros((0, len(mcols)))
        halo = np.vstack(halo_rows) if halo_rows else np.zeros((0, len(mcols)))
        valid = np.isfinite(nuc).all(axis=1) if nuc.size else np.array([], bool)
        nuc_v = nuc[valid] if nuc.size else nuc
        halo_v = halo[valid] if halo.size else halo
        n = sub.height
        sqrt_n = max(np.sqrt(n), 1.0)
        out[cond] = {
            "n": n,
            "n_nuc": nuc_v.shape[0],
            "mito_mean": np.nanmean(mito, axis=0),
            "mito_sem": np.nanstd(mito, axis=0, ddof=1) / sqrt_n,
            "nuc_mean": np.nanmean(nuc_v, axis=0) if nuc_v.size else None,
            "nuc_sem": (np.nanstd(nuc_v, axis=0, ddof=1) /
                        max(np.sqrt(nuc_v.shape[0]), 1.0)
                        if nuc_v.size else None),
            "halo_mean": np.nanmean(halo_v, axis=0) if halo_v.size else None,
            "halo_sem": (np.nanstd(halo_v, axis=0, ddof=1) /
                         max(np.sqrt(halo_v.shape[0]), 1.0)
                         if halo_v.size else None),
        }
    return out


def render_sheet(sheet: str, sheet_df: pl.DataFrame, mcols: list[str],
                 out_path: pathlib.Path, cache: dict):
    cfg = SHEET_CONFIG[sheet]
    cond_order = cfg["conditions"]
    color_map = {c: CONDITION_COLORS[i % len(CONDITION_COLORS)]
                 for i, c in enumerate(cond_order)}
    n_bins = len(mcols)
    centers = np.array([i + 0.5 for i in range(n_bins)])

    print(f"\n=== {sheet}: {sheet_df.height} cells, "
          f"{len(cond_order)} conditions ===")
    curves = per_condition(sheet_df, mcols, cache)

    n_panels = len(cond_order)
    width = max(5.5 * n_panels, 12)
    fig, axes = plt.subplots(1, n_panels, figsize=(width, 6),
                             sharey=True, sharex=True)
    if n_panels == 1:
        axes = [axes]

    # Determine y-limit so every condition shares the same scale.
    ymax = 0.0
    for cond in cond_order:
        c = curves.get(cond)
        if c is None:
            continue
        for key in ("mito_mean", "nuc_mean", "halo_mean"):
            mean = c.get(key)
            sem_key = key.replace("_mean", "_sem")
            sem = c.get(sem_key)
            if mean is not None and sem is not None:
                ymax = max(ymax, float(np.nanmax(mean + sem)))
    ymax *= 1.10

    for ax, cond in zip(axes, cond_order):
        if cond not in curves:
            ax.text(0.5, 0.5, f"no rows: {cond}", ha="center", va="center",
                    transform=ax.transAxes)
            continue
        c = curves[cond]
        col = color_map[cond]

        for lo, hi in (INNER_BAND, OUTER_BAND):
            ax.axvspan(lo, hi, color="0.85", zorder=0, linewidth=0)

        ax.plot(centers, c["mito_mean"], color=col, lw=2.0,
                label=f"488 mito (n={c['n']})")
        ax.fill_between(centers,
                        c["mito_mean"] - c["mito_sem"],
                        c["mito_mean"] + c["mito_sem"],
                        color=col, alpha=0.20, linewidth=0)

        if c["nuc_mean"] is not None:
            ax.plot(centers, c["nuc_mean"], color="#1a1a1a", lw=1.8,
                    linestyle="--", label=f"405 nuclear mask (n={c['n_nuc']})")
            ax.fill_between(centers,
                            c["nuc_mean"] - c["nuc_sem"],
                            c["nuc_mean"] + c["nuc_sem"],
                            color="#1a1a1a", alpha=0.12, linewidth=0)
        if c["halo_mean"] is not None:
            # 5 µm dilation halo around the nucleus — extra pixels the
            # perinuclear_5um metric includes beyond the bare nucleus.
            ax.plot(centers, c["halo_mean"], color="#a050c0", lw=1.6,
                    linestyle=":", label="5 µm perinuclear halo")
            ax.fill_between(centers,
                            c["halo_mean"] - c["halo_sem"],
                            c["halo_mean"] + c["halo_sem"],
                            color="#a050c0", alpha=0.15, linewidth=0)

        ax.set_title(cond, fontsize=14)
        ax.set_xlabel("wedge-r (µm from apex)", fontsize=11)
        ax.legend(fontsize=10, loc="upper right")
        ax.grid(alpha=0.3)
        ax.set_xlim(0, n_bins)
        ax.set_ylim(0, ymax)

    axes[0].set_ylabel("mean per-bin % (of wedge total / nuclear-mask area)",
                       fontsize=11)
    fig.suptitle(
        f"{sheet} — 488 mito (solid), 405 nuclear mask (dashed), "
        f"5 µm perinuclear halo (dotted; what the perinuclear_5um metric adds)\n"
        f"shaded: [{INNER_BAND[0]}, {INNER_BAND[1]}) µm centrosomal  ·  "
        f"[{OUTER_BAND[0]}, {OUTER_BAND[1]}) µm peripheral   "
        f"(all curves: mean ± SEM)",
        fontsize=11)
    plt.tight_layout(rect=[0, 0, 1, 0.93])
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=140)
    plt.close(fig)
    print(f"  wrote {out_path}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out-dir",
                    default="analysis/figures_wedge_r_ks")
    args = ap.parse_args()

    df = load_template_matching(pathlib.Path(
        "analysis/wedge_r_ks_out_all_denoised/by_well"))
    df = join_with_metadata(df, REPO / "config/Comparisons_table_v3.xlsx")
    mcols = mito_cols(df)

    cache: dict[str, np.ndarray] = {}
    out_dir = pathlib.Path(args.out_dir)
    for sheet, fname in SHEET_OUT.items():
        sheet_df = df.filter(pl.col("sheet") == sheet)
        if sheet_df.height == 0:
            print(f"\n=== {sheet}: no rows, skipping ===")
            continue
        render_sheet(sheet, sheet_df, mcols, out_dir / fname, cache)

    print(f"\nunique cells processed for nuclear mask: {len(cache)}")


if __name__ == "__main__":
    main()
