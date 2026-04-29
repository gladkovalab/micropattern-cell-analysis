"""Three-panel figure (no TRAK / TRAK1 / TRAK2) of the 60mer wedge-r profile,
with the per-condition nuclear (405) wedge-r profile overlaid as a dashed
line on each panel.

Mito (488) profile comes from the existing wedge_r_NN_NN+1um_pct columns
in the per-well CSVs. Nuclear (405) profile is recomputed offline from the
saved Z-sum projection netCDFs in `projections/...` using the same wedge
geometry (no pipeline rerun required).
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
import xarray as xr

REPO = pathlib.Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / "replication"))
from plot_metrics import (  # noqa: E402
    CONDITION_COLORS, load_template_matching, join_with_metadata,
)
import template_matching_bulk as tmb  # noqa: E402

INNER_BAND = (18, 33)
OUTER_BAND = (41, 56)
SHEET = "TRAK isoform (60mer)"
PITCH_UM = 0.065  # constant across the dataset (verified earlier)


def proj_path(nd2_path: str) -> pathlib.Path:
    parts = pathlib.Path(nd2_path).parts
    idx = parts.index("patterned_data")
    rel = pathlib.Path(*parts[idx + 1:]).with_suffix(".nc")
    return REPO / "projections" / rel


def nuclear_profile(nd2_path: str) -> np.ndarray:
    """Wedge-r profile (% per 1µm bin) of the *binary nuclear mask*.

    The Z-sum 405 channel is dominated by integrated background (the
    per-bin % then just tracks wedge area, peaking at the rim). We use
    the same Otsu-thresholded + largest-connected-component segmentation
    Mark's pipeline uses for nuclear-edge metrics, so what we plot is
    "fraction of nuclear *area* in this 1 µm shell" — a clean readout
    of where the nucleus sits along the wedge-r axis.
    """
    import skimage
    nc = proj_path(nd2_path)
    da = xr.open_dataset(nc)
    nuc = da.sel(C="405").to_array().squeeze().values.reshape((1024, 1024))
    thresh = skimage.filters.threshold_otsu(nuc)
    mask = nuc > thresh
    label = skimage.measure.label(mask)
    props = skimage.measure.regionprops(label)
    if not props:
        return np.full(60, np.nan)
    largest = int(np.argmax([p.area for p in props])) + 1
    nuc_mask = (label == largest).astype(np.float64)
    geom = tmb._get_wedge_geometry(nuc_mask.shape, PITCH_UM)
    return tmb.wedge_r_profile(nuc_mask, wedge_geom=geom)


def mito_profile_columns(df: pl.DataFrame) -> list[str]:
    cols = [f"wedge_r_{i:02d}_{i+1:02d}um_pct" for i in range(60)]
    return [c for c in cols if c in df.columns]


def per_condition_curves(sheet_df: pl.DataFrame, mito_cols: list[str]):
    """Returns dict: condition -> (mito_mean, mito_sem, nuc_mean, nuc_sem, n)."""
    out = {}
    for cond in sheet_df["condition"].unique().to_list():
        sub = sheet_df.filter(pl.col("condition") == cond)
        mito = sub.select(mito_cols).to_numpy().astype(float)
        nuc_rows = []
        for path in sub["path"].to_list():
            try:
                nuc_rows.append(nuclear_profile(path))
            except Exception as e:
                print(f"  WARN: skip {path}: {e}")
        nuc = np.vstack(nuc_rows) if nuc_rows else np.zeros((0, len(mito_cols)))
        n = sub.height
        sqrt_n = max(np.sqrt(n), 1.0)
        out[cond] = {
            "n": n,
            "mito_mean": np.nanmean(mito, axis=0),
            "mito_sem": np.nanstd(mito, axis=0, ddof=1) / sqrt_n,
            "nuc_mean": np.nanmean(nuc, axis=0) if nuc.size else None,
            "nuc_sem": np.nanstd(nuc, axis=0, ddof=1) / sqrt_n if nuc.size else None,
        }
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out",
                    default="replication/figures_wedge_r_ks/profiles_60mer_with_nuclear.png")
    args = ap.parse_args()

    df = load_template_matching(pathlib.Path(
        "replication/wedge_r_ks_out_all_denoised/by_well"))
    df = join_with_metadata(df, REPO / "config/Comparisons_table_v3.xlsx")
    sheet_df = df.filter(pl.col("sheet") == SHEET)
    print(f"60mer cells: {sheet_df.height}")
    print(f"computing nuclear profiles ({sheet_df.height} cells)...")

    mito_cols = mito_profile_columns(sheet_df)
    n_bins = len(mito_cols)
    centers = np.array([i + 0.5 for i in range(n_bins)])
    curves = per_condition_curves(sheet_df, mito_cols)

    # 3-panel figure, one per condition
    cond_order = ["no TRAK", "TRAK1", "TRAK2"]
    color_map = {c: CONDITION_COLORS[i] for i, c in enumerate(cond_order)}

    fig, axes = plt.subplots(1, 3, figsize=(18, 6), sharey=True, sharex=True)
    for ax, cond in zip(axes, cond_order):
        if cond not in curves:
            ax.text(0.5, 0.5, f"no {cond}", ha="center", va="center",
                    transform=ax.transAxes)
            continue
        c = curves[cond]
        col = color_map[cond]

        # Grey slabs underneath
        for lo, hi in (INNER_BAND, OUTER_BAND):
            ax.axvspan(lo, hi, color="0.85", zorder=0, linewidth=0)

        # 488 mito profile (solid + SEM band)
        ax.plot(centers, c["mito_mean"], color=col, lw=2.0,
                label=f"488 mitochondria (n={c['n']})")
        ax.fill_between(centers, c["mito_mean"] - c["mito_sem"],
                        c["mito_mean"] + c["mito_sem"],
                        color=col, alpha=0.20, linewidth=0)

        # 405 nuclear profile (dashed, darker)
        if c["nuc_mean"] is not None:
            ax.plot(centers, c["nuc_mean"], color="#1a1a1a", lw=1.8,
                    linestyle="--", label=f"405 nuclear mask (n={c['n']})")
            ax.fill_between(centers, c["nuc_mean"] - c["nuc_sem"],
                            c["nuc_mean"] + c["nuc_sem"],
                            color="#1a1a1a", alpha=0.12, linewidth=0)

        ax.set_title(f"{cond}", fontsize=14)
        ax.set_xlabel("wedge-r (µm from apex)", fontsize=12)
        ax.legend(fontsize=10, loc="upper right")
        ax.grid(alpha=0.3)
        ax.set_xlim(0, n_bins)

    axes[0].set_ylabel("mean intensity per bin (% of wedge total)", fontsize=12)
    fig.suptitle(
        f"60mer wedge-r profiles per condition — solid = 488 mitochondria (MaxIP), "
        f"dashed = 405 nuclear mask (Otsu) (mean ± SEM)\n"
        f"shaded slabs: [{INNER_BAND[0]}, {INNER_BAND[1]}) µm centrosomal  ·  "
        f"[{OUTER_BAND[0]}, {OUTER_BAND[1]}) µm peripheral",
        fontsize=13)
    plt.tight_layout(rect=[0, 0, 1, 0.93])
    out = pathlib.Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=140)
    plt.close(fig)
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
