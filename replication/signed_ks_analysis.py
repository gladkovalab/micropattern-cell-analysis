"""One-off: explore signed KS for the prior mito-only run.

For each cell:
  F(r) = cumsum(wedge_r_NN_NN+1um_pct) / 100
  G_uniform(r) = (r/R)^2  (analytical area-uniform — matches plot_metrics.py)
  G_60mer(r)   = _REF_CDF_60MER_NOTRAK from template_matching_bulk

  D+ = sup(F - G)+
  D- = sup(G - F)+
  KS_signed = D+ - D-

Convention: positive KS_signed means F rises faster than reference at small r
i.e. cell is more concentrated near the wedge apex than the reference
(perinuclear-leaning). Negative means peripheral-leaning.
"""
from __future__ import annotations

import pathlib
import re
import sys

import fastexcel
import numpy as np
import polars as pl

REPO = pathlib.Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))
from template_matching_bulk import _REF_CDF_60MER_NOTRAK  # noqa: E402


def wedge_r_columns(df: pl.DataFrame) -> list[str]:
    pat = re.compile(r"^wedge_r_(\d{2})_(\d{2})um_pct$")
    matches = []
    for c in df.columns:
        m = pat.match(c)
        if m:
            matches.append((int(m.group(1)), c))
    matches.sort()
    return [c for _, c in matches]


def signed_ks(F: np.ndarray, G: np.ndarray) -> tuple[float, float, float]:
    """Return (D_plus, D_minus, KS_signed) for one CDF F against reference G."""
    diff = F - G
    d_plus = float(np.max(np.maximum(diff, 0.0)))
    d_minus = float(np.max(np.maximum(-diff, 0.0)))
    return d_plus, d_minus, d_plus - d_minus


def load_metadata(xlsx: pathlib.Path) -> pl.DataFrame:
    fe = fastexcel.read_excel(str(xlsx))
    rows = []
    for sheet_name in fe.sheet_names:
        df = fe.load_sheet_by_name(sheet_name).to_polars()
        plate_col = df.columns[0]
        for record in df.iter_rows(named=True):
            plate = record[plate_col]
            if not plate:
                continue
            for cond in df.columns[1:]:
                well = record[cond]
                if well:
                    rows.append({"plate": plate, "well": well,
                                 "sheet": sheet_name, "condition": cond})
    return pl.from_dicts(rows)


def main():
    csv_path = REPO / "replication/wedge_r_ks_out_mito_denoised/combined.csv"
    df = pl.read_csv(csv_path)

    # Join to condition metadata via (plate, well) extracted from path.
    plate_re = r"patterned_data/([^/]+)/"
    well_re = r"patterned_data/[^/]+/([A-Z]\d+)_"
    df = df.with_columns([
        pl.col("path").str.extract(plate_re, 1).alias("plate"),
        pl.col("path").str.extract(well_re, 1).alias("well"),
    ])
    meta = load_metadata(REPO / "config/Comparisons_table_v3.xlsx")
    meta = meta.unique(subset=["plate", "well", "sheet"]).filter(
        pl.col("sheet") == "TRAK isoform (mito)")
    df = df.join(meta, on=["plate", "well"], how="left")
    df = df.filter(pl.col("condition").is_not_null())
    print(f"loaded {df.height} cells from sheet 'TRAK isoform (mito)'")

    cols = wedge_r_columns(df)
    n_bins = len(cols)
    profile = df.select(cols).to_numpy()  # (n_cells, n_bins)
    # Renormalize per-cell to match wedge_r_cdf in the pipeline. The prior
    # mito run was generated with the old normalization (sums to <=100), so
    # dividing by literal 100 would bias F. cumsum / sum reaches 1.0 at the
    # last bin regardless of whether the per-bin pct sums to 100 or 99.7.
    profile = np.where(np.isnan(profile), 0.0, profile)
    row_sum = profile.sum(axis=1, keepdims=True)
    F = np.cumsum(profile, axis=1) / np.where(row_sum > 0, row_sum, 1.0)
    G_uni = (np.arange(1, n_bins + 1) ** 2) / (n_bins ** 2)
    G_60 = np.asarray(_REF_CDF_60MER_NOTRAK)

    rows = []
    for i in range(F.shape[0]):
        dpu, dmu, ksu = signed_ks(F[i], G_uni)
        dp6, dm6, ks6 = signed_ks(F[i], G_60)
        rows.append({
            "path": df["path"][i],
            "condition": df["condition"][i],
            "plate": df["plate"][i],
            "ks_uniform_unsigned": float(df["wedge_r_ks_vs_uniform"][i]),
            "ks_signed_vs_uniform": ksu,
            "d_plus_uniform": dpu,
            "d_minus_uniform": dmu,
            "ks_60mer_unsigned": float(df["wedge_r_ks_vs_60merNoTRAK"][i]),
            "ks_signed_vs_60mer": ks6,
            "d_plus_60mer": dp6,
            "d_minus_60mer": dm6,
        })
    out = pl.from_dicts(rows)

    # Sanity check: signed magnitude vs original unsigned (should be very close)
    delta = (out["ks_uniform_unsigned"] -
             np.maximum(out["d_plus_uniform"].to_numpy(),
                        out["d_minus_uniform"].to_numpy())).abs()
    print(f"\nsanity: max |unsigned - max(D+,D-)| vs uniform = {float(delta.max()):.6f} "
          f"(should be ~0 — confirms F/G computation matches the pipeline)")

    print("\n=== sign distribution: KS_signed vs area-uniform ===")
    print(f"  positive (perinuclear-leaning): {(out['ks_signed_vs_uniform'] > 0).sum()}")
    print(f"  negative (peripheral-leaning):  {(out['ks_signed_vs_uniform'] < 0).sum()}")
    print(f"  exactly zero:                   {(out['ks_signed_vs_uniform'] == 0).sum()}")

    print("\n=== sign distribution: KS_signed vs 60mer no-TRAK ===")
    print(f"  positive: {(out['ks_signed_vs_60mer'] > 0).sum()}")
    print(f"  negative: {(out['ks_signed_vs_60mer'] < 0).sum()}")
    print(f"  zero:     {(out['ks_signed_vs_60mer'] == 0).sum()}")

    print("\n=== per-condition breakdown ===")
    summary = out.group_by("condition").agg([
        pl.len().alias("n"),
        pl.col("ks_signed_vs_uniform").mean().alias("mean_signed_vs_uniform"),
        pl.col("ks_signed_vs_uniform").min().alias("min_signed_vs_uniform"),
        pl.col("ks_signed_vs_uniform").max().alias("max_signed_vs_uniform"),
        (pl.col("ks_signed_vs_uniform") > 0).sum().alias("n_pos_uniform"),
        (pl.col("ks_signed_vs_uniform") < 0).sum().alias("n_neg_uniform"),
        pl.col("ks_signed_vs_60mer").mean().alias("mean_signed_vs_60mer"),
        pl.col("ks_signed_vs_60mer").min().alias("min_signed_vs_60mer"),
        pl.col("ks_signed_vs_60mer").max().alias("max_signed_vs_60mer"),
        (pl.col("ks_signed_vs_60mer") > 0).sum().alias("n_pos_60mer"),
        (pl.col("ks_signed_vs_60mer") < 0).sum().alias("n_neg_60mer"),
    ]).sort("condition")
    with pl.Config(tbl_cols=-1, tbl_width_chars=200, set_fmt_float="full"):
        print(summary)


if __name__ == "__main__":
    main()
