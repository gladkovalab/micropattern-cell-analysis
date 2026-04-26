"""Render the MAPK9-siRNA scalar / wedge-r-profile / wedge-r-CDF figures
on the DENOISED dataset, with the same Šídák m=3 statistics overlay as the
canonical raw-data figures.

Reads:
  - replication/overnight_final_mapk9_denoised/combined_raw.csv
  - replication/overnight_final_out/mapk9_denoised_metadata.csv

Writes 3 figures into
  replication/overnight_final_mapk9_denoised/figures/:
    MAPK9_siRNA_scalars_denoised.png
    MAPK9_siRNA_wedge_r_profile_denoised.png
    MAPK9_siRNA_wedge_r_cdf_denoised.png

Side effect: writes the Šídák stats CSV to
  replication/overnight_final_mapk9_denoised/evaluation_summary_mapk9_denoised.csv
"""
from __future__ import annotations
import pathlib, sys
import polars as pl

REPO = pathlib.Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / "replication"))

from plot_final import (  # noqa: E402
    plot_scalars, plot_profile, plot_wedge_r_cdf,
    SHEET_CONFIG, EVAL_FAMILY_M, slug,
)
from evaluate_final import test_pair, KEEPER_SCALARS, SHEET_PAIRS  # noqa: E402

DENOISED_DIR = REPO / "replication" / "overnight_final_mapk9_denoised"
NEW_CSV = DENOISED_DIR / "combined_raw.csv"
META_CSV = REPO / "replication" / "overnight_final_out" / "mapk9_denoised_metadata.csv"
OUT_DIR = DENOISED_DIR / "figures"
EVAL_CSV = DENOISED_DIR / "evaluation_summary_mapk9_denoised.csv"

SHEET = "MAPK9 siRNA"


def main():
    if not NEW_CSV.exists():
        print(f"Not found: {NEW_CSV}\n(Pipeline still running?)")
        return 1
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    # Load + join metadata
    new = pl.read_csv(NEW_CSV)
    meta = pl.read_csv(META_CSV).select(["path","plate","well","sheet","condition"])
    df = new.join(meta, on="path", how="left").filter(
        pl.col("condition").is_not_null())
    print(f"Loaded {df.height} denoised MAPK9 cells, "
          f"conditions: {sorted(df['condition'].unique().to_list())}")

    cfg = SHEET_CONFIG[SHEET]
    pairs = cfg["pairs"]
    family_m = cfg["family_m"]
    conditions = cfg["conditions"]

    # Compute stats over the user-specified pair list (m=3)
    rows = []
    keeper = []
    for tmpl in KEEPER_SCALARS:
        for proj in ("zsum","maxip"):
            keeper.append(tmpl.format(p=proj))
    keeper = [m for m in keeper if m in df.columns]
    for pair in pairs:
        for m in keeper:
            r = test_pair(df, m, pair, family_m)
            if r is None:
                continue
            rows.append({"sheet": SHEET, "pair": f"{pair[0]} vs {pair[1]}",
                         "metric": m, **r})
    eval_df = pl.from_dicts(rows)
    eval_df.write_csv(EVAL_CSV)
    print(f"Wrote {EVAL_CSV} ({eval_df.height} rows)")

    # p_lookup that uses our just-computed eval table
    EVAL_M = family_m  # we just computed at the requested family

    def p_lookup(sheet, metric, cond_a, cond_b, family):
        for label in (f"{cond_a} vs {cond_b}", f"{cond_b} vs {cond_a}"):
            r = eval_df.filter((pl.col("sheet") == sheet) &
                               (pl.col("metric") == metric) &
                               (pl.col("pair") == label))
            if r.height > 0:
                p = r["p"][0]
                if family == EVAL_M:
                    return float(p) if p == p else None  # NaN check
                return float(1 - (1 - p) ** (family / EVAL_M))
        return None

    title_pre = f"{SHEET} (n={df.height}) DENOISED"
    plot_scalars(df, SHEET, conditions, pairs, family_m, p_lookup,
                 f"{title_pre} · per-cell scalars",
                 OUT_DIR / "MAPK9_siRNA_scalars_denoised.png")
    plot_profile(df, conditions, "wedge_r",
                 f"{title_pre} · wedge-r profile",
                 OUT_DIR / "MAPK9_siRNA_wedge_r_profile_denoised.png")
    plot_wedge_r_cdf(df, conditions, f"{title_pre} · wedge-r CDF",
                     OUT_DIR / "MAPK9_siRNA_wedge_r_cdf_denoised.png")

    print(f"\nFigures written to {OUT_DIR}/")
    for p in sorted(OUT_DIR.glob("*.png")):
        print(f"  {p.name}")


if __name__ == "__main__":
    main()
