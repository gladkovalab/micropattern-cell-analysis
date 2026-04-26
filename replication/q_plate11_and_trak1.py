"""Three probes to answer the user's questions:

1. Pull per-cell Y-Gini / entropy / SD values for the specific plate-11 cells
   flagged by the user (peripheral-centrosome cells): noTRAK 5/11/13 and
   TRAK2 9/13. Do they stick out from the rest of their condition?

2. Rank metrics by p-value for the no-TRAK vs TRAK1 comparison to see what
   (other than Y-Gini) captures the subtle less-clustered phenotype.

3. Cross-check those metrics against no-TRAK vs TRAK2 — metrics that separate
   BOTH pairs are the most robust clustering indicators.
"""
from __future__ import annotations

import pathlib
import sys
import warnings

import numpy as np
import polars as pl

REPO = pathlib.Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / "replication"))

from replicate_stats import ConditionData, nested_oneway_anova, sidak_pairwise  # noqa: E402
from evaluate_metrics import cohens_d  # noqa: E402

warnings.filterwarnings("ignore")

NEW_CSV = REPO / "replication" / "overnight_fig4b_out" / "combined_raw.csv"
OLD_CSV = REPO / "replication" / "overnight_out" / "combined.csv"


def load():
    new = pl.read_csv(NEW_CSV)
    old = pl.read_csv(OLD_CSV).filter(pl.col("sheet") == "TRAK isoform (mito)").select(
        ["path", "plate", "well", "sheet", "condition"])
    return new.join(old, on="path", how="left").filter(pl.col("condition").is_not_null())


def collect(df, metric):
    sub = df.filter(pl.col(metric).is_not_null() & pl.col(metric).is_not_nan())
    conds = []
    for cn in sorted(sub["condition"].unique().to_list()):
        g = sub.filter(pl.col("condition") == cn)
        pc = {}
        for plate, grp in g.group_by("plate"):
            key = plate[0] if isinstance(plate, tuple) else plate
            pc[key] = grp[metric].to_numpy().astype(float)
        conds.append(ConditionData(name=cn, plate_cells=pc))
    return conds


def pair_test(df, metric, pair, family=3):
    conds = collect(df, metric)
    n2i = {c.name: i for i, c in enumerate(conds)}
    if pair[0] not in n2i or pair[1] not in n2i:
        return None
    i, j = n2i[pair[0]], n2i[pair[1]]
    try:
        a = nested_oneway_anova(conds)
        r = sidak_pairwise(conds, a, pairs=[(i, j)])[0]
        p = 1 - (1 - r["p_sidak"]) ** family if np.isfinite(r["p_sidak"]) else np.nan
    except Exception:
        return None
    return {"d": cohens_d(conds[i].all_cells, conds[j].all_cells), "p": p}


def main():
    df = load()

    # -------- Q1: plate 11 flagged cells on Y-metrics --------
    print("=" * 80)
    print("Q1: plate 11 flagged cells (peripheral centrosome)")
    print("=" * 80)
    flagged = [
        ("no TRAK", "D06", "Cell5"), ("no TRAK", "D06", "Cell11"), ("no TRAK", "D06", "Cell13"),
        ("TRAK2",   "F05", "Cell9"), ("TRAK2",   "F05", "Cell13"),
    ]
    metrics_q1 = ["zsum_y_gini", "zsum_y_entropy", "zsum_y_sd_u",
                  "maxip_y_gini", "maxip_y_entropy", "maxip_y_sd_u"]

    # Get plate-11 population stats per condition to provide reference
    p11 = df.filter(pl.col("plate") == "250731_patterned_plate_11_good")
    other_plates = df.filter(pl.col("plate") != "250731_patterned_plate_11_good")

    print(f"\n{'Cell':<20} {'Cond':<10} {'Plate':<5} " + " ".join(f"{m.split('_',1)[1]:>14}" for m in metrics_q1))
    print("-" * (40 + 16 * len(metrics_q1)))

    def fmt(r, m):
        v = r[m] if m in r else float("nan")
        return f"{v:>14.3f}" if v is not None and np.isfinite(v) else f"{'nan':>14}"

    # Flagged cells
    for cond, well_pfx, cell in flagged:
        sub = df.filter(
            (pl.col("plate") == "250731_patterned_plate_11_good") &
            (pl.col("well") == well_pfx) &
            (pl.col("condition") == cond) &
            pl.col("path").str.contains(f"/{cell}.nd2")
        )
        if sub.height == 0:
            # try lowercase
            sub = df.filter(
                (pl.col("plate") == "250731_patterned_plate_11_good") &
                (pl.col("well") == well_pfx) &
                (pl.col("condition") == cond) &
                pl.col("path").str.contains(f"/{cell.lower()}.nd2")
            )
        if sub.height == 0:
            print(f"  {cell:<20} {cond:<10} p11   (NOT FOUND in CSV)")
            continue
        r = sub.row(0, named=True)
        print(f"  {cell:<20} {cond:<10} p11 " + " ".join(fmt(r, m) for m in metrics_q1))

    # Plate-11 condition means for reference
    print("\nPlate-11 condition means ± SD (all cells in that condition on plate 11):")
    print(f"  {'Cond':<30} " + " ".join(f"{m.split('_',1)[1]:>14}" for m in metrics_q1))
    for cond in ("no TRAK", "TRAK2"):
        sub = p11.filter(pl.col("condition") == cond)
        if sub.height == 0:
            continue
        line = f"  {cond + f' (n={sub.height})':<30} "
        for m in metrics_q1:
            v = sub[m].to_numpy()
            v = v[~np.isnan(v)]
            line += f"{v.mean():>7.3f}±{v.std():<6.3f}" if v.size else f"{'nan':>14}"
        print(line)

    print("\nOther-plates (3,9) condition means for comparison:")
    for cond in ("no TRAK", "TRAK2"):
        sub = other_plates.filter(pl.col("condition") == cond)
        if sub.height == 0:
            continue
        line = f"  {cond + f' (n={sub.height})':<30} "
        for m in metrics_q1:
            v = sub[m].to_numpy()
            v = v[~np.isnan(v)]
            line += f"{v.mean():>7.3f}±{v.std():<6.3f}" if v.size else f"{'nan':>14}"
        print(line)

    # -------- Q2: top metrics for no-TRAK vs TRAK1 --------
    print("\n\n" + "=" * 80)
    print("Q2: top metrics separating no-TRAK vs TRAK1 (biology: TRAK1 slightly less clustered)")
    print("=" * 80)
    META = {"path", "plate", "well", "sheet", "condition",
            "template_matching_score", "lateral_pixel_pitch_um",
            "zsum_bg_threshold", "maxip_bg_threshold",
            "nuc_area_um2", "nuc_solidity", "nuc_eccentricity", "nuc_euler_number",
            "nuc_n_components", "nuc_largest_area_frac",
            "zsum_total_signal", "maxip_total_signal"}
    metric_cols = [c for c in df.columns if c not in META and
                   df[c].dtype in (pl.Float64, pl.Float32, pl.Int64)]

    rows_1 = []
    rows_2 = []
    for m in metric_cols:
        r1 = pair_test(df, m, ("no TRAK", "TRAK1"), family=3)
        r2 = pair_test(df, m, ("no TRAK", "TRAK2"), family=3)
        if r1 is not None and np.isfinite(r1["p"]):
            rows_1.append({"metric": m, "d_vs_T1": r1["d"], "p_vs_T1": r1["p"],
                           "d_vs_T2": r2["d"] if r2 else np.nan,
                           "p_vs_T2": r2["p"] if r2 else np.nan})

    s1 = pl.from_dicts(rows_1).sort("p_vs_T1")
    print(f"\nTop 20 metrics for no TRAK vs TRAK1 (Šídák m=3), plus their no-TRAK-vs-TRAK2 p:\n")
    print(f"{'Metric':<50}  {'d(T1)':>6} {'p(T1)':>8}  {'d(T2)':>6} {'p(T2)':>8}  same dir?")
    print("-" * 100)
    for r in s1.head(25).iter_rows(named=True):
        sig1 = "***" if r["p_vs_T1"] < 0.001 else "**" if r["p_vs_T1"] < 0.01 else "*" if r["p_vs_T1"] < 0.05 else "ns"
        sig2 = "***" if r["p_vs_T2"] < 0.001 else "**" if r["p_vs_T2"] < 0.01 else "*" if r["p_vs_T2"] < 0.05 else "ns"
        dirmatch = "Y" if np.sign(r["d_vs_T1"]) == np.sign(r["d_vs_T2"]) else "FLIP"
        print(f"{r['metric']:<50}  {r['d_vs_T1']:+6.2f} {r['p_vs_T1']:8.4f}{sig1:<3}  "
              f"{r['d_vs_T2']:+6.2f} {r['p_vs_T2']:8.4f}{sig2:<3}  {dirmatch}")

    # -------- Q3: metrics significant on BOTH pairs (robust "clustering" indicators) --------
    print("\n\n" + "=" * 80)
    print("Q3: metrics significant on BOTH no-TRAK→TRAK1 AND no-TRAK→TRAK2 (same direction)")
    print("=" * 80)
    dual = s1.filter(
        (pl.col("p_vs_T1") < 0.05) & (pl.col("p_vs_T2") < 0.05) &
        (pl.col("d_vs_T1").sign() == pl.col("d_vs_T2").sign())
    ).sort("p_vs_T1")
    print(f"\n{dual.height} metrics:\n")
    print(f"{'Metric':<50}  {'d(T1)':>6} {'p(T1)':>8}  {'d(T2)':>6} {'p(T2)':>8}")
    print("-" * 90)
    for r in dual.iter_rows(named=True):
        print(f"{r['metric']:<50}  {r['d_vs_T1']:+6.2f} {r['p_vs_T1']:8.4f}  "
              f"{r['d_vs_T2']:+6.2f} {r['p_vs_T2']:8.4f}")

    # borderline: p<0.1 on both, same direction
    near = s1.filter(
        (pl.col("p_vs_T1") < 0.10) & (pl.col("p_vs_T2") < 0.10) &
        (pl.col("d_vs_T1").sign() == pl.col("d_vs_T2").sign())
    ).sort("p_vs_T1")
    print(f"\n\nBorderline both-p<0.10 same-direction ({near.height} metrics):")
    print(f"{'Metric':<50}  {'d(T1)':>6} {'p(T1)':>8}  {'d(T2)':>6} {'p(T2)':>8}")
    print("-" * 90)
    for r in near.head(30).iter_rows(named=True):
        print(f"{r['metric']:<50}  {r['d_vs_T1']:+6.2f} {r['p_vs_T1']:8.4f}  "
              f"{r['d_vs_T2']:+6.2f} {r['p_vs_T2']:8.4f}")


if __name__ == "__main__":
    main()
