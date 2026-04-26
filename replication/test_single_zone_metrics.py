"""Test single-zone perinuclear-clustering metrics on Fig 4B (TRAK isoform mito).

The reviewer's test case: 'no TRAK' vs 'TRAK2'. Mark's current perinuclear and
peripheral 5 µm metrics are ns (peri p=0.41). Can a single-zone metric that uses
the full intensity distribution (no peri/nuc scaling) catch the clustering that
is visible by eye?

Metrics tested, in two variants (zsum and maxip):
  - Mark's baselines: perinuclear_5um_pct, peripheral_5um_pct
  - Intensity-weighted mean / median / Q90 distance to nucleus
  - Perinuclear concentration = radial_0_2 + radial_2_5 (equals Mark's when 5 µm)
  - Radial Gini (Gini over 5 radial bin fractions)
  - Radial Shannon entropy (base-e; lower = more clustered)
  - Radial CoV (stddev / mean of 5 bin fractions)
  - 2nd-moment radius estimate (sqrt of variance of r over radial bins)
  - CoM offset from nucleus CoM (already directional, included for completeness)
  - CoM offset from pattern CoM

Reports Cohen's d and Šídák-corrected p-values (family=3 to match Fig 4B).
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

COMBINED = REPO / "replication" / "overnight_out" / "combined.csv"
OUT_DIR = REPO / "replication" / "single_zone_test_out"

SHEET = "TRAK isoform (mito)"
PAIRS = [("no TRAK", "TRAK1"), ("no TRAK", "TRAK2"), ("TRAK1", "TRAK2")]
FAMILY = 3

RADIAL_BIN_CENTERS_UM = np.array([1.0, 3.5, 7.5, 12.5, 20.0])  # approximate bin centers
RADIAL_BIN_EDGES_UM = np.array([0, 2, 5, 10, 15, 25])  # last bin is open-ended; 25 as cap estimate


def add_radial_derived(df: pl.DataFrame, proj: str, mask: str) -> pl.DataFrame:
    """Add radial-profile-derived scalar metrics for one (projection, mask)."""
    prefix = f"{proj}_{mask}"
    cols = [f"{prefix}_radial_0_2um_pct",
            f"{prefix}_radial_2_5um_pct",
            f"{prefix}_radial_5_10um_pct",
            f"{prefix}_radial_10_15um_pct",
            f"{prefix}_radial_ge15um_pct"]
    if not all(c in df.columns for c in cols):
        return df

    # Pull as numpy for per-row computation
    bins = df.select(cols).to_numpy()  # shape (N, 5), in percent
    fracs = bins / 100.0  # shape (N, 5), sums ≈ 1 per row

    # --- radial Gini ---
    # Gini of a length-k discrete distribution: same formula as intensity Gini.
    def _row_gini(p):
        p = np.sort(p)
        n = p.size
        s = p.sum()
        if s <= 0:
            return np.nan
        idx = np.arange(1, n + 1)
        return (2 * (idx * p).sum() - (n + 1) * s) / (n * s)
    gini = np.array([_row_gini(r) for r in fracs])

    # --- radial entropy (nats) ---
    # High = uniform (log 5 ≈ 1.609), low = concentrated
    with np.errstate(divide="ignore", invalid="ignore"):
        logp = np.where(fracs > 0, np.log(fracs), 0.0)
    entropy = -(fracs * logp).sum(axis=1)

    # --- radial CoV ---
    mean_f = fracs.mean(axis=1)
    std_f = fracs.std(axis=1)
    cov = np.where(mean_f > 0, std_f / mean_f, np.nan)

    # --- perinuclear concentration (0-5 µm) ---
    perinuc_0_5 = bins[:, 0] + bins[:, 1]  # in percent

    # --- 2nd-moment radius estimate from radial bins ---
    # E[r] = Σ p_i * c_i; Var[r] = Σ p_i * (c_i - E[r])^2
    centers = RADIAL_BIN_CENTERS_UM[None, :]  # (1, 5)
    mean_r = (fracs * centers).sum(axis=1)
    var_r = (fracs * (centers - mean_r[:, None]) ** 2).sum(axis=1)
    sd_r = np.sqrt(np.clip(var_r, 0, None))

    return df.with_columns([
        pl.Series(f"{prefix}_radial_gini", gini),
        pl.Series(f"{prefix}_radial_entropy", entropy),
        pl.Series(f"{prefix}_radial_cov", cov),
        pl.Series(f"{prefix}_perinuc_0_5um_pct", perinuc_0_5),
        pl.Series(f"{prefix}_radial_mean_r_um", mean_r),
        pl.Series(f"{prefix}_radial_sd_r_um", sd_r),
    ])


def collect(df: pl.DataFrame, metric: str) -> list[ConditionData]:
    sub = df.filter(pl.col(metric).is_not_null() & pl.col(metric).is_not_nan())
    conds: list[ConditionData] = []
    for cond_name in sorted(sub["condition"].unique().to_list()):
        g = sub.filter(pl.col("condition") == cond_name)
        plate_cells: dict[str, np.ndarray] = {}
        for plate, grp in g.group_by("plate"):
            key = plate[0] if isinstance(plate, tuple) else plate
            plate_cells[key] = grp[metric].to_numpy().astype(float)
        conds.append(ConditionData(name=cond_name, plate_cells=plate_cells))
    return conds


def test_metric(df: pl.DataFrame, metric: str) -> dict | None:
    conds = collect(df, metric)
    if len(conds) < 2:
        return None
    name_to_idx = {c.name: i for i, c in enumerate(conds)}
    a = nested_oneway_anova(conds)
    out: dict = {"metric": metric}
    for (p, q) in PAIRS:
        if p not in name_to_idx or q not in name_to_idx:
            continue
        i, j = name_to_idx[p], name_to_idx[q]
        r = sidak_pairwise(conds, a, pairs=[(i, j)])[0]
        p_raw = r["p_sidak"]  # m=1 at this call
        p_sidak = 1 - (1 - p_raw) ** FAMILY if np.isfinite(p_raw) else np.nan
        d = cohens_d(conds[i].all_cells, conds[j].all_cells)
        out[f"{p} vs {q} d"] = d
        out[f"{p} vs {q} p"] = p_sidak
        out[f"{p} vs {q} n_i"] = int(conds[i].all_cells.size)
        out[f"{p} vs {q} n_j"] = int(conds[j].all_cells.size)
        out[f"{p} vs {q} mean_i"] = float(conds[i].all_cells.mean())
        out[f"{p} vs {q} mean_j"] = float(conds[j].all_cells.mean())
    return out


def main():
    df = pl.read_csv(COMBINED)
    df = df.filter(pl.col("sheet") == SHEET)
    print(f"{df.height} cells on sheet '{SHEET}'")
    print("Conditions:", sorted(df["condition"].unique().to_list()))
    print("Plates present:", sorted(df["plate"].unique().to_list()))
    print()

    # Add derived radial metrics for each (proj, mask)
    for proj in ("zsum", "maxip"):
        for mask in ("crop", "pattern"):
            df = add_radial_derived(df, proj, mask)

    # Build the candidate list. Focus zsum + crop (default pitch), then
    # maxip + crop, plus pattern-mask variants as sensitivity checks.
    def variants(leaf_names):
        for proj in ("zsum", "maxip"):
            for mask in ("crop",):  # crop only for first pass; add pattern later
                for leaf in leaf_names:
                    yield f"{proj}_{mask}_{leaf}"

    candidates = list(variants([
        "perinuclear_5um_pct",  # Mark baseline (perinuc zone)
        "peripheral_5um_pct",   # Mark baseline (peri zone)
        "perinuc_0_5um_pct",    # = 0-2 + 2-5; should match perinuclear_5um_pct up to definitional diff
        "mean_dist_to_nucleus_um",
        "median_dist_to_nucleus_um",
        "q90_dist_to_nucleus_um",
        "com_offset_um",
        "com_vs_pattern_offset_um",
        "apical_fraction_pct",
        "radial_gini",
        "radial_entropy",
        "radial_cov",
        "radial_mean_r_um",
        "radial_sd_r_um",
    ]))

    rows = []
    for m in candidates:
        if m not in df.columns:
            continue
        r = test_metric(df, m)
        if r is not None:
            rows.append(r)
    res = pl.from_dicts(rows)

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    res.write_csv(OUT_DIR / "fig4B_single_zone_metrics.csv")

    # Pretty-print headline: focus on no-TRAK vs TRAK2
    print(f"=== Fig 4B: no TRAK vs TRAK2 (family size m={FAMILY}) ===\n")
    print(f"{'Metric':<55}  {'d':>7}  {'p (Šídák)':>10}  {'sig':>4}")
    print("-" * 82)
    key = "no TRAK vs TRAK2"
    # sort by p for readability, zsum first
    zrows = [r for r in rows if r["metric"].startswith("zsum_")]
    mrows = [r for r in rows if r["metric"].startswith("maxip_")]
    for group_name, group_rows in (("z-sum (SumIP)", zrows), ("MaxIP", mrows)):
        print(f"\n--- {group_name} ---")
        # order: put Mark baseline first, then the rest sorted by p
        mark_first = [r for r in group_rows if "perinuclear_5um_pct" in r["metric"] or "peripheral_5um_pct" in r["metric"]]
        rest = sorted([r for r in group_rows if r not in mark_first], key=lambda r: r.get(f"{key} p", 1.0) or 1.0)
        for r in mark_first + rest:
            d = r.get(f"{key} d")
            p = r.get(f"{key} p")
            if d is None or p is None:
                continue
            sig = "***" if p < 0.001 else "**" if p < 0.01 else "*" if p < 0.05 else "ns"
            print(f"{r['metric']:<55}  {d:+7.3f}  {p:10.4f}  {sig:>4}")

    # Same for the other two pairs
    for key in ("no TRAK vs TRAK1", "TRAK1 vs TRAK2"):
        print(f"\n\n=== Fig 4B: {key} (m={FAMILY}) ===\n")
        print(f"{'Metric':<55}  {'d':>7}  {'p (Šídák)':>10}  {'sig':>4}")
        print("-" * 82)
        for group_name, group_rows in (("z-sum (SumIP)", zrows), ("MaxIP", mrows)):
            print(f"\n--- {group_name} ---")
            mark_first = [r for r in group_rows if "perinuclear_5um_pct" in r["metric"] or "peripheral_5um_pct" in r["metric"]]
            rest = sorted([r for r in group_rows if r not in mark_first], key=lambda r: r.get(f"{key} p", 1.0) or 1.0)
            for r in mark_first + rest:
                d = r.get(f"{key} d")
                p = r.get(f"{key} p")
                if d is None or p is None:
                    continue
                sig = "***" if p < 0.001 else "**" if p < 0.01 else "*" if p < 0.05 else "ns"
                print(f"{r['metric']:<55}  {d:+7.3f}  {p:10.4f}  {sig:>4}")

    print(f"\n\nFull CSV: {OUT_DIR / 'fig4B_single_zone_metrics.csv'}")


if __name__ == "__main__":
    main()
