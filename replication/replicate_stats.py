"""Replicate Mark's Fig S11 / Fig 4 nested one-way ANOVA + Šídák stats.

Aggregates per-cell metrics from the existing analysis CSVs via the
Comparisons_table_v3.xlsx sheet/condition/plate/well mapping, then runs a
Prism-equivalent nested one-way ANOVA with Šídák pairwise corrections.

Reads mark_data/ only; writes output to replication/stats_out/.
"""
from __future__ import annotations

import pathlib
import re
import sys
from dataclasses import dataclass
from itertools import combinations

import warnings

import fastexcel
import numpy as np
import pandas as pd
import polars as pl
import statsmodels.formula.api as smf
from scipy import stats

warnings.filterwarnings(
    "ignore",
    message=".*Random effects covariance is singular|Maximum Likelihood optimization|"
    "MixedLM optimization failed|Gradient optimization failed|Hessian matrix|MLE may be on the boundary.*",
)


REPO = pathlib.Path(__file__).resolve().parent.parent
DATA = REPO / "mark_data"
COMPARISONS_XLSX = REPO / "config" / "Comparisons_table_v3.xlsx"

# Analysis roots in priority order (first = preferred). Plate 9 was dropped
# from the 260224 snapshot, so we fall back to older dirs for that plate.
ANALYSIS_ROOTS = [
    DATA / "analysis" / d
    for d in ("260224", "260124", "260116", "260113", "251229")
]

# Column-name variants across pipeline versions — first match wins.
METRIC_ALIASES = {
    "perinuclear_5um_percent_total": [
        "perinuclear_5um_percent_total",
        "perinuclear_percent_total",
    ],
    "peripheral_5um_simple_percent_total": [
        "peripheral_5um_simple_percent_total",
    ],
}


def _find_well_dir(plate: str, well: str) -> pathlib.Path | None:
    """Return the analysis well directory for `{plate}/{well}_*` using the
    newest analysis root that has it. Returns the top-level well dir; the
    caller picks `template_matching.csv` or `denoised/template_matching.csv`."""
    for root in ANALYSIS_ROOTS:
        plate_dir = root / plate
        if not plate_dir.is_dir():
            continue
        for sub in plate_dir.iterdir():
            if sub.is_dir() and sub.name.startswith(well + "_"):
                return sub
    return None


def _read_metric(csv: pathlib.Path, metric: str) -> np.ndarray | None:
    if not csv.exists():
        return None
    df = pl.read_csv(csv)
    for name in METRIC_ALIASES[metric]:
        if name in df.columns:
            arr = df[name].drop_nans().drop_nulls().to_numpy()
            return arr.astype(float)
    return None


@dataclass
class ConditionData:
    """Per-condition, per-plate cell-level values."""
    name: str
    # plate_name -> 1-D array of per-cell values
    plate_cells: dict[str, np.ndarray]

    @property
    def all_cells(self) -> np.ndarray:
        return (
            np.concatenate(list(self.plate_cells.values()))
            if self.plate_cells
            else np.array([])
        )

    @property
    def plate_means(self) -> np.ndarray:
        return np.array([v.mean() for v in self.plate_cells.values() if v.size > 0])

    @property
    def plate_counts(self) -> np.ndarray:
        return np.array([v.size for v in self.plate_cells.values() if v.size > 0])


def collect_sheet(
    sheet_df: pl.DataFrame, metric: str, *, use_denoised: bool
) -> list[ConditionData]:
    """Gather per-cell values per condition from Mark's analysis CSVs."""
    plate_col = sheet_df.columns[0]
    out: list[ConditionData] = []
    for cond in sheet_df.columns[1:]:
        plates: dict[str, np.ndarray] = {}
        for row in sheet_df.iter_rows(named=True):
            plate = row[plate_col]
            well = row[cond]
            if not well:
                continue
            well_dir = _find_well_dir(plate, well)
            if well_dir is None:
                print(f"  [miss] {cond} {plate}/{well}: no analysis dir", file=sys.stderr)
                continue
            csv = well_dir / ("denoised/template_matching.csv" if use_denoised else "template_matching.csv")
            vals = _read_metric(csv, metric)
            if vals is None or vals.size == 0:
                print(f"  [miss] {cond} {plate}/{well}: no {metric} in {csv}", file=sys.stderr)
                continue
            plates[plate] = vals
        out.append(ConditionData(name=cond, plate_cells=plates))
    return out


def nested_oneway_anova(conds: list[ConditionData]) -> dict:
    """Prism-equivalent nested one-way ANOVA: F-test of condition fixed effect
    using plate-within-condition as the error term. Uses cell-count-weighted
    means (matches Prism's descriptive output where the condition mean = mean
    over all cells in that condition).

    With k conditions and r_c plates in condition c:
      SS_between  = Σ_c n_c (mean_c - grand_mean)^2       # n_c = cells in c
      SS_plate_in = Σ_c Σ_p n_cp (plate_mean_cp - mean_c)^2
      MS_between  = SS_between / (k - 1)
      MS_plate_in = SS_plate_in / Σ_c (r_c - 1)
      F           = MS_between / MS_plate_in              # df = (k-1, Σ(r_c-1))
    """
    k = len(conds)
    all_cells = np.concatenate([c.all_cells for c in conds])
    grand_mean = all_cells.mean()

    ss_between = 0.0
    ss_plate_in = 0.0
    df_plate_in = 0
    cond_means: list[float] = []
    for c in conds:
        if not c.plate_cells:
            cond_means.append(float("nan"))
            continue
        cells = c.all_cells
        n_c = cells.size
        mean_c = cells.mean()
        cond_means.append(mean_c)
        ss_between += n_c * (mean_c - grand_mean) ** 2
        for plate_vals in c.plate_cells.values():
            n_cp = plate_vals.size
            plate_mean = plate_vals.mean()
            ss_plate_in += n_cp * (plate_mean - mean_c) ** 2
        df_plate_in += len(c.plate_cells) - 1

    df_between = k - 1
    ms_between = ss_between / df_between if df_between else float("nan")
    ms_plate_in = ss_plate_in / df_plate_in if df_plate_in else float("nan")
    F = ms_between / ms_plate_in if ms_plate_in else float("nan")
    p = float(stats.f.sf(F, df_between, df_plate_in)) if df_plate_in else float("nan")

    return {
        "k": k,
        "df_between": df_between,
        "df_within": df_plate_in,
        "F": F,
        "p": p,
        "ms_between": ms_between,
        "ms_within": ms_plate_in,
        "grand_mean": grand_mean,
        "cond_means": cond_means,
    }


def sidak_pairwise(
    conds: list[ConditionData],
    anova: dict,
    *,
    pairs: list[tuple[int, int]] | None = None,
) -> list[dict]:
    """Pairwise t-tests using the nested-ANOVA error term (MS_plate_within),
    with Šídák correction over the specified pair family.

    Test statistic: t = (mean_i - mean_j) / sqrt(MS_plate_in * (1/n_i + 1/n_j))
    where n_i = total cells in condition i and df = df_plate_within.
    """
    if pairs is None:
        pairs = list(combinations(range(len(conds)), 2))
    m = len(pairs)
    ms_within = anova["ms_within"]
    df = anova["df_within"]
    results = []
    for i, j in pairs:
        ci, cj = conds[i], conds[j]
        xi, xj = ci.all_cells, cj.all_cells
        ni, nj = xi.size, xj.size
        mean_i, mean_j = xi.mean(), xj.mean()
        diff = mean_i - mean_j
        se = np.sqrt(ms_within * (1 / ni + 1 / nj))
        t = diff / se if se else float("nan")
        p_raw = 2 * stats.t.sf(abs(t), df) if np.isfinite(t) else float("nan")
        p_sidak = 1 - (1 - p_raw) ** m
        results.append(
            {
                "i": ci.name, "j": cj.name,
                "mean_i": mean_i, "mean_j": mean_j,
                "n_i": ni, "n_j": nj,
                "mean_diff": diff, "se_diff": se,
                "t": t, "df": df,
                "p_raw": p_raw, "p_sidak": p_sidak,
            }
        )
    return results


def _long_dataframe(conds: list[ConditionData]) -> pd.DataFrame:
    rows = []
    for c in conds:
        for plate, vals in c.plate_cells.items():
            for v in vals:
                rows.append({"value": float(v), "condition": c.name, "plate": plate})
    return pd.DataFrame(rows)


def mixedlm_condition_test(
    conds: list[ConditionData], *, pairs: list[tuple[int, int]] | None = None
) -> dict | None:
    """Fit a mixed-effects model with random intercept per plate nested within
    condition, and test the fixed condition effect via Wald F on the
    k−1 dummy-encoded coefficients. Pairwise contrasts use Wald z tests with
    Šídák correction.

    Returns None if the design is too thin for MixedLM to fit.
    """
    df = _long_dataframe(conds)
    if df.empty or df["condition"].nunique() < 2:
        return None
    # Plate IDs must be globally unique across conditions — same plate processed
    # under different conditions is still a different cluster of cells, so
    # cluster by (condition, plate).
    df["cluster"] = df["condition"].astype(str) + "/" + df["plate"].astype(str)

    # Drop conditions with no data
    df = df[df.groupby("condition")["value"].transform("size") > 0]
    if df["condition"].nunique() < 2:
        return None

    ref = df["condition"].iloc[0]
    levels = [ref] + [c for c in df["condition"].unique() if c != ref]
    df["condition"] = pd.Categorical(df["condition"], categories=levels)

    try:
        model = smf.mixedlm("value ~ C(condition)", data=df, groups=df["cluster"])
        fit = model.fit(method="lbfgs", reml=True)
    except Exception as exc:
        return {"error": str(exc)}

    # Joint Wald F for the k-1 condition dummies
    dummy_names = [p for p in fit.params.index if p.startswith("C(condition)[T.")]
    if not dummy_names:
        return {"fit": fit, "error": "no condition dummies fit"}
    L = np.zeros((len(dummy_names), len(fit.params)))
    for i, name in enumerate(dummy_names):
        L[i, list(fit.params.index).index(name)] = 1.0
    try:
        wald = fit.wald_test(L, scalar=True)
        # wald.fvalue / pvalue come with df_denom from the profile likelihood —
        # for MixedLM statsmodels falls back to chi² (df_num=len(dummy_names)).
        p_joint = float(wald.pvalue)
        f_like = float(wald.statistic) / len(dummy_names)
    except Exception as exc:
        return {"fit": fit, "error": str(exc)}

    # Pairwise contrasts
    if pairs is None:
        pairs = list(combinations(range(len(conds)), 2))
    m = len(pairs)

    # Build contrast matrix: each row picks (mean_i - mean_j) by subtracting
    # the dummy columns. For the reference level, its dummy is 0.
    cond_names = [c.name for c in conds]
    param_names = list(fit.params.index)
    dummy_idx = {name: param_names.index(f"C(condition)[T.{name}]") for name in cond_names if name != ref}

    pairwise = []
    cov = fit.cov_params().to_numpy()
    for i, j in pairs:
        if i >= len(conds) or j >= len(conds):
            continue
        ci, cj = conds[i], conds[j]
        if not ci.plate_cells or not cj.plate_cells:
            continue
        contrast = np.zeros(len(param_names))
        if ci.name != ref:
            contrast[dummy_idx[ci.name]] += 1.0
        if cj.name != ref:
            contrast[dummy_idx[cj.name]] -= 1.0
        est = float(contrast @ fit.params.to_numpy())
        se = float(np.sqrt(contrast @ cov @ contrast))
        z = est / se if se else float("nan")
        p_raw = 2 * stats.norm.sf(abs(z)) if np.isfinite(z) else float("nan")
        p_sidak = 1 - (1 - p_raw) ** m
        pairwise.append({
            "i": ci.name, "j": cj.name,
            "mean_diff": est, "se": se, "z": z,
            "p_raw": p_raw, "p_sidak": p_sidak,
        })

    return {
        "fit": fit,
        "p_joint": p_joint,
        "wald_stat": float(wald.statistic),
        "f_like": f_like,
        "df_num": len(dummy_names),
        "var_plate": float(fit.cov_re.iloc[0, 0]) if hasattr(fit, "cov_re") and fit.cov_re.size else float("nan"),
        "var_resid": float(fit.scale),
        "pairwise": pairwise,
    }


def fmt_p(p: float) -> str:
    if p < 0.0001: return "<0.0001"
    if p < 0.001: return f"{p:.4f} ***"
    if p < 0.01: return f"{p:.4f} **"
    if p < 0.05: return f"{p:.4f} *"
    return f"{p:.4f} ns"


def run_sheet(
    sheet_name: str,
    sheet_df: pl.DataFrame,
    metric: str,
    *,
    use_denoised: bool,
    pairs_name: str = "all",
) -> None:
    label = f"{sheet_name}  |  {metric}  |  {'denoised' if use_denoised else 'raw'}"
    print(f"\n=== {label} ===")
    conds = collect_sheet(sheet_df, metric, use_denoised=use_denoised)

    header = f"{'condition':<25} {'n cells':>7} {'n plates':>9} {'mean':>8} {'plate means':<}"
    print(header)
    for c in conds:
        pm = ", ".join(f"{v:.2f}" for v in c.plate_means)
        print(f"{c.name:<25} {c.all_cells.size:>7d} {len(c.plate_cells):>9d} "
              f"{(c.all_cells.mean() if c.all_cells.size else float('nan')):>8.2f}  [{pm}]")

    if sum(1 for c in conds if c.plate_cells) < 2:
        print("  (skipping — fewer than 2 populated conditions)")
        return

    a = nested_oneway_anova(conds)
    print(f"\nClassical nested 1-way ANOVA  F({a['df_between']},{a['df_within']}) = "
          f"{a['F']:.3f}  p = {fmt_p(a['p'])}")

    if pairs_name == "adjacent":
        pairs = [(i, i + 1) for i in range(len(conds) - 1)]
    else:
        pairs = list(combinations(range(len(conds)), 2))
    rs = sidak_pairwise(conds, a, pairs=pairs)
    print(f"  Šídák pairwise ({len(pairs)} comparisons, df={a['df_within']}):")
    for r in rs:
        print(f"    {r['i']:<20} vs {r['j']:<20}  diff={r['mean_diff']:+7.3f}  "
              f"t={r['t']:+.3f}  p_raw={r['p_raw']:.4f}  p_šídák={fmt_p(r['p_sidak'])}")

    mm = mixedlm_condition_test(conds, pairs=pairs)
    if mm is None:
        print("  MixedLM: skipped (insufficient data)")
    elif "error" in mm:
        print(f"  MixedLM: FAILED — {mm['error']}")
    else:
        print(f"\nMixedLM (random intercept per plate×condition, REML)")
        print(f"  σ²_plate = {mm['var_plate']:.3f}   σ²_resid = {mm['var_resid']:.3f}")
        print(f"  joint Wald χ²({mm['df_num']}) = {mm['wald_stat']:.2f}   p = {fmt_p(mm['p_joint'])}")
        print(f"  Šídák pairwise (Wald z, {len(mm['pairwise'])} comparisons):")
        for r in mm["pairwise"]:
            print(f"    {r['i']:<20} vs {r['j']:<20}  diff={r['mean_diff']:+7.3f}  "
                  f"z={r['z']:+.3f}  p_raw={r['p_raw']:.4f}  p_šídák={fmt_p(r['p_sidak'])}")


def load_comparisons() -> dict[str, pl.DataFrame]:
    r = fastexcel.read_excel(COMPARISONS_XLSX)
    return {
        name: pl.read_excel(COMPARISONS_XLSX, sheet_name=name)
        for name in r.sheet_names
    }


def main() -> int:
    sheets = load_comparisons()
    # (sheet, metric, use_denoised, pairs_mode)
    # Mark's Prism templates pair the perinuclear metric with raw z-sum
    # projections, and the peripheral metric with denoised.
    plan = [
        ("TRAK1 helix muts", "perinuclear_5um_percent_total", False, "adjacent"),
        ("TRAK1 helix muts", "peripheral_5um_simple_percent_total", True, "adjacent"),
        ("TRAK2 helix muts", "perinuclear_5um_percent_total", False, "adjacent"),
        ("TRAK2 helix muts", "peripheral_5um_simple_percent_total", True, "adjacent"),
        ("TRAK isoform (mito)", "perinuclear_5um_percent_total", False, "all"),
        ("TRAK isoform (mito)", "peripheral_5um_simple_percent_total", True, "all"),
        ("TRAK isoform (peroxisome)", "perinuclear_5um_percent_total", False, "all"),
        ("TRAK isoform (peroxisome)", "peripheral_5um_simple_percent_total", True, "all"),
        ("MAPK9 siRNA", "perinuclear_5um_percent_total", False, "all"),
        ("MAPK9 siRNA", "peripheral_5um_simple_percent_total", True, "all"),
        ("TRAK2 S84 Ars", "perinuclear_5um_percent_total", False, "all"),
        ("TRAK2 S84 Ars", "peripheral_5um_simple_percent_total", True, "all"),
    ]
    for sheet, metric, denoised, pairs in plan:
        if sheet not in sheets:
            print(f"\n=== {sheet} — MISSING from comparisons table ===")
            continue
        run_sheet(sheet, sheets[sheet], metric, use_denoised=denoised, pairs_name=pairs)
    return 0


if __name__ == "__main__":
    sys.exit(main())
