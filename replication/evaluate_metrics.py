"""Evaluate candidate metrics from metric_pipeline output against the same
nested ANOVA framework used to replicate Mark's Fig S11 stats, ranking metrics
by effect size × significance for a target comparison (e.g. TRAK1_mDRH vs wt).

Reads replication/metrics_out/by_well/{plate}/{well}/metrics.csv produced by
metric_pipeline.py, plus config/Comparisons_table_v3.xlsx.
"""
from __future__ import annotations

import pathlib
import sys
import warnings
from dataclasses import dataclass
from itertools import combinations

import fastexcel
import numpy as np
import pandas as pd
import polars as pl
from scipy import stats

REPO = pathlib.Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / "replication"))

from replicate_stats import (  # noqa: E402
    ConditionData,
    mixedlm_condition_test,
    nested_oneway_anova,
    sidak_pairwise,
)

warnings.filterwarnings("ignore")

DEFAULT_METRICS_ROOT = REPO / "replication" / "metrics_out" / "by_well"


def _find_metrics_csv(root: pathlib.Path, plate: str, well: str) -> pathlib.Path | None:
    plate_dir = root / plate
    if not plate_dir.is_dir():
        return None
    for sub in plate_dir.iterdir():
        if sub.is_dir() and sub.name.startswith(well + "_"):
            csv = sub / "metrics.csv"
            if csv.exists():
                return csv
    return None


def collect(sheet_df: pl.DataFrame, root: pathlib.Path, metric: str) -> list[ConditionData]:
    plate_col = sheet_df.columns[0]
    out: list[ConditionData] = []
    for cond in sheet_df.columns[1:]:
        plates: dict[str, np.ndarray] = {}
        for row in sheet_df.iter_rows(named=True):
            plate = row[plate_col]
            well = row[cond]
            if not well:
                continue
            csv = _find_metrics_csv(root, plate, well)
            if csv is None:
                continue
            df = pl.read_csv(csv)
            if metric not in df.columns:
                continue
            vals = df[metric].drop_nans().drop_nulls().to_numpy().astype(float)
            if vals.size:
                plates[plate] = vals
        out.append(ConditionData(name=cond, plate_cells=plates))
    return out


def cohens_d(a: np.ndarray, b: np.ndarray) -> float:
    sa = a.std(ddof=1)
    sb = b.std(ddof=1)
    na, nb = a.size, b.size
    if na < 2 or nb < 2:
        return float("nan")
    pooled = np.sqrt(((na - 1) * sa ** 2 + (nb - 1) * sb ** 2) / (na + nb - 2))
    if pooled == 0:
        return float("nan")
    return float((a.mean() - b.mean()) / pooled)


@dataclass
class MetricResult:
    sheet: str
    metric: str
    pair: tuple[str, str]
    n_i: int
    n_j: int
    mean_i: float
    mean_j: float
    cohens_d: float
    classical_t: float
    classical_p_sidak: float
    mixedlm_z: float
    mixedlm_p_sidak: float
    f_classical: float
    f_classical_p: float
    mixedlm_joint_p: float
    var_plate: float
    var_resid: float


def evaluate(
    sheet_name: str, sheet_df: pl.DataFrame, metric: str,
    root: pathlib.Path, target_pair: tuple[str, str] | None = None,
) -> list[MetricResult]:
    conds = collect(sheet_df, root, metric)
    if sum(1 for c in conds if c.plate_cells) < 2:
        return []
    if target_pair is None:
        pairs = list(combinations(range(len(conds)), 2))
    else:
        cond_names = [c.name for c in conds]
        if target_pair[0] not in cond_names or target_pair[1] not in cond_names:
            return []
        pairs = [(cond_names.index(target_pair[0]), cond_names.index(target_pair[1]))]

    a = nested_oneway_anova(conds)
    rs_classical = sidak_pairwise(conds, a, pairs=pairs)
    mm = mixedlm_condition_test(conds, pairs=pairs)

    results: list[MetricResult] = []
    for idx, (i, j) in enumerate(pairs):
        ci, cj = conds[i], conds[j]
        if not ci.plate_cells or not cj.plate_cells:
            continue
        clr = rs_classical[idx]
        mmr = next((r for r in (mm["pairwise"] if mm and "pairwise" in mm else []) if r["i"] == ci.name and r["j"] == cj.name), None)
        results.append(MetricResult(
            sheet=sheet_name,
            metric=metric,
            pair=(ci.name, cj.name),
            n_i=ci.all_cells.size,
            n_j=cj.all_cells.size,
            mean_i=float(ci.all_cells.mean()),
            mean_j=float(cj.all_cells.mean()),
            cohens_d=cohens_d(ci.all_cells, cj.all_cells),
            classical_t=clr["t"], classical_p_sidak=clr["p_sidak"],
            mixedlm_z=mmr["z"] if mmr else float("nan"),
            mixedlm_p_sidak=mmr["p_sidak"] if mmr else float("nan"),
            f_classical=a["F"], f_classical_p=a["p"],
            mixedlm_joint_p=mm["p_joint"] if mm and "p_joint" in mm else float("nan"),
            var_plate=mm["var_plate"] if mm and "var_plate" in mm else float("nan"),
            var_resid=mm["var_resid"] if mm and "var_resid" in mm else float("nan"),
        ))
    return results


def load_comparisons(xlsx: pathlib.Path) -> dict[str, pl.DataFrame]:
    r = fastexcel.read_excel(xlsx)
    return {name: pl.read_excel(xlsx, sheet_name=name) for name in r.sheet_names}


def main():
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", default=str(DEFAULT_METRICS_ROOT))
    ap.add_argument("--sheet", required=True, help="Sheet name in Comparisons_table_v3.xlsx")
    ap.add_argument("--pair", nargs=2, default=None,
                    help="Two condition names to focus on; default = all pairs")
    ap.add_argument("--metrics", nargs="*", default=None,
                    help="Specific metric names; default = all numeric columns")
    ap.add_argument("--exclude", nargs="*", default=["path", "template_matching_score",
                                                     "lateral_pixel_pitch_um"],
                    help="Columns to skip when auto-discovering metrics")
    ap.add_argument("--out", default=None)
    args = ap.parse_args()

    root = pathlib.Path(args.root).resolve()
    sheets = load_comparisons(REPO / "config" / "Comparisons_table_v3.xlsx")
    if args.sheet not in sheets:
        print(f"Sheet not found: {args.sheet}. Available: {list(sheets)}")
        return 1
    sheet_df = sheets[args.sheet]

    # Discover metrics from any one CSV
    sample_csv = None
    for p in root.rglob("metrics.csv"):
        sample_csv = p
        break
    if sample_csv is None:
        print(f"No metrics.csv found under {root}")
        return 1
    all_cols = pl.read_csv(sample_csv).columns
    if args.metrics:
        metrics = args.metrics
    else:
        metrics = [c for c in all_cols if c not in args.exclude]

    target_pair = tuple(args.pair) if args.pair else None
    all_results: list[MetricResult] = []
    for m in metrics:
        try:
            rs = evaluate(args.sheet, sheet_df, m, root, target_pair=target_pair)
            all_results.extend(rs)
        except Exception as e:
            print(f"[skip] {m}: {e}")

    if not all_results:
        print("No results.")
        return 1

    rows = [
        {
            "metric": r.metric,
            "pair": f"{r.pair[0]} vs {r.pair[1]}",
            "n_i": r.n_i, "n_j": r.n_j,
            "mean_i": r.mean_i, "mean_j": r.mean_j,
            "diff": r.mean_i - r.mean_j,
            "cohens_d": r.cohens_d,
            "classical_t": r.classical_t,
            "classical_p_sidak": r.classical_p_sidak,
            "mixedlm_z": r.mixedlm_z,
            "mixedlm_p_sidak": r.mixedlm_p_sidak,
            "f_classical": r.f_classical,
            "f_classical_p": r.f_classical_p,
            "mixedlm_joint_p": r.mixedlm_joint_p,
            "var_plate": r.var_plate,
            "var_resid": r.var_resid,
        }
        for r in all_results
    ]
    df = pl.DataFrame(rows).sort("mixedlm_p_sidak")
    pl.Config.set_tbl_cols(20)
    pl.Config.set_tbl_rows(200)

    print(f"\n[{args.sheet}]  {len(rows)} metric×pair evaluations\n")
    print(df.select(["metric", "pair", "n_i", "n_j", "mean_i", "mean_j",
                     "cohens_d", "classical_p_sidak", "mixedlm_p_sidak",
                     "f_classical", "mixedlm_joint_p"]))

    if args.out:
        df.write_csv(args.out)
        print(f"\nWrote {args.out}")


if __name__ == "__main__":
    sys.exit(main())
