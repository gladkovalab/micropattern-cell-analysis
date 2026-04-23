"""Derive richer distribution metrics from Mark's existing per-cell CSVs.

Mark's CSVs already contain the radial signal profile at 5 discrete bins
(perinuclear 1..5 µm, peripheral 1..5 µm). That's enough to define shape-of-
distribution metrics without re-processing the raw ND2s. Each cell's row in
the CSV is converted into ~25 derived metrics, which are then evaluated with
nested ANOVA + Šídák per sheet and pair.

Outputs:
    replication/derived_metrics_out/per_cell.csv — per-cell metrics with sheet/condition/plate tags
    replication/derived_metrics_out/per_metric_summary.csv — one row per (sheet, pair, metric)
"""
from __future__ import annotations

import pathlib
import sys
import warnings
from itertools import combinations

import fastexcel
import numpy as np
import polars as pl

REPO = pathlib.Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / "replication"))
from replicate_stats import (  # noqa: E402
    ConditionData,
    mixedlm_condition_test,
    nested_oneway_anova,
    sidak_pairwise,
)
from evaluate_metrics import cohens_d  # noqa: E402

warnings.filterwarnings("ignore")

DATA = REPO / "mark_data"
COMPARISONS_XLSX = REPO / "config" / "Comparisons_table_v3.xlsx"
OUT = REPO / "replication" / "derived_metrics_out"

ANALYSIS_ROOTS = [
    DATA / "analysis" / d
    for d in ("260224", "260124", "260116", "260113", "251229")
]

DISTANCES = (1, 2, 3, 4, 5)


def _find_well_dir(plate: str, well: str) -> pathlib.Path | None:
    for root in ANALYSIS_ROOTS:
        plate_dir = root / plate
        if not plate_dir.is_dir():
            continue
        for sub in plate_dir.iterdir():
            if sub.is_dir() and sub.name.startswith(well + "_"):
                return sub
    return None


def _read_csv_values(csv: pathlib.Path) -> dict[str, np.ndarray] | None:
    """Read per-cell values for all distance bins (raw or denoised CSV).
    Returns a dict keyed by canonicalized column name.
    """
    if not csv.exists():
        return None
    df = pl.read_csv(csv)
    out: dict[str, np.ndarray] = {}
    # canonicalize: both old (`perinuclear_percent_total`) and new
    # (`perinuclear_5um_percent_total`) schemas. We only harvest 1..5 µm variants.
    for d in DISTANCES:
        for metric in ("perinuclear", "peripheral"):
            # new schema
            suffix = "percent_total" if metric == "perinuclear" else "simple_percent_total"
            col_new = f"{metric}_{d}um_{suffix}"
            col_old = None
            if d == 5:
                col_old = f"{metric}_percent_total" if metric == "perinuclear" else f"{metric}_5um_simple_percent_total"
            if col_new in df.columns:
                out[f"{metric}_{d}um"] = df[col_new].to_numpy().astype(float)
            elif col_old and col_old in df.columns:
                out[f"{metric}_{d}um"] = df[col_old].to_numpy().astype(float)
    return out or None


def _cell_level_derived(raw: dict[str, np.ndarray], den: dict[str, np.ndarray] | None) -> dict[str, np.ndarray]:
    """Compute per-cell derived metrics. Each entry maps metric_name -> 1D array."""
    out: dict[str, np.ndarray] = {}
    # Pass-through of Mark's primary metrics for side-by-side comparison.
    for d in DISTANCES:
        out[f"raw_perinuclear_{d}um"] = raw.get(f"perinuclear_{d}um", np.array([]))
        out[f"raw_peripheral_{d}um"] = raw.get(f"peripheral_{d}um", np.array([]))
        if den is not None:
            out[f"den_perinuclear_{d}um"] = den.get(f"perinuclear_{d}um", np.array([]))
            out[f"den_peripheral_{d}um"] = den.get(f"peripheral_{d}um", np.array([]))

    def _get(src: dict[str, np.ndarray], metric: str, d: int) -> np.ndarray:
        return src.get(f"{metric}_{d}um", np.array([]))

    def _radial_bins(src: dict[str, np.ndarray], metric: str, prefix: str):
        """Non-cumulative annular bins from cumulative-from-boundary percentages.
        Only emits bin d-1_d when both endpoints are available; older CSVs only
        carry the 5 µm value, so intermediate bins will be absent for those wells."""
        by_d = {d: _get(src, metric, d) for d in DISTANCES}
        ref_n = next((a.size for a in by_d.values() if a.size), 0)
        if ref_n == 0:
            return
        # bin 0-1 = metric_1um
        if by_d[1].size == ref_n:
            out[f"{prefix}_{metric}_bin_0_1"] = by_d[1]
        for lo, hi in [(1, 2), (2, 3), (3, 4), (4, 5)]:
            a_lo = by_d[lo]
            a_hi = by_d[hi]
            if a_lo.size == ref_n and a_hi.size == ref_n:
                out[f"{prefix}_{metric}_bin_{lo}_{hi}"] = a_hi - a_lo

    _radial_bins(raw, "perinuclear", "raw")
    _radial_bins(raw, "peripheral", "raw")
    if den is not None:
        _radial_bins(den, "perinuclear", "den")
        _radial_bins(den, "peripheral", "den")

    def _same_shape(*arrs):
        if not arrs or not arrs[0].size:
            return False
        n = arrs[0].size
        return all(a.size == n for a in arrs)

    for src_name, src in ((("raw", raw),) if den is None else (("raw", raw), ("den", den))):
        peri_5 = _get(src, "perinuclear", 5)
        peri_1 = _get(src, "perinuclear", 1)
        per_5 = _get(src, "peripheral", 5)
        per_1 = _get(src, "peripheral", 1)
        if _same_shape(peri_5, per_5):
            out[f"{src_name}_midcyto_pct"] = np.clip(100 - peri_5 - per_5, 0, None)
            out[f"{src_name}_peri_over_nuc"] = per_5 / (peri_5 + 1e-9)
            out[f"{src_name}_peri_minus_nuc"] = per_5 - peri_5
            # log-ratio: symmetric, compresses tails. log2 so "+1" == "2× more peripheral".
            out[f"{src_name}_log2_peri_over_nuc"] = np.log2((per_5 + 1e-3) / (peri_5 + 1e-3))
            # "peripheral share of zoned signal": bounded in [0, 1], avoids the unbounded
            # ratio's sensitivity when peri_5 is tiny. Reviewers usually prefer bounded %.
            out[f"{src_name}_peripheral_share_pct"] = 100 * per_5 / (per_5 + peri_5 + 1e-9)
            # angular representation: atan2 in degrees, 0° = all perinuclear, 90° = all peripheral
            out[f"{src_name}_polarization_angle_deg"] = np.degrees(np.arctan2(per_5, peri_5 + 1e-9))
        if _same_shape(peri_1, peri_5):
            out[f"{src_name}_perinuclear_concentration"] = peri_1 / (peri_5 + 1e-9)
        if _same_shape(per_1, per_5):
            out[f"{src_name}_peripheral_concentration"] = per_1 / (per_5 + 1e-9)
        b01 = out.get(f"{src_name}_perinuclear_bin_0_1")
        b45 = out.get(f"{src_name}_perinuclear_bin_4_5")
        if b01 is not None and b45 is not None and _same_shape(b01, b45):
            out[f"{src_name}_perinuclear_inner_outer_ratio"] = b01 / (b45 + 1e-9)

    return out


def collect_all(sheets: dict[str, pl.DataFrame]) -> pl.DataFrame:
    """Build a long per-cell DataFrame with all derived metrics + metadata."""
    rows = []
    for sheet_name, df in sheets.items():
        plate_col = df.columns[0]
        for cond in df.columns[1:]:
            for row in df.iter_rows(named=True):
                plate = row[plate_col]
                well = row[cond]
                if not well:
                    continue
                well_dir = _find_well_dir(plate, well)
                if well_dir is None:
                    continue
                raw_csv = well_dir / "template_matching.csv"
                den_csv = well_dir / "denoised" / "template_matching.csv"
                raw = _read_csv_values(raw_csv)
                den = _read_csv_values(den_csv)
                if raw is None:
                    continue
                derived = _cell_level_derived(raw, den)
                n_cells = next(iter(raw.values())).size
                for i in range(n_cells):
                    row_out = {
                        "sheet": sheet_name,
                        "condition": cond,
                        "plate": plate,
                        "well": well,
                        "cell_idx": i,
                    }
                    for k, arr in derived.items():
                        row_out[k] = float(arr[i]) if i < arr.size else None
                    rows.append(row_out)
    return pl.from_dicts(rows)


def _collect_all_conditions(long: pl.DataFrame, sheet: str, metric: str) -> list[ConditionData]:
    sub = long.filter(
        (pl.col("sheet") == sheet) &
        pl.col(metric).is_not_null() &
        pl.col(metric).is_not_nan()
    )
    conds: list[ConditionData] = []
    for cond_name in sub["condition"].unique().to_list():
        g = sub.filter(pl.col("condition") == cond_name)
        plate_cells: dict[str, np.ndarray] = {}
        for plate, grp in g.group_by("plate"):
            plate_key = plate[0] if isinstance(plate, tuple) else plate
            plate_cells[plate_key] = grp[metric].to_numpy().astype(float)
        conds.append(ConditionData(name=cond_name, plate_cells=plate_cells))
    return conds


def evaluate_one(long: pl.DataFrame, sheet: str, metric: str,
                 pair: tuple[str, str]) -> dict | None:
    """Evaluate a single pair against the FULL-sheet nested ANOVA, matching
    Prism's behaviour: the error term (MS_plate_within_condition) is pooled
    across every condition in the sheet, not just the two being contrasted.
    Šídák correction is applied upstream via the family-size passed by the
    caller; here we report the raw p-value — the plot layer then applies
    `1-(1-p)^m` with the panel's family size."""
    all_conds = _collect_all_conditions(long, sheet, metric)
    if len(all_conds) < 2:
        return None
    name_to_idx = {c.name: i for i, c in enumerate(all_conds)}
    if pair[0] not in name_to_idx or pair[1] not in name_to_idx:
        return None
    i, j = name_to_idx[pair[0]], name_to_idx[pair[1]]
    ci, cj = all_conds[i], all_conds[j]
    if not ci.plate_cells or not cj.plate_cells:
        return None

    a = nested_oneway_anova(all_conds)
    # Pairwise t using the sheet-wide pooled error term.
    # Note: passing pairs=[(i,j)] makes sidak_pairwise return m=1 (= raw p);
    # we store raw p, and the plotter applies the panel-specific Šídák.
    clr = sidak_pairwise(all_conds, a, pairs=[(i, j)])[0]
    mm = mixedlm_condition_test(all_conds, pairs=[(i, j)])
    mmr = None
    if mm and "pairwise" in mm:
        mmr = mm["pairwise"][0] if mm["pairwise"] else None
    return {
        "sheet": sheet,
        "metric": metric,
        "pair": f"{pair[0]} vs {pair[1]}",
        "n_i": ci.all_cells.size,
        "n_j": cj.all_cells.size,
        "mean_i": float(ci.all_cells.mean()),
        "mean_j": float(cj.all_cells.mean()),
        "diff": float(ci.all_cells.mean() - cj.all_cells.mean()),
        "cohens_d": cohens_d(ci.all_cells, cj.all_cells),
        "f_classical": a["F"],
        "p_classical_anova": a["p"],
        "t_classical_pair": clr["t"],
        "p_classical_sidak": clr["p_sidak"],  # raw p; plotter Šídák-corrects
        "p_mixedlm_joint": mm["p_joint"] if mm and "p_joint" in mm else float("nan"),
        "p_mixedlm_sidak": mmr["p_sidak"] if mmr else float("nan"),
        "var_plate": mm["var_plate"] if mm and "var_plate" in mm else float("nan"),
        "var_resid": mm["var_resid"] if mm and "var_resid" in mm else float("nan"),
    }


def main():
    sheets = {
        name: pl.read_excel(COMPARISONS_XLSX, sheet_name=name)
        for name in fastexcel.read_excel(COMPARISONS_XLSX).sheet_names
    }
    print("Collecting per-cell metrics…", flush=True)
    long = collect_all(sheets)
    OUT.mkdir(parents=True, exist_ok=True)
    long.write_csv(OUT / "per_cell.csv")
    print(f"  wrote {OUT / 'per_cell.csv'}  ({long.height} cells, {len(long.columns)-5} metrics)")

    # Evaluate every metric × every adjacent pair × every sheet.
    # For Fig S11 / Fig 4 baseline pairs we focus on the comparisons
    # Mark chose in his Prism files (adjacent columns).
    metric_cols = [c for c in long.columns
                   if c not in ("sheet", "condition", "plate", "well", "cell_idx")]
    # Pair choices mirror Mark's Prism selections (see Prism files under
    # mark_data/analysis/260224/prism_plots/). For MAPK9 the perinuclear and
    # peripheral panels used DIFFERENT pair sets — we evaluate the union so
    # both can be annotated per metric.
    sheet_pairs = {
        "TRAK1 helix muts": [("T1 wt", "T1 mDRH"), ("T1 mDRH", "T1 mDRH / dSp")],
        "TRAK2 helix muts": [("TRAK2", "TRAK2 mDRH"), ("TRAK2 mDRH", "TRAK2 mDRH mSpindly")],
        "TRAK isoform (mito)": [("no TRAK", "TRAK1"), ("no TRAK", "TRAK2"), ("TRAK1", "TRAK2")],
        "TRAK isoform (peroxisome)": [("no TRAK", "TRAK1"), ("no TRAK", "TRAK2"), ("TRAK1", "TRAK2")],
        "MAPK9 siRNA": [
            # Mark's Prism column order (different from the comparisons sheet):
            #   A = ctrl ctrl   B = MAPK9 ctrl   C = ctrl Ars   D = MAPK9 Ars
            # PUBLISHED Fig 4E peripheral uses A-B, A-C, A-D
            ("ctrl ctrl", "MAPK9 ctrl"), ("ctrl ctrl", "ctrl Ars"), ("ctrl ctrl", "MAPK9 Ars"),
            # PUBLISHED Fig S11 F perinuclear uses A-C, B-D, C-D
            #   = (ctrl ctrl, ctrl Ars), (MAPK9 ctrl, MAPK9 Ars), (ctrl Ars, MAPK9 Ars)
            ("MAPK9 ctrl", "MAPK9 Ars"), ("ctrl Ars", "MAPK9 Ars"),
        ],
        "TRAK2 S84 Ars": [("wt ctrl", "wt Ars"), ("wt Ars", "S84A Ars"),
                          ("S84A ctrl", "S84A Ars")],
    }

    results: list[dict] = []
    total = sum(len(pairs) for pairs in sheet_pairs.values()) * len(metric_cols)
    print(f"Evaluating {total} (sheet × pair × metric) combos…", flush=True)
    done = 0
    for sheet, pairs in sheet_pairs.items():
        for pair in pairs:
            for m in metric_cols:
                try:
                    r = evaluate_one(long, sheet, m, pair)
                    if r is not None:
                        results.append(r)
                except Exception as e:
                    pass
                done += 1
    df = pl.from_dicts(results)
    df.write_csv(OUT / "per_metric_summary.csv")
    print(f"  wrote {OUT / 'per_metric_summary.csv'}  ({df.height} rows)")


if __name__ == "__main__":
    sys.exit(main())
