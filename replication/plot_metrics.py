"""Plot wedge-r profile / CDF / scalars from the patched template_matching_bulk
pipeline output.

Reads per-well CSVs from `template_matching/{plate}/{well}/template_matching.csv`
(written by template_matching_bulk.main), joins them against
`config/Comparisons_table_v3.xlsx` for sheet/condition metadata, and
produces a single per-sheet figure with three panels:

  1. wedge-r intensity profile (mean ± SEM per condition)
  2. wedge-r CDF (per-cell CDFs as faint traces; mean per condition bold)
  3. strip plots of perinuclear_5um_percent_total, peripheral_5um_percent_total,
     wedge_r_ks_vs_uniform, wedge_r_ks_vs_60merNoTRAK, with nested-ANOVA +
     Šídák pairwise brackets (Welch fallback for single-plate sheets).

CLI:
    pixi run python replication/plot_metrics.py \\
        --template-matching template_matching \\
        --sheet "TRAK isoform (mito)" \\
        --out replication/figures_wedge_r_ks/trak_isoform_mito.png
"""
from __future__ import annotations

import argparse
import colorsys
import math
import pathlib
import re
import sys
from dataclasses import dataclass

import fastexcel
import matplotlib
matplotlib.use("Agg")
import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import numpy as np
import polars as pl
from scipy import stats

REPO = pathlib.Path(__file__).resolve().parent.parent

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
        "pairs": [("TRAK2", "TRAK2 mDRH"), ("TRAK2 mDRH", "TRAK2 mDRH mSpindly")],
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

CONDITION_COLORS = ["#4c78a8", "#59a14f", "#e15759", "#f28e2b", "#b07aa1"]


def _shade(base, factor: float) -> tuple:
    """Lighter (factor>1) or darker (factor<1) HLS-shifted variant of `base`."""
    rgb = mcolors.to_rgb(base)
    h, l, s = colorsys.rgb_to_hls(*rgb)
    return colorsys.hls_to_rgb(h, max(0.05, min(0.95, l * factor)), s)


def _plate_shades(base, n: int) -> list[tuple]:
    """`n` lightness variants of `base`, evenly spaced across [0.7×, 1.3×];
    the middle index reproduces the base color when n is odd."""
    if n <= 1:
        return [mcolors.to_rgb(base)]
    factors = np.linspace(0.7, 1.3, n)
    return [_shade(base, f) for f in factors]


def _wedge_r_columns(df: pl.DataFrame) -> list[str]:
    pat = re.compile(r"^wedge_r_(\d{2})_(\d{2})um_pct$")
    matches = []
    for c in df.columns:
        m = pat.match(c)
        if m:
            matches.append((int(m.group(1)), c))
    matches.sort()
    return [c for _, c in matches]


def _per_cell_cdf(profile: np.ndarray) -> np.ndarray:
    row_sums = np.nansum(profile, axis=1, keepdims=True)
    out = np.full_like(profile, np.nan, dtype=float)
    ok = (row_sums.flatten() > 0)
    out[ok] = np.nancumsum(profile[ok] / row_sums[ok], axis=1)
    return out


def _sidak(p_raw: float, m: int) -> float:
    if not np.isfinite(p_raw):
        return p_raw
    return 1 - (1 - min(max(p_raw, 0.0), 1.0)) ** m


# ---- Stats: nested one-way ANOVA + Šídák pairwise (matches v3 evaluator) ----
# Ported from replication/replicate_stats.py (on wpg/alt-metrics) so this
# branch stays self-contained. Cell-level Welch t-tests treat every cell as
# independent and inflate significance for clustered designs; nested ANOVA
# uses plate-within-condition as the error term to give Prism-equivalent
# p-values.

@dataclass
class _ConditionData:
    name: str
    plate_cells: dict[str, np.ndarray]

    @property
    def all_cells(self) -> np.ndarray:
        return (np.concatenate(list(self.plate_cells.values()))
                if self.plate_cells else np.array([]))


def _collect(df: pl.DataFrame, metric: str) -> list[_ConditionData]:
    sub = df.filter(pl.col(metric).is_not_null() & pl.col(metric).is_not_nan())
    out = []
    for cn in sorted(sub["condition"].unique().to_list()):
        g = sub.filter(pl.col("condition") == cn)
        plates = {}
        for plate, grp in g.group_by("plate"):
            key = plate[0] if isinstance(plate, tuple) else plate
            plates[key] = grp[metric].to_numpy().astype(float)
        out.append(_ConditionData(name=cn, plate_cells=plates))
    return out


def _nested_oneway_anova(conds: list[_ConditionData]) -> dict:
    """SS_between (n_c-weighted) over MS_plate_within. Returns {ms_within,
    df_within, p, F}. NaN if too thin."""
    k = len(conds)
    all_cells = np.concatenate([c.all_cells for c in conds])
    if all_cells.size == 0:
        return {"k": k, "ms_within": float("nan"), "df_within": 0,
                "F": float("nan"), "p": float("nan")}
    grand = all_cells.mean()
    ss_b, ss_w, df_w = 0.0, 0.0, 0
    for c in conds:
        if not c.plate_cells:
            continue
        cells = c.all_cells
        n_c = cells.size
        m_c = cells.mean()
        ss_b += n_c * (m_c - grand) ** 2
        for vals in c.plate_cells.values():
            ss_w += vals.size * (vals.mean() - m_c) ** 2
        df_w += len(c.plate_cells) - 1
    df_b = k - 1
    ms_b = ss_b / df_b if df_b else float("nan")
    ms_w = ss_w / df_w if df_w else float("nan")
    F = ms_b / ms_w if ms_w else float("nan")
    p = float(stats.f.sf(F, df_b, df_w)) if df_w else float("nan")
    return {"k": k, "ms_within": ms_w, "df_within": df_w, "F": F, "p": p}


def _nested_pair_p(conds: list[_ConditionData], anova: dict,
                   i: int, j: int) -> float:
    """Pairwise t using the nested-ANOVA error term."""
    ms_w = anova["ms_within"]
    df = anova["df_within"]
    if not np.isfinite(ms_w) or df <= 0:
        return float("nan")
    xi, xj = conds[i].all_cells, conds[j].all_cells
    if xi.size < 2 or xj.size < 2:
        return float("nan")
    se = np.sqrt(ms_w * (1 / xi.size + 1 / xj.size))
    if se == 0:
        return float("nan")
    t = (xi.mean() - xj.mean()) / se
    return float(2 * stats.t.sf(abs(t), df))


def _test_pair(df: pl.DataFrame, metric: str, pair: tuple[str, str],
               family: int) -> float:
    """Return Šídák-corrected p for a single (a, b) pair on `metric`. Tries
    nested ANOVA first; falls back to Welch on pooled cells when the design
    has only one plate per condition."""
    conds = _collect(df, metric)
    n2i = {c.name: i for i, c in enumerate(conds)}
    if pair[0] not in n2i or pair[1] not in n2i:
        return float("nan")
    i, j = n2i[pair[0]], n2i[pair[1]]
    p_raw = float("nan")
    try:
        a = _nested_oneway_anova(conds)
        if np.isfinite(a["p"]):
            p_raw = _nested_pair_p(conds, a, i, j)
    except Exception:
        pass
    if not np.isfinite(p_raw):
        a_vals = conds[i].all_cells
        b_vals = conds[j].all_cells
        if a_vals.size > 1 and b_vals.size > 1:
            p_raw = float(stats.ttest_ind(a_vals, b_vals, equal_var=False).pvalue)
    if not np.isfinite(p_raw):
        return float("nan")
    return _sidak(p_raw, family)


def _format_p(p: float) -> str:
    if not np.isfinite(p):
        return "ns"
    if p < 0.001:
        return f"*** p={p:.1e}"
    if p < 0.01:
        return f"** p={p:.3f}"
    if p < 0.05:
        return f"* p={p:.3f}"
    return f"ns p={p:.3f}"


def load_template_matching(tm_root: pathlib.Path) -> pl.DataFrame:
    """Concatenate all per-well CSVs under {tm_root}/**/template_matching.csv."""
    csvs = sorted(tm_root.rglob("template_matching.csv"))
    if not csvs:
        raise FileNotFoundError(f"No template_matching.csv files under {tm_root}")
    parts = [pl.read_csv(c) for c in csvs]
    common = set(parts[0].columns)
    for p in parts[1:]:
        common &= set(p.columns)
    common = list(common)
    return pl.concat([p.select(common) for p in parts], how="vertical")


def load_comparisons_table(xlsx: pathlib.Path) -> pl.DataFrame:
    """Long-form (plate, well, sheet, condition) view of Mark's
    `config/Comparisons_table_v3.xlsx`. See `run_pipeline_paths.load_comparisons_table`
    for shape details — duplicated here to keep both scripts independently runnable."""
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


def join_with_metadata(df: pl.DataFrame, comparisons_xlsx: pathlib.Path) -> pl.DataFrame:
    """Join the pipeline output to (sheet, condition) metadata via
    (plate, well_short). Pipeline-side keys come from the path; the
    metadata side comes from Mark's Comparisons table."""
    plate_re = r"patterned_data/([^/]+)/"
    well_re = r"patterned_data/[^/]+/([A-Z]\d+)_"

    df = df.with_columns([
        pl.col("path").str.extract(plate_re, 1).alias("_plate"),
        pl.col("path").str.extract(well_re, 1).alias("_well"),
    ])
    meta = load_comparisons_table(comparisons_xlsx).rename(
        {"plate": "_plate", "well": "_well"})
    # Each (plate, well) maps to one (sheet, condition) for the sheet of
    # interest; collapse duplicates so the join doesn't fan out.
    meta = meta.unique(subset=["_plate", "_well", "sheet"])

    out = df.join(meta, on=["_plate", "_well"], how="left")
    out = out.filter(pl.col("condition").is_not_null())
    out = out.rename({"_plate": "plate", "_well": "well"})
    return out


def make_figure(df: pl.DataFrame, sheet: str, out_path: pathlib.Path):
    cfg = SHEET_CONFIG[sheet]
    conditions = cfg["conditions"]
    pairs = cfg["pairs"]
    family_m = cfg["family_m"]

    sheet_df = df.filter(pl.col("sheet") == sheet)
    if sheet_df.height == 0:
        print(f"[plot_metrics] no rows for sheet {sheet!r}, skipping")
        return
    print(f"[plot_metrics] {sheet}: {sheet_df.height} cells")

    wedge_cols = _wedge_r_columns(sheet_df)
    n_bins = len(wedge_cols)
    centers = np.array([i + 0.5 for i in range(n_bins)])
    color_map = {c: CONDITION_COLORS[i % len(CONDITION_COLORS)]
                 for i, c in enumerate(conditions)}

    rng = np.random.default_rng(0)
    fig = plt.figure(figsize=(15, 11))
    gs = fig.add_gridspec(2, 4, height_ratios=[1, 1])
    ax_prof = fig.add_subplot(gs[0, 0:2])
    ax_cdf = fig.add_subplot(gs[0, 2:4])

    # --- Panels 1+2: wedge-r profile + CDF
    for cond in conditions:
        sub = sheet_df.filter(pl.col("condition") == cond)
        if sub.height == 0:
            continue
        prof = sub.select(wedge_cols).to_numpy()
        mean = np.nanmean(prof, axis=0)
        sem = np.nanstd(prof, axis=0, ddof=1) / np.sqrt(np.maximum(prof.shape[0], 1))
        col = color_map[cond]
        ax_prof.plot(centers, mean, color=col, lw=1.6,
                     label=f"{cond} (n={sub.height})")
        ax_prof.fill_between(centers, mean - sem, mean + sem,
                             color=col, alpha=0.18, linewidth=0)

        cdf_per_cell = _per_cell_cdf(prof)
        cdf_mean = np.nanmean(cdf_per_cell, axis=0)
        cdf_sem = (np.nanstd(cdf_per_cell, axis=0, ddof=1) /
                   np.sqrt(np.maximum(cdf_per_cell.shape[0], 1)))
        ax_cdf.plot(centers, cdf_mean, color=col, lw=2.0,
                    label=f"{cond} (n={sub.height})")
        ax_cdf.fill_between(centers, cdf_mean - cdf_sem, cdf_mean + cdf_sem,
                            color=col, alpha=0.18, linewidth=0)
    cdf_uni = (np.arange(1, n_bins + 1) ** 2) / (n_bins ** 2)
    ax_cdf.plot(centers, cdf_uni, color="black", linestyle="--",
                linewidth=1, label="area-uniform")
    ax_prof.set_xlabel("wedge-r (µm from apex)")
    ax_prof.set_ylabel("mean intensity per bin (% of wedge total)")
    ax_prof.set_title("Wedge-r intensity profile")
    ax_prof.legend(fontsize=8, loc="upper right")
    ax_prof.grid(alpha=0.3)
    ax_cdf.set_xlabel("wedge-r (µm from apex)")
    ax_cdf.set_ylabel("intensity-weighted CDF")
    ax_cdf.set_title("Wedge-r CDF (mean ± SEM per condition)")
    ax_cdf.legend(fontsize=8, loc="lower right")
    ax_cdf.grid(alpha=0.3)
    ax_cdf.set_ylim(0, 1.02)

    # Bottom row: 4 columns, one per scalar metric
    metric_specs = [
        ("peripheral_5um_percent_total", "peripheral 5µm\n(% of cell signal)"),
        ("perinuclear_5um_percent_total", "perinuclear 5µm\n(% of cell signal)"),
        ("wedge_r_ks_vs_uniform", "wedge-r KS\nvs area-uniform"),
        ("wedge_r_ks_vs_60merNoTRAK", "wedge-r KS\nvs 60mer no-TRAK"),
    ]
    plates = sorted(sheet_df["plate"].unique().to_list())
    plate_idx = {p: i for i, p in enumerate(plates)}
    plate_palette = {cond: _plate_shades(color_map[cond], len(plates))
                     for cond in conditions}

    for col_i, (metric, ylab) in enumerate(metric_specs):
        ax = fig.add_subplot(gs[1, col_i])
        if metric not in sheet_df.columns:
            ax.text(0.5, 0.5, f"missing\n{metric}", ha="center", va="center",
                    transform=ax.transAxes)
            ax.set_xticks([])
            continue
        all_vals = []
        for j, cond in enumerate(conditions):
            cond_sub = sheet_df.filter(pl.col("condition") == cond)
            plate_means = []  # (x, mean, shade) tuples
            for plate in plates:
                psub = cond_sub.filter(pl.col("plate") == plate)
                vals = psub[metric].to_numpy()
                vals = vals[np.isfinite(vals)]
                if vals.size == 0:
                    continue
                shade = plate_palette[cond][plate_idx[plate]]
                xs = j + rng.uniform(-0.18, 0.18, size=vals.size)
                ax.scatter(xs, vals, s=28, color=shade,
                           edgecolor="black", linewidth=0.3, alpha=0.85)
                plate_means.append((j + rng.uniform(-0.08, 0.08),
                                    float(np.mean(vals)), shade))
                all_vals.extend(vals.tolist())
            # Per-plate (biological-replicate) means as larger diamonds.
            for mx, mv, sh in plate_means:
                ax.scatter(mx, mv, s=110, color=sh, marker="D",
                           edgecolor="black", linewidth=0.9, zorder=6)
            # Overall condition mean line, drawn last so it sits on top.
            cond_vals = cond_sub[metric].to_numpy()
            cond_vals = cond_vals[np.isfinite(cond_vals)]
            if cond_vals.size > 0:
                ax.hlines(float(np.mean(cond_vals)), j - 0.32, j + 0.32,
                          color="black", linewidth=2.2, zorder=7)
        ax.set_xticks(range(len(conditions)))
        ax.set_xticklabels(conditions, rotation=20, ha="right", fontsize=8)
        ax.set_title(ylab, fontsize=10)
        ax.grid(axis="y", alpha=0.3)

        if not all_vals:
            continue
        data_max = max(all_vals)
        data_min = min(all_vals)
        span = data_max - data_min
        bracket_base = data_max + span * 0.05
        bracket_step = span * 0.10
        ax.set_ylim(data_min - span * 0.05, bracket_base + bracket_step * (len(pairs) + 0.5))

        for k, (a, b) in enumerate(pairs):
            p = _test_pair(sheet_df, metric, (a, b), family_m)
            if not np.isfinite(p):
                continue
            xa = conditions.index(a)
            xb = conditions.index(b)
            y = bracket_base + bracket_step * (k + 0.3)
            tick = bracket_step * 0.18
            ax.plot([xa, xa, xb, xb], [y - tick, y, y, y - tick],
                    color="black", linewidth=0.9)
            ax.text((xa + xb) / 2, y + bracket_step * 0.05, _format_p(p),
                    ha="center", va="bottom", fontsize=7.5)

    fig.suptitle(f"{sheet}  ·  wedge-r KS  ·  nested ANOVA + Šídák m={family_m}",
                 fontsize=12)
    plt.tight_layout(rect=[0, 0.0, 1, 0.96])
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=140)
    plt.close(fig)
    print(f"[plot_metrics] wrote {out_path}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--template-matching", default="template_matching",
                    help="Root directory of per-well template_matching.csv files "
                         "produced by the patched pipeline.")
    ap.add_argument("--comparisons-xlsx",
                    default=str(REPO / "config/Comparisons_table_v3.xlsx"),
                    help="Authoritative sheet/condition/plate/well map.")
    ap.add_argument("--sheet", required=True,
                    help='Sheet to plot, e.g. "TRAK isoform (mito)".')
    ap.add_argument("--out", required=True,
                    help="Output figure path (PNG).")
    args = ap.parse_args()

    tm = pathlib.Path(args.template_matching).resolve()
    df = load_template_matching(tm)
    df = join_with_metadata(df, pathlib.Path(args.comparisons_xlsx).resolve())
    print(f"[plot_metrics] joined: {df.height} cells")

    if args.sheet not in SHEET_CONFIG:
        sys.exit(f"unknown sheet {args.sheet!r}; "
                 f"choices: {list(SHEET_CONFIG.keys())}")

    make_figure(df, args.sheet, pathlib.Path(args.out).resolve())


if __name__ == "__main__":
    main()
