"""PoseBusters heatmaps with All targets summary row/column."""

from __future__ import annotations

import numpy as np
import pandas as pd

from ..config import AnalysisConfig
from ..constants import POSEBUSTERS_EXCLUDE
from .style import apply_style, plot_heatmap_with_summary_col, plot_heatmap_with_summary_row, sequential_cmap


def _pretty_check(name: str) -> str:
    return (
        name.replace("minimum_distance_to_", "min_dist_")
        .replace("volume_overlap_with_", "overlap_")
        .replace("protein-ligand_", "prot_lig_")
    )


def plot_posebusters_summary(summary: pd.DataFrame, cfg: AnalysisConfig) -> None:
    apply_style(cfg)
    cmap = sequential_cmap(cfg)
    out_dir = cfg.figures_dir / "posebusters" / "heatmaps"
    methods = [m for m in cfg.methods if m not in POSEBUSTERS_EXCLUDE]

    for rate_col, title in (
        ("pass_rate_all", "Pass all checks (100%)"),
        ("pass_rate_90pct", "Pass >=90% checks"),
        ("pass_rate_50pct", "Pass >=50% checks"),
    ):
        labels = [cfg.method_label(m) for m in methods]
        idx = [t.upper() for t in cfg.targets]
        body = pd.DataFrame(index=idx, columns=labels, dtype=float)
        for _, row in summary.iterrows():
            if row["method_id"] not in methods:
                continue
            body.loc[row["target"].upper(), cfg.method_label(row["method_id"])] = 100.0 * row[rate_col]

        if body.isna().all().all():
            continue

        sub = summary[summary["method_id"].isin(methods)].copy()
        sub["weighted"] = sub[rate_col] * sub["n_poses"]
        summary_row = pd.Series(dtype=float)
        for method_id in methods:
            label = cfg.method_label(method_id)
            part = sub[sub["method_id"] == method_id]
            if part.empty or part["n_poses"].sum() == 0:
                summary_row[label] = np.nan
            else:
                summary_row[label] = 100.0 * part["weighted"].sum() / part["n_poses"].sum()

        plot_heatmap_with_summary_row(
            summary_row,
            body.astype(float),
            title=title,
            cbar_label="Pass rate (%)",
            cmap=cmap,
            vmin=0.0,
            vmax=100.0,
            center=None,
            fmt=".1f",
            path_base=out_dir / f"heatmap_{rate_col}",
        )


def plot_posebusters_per_check(checks: pd.DataFrame, cfg: AnalysisConfig) -> None:
    apply_style(cfg)
    cmap = sequential_cmap(cfg)
    out_dir = cfg.figures_dir / "posebusters" / "heatmaps"
    methods = [m for m in cfg.methods if m not in POSEBUSTERS_EXCLUDE]

    for method_id in methods:
        sub = checks[checks.method_id == method_id].copy()
        if sub.empty:
            continue
        label = cfg.method_label(method_id)
        body = sub.pivot(index="check", columns="target", values="pass_rate")
        body.index = [_pretty_check(c) for c in body.index]
        body.columns = [c.upper() for c in body.columns]
        body = 100.0 * body

        summary_col = body.mean(axis=1, skipna=True)

        plot_heatmap_with_summary_col(
            summary_col,
            body.astype(float),
            title=f"PoseBusters per-check pass rate — {label}",
            cbar_label="Pass rate (%)",
            cmap=cmap,
            vmin=0.0,
            vmax=100.0,
            fmt="",
            path_base=out_dir / f"heatmap_checks_{method_id}",
        )
