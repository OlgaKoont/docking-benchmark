"""Correlation heatmaps (article palette, All targets summary row)."""

from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.stats import kendalltau, pearsonr, spearmanr

from ..config import AnalysisConfig
from ..data import load_merged, prepare_xy
from .style import apply_style, diverging_cmap, plot_heatmap_with_summary_row


def _pivot(corr: pd.DataFrame, value_col: str, cfg: AnalysisConfig) -> pd.DataFrame:
    labels = [cfg.method_label(m) for m in cfg.methods]
    idx = [t.upper() for t in cfg.targets]
    mat = pd.DataFrame(index=idx, columns=labels, dtype=float)
    for _, row in corr.iterrows():
        if row["method_id"] not in cfg.methods:
            continue
        mat.loc[row["target"].upper(), cfg.method_label(row["method_id"])] = row[value_col]
    return mat


def _pooled_correlations(cfg: AnalysisConfig) -> dict[str, dict[str, float]]:
    """Pooled pKi vs score across all targets (All targets row)."""
    pooled: dict[str, dict[str, float]] = {}
    frames = []
    for target in cfg.targets:
        try:
            frames.append(load_merged(cfg.merged_dir, target))
        except FileNotFoundError:
            continue
    if not frames:
        return pooled

    merged = pd.concat(frames, ignore_index=True)
    for method_id in cfg.methods:
        metric = cfg.primary_metric(method_id)
        if metric not in merged.columns:
            continue
        x, y = prepare_xy(merged, metric, cfg.exp_col)
        label = cfg.method_label(method_id)
        if len(x) < 3:
            pooled[label] = {"pearson_r": np.nan, "spearman_rho": np.nan, "kendall_tau": np.nan}
            continue
        pr, _ = pearsonr(x, y)
        sr, _ = spearmanr(x, y)
        kt, _ = kendalltau(x, y)
        pooled[label] = {
            "pearson_r": float(pr),
            "spearman_rho": float(sr),
            "kendall_tau": float(kt),
        }
    return pooled


def plot_correlation_heatmaps(corr: pd.DataFrame, cfg: AnalysisConfig) -> None:
    apply_style(cfg)
    cmap = diverging_cmap(cfg)
    out_dir = cfg.figures_dir / "correlations" / "heatmaps"
    pooled = _pooled_correlations(cfg)

    specs = [
        ("pearson_r", "Pearson r"),
        ("spearman_rho", "Spearman rho"),
        ("kendall_tau", "Kendall tau"),
    ]
    for col, title in specs:
        body = _pivot(corr, col, cfg)
        if body.isna().all().all():
            continue

        summary = pd.Series(
            {m: pooled.get(m, {}).get(col, np.nan) for m in body.columns},
            dtype=float,
        )
        plot_heatmap_with_summary_row(
            summary,
            body.astype(float),
            title=f"{title} by target and method",
            cbar_label="Correlation coefficient",
            cmap=cmap,
            vmin=-1.0,
            vmax=1.0,
            center=0.0,
            fmt=".2f",
            path_base=out_dir / f"heatmap_{col}",
        )
