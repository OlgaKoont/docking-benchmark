"""Figure styling: palette, DPI, save helpers, heatmap layout with summary row/column."""

from __future__ import annotations

from decimal import Decimal, ROUND_HALF_UP
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.gridspec import GridSpec

from ..config import AnalysisConfig

ALL_TARGETS_LABEL = "All targets"
GAP_RATIO = 0.10

# Article diverging palette: blue (negative) → gray/white (0) → salmon (positive).
# Positions are in colormap space [0, 1] matching data range [-1, +1] with center=0.
CORRELATION_DIVERGING_STOPS: list[tuple[float, str]] = [
    (0.0, "#6584E1"),    # -1.0  saturated blue
    (0.25, "#A2BFFF"),   # -0.5  periwinkle
    (0.375, "#B7CAEA"),  # -0.25 light blue
    (0.5, "#DEDAD7"),    #  0.0  light gray / near-white
    (0.625, "#F0D1BF"),  # +0.25 peach
    (0.75, "#FDBBAB"),   # +0.5  light salmon
    (0.875, "#F18C6E"),  # +0.75 salmon-orange
    (1.0, "#E5885F"),    # +1.0  saturated salmon
]

# Sequential variant for bounded positive metrics (nEF 0–1, pass rates 0–100).
POSITIVE_SEQUENTIAL_STOPS: list[tuple[float, str]] = [
    (0.0, "#6584E1"),
    (0.35, "#B7CAEA"),
    (0.55, "#DEDAD7"),
    (0.75, "#FDBBAB"),
    (1.0, "#E5885F"),
]


def round_half_up(value: float, digits: int = 2) -> float:
    quant = Decimal("1").scaleb(-digits)
    return float(Decimal(str(value)).quantize(quant, rounding=ROUND_HALF_UP))


def correlation_diverging_cmap(
    cfg: AnalysisConfig | None = None, name: str = "article_correlation"
) -> LinearSegmentedColormap:
    """Salmon (positive) — gray (0) — blue (negative), article heatmap palette."""
    _ = cfg  # reserved for future per-run overrides
    return LinearSegmentedColormap.from_list(name, CORRELATION_DIVERGING_STOPS, N=256)


def positive_sequential_cmap(
    cfg: AnalysisConfig | None = None, name: str = "article_positive"
) -> LinearSegmentedColormap:
    """Low values blue → high values salmon (nEF, PoseBusters pass rates)."""
    _ = cfg
    return LinearSegmentedColormap.from_list(name, POSITIVE_SEQUENTIAL_STOPS, N=256)


def diverging_cmap(cfg: AnalysisConfig, name: str = "article_diverging") -> LinearSegmentedColormap:
    return correlation_diverging_cmap(cfg, name=name)


def sequential_cmap(cfg: AnalysisConfig, name: str = "article_sequential") -> LinearSegmentedColormap:
    return positive_sequential_cmap(cfg, name=name)


def apply_style(cfg: AnalysisConfig) -> None:
    plt.rcParams.update(
        {
            "font.family": cfg.font_family,
            "font.size": 9,
            "axes.labelsize": 10,
            "axes.titlesize": 10,
            "xtick.labelsize": 9,
            "ytick.labelsize": 9,
            "figure.dpi": cfg.figure_dpi,
            "savefig.dpi": cfg.figure_dpi,
            "svg.fonttype": "none",
        }
    )


def save_figure(fig: plt.Figure, path_base: Path) -> None:
    path_base.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path_base.with_suffix(".png"), bbox_inches="tight")
    fig.savefig(path_base.with_suffix(".svg"), bbox_inches="tight")
    plt.close(fig)


def method_color(cfg: AnalysisConfig, method_id: str, default: str = "#555555") -> str:
    return cfg.method_colors.get(method_id, default)


def _format_annot(mat: pd.DataFrame, fmt: str) -> np.ndarray:
    annot = np.empty(mat.shape, dtype=object)
    for i in range(mat.shape[0]):
        for j in range(mat.shape[1]):
            val = mat.iat[i, j]
            if pd.isna(val):
                annot[i, j] = ""
            elif fmt == ".2f":
                annot[i, j] = f"{round_half_up(float(val), 2):.2f}"
            elif fmt == ".1f":
                annot[i, j] = f"{round_half_up(float(val), 1):.1f}"
            else:
                annot[i, j] = f"{float(val):{fmt}}"
    return annot


def plot_heatmap_with_summary_row(
    summary_row: pd.Series,
    body: pd.DataFrame,
    *,
    title: str,
    cbar_label: str,
    cmap: LinearSegmentedColormap,
    vmin: float | None,
    vmax: float | None,
    center: float | None,
    fmt: str,
    path_base: Path,
    figsize_scale: tuple[float, float] = (0.55, 0.35),
) -> None:
    """Heatmap: one summary row (All targets) + gap + per-target rows."""
    summary_df = summary_row.to_frame().T
    summary_df.index = [ALL_TARGETS_LABEL]
    n_cols = len(body.columns)
    n_body = len(body)

    fig_w = max(6.5, figsize_scale[0] * n_cols)
    fig_h = max(4.5, figsize_scale[1] * (1 + n_body) + 0.8)
    fig = plt.figure(figsize=(fig_w, fig_h))
    gs = GridSpec(
        2, 2,
        figure=fig,
        height_ratios=[1.0, max(n_body, 1)],
        width_ratios=[1.0, 0.045],
        hspace=GAP_RATIO,
        wspace=0.05,
    )
    ax_top = fig.add_subplot(gs[0, 0])
    ax_bot = fig.add_subplot(gs[1, 0])
    cbar_ax = fig.add_subplot(gs[:, 1])

    heatmap_kw = dict(
        cmap=cmap,
        vmin=vmin,
        vmax=vmax,
        center=center,
        linewidths=0.5,
        linecolor="white",
        annot_kws={"size": 9},
        cbar_ax=cbar_ax,
        cbar_kws={"label": cbar_label},
    )

    sns.heatmap(
        summary_df,
        ax=ax_top,
        annot=_format_annot(summary_df, fmt),
        fmt="",
        cbar=False,
        **heatmap_kw,
    )
    sns.heatmap(
        body,
        ax=ax_bot,
        annot=_format_annot(body, fmt) if fmt else None,
        fmt="" if fmt else "",
        **heatmap_kw,
    )

    ax_top.set_ylabel("")
    ax_top.set_xlabel("")
    ax_top.tick_params(axis="x", labelbottom=False)
    ax_bot.set_xlabel("Docking method")
    ax_bot.set_ylabel("Target")
    ax_bot.tick_params(axis="x", rotation=45)
    fig.suptitle(title, y=1.02, fontsize=11)
    save_figure(fig, path_base)


def plot_heatmap_with_summary_col(
    summary_col: pd.Series,
    body: pd.DataFrame,
    *,
    title: str,
    cbar_label: str,
    cmap: LinearSegmentedColormap,
    vmin: float,
    vmax: float,
    fmt: str,
    path_base: Path,
    figsize_scale: tuple[float, float] = (0.45, 0.22),
) -> None:
    """Heatmap: one summary column (All targets) + gap + per-target columns."""
    summary_df = summary_col.to_frame()
    summary_df.columns = [ALL_TARGETS_LABEL]
    n_rows = len(body.index)
    n_cols = len(body.columns)

    fig_w = max(8.0, figsize_scale[0] * (1 + n_cols) + 1.2)
    fig_h = max(6.0, figsize_scale[1] * n_rows)
    fig = plt.figure(figsize=(fig_w, fig_h))
    gs = GridSpec(
        1, 2,
        figure=fig,
        width_ratios=[1.35, max(n_cols, 1)],
        wspace=GAP_RATIO,
    )
    gs_inner = gs[1].subgridspec(1, 2, width_ratios=[1.0, 0.045], wspace=0.05)
    ax_left = fig.add_subplot(gs[0, 0])
    ax_mid = fig.add_subplot(gs_inner[0, 0])
    cbar_ax = fig.add_subplot(gs_inner[0, 1])

    heatmap_kw = dict(
        cmap=cmap,
        vmin=vmin,
        vmax=vmax,
        linewidths=0.3,
        linecolor="white",
        cbar_ax=cbar_ax,
        cbar_kws={"label": cbar_label},
    )

    sns.heatmap(
        summary_df,
        ax=ax_left,
        annot=_format_annot(summary_df, fmt) if fmt else None,
        fmt="" if fmt else "",
        cbar=False,
        **heatmap_kw,
    )
    sns.heatmap(
        body,
        ax=ax_mid,
        annot=False,
        **heatmap_kw,
    )

    ax_left.set_ylabel("Check")
    ax_left.set_xlabel("")
    ax_left.tick_params(axis="y", labelleft=True, pad=2)
    for tick in ax_left.get_yticklabels():
        tick.set_fontsize(8)
    ax_mid.set_ylabel("")
    ax_mid.set_xlabel("Target")
    ax_mid.tick_params(axis="x", rotation=45)
    # Keep check labels only on the left summary panel.
    ax_mid.tick_params(axis="y", left=False, labelleft=False)
    fig.suptitle(title, y=1.02, fontsize=11)
    save_figure(fig, path_base)
