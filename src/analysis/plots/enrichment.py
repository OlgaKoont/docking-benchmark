"""nEF violin figures (article style, active + inactive panels)."""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from ..config import AnalysisConfig
from .style import apply_style, positive_sequential_cmap, save_figure


def _panel_plot(
    ax: plt.Axes,
    df: pd.DataFrame,
    col: str,
    title: str,
    panel_letter: str,
    order: list[str],
    cmap,
    seed: int,
    frac: int,
) -> None:
    rng = np.random.default_rng(seed)
    x_positions = np.arange(len(order))
    values_by_method = [
        df.loc[df["method_label"] == method, col].astype(float).dropna().to_numpy()
        for method in order
    ]

    violin = ax.violinplot(
        values_by_method,
        positions=x_positions,
        widths=0.82,
        showmeans=False,
        showmedians=False,
        showextrema=False,
    )
    for body in violin["bodies"]:
        body.set_facecolor("#DEDAD7")
        body.set_edgecolor("#2A2A2A")
        body.set_linewidth(0.6)
        body.set_alpha(0.22)

    label_thr = 0.6
    for i, vals in enumerate(values_by_method):
        if len(vals) == 0:
            continue
        method = order[i]
        df_method = (
            df.loc[df["method_label"] == method, ["target", col]]
            .dropna()
            .copy()
            .reset_index(drop=True)
        )
        proteins = df_method["target"].astype(str).str.upper().to_numpy()
        jitter_x = x_positions[i] + rng.uniform(-0.18, 0.18, size=len(vals))
        jitter_y = vals + rng.uniform(-0.015, 0.015, size=len(vals))
        jitter_y = np.clip(jitter_y, 0.0, 1.0)
        low_mask = vals < label_thr
        jitter_x[low_mask] = x_positions[i]
        jitter_y[low_mask] = vals[low_mask]
        ax.scatter(
            jitter_x,
            jitter_y,
            s=34,
            c=[cmap(float(np.clip(v, 0.0, 1.0))) for v in vals],
            alpha=0.82,
            edgecolors="white",
            linewidths=0.55,
            zorder=3,
        )
        label_idx = sorted([j for j, yi in enumerate(vals) if yi < label_thr], key=lambda j: vals[j])
        last_y_by_side = {"right": -10.0, "left": -10.0}
        min_gap = 0.03
        for rank, j in enumerate(label_idx):
            side = "right" if rank % 2 == 0 else "left"
            ha = "left" if side == "right" else "right"
            x_text = (jitter_x[j] + 0.08) if side == "right" else (jitter_x[j] - 0.08)
            y_text = float(jitter_y[j])
            if (y_text - last_y_by_side[side]) < min_gap:
                y_text = last_y_by_side[side] + min_gap
            y_text = float(np.clip(y_text, 0.0, 1.0))
            last_y_by_side[side] = y_text
            ax.annotate(
                proteins[j],
                (jitter_x[j], jitter_y[j]),
                textcoords="data",
                xytext=(x_text, y_text),
                ha=ha,
                va="center",
                fontsize=7,
                alpha=0.9,
                zorder=5,
            )
        ax.scatter(
            [x_positions[i]],
            [float(np.median(vals))],
            s=50,
            marker="D",
            c=["#222222"],
            edgecolors="white",
            linewidths=0.7,
            zorder=4,
        )

    ax.set_xticks(x_positions)
    ax.set_xticklabels(order, fontsize=9)
    ax.set_ylim(-0.03, 1.03)
    ax.set_ylabel(f"nEF{frac}", fontsize=11)
    ax.grid(axis="y", linestyle="--", alpha=0.35, linewidth=0.7)
    ax.axhline(0.5, color="#666666", linestyle=":", linewidth=1.0, alpha=0.7)
    ax.set_title(title, fontsize=12, pad=10)
    ax.text(
        -0.1, 1.02, panel_letter,
        transform=ax.transAxes,
        fontsize=18, fontweight="bold", va="bottom", ha="right",
    )


def plot_nef_violins(ef: pd.DataFrame, cfg: AnalysisConfig) -> None:
    apply_style(cfg)
    out_dir = cfg.figures_dir / "enrichment" / "figures"
    out_dir.mkdir(parents=True, exist_ok=True)

    for old in out_dir.glob("heatmap_nEF*"):
        old.unlink()

    order = [cfg.method_label(m) for m in cfg.methods if m in set(ef["method_id"])]
    point_cmap = positive_sequential_cmap(cfg, name="nef_point_gradient")

    for frac in ("1", "5", "10"):
        active_col = f"nEF{frac}_active"
        low_col = f"nEF{frac}_low"
        if active_col not in ef.columns or low_col not in ef.columns:
            continue

        fig, axes = plt.subplots(1, 2, figsize=(12.8, 5.6), sharey=True)
        _panel_plot(
            axes[0], ef, active_col,
            f"nEF{frac} (actives: Ki < 1000 nM)", "A", order, point_cmap,
            seed=cfg.random_seed, frac=int(frac),
        )
        _panel_plot(
            axes[1], ef, low_col,
            f"nEF{frac} (inactives: Ki ≥ 1000 nM)", "B", order, point_cmap,
            seed=cfg.random_seed + 1, frac=int(frac),
        )
        plt.tight_layout()
        save_figure(fig, out_dir / f"summary_nEF{frac}_combined")
