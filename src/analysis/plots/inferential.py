"""Inferential statistics plots for pairwise Wilcoxon comparisons."""

from __future__ import annotations

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

from ..config import AnalysisConfig
from .style import apply_style, save_figure, diverging_cmap, sequential_cmap


FAMILY_ORDER: tuple[str, ...] = (
    "Pearson r",
    "nEF active (top 10%)",
    "Spearman rho",
    "Kendall tau",
    "nEF inactive (top 10%)",
)
CORR_FAMILIES: tuple[str, ...] = (
    "Pearson r",
    "Spearman rho",
    "Kendall tau",
)
NEF_FAMILIES: tuple[str, ...] = (
    "nEF active (top 10%)",
    "nEF inactive (top 10%)",
)
PER_TARGET_FAMILY_ORDER: tuple[str, ...] = (
    "Pearson r (per-target)",
    "Spearman rho (per-target)",
    "nEF10 active (per-target)",
    "nEF10 inactive (per-target)",
)


def _method_order(tests: pd.DataFrame, cfg: AnalysisConfig) -> list[str]:
    ordered = [m for m in cfg.methods if m in set(tests["method_a"]) | set(tests["method_b"])]
    extras = sorted((set(tests["method_a"]) | set(tests["method_b"])) - set(ordered))
    return ordered + extras


def _family_matrix(
    tests: pd.DataFrame,
    family: str,
    methods: list[str],
    value_col: str,
) -> pd.DataFrame:
    mat = pd.DataFrame(np.nan, index=methods, columns=methods, dtype=float)
    np.fill_diagonal(mat.values, 0.0)

    sub = tests[tests["family"] == family]
    for _, row in sub.iterrows():
        a = row["method_a"]
        b = row["method_b"]
        v = float(row[value_col])
        mat.loc[a, b] = v
        mat.loc[b, a] = -v if value_col == "median_diff_a_minus_b" else v
    return mat


def _annot_with_sig(effect: pd.DataFrame, pvals: pd.DataFrame, alpha: float = 0.05) -> np.ndarray:
    annot = np.empty(effect.shape, dtype=object)
    for i in range(effect.shape[0]):
        for j in range(effect.shape[1]):
            val = effect.iat[i, j]
            if pd.isna(val):
                annot[i, j] = ""
                continue
            if i == j:
                annot[i, j] = "0.00"
                continue
            star = "*" if pd.notna(pvals.iat[i, j]) and pvals.iat[i, j] < alpha else ""
            annot[i, j] = f"{val:.2f}{star}"
    return annot


def plot_inferential_heatmaps(tests: pd.DataFrame, cfg: AnalysisConfig) -> None:
    """Create SI-ready heatmaps of pairwise effect sizes and Holm-adjusted p-values."""
    if tests.empty:
        return

    apply_style(cfg)
    out_dir = cfg.figures_dir / "inferential"
    out_dir.mkdir(parents=True, exist_ok=True)

    methods = _method_order(tests, cfg)
    labels = {m: cfg.method_label(m) for m in methods}
    families = [f for f in FAMILY_ORDER if f in set(tests["family"])]
    if not families:
        families = sorted(set(tests["family"]))

    vmax = 0.0
    effect_mats: dict[str, pd.DataFrame] = {}
    pval_mats: dict[str, pd.DataFrame] = {}
    for family in families:
        em = _family_matrix(tests, family, methods, "median_diff_a_minus_b")
        pm = _family_matrix(tests, family, methods, "p_holm")
        effect_mats[family] = em
        pval_mats[family] = pm
        finite_vals = np.abs(em.values[np.isfinite(em.values)])
        if finite_vals.size:
            vmax = max(vmax, float(np.nanmax(finite_vals)))
    if vmax == 0.0:
        vmax = 0.5

    corr_families = [f for f in CORR_FAMILIES if f in set(tests["family"])]
    nef_families = [f for f in NEF_FAMILIES if f in set(tests["family"])]
    if len(corr_families) < 3:
        corr_families = [f for f in families if "nEF" not in f][:3]
    if len(nef_families) < 2:
        nef_families = [f for f in families if "nEF" in f][:2]

    # (1) Correlation families in one horizontal row, single scale after third panel.
    fig1 = plt.figure(figsize=(16.0, 4.8))
    gs1 = fig1.add_gridspec(1, 4, width_ratios=[1.0, 1.0, 1.0, 0.06], wspace=0.25)
    cmap_eff = diverging_cmap(cfg, name="inferential_effect")
    for idx, family in enumerate(corr_families[:3]):
        ax = fig1.add_subplot(gs1[0, idx])
        em = effect_mats[family].rename(index=labels, columns=labels)
        pm = pval_mats[family].rename(index=labels, columns=labels)
        cbar_ax = fig1.add_subplot(gs1[0, 3]) if idx == 2 else None
        sns.heatmap(
            em,
            ax=ax,
            cmap=cmap_eff,
            center=0.0,
            vmin=-vmax,
            vmax=vmax,
            annot=_annot_with_sig(em, pm),
            fmt="",
            linewidths=0.5,
            linecolor="white",
            cbar=(idx == 2),
            cbar_ax=cbar_ax,
            cbar_kws={"label": "Median paired difference (A - B)"},
        )
        ax.set_title(family, fontsize=10)
        ax.tick_params(axis="x", rotation=45)
        ax.tick_params(axis="y", rotation=0)
    fig1.suptitle("Pairwise Wilcoxon effects: Pearson / Spearman / Kendall", y=1.02)
    save_figure(fig1, out_dir / "heatmap_wilcoxon_effects")

    # (2) nEF families in one horizontal row, scale after first panel.
    fig_nef = plt.figure(figsize=(11.5, 4.8))
    gs_nef = fig_nef.add_gridspec(1, 3, width_ratios=[1.0, 0.06, 1.0], wspace=0.30)
    for idx, family in enumerate(nef_families[:2]):
        ax = fig_nef.add_subplot(gs_nef[0, 0 if idx == 0 else 2])
        em = effect_mats[family].rename(index=labels, columns=labels)
        pm = pval_mats[family].rename(index=labels, columns=labels)
        cbar_ax = fig_nef.add_subplot(gs_nef[0, 1]) if idx == 0 else None
        sns.heatmap(
            em,
            ax=ax,
            cmap=cmap_eff,
            center=0.0,
            vmin=-vmax,
            vmax=vmax,
            annot=_annot_with_sig(em, pm),
            fmt="",
            linewidths=0.5,
            linecolor="white",
            cbar=(idx == 0),
            cbar_ax=cbar_ax,
            cbar_kws={"label": "Median paired difference (A - B)"},
        )
        ax.set_title(family, fontsize=10)
        ax.tick_params(axis="x", rotation=45)
        ax.tick_params(axis="y", rotation=0)
    fig_nef.suptitle(r"Pairwise Wilcoxon effects: nEF$_{10,\mathrm{active}}$ / nEF$_{10,\mathrm{inactive}}$", y=1.02)
    save_figure(fig_nef, out_dir / "heatmap_wilcoxon_effects_nef10")

    n = len(families)
    ncols = 2 if n > 1 else 1
    nrows = int(np.ceil(n / ncols))
    fig2, axes2 = plt.subplots(nrows, ncols, figsize=(6.0 * ncols, 4.8 * nrows), squeeze=False)
    cmap_p = sequential_cmap(cfg, name="inferential_logp")
    for idx, family in enumerate(families):
        r, c = divmod(idx, ncols)
        ax = axes2[r][c]
        pm = pval_mats[family]
        logp = -np.log10(pm.clip(lower=1e-12)).rename(index=labels, columns=labels)
        sns.heatmap(
            logp,
            ax=ax,
            cmap=cmap_p,
            vmin=0.0,
            vmax=max(3.0, float(np.nanmax(logp.values[np.isfinite(logp.values)])) if np.isfinite(logp.values).any() else 3.0),
            annot=True,
            fmt=".2f",
            linewidths=0.5,
            linecolor="white",
            cbar=(idx == 0),
            cbar_kws={"label": r"$-\log_{10}(p_{\mathrm{Holm}})$"},
        )
        ax.set_title(family, fontsize=10)
        ax.tick_params(axis="x", rotation=45)
        ax.tick_params(axis="y", rotation=0)

    for idx in range(len(families), nrows * ncols):
        r, c = divmod(idx, ncols)
        axes2[r][c].axis("off")
    fig2.suptitle(r"Pairwise Wilcoxon significance maps ($-\log_{10}(p_{\mathrm{Holm}})$)", y=1.01)
    save_figure(fig2, out_dir / "heatmap_wilcoxon_logp_holm")


def plot_per_target_inferential_heatmaps(per_target_tests: pd.DataFrame, cfg: AnalysisConfig) -> None:
    """Create per-target pairwise effect and significance heatmaps."""
    if per_target_tests.empty:
        return

    apply_style(cfg)
    out_dir = cfg.figures_dir / "inferential"
    out_dir.mkdir(parents=True, exist_ok=True)

    targets = [t.lower() for t in cfg.targets]
    fams = [f for f in PER_TARGET_FAMILY_ORDER if f in set(per_target_tests["family"])]
    if not fams:
        fams = sorted(per_target_tests["family"].unique())
    all_comparisons = sorted(per_target_tests["comparison"].unique())

    vmax = float(np.nanmax(np.abs(per_target_tests["effect_obs"].to_numpy(dtype=float))))
    if not np.isfinite(vmax) or vmax == 0.0:
        vmax = 0.5

    cmap_eff = diverging_cmap(cfg, name="per_target_effect")
    cmap_logp = sequential_cmap(cfg, name="per_target_logp")
    vmax_logp = float(
        np.nanmax(
            -np.log10(per_target_tests["perm_p_bh_target"].clip(lower=1e-12).to_numpy(dtype=float))
        )
    )
    if not np.isfinite(vmax_logp):
        vmax_logp = 3.0
    vmax_logp = max(3.0, vmax_logp)

    for target in targets:
        sub_t = per_target_tests[per_target_tests["target"] == target].copy()
        if sub_t.empty:
            continue

        eff = sub_t.pivot(index="family", columns="comparison", values="effect_obs")
        pbh = sub_t.pivot(index="family", columns="comparison", values="perm_p_bh_target")
        eff = eff.reindex(index=fams, columns=all_comparisons)
        pbh = pbh.reindex(index=fams, columns=all_comparisons)
        logp = -np.log10(pbh.clip(lower=1e-12))

        fig_e, ax_e = plt.subplots(1, 1, figsize=(max(9.0, 0.40 * len(all_comparisons)), 3.8))
        sns.heatmap(
            eff,
            ax=ax_e,
            cmap=cmap_eff,
            center=0.0,
            vmin=-vmax,
            vmax=vmax,
            annot=True,
            fmt=".2f",
            linewidths=0.4,
            linecolor="white",
            cbar=True,
            cbar_kws={"label": "Effect size (A - B)"},
        )
        ax_e.set_title(f"{target.upper()}: per-target pairwise effects", fontsize=10)
        ax_e.set_xlabel("Method pair")
        ax_e.set_ylabel("Family")
        ax_e.tick_params(axis="x", rotation=45)
        ax_e.tick_params(axis="y", rotation=0)
        save_figure(fig_e, out_dir / f"per_target_{target}_wilcoxon_effects")

        fig_p, ax_p = plt.subplots(1, 1, figsize=(max(9.0, 0.40 * len(all_comparisons)), 3.8))
        sns.heatmap(
            logp,
            ax=ax_p,
            cmap=cmap_logp,
            vmin=0.0,
            vmax=vmax_logp,
            annot=True,
            fmt=".2f",
            linewidths=0.4,
            linecolor="white",
            cbar=True,
            cbar_kws={"label": r"$-\log_{10}(p_{\mathrm{BH,target}})$"},
        )
        ax_p.set_title(f"{target.upper()}: per-target pairwise significance", fontsize=10)
        ax_p.set_xlabel("Method pair")
        ax_p.set_ylabel("Family")
        ax_p.tick_params(axis="x", rotation=45)
        ax_p.tick_params(axis="y", rotation=0)
        save_figure(fig_p, out_dir / f"per_target_{target}_wilcoxon_logp")
