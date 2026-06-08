"""Inferential statistics: bootstrap, permutation, Wilcoxon, multiple testing."""

from __future__ import annotations

import warnings
from itertools import combinations

import numpy as np
import pandas as pd
from scipy.stats import ConstantInputWarning, pearsonr, spearmanr, wilcoxon

from .config import AnalysisConfig
from .data import ensure_pki, load_merged
from .enrichment import compute_ef, detect_score_direction


def _correlation_coeff(x: np.ndarray, y: np.ndarray, method: str) -> float:
    """Pearson or Spearman r; NaN if undefined (constant bootstrap sample)."""
    if len(x) < 3:
        return np.nan
    if np.std(x, ddof=1) == 0 or np.std(y, ddof=1) == 0:
        return np.nan
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", ConstantInputWarning)
        if method == "pearson":
            r = pearsonr(x, y)[0]
        else:
            r = spearmanr(x, y)[0]
    return float(r) if np.isfinite(r) else np.nan


def benjamini_hochberg(pvals: np.ndarray) -> np.ndarray:
    p = np.asarray(pvals, dtype=float)
    m = len(p)
    if m == 0:
        return p
    order = np.argsort(p)
    ranked = p[order]
    out = np.empty(m)
    prev = 1.0
    for i in range(m - 1, -1, -1):
        rank = i + 1
        v = min(ranked[i] * m / rank, prev)
        prev = v
        out[order[i]] = min(v, 1.0)
    return out


def holm(pvals: np.ndarray) -> np.ndarray:
    p = np.asarray(pvals, dtype=float)
    m = len(p)
    if m == 0:
        return p
    order = np.argsort(p)
    out = np.empty(m)
    prev = 0.0
    for i, idx in enumerate(order):
        v = (m - i) * p[idx]
        v = max(v, prev)
        prev = v
        out[idx] = min(v, 1.0)
    return out


def bootstrap_ci(
    x: np.ndarray,
    y: np.ndarray,
    method: str = "pearson",
    n_boot: int = 10_000,
    seed: int = 42,
) -> tuple[float, float, float]:
    rng = np.random.default_rng(seed)
    n = len(x)
    vals: list[float] = []
    for _ in range(n_boot):
        idx = rng.integers(0, n, n)
        xb, yb = x[idx], y[idx]
        r = _correlation_coeff(xb, yb, method)
        if np.isfinite(r):
            vals.append(r)
    if not vals:
        return np.nan, np.nan, np.nan
    arr = np.array(vals)
    return float(np.percentile(arr, 2.5)), float(np.percentile(arr, 97.5)), float(np.mean(arr))


def permutation_p(
    x: np.ndarray,
    y: np.ndarray,
    method: str = "pearson",
    n_perm: int = 10_000,
    seed: int = 43,
) -> float:
    r_obs_val = _correlation_coeff(x, y, method)
    if not np.isfinite(r_obs_val):
        return np.nan
    r_obs = abs(r_obs_val)
    rng = np.random.default_rng(seed)
    n = len(x)
    count = 0
    for _ in range(n_perm):
        yp = y[rng.permutation(n)]
        r = _correlation_coeff(x, yp, method)
        if np.isfinite(r) and abs(r) >= r_obs:
            count += 1
    return (count + 1) / (n_perm + 1)


def wilcoxon_paired(
    x: list[float], y: list[float]
) -> tuple[float, int, int, float, float]:
    """
    Paired Wilcoxon signed-rank (two-sided).

    Returns:
        p_raw, n_nonzero_diff, n_targets, median_diff (A - B), mean_diff (A - B)

    n_nonzero_diff excludes target pairs with |diff| <= 1e-12 (Wilcoxon convention).
    For nEF this often happens when both methods reach the ceiling nEF = 1.0.
    """
    xa = np.array(x, dtype=float)
    ya = np.array(y, dtype=float)
    n_targets = len(xa)
    if n_targets == 0:
        return np.nan, 0, 0, np.nan, np.nan

    diff = xa - ya
    median_diff = float(np.median(diff))
    mean_diff = float(np.mean(diff))

    nz = diff[np.abs(diff) > 1e-12]
    n_nonzero = len(nz)
    if n_nonzero < 1:
        return np.nan, 0, n_targets, median_diff, mean_diff

    p = wilcoxon(
        nz, zero_method="wilcox", alternative="two-sided", correction=False, method="auto"
    ).pvalue
    return float(p), n_nonzero, n_targets, median_diff, mean_diff


def _append_wilcoxon_row(
    rows: list[dict],
    *,
    family: str,
    label1: str,
    label2: str,
    m1: str,
    m2: str,
    values_a: list[float],
    values_b: list[float],
) -> None:
    p, n_nz, n_tgt, med, mean = wilcoxon_paired(values_a, values_b)
    rows.append(
        {
            "family": family,
            "comparison": f"{label1} vs {label2}",
            "method_a": m1,
            "method_b": m2,
            "n_targets": n_tgt,
            "n_nonzero_diff": n_nz,
            "median_diff_a_minus_b": med,
            "mean_diff_a_minus_b": mean,
            "p_raw": p,
        }
    )


def run_pairwise_method_tests(
    corr: pd.DataFrame,
    ef: pd.DataFrame,
    methods: list[str],
    method_labels: dict[str, str],
    primary_metrics: dict[str, str],
) -> pd.DataFrame:
    rows: list[dict] = []

    corr_cols = [
        ("Pearson r", "pearson_r"),
        ("Spearman rho", "spearman_rho"),
        ("Kendall tau", "kendall_tau"),
    ]

    for m1, m2 in combinations(methods, 2):
        metric1 = primary_metrics[m1]
        metric2 = primary_metrics[m2]
        label1 = method_labels[m1]
        label2 = method_labels[m2]

        for family, col in corr_cols:
            a, b = [], []
            for target in sorted(corr["target"].unique()):
                v1 = corr[(corr.target == target) & (corr.metric == metric1)][col]
                v2 = corr[(corr.target == target) & (corr.metric == metric2)][col]
                if len(v1) and len(v2) and pd.notna(v1.iloc[0]) and pd.notna(v2.iloc[0]):
                    a.append(float(v1.iloc[0]))
                    b.append(float(v2.iloc[0]))
            _append_wilcoxon_row(
                rows,
                family=family,
                label1=label1,
                label2=label2,
                m1=m1,
                m2=m2,
                values_a=a,
                values_b=b,
            )

        for frac in ("1", "5", "10"):
            for set_name, family in (("active", "nEF active"), ("low", "nEF inactive")):
                col = f"nEF{frac}_{set_name}"
                a, b = [], []
                for target in sorted(ef["target"].unique()):
                    v1 = ef[(ef.target == target) & (ef.method_id == m1)][col]
                    v2 = ef[(ef.target == target) & (ef.method_id == m2)][col]
                    if len(v1) and len(v2) and pd.notna(v1.iloc[0]) and pd.notna(v2.iloc[0]):
                        a.append(float(v1.iloc[0]))
                        b.append(float(v2.iloc[0]))
                _append_wilcoxon_row(
                    rows,
                    family=f"{family} (top {frac}%)",
                    label1=label1,
                    label2=label2,
                    m1=m1,
                    m2=m2,
                    values_a=a,
                    values_b=b,
                )

    result = pd.DataFrame(rows)
    if result.empty:
        return result

    for family in result["family"].unique():
        mask = result["family"] == family
        pvals = result.loc[mask, "p_raw"].fillna(1.0).values
        result.loc[mask, "p_holm"] = holm(pvals)
        result.loc[mask, "p_bh"] = benjamini_hochberg(pvals)
    return result


def _classes_from_df(df: pd.DataFrame) -> pd.Series:
    if "pValue" in df.columns:
        pki_vals = pd.to_numeric(df["pValue"], errors="coerce")
        return pd.Series(
            np.where(
                pki_vals >= 7.0,
                "high",
                np.where(pki_vals > 6.0, "medium", "low"),
            ),
            index=df.index,
        ).astype(str)
    ki_vals = pd.to_numeric(df.get("standard_value"), errors="coerce")
    return pd.Series(
        np.where(
            ki_vals <= 100.0,
            "high",
            np.where(ki_vals < 1000.0, "medium", "low"),
        ),
        index=df.index,
    ).astype(str)


def _paired_target_frame(
    df: pd.DataFrame,
    metric_a: str,
    metric_b: str,
    exp_col: str,
) -> pd.DataFrame:
    pki = ensure_pki(df, exp_col)
    sa = pd.to_numeric(df[metric_a], errors="coerce")
    sb = pd.to_numeric(df[metric_b], errors="coerce")
    classes = _classes_from_df(df)
    out = pd.DataFrame({"pki": pki, "score_a": sa, "score_b": sb, "class": classes})
    out = out[np.isfinite(out["pki"]) & np.isfinite(out["score_a"]) & np.isfinite(out["score_b"])]
    return out.reset_index(drop=True)


def _corr_diff(
    pki: np.ndarray,
    s1: np.ndarray,
    s2: np.ndarray,
    method: str,
) -> float:
    r1 = _correlation_coeff(pki, s1, method)
    r2 = _correlation_coeff(pki, s2, method)
    if not np.isfinite(r1) or not np.isfinite(r2):
        return np.nan
    return float(r1 - r2)


def _nef_diff(
    classes: np.ndarray,
    scores_a: np.ndarray,
    scores_b: np.ndarray,
    metric_a: str,
    metric_b: str,
    set_name: str,
) -> float:
    is_active = classes == ("low" if set_name == "low" else "high")
    if set_name == "active":
        is_active = np.isin(classes, ["high", "medium"])
    direction_a = detect_score_direction(metric_a, set_name)
    direction_b = detect_score_direction(metric_b, set_name)
    _, n_ef_a, *_ = compute_ef(is_active, scores_a, 0.10, direction_a)
    _, n_ef_b, *_ = compute_ef(is_active, scores_b, 0.10, direction_b)
    if not np.isfinite(n_ef_a) or not np.isfinite(n_ef_b):
        return np.nan
    return float(n_ef_a - n_ef_b)


def run_per_target_pairwise_tests(
    cfg: AnalysisConfig,
    methods: list[str],
    method_labels: dict[str, str],
    primary_metrics: dict[str, str],
) -> pd.DataFrame:
    """
    Per-target pairwise method comparisons on shared ligand subsets.

    Endpoints:
      - Pearson difference: r(method_a) - r(method_b)
      - Spearman difference: rho(method_a) - rho(method_b)
      - nEF10 active difference
      - nEF10 inactive difference
    """
    rows: list[dict] = []
    n_boot = min(int(cfg.n_bootstrap), 2000)
    n_perm = min(int(cfg.n_permutation), 2000)

    for target in cfg.targets:
        df = load_merged(cfg.merged_dir, target)
        for m1, m2 in combinations(methods, 2):
            metric1 = primary_metrics[m1]
            metric2 = primary_metrics[m2]
            if metric1 not in df.columns or metric2 not in df.columns:
                continue

            pf = _paired_target_frame(df, metric1, metric2, cfg.exp_col)
            if len(pf) < 20:
                continue

            pki = pf["pki"].to_numpy(dtype=float)
            s1 = pf["score_a"].to_numpy(dtype=float)
            s2 = pf["score_b"].to_numpy(dtype=float)
            classes = pf["class"].to_numpy(dtype=str)

            corr_specs = (
                ("Pearson r (per-target)", "pearson"),
                ("Spearman rho (per-target)", "spearman"),
            )
            for family, corr_kind in corr_specs:
                d_obs = _corr_diff(pki, s1, s2, corr_kind)
                rng_boot = np.random.default_rng(cfg.random_seed)
                boot_vals: list[float] = []
                n = len(pf)
                for _ in range(n_boot):
                    idx = rng_boot.integers(0, n, n)
                    d = _corr_diff(pki[idx], s1[idx], s2[idx], corr_kind)
                    if np.isfinite(d):
                        boot_vals.append(float(d))
                if boot_vals:
                    arr = np.asarray(boot_vals, dtype=float)
                    ci_lo = float(np.percentile(arr, 2.5))
                    ci_hi = float(np.percentile(arr, 97.5))
                    d_boot_mean = float(np.mean(arr))
                else:
                    ci_lo = np.nan
                    ci_hi = np.nan
                    d_boot_mean = np.nan

                if np.isfinite(d_obs):
                    rng_perm = np.random.default_rng(cfg.random_seed + 7)
                    obs_abs = abs(float(d_obs))
                    count = 0
                    for _ in range(n_perm):
                        swap = rng_perm.random(n) < 0.5
                        pa = s1.copy()
                        pb = s2.copy()
                        pa[swap], pb[swap] = pb[swap], pa[swap]
                        d = _corr_diff(pki, pa, pb, corr_kind)
                        if np.isfinite(d) and abs(float(d)) >= obs_abs:
                            count += 1
                    p_perm = (count + 1) / (n_perm + 1)
                else:
                    p_perm = np.nan
                rows.append(
                    {
                        "target": target.lower(),
                        "family": family,
                        "comparison": f"{method_labels[m1]} vs {method_labels[m2]}",
                        "method_a": m1,
                        "method_b": m2,
                        "n_points": int(len(pf)),
                        "effect_obs": d_obs,
                        "effect_boot_mean": d_boot_mean,
                        "ci95_low": ci_lo,
                        "ci95_high": ci_hi,
                        "perm_p_raw": p_perm,
                    }
                )

            nef_specs = (
                ("nEF10 active (per-target)", "active"),
                ("nEF10 inactive (per-target)", "low"),
            )
            for family, set_name in nef_specs:
                d_obs = _nef_diff(classes, s1, s2, metric1, metric2, set_name)
                rng_boot = np.random.default_rng(cfg.random_seed + 1)
                boot_vals: list[float] = []
                n = len(pf)
                for _ in range(n_boot):
                    idx = rng_boot.integers(0, n, n)
                    d = _nef_diff(classes[idx], s1[idx], s2[idx], metric1, metric2, set_name)
                    if np.isfinite(d):
                        boot_vals.append(float(d))
                if boot_vals:
                    arr = np.asarray(boot_vals, dtype=float)
                    ci_lo = float(np.percentile(arr, 2.5))
                    ci_hi = float(np.percentile(arr, 97.5))
                    d_boot_mean = float(np.mean(arr))
                else:
                    ci_lo = np.nan
                    ci_hi = np.nan
                    d_boot_mean = np.nan

                if np.isfinite(d_obs):
                    rng_perm = np.random.default_rng(cfg.random_seed + 11)
                    obs_abs = abs(float(d_obs))
                    count = 0
                    for _ in range(n_perm):
                        swap = rng_perm.random(n) < 0.5
                        pa = s1.copy()
                        pb = s2.copy()
                        pa[swap], pb[swap] = pb[swap], pa[swap]
                        d = _nef_diff(classes, pa, pb, metric1, metric2, set_name)
                        if np.isfinite(d) and abs(float(d)) >= obs_abs:
                            count += 1
                    p_perm = (count + 1) / (n_perm + 1)
                else:
                    p_perm = np.nan
                rows.append(
                    {
                        "target": target.lower(),
                        "family": family,
                        "comparison": f"{method_labels[m1]} vs {method_labels[m2]}",
                        "method_a": m1,
                        "method_b": m2,
                        "n_points": int(len(pf)),
                        "effect_obs": d_obs,
                        "effect_boot_mean": d_boot_mean,
                        "ci95_low": ci_lo,
                        "ci95_high": ci_hi,
                        "perm_p_raw": p_perm,
                    }
                )

    result = pd.DataFrame(rows)
    if result.empty:
        return result

    result["perm_p_bh_target"] = np.nan
    for (target, family), sub in result.groupby(["target", "family"]):
        idx = sub.index
        result.loc[idx, "perm_p_bh_target"] = benjamini_hochberg(
            sub["perm_p_raw"].fillna(1.0).to_numpy(dtype=float)
        )
    return result
