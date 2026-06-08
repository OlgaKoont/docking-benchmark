"""
EF / nEF — canonical implementation aligned with
ToxAffinity/docking-benchmark-2/scripts/compute_enrichment_from_merged.py

Do not sign-flip docking scores here; ranking direction is handled only via
detect_score_direction (asc/desc).
"""

from __future__ import annotations

from typing import Optional

import numpy as np
import pandas as pd

from .config import AnalysisConfig
from .constants import TOP_FRACS
from .data import load_merged

# ---------------------------------------------------------------------------
# Canonical EF logic (keep in sync with compute_enrichment_from_merged.py)
# ---------------------------------------------------------------------------


def classify_activity(standard_value: float) -> Optional[str]:
    if standard_value is None or (isinstance(standard_value, float) and np.isnan(standard_value)):
        return None
    try:
        v = float(standard_value)
    except Exception:
        return None
    if v <= 100.0:
        return "high"
    if v < 1000.0:
        return "medium"
    return "low"


def detect_score_direction(metric: str, set_name: str) -> str:
    metric_lower = metric.lower()

    if metric_lower == "pki_oracle":
        base_direction = "desc"
    elif "boltz2" in metric_lower:
        base_direction = "asc"
    elif "dynamicbind" in metric_lower:
        base_direction = "desc"
    elif "plapt_affinity" in metric_lower and "um" not in metric_lower:
        base_direction = "desc"
    elif "plapt_affinity_um" in metric_lower:
        base_direction = "asc"
    elif "qvina" in metric_lower:
        base_direction = "asc"
    elif "gnina" in metric_lower:
        if "cnn_affinity" in metric_lower:
            base_direction = "desc"
        else:
            base_direction = "asc"
    else:
        base_direction = "asc"

    if set_name == "low":
        return "desc" if base_direction == "asc" else "asc"
    return base_direction


def compute_ef(
    is_active: np.ndarray,
    scores: np.ndarray,
    top_frac: float,
    direction: str,
) -> tuple[float, float, int, int, int]:
    assert direction in ("asc", "desc")

    mask = np.isfinite(scores) & np.isfinite(is_active.astype(float))
    scores = scores[mask]
    is_active = is_active[mask]

    n = int(scores.shape[0])
    if n == 0:
        return np.nan, np.nan, 0, 0, 0

    a = int(is_active.sum())
    if a == 0:
        return np.nan, np.nan, n, 0, 0

    top_n = max(1, int(np.ceil(top_frac * n)))
    order = np.argsort(scores)
    if direction == "desc":
        order = order[::-1]
    top_idx = order[:top_n]
    h_top = int(is_active[top_idx].sum())

    ef = (h_top / top_n) / (a / n) if a > 0 else np.nan
    if a == 0:
        n_ef = np.nan
    else:
        ef_max = (min(a, top_n) / top_n) / (a / n)
        n_ef = ef / ef_max if ef_max > 0 else np.nan

    return (
        float(ef) if ef == ef else np.nan,
        float(n_ef) if n_ef == n_ef else np.nan,
        n,
        a,
        h_top,
    )


def compute_ef_for_metric(df: pd.DataFrame, metric: str) -> Optional[dict[str, float]]:
    if metric != "pki_oracle" and metric not in df.columns:
        return None

    if "standard_value" not in df.columns and "pValue" not in df.columns:
        return None

    if metric == "pki_oracle":
        if "pValue" in df.columns:
            score_series = pd.to_numeric(df["pValue"], errors="coerce")
        else:
            ki = pd.to_numeric(df["standard_value"], errors="coerce")
            score_series = 9.0 - np.log10(ki.clip(lower=1e-12))
        mask = score_series.notna()
    else:
        mask = df["standard_value"].notna() & df[metric].notna()
        if "pValue" in df.columns:
            mask = mask & df["pValue"].notna()

    dfm = df.loc[mask].copy()
    if dfm.empty:
        return None

    if "pValue" in dfm.columns:
        pki_vals = pd.to_numeric(dfm["pValue"], errors="coerce")
        classes = pd.Series(
            np.where(
                pki_vals >= 7.0,
                "high",
                np.where(pki_vals > 6.0, "medium", "low"),
            ),
            index=dfm.index,
        ).astype(str)
    elif "activity_class" in dfm.columns:
        classes = dfm["activity_class"].astype(str)
    else:
        classes = dfm["standard_value"].apply(classify_activity).astype(str)

    if metric == "pki_oracle":
        scores = score_series.loc[dfm.index].astype(float).to_numpy()
    else:
        scores = dfm[metric].astype(float).to_numpy()

    result: dict[str, float] = {}
    sets = {
        "high_active": (classes == "high"),
        "active": (classes.isin(["high", "medium"])),
        "low": (classes == "low"),
    }

    for frac_name, frac in TOP_FRACS.items():
        for set_name, is_act_mask in sets.items():
            direction = detect_score_direction(metric, set_name)
            is_active = is_act_mask.to_numpy()
            ef, n_ef, n, a, h_top = compute_ef(is_active, scores, frac, direction)
            result[f"EF{frac_name}_{set_name}"] = ef
            result[f"nEF{frac_name}_{set_name}"] = n_ef
            result[f"N_{frac_name}_{set_name}"] = n
            result[f"A_{frac_name}_{set_name}"] = a
            result[f"H_top{frac_name}_{set_name}"] = h_top
            result[f"direction_{frac_name}_{set_name}"] = direction

    return result


# ---------------------------------------------------------------------------
# Pipeline entry
# ---------------------------------------------------------------------------


def run_enrichment(cfg: AnalysisConfig) -> pd.DataFrame:
    all_rows: list[dict] = []
    for target in cfg.targets:
        df = load_merged(cfg.merged_dir, target)
        target_rows: list[dict] = []
        for method_id in cfg.methods:
            metric = cfg.primary_metric(method_id)
            res = compute_ef_for_metric(df, metric)
            if not res:
                continue
            row = {"metric": metric, **res}
            row.update(
                target=target.lower(),
                protein=target.upper(),
                method_id=method_id,
                method_label=cfg.method_label(method_id),
            )
            target_rows.append(row)
            all_rows.append(row)

        if target_rows:
            pd.DataFrame(target_rows).to_csv(
                cfg.tables_dir / "enrichment" / "per_target" / f"ef_summary_{target}.csv",
                index=False,
            )

    summary = pd.DataFrame(all_rows)
    out = cfg.tables_dir / "enrichment" / "ef_summary_all_proteins.csv"
    summary.to_csv(out, index=False)
    return summary
