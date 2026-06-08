"""PoseBusters pass rates: 100%, 90%, 50% and per-check breakdown."""

from __future__ import annotations

import pandas as pd

from .config import AnalysisConfig
from .constants import POSEBUSTERS_EXCLUDE
from .data import resolve_posebusters_csv


def _boolean_check_columns(df: pd.DataFrame) -> list[str]:
    skip = {"file", "molecule", "position"}
    cols: list[str] = []
    for c in df.columns:
        if c in skip:
            continue
        if pd.api.types.is_bool_dtype(df[c]):
            cols.append(c)
        else:
            s = df[c].dropna().astype(str).str.lower()
            if len(s) > 0 and s.isin(["true", "false"]).all():
                cols.append(c)
    return cols


def _to_bool(series: pd.Series) -> pd.Series:
    if pd.api.types.is_bool_dtype(series):
        return series
    return series.astype(str).str.lower().map({"true": True, "false": False})


def run_posebusters(cfg: AnalysisConfig) -> tuple[pd.DataFrame, pd.DataFrame]:
    summary_rows: list[dict] = []
    check_rows: list[dict] = []

    for target in cfg.targets:
        for method_id in cfg.methods:
            if method_id in POSEBUSTERS_EXCLUDE:
                continue
            csv_path = resolve_posebusters_csv(
                cfg.posebusters_dir,
                cfg.posebusters_boltz2_dir,
                target,
                method_id,
                cfg.dynamicbind_posebusters_label,
            )
            if csv_path is None:
                continue

            df = pd.read_csv(csv_path)
            bool_cols = _boolean_check_columns(df)
            if not bool_cols:
                continue

            bool_df = pd.DataFrame({c: _to_bool(df[c]) for c in bool_cols})
            frac_pass = bool_df.mean(axis=1)

            summary_rows.append(
                {
                    "target": target.lower(),
                    "method_id": method_id,
                    "method_label": cfg.method_label(method_id),
                    "n_poses": len(bool_df),
                    "n_checks": len(bool_cols),
                    "pass_rate_all": float((frac_pass >= 1.0).mean()),
                    "pass_rate_90pct": float((frac_pass >= 0.90).mean()),
                    "pass_rate_50pct": float((frac_pass >= 0.50).mean()),
                    "mean_frac_pass": float(frac_pass.mean()),
                }
            )

            for check in bool_cols:
                check_rows.append(
                    {
                        "target": target.lower(),
                        "method_id": method_id,
                        "method_label": cfg.method_label(method_id),
                        "check": check,
                        "pass_rate": float(_to_bool(df[check]).mean()),
                    }
                )

    summary = pd.DataFrame(summary_rows)
    checks = pd.DataFrame(check_rows)

    summary.to_csv(cfg.tables_dir / "posebusters" / "pass_rates_by_target_method.csv", index=False)
    checks.to_csv(cfg.tables_dir / "posebusters" / "pass_rates_by_check_target_method.csv", index=False)
    return summary, checks
