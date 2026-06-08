"""Statistical analysis - per protein."""

import pandas as pd
import numpy as np
import json
from pathlib import Path
from typing import Dict, List, Optional


def calculate_protein_statistics(
    output_dir: Path,
    protein_name: str,
    combined_df: Optional[pd.DataFrame] = None,
) -> Dict:
    """
    Calculate statistics for a protein's docking results.
    
    Args:
        output_dir: Directory with results.
        protein_name: Name of the protein.
        combined_df: Optional pre-loaded combined DataFrame.
    
    Returns:
        Dictionary with statistics.
    """
    if combined_df is None:
        metrics_file = output_dir / protein_name / "metrics" / "combined_metrics.csv"
        if not metrics_file.exists():
            return {}
        combined_df = pd.read_csv(metrics_file)
    
    if combined_df.empty:
        return {}
    
    stats_dir = output_dir / protein_name / "statistics"
    stats_dir.mkdir(parents=True, exist_ok=True)
    
    # Calculate statistics for numeric columns
    numeric_cols = combined_df.select_dtypes(include=[np.number]).columns
    stats = {}
    
    for col in numeric_cols:
        stats[col] = {
            'mean': float(combined_df[col].mean()) if not combined_df[col].isna().all() else None,
            'std': float(combined_df[col].std()) if not combined_df[col].isna().all() else None,
            'min': float(combined_df[col].min()) if not combined_df[col].isna().all() else None,
            'max': float(combined_df[col].max()) if not combined_df[col].isna().all() else None,
            'median': float(combined_df[col].median()) if not combined_df[col].isna().all() else None,
            'count': int(combined_df[col].notna().sum()),
        }
    
    # Save statistics
    stats_file = stats_dir / "summary_stats.json"
    with open(stats_file, 'w') as f:
        json.dump(stats, f, indent=2)
    
    print(f"    Calculated statistics for {protein_name}")
    return stats


def compare_methods(
    output_dir: Path,
    protein_name: str,
    combined_df: Optional[pd.DataFrame] = None,
    metric: str = 'affinity',
) -> pd.DataFrame:
    """
    Compare methods for a protein.
    
    Args:
        output_dir: Directory with results.
        protein_name: Name of the protein.
        combined_df: Optional pre-loaded combined DataFrame.
        metric: Metric to compare (e.g., 'affinity').
    
    Returns:
        DataFrame with method comparison.
    """
    if combined_df is None:
        metrics_file = output_dir / protein_name / "metrics" / "combined_metrics.csv"
        if not metrics_file.exists():
            return pd.DataFrame()
        combined_df = pd.read_csv(metrics_file)
    
    if combined_df.empty:
        return pd.DataFrame()
    
    stats_dir = output_dir / protein_name / "statistics"
    stats_dir.mkdir(parents=True, exist_ok=True)
    
    # Find columns with the metric (e.g., gnina_affinity, qvina_affinity)
    metric_cols = [col for col in combined_df.columns if col.endswith(f'_{metric}')]
    
    if not metric_cols:
        print(f"    Warning: No {metric} columns found for {protein_name}")
        return pd.DataFrame()
    
    # Calculate statistics per method
    comparison_data = []
    for col in metric_cols:
        method = col.replace(f'_{metric}', '')
        values = combined_df[col].dropna()
        
        if len(values) > 0:
            comparison_data.append({
                'method': method,
                'metric': metric,
                'mean': float(values.mean()),
                'std': float(values.std()),
                'min': float(values.min()),
                'max': float(values.max()),
                'median': float(values.median()),
                'count': int(len(values)),
            })
    
    comparison_df = pd.DataFrame(comparison_data)
    
    # Save comparison
    comparison_file = stats_dir / "method_comparison.csv"
    comparison_df.to_csv(comparison_file, index=False)
    
    print(f"    Compared methods for {protein_name}: {len(comparison_df)} methods")
    return comparison_df

