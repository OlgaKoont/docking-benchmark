"""Results aggregation - combine metrics from all methods per protein."""

import pandas as pd
from pathlib import Path
from typing import Dict, List, Optional


def aggregate_protein_results(
    output_dir: Path,
    protein_name: str,
    methods: List[str],
) -> pd.DataFrame:
    """
    Aggregate results from all methods for a single protein.
    
    Creates a combined table: ligand_id | method1_affinity | method2_affinity | ...
    
    Args:
        output_dir: Directory with results organized by protein.
        protein_name: Name of the protein.
        methods: List of method names to aggregate.
    
    Returns:
        DataFrame with combined metrics.
    """
    protein_dir = output_dir / protein_name
    metrics_dir = protein_dir / "metrics"
    metrics_dir.mkdir(parents=True, exist_ok=True)
    
    # Collect metrics from each method
    all_metrics = {}
    
    for method in methods:
        method_metrics_file = metrics_dir / f"{method}_metrics.csv"
        if method_metrics_file.exists():
            try:
                df = pd.read_csv(method_metrics_file)
                all_metrics[method] = df
            except Exception as e:
                print(f"    Warning: Could not load metrics for {method}: {e}")
                continue
    
    if not all_metrics:
        print(f"    Warning: No metrics found for protein {protein_name}")
        return pd.DataFrame()
    
    # Combine metrics by ligand
    combined_data = {}
    
    for method, df in all_metrics.items():
        for _, row in df.iterrows():
            ligand = row.get('ligand', row.get('ligand_id', 'unknown'))
            if ligand not in combined_data:
                combined_data[ligand] = {'ligand': ligand}
            
            # Add method-specific metrics
            for col in df.columns:
                if col not in ['method', 'protein', 'ligand', 'ligand_id', 'output_file']:
                    col_name = f"{method}_{col}"
                    value = row[col]
                    # Handle NaN values
                    if pd.notna(value):
                        combined_data[ligand][col_name] = value
    
    # Create combined DataFrame
    combined_df = pd.DataFrame(list(combined_data.values()))
    
    # Save combined metrics
    combined_file = metrics_dir / "combined_metrics.csv"
    combined_df.to_csv(combined_file, index=False)
    print(f"    Saved combined metrics for {protein_name}: {len(combined_df)} ligands, {len(methods)} methods")
    
    return combined_df


def aggregate_all_proteins(
    output_dir: Path,
    methods: List[str],
) -> Dict[str, pd.DataFrame]:
    """
    Aggregate results for all proteins.
    
    Args:
        output_dir: Directory with results organized by protein.
        methods: List of method names.
    
    Returns:
        Dictionary mapping protein names to combined DataFrames.
    """
    results = {}
    
    # Find all protein directories
    for protein_dir in output_dir.iterdir():
        if not protein_dir.is_dir() or protein_dir.name == 'global':
            continue
        
        protein_name = protein_dir.name
        combined_df = aggregate_protein_results(output_dir, protein_name, methods)
        if not combined_df.empty:
            results[protein_name] = combined_df
    
    print(f"  Aggregated results for {len(results)} proteins")
    return results

