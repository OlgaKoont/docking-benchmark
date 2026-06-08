"""Plot generation - per protein."""

import pandas as pd
import matplotlib
matplotlib.use('Agg')  # Non-interactive backend
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
from typing import Optional, List

# Set style
sns.set_style("whitegrid")
plt.rcParams['figure.figsize'] = (10, 6)


def generate_protein_plots(
    output_dir: Path,
    protein_name: str,
    combined_df: Optional[pd.DataFrame] = None,
) -> List[Path]:
    """
    Generate plots for a protein's docking results.
    
    Args:
        output_dir: Directory with results.
        protein_name: Name of the protein.
        combined_df: Optional pre-loaded combined DataFrame.
    
    Returns:
        List of generated plot file paths.
    """
    if combined_df is None:
        metrics_file = output_dir / protein_name / "metrics" / "combined_metrics.csv"
        if not metrics_file.exists():
            return []
        combined_df = pd.read_csv(metrics_file)
    
    if combined_df.empty:
        return []
    
    plots_dir = output_dir / protein_name / "plots"
    plots_dir.mkdir(parents=True, exist_ok=True)
    
    generated_plots = []
    
    # Find affinity columns
    affinity_cols = [col for col in combined_df.columns if 'affinity' in col.lower()]
    
    if not affinity_cols:
        print(f"    Warning: No affinity columns found for {protein_name}")
        return []
    
    # 1. Affinity comparison boxplot
    try:
        plot_data = []
        for col in affinity_cols:
            method = col.replace('_affinity', '').replace('affinity', '')
            for val in combined_df[col].dropna():
                plot_data.append({'method': method, 'affinity': val})
        
        if plot_data:
            plot_df = pd.DataFrame(plot_data)
            plt.figure(figsize=(10, 6))
            sns.boxplot(data=plot_df, x='method', y='affinity')
            plt.title(f'Affinity Comparison - {protein_name}')
            plt.xlabel('Method')
            plt.ylabel('Affinity (kcal/mol)')
            plt.xticks(rotation=45)
            plt.tight_layout()
            
            plot_file = plots_dir / "affinity_comparison.png"
            plt.savefig(plot_file, dpi=300, bbox_inches='tight')
            plt.close()
            generated_plots.append(plot_file)
    except Exception as e:
        print(f"    Warning: Could not generate affinity comparison plot: {e}")
    
    # 2. Method correlation heatmap
    try:
        if len(affinity_cols) > 1:
            corr_data = combined_df[affinity_cols].corr()
            plt.figure(figsize=(8, 6))
            sns.heatmap(corr_data, annot=True, fmt='.2f', cmap='coolwarm', center=0)
            plt.title(f'Method Correlation - {protein_name}')
            plt.tight_layout()
            
            plot_file = plots_dir / "method_correlation.png"
            plt.savefig(plot_file, dpi=300, bbox_inches='tight')
            plt.close()
            generated_plots.append(plot_file)
    except Exception as e:
        print(f"    Warning: Could not generate correlation plot: {e}")
    
    # 3. Scatter plot comparing methods (if 2 methods)
    try:
        if len(affinity_cols) == 2:
            col1, col2 = affinity_cols
            method1 = col1.replace('_affinity', '').replace('affinity', '')
            method2 = col2.replace('_affinity', '').replace('affinity', '')
            
            plt.figure(figsize=(8, 8))
            plt.scatter(combined_df[col1], combined_df[col2], alpha=0.5)
            plt.xlabel(f'{method1} Affinity')
            plt.ylabel(f'{method2} Affinity')
            plt.title(f'Affinity Comparison: {method1} vs {method2} - {protein_name}')
            
            # Add diagonal line
            min_val = min(combined_df[col1].min(), combined_df[col2].min())
            max_val = max(combined_df[col1].max(), combined_df[col2].max())
            plt.plot([min_val, max_val], [min_val, max_val], 'r--', alpha=0.5)
            
            plt.tight_layout()
            plot_file = plots_dir / f"affinity_scatter_{method1}_vs_{method2}.png"
            plt.savefig(plot_file, dpi=300, bbox_inches='tight')
            plt.close()
            generated_plots.append(plot_file)
    except Exception as e:
        print(f"    Warning: Could not generate scatter plot: {e}")
    
    if generated_plots:
        print(f"    Generated {len(generated_plots)} plots for {protein_name}")
    
    return generated_plots

