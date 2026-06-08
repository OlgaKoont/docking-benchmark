"""Main pipeline - simplified function-based approach."""

from pathlib import Path
from typing import Dict, List, Optional

from .config import load_config, load_methods_config
from .docking import METHODS
from .aggregation import aggregate_all_proteins
# Lazy import analysis to avoid dependency issues in method-specific environments
# from .analysis import calculate_protein_statistics, compare_methods, generate_protein_plots
from .utils.settings import load_interaction_config, load_protein_settings, load_box_settings


def run_pipeline(
    config: Optional[Dict] = None,
    methods_config: Optional[Dict] = None,
    stage: str = 'all',
) -> None:
    """
    Run the complete docking benchmark pipeline.
    
    Pipeline stages:
    1. PREPARATION: Prepare proteins, ligands, and boxes (once for all methods)
    2. DOCKING: Run docking for each method
    3. AGGREGATION: Combine results from all methods per protein
    4. ANALYSIS: Calculate statistics and generate plots per protein
    
    Args:
        config: Main configuration dictionary. If None, loads from default.
        methods_config: Methods configuration. If None, loads from default.
        stage: Which stage to run: 'preparation', 'docking', 'aggregation', 'analysis', or 'all'.
    """
    # Load configs
    if config is None:
        config = load_config()
    if methods_config is None:
        methods_config = load_methods_config()
    
    # Setup directories
    base_dir = Path(config['base_dir'])
    protein_dir = Path(config['protein_dir'])
    ligand_dir = Path(config['ligand_dir'])
    processed_dir = base_dir / config.get('processed_dir', 'processed')
    output_dir = base_dir / config.get('output_dir', 'results')
    
    # Create directories
    processed_dir.mkdir(parents=True, exist_ok=True)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Get methods to run
    methods = config.get('methods', ['gnina'])
    
    # Load settings
    protein_settings = load_protein_settings(config.get('protein_settings_file'))
    box_settings = load_box_settings(config.get('box_settings_file'))
    labox_config = config.get('labox', {})
    interaction_config = load_interaction_config(config.get('interaction_config_file'))
    
    print("=" * 80)
    print("DOCKING BENCHMARK PIPELINE 2.0")
    print("=" * 80)
    print(f"Methods: {', '.join(methods)}")
    print(f"Protein directory: {protein_dir}")
    print(f"Ligand directory: {ligand_dir}")
    print(f"Processed directory: {processed_dir}")
    print(f"Output directory: {output_dir}")
    print("=" * 80)
    
    # Parse stage - can be comma-separated list or 'all'
    if isinstance(stage, str):
        stage_parts = [s.strip() for s in stage.split(',')]
    else:
        stage_parts = [stage]
    
    run_preparation = 'preparation' in stage_parts or 'all' in stage_parts
    run_docking = 'docking' in stage_parts or 'all' in stage_parts
    run_aggregation = 'aggregation' in stage_parts or 'all' in stage_parts
    run_analysis = 'analysis' in stage_parts or 'all' in stage_parts
    
    # STAGE 1: PREPARATION
    if run_preparation:
        print("\n" + "=" * 80)
        print("STAGE 1: PREPARATION")
        print("=" * 80)
        
        # Lazy import preprocessing functions to avoid loading meeko/gemmi when not needed
        from .preprocessing.proteins import prepare_proteins
        from .preprocessing.ligands import prepare_ligands
        from .preprocessing.boxes import prepare_boxes
        
        # Prepare proteins (once for all methods)
        print("\n[1.1] Preparing proteins...")
        prepare_proteins(
            protein_dir=protein_dir,
            processed_dir=processed_dir,
            settings=protein_settings,
            interaction_config=interaction_config,
        )
        
        # Prepare ligands (once for all methods)
        print("\n[1.2] Preparing ligands...")
        prepare_ligands(
            ligand_dir=ligand_dir,
            processed_dir=processed_dir,
            interaction_config=interaction_config,
        )
        
        # Prepare boxes (once for all methods)
        print("\n[1.3] Preparing boxes...")
        prepare_boxes(
            protein_dir=protein_dir,
            processed_dir=processed_dir,
            settings=box_settings,
            labox_config=labox_config,
        )
        
        # Prepare PLAPT sequences if PLAPT is in methods
        if 'plapt' in methods:
            print("\n[1.4] Preparing PLAPT sequences...")
            try:
                from .docking.plapt import preprocess_plapt
                preprocess_plapt(
                    protein_dir=protein_dir,
                    processed_dir=processed_dir,
                    interaction_config=interaction_config,
                )
            except Exception as e:
                print(f"  Warning: PLAPT preprocessing failed: {e}")
        
        # Prepare DynamicBind if DynamicBind is in methods
        if 'dynamicbind' in methods:
            print("\n[1.5] Preparing DynamicBind...")
            try:
                from .docking.dynamicbind import preprocess_dynamicbind
                preprocess_dynamicbind(
                    protein_dir=protein_dir,
                    processed_dir=processed_dir,
                    interaction_config=interaction_config,
                )
            except Exception as e:
                print(f"  Warning: DynamicBind preprocessing failed: {e}")
        
        # Prepare Interformer pockets/ligands
        if 'interformer' in methods:
            print("\n[1.6] Preparing Interformer pockets and ligands...")
            try:
                from .docking.interformer import preprocess_interformer
                interformer_config = methods_config.get('interformer', {})
                preprocess_interformer(
                    protein_dir=protein_dir,
                    processed_dir=processed_dir,
                    interaction_config=interaction_config,
                    config=interformer_config,
                )
            except Exception as e:
                print(f"  Warning: Interformer preprocessing failed: {e}")
        
        # NOTE: Method integration removed by request.
        
    # STAGE 2: DOCKING
    if run_docking:
        print("\n" + "=" * 80)
        print("STAGE 2: DOCKING")
        print("=" * 80)
        
        for method in methods:
            if method not in METHODS:
                print(f"\n  Warning: Method '{method}' not found, skipping")
                continue
            
            print(f"\n[2.{methods.index(method)+1}] Running {method.upper()} docking...")
            dock_func = METHODS[method]['dock']
            method_config = methods_config.get(method, {})
            method_config['labox'] = labox_config
            
            try:
                # PLAPT and DynamicBind require ligand_dir parameter
                if method == 'plapt':
                    dock_func(
                        processed_dir=processed_dir,
                        output_dir=output_dir,
                        config=method_config,
                        interaction_config=interaction_config,
                        ligand_dir=ligand_dir,
                    )
                elif method == 'dynamicbind':
                    dock_func(
                        processed_dir=processed_dir,
                        output_dir=output_dir,
                        config=method_config,
                        interaction_config=interaction_config,
                        ligand_dir=ligand_dir,
                        protein_dir=protein_dir,
                    )
                elif method == 'interformer':
                    dock_func(
                        processed_dir=processed_dir,
                        output_dir=output_dir,
                        config=method_config,
                        interaction_config=interaction_config,
                        ligand_dir=ligand_dir,
                        protein_dir=protein_dir,
                    )
                else:
                    dock_func(
                        processed_dir=processed_dir,
                        output_dir=output_dir,
                        config=method_config,
                        interaction_config=interaction_config,
                    )
            except Exception as e:
                print(f"  Error in {method} docking: {e}")
                import traceback
                traceback.print_exc()
                continue
    
    # STAGE 3: EXTRACT METRICS
    if run_docking or run_aggregation or run_analysis:
        print("\n" + "=" * 80)
        print("STAGE 3: EXTRACTING METRICS")
        print("=" * 80)
        
        for method in methods:
            if method not in METHODS:
                continue
            
            print(f"\n[3.{methods.index(method)+1}] Extracting {method.upper()} metrics...")
            extract_func = METHODS[method]['extract_metrics']
            
            try:
                metrics = extract_func(output_dir)
                
                # Save metrics per protein
                import pandas as pd
                if metrics:
                    df = pd.DataFrame(metrics)
                    for protein_name in df['protein'].unique():
                        protein_metrics = df[df['protein'] == protein_name]
                        metrics_dir = output_dir / protein_name / "metrics"
                        metrics_dir.mkdir(parents=True, exist_ok=True)
                        metrics_file = metrics_dir / f"{method}_metrics.csv"
                        protein_metrics.to_csv(metrics_file, index=False)
                        print(f"    Saved metrics for {protein_name}: {len(protein_metrics)} ligands")
            except Exception as e:
                print(f"  Error extracting {method} metrics: {e}")
                import traceback
                traceback.print_exc()
                continue
    
    # STAGE 4: AGGREGATION
    if run_aggregation or run_analysis:
        print("\n" + "=" * 80)
        print("STAGE 4: AGGREGATING RESULTS")
        print("=" * 80)
        
        # Extract metrics if not already done (for aggregation/analysis stages)
        # Metrics are already extracted in stage 3 for 'all' stage
        if 'aggregation' in stage_parts or 'analysis' in stage_parts:
            for method in methods:
                if method not in METHODS:
                    continue
                print(f"  Extracting {method.upper()} metrics...")
                extract_func = METHODS[method]['extract_metrics']
                try:
                    metrics = extract_func(output_dir)
                    import pandas as pd
                    if metrics:
                        df = pd.DataFrame(metrics)
                        for protein_name in df['protein'].unique():
                            protein_metrics = df[df['protein'] == protein_name]
                            metrics_dir = output_dir / protein_name / "metrics"
                            metrics_dir.mkdir(parents=True, exist_ok=True)
                            metrics_file = metrics_dir / f"{method}_metrics.csv"
                            protein_metrics.to_csv(metrics_file, index=False)
                except Exception as e:
                    print(f"  Warning: Could not extract {method} metrics: {e}")
        
        # Aggregate results from all methods
        try:
            aggregate_all_proteins(output_dir, methods)
        except Exception as e:
            print(f"  Error aggregating results: {e}")
            import traceback
            traceback.print_exc()
    
    # STAGE 5: ANALYSIS
    if run_analysis:
        print("\n" + "=" * 80)
        print("STAGE 5: ANALYSIS (STATISTICS & PLOTS)")
        print("=" * 80)
        
        # Lazy import analysis to avoid dependency issues
        try:
            from .analysis import calculate_protein_statistics, compare_methods, generate_protein_plots
        except ImportError as e:
            print(f"  Warning: Analysis dependencies not available: {e}")
            print("  Skipping analysis stage. Install seaborn, scipy, and matplotlib for analysis.")
            return
        
        # Process each protein
        for protein_dir in output_dir.iterdir():
            if not protein_dir.is_dir() or protein_dir.name == 'global':
                continue
            
            protein_name = protein_dir.name
            print(f"\n[5.{list(output_dir.iterdir()).index(protein_dir)+1}] Analyzing {protein_name}...")
            
            # Load combined metrics
            combined_file = protein_dir / "metrics" / "combined_metrics.csv"
            if not combined_file.exists():
                print(f"    Warning: No combined metrics for {protein_name}, skipping")
                continue
            
            import pandas as pd
            combined_df = pd.read_csv(combined_file)
            
            if combined_df.empty:
                print(f"    Warning: Empty metrics for {protein_name}, skipping")
                continue
            
            # Calculate statistics
            try:
                calculate_protein_statistics(output_dir, protein_name, combined_df)
                compare_methods(output_dir, protein_name, combined_df)
            except Exception as e:
                print(f"    Error calculating statistics for {protein_name}: {e}")
            
            # Generate plots
            try:
                generate_protein_plots(output_dir, protein_name, combined_df)
            except Exception as e:
                print(f"    Error generating plots for {protein_name}: {e}")
    
    print("\n" + "=" * 80)
    print("PIPELINE COMPLETE!")
    print("=" * 80)

