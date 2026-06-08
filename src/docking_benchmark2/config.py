"""Configuration loading and management."""

import yaml
from pathlib import Path
from typing import Dict, Any, Optional


def load_config(config_path: Optional[Path] = None) -> Dict[str, Any]:
    """
    Load configuration from YAML file.
    
    Args:
        config_path: Path to config file. If None, uses default.
    
    Returns:
        Dictionary with configuration.
    """
    if config_path is None:
        config_path = Path(__file__).parent.parent.parent / "config" / "default_config.yaml"
    
    if not config_path.exists():
        # Return default config if file doesn't exist
        return {
            'base_dir': 'data',
            'protein_dir': 'data/input/proteins',
            'ligand_dir': 'data/input/ligands',
            'processed_dir': 'data/processed',
            'output_dir': 'data/results',
            'random_state': 42,
            'methods': ['gnina'],
        }
    
    with open(config_path, 'r') as f:
        config = yaml.safe_load(f)
    
    return config


def load_methods_config(config_path: Optional[Path] = None) -> Dict[str, Any]:
    """
    Load methods configuration.
    
    Args:
        config_path: Path to methods config file. If None, uses default.
    
    Returns:
        Dictionary with methods configuration.
    """
    if config_path is None:
        config_path = Path(__file__).parent.parent.parent / "config" / "methods_config.yaml"
    
    if not config_path.exists():
        # Return default methods config
        return {
            'gnina': {
                'binary': 'gnina',
                'exhaustiveness': 8,
                'num_modes': 9,
                'use_cnn': False,
            },
            'qvina': {
                'binary': 'qvina02',
                'exhaustiveness': 8,
            },
        }
    
    with open(config_path, 'r') as f:
        config = yaml.safe_load(f)
    
    return config

