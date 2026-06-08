"""Docking box preparation - simplified function-based approach."""

import json
from pathlib import Path
from typing import Dict, List, Optional, Any, Tuple


def _normalize_protein_id(protein_id: str) -> str:
    """Normalize protein ID to lowercase for case-insensitive comparison."""
    return protein_id.lower()


def _get_setting(protein_id: str, settings: Dict[str, Any], key: str, default: Any = None) -> Any:
    """Get setting for protein with case-insensitive lookup."""
    normalized_id = _normalize_protein_id(protein_id)
    if protein_id in settings:
        return settings[protein_id].get(key, default)
    for setting_key, setting_value in settings.items():
        if _normalize_protein_id(setting_key) == normalized_id:
            return setting_value.get(key, default)
    return default


def _extract_coordinates(pdb_path: Path) -> List[Tuple[float, float, float]]:
    """Extract coordinates from PDB file."""
    coords = []
    with pdb_path.open("r", encoding="utf-8") as handle:
        for line in handle:
            if line.startswith(("ATOM", "HETATM")):
                try:
                    x = float(line[30:38].strip())
                    y = float(line[38:46].strip())
                    z = float(line[46:54].strip())
                    coords.append((x, y, z))
                except (ValueError, IndexError):
                    continue
    return coords


def _locate_ligand_pdb(protein_path: Path, processed_dir: Path) -> Optional[Path]:
    """Try to locate cleaned ligand PDB file."""
    protein_id = protein_path.stem
    normalized_id = _normalize_protein_id(protein_id)
    
    # Try cleaned ligand PDB first (extracted from complex)
    ligand_pdb = processed_dir / "proteins" / "cleaned" / f"{normalized_id}_ligand.pdb"
    if ligand_pdb.exists():
        return ligand_pdb
    
    return None


def _calculate_labox(
    protein_path: Path,
    processed_dir: Path,
    settings: Dict[str, Any],
    labox_config: Dict[str, Any],
) -> Dict[str, List[float]]:
    """
    Calculate docking box using LaBOX algorithm.
    Uses ligand coordinates from complex if available, otherwise falls back to protein.
    """
    protein_id = protein_path.stem
    normalized_id = _normalize_protein_id(protein_id)
    
    # Get settings
    overrides = {}
    if protein_id in settings:
        overrides = settings[protein_id]
    else:
        for setting_key, setting_value in settings.items():
            if _normalize_protein_id(setting_key) == normalized_id:
                overrides = setting_value
                break
    
    scale = overrides.get("scale", labox_config.get("scale", 2.0))
    min_size = overrides.get("min_size", labox_config.get("min_size", 4.0))
    
    # Try to use cleaned ligand PDB first (extracted from complex)
    source_path = _locate_ligand_pdb(protein_path, processed_dir)
    use_ligand = source_path is not None
    
    # Fallback to protein if ligand not found
    if source_path is None or not source_path.exists():
        # Try cleaned protein PDB
        cleaned_pdb = processed_dir / "proteins" / "cleaned" / f"{normalized_id}_chainA.pdb"
        if cleaned_pdb.exists():
            source_path = cleaned_pdb
        else:
            # Try processed PDB
            processed_pdb = processed_dir / "proteins" / f"{normalized_id}.pdb"
            if processed_pdb.exists():
                source_path = processed_pdb
            else:
                # Fallback to original
                source_path = protein_path
        use_ligand = False
    
    coords = _extract_coordinates(source_path)
    if not coords:
        raise RuntimeError(f"No coordinates found for box calculation in {source_path}")
    
    x_vals, y_vals, z_vals = zip(*coords)
    min_x, max_x = min(x_vals), max(x_vals)
    min_y, max_y = min(y_vals), max(y_vals)
    min_z, max_z = min(z_vals), max(z_vals)
    
    center = [
        round((min_x + max_x) / 2.0, 3),
        round((min_y + max_y) / 2.0, 3),
        round((min_z + max_z) / 2.0, 3)
    ]
    
    # Calculate size for each dimension
    size_x = round(abs(max_x - min_x) * scale, 3)
    size_y = round(abs(max_y - min_y) * scale, 3)
    size_z = round(abs(max_z - min_z) * scale, 3)
    
    # Ensure minimal size of 20 Å for each dimension separately
    # If size < 20, set to 20; if size >= 20, keep calculated size
    min_dimension_size = 20.0
    if size_x < min_dimension_size:
        size_x = min_dimension_size
    if size_y < min_dimension_size:
        size_y = min_dimension_size
    if size_z < min_dimension_size:
        size_z = min_dimension_size
    
    size = [size_x, size_y, size_z]
    
    source_type = "ligand" if use_ligand else "protein"
    print(f"    Using {source_type} coordinates from: {source_path.name}")
    
    return {"center": center, "size": size}


def _default_box() -> Dict[str, List[float]]:
    """Return default box."""
    return {"center": [0.0, 0.0, 0.0], "size": [20.0, 20.0, 20.0]}


def prepare_boxes(
    protein_dir: Path,
    processed_dir: Path,
    settings: Optional[Dict[str, Any]] = None,
    labox_config: Optional[Dict[str, Any]] = None,
    overwrite: bool = False,
) -> Dict[str, Dict[str, List[float]]]:
    """
    Prepare docking boxes for all proteins.
    
    Args:
        protein_dir: Directory with input PDB files.
        processed_dir: Directory with processed proteins and to save boxes.
        settings: Optional per-protein box settings.
        labox_config: Optional LaBOX configuration.
        overwrite: Whether to overwrite existing boxes.
    
    Returns:
        Dictionary mapping protein names to box info.
    """
    settings = settings or {}
    labox_config = labox_config or {}
    box_dir = processed_dir / "boxes"
    box_dir.mkdir(parents=True, exist_ok=True)
    
    results = {}
    
    # Find all PDB files
    pdb_files = list(protein_dir.glob("*.pdb"))
    if not pdb_files:
        print(f"  Warning: No PDB files found in {protein_dir}")
        return results
    
    print(f"  Preparing boxes for {len(pdb_files)} proteins...")
    
    for pdb_path in pdb_files:
        protein_id = pdb_path.stem
        normalized_id = _normalize_protein_id(protein_id)
        
        # Check if box already exists
        box_file = box_dir / f"{normalized_id}.json"
        if box_file.exists() and not overwrite:
            with open(box_file, 'r') as f:
                results[protein_id] = json.load(f)
            continue
        
        try:
            # Calculate box using LaBOX
            box_info = _calculate_labox(pdb_path, processed_dir, settings, labox_config)
            
            # Save box
            with open(box_file, 'w') as f:
                json.dump(box_info, f, indent=2)
            
            results[protein_id] = box_info
            print(f"    Prepared box for {protein_id}: center={box_info['center']}, size={box_info['size']}")
            
        except Exception as e:
            print(f"    Error preparing box for {protein_id}: {e}")
            # Use default box as fallback
            box_info = _default_box()
            results[protein_id] = box_info
            with open(box_file, 'w') as f:
                json.dump(box_info, f, indent=2)
            continue
    
    print(f"  Successfully prepared boxes for {len(results)}/{len(pdb_files)} proteins")
    return results

