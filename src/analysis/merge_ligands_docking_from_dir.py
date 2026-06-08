#!/usr/bin/env python3
"""
Скрипт для объединения данных лигандов из ligands_nodubl с метриками докинга,
собранными из директории docking.

Берет CSV файлы из ligands_nodubl, извлекает нужные столбцы и добавляет
метрики докинга (boltz2 и exp_*) из директории data/results_bindingdb/{protein}/docking
или из собранных файлов, объединяя по canonical_smiles.
"""

import argparse
import json
import os
import re
import sys
from pathlib import Path
from typing import Dict, Optional, List, Tuple

import pandas as pd


# Маппинг PDB кодов на имена CSV файлов в ligands_nodubl (без _nodubl.csv)
PDB_TO_CSV_MAPPING = {
    "1ere": "Estrogen_receptor_alpha_Ki_WT_ChEMBL_494",
    "1g5m": "BCL2_Ki_WT_ChEMBL_252",
    "1g5m_wrong": "BCL2L1_BCL-XL_Ki_WT_ChEMBL_287",
    "1tqn": "CYP3A4_Ki_WT_ChEMBL_1",
    "2z5x": "MAO-B_Ki_WT_ChEMBL_246",
    "3eyg": "JAK1_Ki_WT_ChEMBL_2255",
    "3jy9": "JAK2_Ki_WT_ChEMBL_2027",
    "3lxk": "JAK3_Ki_WT_ChEMBL_786",
    "3mjg": "PDGFRB_Ki_WT_ChEMBL_275",
    "4ase": "VEGFR2_Ki_WT_ChEMBL_875",
    "4f65": "FGFR1_Ki_WT_ChEMBL_134",
    "4k7a": "Androgen_receptor_Ki_WT_ChEMBL_516",
    "4tz4": "CRBN_Ki_WT_ChEMBL_127",
    "4wnv": "HTR2B_5HT2B_IC50_WT_ChEMBL_61",
    "4zau": "EGFR_Ki_WT_curated_251",
    "5jkv": "CYP19A1_Aromatase_Ki_WT_ChEMBL_548",
    "5mo4": "ABL1_BCR-ABL_Ki_WT_ChEMBL_693",
    "6gqj": "KIT_Ki_WT_curated_1298",
    "6jok": "PDGFRA_Ki_WT_curated_250",
    "7awe": "PSMB5_Ki_WT_ChEMBL_88",
    "7kk3": "PARP1_Ki_WT_ChEMBL_1075",
    "8zyq": "hERG_Ki_WT_curated_417",
}

# Столбцы для извлечения из ligands_nodubl CSV
LIGANDS_COLUMNS = [
    "canonical_smiles",
    "assay_chembl_id",
    "standard_value",
    "value",
    "pchembl_value",
    "assay_type",
    "document_chembl_id",
    "type",
    "units",
    "uo_units",
]

# Столбцы boltz для добавления
BOLTZ_COLUMNS = [
    "boltz2_affinity_pred_value",
    "boltz2_affinity_probability_binary",
    "boltz2_affinity_pred_value1",
    "boltz2_affinity_probability_binary1",
    "boltz2_affinity_pred_value2",
    "boltz2_affinity_probability_binary2",
]

# Столбцы экспериментальных данных для добавления
EXP_COLUMNS = [
    "exp_value",
    "exp_standard_value",
    "exp_pchembl_value",
]

# Path to Boltz affinity JSON outputs (override via BOLTZ_RESULTS_DIR env or --boltz-results-dir)
BOLTZ_RESULTS_DIR = os.environ.get("BOLTZ_RESULTS_DIR", "")


def extract_chembl_id_from_path(json_path: str) -> Optional[str]:
    """Извлекает CHEMBL ID из пути к JSON файлу."""
    match = re.search(r'/(CHEMBL\d+)/', json_path)
    if match:
        return match.group(1)
    match = re.search(r'affinity_(CHEMBL\d+)\.json', json_path)
    if match:
        return match.group(1)
    return None


def parse_boltz_affinity_json(json_path: str) -> Optional[Dict[str, float]]:
    """Разбор одного Boltz affinity JSON файла."""
    try:
        with open(json_path, "r") as f:
            data = json.load(f)
    except (OSError, json.JSONDecodeError):
        return None

    result: Dict[str, float] = {}
    fields = [
        "affinity_pred_value",
        "affinity_probability_binary",
        "affinity_pred_value1",
        "affinity_probability_binary1",
        "affinity_pred_value2",
        "affinity_probability_binary2",
    ]
    
    for field in fields:
        if field in data:
            try:
                result[f"boltz2_{field}"] = float(data[field])
            except (ValueError, TypeError):
                pass
    
    return result if result else None


def collect_boltz_metrics(
    protein: str, boltz_results_dir: str = BOLTZ_RESULTS_DIR
) -> Dict[str, Dict[str, float]]:
    """
    Собирает метрики Boltz для данного белка.
    Возвращает словарь: {molecule_chembl_id: {boltz2_...: value, ...}, ...}
    
    Поддерживает две структуры директорий:
    1. Старая: {boltz_results_dir}/{protein}/docking/boltz2_*/predictions/{ligand_dataset}/boltz_results_{ligand_dataset}/predictions/{CHEMBL_ID}/affinity_{CHEMBL_ID}.json
    2. Новая: {boltz_results_dir}/{protein}/docking/boltz2_*/boltz_results_{CHEMBL_ID}/predictions/{CHEMBL_ID}/affinity_{CHEMBL_ID}.json
    """
    result: Dict[str, Dict[str, float]] = {}
    
    protein_dir = os.path.join(boltz_results_dir, protein, "docking")
    if not os.path.isdir(protein_dir):
        return result
    
    boltz_docking_dir = None
    for item in os.listdir(protein_dir):
        if item.startswith("boltz2_") and "aff" in item:
            boltz_docking_dir = os.path.join(protein_dir, item)
            break
    
    if not boltz_docking_dir or not os.path.isdir(boltz_docking_dir):
        return result
    
    # НОВАЯ СТРУКТУРА: boltz_results_{CHEMBL_ID}/predictions/{CHEMBL_ID}/affinity_{CHEMBL_ID}.json
    for item in os.listdir(boltz_docking_dir):
        if item.startswith("boltz_results_") and item.startswith("boltz_results_CHEMBL"):
            boltz_results_chembl_dir = os.path.join(boltz_docking_dir, item)
            if not os.path.isdir(boltz_results_chembl_dir):
                continue
            
            # Извлекаем CHEMBL ID из имени директории
            match = re.search(r'boltz_results_(CHEMBL\d+)', item)
            if not match:
                continue
            chembl_id = match.group(1)
            
            # Ищем файл affinity_{CHEMBL_ID}.json
            predictions_dir = os.path.join(boltz_results_chembl_dir, "predictions", chembl_id)
            if os.path.isdir(predictions_dir):
                json_file = os.path.join(predictions_dir, f"affinity_{chembl_id}.json")
                if os.path.isfile(json_file):
                    metrics = parse_boltz_affinity_json(json_file)
                    if metrics:
                        result[chembl_id] = metrics
    
    # СТАРАЯ СТРУКТУРА: predictions/{ligand_dataset}/boltz_results_{ligand_dataset}/predictions/{CHEMBL_ID}/affinity_{CHEMBL_ID}.json
    predictions_dir = os.path.join(boltz_docking_dir, "predictions")
    if os.path.isdir(predictions_dir):
        for ligand_dataset in os.listdir(predictions_dir):
            ligand_dataset_dir = os.path.join(predictions_dir, ligand_dataset)
            if not os.path.isdir(ligand_dataset_dir):
                continue
            
            boltz_results_subdir = os.path.join(
                ligand_dataset_dir, f"boltz_results_{ligand_dataset}"
            )
            if not os.path.isdir(boltz_results_subdir):
                continue
            
            # Ищем JSON файлы с метриками в поддиректории predictions
            predictions_subdir = os.path.join(boltz_results_subdir, "predictions")
            if os.path.isdir(predictions_subdir):
                for chembl_dir in os.listdir(predictions_subdir):
                    chembl_path = os.path.join(predictions_subdir, chembl_dir)
                    if not os.path.isdir(chembl_path):
                        continue
                    
                    # Ищем файл affinity_{CHEMBL_ID}.json
                    json_file = os.path.join(chembl_path, f"affinity_{chembl_dir}.json")
                    if os.path.isfile(json_file):
                        chembl_id = extract_chembl_id_from_path(json_file) or chembl_dir
                        # Пропускаем, если уже добавили из новой структуры
                        if chembl_id in result:
                            continue
                        metrics = parse_boltz_affinity_json(json_file)
                        if metrics:
                            result[chembl_id] = metrics
    
    return result


def group_by_smiles(df: pd.DataFrame, smiles_col: str = "canonical_smiles") -> pd.DataFrame:
    """
    Группирует DataFrame по SMILES, усредняя числовые метрики.
    """
    if smiles_col not in df.columns:
        return df
    
    # Фильтруем строки с валидными SMILES
    mask = df[smiles_col].notna() & (
        df[smiles_col].astype(str).str.strip() != ''
    ) & (
        df[smiles_col].astype(str).str.lower() != 'nan'
    )
    df_filtered = df[mask].copy()
    
    if len(df_filtered) == 0:
        return df
    
    # Определяем колонки для группировки
    group_cols = [smiles_col]
    if "molecule_chembl_id" in df_filtered.columns:
        group_cols.append("molecule_chembl_id")
    
    # Определяем колонки с метриками (исключаем служебные)
    exclude_cols = set(group_cols)
    metric_cols = [c for c in df_filtered.columns if c not in exclude_cols]
    
    # Группируем по SMILES и агрегируем метрики
    agg_dict = {}
    for col in metric_cols:
        if pd.api.types.is_numeric_dtype(df_filtered[col]):
            agg_dict[col] = 'mean'
        else:
            agg_dict[col] = 'first'
    
    df_grouped = df_filtered.groupby(group_cols, as_index=False).agg(agg_dict)
    
    return df_grouped


# Импортируем функции парсинга из extract_raw_docking_metrics.py
# (копируем их сюда, чтобы не зависеть от другого скрипта)

def parse_qvina_log(log_path: str) -> Optional[Dict[str, float]]:
    """Разбор одного qvina .log файла."""
    try:
        with open(log_path, "r") as f:
            lines = f.readlines()
    except OSError:
        return None

    data_rows: List[Tuple[int, float, float, float]] = []
    in_table = False

    for line in lines:
        line = line.strip()
        if "mode |" in line and "affinity" in line:
            in_table = True
            continue
        if not in_table:
            continue
        if line.startswith("-----"):
            continue
        if re.match(r"^\s*\d+\s+[-\d\.]+\s+[-\d\.]+\s+[-\d\.]+$", line):
            parts = re.split(r"\s+", line)
            if len(parts) >= 4:
                try:
                    mode = int(parts[0])
                    affinity = float(parts[1])
                    rmsd_lb = float(parts[2])
                    rmsd_ub = float(parts[3])
                    data_rows.append((mode, affinity, rmsd_lb, rmsd_ub))
                except (ValueError, IndexError):
                    continue

    if not data_rows:
        return None

    min_aff_row = min(data_rows, key=lambda x: x[1])
    max_aff_row = max(data_rows, key=lambda x: x[1])
    bestpose_row = data_rows[0] if data_rows else None

    if bestpose_row is None:
        return None

    _, aff_min, rmsd_lb_min, rmsd_ub_min = min_aff_row
    _, aff_max, _, _ = max_aff_row
    _, aff_bp, rmsd_lb_bp, rmsd_ub_bp = bestpose_row

    lb_candidates = [row for row in data_rows if row[2] > 0]
    lb_aff = None
    if lb_candidates:
        lb_best = min(lb_candidates, key=lambda x: x[2])
        lb_aff = lb_best[1]

    ub_candidates = [row for row in data_rows if row[3] > 0]
    ub_aff = None
    if ub_candidates:
        ub_best = min(ub_candidates, key=lambda x: x[3])
        ub_aff = ub_best[1]

    return {
        "affinity_min": aff_min,
        "rmsd_lb_min": rmsd_lb_min,
        "rmsd_ub_min": rmsd_ub_min,
        "affinity_bestpose": aff_bp,
        "rmsd_lb_bestpose": rmsd_lb_bp,
        "rmsd_ub_bestpose": rmsd_ub_bp,
        "max_affinity": aff_max,
        "lb_affinity": lb_aff,
        "ub_affinity": ub_aff,
    }


def parse_gnina_log(log_path: str) -> Optional[Dict[str, float]]:
    """Разбор одного gnina .log файла."""
    try:
        with open(log_path, "r") as f:
            lines = f.readlines()
    except OSError:
        return None

    data_rows: List[Tuple[int, float, float, float, float]] = []

    for line in lines:
        line = line.strip()
        if re.match(r"^\d+\s+[-\d\.]+\s+[-\d\.]+\s+[-\d\.]+\s+[-\d\.]+$", line):
            parts = re.split(r"\s+", line)
            if len(parts) != 5:
                continue
            try:
                mode = int(parts[0])
                affinity = float(parts[1])
                intramol = float(parts[2])
                cnn_pose = float(parts[3])
                cnn_affinity = float(parts[4])
            except ValueError:
                continue
            data_rows.append((mode, affinity, intramol, cnn_pose, cnn_affinity))

    if not data_rows:
        return None

    min_aff_row = min(data_rows, key=lambda x: x[1])
    bestpose_row = max(data_rows, key=lambda x: x[3])

    _, aff_min, intr_min, cnn_pose_min, cnn_aff_min = min_aff_row
    _, aff_bp, intr_bp, cnn_pose_bp, cnn_aff_bp = bestpose_row

    return {
        "affinity_min": aff_min,
        "intramol_min": intr_min,
        "cnn_pose_score_min": cnn_pose_min,
        "cnn_affinity_min": cnn_aff_min,
        "affinity_bestpose": aff_bp,
        "intramol_bestpose": intr_bp,
        "cnn_pose_score_bestpose": cnn_pose_bp,
        "cnn_affinity_bestpose": cnn_aff_bp,
    }


def parse_plapt_json(json_path: str) -> Optional[Dict[str, float]]:
    """Разбор одного plapt .json файла."""
    try:
        with open(json_path, "r") as f:
            data = json.load(f)
    except (OSError, json.JSONDecodeError):
        return None

    result: Dict[str, float] = {}
    if "affinity" in data:
        result["affinity"] = float(data["affinity"])
    if "affinity_uM" in data:
        result["affinity_uM"] = float(data["affinity_uM"])
    return result if result else None


def parse_dynamicbind_filename(fname: str) -> Optional[Tuple[float, float]]:
    """Разбор имени sdf-файла DynamicBind."""
    m = re.search(r"_lddt([0-9.]+)_affinity([0-9.\-]+)\.sdf$", fname)
    if not m:
        return None
    try:
        lddt = float(m.group(1))
        aff = float(m.group(2))
    except ValueError:
        return None
    return lddt, aff


def dynamicbind_ligand_from_folder(folder_name: str) -> str:
    """По имени папки index470_idx_470 делаем идентификатор лиганда idx_470."""
    m = re.search(r"index(\d+)_idx_(\d+)", folder_name)
    if m:
        return f"idx_{m.group(2)}"
    return folder_name


def collect_metrics_from_docking_dir(
    protein: str, docking_dir: Path, dynamicbind_subdir: str = "dynamicbind"
) -> Dict[str, Dict[str, float]]:
    """
    Собирает метрики всех методов докинга из директории docking.
    
    Returns:
        Словарь: {ligand_id: {method_metric: value, ...}, ...}
        где ligand_id = "ligand_1", "ligand_2", ... или "idx_0", "idx_1", ...
    """
    from typing import Tuple
    
    all_metrics: Dict[str, Dict[str, float]] = {}
    
    # QVina
    qvina_dir = docking_dir / "qvina"
    if qvina_dir.exists():
        print(f"   Собираю метрики qvina...")
        count = 0
        for fname in os.listdir(qvina_dir):
            if not fname.endswith(".log") or not fname.startswith("ligand_"):
                continue
            ligand = os.path.splitext(fname)[0]  # ligand_1
            log_path = qvina_dir / fname
            parsed = parse_qvina_log(str(log_path))
            if parsed:
                if ligand not in all_metrics:
                    all_metrics[ligand] = {}
                for k, v in parsed.items():
                    all_metrics[ligand][f"qvina_{k}"] = v
                count += 1
        print(f"     Найдено {count} метрик qvina")
    
    # GNina
    gnina_dir = docking_dir / "gnina"
    if gnina_dir.exists():
        print(f"   Собираю метрики gnina...")
        count = 0
        for fname in os.listdir(gnina_dir):
            if not fname.endswith(".log") or not fname.startswith("ligand_"):
                continue
            ligand = os.path.splitext(fname)[0]  # ligand_1
            log_path = gnina_dir / fname
            parsed = parse_gnina_log(str(log_path))
            if parsed:
                if ligand not in all_metrics:
                    all_metrics[ligand] = {}
                for k, v in parsed.items():
                    all_metrics[ligand][f"gnina_{k}"] = v
                count += 1
        print(f"     Найдено {count} метрик gnina")
    
    # PLAPT
    plapt_dir = docking_dir / "plapt"
    if plapt_dir.exists():
        print(f"   Собираю метрики plapt...")
        count = 0
        for sub in os.listdir(plapt_dir):
            sub_dir = plapt_dir / sub
            if not sub_dir.is_dir():
                continue
            for fname in os.listdir(sub_dir):
                if not fname.endswith(".json") or not fname.startswith("ligand_"):
                    continue
                ligand = os.path.splitext(fname)[0]  # ligand_1
                jpath = sub_dir / fname
                parsed = parse_plapt_json(str(jpath))
                if parsed:
                    if ligand not in all_metrics:
                        all_metrics[ligand] = {}
                    for k, v in parsed.items():
                        all_metrics[ligand][f"plapt_{k}"] = v
                    count += 1
        print(f"     Найдено {count} метрик plapt")
    
    # DynamicBind
    dynamicbind_dir = docking_dir / dynamicbind_subdir
    if dynamicbind_dir.exists():
        print(f"   Собираю метрики dynamicbind...")
        count = 0
        for dataset in os.listdir(dynamicbind_dir):
            ds_dir = dynamicbind_dir / dataset
            if not ds_dir.is_dir():
                continue
            for pd_name in os.listdir(ds_dir):
                pd_dir = ds_dir / pd_name
                if not pd_dir.is_dir():
                    continue
                for folder in os.listdir(pd_dir):
                    idx_dir = pd_dir / folder
                    if not idx_dir.is_dir() or not folder.startswith("index"):
                        continue
                    
                    ligand = dynamicbind_ligand_from_folder(folder)  # idx_0, idx_1, ...
                    
                    lddt_vals: List[float] = []
                    aff_vals: List[float] = []
                    
                    for fname in os.listdir(idx_dir):
                        if not fname.endswith(".sdf") or "lddt" not in fname or "affinity" not in fname:
                            continue
                        parsed = parse_dynamicbind_filename(fname)
                        if parsed:
                            lddt, aff = parsed
                            lddt_vals.append(lddt)
                            aff_vals.append(aff)
                    
                    if lddt_vals:
                        max_aff = max(aff_vals)
                        idx_max_aff = aff_vals.index(max_aff)
                        lddt_at_max_aff = lddt_vals[idx_max_aff]
                        max_lddt = max(lddt_vals)
                        idx_max_lddt = lddt_vals.index(max_lddt)
                        aff_at_max_lddt = aff_vals[idx_max_lddt]
                        
                        if ligand not in all_metrics:
                            all_metrics[ligand] = {}
                        all_metrics[ligand]["dynamicbind_affinity_maxaff"] = max_aff
                        all_metrics[ligand]["dynamicbind_lddt_maxaff"] = lddt_at_max_aff
                        all_metrics[ligand]["dynamicbind_lddt_bestpose"] = max_lddt
                        all_metrics[ligand]["dynamicbind_affinity_bestpose"] = aff_at_max_lddt
                        count += 1
                # Fallback: если нет папок index*_idx_*, читаем affinity из CSV в pd_dir
                for csv_name in ("affinity_prediction.csv", "complete_affinity_prediction.csv"):
                    csv_path = pd_dir / csv_name
                    if not csv_path.exists():
                        continue
                    try:
                        df_csv = pd.read_csv(csv_path)
                    except Exception:
                        continue
                    aff_col = None
                    for c in ("affinity", "predicted_affinity", "affinity_pred"):
                        if c in df_csv.columns:
                            aff_col = c
                            break
                    if aff_col is None:
                        continue
                    for i, row in df_csv.iterrows():
                        try:
                            aff = float(row[aff_col])
                        except (TypeError, ValueError):
                            continue
                        ligand = f"idx_{i}"
                        if ligand not in all_metrics:
                            all_metrics[ligand] = {}
                        all_metrics[ligand]["dynamicbind_affinity_maxaff"] = aff
                        all_metrics[ligand]["dynamicbind_affinity_bestpose"] = aff
                        all_metrics[ligand]["dynamicbind_lddt_maxaff"] = float("nan")
                        all_metrics[ligand]["dynamicbind_lddt_bestpose"] = float("nan")
                        count += 1
                    break  # один CSV на pd_dir достаточно
        print(f"     Найдено {count} метрик dynamicbind")
    
    return all_metrics


def merge_ligands_with_docking(
    ligands_csv: Path,
    protein: str,
    base_dir: Path,
    boltz_results_dir: str = BOLTZ_RESULTS_DIR,
    use_metrics_csv: bool = True,
    dynamicbind_subdir: str = "dynamicbind",
    results_dir: Path | None = None,
) -> pd.DataFrame:
    """
    Объединяет данные лигандов с метриками докинга (boltz2 и exp_*).
    
    Args:
        ligands_csv: Путь к CSV файлу из ligands_nodubl
        protein: PDB ID белка
        base_dir: Базовая директория проекта
        boltz_results_dir: Директория с результатами Boltz
    
    Returns:
        Объединенный DataFrame
    """
    print(f"📖 Читаю лиганды из: {ligands_csv.name}")
    
    # Читаем ligands CSV (разделитель ;)
    try:
        ligands_df = pd.read_csv(ligands_csv, sep=";", low_memory=False)
    except Exception as e:
        print(f"❌ Ошибка при чтении {ligands_csv}: {e}")
        sys.exit(1)
    
    # Проверяем наличие нужных столбцов
    missing_cols = [col for col in LIGANDS_COLUMNS if col not in ligands_df.columns]
    if missing_cols:
        print(f"⚠️  Отсутствуют столбцы в ligands CSV: {missing_cols}")
        available_cols = [col for col in LIGANDS_COLUMNS if col in ligands_df.columns]
        ligands_selected = ligands_df[available_cols].copy()
    else:
        ligands_selected = ligands_df[LIGANDS_COLUMNS].copy()
    
    # Добавляем molecule_chembl_id, если его нет, но он есть в исходном файле
    if "molecule_chembl_id" not in ligands_selected.columns and "molecule_chembl_id" in ligands_df.columns:
        ligands_selected["molecule_chembl_id"] = ligands_df["molecule_chembl_id"]
    
    print(f"   Найдено {len(ligands_selected)} записей лигандов")
    
    # Collect docking metrics from per-target raw results
    if results_dir is None:
        results_dir = base_dir / "results"
    docking_dir = results_dir / protein / "docking"
    print(f"🔍 Собираю метрики докинга из исходных файлов...")
    all_metrics_by_ligand = collect_metrics_from_docking_dir(
        protein, docking_dir, dynamicbind_subdir=dynamicbind_subdir
    )
    print(f"   Найдено метрик для {len(all_metrics_by_ligand)} лигандов")
    
    # Собираем метрики Boltz из JSON файлов
    print(f"🔍 Собираю метрики Boltz для {protein}...")
    boltz_metrics = collect_boltz_metrics(protein, boltz_results_dir)
    print(f"   Найдено {len(boltz_metrics)} метрик Boltz")
    
    # Создаем DataFrame с метриками Boltz
    if boltz_metrics:
        boltz_data = []
        for chembl_id, metrics in boltz_metrics.items():
            boltz_data.append({"molecule_chembl_id": chembl_id, **metrics})
        boltz_df = pd.DataFrame(boltz_data)
    else:
        # Создаем пустой DataFrame с нужными колонками
        boltz_df = pd.DataFrame(columns=["molecule_chembl_id"] + BOLTZ_COLUMNS)
    
    # Преобразуем метрики в DataFrame и маппим к molecule_chembl_id
    # Маппинг: ligand_1 -> индекс 0, ligand_2 -> индекс 1, ... (для qvina, gnina, plapt)
    # Маппинг: idx_0 -> индекс 0, idx_1 -> индекс 1, ... (для dynamicbind)
    
    all_metrics_dfs = {}
    if all_metrics_by_ligand and "molecule_chembl_id" in ligands_selected.columns:
        # Создаем маппинг ligand_id -> csv_index -> molecule_chembl_id
        ligands_for_merge = ligands_selected[["molecule_chembl_id"]].copy()
        ligands_for_merge["csv_index"] = ligands_for_merge.index
        
        # Преобразуем метрики в список словарей
        metrics_rows = []
        for ligand_id, metrics in all_metrics_by_ligand.items():
            # Определяем индекс в CSV
            if ligand_id.startswith("ligand_"):
                # ligand_1 -> индекс 0, ligand_2 -> индекс 1, ...
                match = re.search(r'ligand_(\d+)$', ligand_id)
                if match:
                    ligand_num = int(match.group(1))
                    csv_index = ligand_num - 1
                else:
                    continue
            elif ligand_id.startswith("idx_"):
                # idx_0 -> индекс 0, idx_1 -> индекс 1, ...
                match = re.search(r'idx_(\d+)$', ligand_id)
                if match:
                    csv_index = int(match.group(1))
                else:
                    continue
            else:
                continue
            
            if csv_index >= 0 and csv_index < len(ligands_for_merge):
                row = {"csv_index": csv_index, **metrics}
                metrics_rows.append(row)
        
        if metrics_rows:
            metrics_df = pd.DataFrame(metrics_rows)
            # Мерджим molecule_chembl_id
            metrics_df = metrics_df.merge(
                ligands_for_merge,
                on="csv_index",
                how="left",
            )
            metrics_df = metrics_df.drop(columns=["csv_index"])
            all_metrics_dfs = {"all_methods": metrics_df}
    
    # Добавляем экспериментальные данные из ligands_df
    # exp_value = value, exp_standard_value = standard_value, exp_pchembl_value = pchembl_value
    ligands_selected["exp_value"] = ligands_selected.get("value", None)
    ligands_selected["exp_standard_value"] = ligands_selected.get("standard_value", None)
    ligands_selected["exp_pchembl_value"] = ligands_selected.get("pchembl_value", None)
    
    # Объединяем все метрики
    print(f"🔗 Объединяю данные...")
    
    merged_df = ligands_selected.copy()
    
    # Объединяем метрики других методов
    for method, method_df in all_metrics_dfs.items():
        if "molecule_chembl_id" in method_df.columns and "molecule_chembl_id" in merged_df.columns:
            merged_df = merged_df.merge(
                method_df,
                on="molecule_chembl_id",
                how="left",
                suffixes=("", f"_{method}"),
            )
    
    # Объединяем с метриками Boltz по molecule_chembl_id
    if "molecule_chembl_id" in merged_df.columns and "molecule_chembl_id" in boltz_df.columns:
        merged_df = merged_df.merge(
            boltz_df,
            on="molecule_chembl_id",
            how="left",
            suffixes=("", "_boltz"),
        )
    elif "canonical_smiles" in merged_df.columns:
        # Если нет molecule_chembl_id, объединяем по canonical_smiles
        if "canonical_smiles" in ligands_selected.columns:
            # Нужно сначала добавить canonical_smiles к boltz_df через molecule_chembl_id
            if "molecule_chembl_id" in ligands_selected.columns and "molecule_chembl_id" in boltz_df.columns:
                # Создаем маппинг molecule_chembl_id -> canonical_smiles
                smiles_mapping = ligands_selected[["molecule_chembl_id", "canonical_smiles"]].drop_duplicates()
                boltz_df_with_smiles = boltz_df.merge(
                    smiles_mapping,
                    on="molecule_chembl_id",
                    how="left",
                )
                # Группируем boltz метрики по SMILES (если один SMILES имеет несколько molecule_chembl_id)
                boltz_df_grouped = group_by_smiles(boltz_df_with_smiles, smiles_col="canonical_smiles")
                merged_df = ligands_selected.merge(
                    boltz_df_grouped.drop(columns=["molecule_chembl_id"]),
                    on="canonical_smiles",
                    how="left",
                    suffixes=("", "_boltz"),
                )
            else:
                merged_df = ligands_selected.copy()
        else:
            merged_df = ligands_selected.copy()
    
    # Группируем по canonical_smiles, если есть дубликаты
    if "canonical_smiles" in merged_df.columns:
        initial_len = len(merged_df)
        merged_df = group_by_smiles(merged_df, smiles_col="canonical_smiles")
        if len(merged_df) < initial_len:
            print(f"   После группировки по SMILES: {len(merged_df)} уникальных SMILES (было {initial_len})")
    
    # Добавляем PDB ID как столбец
    merged_df["protein_pdb_id"] = protein.upper()
    
    # Статистика объединения
    total_ligands = len(merged_df)
    boltz_cols_present = [col for col in BOLTZ_COLUMNS if col in merged_df.columns]
    if boltz_cols_present:
        has_boltz_data = merged_df[boltz_cols_present[0]].notna().sum()
    else:
        has_boltz_data = 0
    
    print(f"   Всего лигандов: {total_ligands}")
    print(f"   Лигандов с данными Boltz: {has_boltz_data}")
    print(f"   Лигандов без данных Boltz: {total_ligands - has_boltz_data}")
    
    return merged_df


def process_single_protein(
    pdb_id: str,
    base_dir: Path,
    output_dir: Optional[Path] = None,
    boltz_results_dir: str = BOLTZ_RESULTS_DIR,
    dynamicbind_subdir: str = "dynamicbind_new",
    results_dir: Path | None = None,
) -> bool:
    """
    Обрабатывает один белок: объединяет лиганды с метриками докинга.
    
    Args:
        pdb_id: PDB ID белка (например, "3eyg")
        base_dir: Базовая директория проекта
        output_dir: Директория для сохранения результатов
        boltz_results_dir: Директория с результатами Boltz
    
    Returns:
        True если успешно, False иначе
    """
    pdb_id_lower = pdb_id.lower()
    
    # Проверяем маппинг
    if pdb_id_lower not in PDB_TO_CSV_MAPPING:
        print(f"❌ PDB ID {pdb_id} не найден в маппинге")
        return False
    
    csv_name = PDB_TO_CSV_MAPPING[pdb_id_lower]
    ligands_csv = base_dir / "input" / "ligands_nodubl" / f"{csv_name}_nodubl.csv"
    
    if not ligands_csv.exists():
        print(f"❌ Файл лигандов не найден: {ligands_csv}")
        return False
    
    # Output directory for merged tables
    if output_dir is None:
        output_dir = base_dir / "analysis" / "tables"
    
    # Создаем выходной файл
    output_csv = output_dir / f"merged_ligands_docking_{pdb_id_lower}.csv"
    
    # Объединяем данные
    try:
        merged_df = merge_ligands_with_docking(
            ligands_csv,
            pdb_id_lower,
            base_dir,
            boltz_results_dir,
            dynamicbind_subdir=dynamicbind_subdir,
            results_dir=results_dir,
        )
        
        # Сохраняем результат
        print(f"💾 Сохраняю результат в: {output_csv}")
        merged_df.to_csv(output_csv, index=False, sep=",")
        print(f"✅ Готово! Сохранено {len(merged_df)} записей")
        return True
    except Exception as e:
        print(f"❌ Ошибка при обработке {pdb_id}: {e}")
        import traceback
        traceback.print_exc()
        return False


def process_all_proteins(
    base_dir: Path,
    output_dir: Optional[Path] = None,
    boltz_results_dir: str = BOLTZ_RESULTS_DIR,
    dynamicbind_subdir: str = "dynamicbind_new",
    results_dir: Path | None = None,
) -> None:
    """
    Обрабатывает все белки из маппинга.
    
    Args:
        base_dir: Базовая директория проекта
        output_dir: Директория для сохранения результатов
        boltz_results_dir: Директория с результатами Boltz
    """
    print("=" * 80)
    print("ОБЪЕДИНЕНИЕ ДАННЫХ ЛИГАНДОВ С МЕТРИКАМИ ДОКИНГА (ИЗ ДИРЕКТОРИИ DOCKING)")
    print("=" * 80)
    print()
    
    if output_dir is None:
        output_dir = base_dir / "analysis" / "tables"
    
    success_count = 0
    fail_count = 0
    
    for pdb_id in sorted(PDB_TO_CSV_MAPPING.keys()):
        print(f"\n{'=' * 80}")
        print(f"Обработка белка: {pdb_id.upper()}")
        print(f"{'=' * 80}")
        
        if process_single_protein(
            pdb_id,
            base_dir,
            output_dir,
            boltz_results_dir,
            dynamicbind_subdir=dynamicbind_subdir,
            results_dir=results_dir,
        ):
            success_count += 1
        else:
            fail_count += 1
    
    print(f"\n{'=' * 80}")
    print("ИТОГИ")
    print(f"{'=' * 80}")
    print(f"✅ Успешно обработано: {success_count}")
    print(f"❌ Ошибок: {fail_count}")
    print(f"📁 Результаты сохранены в: {output_dir}")


def main():
    parser = argparse.ArgumentParser(
        description="Объединяет данные лигандов из ligands_nodubl с метриками докинга из директории docking"
    )
    parser.add_argument(
        "pdb_id",
        nargs="?",
        type=str,
        help="PDB ID белка для обработки (например, 3eyg). Если не указан, обрабатываются все белки.",
    )
    parser.add_argument(
        "--base-dir",
        type=str,
        default=None,
        help="Базовая директория проекта (по умолчанию: директория скрипта)",
    )
    parser.add_argument(
        "--output-dir",
        type=str,
        default=None,
        help="Директория для сохранения результатов (по умолчанию: analysis/tables)",
    )
    parser.add_argument(
        "--results-dir",
        type=str,
        default=None,
        help="Корневая директория сырых результатов докинга (по умолчанию: <base-dir>/results)",
    )
    parser.add_argument(
        "--boltz-results-dir",
        type=str,
        default=BOLTZ_RESULTS_DIR or None,
        help="Директория с результатами Boltz (env: BOLTZ_RESULTS_DIR)",
    )
    parser.add_argument(
        "--dynamicbind-subdir",
        type=str,
        default="dynamicbind_new",
        help="Подпапка DynamicBind внутри docking (по умолчанию: dynamicbind_new)",
    )
    
    args = parser.parse_args()
    
    if args.base_dir:
        base_dir = Path(args.base_dir)
    else:
        base_dir = Path(__file__).resolve().parent.parent.parent
    
    output_dir = Path(args.output_dir) if args.output_dir else None
    results_dir = Path(args.results_dir) if args.results_dir else None
    boltz_dir = args.boltz_results_dir or os.environ.get("BOLTZ_RESULTS_DIR", "")
    
    if args.pdb_id:
        success = process_single_protein(
            args.pdb_id,
            base_dir,
            output_dir,
            boltz_dir,
            dynamicbind_subdir=args.dynamicbind_subdir,
            results_dir=results_dir,
        )
        sys.exit(0 if success else 1)
    else:
        process_all_proteins(
            base_dir,
            output_dir,
            boltz_dir,
            dynamicbind_subdir=args.dynamicbind_subdir,
            results_dir=results_dir,
        )


if __name__ == "__main__":
    main()

