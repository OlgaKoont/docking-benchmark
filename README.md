# ToxAffinity

ToxDock-Bench is a toxicity-oriented docking benchmark built on a curated 16-target BindingDB panel.  
The repository provides a reproducible pipeline from structure preparation to docking, postprocessing, pose validation, and manuscript-grade statistical analysis.

## What This Repository Delivers

- A fixed 16-target toxicity panel with harmonized affinity labels.
- Multi-method docking benchmark across classical, hybrid, and deep-learning methods.
- Unified postprocessing outputs (`merged` tables + PoseBusters checks).
- Publication-ready statistics and figures for scoring, ranking, screening, and pose plausibility.

## Quick Start (Analysis Only)

Use this mode when `analysis/tables/` and `analysis/tables/posebuster/` are already populated.

```bash
export TOXAFFINITY_ROOT="$(pwd)"
export PYTHON=python3   # or your conda env python

bash pipeline/postprocess/run_article_analysis.sh
# Full SI profile (10k bootstrap + 10k permutation):
bash pipeline/postprocess/run_si_full.sh
```

Main outputs:
- tables: `analysis/tables/`
- figures: `analysis/figures/`

## Full End-to-End Reproduction

```bash
source config/project.env.sh
export BOLTZ_RESULTS_DIR=/path/to/boltz/data/results   # optional, for Boltz-2 affinity outputs

bash run_pipeline.sh install
bash run_pipeline.sh prepare
bash run_pipeline.sh dock
bash run_pipeline.sh merge
bash run_pipeline.sh posebusters
bash run_pipeline.sh analysis
```

Or run everything in one command:

```bash
bash run_pipeline.sh all
```

## Pipeline Stages

- `install`: environment and package setup.
- `prepare`: protein/ligand/box preparation.
- `dock`: method-specific docking/inference.
- `merge`: docking output merge + `pValue` augmentation.
- `posebusters`: pose validity tables (all supported pose-generating methods).
- `analysis`: statistics and figures via `src/analysis/run_pipeline.py`.

## Repository layout

```text
ToxAffinity/
├── config/                    path contract + YAML configs
├── input/                     raw proteins and ligand tables
├── processed/                 prepared structures and intermediate assets
├── results/                   raw docking outputs ({pdb}/docking/{method}/...)
├── src/
│   ├── docking_benchmark2/    core benchmark package
│   └── analysis/              merge/posebusters/statistics/plotting modules
├── pipeline/
│   ├── install/
│   ├── prepare/
│   ├── dock/
│   └── postprocess/           shell entrypoints calling src/analysis modules
├── analysis/
│   ├── tables/                merged ligand tables + stats CSVs
│   │   └── posebuster/        unified PoseBusters CSVs (all methods incl. Boltz-2)
│   ├── figures/               generated figures (PNG/SVG)
└── run_pipeline.sh            top-level stage driver
```

## Data and Output Contract

| Path | Role |
|---|---|
| `analysis/tables/merged_ligands_docking_<pdb>.csv` | merged ligand + docking + pValue table per target |
| `analysis/tables/posebuster/posebusters_results_<pdb>_<method>.csv` | pose-level quality checks |
| `analysis/tables/correlations/` | Pearson/Spearman/Kendall summaries, CI/permutation |
| `analysis/tables/enrichment/` | EF/nEF summaries |
| `analysis/tables/inferential/` | cross-target and per-target pairwise tests |
| `analysis/figures/` | all generated manuscript figures |

## Key Environment Variables

| Variable | Default |
|---|---|
| `TOXAFFINITY_ROOT` | repo root |
| `RESULTS_DIR` | `$ROOT/results` |
| `ANALYSIS_ROOT` | `$ROOT/analysis` |
| `MERGED_DATA_DIR` | `$ROOT/analysis/tables` |
| `POSEBUSTERS_DIR` | `$ROOT/analysis/tables/posebuster` |
| `POSEBUSTERS_BOLTZ2_DIR` | `$ROOT/analysis/tables/posebuster` |
| `BOLTZ_RESULTS_DIR` | unset (optional external Boltz results) |

Full path contract: `config/project.env.sh`.

## Curated Benchmark Panel

- Targets (16): `1g5m 2z5x 3eyg 3jy9 3lxk 3mjg 4ase 4f65 4tz4 4zau 5jkv 5mo4 6gqj 6jok 7awe 7kk3`
- Methods (5): Boltz-2, DynamicBind, GNINA 1.3, PLAPT, QVina2

## Zenodo Artifacts

Use Zenodo records for heavyweight or archival assets that are not tracked in this Git repository.

| Zenodo link | What it contains | How to use |
|---|---|---|
| `https://zenodo.org/records/<INPUT_ARCHIVE_RECORD_ID>` | archived raw inputs (for example BindingDB snapshots) | optional provenance rebuilds |
| `https://zenodo.org/records/<RESULTS_ARCHIVE_RECORD_ID>` | large raw docking outputs (`results/`) | reproduce merged/postprocess steps without rerunning docking |
| `https://zenodo.org/records/<ANALYSIS_ARCHIVE_RECORD_ID>` | generated analysis artifacts (`analysis/tables`, `analysis/figures`) | manuscript/SI reproduction and audit |

Replace placeholder record IDs with your published Zenodo links.

## Citation

Manuscript in preparation (JCIM).
