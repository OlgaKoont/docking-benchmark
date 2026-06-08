# Pipeline Prepare Stage

This directory contains wrappers for converting raw `input/` assets into method-ready artifacts (`processed/`).  
It is intended for users who need reproducible receptor preparation, ligand preparation, and docking box generation before docking.

## Repository layout

```text
pipeline/prepare/
└── run_prepare.sh            preparation stage wrapper
```

## Minimal run

```bash
source config/project.env.sh
bash pipeline/prepare/run_prepare.sh
```

Pass-through CLI options to underlying runner:

```bash
bash pipeline/prepare/run_prepare.sh --config config/toxdock_config.yaml --stage preparation
```

## Output contract

| Output path | Meaning | Use case |
|---|---|---|
| `processed/proteins/*.pdbqt` | prepared receptors | GNINA/QVina docking input |
| `processed/ligands/<target>/<dataset>/*.pdbqt` | prepared ligand sets | method input ligands |
| `processed/boxes/*.json` | docking search boxes | docking grid definition |
| `processed/plapt_sequences/*.txt` | optional sequence artifacts | PLAPT stage support |

## Interpretation scope

- Preparation artifacts are operational intermediates; do not treat them as performance endpoints.
- Different preprocessing settings can materially change downstream docking scores and ranking behavior.

## Reproducibility

Primary parameter files:
- `TOXDOCK_CONFIG` (default `config/toxdock_config.yaml`)
- `METHODS_CONFIG` (default `config/methods_config.yaml`)

Important parameters and effect on results:
- `protein_settings_file` (in `TOXDOCK_CONFIG`)  
  - Changes chain/cofactor/protonation behavior; impacts pocket geometry and docking feasibility.
- `ligand_dir` / ligand curation set  
  - Changes chemical panel composition and all downstream statistics.
- `labox.scale` / `labox.min_size`  
  - Alters search space volume; larger boxes can increase runtime and false positives.
- `random_state`  
  - Affects stochastic components (e.g., ligand embedding); different seeds may shift poses/scores.

Examples:

```bash
# Use alternative preparation config
TOXDOCK_CONFIG=config/toxdock_config.yaml \
bash pipeline/prepare/run_prepare.sh
```

```bash
# Through top-level driver
bash run_pipeline.sh prepare
```

## Evaluation status

- Implemented now:
  - receptor/ligand/box preparation and optional method-specific helper preprocessing.
- Planned / evolving:
  - expanded preprocessing hooks for additional methods.

## Citation, license, contribution

### Citation

Cite method and data sources used in downstream stages; this stage itself is preprocessing orchestration.

### License

Respect licenses of preprocessing dependencies (RDKit, Meeko, etc.).

### Contribution

- Any new preprocessing parameter must be documented with effect direction (runtime, geometry, reproducibility).
