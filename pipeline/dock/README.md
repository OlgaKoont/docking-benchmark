# Pipeline Dock Stage

This directory contains wrappers for executing configured docking/prediction methods on prepared assets.  
It is intended for users who need comparable per-method outputs under a common directory contract (`results/<target>/docking/<method>/...`).

## Repository layout

```text
pipeline/dock/
└── run_docking.sh            docking stage wrapper
```

## Minimal run

```bash
source config/project.env.sh
bash pipeline/dock/run_docking.sh
```

Equivalent top-level call:

```bash
bash run_pipeline.sh dock
```

## Output contract

| Output path | Meaning | Use case |
|---|---|---|
| `results/<target>/docking/gnina/*` | GNINA outputs/logs | structure-based docking benchmark |
| `results/<target>/docking/qvina/*` | QVina outputs/logs | classical docking benchmark |
| `results/<target>/docking/plapt/*` | PLAPT predictions | sequence+SMILES affinity prediction |
| `results/<target>/docking/dynamicbind*/*` | DynamicBind outputs/logs | DL pose/affinity inference |
| `results/<target>/metrics/*_metrics.csv` | extracted per-method metrics | standardized downstream merge/analysis |

## Interpretation scope

- Raw docking scores are method-specific and not directly comparable without normalization/statistical analysis.
- Missing outputs for some ligands/methods are possible; confirm completeness before interpretation.

## Reproducibility

Method behavior is controlled primarily by `config/methods_config.yaml`.

Common parameters and effect on results:
- `gnina.exhaustiveness`, `qvina.exhaustiveness`  
  - Higher values increase search depth/runtime and can change best poses/scores.
- `qvina.num_modes`, `qvina.energy_range`  
  - Controls number/range of retained poses; impacts score extraction and pose diversity.
- `*_random_seed`  
  - Improves run-to-run stability; changing seed can alter stochastic outputs.
- `plapt.device` (`cpu`/`cuda`) and `batch_size`  
  - Mainly runtime/performance stability; numeric drift may occur across hardware backends.
- `dynamicbind` sampling/steps parameters  
  - Directly affect generated pose distribution, affinity predictions, and runtime.
- `docking_timeout` per method  
  - Prevents hangs; too strict timeout increases failed/missing ligand outputs.

Examples:

```bash
# Run only selected methods via project env override
METHODS="gnina qvina" bash pipeline/dock/run_docking.sh
```

```bash
# Keep default methods, but swap methods config
METHODS_CONFIG=config/methods_config.yaml bash pipeline/dock/run_docking.sh
```

## Evaluation status

- Implemented now:
  - wrappers for GNINA, QVina, PLAPT, DynamicBind integration in one stage contract.
- Planned / evolving:
  - extended adapters and richer confidence extraction fields.

## Citation, license, contribution

### Citation

Cite each method you execute (e.g., GNINA, Vina/QVina, PLAPT, DynamicBind, Boltz-2 where applicable).

### License

Respect each method's upstream license terms when running/redistributing outputs.

### Contribution

- Add new method wrappers with clear output schema and metric extraction mapping.
- Document every parameter that can change ranking, score scale, or runtime-completeness tradeoff.
