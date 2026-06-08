# Pipeline Postprocess Stage

This directory contains wrappers that transform method outputs into analysis-ready tables and publication figures.  
It is intended for users who need reproducible merge, PoseBusters validation, and inferential analysis execution from `results/` to `analysis/`.

## Repository layout

```text
pipeline/postprocess/
├── run_merge.sh              merge raw method outputs + add pValue column
├── run_posebusters.sh        PoseBusters for generated poses
├── run_article_analysis.sh   analysis pipeline entrypoint
└── run_si_full.sh            full SI profile (10k bootstrap/permutation)
```

## Minimal run

```bash
source config/project.env.sh
bash pipeline/postprocess/run_merge.sh
bash pipeline/postprocess/run_posebusters.sh
bash pipeline/postprocess/run_article_analysis.sh
```

Full SI profile:

```bash
bash pipeline/postprocess/run_si_full.sh
```

## Output contract

| Script | Main output | Use case |
|---|---|---|
| `run_merge.sh` | `analysis/tables/merged_ligands_docking_<pdb>.csv` | unified panel-level merge |
| `run_posebusters.sh` | `analysis/tables/posebuster/*.csv` | pose plausibility summary tables |
| `run_article_analysis.sh` | `analysis/tables/*`, `analysis/figures/*` | standard manuscript analysis run |
| `run_si_full.sh` | same as article analysis + SI log profile | high-precision SI reproducibility run |

## Interpretation scope

- Postprocess outputs become the canonical basis for scientific interpretation.
- Statistical significance and effect sizes must be interpreted jointly (not isolated p-values).
- SI/full profile and dev profile should not be mixed in comparative reporting.

## Reproducibility

### Parameter variants and expected effect

Merge / PoseBusters level:
- `BOLTZ_RESULTS_DIR`  
  - If unset, Boltz-2 outputs are excluded from merged and PoseBusters Boltz paths.
- `DYNAMICBIND_POSEBUSTERS_LABEL`  
  - Must match actual DynamicBind subfolder naming; mismatch yields missing rows.
- `RESULTS_DIR`, `PROCESSED_DIR`, `MERGED_DATA_DIR`, `POSEBUSTERS_DIR`  
  - Path contract overrides; wrong values cause partial/empty tables.

Analysis level (`analysis/config/defaults.sh`):
- `TARGETS`, `METHODS`  
  - Subsetting changes sample size and all comparative statistics/plots.
- `N_BOOTSTRAP`, `N_PERMUTATION`  
  - Higher values increase stability of uncertainty/p-values, but increase runtime.
- `RANDOM_SEED`  
  - Controls bootstrap/permutation reproducibility.
- `RUN_CORRELATIONS`, `RUN_ENRICHMENT`, `RUN_POSEBUSTERS`, `RUN_INFERENTIAL`, `RUN_FIGURES`  
  - Stage switches; disabling parts produces incomplete manuscript output sets.
- `FIGURE_DPI`, `METHOD_COLORS`, `COLOR_*`, `FONT_FAMILY`  
  - Visual rendering changes only; statistical tables remain unaffected.

Examples:

```bash
# Fast development run (reduced resampling)
N_BOOTSTRAP=300 N_PERMUTATION=300 bash pipeline/postprocess/run_article_analysis.sh
```

```bash
# Disable figure regeneration, keep tables
RUN_FIGURES=0 bash pipeline/postprocess/run_article_analysis.sh
```

```bash
# Target/method subset run
TARGETS="1g5m 3eyg 6gqj" METHODS="gnina qvina" \
bash pipeline/postprocess/run_article_analysis.sh
```

## Evaluation status

- Implemented now:
  - complete postprocess chain from merge to SI-grade inferential outputs.
- Planned / evolving:
  - additional table harmonization and plotting refinements as endpoint families expand.

## Citation, license, contribution

### Citation

Use the full citation block from `analysis/README.md` for metrics/statistics and method context.

### License

Respect repository and third-party tool/data licensing for generated artifacts.

### Contribution

- Any new postprocess variable must include default, allowed range/options, and impact statement.
- Keep script outputs path-stable for downstream manuscript tooling.
