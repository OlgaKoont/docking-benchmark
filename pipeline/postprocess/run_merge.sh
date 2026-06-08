#!/usr/bin/env bash
# Stage 3: merge ligand tables, add pValue, optional PoseBusters.
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/../../config/project.env.sh"

POST_PY="${TOXAFFINITY_ROOT}/src/analysis"
mkdir -p "${MERGED_DATA_DIR}" "${POSEBUSTERS_DIR}" "${POSEBUSTERS_BOLTZ2_DIR}"

cd "${TOXAFFINITY_ROOT}"

echo "=== Merge docking metrics -> ${MERGED_DATA_DIR} ==="
"${PYTHON}" "${POST_PY}/merge_ligands_docking_from_dir.py" \
  --base-dir "${TOXAFFINITY_ROOT}" \
  --output-dir "${MERGED_DATA_DIR}" \
  --results-dir "${RESULTS_DIR}" \
  --dynamicbind-subdir dynamicbind_new \
  ${BOLTZ_RESULTS_DIR:+--boltz-results-dir "${BOLTZ_RESULTS_DIR}"} \
  "$@"

echo "=== Add pValue column ==="
"${PYTHON}" "${POST_PY}/add_pvalue_column.py" \
  --input-dir "${MERGED_DATA_DIR}" \
  --output-dir "${MERGED_DATA_DIR}"

echo "Postprocess merge complete."
