#!/usr/bin/env bash
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/../../config/project.env.sh"

POST_PY="${TOXAFFINITY_ROOT}/src/analysis"
mkdir -p "${POSEBUSTERS_DIR}" "${POSEBUSTERS_BOLTZ2_DIR}"

cd "${TOXAFFINITY_ROOT}"

echo "=== PoseBusters (gnina, qvina, dynamicbind) ==="
"${PYTHON}" "${POST_PY}/prepare_and_run_posebusters.py" \
  --poses-dir "${RESULTS_DIR}" \
  --proteins-dir "${PROCESSED_DIR}/proteins" \
  --output-dir "${ANALYSIS_ROOT}/tables/posebuster" \
  "$@"
if [[ -d "${POSEBUSTERS_DIR}/results_by_method" ]]; then
  mv "${POSEBUSTERS_DIR}/results_by_method/"*.csv "${POSEBUSTERS_DIR}/" 2>/dev/null || true
fi

if [[ -n "${BOLTZ_RESULTS_DIR}" ]]; then
  echo "=== PoseBusters (Boltz-2) ==="
  "${PYTHON}" "${POST_PY}/prepare_and_run_posebusters_boltz2.py" \
    --boltz-results-dir "${BOLTZ_RESULTS_DIR}" \
    --output-dir "${ANALYSIS_ROOT}/tables/posebuster" \
    "$@"
  if [[ -d "${POSEBUSTERS_BOLTZ2_DIR}/results_by_method" ]]; then
    mv "${POSEBUSTERS_BOLTZ2_DIR}/results_by_method/"*.csv "${POSEBUSTERS_BOLTZ2_DIR}/" 2>/dev/null || true
  fi
fi
