#!/usr/bin/env bash
# run.sh — launch the standalone mix-generation pipeline
#
# Usage:
#   bash run.sh [-resume] [-params-file params_test.yaml] [--param value ...]
#
# Defaults to params.yaml; pass -params-file to override (e.g. for params_test.yaml).

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

LOG_DIR="${SCRIPT_DIR}/logs"
mkdir -p "${LOG_DIR}"
TIMESTAMP="$(date +%Y%m%d_%H%M%S)"

PARAMS_FILE="${SCRIPT_DIR}/params.yaml"
if [[ "$*" == *"-params-file"* ]]; then
    PARAMS_FILE=""
fi

cd "${SCRIPT_DIR}"

nextflow -log "${LOG_DIR}/nextflow_${TIMESTAMP}.log" run "${SCRIPT_DIR}/main.nf" \
    -c "${SCRIPT_DIR}/nextflow.config" \
    ${PARAMS_FILE:+-params-file "${PARAMS_FILE}"} \
    "$@"
