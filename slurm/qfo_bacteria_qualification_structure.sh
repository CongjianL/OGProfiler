#!/usr/bin/env bash
#SBATCH --job-name=ogp-qfo-structure
#SBATCH --partition=cu
#SBATCH --account=mselab
#SBATCH --qos=normal
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=02:00:00
set -euo pipefail
: "${DEV_SOURCE_DIR:?}" "${DEV_RUN_DIR:?}" "${DEV_CONDA:?}" "${DEV_CONDA_ENV:?}"
V2_ROOT=${1:?}; OF3_ROOT=${2:?}; AUDIT_ROOT=${3:?}
if [[ ${QFO_FIXED_READY:-0} != 1 ]]; then
    exec "$DEV_CONDA" run -n "$DEV_CONDA_ENV" env QFO_FIXED_READY=1 bash "$0" "$@"
fi
export PYTHONPATH="$DEV_SOURCE_DIR/src:$DEV_SOURCE_DIR"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
export BENCHMARK_THREADS=4 BENCHMARK_SEED=42
REFERENCE=$(dirname "$(cat "$OF3_ROOT/campaign/of3-assigned-primary.txt")")
"$DEV_SOURCE_DIR/OGProfiler2_benchmark/workflows/run_timed.sh" Qualification_structure QFO_bacteria 1 "$DEV_RUN_DIR/diagnostic-timing" -- \
    python -m benchmarks.qfo.qualification_structure --run "$V2_ROOT/campaign/v2" \
    --reference-dir "$REFERENCE" --baseline "$AUDIT_ROOT/fixed-tree-audit/report.json" --out "$DEV_RUN_DIR/qualification-structure"
