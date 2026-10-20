#!/usr/bin/env bash
#SBATCH --job-name=ogp-qfo-bacteria-evaluate
#SBATCH --partition=cu
#SBATCH --account=mselab
#SBATCH --qos=normal
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=02:00:00
set -euo pipefail
: "${DEV_SOURCE_DIR:?}" "${DEV_RUN_DIR:?}" "${DEV_CONDA:?}" "${DEV_CONDA_ENV:?}"
if [[ ${QFO_EVAL_READY:-0} != 1 ]]; then
    exec "$DEV_CONDA" run -n "$DEV_CONDA_ENV" env QFO_EVAL_READY=1 bash "$0" "$@"
fi
export PYTHONPATH="$DEV_SOURCE_DIR/src:$DEV_SOURCE_DIR"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
export BENCHMARK_THREADS=4 BENCHMARK_SEED=42
"$DEV_SOURCE_DIR/OGProfiler2_benchmark/workflows/run_timed.sh" OG_concordance QFO_bacteria 1 "$DEV_RUN_DIR/evaluation-timing" -- \
    python -m benchmarks.qfo.evaluate --v2 "${1:?}" --of3 "${2:?}" --v1 "${3:?}" --out "$DEV_RUN_DIR/metrics"
