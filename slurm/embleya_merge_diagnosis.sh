#!/usr/bin/env bash
#SBATCH --job-name=ogp-graph-merge-diagnosis
#SBATCH --partition=cu
#SBATCH --account=mselab
#SBATCH --qos=normal
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=02:00:00
set -euo pipefail
: "${DEV_SOURCE_DIR:?}" "${DEV_RUN_DIR:?}" "${DEV_CONDA:?}" "${DEV_CONDA_ENV:?}"
RUN=${1:?Frozen run required}
REFERENCE=${2:?Reference directory required}
BATCH=${3:?Frozen baseline batch required}
if [[ ${OGP_GRAPH_READY:-0} != 1 ]]; then
    exec "$DEV_CONDA" run -n "$DEV_CONDA_ENV" env OGP_GRAPH_READY=1 bash "$0" "$RUN" "$REFERENCE" "$BATCH"
fi
export PYTHONPATH="$DEV_SOURCE_DIR/src:$DEV_SOURCE_DIR"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
python -m pytest -q "$DEV_SOURCE_DIR/tests/test_graph_merge_diagnosis.py" > "$DEV_RUN_DIR/smoke.txt" 2>&1
"$DEV_CONDA" list -n "$DEV_CONDA_ENV" > "$DEV_RUN_DIR/environment.txt"
/usr/bin/time -v -o "$DEV_RUN_DIR/diagnosis.time" \
    python -m benchmarks.og_extraction.graph_merge_diagnosis \
    --run "$RUN" --reference-dir "$REFERENCE" --out "$DEV_RUN_DIR/merge-diagnosis" \
    --batch "$BATCH" > "$DEV_RUN_DIR/diagnosis.stdout" 2> "$DEV_RUN_DIR/diagnosis.stderr"
