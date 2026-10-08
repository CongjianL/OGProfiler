#!/usr/bin/env bash
# Read-only fixed hierarchy audit, no Leiden or scientific parameter changes.
#SBATCH --job-name=ogp-embleya-representation
#SBATCH --partition=cu
#SBATCH --account=mselab
#SBATCH --qos=normal
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=02:00:00
set -euo pipefail
: "${DEV_SOURCE_DIR:?}"
: "${DEV_RUN_DIR:?}"
: "${DEV_CONDA:?}"
: "${DEV_CONDA_ENV:?}"
RUN=${1:?Usage: embleya_representation_audit.sh FROZEN_SOFT42_RUN OF_REFERENCE_DIR}
REFERENCE=${2:?}
if [[ ${OGP_REPRESENTATION_READY:-0} != 1 ]]; then
    exec "$DEV_CONDA" run -n "$DEV_CONDA_ENV" env OGP_REPRESENTATION_READY=1 bash "$0" "$@"
fi
export PYTHONPATH="$DEV_SOURCE_DIR/src:$DEV_SOURCE_DIR"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
python -m pytest -q "$DEV_SOURCE_DIR/tests/test_tree_cut_audit.py" \
    "$DEV_SOURCE_DIR/tests/test_embleya_reference.py" > "$DEV_RUN_DIR/smoke.txt" 2>&1
"$DEV_CONDA" list -n "$DEV_CONDA_ENV" > "$DEV_RUN_DIR/environment.txt"
/usr/bin/time -v -o "$DEV_RUN_DIR/audit.time" \
    python -m benchmarks.og_extraction.embleya_representation \
    --run "$RUN" --reference-dir "$REFERENCE" --out "$DEV_RUN_DIR/representation-audit" \
    > "$DEV_RUN_DIR/audit.stdout" 2> "$DEV_RUN_DIR/audit.stderr"
