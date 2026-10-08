#!/usr/bin/env bash
#SBATCH --job-name=ogp-ref019-fixed
#SBATCH --partition=cu
#SBATCH --account=mselab
#SBATCH --qos=normal
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=01:00:00
# Same small-diagnostic envelope as v1_original_mean_ssn_targets.sh; no array.
set -euo pipefail
: "${DEV_SOURCE_DIR:?Missing immutable source}"
: "${DEV_RUN_DIR:?Missing independent run directory}"
: "${DEV_CONDA:?Missing environment manager}"
: "${DEV_CONDA_ENV:?Missing environment}"
ORIGIN=${1:?Provide frozen job1410868 root}
V1_RUN=${2:?Provide completed job1411083 root}
if [[ ${OGP_FIXED_READY:-0} != 1 ]]; then
    exec "$DEV_CONDA" run -n "$DEV_CONDA_ENV" env OGP_FIXED_READY=1 bash "$0" "$@"
fi
export PYTHONPATH="$DEV_SOURCE_DIR/src:$DEV_SOURCE_DIR"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
"$DEV_CONDA" list -n "$DEV_CONDA_ENV" > "$DEV_RUN_DIR/environment.txt"
/usr/bin/time -v -o "$DEV_RUN_DIR/fixed-point.time" \
    python -m benchmarks.og_extraction.refog019_fixed_gamma \
    --origin "$ORIGIN" --v1-run "$V1_RUN" --out "$DEV_RUN_DIR/fixed-point" \
    > "$DEV_RUN_DIR/fixed-point.stdout" 2> "$DEV_RUN_DIR/fixed-point.stderr"
