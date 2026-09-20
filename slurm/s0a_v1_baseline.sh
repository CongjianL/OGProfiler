#!/usr/bin/env bash
# S0a: V1 full baseline (search + SSN edges).
#
# Part of the 5-submission S0 workflow: s0p (Prodigal) -> s0a/s0b/s0c -> s0d.
# All share one results root: $(dirname "$DEV_RUN_DIR")/s0_v1_v2_of_regression

#SBATCH --job-name=ogp-s0a-v1
#SBATCH --partition=cu
#SBATCH --account=mselab
#SBATCH --qos=normal
#SBATCH --cpus-per-task=56
#SBATCH --mem=128G
#SBATCH --time=24:00:00

set -euo pipefail

: "${DEV_SOURCE_DIR:?Missing immutable source snapshot}"
: "${DEV_RUN_DIR:?Missing immutable run directory}"
: "${DEV_CONDA:?Missing micromamba executable from project remote config}"

ENV="ogprofiler-v1-first"
DATA_NAME="real_embleya"
SHARED_ROOT="$(dirname "$DEV_RUN_DIR")/s0_v1_v2_of_regression"
INPUT_DIR="$SHARED_ROOT/proteomes"
TASK_ROOT="$SHARED_ROOT/results/$DATA_NAME"
V1_ROOT="$TASK_ROOT/v1_full"
mkdir -p "$V1_ROOT"

"$DEV_CONDA" list -n "$ENV" > "$TASK_ROOT/v1-conda-list.txt" 2>&1 || true
{
  "$DEV_CONDA" run -n "$ENV" python --version || true
  "$DEV_CONDA" run -n "$ENV" diamond version || true
} > "$TASK_ROOT/v1-tool-versions.txt" 2>&1 || true

echo "==> [$DATA_NAME] V1 full baseline (env=$ENV)"
/usr/bin/time -v -o "$TASK_ROOT/v1.time" \
  "$DEV_CONDA" run -n "$ENV" python "$DEV_SOURCE_DIR/legacy/OGProfiler_v1.py" \
    --in "$INPUT_DIR" \
    --out "$V1_ROOT" \
    --extension faa \
    --search_method diamond \
    --evalue 0.001 \
    --threads "$SLURM_CPUS_PER_TASK" \
    --distance lrb \
    --weight NBS \
    --community_method rber \
    --network_threads "$SLURM_CPUS_PER_TASK" \
    --gamma_coefficient 1.0 \
    --species_overlap 0 \
  > "$TASK_ROOT/v1.stdout" 2> "$TASK_ROOT/v1.stderr"

echo "==> V1 baseline complete: $V1_ROOT/WorkingDirectory"
