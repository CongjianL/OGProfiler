#!/usr/bin/env bash
# S0c: OrthoFinder3 full run (search + orthogroups).
#
# OrthoFinder 3.1.5 handles its own clustering (no separate MCL invocation).
# Part of the 5-submission S0 workflow (s0p -> s0a/s0b/s0c -> s0d); writes into the shared results root.

#SBATCH --job-name=ogp-s0c-orthofinder
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

ENV="orthofinder-3.1.5"
DATA_NAME="real_embleya"
SHARED_ROOT="$(dirname "$DEV_RUN_DIR")/s0_v1_v2_of_regression"
INPUT_DIR="$SHARED_ROOT/proteomes"
TASK_ROOT="$SHARED_ROOT/results/$DATA_NAME"
OF_ROOT="$TASK_ROOT/orthofinder"
mkdir -p "$TASK_ROOT"

"$DEV_CONDA" list -n "$ENV" > "$TASK_ROOT/of-conda-list.txt" 2>&1 || true
{
  "$DEV_CONDA" run -n "$ENV" orthofinder --help 2>&1 | head -1 || true
  "$DEV_CONDA" run -n "$ENV" diamond version || true
} > "$TASK_ROOT/of-tool-versions.txt" 2>&1 || true

echo "==> [$DATA_NAME] OrthoFinder3 full run (env=$ENV)"
/usr/bin/time -v -o "$TASK_ROOT/of.time" \
  "$DEV_CONDA" run -n "$ENV" orthofinder \
    -f "$INPUT_DIR" \
    -o "$OF_ROOT" \
    -S diamond \
    -t "$SLURM_CPUS_PER_TASK" \
    -a "$SLURM_CPUS_PER_TASK" \
  > "$TASK_ROOT/of.stdout" 2> "$TASK_ROOT/of.stderr"

OF_WORKDIR="$(find "$OF_ROOT" -maxdepth 2 -type d -name WorkingDirectory | head -1)"
[[ -n "$OF_WORKDIR" ]] || { echo "OrthoFinder WorkingDirectory not found under $OF_ROOT" >&2; exit 3; }
printf '%s\n' "$OF_WORKDIR" > "$TASK_ROOT/of_workdir.txt"

echo "==> OrthoFinder3 complete; WorkingDirectory: $OF_WORKDIR"
