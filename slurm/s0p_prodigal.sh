#!/usr/bin/env bash
# S0p: predict proteins from real genome assemblies with Prodigal (Slurm array).
#
# One array task per genome; all genomes are annotated simultaneously.
# Input:  13 Embleya genome assemblies (.fna) under the remote project_data dir.
# Output: 13 protein FASTA (.faa) files in the shared S0 proteomes directory.
# Part of the 5-submission S0 workflow: s0p -> s0a/s0b/s0c -> s0d.

#SBATCH --job-name=ogp-s0p-prodigal
#SBATCH --partition=cu
#SBATCH --account=mselab
#SBATCH --qos=normal
#SBATCH --array=0-12
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=04:00:00

set -euo pipefail

: "${DEV_SOURCE_DIR:?Missing immutable source snapshot}"
: "${DEV_RUN_DIR:?Missing immutable run directory}"
: "${DEV_CONDA:?Missing micromamba executable from project remote config}"
: "${SLURM_ARRAY_TASK_ID:?This script must run as a Slurm array}"

# NOTE: avoid the shell variable name ENV; it collides with the POSIX/cluster
# module-profile variable and gets overwritten inside subshells.
PRODIGAL_ENV="prodigal-2.6.3"
GENOME_DIR="/home/mselab/licj/project_data/OGProfiler/OGProfiler_test_data"
SHARED_ROOT="$(dirname "$DEV_RUN_DIR")/s0_v1_v2_of_regression"
PROT_DIR="$SHARED_ROOT/proteomes"
mkdir -p "$PROT_DIR"

mapfile -t GENOMES < <(find "$GENOME_DIR" -maxdepth 1 -name "*.fna" | sort)
if (( SLURM_ARRAY_TASK_ID < 0 || SLURM_ARRAY_TASK_ID >= ${#GENOMES[@]} )); then
  echo "Invalid array task: $SLURM_ARRAY_TASK_ID" >&2
  exit 2
fi

genome="${GENOMES[$SLURM_ARRAY_TASK_ID]}"
base="$(basename "$genome" .fna)"

echo "==> [array=$SLURM_ARRAY_TASK_ID] Prodigal on $base (env=$PRODIGAL_ENV)"
"$DEV_CONDA" run -n "$PRODIGAL_ENV" prodigal \
  -i "$genome" -a "$PROT_DIR/${base}.faa" -o /dev/null -p single -q

echo "==> [array=$SLURM_ARRAY_TASK_ID] done: $PROT_DIR/${base}.faa"
