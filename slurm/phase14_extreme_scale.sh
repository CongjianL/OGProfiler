#!/usr/bin/env bash
#SBATCH --job-name=ogp-phase14
#SBATCH --partition=cu
#SBATCH --account=mselab
#SBATCH --qos=normal
#SBATCH --array=0-1%2
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=00:30:00

set -euo pipefail

: "${DEV_SOURCE_DIR:?Missing immutable source snapshot}"
: "${DEV_RUN_DIR:?Missing immutable run directory}"
: "${DEV_CONDA:?Missing Conda/Micromamba executable}"
: "${DEV_CONDA_ENV:?Missing project environment}"
: "${SLURM_ARRAY_TASK_ID:?This script must run as a Slurm array}"

DATASETS=(C_gene_family_expansion E_large_connected_component)
DATASET="${DATASETS[$SLURM_ARRAY_TASK_ID]}"
TASK_ROOT="$DEV_RUN_DIR/results/$DATASET"
RUN_ROOT="$TASK_ROOT/run"
REFERENCE="$DEV_SOURCE_DIR/benchmarks/datasets/$DATASET/reference_ssn.gml"
TASK_MANIFEST="$TASK_ROOT/component-tasks.parquet"
mkdir -p "$TASK_ROOT"
cd "$DEV_SOURCE_DIR"
export PYTHONPATH="$DEV_SOURCE_DIR/src"

"$DEV_CONDA" run -n "$DEV_CONDA_ENV" python -m ogprofiler.cli prepare \
  --proteomes "$DEV_SOURCE_DIR/benchmarks/datasets/$DATASET/proteomes" --out "$RUN_ROOT"
"$DEV_CONDA" run -n "$DEV_CONDA_ENV" python \
  benchmarks/prepare_phase5_reference_edges.py --run "$RUN_ROOT" --ssn "$REFERENCE"
"$DEV_CONDA" run -n "$DEV_CONDA_ENV" python -m ogprofiler.cli components \
  --run "$RUN_ROOT" --set components.edge_batch_size=128
"$DEV_CONDA" run -n "$DEV_CONDA_ENV" python \
  benchmarks/build_phase14_task_manifest.py --run "$RUN_ROOT" --out "$TASK_MANIFEST"

/usr/bin/time -v -o "$TASK_ROOT/parallel.time" \
  "$DEV_CONDA" run -n "$DEV_CONDA_ENV" python \
  benchmarks/run_phase14_component_task.py \
  --run "$RUN_ROOT" --manifest "$TASK_MANIFEST" --shard 0 --shards 1 \
  --subtree-workers "$SLURM_CPUS_PER_TASK" --subtree-release-size 1

"$DEV_CONDA" run -n "$DEV_CONDA_ENV" python \
  benchmarks/aggregate_phase14_components.py --run "$RUN_ROOT" \
  --manifest "$TASK_MANIFEST" --out "$TASK_ROOT/component-results.parquet"

/usr/bin/time -v -o "$TASK_ROOT/resume.time" \
  "$DEV_CONDA" run -n "$DEV_CONDA_ENV" python \
  benchmarks/run_phase14_component_task.py \
  --run "$RUN_ROOT" --manifest "$TASK_MANIFEST" --shard 0 --shards 1 \
  --subtree-workers "$SLURM_CPUS_PER_TASK" --subtree-release-size 1

"$DEV_CONDA" run -n "$DEV_CONDA_ENV" python benchmarks/run_phase14_io.py \
  --out "$TASK_ROOT/io" --rows 2000000 --seed "$((42 + SLURM_ARRAY_TASK_ID))"

echo "Phase 14 $DATASET complete"
