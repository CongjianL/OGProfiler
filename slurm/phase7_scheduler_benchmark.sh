#!/usr/bin/env bash
#SBATCH --job-name=ogp-p7-scheduler
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
: "${DEV_CONDA:?Missing micromamba executable}"
: "${DEV_CONDA_ENV:?Missing micromamba environment}"
: "${SLURM_ARRAY_TASK_ID:?This script must run as a Slurm array}"

DATASETS=(C_gene_family_expansion E_large_connected_component)
DATASET="${DATASETS[$SLURM_ARRAY_TASK_ID]}"
TASK_ROOT="$DEV_RUN_DIR/results/$DATASET"
RUN_ROOT="$TASK_ROOT/run"
REFERENCE="$DEV_SOURCE_DIR/benchmarks/datasets/$DATASET/reference_ssn.gml"
mkdir -p "$TASK_ROOT"
export PYTHONPATH="$DEV_SOURCE_DIR/src"

"$DEV_CONDA" run -n "$DEV_CONDA_ENV" python -m ogprofiler.cli prepare \
  --proteomes "$DEV_SOURCE_DIR/benchmarks/datasets/$DATASET/proteomes" --out "$RUN_ROOT"
"$DEV_CONDA" run -n "$DEV_CONDA_ENV" python \
  "$DEV_SOURCE_DIR/benchmarks/prepare_phase5_reference_edges.py" \
  --run "$RUN_ROOT" --ssn "$REFERENCE"
"$DEV_CONDA" run -n "$DEV_CONDA_ENV" python -m ogprofiler.cli components \
  --run "$RUN_ROOT" --set components.edge_batch_size=128

/usr/bin/time -v -o "$TASK_ROOT/initial.time" \
  "$DEV_CONDA" run -n "$DEV_CONDA_ENV" python -m ogprofiler.cli hierarchy-all \
  --run "$RUN_ROOT" --set runtime.workers=4 \
  --set hierarchy.stability_mode=robust \
  > "$TASK_ROOT/initial.stdout" 2> "$TASK_ROOT/initial.stderr"
cp "$RUN_ROOT/hierarchy/scheduler-manifest.json" "$TASK_ROOT/scheduler-initial.json"

/usr/bin/time -v -o "$TASK_ROOT/resume.time" \
  "$DEV_CONDA" run -n "$DEV_CONDA_ENV" python -m ogprofiler.cli hierarchy-all \
  --run "$RUN_ROOT" --set runtime.workers=4 \
  --set hierarchy.stability_mode=robust \
  > "$TASK_ROOT/resume.stdout" 2> "$TASK_ROOT/resume.stderr"
cp "$RUN_ROOT/hierarchy/scheduler-manifest.json" "$TASK_ROOT/scheduler-resume.json"

"$DEV_CONDA" run -n "$DEV_CONDA_ENV" python \
  "$DEV_SOURCE_DIR/benchmarks/summarize_phase7_scheduler.py" \
  --run "$RUN_ROOT" --initial-time "$TASK_ROOT/initial.time" \
  --resume-time "$TASK_ROOT/resume.time" --out "$TASK_ROOT/summary.json"

echo "Phase 7 $DATASET complete"
