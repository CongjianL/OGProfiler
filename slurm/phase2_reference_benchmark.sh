#!/usr/bin/env bash
#SBATCH --job-name=ogp-p2-reference
#SBATCH --partition=cu
#SBATCH --account=mselab
#SBATCH --qos=normal
#SBATCH --array=0-5%2
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=01:00:00

set -euo pipefail

: "${DEV_SOURCE_DIR:?Missing immutable source snapshot}"
: "${DEV_RUN_DIR:?Missing immutable run directory}"
: "${DEV_CONDA:?Missing micromamba executable from project remote config}"
: "${DEV_CONDA_ENV:?Missing micromamba environment from project remote config}"
: "${SLURM_ARRAY_TASK_ID:?This script must run as a Slurm array}"

DATASETS=(C_gene_family_expansion E_large_connected_component)
REPEATS=3
DATASET_INDEX=$((SLURM_ARRAY_TASK_ID / REPEATS))
REPEAT_INDEX=$((SLURM_ARRAY_TASK_ID % REPEATS))
if (( DATASET_INDEX < 0 || DATASET_INDEX >= ${#DATASETS[@]} )); then
  echo "Invalid array task: $SLURM_ARRAY_TASK_ID" >&2
  exit 2
fi

DATASET="${DATASETS[$DATASET_INDEX]}"
REPEAT="$(printf '%02d' "$REPEAT_INDEX")"
REFERENCE_SSN="$DEV_SOURCE_DIR/benchmarks/datasets/$DATASET/reference_ssn.gml"
TASK_ROOT="$DEV_RUN_DIR/results/$DATASET/repeat_$REPEAT"
V1_ROOT="$TASK_ROOT/v1_hierarchy"
V1_SUMMARY="$TASK_ROOT/v1_summary"
V2_ROOT="$TASK_ROOT/v2_hierarchy"
COMPARE_ROOT="$TASK_ROOT/comparison"
mkdir -p "$V1_ROOT" "$V1_SUMMARY" "$V2_ROOT" "$COMPARE_ROOT"

cp "$DEV_SOURCE_DIR/benchmarks/datasets/$DATASET/manifest.json" \
  "$TASK_ROOT/dataset_manifest.json"
"$DEV_CONDA" run -n "$DEV_CONDA_ENV" python \
  "$DEV_SOURCE_DIR/benchmarks/prepare_reference_ssn.py" \
  --ssn "$REFERENCE_SSN" --out "$V1_ROOT"

if [[ "$REPEAT_INDEX" == "0" ]]; then
  "$DEV_CONDA" list -n "$DEV_CONDA_ENV" > "$DEV_RUN_DIR/conda-list.txt"
  {
    "$DEV_CONDA" run -n "$DEV_CONDA_ENV" python --version || true
    "$DEV_CONDA" run -n "$DEV_CONDA_ENV" python -c \
      'import igraph, leidenalg; print("igraph", igraph.__version__); print("leidenalg", leidenalg.__version__)' || true
  } > "$TASK_ROOT/tool-versions.txt"
fi

echo "==> [$DATASET repeat=$REPEAT] frozen V1 hierarchy"
/usr/bin/time -v -o "$TASK_ROOT/v1.time" \
  "$DEV_CONDA" run -n "$DEV_CONDA_ENV" python \
  "$DEV_SOURCE_DIR/benchmarks/run_v1_hierarchy.py" \
  --v1 "$DEV_SOURCE_DIR/legacy/OGProfiler_v1.py" \
  --ssn "$V1_ROOT/ssn.gml" \
  --out "$V1_ROOT" \
  --method rber --weight NBS --threads "$SLURM_CPUS_PER_TASK" --gamma-coefficient 1.0 \
  > "$TASK_ROOT/v1.stdout" 2> "$TASK_ROOT/v1.stderr"

"$DEV_CONDA" run -n "$DEV_CONDA_ENV" python \
  "$DEV_SOURCE_DIR/benchmarks/summarize_v1_run.py" \
  --working-dir "$V1_ROOT" --time-report "$TASK_ROOT/v1.time" \
  --stdout "$TASK_ROOT/v1.stdout" --out "$V1_SUMMARY"

echo "==> [$DATASET repeat=$REPEAT] V2 hierarchy"
/usr/bin/time -v -o "$TASK_ROOT/v2.time" \
  "$DEV_CONDA" run -n "$DEV_CONDA_ENV" bash -lc \
  "cd '$DEV_SOURCE_DIR' && PYTHONPATH=src python -m ogprofiler prototype-hierarchy \
    --ssn '$REFERENCE_SSN' --out '$V2_ROOT' --weight-attribute NBS" \
  > "$TASK_ROOT/v2.stdout" 2> "$TASK_ROOT/v2.stderr"

"$DEV_CONDA" run -n "$DEV_CONDA_ENV" python \
  "$DEV_SOURCE_DIR/benchmarks/compare_hierarchies.py" \
  --v1-normalized "$V1_SUMMARY/normalized_hierarchy.json" \
  --v1-hierarchy-metrics "$V1_ROOT/metrics.json" \
  --v2-run "$V2_ROOT" \
  --sequence-ids "$V1_ROOT/SequenceIDs.txt" \
  --out "$COMPARE_ROOT"

echo "==> [$DATASET repeat=$REPEAT] complete"

