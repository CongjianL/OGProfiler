#!/usr/bin/env bash
#SBATCH --job-name=ogp-v1-baseline
#SBATCH --partition=cu
#SBATCH --account=mselab
#SBATCH --qos=normal
#SBATCH --array=0-4%2
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=02:00:00

set -euo pipefail

: "${DEV_SOURCE_DIR:?Missing immutable source snapshot}"
: "${DEV_RUN_DIR:?Missing immutable run directory}"
: "${DEV_CONDA:?Missing micromamba executable from project remote config}"
: "${DEV_CONDA_ENV:?Missing micromamba environment from project remote config}"
: "${SLURM_ARRAY_TASK_ID:?This script must run as a Slurm array}"

DATASETS=(
  A_small_sanity
  B_paralog
  C_gene_family_expansion
  D_fusion_multidomain
  E_large_connected_component
)

if (( SLURM_ARRAY_TASK_ID < 0 || SLURM_ARRAY_TASK_ID >= ${#DATASETS[@]} )); then
  echo "Invalid array task: $SLURM_ARRAY_TASK_ID" >&2
  exit 2
fi

DATASET="${DATASETS[$SLURM_ARRAY_TASK_ID]}"
INPUT_DIR="$DEV_SOURCE_DIR/benchmarks/datasets/$DATASET/proteomes"
DATASET_MANIFEST="$DEV_SOURCE_DIR/benchmarks/datasets/$DATASET/manifest.json"
TASK_ROOT="$DEV_RUN_DIR/results/$DATASET"
FULL_ROOT="$TASK_ROOT/v1_full"
WORKING_DIR="$FULL_ROOT/WorkingDirectory"
V1_HIERARCHY_ROOT="$TASK_ROOT/v1_hierarchy_only"
V2_ROOT="$TASK_ROOT/v2_hierarchy"
COMPARE_ROOT="$TASK_ROOT/comparison"

mkdir -p "$FULL_ROOT" "$V1_HIERARCHY_ROOT" "$V2_ROOT" "$COMPARE_ROOT"

COMMAND=(
  python "$DEV_SOURCE_DIR/legacy/OGProfiler_v1.py"
  --in "$INPUT_DIR"
  --out "$FULL_ROOT"
  --extension faa
  --search_method diamond
  --evalue 0.001
  --threads "$SLURM_CPUS_PER_TASK"
  --distance lrb
  --weight NBS
  --community_method rber
  --network_threads "$SLURM_CPUS_PER_TASK"
  --gamma_coefficient 1.0
  --species_overlap 0
)

printf '%q ' "${COMMAND[@]}" > "$TASK_ROOT/v1_full_command.sh"
printf '\n' >> "$TASK_ROOT/v1_full_command.sh"
cp "$DATASET_MANIFEST" "$TASK_ROOT/dataset_manifest.json"
"$DEV_CONDA" list -n "$DEV_CONDA_ENV" > "$TASK_ROOT/conda-list.txt"
{
  "$DEV_CONDA" run -n "$DEV_CONDA_ENV" python --version || true
  "$DEV_CONDA" run -n "$DEV_CONDA_ENV" diamond version || true
  "$DEV_CONDA" run -n "$DEV_CONDA_ENV" mmseqs version || true
  "$DEV_CONDA" run -n "$DEV_CONDA_ENV" blastp -version || true
  "$DEV_CONDA" run -n "$DEV_CONDA_ENV" FastTree 2>&1 || true
  "$DEV_CONDA" run -n "$DEV_CONDA_ENV" mafft --version 2>&1 || true
} > "$TASK_ROOT/tool-versions.txt"

echo "==> [$DATASET] frozen V1 full baseline"
/usr/bin/time -v -o "$TASK_ROOT/v1_full.time" \
  "$DEV_CONDA" run -n "$DEV_CONDA_ENV" "${COMMAND[@]}" \
  > "$TASK_ROOT/v1_full.stdout" 2> "$TASK_ROOT/v1_full.stderr"

"$DEV_CONDA" run -n "$DEV_CONDA_ENV" python \
  "$DEV_SOURCE_DIR/benchmarks/summarize_v1_run.py" \
  --working-dir "$WORKING_DIR" \
  --time-report "$TASK_ROOT/v1_full.time" \
  --stdout "$TASK_ROOT/v1_full.stdout" \
  --out "$TASK_ROOT/v1_full_summary"

cp "$WORKING_DIR/ssn.gml" "$V1_HIERARCHY_ROOT/ssn.gml"
cp "$WORKING_DIR/SequenceIDs.txt" "$V1_HIERARCHY_ROOT/SequenceIDs.txt"

echo "==> [$DATASET] frozen V1 hierarchy-only baseline"
/usr/bin/time -v -o "$TASK_ROOT/v1_hierarchy.time" \
  "$DEV_CONDA" run -n "$DEV_CONDA_ENV" python \
  "$DEV_SOURCE_DIR/benchmarks/run_v1_hierarchy.py" \
  --v1 "$DEV_SOURCE_DIR/legacy/OGProfiler_v1.py" \
  --ssn "$WORKING_DIR/ssn.gml" \
  --out "$V1_HIERARCHY_ROOT" \
  --method rber \
  --weight NBS \
  --threads "$SLURM_CPUS_PER_TASK" \
  --gamma-coefficient 1.0 \
  > "$TASK_ROOT/v1_hierarchy.stdout" 2> "$TASK_ROOT/v1_hierarchy.stderr"

"$DEV_CONDA" run -n "$DEV_CONDA_ENV" python \
  "$DEV_SOURCE_DIR/benchmarks/summarize_v1_run.py" \
  --working-dir "$V1_HIERARCHY_ROOT" \
  --time-report "$TASK_ROOT/v1_hierarchy.time" \
  --stdout "$TASK_ROOT/v1_hierarchy.stdout" \
  --out "$TASK_ROOT/v1_hierarchy_summary"

echo "==> [$DATASET] V2 hierarchy on the identical V1 SSN"
/usr/bin/time -v -o "$TASK_ROOT/v2_hierarchy.time" \
  "$DEV_CONDA" run -n "$DEV_CONDA_ENV" bash -lc \
  "cd '$DEV_SOURCE_DIR' && PYTHONPATH=src python -m ogprofiler prototype-hierarchy \
    --ssn '$WORKING_DIR/ssn.gml' --out '$V2_ROOT' --weight-attribute NBS" \
  > "$TASK_ROOT/v2_hierarchy.stdout" 2> "$TASK_ROOT/v2_hierarchy.stderr"

"$DEV_CONDA" run -n "$DEV_CONDA_ENV" python \
  "$DEV_SOURCE_DIR/benchmarks/compare_hierarchies.py" \
  --v1-normalized "$TASK_ROOT/v1_hierarchy_summary/normalized_hierarchy.json" \
  --v1-hierarchy-metrics "$V1_HIERARCHY_ROOT/metrics.json" \
  --v2-run "$V2_ROOT" \
  --sequence-ids "$WORKING_DIR/SequenceIDs.txt" \
  --out "$COMPARE_ROOT"

echo "==> [$DATASET] baseline and comparison complete"
