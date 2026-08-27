#!/usr/bin/env bash
#SBATCH --job-name=ogp2-b1-stage1c
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=04:00:00
set -euo pipefail

# Compute nodes start with a minimal PATH.  Re-enter the frozen project environment.
if [[ "${OGP2_ENV_READY:-0}" != 1 ]]; then
  exec "${DEV_CONDA:?}" run -n "${DEV_CONDA_ENV:?}" env OGP2_ENV_READY=1 bash "$0"
fi

SOURCE_ROOT=${DEV_SOURCE_DIR:?Submitted immutable source snapshot is required}
ROOT=/home/mselab/licj/project_data/OGProfiler/OGProfiler2_benchmark
DATA_ROOT=/home/mselab/licj/project_data/OGProfiler/Open_Orthobench-master
RUN="$ROOT/03_runs/B1_orthobench/OGProfiler2/seed42_rep2/run"
GROUPS="$ROOT/04_standardized/groups/OGProfiler2/OpenOrthobench/groups_seed42_rep2.tsv"
REFOGS="$DATA_ROOT/BENCHMARKS/RefOGs"
METRICS="$ROOT/05_metrics/orthobench"
DIAG_RUN="$ROOT/03_runs/B1_orthobench/OGProfiler2/seed42_rep2/stage1c_diagnostics"
REF_METRICS="$METRICS/OGProfiler2_seed42_refog_metrics.tsv"

for attempt in 1 2 3 4 5; do
  missing=0
  [[ -d "$RUN" ]] || { echo "Missing directory: $RUN (attempt $attempt)" >&2; missing=1; }
  [[ -f "$GROUPS" ]] || { echo "Missing file: $GROUPS (attempt $attempt)" >&2; missing=1; }
  [[ -d "$REFOGS" ]] || { echo "Missing directory: $REFOGS (attempt $attempt)" >&2; missing=1; }
  [[ -f "$REF_METRICS" ]] || { echo "Missing file: $REF_METRICS (attempt $attempt)" >&2; missing=1; }
  [[ "$missing" == 0 ]] && break
  sleep 5
done
[[ "$missing" == 0 ]] || exit 3
[[ ! -e "$DIAG_RUN/STAGE1C_COMPLETE" ]] || {
  echo "Stage 1C diagnostics already complete: $DIAG_RUN" >&2
  exit 4
}
mkdir -p "$DIAG_RUN" "$METRICS"

cat > "$DIAG_RUN/STAGE1C_RUN_MANIFEST.tsv" <<EOF
diagnostic_run_id	B1_orthobench_OGProfiler2_seed42_rep2_stage1c
source_run_id	B1_orthobench_OGProfiler2_seed42_rep2
algorithm_commit	a7838963de4882d867203e60ce641b5689bde45e
formal_benchmark_commit	1375dc1b7a4b053113432db69203f5f7aa9828cd
postprocessing_commit	58531c1c3fb7df097f1ab600d3b00cd85b729eef
diagnostics_commit	${DEV_GIT_COMMIT:?}
diagnostics_source_sha256	${DEV_SOURCE_HASH:?}
slurm_job_id	${SLURM_JOB_ID:?}
started_utc	$(date -u +%Y-%m-%dT%H:%M:%SZ)
mode	read-only attribution of existing artifacts; no prediction rerun
EOF

export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
python "$SOURCE_ROOT/OGProfiler2_benchmark/scripts/diagnostics/orthobench_pairwise_attribution.py" \
  --groups "$GROUPS" --refogs "$REFOGS" --outdir "$METRICS"

/usr/bin/time -v -o "$DIAG_RUN/stage_attribution.time.txt" \
python "$SOURCE_ROOT/OGProfiler2_benchmark/scripts/diagnostics/refog_stage_attribution.py" \
  --run "$RUN" \
  --refogs "$REFOGS" \
  --metrics "$REF_METRICS" \
  --burden "$METRICS/predicted_family_pairwise_burden.tsv" \
  --outdir "$METRICS"

printf 'completed_utc\t%s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)" >> "$DIAG_RUN/STAGE1C_RUN_MANIFEST.tsv"
touch "$DIAG_RUN/STAGE1C_COMPLETE"
