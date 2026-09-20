#!/usr/bin/env bash
#SBATCH --job-name=ogp2-ob-v2mod
#SBATCH --cpus-per-task=32
#SBATCH --mem=128G
#SBATCH --time=72:00:00
#SBATCH --output=/home/mselab/licj/project_data/OGProfiler/OGProfiler2_benchmark/03_runs/B1_orthobench/OGProfiler2/v2_modified_rep1/slurm-%j.out
#SBATCH --error=/home/mselab/licj/project_data/OGProfiler/OGProfiler2_benchmark/03_runs/B1_orthobench/OGProfiler2/v2_modified_rep1/slurm-%j.err
set -euo pipefail
# Slurm compute nodes start with a minimal PATH; re-exec in the configured immutable environment.
if [[ "${OGP2_ENV_READY:-0}" != 1 ]]; then
  exec "${DEV_CONDA:?}" run -n "${DEV_CONDA_ENV:?}" env OGP2_ENV_READY=1 bash "$0"
fi
RUN_ID=B1_orthobench_OGProfiler2_v2modified_rep1
ALGORITHM_COMMIT=32a47647760f97ad9097b094f1e507de08f35683
DATA=/home/mselab/licj/project_data/OGProfiler/Open_Orthobench-master/BENCHMARKS/Input
ROOT=/home/mselab/licj/project_data/OGProfiler/OGProfiler2_benchmark
OUT="$ROOT/03_runs/B1_orthobench/OGProfiler2/v2_modified_rep1"
SOURCE_ROOT=${DEV_SOURCE_DIR:?Submitted immutable source snapshot is required}
[[ ! -e "$OUT" ]] || { echo "Refusing existing run directory: $OUT" >&2; exit 3; }
mkdir -p "$OUT"
python "$SOURCE_ROOT/OGProfiler2_benchmark/scripts/audit/audit_open_orthobench_input.py" --input "$DATA" --refogs /home/mselab/licj/project_data/OGProfiler/Open_Orthobench-master/BENCHMARKS/RefOGs --out "$ROOT/manifest/open_orthobench_run_input.tsv"
DATASET_DIGEST=$(cat "$ROOT/manifest/OPEN_ORTHOBENCH_INPUT_DIGEST")
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
export BENCHMARK_THREADS=32 BENCHMARK_SEED=42 BENCHMARK_GIT_COMMIT=${DEV_GIT_COMMIT:?}
command=(env "PYTHONPATH=$SOURCE_ROOT/src" python -m ogprofiler run --proteomes "$DATA" --out "$OUT/run" --set search.threads=32 --set runtime.workers=32 --set hierarchy.seed=42)
printf 'run_id\t%s\ntool\tOGProfiler2\ndataset\tOpenOrthobench\nseed\t42\nthreads\t32\nworkers\t32\nalgorithm_commit\t%s\nbenchmark_commit\t%s\nstart_time_utc\t%s\nhostname\t%s\nslurm_job_id\t%s\ndataset_digest\t%s\ncommand\t' "$RUN_ID" "$ALGORITHM_COMMIT" "$DEV_GIT_COMMIT" "$(date -u +%Y-%m-%dT%H:%M:%SZ)" "$(hostname)" "$SLURM_JOB_ID" "$DATASET_DIGEST" > "$OUT/RUN_MANIFEST.tsv"
printf '%q ' "${command[@]}" >> "$OUT/RUN_MANIFEST.tsv"; printf '\n' >> "$OUT/RUN_MANIFEST.tsv"
"$SOURCE_ROOT/OGProfiler2_benchmark/workflows/run_timed.sh" OGProfiler2 OpenOrthobench 1 "$OUT/raw_ogprofiler" -- "${command[@]}"
