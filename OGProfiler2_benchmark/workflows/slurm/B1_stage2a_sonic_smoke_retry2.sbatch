#!/usr/bin/env bash
#SBATCH --job-name=b1-sonic-smoke-r2
#SBATCH --cpus-per-task=32
#SBATCH --mem=128G
#SBATCH --time=24:00:00
set -euo pipefail
SRC=${DEV_SOURCE_DIR:?}/OGProfiler2_benchmark
ROOT=/home/mselab/licj/project_data/OGProfiler/OGProfiler2_benchmark
SMOKE=$ROOT/01_data/processed/smoke_proteomes
MM=/home/mselab/licj/micromamba/bin/micromamba
OUT=$ROOT/03_runs/B1_orthobench/competitor_smoke_stage2a/sonicparanoid-2.0.9_slurm3
[[ ! -e "$OUT" ]] || { echo "Refusing existing smoke retry directory: $OUT" >&2; exit 3; }
mkdir -p "$OUT/input"
cp "$SMOKE"/*.fa "$OUT/input/"
cat > "$OUT/RETRY_MANIFEST.tsv" <<EOF
retry_of_slurm_job_id\t1404276
failure\tPfam_index_OOM_plus_non_FASTA_files_in_shared_smoke_directory
resource_change\t32_CPU_128_GB_for_one-time_profile_index_and_smoke
input_change\trun-local_copy_of_the_same_two_FASTA_files_only
slurm_job_id\t${SLURM_JOB_ID:?}
EOF
export BENCHMARK_THREADS=32
"$SRC/workflows/run_timed.sh" SonicParanoid2 OpenOrthobenchSmoke 1 "$OUT/timing" -- \
 "$MM" run -n sonicparanoid-2.0.9 sonicparanoid -i "$OUT/input" -o "$OUT/result" -p smoke -t 32 -m default
PRIMARY=$(find "$OUT/result/runs/smoke" -type f -name 'ortholog_groups.tsv' -print -quit)
[[ -f "$PRIMARY" ]] || { echo "Observed multi-species primary output not found" >&2; find "$OUT/result" -type f; exit 4; }
"$MM" run -n ogprofiler python "$SRC/scripts/converters/convert_competitor_groups.py" --tool sonicparanoid --input "$PRIMARY" --fasta "$SMOKE" --out "$OUT/groups.tsv"
"$MM" run -n ogprofiler python "$SRC/scripts/converters/validate_groups.py" "$OUT/groups.tsv" --fasta "$SMOKE"
printf 'primary_output\t%s\nstatus\tPASS\n' "$PRIMARY" > "$OUT/SMOKE_STATUS.tsv"
