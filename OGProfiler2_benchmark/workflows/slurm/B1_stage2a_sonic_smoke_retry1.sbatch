#!/usr/bin/env bash
#SBATCH --job-name=b1-sonic-smoke-r1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=24:00:00
set -euo pipefail
SRC=${DEV_SOURCE_DIR:?}/OGProfiler2_benchmark
ROOT=/home/mselab/licj/project_data/OGProfiler/OGProfiler2_benchmark
SMOKE=$ROOT/01_data/processed/smoke_proteomes
MM=/home/mselab/licj/micromamba/bin/micromamba
OUT=$ROOT/03_runs/B1_orthobench/competitor_smoke_stage2a/sonicparanoid-2.0.9_slurm2
[[ ! -e "$OUT" ]] || { echo "Refusing existing smoke retry directory: $OUT" >&2; exit 3; }
mkdir -p "$OUT"
cat > "$OUT/RETRY_MANIFEST.tsv" <<EOF
retry_of_slurm_job_id\t1404264
failure\tfirst-run_Pfam_profile_indexing_exceeded_4h_smoke_limit
resource_change\ttime_04h_to_24h_smoke_only
slurm_job_id\t${SLURM_JOB_ID:?}
EOF
export BENCHMARK_THREADS=8
"$SRC/workflows/run_timed.sh" SonicParanoid2 OpenOrthobenchSmoke 1 "$OUT/timing" -- \
 "$MM" run -n sonicparanoid-2.0.9 sonicparanoid -i "$SMOKE" -o "$OUT/result" -p smoke -t 8 -m default
PRIMARY=$(find "$OUT/result/runs/smoke" -type f -name 'ortholog_groups.tsv' -print -quit)
[[ -f "$PRIMARY" ]] || { echo "Observed multi-species primary output not found" >&2; find "$OUT/result" -type f; exit 4; }
"$MM" run -n ogprofiler python "$SRC/scripts/converters/convert_competitor_groups.py" --tool sonicparanoid --input "$PRIMARY" --fasta "$SMOKE" --out "$OUT/groups.tsv"
"$MM" run -n ogprofiler python "$SRC/scripts/converters/validate_groups.py" "$OUT/groups.tsv" --fasta "$SMOKE"
printf 'primary_output\t%s\nstatus\tPASS\n' "$PRIMARY" > "$OUT/SMOKE_STATUS.tsv"
