#!/usr/bin/env bash
#SBATCH --job-name=b1-proteinortho6
#SBATCH --cpus-per-task=32
#SBATCH --mem=128G
#SBATCH --time=72:00:00
set -euo pipefail
SRC=${DEV_SOURCE_DIR:?}/OGProfiler2_benchmark
ROOT=/home/mselab/licj/project_data/OGProfiler/OGProfiler2_benchmark
OB=/home/mselab/licj/project_data/OGProfiler/Open_Orthobench-master
DATA=$OB/BENCHMARKS/Input
MM=/home/mselab/licj/micromamba/bin/micromamba
OUT=$ROOT/03_runs/B1_orthobench/Proteinortho6/rep1
GROUPS_TSV=$ROOT/04_standardized/groups/Proteinortho6/OpenOrthobench/groups_rep1.tsv
METRICS=$ROOT/05_metrics/orthobench
EXPECTED=380c85d9548c607f5df6daf656d74b21597fe215ef78d4f8ff296776bb1fdd07
[[ ! -e "$OUT" ]] || { echo "Refusing existing run directory: $OUT" >&2; exit 3; }
mkdir -p "$OUT/work"
"$MM" run -n ogprofiler python "$SRC/scripts/audit/audit_open_orthobench_input.py" --input "$DATA" --refogs "$OB/BENCHMARKS/RefOGs" --out "$OUT/input_audit.tsv"
OBSERVED=$(cat "$OUT/OPEN_ORTHOBENCH_INPUT_DIGEST"); [[ "$OBSERVED" == "$EXPECTED" ]] || { echo "Dataset digest mismatch" >&2; exit 4; }
cat > "$OUT/RUN_MANIFEST.tsv" <<EOF
run_id\tB1_orthobench_Proteinortho6_rep1
tool\tProteinortho
version\t6.3.6
dataset_digest\t$OBSERVED
git_commit\t${DEV_GIT_COMMIT:?}
source_sha256\t${DEV_SOURCE_HASH:?}
slurm_job_id\t${SLURM_JOB_ID:?}
primary_output\tB1_orthobench_Proteinortho6_rep1.proteinortho.tsv
synteny\tOFF
EOF
export BENCHMARK_THREADS=32
CMD="cd '$OUT/work' && proteinortho6.pl -project=B1_orthobench_Proteinortho6_rep1 -cpus=32 -singles '$DATA'/*.fa"
"$SRC/workflows/run_timed.sh" Proteinortho6 OpenOrthobench 1 "$OUT/raw_tool" -- "$MM" run -n proteinortho-6.3.6 bash -lc "$CMD"
PRIMARY=$OUT/work/B1_orthobench_Proteinortho6_rep1.proteinortho.tsv; [[ -f "$PRIMARY" ]]
"$MM" run -n ogprofiler python "$SRC/scripts/converters/convert_competitor_groups.py" --tool proteinortho --input "$PRIMARY" --fasta "$DATA" --out "$GROUPS_TSV"
"$MM" run -n ogprofiler python "$SRC/scripts/converters/validate_groups.py" "$GROUPS_TSV" --fasta "$DATA" > "$OUT/validation.log"
awk -F '\t' 'NR>1 {sub(/\r$/, "", $2); a[$1]=a[$1] " " $2} END {for (x in a) print x ":" a[x]}' "$GROUPS_TSV" | LC_ALL=C sort > "$OUT/prediction.txt"
"$SRC/workflows/run_timed.sh" OpenOrthobench Proteinortho6 1 "$OUT/raw_official" -- "$MM" run -n ogprofiler python "$OB/BENCHMARKS/benchmark.py" "$OUT/prediction.txt"
mkdir -p "$METRICS/by_method/Proteinortho6"
"$MM" run -n ogprofiler python "$SRC/scripts/metrics/parse_open_orthobench_official.py" "$OUT/raw_official/stdout.log" --out "$METRICS/by_method/Proteinortho6/official.tsv"
"$MM" run -n ogprofiler python "$SRC/scripts/metrics/orthobench_metrics.py" --refogs "$OB/BENCHMARKS/RefOGs" --groups "$GROUPS_TSV" --out "$METRICS/by_method/Proteinortho6/refog_metrics.tsv" --summary "$METRICS/by_method/Proteinortho6/extended_summary_raw.tsv"
"$MM" run -n ogprofiler python "$SRC/scripts/diagnostics/orthobench_pairwise_attribution.py" --groups "$GROUPS_TSV" --refogs "$OB/BENCHMARKS/RefOGs" --outdir "$METRICS/by_method/Proteinortho6"
