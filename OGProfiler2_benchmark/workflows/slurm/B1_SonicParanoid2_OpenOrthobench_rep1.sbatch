#!/usr/bin/env bash
#SBATCH --job-name=b1-sonicparanoid2
#SBATCH --cpus-per-task=32
#SBATCH --mem=128G
#SBATCH --time=72:00:00
set -euo pipefail

SRC=${DEV_SOURCE_DIR:?}/OGProfiler2_benchmark
ROOT=/home/mselab/licj/project_data/OGProfiler/OGProfiler2_benchmark
OB=/home/mselab/licj/project_data/OGProfiler/Open_Orthobench-master
DATA=$OB/BENCHMARKS/Input
MM=/home/mselab/licj/micromamba/bin/micromamba
SONIC_ENV=sonicparanoid-2.0.9-mmseqs13
ENV=/home/mselab/licj/micromamba/envs/$SONIC_ENV
PFAM=$ENV/lib/python3.12/site-packages/sonicparanoid/pfam_files/profile_db/pfama.mmseqs.idx
PREP=$ROOT/03_runs/B1_orthobench/environment_preparation/sonicparanoid-2.0.9_mmseqs13_pfam_index_unsplit5
OUT=$ROOT/03_runs/B1_orthobench/SonicParanoid2/rep1
GROUPS_TSV=$ROOT/04_standardized/groups/SonicParanoid2/OpenOrthobench/groups_rep1.tsv
METRICS=$ROOT/05_metrics/orthobench/by_method/SonicParanoid2
EXPECTED=380c85d9548c607f5df6daf656d74b21597fe215ef78d4f8ff296776bb1fdd07

[[ ! -e "$OUT" ]] || { echo "Refusing existing formal run directory: $OUT" >&2; exit 3; }
[[ ! -e "$GROUPS_TSV" ]] || { echo "Refusing existing standardized output: $GROUPS_TSV" >&2; exit 3; }
[[ ! -e "$METRICS" ]] || { echo "Refusing existing metric output: $METRICS" >&2; exit 3; }
[[ -f "$PREP/PREPARATION_STATUS.tsv" && -f "$PFAM" && -f "$PFAM.index" && -f "$PFAM.dbtype" ]] || {
 echo "Verified SonicParanoid Pfam index is absent" >&2
 exit 4
}
grep -q $'^status\tPASS$' "$PREP/PREPARATION_STATUS.tsv"
sha256sum -c "$PREP/index_files.sha256"

mkdir -p "$OUT/input"
"$MM" run -n ogprofiler python "$SRC/scripts/audit/audit_open_orthobench_input.py" \
 --input "$DATA" --refogs "$OB/BENCHMARKS/RefOGs" --out "$OUT/input_audit.tsv"
OBSERVED=$(cat "$OUT/OPEN_ORTHOBENCH_INPUT_DIGEST")
[[ "$OBSERVED" == "$EXPECTED" ]] || { echo "Dataset digest mismatch" >&2; exit 5; }
cp "$DATA"/*.fa "$OUT/input/"
[[ $(find "$OUT/input" -maxdepth 1 -type f -name '*.fa' | wc -l) -eq 12 ]]
sha256sum "$OUT"/input/*.fa > "$OUT/EXECUTION_INPUT_SHA256SUMS"
"$MM" list -n "$SONIC_ENV" --explicit > "$OUT/conda-explicit.txt"

cat > "$OUT/RUN_MANIFEST.tsv" <<EOF
run_id\tB1_orthobench_SonicParanoid2_rep1
tool\tSonicParanoid2
version\t2.0.9
environment\t$SONIC_ENV
mmseqs_version\t13.45111
dataset_digest\t$OBSERVED
git_commit\t${DEV_GIT_COMMIT:?}
source_sha256\t${DEV_SOURCE_HASH:?}
slurm_job_id\t${SLURM_JOB_ID:?}
environment_preparation_slurm_job_id\t1404430
primary_output\truns/B1_orthobench_SonicParanoid2_rep1/ortholog_groups/ortholog_groups.tsv
tool_mode\tdefault
threads\t32
EOF

export BENCHMARK_THREADS=32
"$SRC/workflows/run_timed.sh" SonicParanoid2 OpenOrthobench 1 "$OUT/raw_tool" -- \
 "$MM" run -n "$SONIC_ENV" sonicparanoid -i "$OUT/input" -o "$OUT/result" \
 -p B1_orthobench_SonicParanoid2_rep1 -t 32 -m default

PRIMARY=$OUT/result/runs/B1_orthobench_SonicParanoid2_rep1/ortholog_groups/ortholog_groups.tsv
[[ -f "$PRIMARY" ]] || { echo "Frozen SonicParanoid primary output is absent" >&2; exit 6; }
"$MM" run -n ogprofiler python "$SRC/scripts/converters/convert_competitor_groups.py" \
 --tool sonicparanoid --input "$PRIMARY" --fasta "$DATA" --out "$GROUPS_TSV"
"$MM" run -n ogprofiler python "$SRC/scripts/converters/validate_groups.py" \
 "$GROUPS_TSV" --fasta "$DATA" > "$OUT/validation.log"
awk -F '\t' 'NR>1 {sub(/\r$/, "", $2); a[$1]=a[$1] " " $2} END {for (x in a) print x ":" a[x]}' \
 "$GROUPS_TSV" | LC_ALL=C sort > "$OUT/prediction.txt"

"$SRC/workflows/run_timed.sh" OpenOrthobench SonicParanoid2 1 "$OUT/raw_official" -- \
 "$MM" run -n ogprofiler python "$OB/BENCHMARKS/benchmark.py" "$OUT/prediction.txt"
mkdir -p "$METRICS"
"$MM" run -n ogprofiler python "$SRC/scripts/metrics/parse_open_orthobench_official.py" \
 "$OUT/raw_official/stdout.log" --out "$METRICS/official.tsv"
"$MM" run -n ogprofiler python "$SRC/scripts/metrics/orthobench_metrics.py" \
 --refogs "$OB/BENCHMARKS/RefOGs" --groups "$GROUPS_TSV" \
 --out "$METRICS/refog_metrics.tsv" --summary "$METRICS/extended_summary_raw.tsv"
"$MM" run -n ogprofiler python "$SRC/scripts/diagnostics/orthobench_pairwise_attribution.py" \
 --groups "$GROUPS_TSV" --refogs "$OB/BENCHMARKS/RefOGs" --outdir "$METRICS"
[[ $(($(wc -l < "$METRICS/refog_metrics.tsv")-1)) -eq 70 ]]
