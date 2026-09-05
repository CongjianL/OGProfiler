#!/usr/bin/env bash
#SBATCH --job-name=b1-ogp1-r1
#SBATCH --cpus-per-task=32
#SBATCH --mem=128G
#SBATCH --time=72:00:00
set -euo pipefail
SRC=${DEV_SOURCE_DIR:?}/OGProfiler2_benchmark
ROOT=/home/mselab/licj/project_data/OGProfiler/OGProfiler2_benchmark
OB=/home/mselab/licj/project_data/OGProfiler/Open_Orthobench-master
DATA=$OB/BENCHMARKS/Input
MM=/home/mselab/licj/micromamba/bin/micromamba
ENV=ogprofiler-v1-first
OUT=$ROOT/03_runs/B1_orthobench/OGProfiler1First/rep1_retry1
GROUPS_TSV=$ROOT/04_standardized/groups/OGProfiler1First/OpenOrthobench/groups_rep1.tsv
METRICS=$ROOT/05_metrics/orthobench
EXPECTED=380c85d9548c607f5df6daf656d74b21597fe215ef78d4f8ff296776bb1fdd07
V1_COMMIT=729675fee0dabc02f63ae5c2cda7571748751deb
V1_SOURCE_SHA256=65ec43d269b410956fadc4f8215a3e0a0772eac8bcc1060ab011b090823cca41
[[ ! -e "$OUT" ]] || { echo "Refusing existing run directory: $OUT" >&2; exit 3; }
[[ ! -e "$GROUPS_TSV" ]] || { echo "Refusing existing standardized output: $GROUPS_TSV" >&2; exit 3; }
mkdir -p "$OUT"
"$MM" run -n ogprofiler python "$SRC/scripts/audit/audit_open_orthobench_input.py" --input "$DATA" --refogs "$OB/BENCHMARKS/RefOGs" --out "$OUT/input_audit.tsv"
OBSERVED=$(cat "$OUT/OPEN_ORTHOBENCH_INPUT_DIGEST"); [[ "$OBSERVED" == "$EXPECTED" ]] || { echo "Dataset digest mismatch" >&2; exit 4; }
OBSERVED_SOURCE=$(sha256sum "$SRC/00_env/ogprofiler_v1_first/OGProfiler.py" | awk '{print $1}')
[[ "$OBSERVED_SOURCE" == "$V1_SOURCE_SHA256" ]] || { echo "OGProfiler v1 source digest mismatch" >&2; exit 5; }
"$MM" run -n ogprofiler python "$SRC/scripts/preprocess/prepare_ogprofiler_v1_input.py" --input "$DATA" --out "$OUT/input" --map "$OUT/id_map.tsv"
"$MM" list -n "$ENV" --explicit > "$OUT/environment.explicit.txt"
cat > "$OUT/RUN_MANIFEST.tsv" <<EOF
run_id	B1_orthobench_OGProfiler1First_rep1_retry1
retry_of_slurm_job	1404591
retry_reason	DIAMOND protein-type heuristic rejected two valid Open Orthobench proteomes during makedb
tool	OGProfiler1First
definition	first historical repository version
version_commit	$V1_COMMIT
version_source_sha256	$OBSERVED_SOURCE
dataset_digest	$OBSERVED
benchmark_git_commit	${DEV_GIT_COMMIT:?}
benchmark_source_sha256	${DEV_SOURCE_HASH:?}
slurm_job_id	${SLURM_JOB_ID:?}
threads	32
random_seed	unavailable_in_v1_cli
diamond_makedb_compatibility	--ignore-warnings for input-type heuristic only
primary_output	input/WorkingDirectory/OGFile_coalescence_SameGenome.txt
input_adapter	bijective OGPV1_NN|ORIGINAL_ID; restored before scoring
EOF
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 BENCHMARK_THREADS=32
export OGPROFILER_V1_DIAMOND_REAL=/home/mselab/licj/micromamba/envs/ogprofiler-v1-first/bin/diamond
export PATH="$SRC/02_configs/ogprofiler_v1_first/bin:$PATH"
cd "$OUT"
"$SRC/workflows/run_timed.sh" OGProfiler1First OpenOrthobench 1 "$OUT/raw_tool" -- \
 "$MM" run -n "$ENV" env PATH="$PATH" OGPROFILER_V1_DIAMOND_REAL="$OGPROFILER_V1_DIAMOND_REAL" \
 python "$SRC/00_env/ogprofiler_v1_first/OGProfiler.py" \
 -i "$OUT/input" -o "$OUT/ogprofiler_v1_first.gml" -s diamond -e 1e-5 -t 32 -d lrb -w NBS -m rber -a 32 -g 1.0 --so 0
PRIMARY="$OUT/input/WorkingDirectory/OGFile_coalescence_SameGenome.txt"
[[ -s "$PRIMARY" ]]
"$MM" run -n ogprofiler python "$SRC/scripts/converters/convert_competitor_groups.py" --tool ogprofiler_v1 --input "$PRIMARY" --id-map "$OUT/id_map.tsv" --fasta "$DATA" --out "$GROUPS_TSV"
"$MM" run -n ogprofiler python "$SRC/scripts/converters/validate_groups.py" "$GROUPS_TSV" --fasta "$DATA" > "$OUT/validation.log"
awk -F '\t' 'NR>1 {sub(/\r$/, "", $2); a[$1]=a[$1] " " $2} END {for (x in a) print x ":" a[x]}' "$GROUPS_TSV" | LC_ALL=C sort > "$OUT/prediction.txt"
"$SRC/workflows/run_timed.sh" OpenOrthobench OGProfiler1First 1 "$OUT/raw_official" -- "$MM" run -n ogprofiler python "$OB/BENCHMARKS/benchmark.py" "$OUT/prediction.txt"
mkdir -p "$METRICS/by_method/OGProfiler1First"
"$MM" run -n ogprofiler python "$SRC/scripts/metrics/parse_open_orthobench_official.py" "$OUT/raw_official/stdout.log" --out "$METRICS/by_method/OGProfiler1First/official.tsv"
"$MM" run -n ogprofiler python "$SRC/scripts/metrics/orthobench_metrics.py" --refogs "$OB/BENCHMARKS/RefOGs" --groups "$GROUPS_TSV" --out "$METRICS/by_method/OGProfiler1First/refog_metrics.tsv" --summary "$METRICS/by_method/OGProfiler1First/extended_summary_raw.tsv"
"$MM" run -n ogprofiler python "$SRC/scripts/diagnostics/orthobench_pairwise_attribution.py" --groups "$GROUPS_TSV" --refogs "$OB/BENCHMARKS/RefOGs" --outdir "$METRICS/by_method/OGProfiler1First"
