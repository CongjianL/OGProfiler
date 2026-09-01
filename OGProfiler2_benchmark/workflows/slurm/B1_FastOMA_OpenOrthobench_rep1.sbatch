#!/usr/bin/env bash
#SBATCH --job-name=b1-fastoma
#SBATCH --cpus-per-task=32
#SBATCH --mem=128G
#SBATCH --time=72:00:00
set -euo pipefail

SRC=${DEV_SOURCE_DIR:?}/OGProfiler2_benchmark
ROOT=/home/mselab/licj/project_data/OGProfiler/OGProfiler2_benchmark
OB=/home/mselab/licj/project_data/OGProfiler/Open_Orthobench-master
DATA=$OB/BENCHMARKS/Input
MM=/home/mselab/licj/micromamba/bin/micromamba
ENVROOT=$ROOT/00_env/fastoma-0.5.1
PIPELINE_ENV=$ENVROOT/pipeline-env
WORKFLOW=$ENVROOT/workflow-v0.5.1
DB=$ENVROOT/LUCA.h5
OUT=$ROOT/03_runs/B1_orthobench/FastOMA/rep1
GROUPS_TSV=$ROOT/04_standardized/groups/FastOMA/OpenOrthobench/groups_rep1.tsv
METRICS=$ROOT/05_metrics/orthobench/by_method/FastOMA
TREE=$SRC/02_configs/fastoma/open_orthobench_ncbi_taxonomy_tree.nwk
EXPECTED_DATASET=380c85d9548c607f5df6daf656d74b21597fe215ef78d4f8ff296776bb1fdd07
EXPECTED_DB_SIZE=9918199250
EXPECTED_DB_MD5=87a322a0a6381a8de0b539d62e277054

[[ ! -e "$OUT" ]] || { echo "Refusing existing formal run directory: $OUT" >&2; exit 3; }
[[ ! -e "$GROUPS_TSV" ]] || { echo "Refusing existing standardized output: $GROUPS_TSV" >&2; exit 3; }
[[ ! -e "$METRICS" ]] || { echo "Refusing existing metric output: $METRICS" >&2; exit 3; }
[[ -f "$DB" && $(stat -c %s "$DB") -eq $EXPECTED_DB_SIZE ]] || { echo "Verified LUCA.h5 is absent" >&2; exit 4; }
echo "$EXPECTED_DB_MD5  $DB" | md5sum -c -
sha256sum -c "$DB.sha256"
[[ -x "$PIPELINE_ENV/bin/nextflow" && -f "$WORKFLOW/FastOMA.nf" && -f "$TREE" ]]
"$PIPELINE_ENV/bin/python" -c 'import FastOMA; assert FastOMA.__version__ == "0.5.1"'

mkdir -p "$OUT/input/proteome"
"$MM" run -n ogprofiler python "$SRC/scripts/audit/audit_open_orthobench_input.py" \
 --input "$DATA" --refogs "$OB/BENCHMARKS/RefOGs" --out "$OUT/frozen_input_audit.tsv"
OBSERVED=$(cat "$OUT/OPEN_ORTHOBENCH_INPUT_DIGEST")
[[ "$OBSERVED" == "$EXPECTED_DATASET" ]] || { echo "Dataset digest mismatch" >&2; exit 5; }
cp "$DATA"/*.fa "$OUT/input/proteome/"
cp "$TREE" "$OUT/input/species_tree.nwk"
[[ $(find "$OUT/input/proteome" -maxdepth 1 -type f -name '*.fa' | wc -l) -eq 12 ]]
sha256sum "$OUT"/input/proteome/*.fa "$OUT/input/species_tree.nwk" > "$OUT/EXECUTION_INPUT_SHA256SUMS"
sha256sum "$ENVROOT/pipeline-conda-explicit.txt" "$ENVROOT/pipeline-pip-freeze.txt" \
 "$ENVROOT/workflow-provenance.tsv" > "$OUT/ENVIRONMENT_SHA256SUMS"

cat > "$OUT/nextflow-extra.config" <<EOF
conda.enabled = false
docker.enabled = false
singularity.enabled = false
apptainer.enabled = false
process.executor = 'local'
process.resourceLimits = [cpus: 32, memory: 120.GB, time: 71.h]
EOF
cat > "$OUT/RUN_MANIFEST.tsv" <<EOF
run_id\tB1_orthobench_FastOMA_rep1
tool\tFastOMA
version\t0.5.1
workflow_tag\tv0.5.1
workflow_commit\t553f860bfd0855251b7e8bf731d226e2b85747b7
environment\t$PIPELINE_ENV
dataset_digest\t$OBSERVED
species_tree\t02_configs/fastoma/open_orthobench_ncbi_taxonomy_tree.nwk
species_tree_source\tNCBI_taxonomy_independent_of_RefOG_labels_and_scores
omamer_database\tLUCA.h5
omamer_database_size_bytes\t$EXPECTED_DB_SIZE
omamer_database_md5\t$EXPECTED_DB_MD5
git_commit\t${DEV_GIT_COMMIT:?}
source_sha256\t${DEV_SOURCE_HASH:?}
slurm_job_id\t${SLURM_JOB_ID:?}
primary_output\tRootHOGs.tsv
tool_mode\tdocumented_default
threads\t32
EOF

export NXF_HOME=$ENVROOT/nxf-home
export NXF_WORK=$OUT/work
export NXF_OFFLINE=true
export PATH=$PIPELINE_ENV/bin:$PATH
export BENCHMARK_THREADS=32
"$SRC/workflows/run_timed.sh" FastOMA OpenOrthobench 1 "$OUT/raw_tool" -- \
 "$PIPELINE_ENV/bin/nextflow" run "$WORKFLOW" -c "$OUT/nextflow-extra.config" \
 --input "$OUT/input" --species_tree "$OUT/input/species_tree.nwk" \
 --omamer_db "$DB" --output_folder "$OUT/result"

PRIMARY=$OUT/result/RootHOGs.tsv
[[ -f "$PRIMARY" ]] || { echo "Frozen FastOMA RootHOGs.tsv primary output is absent" >&2; exit 6; }
"$MM" run -n ogprofiler python "$SRC/scripts/converters/convert_competitor_groups.py" \
 --tool fastoma --input "$PRIMARY" --fasta "$OUT/input/proteome" --out "$GROUPS_TSV"
"$MM" run -n ogprofiler python "$SRC/scripts/converters/validate_groups.py" \
 "$GROUPS_TSV" --fasta "$OUT/input/proteome" > "$OUT/validation.log"
awk -F '\t' 'NR>1 {sub(/\r$/, "", $2); a[$1]=a[$1] " " $2} END {for (x in a) print x ":" a[x]}' \
 "$GROUPS_TSV" | LC_ALL=C sort > "$OUT/prediction.txt"

"$SRC/workflows/run_timed.sh" OpenOrthobench FastOMA 1 "$OUT/raw_official" -- \
 "$MM" run -n ogprofiler python "$OB/BENCHMARKS/benchmark.py" "$OUT/prediction.txt"
mkdir -p "$METRICS"
"$MM" run -n ogprofiler python "$SRC/scripts/metrics/parse_open_orthobench_official.py" \
 "$OUT/raw_official/stdout.log" --out "$METRICS/official.tsv"
"$MM" run -n ogprofiler python "$SRC/scripts/metrics/orthobench_metrics.py" \
 --refogs "$OB/BENCHMARKS/RefOGs" --groups "$GROUPS_TSV" \
 --out "$METRICS/refog_metrics.tsv" --summary "$METRICS/extended_summary_raw.tsv"
"$MM" run -n ogprofiler python "$SRC/scripts/diagnostics/orthobench_pairwise_attribution.py" \
 --groups "$GROUPS_TSV" --refogs "$OB/BENCHMARKS/RefOGs" --outdir "$METRICS"
[[ $(($(wc -l < "$METRICS/refog_metrics.tsv") - 1)) -eq 70 ]]
printf 'primary_output\t%s\ngroups\t%s\nstatus\tPASS\n' "$PRIMARY" "$GROUPS_TSV" > "$OUT/FORMAL_STATUS.tsv"
