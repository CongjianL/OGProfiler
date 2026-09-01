#!/usr/bin/env bash
#SBATCH --job-name=b1-fastoma-smoke-r2
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --time=24:00:00
set -euo pipefail

SRC=${DEV_SOURCE_DIR:?}/OGProfiler2_benchmark
ROOT=/home/mselab/licj/project_data/OGProfiler/OGProfiler2_benchmark
OB=/home/mselab/licj/project_data/OGProfiler/Open_Orthobench-master
DATA=$OB/BENCHMARKS/Input
SMOKE=$ROOT/01_data/processed/smoke_proteomes
MM=/home/mselab/licj/micromamba/bin/micromamba
ENVROOT=$ROOT/00_env/fastoma-0.5.1
PIPELINE_ENV=$ENVROOT/pipeline-env
WORKFLOW=$ENVROOT/workflow-v0.5.1
DB=$ENVROOT/LUCA.h5
OUT=$ROOT/03_runs/B1_orthobench/competitor_smoke_stage2a/fastoma-0.5.1_slurm3
EXPECTED_DATASET=380c85d9548c607f5df6daf656d74b21597fe215ef78d4f8ff296776bb1fdd07
EXPECTED_DB_SIZE=9918199250
EXPECTED_DB_MD5=87a322a0a6381a8de0b539d62e277054

[[ ! -e "$OUT" ]] || { echo "Refusing existing smoke retry directory: $OUT" >&2; exit 3; }
[[ -f "$DB" && $(stat -c %s "$DB") -eq $EXPECTED_DB_SIZE ]] || { echo "Verified LUCA.h5 is absent" >&2; exit 4; }
echo "$EXPECTED_DB_MD5  $DB" | md5sum -c -
sha256sum -c "$DB.sha256"
[[ -x "$PIPELINE_ENV/bin/nextflow" && -f "$WORKFLOW/FastOMA.nf" ]]
"$PIPELINE_ENV/bin/python" -c 'import FastOMA; assert FastOMA.__version__ == "0.5.1"'

mkdir -p "$OUT/input/proteome"
"$MM" run -n ogprofiler python "$SRC/scripts/audit/audit_open_orthobench_input.py" \
 --input "$DATA" --refogs "$OB/BENCHMARKS/RefOGs" --out "$OUT/frozen_input_audit.tsv"
OBSERVED=$(cat "$OUT/OPEN_ORTHOBENCH_INPUT_DIGEST")
[[ "$OBSERVED" == "$EXPECTED_DATASET" ]] || { echo "Dataset digest mismatch" >&2; exit 5; }
cp "$SMOKE"/*.fa "$OUT/input/proteome/"
FLY=Drosophila_melanogaster.BDGP6.28.pep.all.fa
awk '/^>/{n++} n<=10 {print}' "$DATA/$FLY" > "$OUT/input/proteome/$FLY"
[[ $(grep -c '^>' "$OUT/input/proteome/$FLY") -eq 10 ]]
[[ $(find "$OUT/input/proteome" -maxdepth 1 -type f -name '*.fa' | wc -l) -eq 3 ]]
cat > "$OUT/input/species_tree.nwk" <<'EOF'
((Caenorhabditis_elegans.WBcel235.pep.all,Drosophila_melanogaster.BDGP6.28.pep.all)Protostomia,Canis_familiaris.CanFam3.1.pep.all)Bilateria;
EOF
cat > "$OUT/nextflow-extra.config" <<EOF
conda.enabled = false
docker.enabled = false
singularity.enabled = false
apptainer.enabled = false
process.executor = 'local'
process.resourceLimits = [cpus: 8, memory: 40.GB, time: 23.h]
EOF
cat > "$OUT/RETRY_MANIFEST.tsv" <<EOF
retry_of_slurm_job_id\t1404529
failure\tFastOMA_0.5.1_requires_more_than_two_proteome_files
stage1b_smoke_species_retained\tCaenorhabditis_elegans,Canis_familiaris
added_species\tDrosophila_melanogaster
added_species_selection\tfirst_10_FASTA_records_in_frozen_source_order
added_species_selection_independent_of\tRefOG_labels_and_scores
frozen_dataset_digest\t$OBSERVED
database_size_bytes\t$EXPECTED_DB_SIZE
database_md5\t$EXPECTED_DB_MD5
workflow_tag\tv0.5.1
workflow_commit\t553f860bfd0855251b7e8bf731d226e2b85747b7
slurm_job_id\t${SLURM_JOB_ID:?}
git_commit\t${DEV_GIT_COMMIT:?}
source_sha256\t${DEV_SOURCE_HASH:?}
EOF
sha256sum "$OUT"/input/proteome/*.fa "$OUT/input/species_tree.nwk" > "$OUT/INPUT_SHA256SUMS"
sha256sum "$ENVROOT/pipeline-conda-explicit.txt" "$ENVROOT/pipeline-pip-freeze.txt" \
 "$ENVROOT/workflow-provenance.tsv" > "$OUT/ENVIRONMENT_SHA256SUMS"

export NXF_HOME=$ENVROOT/nxf-home
export NXF_WORK=$OUT/work
export NXF_OFFLINE=true
export PATH=$PIPELINE_ENV/bin:$PATH
export BENCHMARK_THREADS=8
"$SRC/workflows/run_timed.sh" FastOMA OpenOrthobenchSmoke 1 "$OUT/timing" -- \
 "$PIPELINE_ENV/bin/nextflow" run "$WORKFLOW" -c "$OUT/nextflow-extra.config" \
 --input "$OUT/input" --species_tree "$OUT/input/species_tree.nwk" \
 --omamer_db "$DB" --output_folder "$OUT/result"

PRIMARY=$OUT/result/RootHOGs.tsv
[[ -f "$PRIMARY" ]] || { echo "Frozen FastOMA RootHOGs.tsv primary output is absent" >&2; exit 6; }
"$MM" run -n ogprofiler python "$SRC/scripts/converters/convert_competitor_groups.py" \
 --tool fastoma --input "$PRIMARY" --fasta "$OUT/input/proteome" --out "$OUT/groups.tsv"
"$MM" run -n ogprofiler python "$SRC/scripts/converters/validate_groups.py" \
 "$OUT/groups.tsv" --fasta "$OUT/input/proteome" > "$OUT/validation.log"
printf 'primary_output\t%s\nstatus\tPASS\n' "$PRIMARY" > "$OUT/SMOKE_STATUS.tsv"
