#!/usr/bin/env bash
#SBATCH --job-name=b1-sonic-smoke-r3
#SBATCH --cpus-per-task=32
#SBATCH --mem=128G
#SBATCH --time=24:00:00
set -euo pipefail

SRC=${DEV_SOURCE_DIR:?}/OGProfiler2_benchmark
ROOT=/home/mselab/licj/project_data/OGProfiler/OGProfiler2_benchmark
OB=/home/mselab/licj/project_data/OGProfiler/Open_Orthobench-master
DATA=$OB/BENCHMARKS/Input
MM=/home/mselab/licj/micromamba/bin/micromamba
SONIC_ENV=sonicparanoid-2.0.9-mmseqs13
PFAM=/home/mselab/licj/micromamba/envs/$SONIC_ENV/lib/python3.12/site-packages/sonicparanoid/pfam_files/profile_db/pfama.mmseqs.idx
PREP=$ROOT/03_runs/B1_orthobench/environment_preparation/sonicparanoid-2.0.9_mmseqs13_pfam_index_unsplit5
OUT=$ROOT/03_runs/B1_orthobench/competitor_smoke_stage2a/sonicparanoid-2.0.9_slurm4

[[ ! -e "$OUT" ]] || { echo "Refusing existing smoke retry directory: $OUT" >&2; exit 3; }
[[ -f "$PREP/PREPARATION_STATUS.tsv" && -f "$PFAM" && -f "$PFAM.index" && -f "$PFAM.dbtype" ]] || {
 echo "Verified split Pfam index is absent" >&2
 exit 4
}

mkdir -p "$OUT/input"
for name in \
 Caenorhabditis_elegans.WBcel235.pep.all.fa \
 Canis_familiaris.CanFam3.1.pep.all.fa
do
 [[ -f "$DATA/$name" ]] || { echo "Missing frozen Open Orthobench proteome: $name" >&2; exit 5; }
 cp "$DATA/$name" "$OUT/input/$name"
done

cat > "$OUT/RETRY_MANIFEST.tsv" <<EOF
retry_of_slurm_job_id	1404284
failure	unsplit_Pfam_index_OOM_plus_empty_alignment_in_tiny_smoke
environment_preparation_slurm_job_id	1404430
environment_change	verified_Sonic-compatible_index_of_the_identical_bundled_Pfam_database_with_MMseqs2_13.45111
input_change	the_same_two_species_as_retry2_using_their_complete_frozen_Open_Orthobench_FASTA_files
selection_policy	explicit_species_names_independent_of_RefOG_labels
tool_mode	default
threads	32
slurm_job_id	${SLURM_JOB_ID:?}
git_commit	${DEV_GIT_COMMIT:?}
source_sha256	${DEV_SOURCE_HASH:?}
EOF
sha256sum "$OUT"/input/*.fa > "$OUT/INPUT_SHA256SUMS"

export BENCHMARK_THREADS=32
"$SRC/workflows/run_timed.sh" SonicParanoid2 OpenOrthobenchSmoke 1 "$OUT/timing" -- \
 "$MM" run -n "$SONIC_ENV" sonicparanoid -i "$OUT/input" -o "$OUT/result" -p smoke -t 32 -m default

PRIMARY=$(find "$OUT/result/runs/smoke" -type f -name 'ortholog_groups.tsv' -print -quit)
[[ -f "$PRIMARY" ]] || { echo "Observed multi-species primary output not found" >&2; find "$OUT/result" -type f; exit 6; }
"$MM" run -n ogprofiler python "$SRC/scripts/converters/convert_competitor_groups.py" \
 --tool sonicparanoid --input "$PRIMARY" --fasta "$OUT/input" --out "$OUT/groups.tsv"
"$MM" run -n ogprofiler python "$SRC/scripts/converters/validate_groups.py" \
 "$OUT/groups.tsv" --fasta "$OUT/input" > "$OUT/validation.log"
printf 'primary_output\t%s\nstatus\tPASS\n' "$PRIMARY" > "$OUT/SMOKE_STATUS.tsv"
