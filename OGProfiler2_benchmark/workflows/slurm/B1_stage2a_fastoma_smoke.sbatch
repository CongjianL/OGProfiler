#!/usr/bin/env bash
#SBATCH --job-name=b1-fastoma-smoke
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --time=24:00:00
set -euo pipefail
SRC=${DEV_SOURCE_DIR:?}/OGProfiler2_benchmark
ROOT=/home/mselab/licj/project_data/OGProfiler/OGProfiler2_benchmark
SMOKE=$ROOT/01_data/processed/smoke_proteomes
MM=/home/mselab/licj/micromamba/bin/micromamba
ENVROOT=$ROOT/00_env/fastoma-0.5.1
OUT=$ROOT/03_runs/B1_orthobench/competitor_smoke_stage2a/fastoma-0.5.1_slurm1
[[ ! -e "$OUT" ]] || { echo "Refusing existing smoke directory: $OUT" >&2; exit 3; }
mkdir -p "$OUT/input/proteome" "$ENVROOT"
cp "$SMOKE"/*.fa "$OUT/input/proteome/"
cat > "$OUT/input/species_tree.nwk" <<'EOF'
(Caenorhabditis_elegans.WBcel235.pep.all,Canis_familiaris.CanFam3.1.pep.all)Bilateria;
EOF
DB=$ENVROOT/LUCA.h5
if [[ ! -s "$DB" ]]; then
 curl -L --fail --retry 5 --continue-at - -o "$DB.part" https://omabrowser.org/All/LUCA.h5
 mv "$DB.part" "$DB"
fi
sha256sum "$DB" > "$ENVROOT/LUCA.h5.sha256"
cat > "$OUT/nextflow-extra.config" <<EOF
conda.useMamba = true
conda.cacheDir = '$ENVROOT/conda-cache'
process.cpus = 8
process.memory = '40 GB'
EOF
export NXF_HOME=$ENVROOT/nxf-home NXF_WORK=$OUT/work NXF_CONDA_CACHEDIR=$ENVROOT/conda-cache BENCHMARK_THREADS=8
"$SRC/workflows/run_timed.sh" FastOMA OpenOrthobenchSmoke 1 "$OUT/timing" -- \
 "$MM" run -n fastoma-0.5.1-launcher nextflow run dessimozlab/FastOMA -r v0.5.1 \
 -profile conda -c "$OUT/nextflow-extra.config" --input "$OUT/input" --species_tree "$OUT/input/species_tree.nwk" \
 --omamer_db "$DB" --output_folder "$OUT/result"
PRIMARY=$OUT/result/RootHOGs.tsv
[[ -f "$PRIMARY" ]] || { echo "RootHOGs.tsv missing" >&2; exit 4; }
"$MM" run -n ogprofiler python "$SRC/scripts/converters/convert_competitor_groups.py" \
 --tool fastoma --input "$PRIMARY" --fasta "$SMOKE" --out "$OUT/groups.tsv"
"$MM" run -n ogprofiler python "$SRC/scripts/converters/validate_groups.py" "$OUT/groups.tsv" --fasta "$SMOKE"
printf 'primary_output\t%s\nstatus\tPASS\n' "$PRIMARY" > "$OUT/SMOKE_STATUS.tsv"
