#!/usr/bin/env bash
# Fresh FASTA -> production edges -> independent OF layer/artifact verification.
# Resources follow the existing S0 search job; all outputs are run-local.
#SBATCH --job-name=ogp-verify-ssn-repair
#SBATCH --partition=cu
#SBATCH --account=mselab
#SBATCH --qos=normal
#SBATCH --cpus-per-task=56
#SBATCH --mem=128G
#SBATCH --time=24:00:00
set -euo pipefail
: "${DEV_SOURCE_DIR:?Missing immutable source snapshot}"
: "${DEV_RUN_DIR:?Missing run directory}"
: "${DEV_CONDA:?Missing configured environment manager}"
: "${DEV_CONDA_ENV:?Missing configured execution environment}"
INPUT_DIR=${1:?Usage: verify_ssn_repair.sh PROTEOMES OF_SOURCE_ROOT}
OF_SOURCE_ROOT=${2:?Provide the actual OF3 source directory}
[[ -d "$INPUT_DIR" && -f "$OF_SOURCE_ROOT/tools/waterfall.py" ]]
if [[ ${OGP_VERIFY_ENV_READY:-0} != 1 ]]; then
    exec "$DEV_CONDA" run -n "$DEV_CONDA_ENV" env OGP_VERIFY_ENV_READY=1 bash "$0" "$@"
fi
export PYTHONPATH="$DEV_SOURCE_DIR/src"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
THREADS=${SLURM_CPUS_PER_TASK:?Missing allocated CPU count}
PROTEOMES="$DEV_RUN_DIR/proteomes"
V2_RUN="$DEV_RUN_DIR/v2-run"
[[ ! -e "$PROTEOMES" && ! -e "$V2_RUN" && ! -e "$DEV_RUN_DIR/ssn-audit" ]]
input_manifest() {
    (cd "$1" && find . -maxdepth 1 -type f -print0 | sort -z | xargs -0 sha256sum)
}
input_manifest "$INPUT_DIR" > "$DEV_RUN_DIR/input-before.sha256"
[[ -s "$DEV_RUN_DIR/input-before.sha256" ]]
mkdir "$PROTEOMES"
cp -a "$INPUT_DIR/." "$PROTEOMES/"
input_manifest "$PROTEOMES" > "$DEV_RUN_DIR/proteomes-snapshot.sha256"
cmp "$DEV_RUN_DIR/input-before.sha256" "$DEV_RUN_DIR/proteomes-snapshot.sha256"
"$DEV_CONDA" list -n "$DEV_CONDA_ENV" > "$DEV_RUN_DIR/environment.txt"
diamond version > "$DEV_RUN_DIR/diamond-version.txt"
printf '%s\n' "$INPUT_DIR" > "$DEV_RUN_DIR/input-origin.txt"
echo "==> Fresh prepare/search/edges: $V2_RUN"
/usr/bin/time -v -o "$DEV_RUN_DIR/v2.time" python -m ogprofiler run \
    --proteomes "$PROTEOMES" --out "$V2_RUN" --until-stage edges \
    --set "search.threads=$THREADS" --set edges.symmetrization=mean \
    > "$DEV_RUN_DIR/v2.stdout" 2> "$DEV_RUN_DIR/v2.stderr"
echo "==> OF directional layers and actual persisted V2 artifact"
/usr/bin/time -v -o "$DEV_RUN_DIR/audit.time" python \
    "$DEV_SOURCE_DIR/benchmarks/compare_fixed_hits_ssn.py" \
    --run "$V2_RUN" --of-source "$OF_SOURCE_ROOT" --out "$DEV_RUN_DIR/ssn-audit" \
    --production-edges "$V2_RUN/edges/retained_edges.parquet" --require-equivalence
input_manifest "$INPUT_DIR" > "$DEV_RUN_DIR/input-after.sha256"
input_manifest "$PROTEOMES" > "$DEV_RUN_DIR/proteomes-after.sha256"
cmp "$DEV_RUN_DIR/input-before.sha256" "$DEV_RUN_DIR/input-after.sha256"
cmp "$DEV_RUN_DIR/proteomes-snapshot.sha256" "$DEV_RUN_DIR/proteomes-after.sha256"
python - <<'PY'
import json, os
from pathlib import Path
root = Path(os.environ['DEV_RUN_DIR'])
report = json.loads((root / 'ssn-audit/report.json').read_text())
assert report['completed'] and report['input_integrity_verified'] and report['validation']['passed']
summary = dict(passed=True, job_id=os.environ['SLURM_JOB_ID'], run_id=os.environ['DEV_RUN_ID'],
               git_commit=os.environ['DEV_GIT_COMMIT'], source_hash=os.environ['DEV_SOURCE_HASH'],
               original_proteomes=(root/'input-origin.txt').read_text().strip(),
               proteomes_snapshot_verified=True, original_proteomes_unchanged=True,
               fixed_hits_layers=report['validation'], production_artifact=report['production_artifact'],
               graph_objects=report['graph_objects'], report=str(root/'ssn-audit/report.json'),
               scope='OF full-run directional SSN semantics and Leiden mean projection; not OG accuracy')
(root/'verification.json').write_text(json.dumps(summary, indent=2))
print(json.dumps(summary, indent=2))
PY
