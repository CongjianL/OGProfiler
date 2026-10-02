#!/usr/bin/env bash
# P5: reuse verified SSN, construct hierarchy once with its existing run.yaml,
# freeze inputs, then compare actual V1 functions and V2 OG artifacts.
# Resource envelope follows verify_ssn_repair.sh; no array / scientific overrides.
#SBATCH --job-name=ogp-p5-fixed-hierarchy
#SBATCH --partition=cu
#SBATCH --account=mselab
#SBATCH --qos=normal
#SBATCH --cpus-per-task=56
#SBATCH --mem=128G
#SBATCH --time=24:00:00
set -euo pipefail
: "${DEV_SOURCE_DIR:?Missing immutable source}"
: "${DEV_RUN_DIR:?Missing independent run directory}"
: "${DEV_CONDA:?Missing configured environment manager}"
: "${DEV_CONDA_ENV:?Missing configured environment}"
ORIGIN=${1:?Usage: p5_fixed_hierarchy.sh VERIFIED_SSN_RUN}
[[ -f "$ORIGIN/../verification.json" && -f "$ORIGIN/edges/edge-manifest.json" ]]
if [[ ${OGP_P5_ENV_READY:-0} != 1 ]]; then
    exec "$DEV_CONDA" run -n "$DEV_CONDA_ENV" env OGP_P5_ENV_READY=1 bash "$0" "$@"
fi
export PYTHONPATH="$DEV_SOURCE_DIR/src:$DEV_SOURCE_DIR"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
RUN="$DEV_RUN_DIR/fixed-run"
[[ ! -e "$RUN" ]]
mkdir "$RUN"
python - "$ORIGIN" "$RUN" <<'PY'
import json, shutil, sys
from pathlib import Path
from ogprofiler.core.manifest import sha256_file
origin, run = map(Path, sys.argv[1:])
verified = json.loads((origin.parent/'verification.json').read_text())
assert verified['passed'] and verified['job_id'] == '1410705'
manifest = json.loads((origin/'edges/edge-manifest.json').read_text())
assert manifest['algorithm_version'] == 'of315-nbs-lrb-directional-mean-v4'
assert manifest['parameters']['symmetrization'] == 'mean'
for name in ('input', 'edges'):
    shutil.copytree(origin/name, run/name)
shutil.copy2(origin/'run.yaml', run/'run.yaml')
for name in ('components', 'hierarchy', 'evolution', 'results', 'search'):
    (run/name).mkdir()
for name, digest in manifest['output_checksums'].items():
    assert sha256_file(run/'edges'/name) == digest == sha256_file(origin/'edges'/name)
for name, digest in manifest['input_checksums'].items():
    if name == 'proteins.parquet':
        assert sha256_file(run/'input'/name) == digest
provenance = dict(origin=str(origin), ssn_verification=verified,
                  edge_manifest=manifest, run_yaml_sha256=sha256_file(origin/'run.yaml'),
                  hierarchy_scientific_overrides=[], search_or_ssn_recomputed=False)
(run.parent/'fixed-origin.json').write_text(json.dumps(provenance, indent=2))
PY
"$DEV_CONDA" list -n "$DEV_CONDA_ENV" > "$DEV_RUN_DIR/environment.txt"
for stage in components hierarchy-all annotate-network; do
    echo "==> P5 $stage"
    /usr/bin/time -v -o "$DEV_RUN_DIR/$stage.time" \
        python -m ogprofiler "$stage" --run "$RUN" \
        > "$DEV_RUN_DIR/$stage.stdout" 2> "$DEV_RUN_DIR/$stage.stderr"
done
echo "==> P5 frozen hierarchy differential and persisted/exported artifact audit"
/usr/bin/time -v -o "$DEV_RUN_DIR/p5.time" \
    python -m benchmarks.og_extraction.fixed_hierarchy \
    --run "$RUN" --out "$DEV_RUN_DIR/p5-audit"
