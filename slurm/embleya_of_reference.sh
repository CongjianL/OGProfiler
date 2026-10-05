#!/usr/bin/env bash
# Frozen Embleya mean SSN -> current default and explicit soft/depth42 candidate.
# Resource envelope unchanged from P5; sequential variants, no array/search.
#SBATCH --job-name=ogp-embleya-of-reference
#SBATCH --partition=cu
#SBATCH --account=mselab
#SBATCH --qos=normal
#SBATCH --cpus-per-task=56
#SBATCH --mem=128G
#SBATCH --time=24:00:00
set -euo pipefail
: "${DEV_SOURCE_DIR:?}"
: "${DEV_RUN_DIR:?}"
: "${DEV_CONDA:?}"
: "${DEV_CONDA_ENV:?}"
ORIGIN=${1:?Usage: embleya_of_reference.sh P5_FIXED_RUN OF_ORTHOGROUPS_DIR}
REFERENCE=${2:?Missing OF reference directory}
if [[ ${OGP_EMBLEYA_READY:-0} != 1 ]]; then
    exec "$DEV_CONDA" run -n "$DEV_CONDA_ENV" env OGP_EMBLEYA_READY=1 bash "$0" "$@"
fi
export PYTHONPATH="$DEV_SOURCE_DIR/src:$DEV_SOURCE_DIR"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
python -m pytest -q "$DEV_SOURCE_DIR/tests/test_embleya_reference.py" \
    > "$DEV_RUN_DIR/evaluator-smoke.txt" 2>&1
python - "$ORIGIN" "$REFERENCE" <<'PY'
import json, os, shutil
from pathlib import Path
import yaml
from ogprofiler.config import DEFAULT_CONFIG, validate_config
from ogprofiler.core.manifest import sha256_file
from benchmarks.og_extraction.embleya_reference import evaluate
origin, reference = map(Path, __import__('sys').argv[1:])
root = Path(os.environ['DEV_RUN_DIR'])
edge = json.loads((origin/'edges/edge-manifest.json').read_text())
assert edge['algorithm_version'] == 'of315-nbs-lrb-directional-mean-v4'
assert edge['parameters']['symmetrization'] == 'mean'
assert json.loads((origin.parent/'p5-audit/report.json').read_text())['passed']
for name, h in edge['output_checksums'].items():
    assert sha256_file(origin/'edges'/name) == h
ref = root/'reference'
ref.mkdir()
for name in ('Orthogroups.tsv', 'Orthogroups_UnassignedGenes.tsv'):
    shutil.copy2(reference/name, ref/name)
    assert sha256_file(reference/name) == sha256_file(ref/name)
evaluate(origin, ref, root/'p5-comparison.json')
frozen = {}
for name, h in edge['output_checksums'].items():
    frozen[str(origin/'edges'/name)] = h
for name in ('input/proteins.parquet', 'input/species.parquet', 'components/index.parquet',
             'edges/edge-manifest.json'):
    frozen[str(origin/name)] = sha256_file(origin/name)
for name in ('Orthogroups.tsv', 'Orthogroups_UnassignedGenes.tsv'):
    frozen[str(reference/name)] = sha256_file(reference/name)
for label, policy, depth in [('current-default','kway_v1',20),
                              ('current-soft42','soft_binary_24_v2',42)]:
    run = root/label
    run.mkdir()
    for name in ('input','edges','components'):
        (run/name).symlink_to(origin/name, target_is_directory=True)
    for name in ('hierarchy','evolution','orthogroups','results','search'):
        (run/name).mkdir()
    config = __import__('copy').deepcopy(DEFAULT_CONFIG)
    config['hierarchy'].update(topology_policy=policy, max_depth=depth)
    validate_config(config)
    (run/'run.yaml').write_text(yaml.safe_dump(config, sort_keys=False))
(root/'frozen-inputs.json').write_text(json.dumps(frozen,indent=2))
(root/'evaluation-plan.json').write_text(json.dumps(dict(
    origin=str(origin), reference=str(reference), search_recomputed=False,
    ssn_recomputed=False, variants=['current-default','current-soft42'],
    soft42_is_explicit_experiment=True, production_defaults_changed=False,
    primary='OF assigned proteins only', secondary='OF unassigned as distinct singletons',
    openbench_used=False), indent=2))
PY
"$DEV_CONDA" list -n "$DEV_CONDA_ENV" > "$DEV_RUN_DIR/environment.txt"
FAILED=0
for LABEL in current-default current-soft42; do
    RUN="$DEV_RUN_DIR/$LABEL"
    SUCCESS=1
    for STAGE in hierarchy-all annotate-network orthogroups export; do
        if /usr/bin/time -v -o "$DEV_RUN_DIR/$LABEL-$STAGE.time" \
            python -m ogprofiler "$STAGE" --run "$RUN" \
            > "$DEV_RUN_DIR/$LABEL-$STAGE.stdout" 2> "$DEV_RUN_DIR/$LABEL-$STAGE.stderr"; then
            echo "$LABEL $STAGE completed"
        else
            echo "$LABEL $STAGE failed; no downstream score" >&2
            printf '%s\n' "$STAGE" > "$DEV_RUN_DIR/$LABEL-failed-stage.txt"
            SUCCESS=0; FAILED=1; break
        fi
    done
    if [[ "$SUCCESS" == 1 ]]; then
        python -m benchmarks.og_extraction.embleya_reference --run "$RUN" \
            --reference-dir "$DEV_RUN_DIR/reference" --out "$DEV_RUN_DIR/$LABEL-comparison.json"
    fi
done
python - <<'PY'
import json, os
from pathlib import Path
from ogprofiler.core.manifest import sha256_file
root=Path(os.environ['DEV_RUN_DIR'])
frozen=json.loads((root/'frozen-inputs.json').read_text())
checks={p:sha256_file(Path(p))==h for p,h in frozen.items()}
(root/'input-integrity.json').write_text(json.dumps(dict(passed=all(checks.values()),checks=checks),indent=2))
assert all(checks.values())
PY
exit "$FAILED"
