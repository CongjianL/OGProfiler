#!/usr/bin/env bash
# Same resource envelope as the existing full Open Orthobench modified run.
# Single campaign, no array, no new similarity search or parameter sweep.
#SBATCH --job-name=ogp-p6-orthobench
#SBATCH --partition=cu
#SBATCH --account=mselab
#SBATCH --qos=normal
#SBATCH --cpus-per-task=56
#SBATCH --mem=250G
#SBATCH --time=72:00:00
set -euo pipefail
: "${DEV_SOURCE_DIR:?Missing immutable source}"
: "${DEV_RUN_DIR:?Missing independent run directory}"
: "${DEV_CONDA:?Missing configured manager}"
: "${DEV_CONDA_ENV:?Missing configured environment}"
ORIGIN=${1:?Usage: p6_orthobench.sh HISTORICAL_RUN ORTHOBENCH_ROOT V1_GROUPS V1_MANIFEST}
OB=${2:?Missing Orthobench root}
V1_GROUPS=${3:?Missing standardized historical V1 groups}
V1_MANIFEST=${4:?Missing historical V1 run manifest}
if [[ ${OGP_P6_ENV_READY:-0} != 1 ]]; then
    exec "$DEV_CONDA" run -n "$DEV_CONDA_ENV" env OGP_P6_ENV_READY=1 bash "$0" "$@"
fi
export PYTHONPATH="$DEV_SOURCE_DIR/src:$DEV_SOURCE_DIR"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
python - "$ORIGIN" "$OB" "$V1_GROUPS" "$V1_MANIFEST" <<'PY'
import json, os, runpy, shutil, sys
from pathlib import Path
import pyarrow.parquet as pq
from ogprofiler.config import load_config
from ogprofiler.core.manifest import sha256_file
from benchmarks.og_extraction.orthobench import OFFICIAL_SCORER_SHA256
import yaml
origin, ob, v1_groups, v1_manifest = map(Path, sys.argv[1:])
root = Path(os.environ['DEV_RUN_DIR'])
assert sha256_file(origin/'search/hits.parquet') == '29bc19b7efa89ae567e047d5c180e27edd916ca407d154b90e4a05e494904c1e'
assert sha256_file(origin/'input/proteins.parquet') == 'e5323530d8f6acbdf6cab152e46869712e1965ff2cdda1b69300f6f245db106f'
metadata = pq.ParquetFile(origin/'input/proteins.parquet').read().to_pylist()
assert len(metadata) == 251378 and len({r['species_id'] for r in metadata}) == 12
assert len({r['original_id'] for r in metadata}) == len(metadata)
assert (origin/'hierarchy/scheduler-manifest.json').is_file()
snapshot = root/'benchmark'
snapshot.mkdir()
shutil.copy2(ob/'BENCHMARKS/benchmark.py', snapshot/'benchmark.py')
assert sha256_file(snapshot/'benchmark.py') == OFFICIAL_SCORER_SHA256
shutil.copytree(ob/'BENCHMARKS/Input', snapshot/'Input')
(snapshot/'RefOGs').mkdir()
for p in (ob/'BENCHMARKS/RefOGs').rglob('*.txt'):
    target = snapshot/'RefOGs'/p.relative_to(ob/'BENCHMARKS/RefOGs')
    target.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(p, target)
official = runpy.run_path(str(snapshot/'benchmark.py'))
assert {r['original_id'] for r in metadata} == official['get_expected_genes']()
assert sum(len(r) for r in official['read_refogs'](str(snapshot/'RefOGs')+'/')) == 1945
preflight = root/'preflight-original-ids.txt'
preflight.write_text('input: '+' '.join(r['original_id'] for r in metadata)+'\n')
parsed = official['read_orthogroups_smart'](str(preflight))
assert len(parsed) == 1 and parsed[0] == {r['original_id'] for r in metadata}
shutil.copy2(v1_groups, root/'historical-v1-groups.tsv')
shutil.copy2(v1_manifest, root/'historical-v1-manifest.tsv')
for arm, dirs in [('historical-frozen', ('input', 'edges', 'components', 'hierarchy', 'evolution')),
                  ('repaired-mean', ('input', 'search'))]:
    run = root/arm
    run.mkdir()
    for name in dirs:
        if name == 'search':
            (run/name).mkdir()
            shutil.copy2(origin/'search/hits.parquet', run/'search/hits.parquet')
        else:
            shutil.copytree(origin/name, run/name)
    for name in ('search', 'edges', 'components', 'hierarchy', 'evolution', 'results'):
        (run/name).mkdir(exist_ok=True)
    config = load_config(str(origin/'run.yaml'), [])
    if arm == 'repaired-mean':
        config['edges']['symmetrization'] = 'mean'
    (run/'run.yaml').write_text(yaml.safe_dump(config, sort_keys=True))
hashes = {str(p.relative_to(snapshot)): sha256_file(p) for p in snapshot.rglob('*') if p.is_file()}
(root/'benchmark-inputs.sha256.json').write_text(json.dumps(hashes, indent=2))
(root/'p6-input-origin.json').write_text(json.dumps(dict(origin=str(origin), orthobench=str(ob),
    v1_groups=str(v1_groups), v1_manifest=str(v1_manifest), new_search=False,
    arms=['frozen historical forward hierarchy', 'repaired mean fixed hits hierarchy'],
    algorithm_commit=os.environ['DEV_GIT_COMMIT'], dirty=os.environ['DEV_GIT_DIRTY'],
    source_hash=os.environ['DEV_SOURCE_HASH']), indent=2))
PY
python "$DEV_SOURCE_DIR/OGProfiler2_benchmark/scripts/audit/audit_open_orthobench_input.py" \
    --input "$OB/BENCHMARKS/Input" --refogs "$OB/BENCHMARKS/RefOGs" \
    --out "$DEV_RUN_DIR/input-audit/dataset.tsv"
[[ $(cat "$DEV_RUN_DIR/input-audit/OPEN_ORTHOBENCH_INPUT_DIGEST") == 380c85d9548c607f5df6daf656d74b21597fe215ef78d4f8ff296776bb1fdd07 ]]
"$DEV_CONDA" list -n "$DEV_CONDA_ENV" > "$DEV_RUN_DIR/environment.txt"
for stage in edges components hierarchy-all annotate-network; do
    echo "==> Repaired mean: $stage (no search)"
    /usr/bin/time -v -o "$DEV_RUN_DIR/repaired-$stage.time" \
        python -m ogprofiler "$stage" --run "$DEV_RUN_DIR/repaired-mean" \
        > "$DEV_RUN_DIR/repaired-$stage.stdout" 2> "$DEV_RUN_DIR/repaired-$stage.stderr"
done
for arm in historical-frozen repaired-mean; do
    echo "==> Fixed $arm: independent V1 parity, OG artifacts and terminal export"
    /usr/bin/time -v -o "$DEV_RUN_DIR/$arm-audit.time" \
        python -m benchmarks.og_extraction.fixed_hierarchy \
        --run "$DEV_RUN_DIR/$arm" --out "$DEV_RUN_DIR/$arm-audit" \
        > "$DEV_RUN_DIR/$arm-audit.stdout" 2> "$DEV_RUN_DIR/$arm-audit.stderr"
    python -m ogprofiler export --run "$DEV_RUN_DIR/$arm" --strategy terminal \
        > "$DEV_RUN_DIR/$arm-terminal.stdout" 2> "$DEV_RUN_DIR/$arm-terminal.stderr"
done
/usr/bin/time -v -o "$DEV_RUN_DIR/scoring.time" \
    python -m benchmarks.og_extraction.orthobench \
    --benchmark "$DEV_RUN_DIR/benchmark" --frozen "$DEV_RUN_DIR/historical-frozen" \
    --repaired "$DEV_RUN_DIR/repaired-mean" --v1-groups "$DEV_RUN_DIR/historical-v1-groups.tsv" \
    --v1-manifest "$DEV_RUN_DIR/historical-v1-manifest.tsv" --out "$DEV_RUN_DIR/p6-metrics"
python - <<'PY'
import json, os
from pathlib import Path
from ogprofiler.core.manifest import sha256_file
root = Path(os.environ['DEV_RUN_DIR'])
hashes = json.loads((root/'benchmark-inputs.sha256.json').read_text())
assert all(sha256_file(root/'benchmark'/p) == h for p, h in hashes.items())
assert all(json.loads((root/f'{arm}-audit/report.json').read_text())['passed']
           for arm in ('historical-frozen', 'repaired-mean'))
report = json.loads((root/'p6-metrics/report.json').read_text())
assert report['evaluation_completed'] and len(report['results']) == 5
(root/'p6-completion.json').write_text(json.dumps(dict(job_id=os.environ['SLURM_JOB_ID'],
    run_id=os.environ['DEV_RUN_ID'], benchmark_unchanged=True, evaluation_completed=True,
    strategy_parity_passed=True, accuracy_improvement_claimed=False), indent=2))
PY
