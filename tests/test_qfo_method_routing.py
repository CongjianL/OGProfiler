import os
import shutil
import subprocess
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]


@pytest.mark.parametrize(
    "method,expected",
    [
        ("v2", ["QFO_stage", "OGProfiler2_soft42"]),
        ("of3", ["QFO_stage", "OrthoFinder3", "OrthoFinder3_conversion"]),
        ("v1", ["QFO_stage", "OGProfiler1First", "OGProfiler1First_conversion"]),
        (
            "all",
            [
                "QFO_stage",
                "OGProfiler2_soft42",
                "OrthoFinder3",
                "OrthoFinder3_conversion",
                "OGProfiler1First",
                "OGProfiler1First_conversion",
            ],
        ),
    ],
)
def test_method_routing_schedules_only_selected_arm(tmp_path, method, expected):
    source = tmp_path / "source"
    b = source / "OGProfiler2_benchmark"
    (b / "workflows").mkdir(parents=True)
    (b / "00_env/ogprofiler_v1_first").mkdir(parents=True)
    shutil.copy2(
        ROOT / "OGProfiler2_benchmark/00_env/ogprofiler_v1_first/OGProfiler.py",
        b / "00_env/ogprofiler_v1_first/OGProfiler.py",
    )
    timed = b / "workflows/run_timed.sh"
    timed.write_text("""#!/usr/bin/env bash
set -eu
echo "$1" >> "$DEV_RUN_DIR/calls.txt"
mkdir -p "$DEV_RUN_DIR/campaign/of3/Results/Orthogroups"
mkdir -p "$DEV_RUN_DIR/campaign/v1-input/WorkingDirectory"
printf 'Orthogroup\\n' > "$DEV_RUN_DIR/campaign/of3/Results/Orthogroups/Orthogroups.tsv"
echo fixture > "$DEV_RUN_DIR/campaign/v1-input/WorkingDirectory/OGFile_coalescence_SameGenome.txt"
""")
    timed.chmod(0o755)
    prefix = tmp_path / "prefix"
    (prefix / "bin").mkdir(parents=True)
    for name in ("python", "diamond"):
        p = prefix / "bin" / name
        p.write_text("#!/bin/sh\nexit 0\n")
        p.chmod(0o755)
    manager = tmp_path / "manager"
    manager.write_text('#!/bin/sh\nif [ "$1" = run ]; then echo "$MOCK_PREFIX"; fi\n')
    manager.chmod(0o755)
    run = tmp_path / "run"
    run.mkdir()
    env = dict(
        os.environ,
        QFO_METHOD=method,
        QFO_METHODS_READY="1",
        DEV_SOURCE_DIR=str(source),
        DEV_RUN_DIR=str(run),
        DEV_CONDA=str(manager),
        DEV_CONDA_ENV="v2-env",
        SLURM_CPUS_PER_TASK="32",
        MOCK_PREFIX=str(prefix),
    )
    subprocess.run(
        [
            "bash",
            str(ROOT / "benchmarks/qfo/run_methods.sh"),
            "full",
            "ORIGIN",
            "DIGEST",
            "of-env",
            "v1-env",
            "1411874",
            "1411889",
        ],
        env=env,
        check=True,
    )
    assert (run / "calls.txt").read_text().splitlines() == expected
    assert "1411874" in (run / "method-validation-provenance.tsv").read_text()
    assert "1411889" in (run / "parallel-provenance.tsv").read_text()
    import json

    assert json.loads((run / "qfo-methods-completion.json").read_text())["method"] == method


def test_unknown_method_fails_before_any_execution(tmp_path):
    result = subprocess.run(
        [
            "bash",
            str(ROOT / "benchmarks/qfo/run_methods.sh"),
            "full",
            "ORIGIN",
            "DIGEST",
            "OF",
            "V1",
        ],
        env=dict(os.environ, QFO_METHOD="typo"),
        capture_output=True,
    )
    assert result.returncode == 2
    assert b"Invalid method" in result.stderr
