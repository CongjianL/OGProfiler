import csv
import os
import subprocess
from pathlib import Path


def test_timed_wrapper_uses_snapshot_commit(tmp_path):
    script = Path(__file__).resolve().parents[1] / "OGProfiler2_benchmark/workflows/run_timed.sh"
    # The production wrapper targets GNU time on the Linux cluster. Use a tiny
    # resource-capture fixture so the provenance regression also runs on macOS.
    timer = tmp_path / "time-fixture.sh"
    timer.write_text(
        "#!/bin/bash\nfile=$3; shift 3\n"
        'printf "Maximum resident set size (kbytes): 1\\n'
        "Elapsed (wall clock) time (h:mm:ss or m:ss): 0:00.01\\n"
        'User time (seconds): 0.01\\nSystem time (seconds): 0.00\\n" > "$file"\n'
        '"$@"\n'
    )
    timer.chmod(0o755)
    fixture_script = tmp_path / "wrapper.sh"
    fixture_script.write_text(script.read_text().replace("/usr/bin/time", str(timer)))
    out = tmp_path / "timed"
    env = dict(os.environ, DEV_GIT_COMMIT="fixture-snapshot-commit")
    subprocess.run(
        ["bash", str(fixture_script), "fixture", "QFO", "1", str(out), "--", "true"],
        env=env,
        check=True,
    )
    with (out / "metadata.tsv").open() as h:
        values = {r["key"]: r["value"] for r in csv.DictReader(h, delimiter="\t")}
    assert values["git_commit"] == "fixture-snapshot-commit"
    assert values["exit_code"] == "0"
