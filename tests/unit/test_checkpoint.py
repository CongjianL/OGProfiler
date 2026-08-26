from __future__ import annotations

import os
import sqlite3
import subprocess
import sys
from pathlib import Path

from ogprofiler.core.checkpoint import CheckpointStore
from ogprofiler.core.workspace import Workspace


def test_checkpoint_recovers_running_and_audits_invalidation(tmp_path: Path) -> None:
    workspace = Workspace.create(tmp_path / "run")
    store = CheckpointStore(workspace.database_path)
    first = store.register("hierarchy", "7", "input-a", "algorithm-a")
    assert first.status == "PENDING"
    store.start("hierarchy", "7")
    assert store.get("hierarchy", "7").status == "RUNNING"  # type: ignore[union-attr]

    recovered = store.register("hierarchy", "7", "input-a", "algorithm-a")
    assert recovered.status == "PENDING"
    assert recovered.error == "recovered stale RUNNING task"
    store.start("hierarchy", "7")
    store.finish("hierarchy", "7", "component=00000007")
    assert store.get("hierarchy", "7").status == "DONE"  # type: ignore[union-attr]

    invalid = store.register("hierarchy", "7", "input-b", "algorithm-a")
    assert invalid.status == "INVALID"
    assert invalid.attempts == 2
    with sqlite3.connect(workspace.database_path) as database:
        statuses = [
            row[0]
            for row in database.execute(
                "SELECT status FROM task_events WHERE stage='hierarchy' AND task_id='7' "
                "ORDER BY event_id"
            )
        ]
    assert statuses == ["PENDING", "RUNNING", "PENDING", "RUNNING", "DONE", "INVALID"]


def test_workspace_migrates_pre_phase7_checkpoint_schema(tmp_path: Path) -> None:
    root = tmp_path / "legacy-run"
    root.mkdir()
    database_path = root / "run.db"
    with sqlite3.connect(database_path) as database:
        database.execute(
            """
            CREATE TABLE tasks (
                stage TEXT NOT NULL, task_id TEXT NOT NULL, status TEXT NOT NULL,
                started TEXT, completed TEXT, input_hash TEXT, output_path TEXT,
                algorithm_version TEXT NOT NULL, PRIMARY KEY(stage, task_id)
            )
            """
        )
        database.execute(
            "INSERT INTO tasks VALUES(?,?,?,?,?,?,?,?)",
            ("hierarchy", "2", "DONE", "start", "end", "hash", "out", "v0"),
        )
    workspace = Workspace.create(root)
    migrated = CheckpointStore(workspace.database_path).get("hierarchy", "2")
    assert migrated is not None
    assert migrated.started_at == "start"
    assert migrated.completed_at == "end"
    assert migrated.attempts == 0


def test_killed_process_leaves_running_task_recoverable(tmp_path: Path) -> None:
    workspace = Workspace.create(tmp_path / "killed-run")
    store = CheckpointStore(workspace.database_path)
    store.register("hierarchy", "3", "identity", "v1")
    script = (
        "import sys,time; "
        "from pathlib import Path; "
        "from ogprofiler.core.checkpoint import CheckpointStore; "
        "s=CheckpointStore(Path(sys.argv[1])); s.start('hierarchy','3'); "
        "print('running',flush=True); time.sleep(60)"
    )
    process = subprocess.Popen(
        [sys.executable, "-c", script, str(workspace.database_path)],
        stdout=subprocess.PIPE,
        text=True,
        env={**os.environ, "PYTHONPATH": str(Path.cwd() / "src")},
    )
    try:
        assert process.stdout is not None
        assert process.stdout.readline().strip() == "running"
        process.terminate()
        process.wait(timeout=5)
    finally:
        if process.poll() is None:
            process.kill()
    recovered = store.register("hierarchy", "3", "identity", "v1")
    assert recovered.status == "PENDING"
    assert recovered.error == "recovered stale RUNNING task"
