"""Run-workspace creation and checkpoint schema initialization."""

from __future__ import annotations

import sqlite3
from dataclasses import dataclass
from pathlib import Path

from ogprofiler.exceptions import CheckpointError

WORKSPACE_DIRECTORIES = (
    "input",
    "search",
    "edges",
    "components",
    "hierarchy",
    "evolution",
    "results",
)


@dataclass(frozen=True, slots=True)
class Workspace:
    root: Path

    @classmethod
    def create(cls, root: Path) -> Workspace:
        root = root.expanduser().resolve()
        if root.exists() and not root.is_dir():
            raise CheckpointError(f"Workspace path is not a directory: {root}")
        root.mkdir(parents=True, exist_ok=True)
        for name in WORKSPACE_DIRECTORIES:
            (root / name).mkdir(exist_ok=True)
        workspace = cls(root)
        workspace.initialize_database()
        return workspace

    @property
    def database_path(self) -> Path:
        return self.root / "run.db"

    def initialize_database(self) -> None:
        try:
            with sqlite3.connect(self.database_path) as database:
                database.execute(
                    """
                    CREATE TABLE IF NOT EXISTS tasks (
                        stage TEXT NOT NULL,
                        task_id TEXT NOT NULL,
                        status TEXT NOT NULL CHECK (
                            status IN ('PENDING', 'RUNNING', 'DONE', 'FAILED', 'INVALID')
                        ),
                        started_at TEXT,
                        completed_at TEXT,
                        input_hash TEXT,
                        output_path TEXT,
                        error TEXT,
                        attempts INTEGER NOT NULL DEFAULT 0,
                        algorithm_version TEXT NOT NULL,
                        PRIMARY KEY (stage, task_id)
                    )
                    """
                )
                columns = {
                    str(row[1]) for row in database.execute("PRAGMA table_info(tasks)").fetchall()
                }
                if "error" not in columns:
                    database.execute("ALTER TABLE tasks ADD COLUMN error TEXT")
                if "attempts" not in columns:
                    database.execute(
                        "ALTER TABLE tasks ADD COLUMN attempts INTEGER NOT NULL DEFAULT 0"
                    )
                if "started_at" not in columns:
                    database.execute("ALTER TABLE tasks ADD COLUMN started_at TEXT")
                    if "started" in columns:
                        database.execute("UPDATE tasks SET started_at=started")
                if "completed_at" not in columns:
                    database.execute("ALTER TABLE tasks ADD COLUMN completed_at TEXT")
                    if "completed" in columns:
                        database.execute("UPDATE tasks SET completed_at=completed")
                database.execute(
                    """
                    CREATE TABLE IF NOT EXISTS task_events (
                        event_id INTEGER PRIMARY KEY AUTOINCREMENT,
                        stage TEXT NOT NULL,
                        task_id TEXT NOT NULL,
                        status TEXT NOT NULL,
                        timestamp TEXT NOT NULL,
                        detail TEXT
                    )
                    """
                )
        except sqlite3.Error as error:
            raise CheckpointError(f"Failed to initialize {self.database_path}: {error}") from error
