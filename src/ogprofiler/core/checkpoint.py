"""SQLite checkpoint transitions owned exclusively by the scheduler process."""

from __future__ import annotations

import sqlite3
from dataclasses import dataclass, replace
from datetime import datetime, timezone
from pathlib import Path

from ogprofiler.exceptions import CheckpointError

TASK_STATUSES = {"PENDING", "RUNNING", "DONE", "FAILED", "INVALID", "UNRESOLVED"}


def _now() -> str:
    return datetime.now(timezone.utc).isoformat()


@dataclass(frozen=True, slots=True)
class TaskRecord:
    stage: str
    task_id: str
    status: str
    started_at: str | None
    completed_at: str | None
    input_hash: str | None
    output_path: str | None
    error: str | None
    attempts: int
    algorithm_version: str


class CheckpointStore:
    def __init__(self, path: Path) -> None:
        self.path = path

    def _connect(self) -> sqlite3.Connection:
        connection = sqlite3.connect(self.path, timeout=30)
        connection.row_factory = sqlite3.Row
        return connection

    def get(self, stage: str, task_id: str) -> TaskRecord | None:
        try:
            with self._connect() as database:
                row = database.execute(
                    """
                    SELECT stage,task_id,status,started_at,completed_at,input_hash,
                           output_path,error,attempts,algorithm_version
                    FROM tasks WHERE stage=? AND task_id=?
                    """,
                    (stage, task_id),
                ).fetchone()
        except sqlite3.Error as error:
            raise CheckpointError(f"Failed reading checkpoint: {error}") from error
        return TaskRecord(**dict(row)) if row is not None else None

    def list(self, stage: str) -> tuple[TaskRecord, ...]:
        try:
            with self._connect() as database:
                rows = database.execute(
                    """
                    SELECT stage,task_id,status,started_at,completed_at,input_hash,
                           output_path,error,attempts,algorithm_version
                    FROM tasks WHERE stage=? ORDER BY CAST(task_id AS INTEGER)
                    """,
                    (stage,),
                ).fetchall()
        except sqlite3.Error as error:
            raise CheckpointError(f"Failed listing checkpoints: {error}") from error
        return tuple(TaskRecord(**dict(row)) for row in rows)

    def register(
        self, stage: str, task_id: str, input_hash: str, algorithm_version: str
    ) -> TaskRecord:
        current = self.get(stage, task_id)
        if current is None:
            self._upsert(
                TaskRecord(
                    stage,
                    task_id,
                    "PENDING",
                    None,
                    None,
                    input_hash,
                    None,
                    None,
                    0,
                    algorithm_version,
                ),
                "registered",
            )
        elif current.input_hash != input_hash or current.algorithm_version != algorithm_version:
            self._upsert(
                TaskRecord(
                    stage,
                    task_id,
                    "INVALID",
                    current.started_at,
                    current.completed_at,
                    input_hash,
                    current.output_path,
                    "algorithm/config/input identity changed",
                    current.attempts,
                    algorithm_version,
                ),
                "identity changed",
            )
        elif current.status == "RUNNING":
            self._upsert(
                TaskRecord(
                    stage,
                    task_id,
                    "PENDING",
                    None,
                    None,
                    input_hash,
                    current.output_path,
                    "recovered stale RUNNING task",
                    current.attempts,
                    algorithm_version,
                ),
                "recovered stale RUNNING task",
            )
        result = self.get(stage, task_id)
        if result is None:
            raise CheckpointError("Checkpoint registration disappeared")
        return result

    def invalidate(self, stage: str, task_id: str, detail: str) -> None:
        current = self._required(stage, task_id)
        self._upsert(
            TaskRecord(
                current.stage,
                current.task_id,
                "INVALID",
                current.started_at,
                current.completed_at,
                current.input_hash,
                current.output_path,
                detail,
                current.attempts,
                current.algorithm_version,
            ),
            detail,
        )

    def start(self, stage: str, task_id: str) -> None:
        current = self._required(stage, task_id)
        if current.status not in {"PENDING", "FAILED", "INVALID"}:
            raise CheckpointError(f"Task {stage}/{task_id} is not runnable: {current.status}")
        self._upsert(
            TaskRecord(
                current.stage,
                current.task_id,
                "RUNNING",
                _now(),
                None,
                current.input_hash,
                current.output_path,
                None,
                current.attempts + 1,
                current.algorithm_version,
            ),
            "worker started",
        )

    def finish(self, stage: str, task_id: str, output_path: str) -> None:
        current = self._required(stage, task_id)
        if current.status not in {"PENDING", "RUNNING", "FAILED", "INVALID"}:
            raise CheckpointError(f"Task {stage}/{task_id} cannot finish from {current.status}")
        self._upsert(
            TaskRecord(
                current.stage,
                current.task_id,
                "DONE",
                current.started_at,
                _now(),
                current.input_hash,
                output_path,
                None,
                current.attempts,
                current.algorithm_version,
            ),
            "worker completed",
        )

    def unresolved(self, stage: str, task_id: str, output_path: str) -> None:
        current = self._required(stage, task_id)
        self._upsert(
            replace(
                current,
                status="UNRESOLVED",
                completed_at=_now(),
                output_path=output_path,
                error="hierarchy search unresolved",
            ),
            "diagnostics persisted; scientific completion pending",
        )

    def fail(self, stage: str, task_id: str, error: str) -> None:
        current = self._required(stage, task_id)
        if current.status != "RUNNING":
            raise CheckpointError(f"Task {stage}/{task_id} cannot fail from {current.status}")
        self._upsert(
            TaskRecord(
                current.stage,
                current.task_id,
                "FAILED",
                current.started_at,
                _now(),
                current.input_hash,
                current.output_path,
                error,
                current.attempts,
                current.algorithm_version,
            ),
            error,
        )

    def _required(self, stage: str, task_id: str) -> TaskRecord:
        value = self.get(stage, task_id)
        if value is None:
            raise CheckpointError(f"Unknown checkpoint task: {stage}/{task_id}")
        return value

    def _upsert(self, record: TaskRecord, detail: str) -> None:
        if record.status not in TASK_STATUSES:
            raise CheckpointError(f"Invalid checkpoint status: {record.status}")
        try:
            with self._connect() as database:
                database.execute(
                    """
                    INSERT INTO tasks (
                        stage, task_id, status, started_at, completed_at, input_hash,
                        output_path, error, attempts, algorithm_version
                    ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
                    ON CONFLICT(stage, task_id) DO UPDATE SET
                        status=excluded.status, started_at=excluded.started_at,
                        completed_at=excluded.completed_at, input_hash=excluded.input_hash,
                        output_path=excluded.output_path, error=excluded.error,
                        attempts=excluded.attempts, algorithm_version=excluded.algorithm_version
                    """,
                    (
                        record.stage,
                        record.task_id,
                        record.status,
                        record.started_at,
                        record.completed_at,
                        record.input_hash,
                        record.output_path,
                        record.error,
                        record.attempts,
                        record.algorithm_version,
                    ),
                )
                database.execute(
                    "INSERT INTO task_events(stage,task_id,status,timestamp,detail) "
                    "VALUES(?,?,?,?,?)",
                    (record.stage, record.task_id, record.status, _now(), detail),
                )
        except sqlite3.Error as error:
            raise CheckpointError(f"Failed updating checkpoint: {error}") from error
