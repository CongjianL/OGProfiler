"""User-facing pipeline planning, status reporting, and run inspection."""

from __future__ import annotations

import json
import platform
import sqlite3
import subprocess
import sys
from collections import Counter
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any

import pyarrow.parquet as pq

from ogprofiler import __version__
from ogprofiler.core.manifest import sha256_file, write_json
from ogprofiler.exceptions import CheckpointError, InputError

PIPELINE_STAGES = (
    "prepare",
    "search",
    "edges",
    "components",
    "hierarchy",
    "annotate-network",
    "export",
)


@dataclass(frozen=True, slots=True)
class StageStatus:
    stage: str
    status: str
    artifact: str


def pipeline_commands(
    run_root: Path,
    proteomes: Path | None,
    *,
    config: Path | None,
    overrides: tuple[str, ...],
    from_stage: str,
    until_stage: str,
) -> tuple[tuple[str, ...], ...]:
    try:
        start = PIPELINE_STAGES.index(from_stage)
        stop = PIPELINE_STAGES.index(until_stage)
    except ValueError as error:
        raise InputError(f"Unknown pipeline stage: {error}") from error
    if start > stop:
        raise InputError("--from-stage must not follow --until-stage")
    if start == 0 and proteomes is None and not (run_root / "manifest.json").is_file():
        raise InputError("run requires --proteomes when creating a workspace")

    common: list[str] = []
    if config is not None:
        common.extend(["--config", str(config)])
    for override in overrides:
        common.extend(["--set", override])
    commands: list[tuple[str, ...]] = []
    for stage in PIPELINE_STAGES[start : stop + 1]:
        if stage == "prepare":
            if (run_root / "manifest.json").is_file() and proteomes is None:
                continue
            if proteomes is None:
                raise InputError("prepare requires --proteomes")
            command = ["prepare", "--proteomes", str(proteomes), "--out", str(run_root)]
            command.extend(common)
        elif stage == "hierarchy":
            command = ["hierarchy-all", "--run", str(run_root), *common]
        elif stage == "export":
            command = ["export", "--run", str(run_root)]
        else:
            command = [stage, "--run", str(run_root), *common]
        commands.append(tuple(command))
    return tuple(commands)


def _artifact_status(run_root: Path) -> tuple[StageStatus, ...]:
    artifacts = (
        ("prepare", run_root / "manifest.json"),
        ("search", run_root / "search/search-manifest.json"),
        ("edges", run_root / "edges/edge-manifest.json"),
        ("components", run_root / "components/component-manifest.json"),
        ("hierarchy", run_root / "hierarchy/scheduler-manifest.json"),
        ("annotate-network", run_root / "evolution/network-event-manifest.json"),
        ("export", run_root / "results/export-manifest.json"),
    )
    first_missing = False
    rows: list[StageStatus] = []
    for stage, path in artifacts:
        if path.is_file():
            status = "DONE"
        elif first_missing:
            status = "PENDING"
        else:
            status = "NEXT"
            first_missing = True
        rows.append(StageStatus(stage, status, path.relative_to(run_root).as_posix()))
    return tuple(rows)


def run_status(run_root: Path) -> dict[str, Any]:
    if not run_root.is_dir():
        raise CheckpointError(f"Run workspace does not exist: {run_root}")
    stages = _artifact_status(run_root)
    tasks: Counter[str] = Counter()
    database = run_root / "run.db"
    if database.is_file():
        try:
            with sqlite3.connect(database) as connection:
                for status, count in connection.execute(
                    "SELECT status, COUNT(*) FROM tasks GROUP BY status"
                ):
                    tasks[str(status)] = int(count)
        except sqlite3.Error as error:
            raise CheckpointError(f"Failed to inspect {database}: {error}") from error
    overall = "FAILED" if tasks["FAILED"] else (
        "COMPLETE" if all(row.status == "DONE" for row in stages) else "IN_PROGRESS"
    )
    return {
        "run_root": str(run_root.resolve()),
        "overall": overall,
        "stages": [asdict(row) for row in stages],
        "task_counts": dict(sorted(tasks.items())),
    }


def inspect_run(run_root: Path, component_id: int | None = None) -> dict[str, Any]:
    report = run_status(run_root)
    manifest_path = run_root / "manifest.json"
    if manifest_path.is_file():
        report["manifest"] = json.loads(manifest_path.read_text(encoding="utf-8"))
    counts: dict[str, int] = {}
    table_specs = {
        "proteins": run_root / "input/proteins.parquet",
        "components": run_root / "components/statistics.parquet",
        "families": run_root / "results/families.tsv",
        "members": run_root / "results/members.tsv",
    }
    for name, path in table_specs.items():
        if not path.is_file():
            continue
        if path.suffix == ".parquet":
            counts[name] = pq.read_metadata(path).num_rows
        else:
            with path.open(encoding="utf-8") as handle:
                counts[name] = max(0, sum(1 for _ in handle) - 1)
    report["counts"] = counts
    if component_id is not None:
        component = run_root / "hierarchy/components" / f"component={component_id:08d}"
        nodes = component / "nodes.parquet"
        members = component / "members.parquet"
        if not nodes.is_file() or not members.is_file():
            raise InputError(f"Hierarchy component is absent: {component_id}")
        node_rows = pq.read_table(nodes).to_pylist()
        report["component"] = {
            "component_id": component_id,
            "nodes": len(node_rows),
            "terminal_families": sum(row["terminal_reason"] is not None for row in node_rows),
            "members": pq.read_metadata(members).num_rows,
            "max_depth": max((int(row["depth"]) for row in node_rows), default=0),
        }
    return report


def _tool_version(
    executable: str, version_args: tuple[str, ...] = ("--version",)
) -> dict[str, Any]:
    try:
        completed = subprocess.run(
            [executable, *version_args],
            check=False,
            capture_output=True,
            text=True,
            timeout=10,
        )
        output = (completed.stdout or completed.stderr).strip().splitlines()
        return {
            "executable": executable,
            "available": completed.returncode == 0,
            "version": output[0] if output else None,
        }
    except (OSError, subprocess.TimeoutExpired):
        return {"executable": executable, "available": False, "version": None}


def write_run_provenance(run_root: Path, config: dict[str, Any], command: list[str]) -> Path:
    """Capture runtime and external tool versions without changing scientific identity."""
    run_yaml = run_root / "run.yaml"
    manifest = run_root / "manifest.json"
    path = run_root / "provenance.json"
    search_backend = str(config["search"]["backend"])
    search_version_args = {
        "diamond": ("version",),
        "mmseqs": ("version",),
        "blastp": ("-version",),
    }[search_backend]
    tools = {
        "search": (
            str(config["search"]["executable"]),
            search_version_args,
        ),
        "alignment": (str(config["phylogeny"]["alignment_executable"]), ("--version",)),
        "tree": (str(config["phylogeny"]["tree_executable"]), ("--version",)),
    }
    write_json(
        path,
        {
            "ogprofiler_version": __version__,
            "python_version": platform.python_version(),
            "platform": platform.platform(),
            "command": command,
            "random_seed": int(config["hierarchy"]["seed"]),
            "run_yaml_sha256": sha256_file(run_yaml) if run_yaml.is_file() else None,
            "manifest_sha256": sha256_file(manifest) if manifest.is_file() else None,
            "external_tools": {
                name: _tool_version(executable, version_args)
                for name, (executable, version_args) in tools.items()
            },
            "python_executable": sys.executable,
        },
    )
    return path
