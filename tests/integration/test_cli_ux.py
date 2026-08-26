from __future__ import annotations

import json
from pathlib import Path

import pyarrow as pa
import pyarrow.parquet as pq

from ogprofiler.cli import main
from ogprofiler.ux import pipeline_commands


def _fake_diamond(path: Path) -> Path:
    executable = path / "diamond-fake"
    executable.write_text(
        """#!/usr/bin/env python3
import pathlib
import sys
args = sys.argv[1:]
if args[0] in {"version", "--version"}:
    print("diamond version fake-2.1")
elif args[0] == "makedb":
    pathlib.Path(args[args.index("--db") + 1]).with_suffix(".dmnd").write_text("db")
elif args[0] == "blastp":
    pathlib.Path(args[args.index("--out") + 1]).write_text(
        "OGP2P000000000000\\tOGP2P000000000001\\t80\\t4\\t4\\t4\\t1e-10\\t50\\n"
        "OGP2P000000000001\\tOGP2P000000000000\\t80\\t4\\t4\\t4\\t1e-10\\t50\\n"
    )
else:
    raise SystemExit(2)
""",
        encoding="utf-8",
    )
    executable.chmod(0o755)
    return executable


def test_run_dry_run_expands_standard_pipeline(tmp_path: Path, capsys: object) -> None:
    proteomes = tmp_path / "proteomes"
    proteomes.mkdir()
    commands = pipeline_commands(
        tmp_path / "run",
        proteomes,
        config=None,
        overrides=("runtime.workers=2",),
        from_stage="prepare",
        until_stage="export",
    )
    assert [command[0] for command in commands] == [
        "prepare",
        "search",
        "edges",
        "components",
        "hierarchy-all",
        "annotate-network",
        "export",
    ]
    assert commands[-1] == ("export", "--run", str(tmp_path / "run"))
    assert (
        main(
            [
                "run",
                "--proteomes",
                str(proteomes),
                "--out",
                str(tmp_path / "run"),
                "--dry-run",
            ]
        )
        == 0
    )
    output = capsys.readouterr().out  # type: ignore[attr-defined]
    assert "ogprofiler hierarchy-all" in output
    assert "ogprofiler export" in output


def test_status_and_inspect_report_workspace_and_component(
    tmp_path: Path, capsys: object
) -> None:
    run = tmp_path / "run"
    (run / "input").mkdir(parents=True)
    (run / "hierarchy/components/component=00000000").mkdir(parents=True)
    (run / "manifest.json").write_text(
        json.dumps({"ogprofiler_version": "test"}), encoding="utf-8"
    )
    pq.write_table(
        pa.table({"protein_id": [0, 1], "species_id": [0, 1]}),
        run / "input/proteins.parquet",
    )
    component = run / "hierarchy/components/component=00000000"
    pq.write_table(
        pa.table(
            {
                "cluster_id": [0, 1],
                "depth": [0, 1],
                "terminal_reason": [None, "MIN_SIZE"],
            }
        ),
        component / "nodes.parquet",
    )
    pq.write_table(pa.table({"protein_id": [0, 1]}), component / "members.parquet")

    assert main(["status", "--run", str(run), "--json"]) == 0
    status = json.loads(capsys.readouterr().out)  # type: ignore[attr-defined]
    assert status["overall"] == "IN_PROGRESS"
    assert status["stages"][0]["status"] == "DONE"
    assert status["stages"][1]["status"] == "NEXT"

    assert main(["inspect", "--run", str(run), "--component", "0", "--json"]) == 0
    report = json.loads(capsys.readouterr().out)  # type: ignore[attr-defined]
    assert report["counts"]["proteins"] == 2
    assert report["component"] == {
        "component_id": 0,
        "max_depth": 1,
        "members": 2,
        "nodes": 2,
        "terminal_families": 1,
    }


def test_run_executes_stage_range_and_captures_provenance(tmp_path: Path) -> None:
    proteomes = tmp_path / "proteomes"
    proteomes.mkdir()
    (proteomes / "a.faa").write_text(">a\nAAAA\n", encoding="utf-8")
    (proteomes / "b.faa").write_text(">b\nAAAA\n", encoding="utf-8")
    executable = _fake_diamond(tmp_path)
    run = tmp_path / "run"
    command = [
        "run",
        "--proteomes",
        str(proteomes),
        "--out",
        str(run),
        "--until-stage",
        "components",
        "--set",
        f"search.executable={executable}",
    ]
    assert main(command) == 0
    assert (run / "components/component-manifest.json").is_file()
    provenance = json.loads((run / "provenance.json").read_text(encoding="utf-8"))
    assert provenance["ogprofiler_version"]
    assert provenance["random_seed"] == 42
    assert provenance["external_tools"]["search"]["available"] is True
