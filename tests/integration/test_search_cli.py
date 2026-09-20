from __future__ import annotations

import json
from collections.abc import Callable
from pathlib import Path

import pyarrow.parquet as pq
import pytest

from ogprofiler.cli import main


def _fake_diamond(path: Path) -> Path:
    executable = path / "diamond-fake"
    executable.write_text(
        """#!/usr/bin/env python3
import pathlib
import sys

args = sys.argv[1:]
root = pathlib.Path(__file__).parent
counter = root / "diamond-invocations.txt"
with counter.open("a", encoding="utf-8") as handle:
    handle.write(args[0] + "\\n")
if args[0] == "version":
    print("diamond version fake-2.1")
elif args[0] == "makedb":
    database = pathlib.Path(args[args.index("--db") + 1]).with_suffix(".dmnd")
    database.parent.mkdir(parents=True, exist_ok=True)
    database.write_text("fake-db", encoding="utf-8")
elif args[0] == "blastp":
    output = pathlib.Path(args[args.index("--out") + 1])
    query = pathlib.Path(args[args.index("--query") + 1])
    db = pathlib.Path(args[args.index("--db") + 1])
    q_species = query.name.split("_")[1].split(".")[0]
    d_species = db.name.split("_")[1]
    fields = []
    if q_species == "0" and d_species == "1":
        fields = ["OGP2P000000000000", "OGP2P000000000001", "80", "4", "4", "4", "1e-10", "50"]
    elif q_species == "1" and d_species == "0":
        fields = ["OGP2P000000000001", "OGP2P000000000000", "80", "4", "4", "4", "1e-10", "50"]
    line = chr(9).join(fields) + chr(10) if fields else ""
    output.write_text(line, encoding="utf-8")
else:
    raise SystemExit(2)
""",
        encoding="utf-8",
    )
    executable.chmod(0o755)
    return executable


def _fake_mmseqs(path: Path) -> Path:
    executable = path / "mmseqs-fake"
    executable.write_text(
        """#!/usr/bin/env python3
import pathlib
import sys
args = sys.argv[1:]
if args[0] == "version":
    print("15.6f452")
elif args[0] == "createdb":
    pathlib.Path(args[2]).write_text("fake-db")
    pathlib.Path(args[2] + ".dbtype").write_text("0")
elif args[0] == "easy-search":
    query = pathlib.Path(args[1])
    db = pathlib.Path(args[2])
    output = pathlib.Path(args[3])
    q_species = query.name.split("_")[1].split(".")[0]
    d_species = db.name.split("_")[1]
    fields = []
    if q_species == "0" and d_species == "1":
        fields = ["OGP2P000000000000", "OGP2P000000000001", "80", "4", "4", "4", "1e-10", "50"]
    elif q_species == "1" and d_species == "0":
        fields = ["OGP2P000000000001", "OGP2P000000000000", "80", "4", "4", "4", "1e-10", "50"]
    line = chr(9).join(fields) + chr(10) if fields else ""
    output.write_text(line)
    pathlib.Path(args[4]).mkdir(parents=True)
else:
    raise SystemExit(2)
""",
        encoding="utf-8",
    )
    executable.chmod(0o755)
    return executable


def _fake_blast(path: Path) -> Path:
    blastp = path / "blastp-fake"
    blastp.write_text(
        """#!/usr/bin/env python3
import pathlib
import sys
args = sys.argv[1:]
if args[0] == "-version":
    print("blastp: 2.16.0+")
else:
    query = pathlib.Path(args[args.index("-query") + 1])
    db = pathlib.Path(args[args.index("-db") + 1])
    output = pathlib.Path(args[args.index("-out") + 1])
    q_species = query.name.split("_")[1].split(".")[0]
    d_species = db.name.split("_")[1]
    fields = []
    if q_species == "0" and d_species == "1":
        fields = ["OGP2P000000000000", "OGP2P000000000001", "80", "4", "4", "4", "1e-10", "50"]
    elif q_species == "1" and d_species == "0":
        fields = ["OGP2P000000000001", "OGP2P000000000000", "80", "4", "4", "4", "1e-10", "50"]
    line = chr(9).join(fields) + chr(10) if fields else ""
    output.write_text(line)
""",
        encoding="utf-8",
    )
    blastp.chmod(0o755)
    makeblastdb = path / "makeblastdb"
    makeblastdb.write_text(
        """#!/usr/bin/env python3
import pathlib
import sys
args = sys.argv[1:]
pathlib.Path(args[args.index("-out") + 1] + ".pin").write_text("fake-db")
""",
        encoding="utf-8",
    )
    makeblastdb.chmod(0o755)
    return blastp


@pytest.mark.parametrize(
    ("backend", "factory"),
    [("mmseqs", _fake_mmseqs), ("blastp", _fake_blast)],
)
def test_alternative_backends_emit_same_standard_hit_schema(
    tmp_path: Path, backend: str, factory: Callable[[Path], Path]
) -> None:
    proteomes = tmp_path / "proteomes"
    proteomes.mkdir()
    (proteomes / "a.faa").write_text(">a\nAAAA\n", encoding="utf-8")
    (proteomes / "b.faa").write_text(">b\nAAAA\n", encoding="utf-8")
    run = tmp_path / "run"
    assert main(["prepare", "--proteomes", str(proteomes), "--out", str(run)]) == 0
    executable = factory(tmp_path)
    command = [
        "search",
        "--run",
        str(run),
        "--backend",
        backend,
        "--set",
        f"search.executable={executable}",
    ]
    assert main(command) == 0
    table = pq.read_table(run / "search/hits.parquet")
    assert table.num_rows == 2
    assert table.column_names == [
        "query_id",
        "target_id",
        "query_species",
        "target_species",
        "bitscore",
        "identity",
        "query_coverage",
        "target_coverage",
        "evalue",
    ]
    manifest = json.loads((run / "search/search-manifest.json").read_text(encoding="utf-8"))
    assert manifest["backend"] == backend
    first_mtime = (run / "search/hits.parquet").stat().st_mtime_ns
    assert main(command) == 0
    assert (run / "search/hits.parquet").stat().st_mtime_ns == first_mtime


def test_search_cli_produces_manifest_and_verified_resume(tmp_path: Path) -> None:
    proteomes = tmp_path / "proteomes"
    proteomes.mkdir()
    (proteomes / "a.faa").write_text(">a\nAAAA\n", encoding="utf-8")
    (proteomes / "b.faa").write_text(">b\nAAAA\n", encoding="utf-8")
    run = tmp_path / "run"
    assert main(["prepare", "--proteomes", str(proteomes), "--out", str(run)]) == 0
    executable = _fake_diamond(tmp_path)
    command = [
        "search",
        "--run",
        str(run),
        "--backend",
        "diamond",
        "--set",
        f"search.executable={executable}",
        "--set",
        "search.max_target_seqs=0",
    ]

    assert main(command) == 0
    hits_path = run / "search" / "hits.parquet"
    manifest_path = run / "search" / "search-manifest.json"
    first_mtime = hits_path.stat().st_mtime_ns
    assert main(command) == 0
    assert hits_path.stat().st_mtime_ns == first_mtime

    hits = pq.read_table(hits_path)
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    assert hits.num_rows == 2
    assert manifest["backend"] == "diamond"
    assert manifest["backend_version"] == "diamond version fake-2.1"
    assert manifest["directional"] is True
    assert manifest["hit_count"] == 2
    assert manifest["parameters"]["max_target_seqs"] == 0
    invocations = (tmp_path / "diamond-invocations.txt").read_text(encoding="utf-8").splitlines()
    assert invocations == [
        "version",
        "makedb",
        "makedb",
        "blastp",
        "blastp",
        "blastp",
        "blastp",
        "version",
    ]

    hits_path.write_bytes(b"corrupt")
    assert main(command) == 0
    assert pq.read_table(hits_path).num_rows == 2
    invocations = (tmp_path / "diamond-invocations.txt").read_text(encoding="utf-8").splitlines()
    assert invocations == [
        "version",
        "makedb",
        "makedb",
        "blastp",
        "blastp",
        "blastp",
        "blastp",
        "version",
        "version",
        "makedb",
        "makedb",
        "blastp",
        "blastp",
        "blastp",
        "blastp",
    ]

    edge_command = [
        "edges",
        "--run",
        str(run),
        "--method",
        "lrb",
        # Minimal fixture has one hit per species-pair group; keep v2_max so the
        # group is not dropped (v1_zero default is covered by unit tests).
        "--set",
        "similarity.nbs_fallback=v2_max",
    ]
    assert main(edge_command) == 0
    edge_path = run / "edges" / "retained_edges.parquet"
    edge_manifest_path = run / "edges" / "edge-manifest.json"
    first_edge_mtime = edge_path.stat().st_mtime_ns
    assert main(edge_command) == 0
    assert edge_path.stat().st_mtime_ns == first_edge_mtime
    edges = pq.read_table(edge_path).to_pylist()
    edge_manifest = json.loads(edge_manifest_path.read_text(encoding="utf-8"))
    assert len(edges) == 1
    assert edges[0]["u"] == 0 and edges[0]["v"] == 1
    assert edges[0]["weight"] == pytest.approx(1.0)
    assert edge_manifest["counts"]["coverage_filtered_hits"] == 2
    assert edge_manifest["counts"]["retained_edges"] == 1

    component_command = [
        "components",
        "--run",
        str(run),
        "--set",
        "components.edge_batch_size=1",
    ]
    assert main(component_command) == 0
    index_path = run / "components" / "index.parquet"
    component_manifest_path = run / "components" / "component-manifest.json"
    first_component_mtime = index_path.stat().st_mtime_ns
    assert main(component_command) == 0
    assert index_path.stat().st_mtime_ns == first_component_mtime
    component_manifest = json.loads(component_manifest_path.read_text(encoding="utf-8"))
    assert component_manifest["counts"] == {
        "components": 1,
        "partition_files": 1,
        "proteins": 2,
        "retained_edges": 1,
        "singletons": 0,
    }
    partition = next((run / "components" / "edges").rglob("*.parquet"))
    partition.write_bytes(b"corrupt")
    assert main(component_command) == 0
    assert pq.read_table(partition).num_rows == 1
