"""Execute selected, unmodified V1 definitions on tiny explicit fixtures."""

from __future__ import annotations

import ast
import hashlib
import os
import tempfile
import time
from pathlib import Path
from types import SimpleNamespace
from typing import Any

import igraph

REFERENCE_PATH = Path(__file__).resolve().parents[2] / "legacy" / "OGProfiler_v1.py"
REFERENCE_SHA256 = "35607fb7cdbd5f8db6c6a3bc3f22df953884b24af2df023f79307bdfcdc4ac3b"


class _QuietProgress:
    def __init__(self, **kwargs: Any):
        pass

    def start(self) -> None:
        pass

    def update(self, value: int) -> None:
        pass

    def finish(self) -> None:
        pass


def load_reference() -> dict[str, Any]:
    source = REFERENCE_PATH.read_bytes()
    digest = hashlib.sha256(source).hexdigest()
    if digest != REFERENCE_SHA256:
        raise ValueError(f"V1 reference changed: {digest}; review and re-freeze the contract")
    tree = ast.parse(source)
    wanted = {"ExtractOGSorted", "GetGenesID", "GetGenesIDs", "GetAttribution"}
    selected = []
    for node in tree.body:
        if isinstance(node, ast.ClassDef) and node.name == "HHN":
            methods = [
                m
                for m in node.body
                if isinstance(m, ast.FunctionDef)
                and m.name
                in {"__init__", "GetEvolutionEvents", "ExtractOG", "RunCommunityDetection"}
            ]
            selected.append(
                ast.ClassDef(
                    name=node.name, bases=node.bases, keywords=[], body=methods, decorator_list=[]
                )
            )
        elif isinstance(node, ast.FunctionDef) and node.name in wanted:
            selected.append(node)
    env: dict[str, Any] = dict(
        __name__="frozen_v1_og_reference",
        igraph=igraph,
        os=os,
        time=time,
        progressbar=SimpleNamespace(
            ProgressBar=_QuietProgress,
            Percentage=lambda: None,
            Bar=lambda value: None,
            Timer=lambda: None,
        ),
    )
    module = ast.fix_missing_locations(ast.Module(body=selected, type_ignores=[]))
    exec(compile(module, str(REFERENCE_PATH), "exec"), env)
    return env


def graph_from_fixture(spec: dict[str, Any]) -> tuple[igraph.Graph, igraph.Graph]:
    vertices = spec["vertices"]
    graph = igraph.Graph(directed=False)
    graph.add_vertices([v["name"] for v in vertices])
    graph.vs["geneIDs"] = [" ".join(v["genes"]) for v in vertices]
    graph.vs["genomeIDs"] = [" ".join(v["species"]) for v in vertices]
    graph.vs["genesNum"] = [len(v["genes"]) for v in vertices]
    graph.vs["genomesNum"] = [len(v["species"]) for v in vertices]
    # A component-local view must initialize Event, as another component's
    # assignment normally creates this column in V1's global igraph.
    if spec.get("initialize_event", True):
        graph.vs["Event"] = [None] * len(vertices)
    graph.add_edges(spec.get("edges", []))
    genes = sorted({g for v in vertices for g in v["genes"]} | set(spec.get("isolates", [])))
    ssn = igraph.Graph(directed=False)
    ssn.add_vertices(genes)
    nonisolates = [g for g in genes if g not in spec.get("isolates", [])]
    ssn.add_edges(list(zip(nonisolates, nonisolates[1:], strict=False)))
    return graph, ssn


def run_reference(spec: dict[str, Any], overlap_count: int | None = None) -> dict[str, Any]:
    if overlap_count is None:
        overlap_count = spec.get("overlap_count", 0)
    env = load_reference()
    graph, ssn = graph_from_fixture(spec)
    obj = env["HHN"](ssn, graph)
    with tempfile.TemporaryDirectory(prefix="v1-og-reference-") as directory:
        obj.GetEvolutionEvents(overlap_count, directory)
    raw_events = dict(zip(graph.vs["name"], graph.vs["Event"], strict=True))
    # Actual unrefined hnn_analysis line, not an inferred terminal policy.
    graph.vs["Event"] = ["None" if event is None else str(event) for event in graph.vs["Event"]]
    remaining, groups = env["ExtractOGSorted"](
        graph, ssn, {"GenomeToUsed": list(range(spec["n_species"]))}, "I", 0, 0
    )
    records = [
        dict(
            source=name,
            level=value[0],
            reported_count=value[1],
            members=sorted(value[2].split(" ")),
            consumed=value[3],
        )
        for name, value in groups.items()
    ]
    all_genes = sorted(
        {g for v in spec["vertices"] for g in v["genes"]} | set(spec.get("isolates", []))
    )
    counts = {gene: sum(gene in row["members"] for row in records) for gene in all_genes}
    return dict(
        reference_sha256=REFERENCE_SHA256,
        raw_events=raw_events,
        groups=records,
        remaining_nodes=remaining.vs["name"],
        remaining_degrees=dict(zip(remaining.vs["name"], remaining.degree(), strict=True)),
        unassigned=[g for g, n in counts.items() if not n],
        duplicate_members=[g for g, n in counts.items() if n > 1],
    )
