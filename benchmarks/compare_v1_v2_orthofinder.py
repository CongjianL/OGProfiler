#!/usr/bin/env python3
"""Three-way artifact comparison: V1 ↔ V2 ↔ OrthoFinder3.

This is the S0 regression harness. It does NOT run the search/clustering tools;
it only compares their already-produced artifacts, remapping every gene back to
its original FASTA identifier before comparing, so that the three different
internal ID schemes (V1 `G0|g0`, V2 integer protein_id, OrthoFinder `species_gene`)
do not leak into the diff.

Two artifact classes are compared:

1. search hits  — (query, target, bitscore) as a multiset;
2. SSN edges    — undirected (u, v) -> weight, with a relative weight tolerance.

Output is a JSON report of pairwise equality, counts, and a bounded list of
first-N differences per pair.

Run (from the repo root):

    python benchmarks/compare_v1_v2_orthofinder.py \
        --v1-working-dir <V1 WorkingDirectory> \
        --v2-run <V2 run root> \
        --of-working-dir <OrthoFinder WorkingDirectory> \
        --out <report.json>
"""

from __future__ import annotations

import argparse
import gzip
import json
from collections import Counter
from pathlib import Path
from typing import Callable, Iterable

import igraph
import pyarrow.parquet as pq


# ---------------------------------------------------------------------------
# ID remapping: every parser returns original-ID-keyed data
# ---------------------------------------------------------------------------


def load_v1_recoded_to_original(seqids_path: Path) -> dict[str, str]:
    """V1 SequenceIDs.txt lines look like ``original_id\\tG0|g0``."""
    mapping: dict[str, str] = {}
    for line in seqids_path.read_text(encoding="utf-8").splitlines():
        parts = line.split("\t")
        if len(parts) >= 2:
            mapping[parts[1].strip()] = parts[0].strip()
    return mapping


def load_v2_id_to_original(proteins_path: Path) -> dict[int, str]:
    rows = pq.read_table(proteins_path, columns=["protein_id", "original_id"]).to_pylist()
    return {int(row["protein_id"]): str(row["original_id"]) for row in rows}


def load_orthofinder_ids(seqids_path: Path) -> tuple[dict[str, str], list[str]]:
    """Parse OrthoFinder SequenceIDs.txt (``OF_ID: original_id``).

    Returns both the ``OF_ID -> original`` map (for Blast rows, whose IDs look
    like ``0_0``) and the originals in file order (for the MCL graph, whose
    nodes are global integer offsets 0..N-1 in SequenceIDs.txt order).

    TODO(remote): confirm the exact delimiter/layout and that line order equals
    global index order from the actual remote run.
    """
    id_to_original: dict[str, str] = {}
    ordered: list[str] = []
    for line in seqids_path.read_text(encoding="utf-8").splitlines():
        if ":" in line:
            of_id, original = line.split(":", 1)
            # OrthoFinder keeps the full FASTA header; V1/V2 keep only the
            # first whitespace token, so normalize to that convention.
            original = original.split()[0]
            id_to_original[of_id.strip()] = original
            ordered.append(original)
    return id_to_original, ordered


# ---------------------------------------------------------------------------
# Search-hit parsers: each returns (query, target, bitscore)
# ---------------------------------------------------------------------------


def _iter_lines(path: Path) -> Iterable[str]:
    if path.suffix == ".gz":
        with gzip.open(path, "rt", encoding="utf-8") as handle:
            yield from handle
    else:
        with path.open(encoding="utf-8") as handle:
            yield from handle


def parse_v1_blast(blast_dir: Path) -> list[tuple[str, str, float]]:
    hits: list[tuple[str, str, float]] = []
    for path in sorted(blast_dir.glob("*.out")):
        for line in _iter_lines(path):
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 12:
                continue
            hits.append((fields[0], fields[1], float(fields[11])))
    return hits


def parse_v2_hits(hits_path: Path) -> list[tuple[int, int, float]]:
    rows = pq.read_table(
        hits_path, columns=["query_id", "target_id", "bitscore"]
    ).to_pylist()
    return [
        (int(row["query_id"]), int(row["target_id"]), float(row["bitscore"]))
        for row in rows
    ]


def parse_orthofinder_blast(workdir: Path) -> list[tuple[str, str, float]]:
    hits: list[tuple[str, str, float]] = []
    for pattern in ("Blast*.txt", "Blast*.txt.gz"):
        for path in sorted(workdir.glob(pattern)):
            for line in _iter_lines(path):
                fields = line.rstrip("\n").split("\t")
                if len(fields) < 12:
                    continue
                hits.append((fields[0], fields[1], float(fields[11])))
    return hits


# ---------------------------------------------------------------------------
# SSN edge parsers: each returns {(u, v): weight} with canonical u < v
# ---------------------------------------------------------------------------


def parse_v1_ssn(gml_path: Path) -> dict[tuple[str, str], float]:
    graph = igraph.Graph.Read_GML(str(gml_path))
    edges: dict[tuple[str, str], float] = {}
    for edge in graph.es:
        left = graph.vs[edge.source]["name"]
        right = graph.vs[edge.target]["name"]
        left, right = sorted((left, right))
        edges[(left, right)] = float(edge["NBS"])
    return edges


def parse_v2_edges(edges_path: Path) -> dict[tuple[int, int], float]:
    rows = pq.read_table(edges_path, columns=["u", "v", "weight"]).to_pylist()
    edges: dict[tuple[int, int], float] = {}
    for row in rows:
        left, right = int(row["u"]), int(row["v"])
        left, right = sorted((left, right))
        edges[(left, right)] = float(row["weight"])
    return edges


def parse_orthofinder_graph(graph_path: Path) -> dict[tuple[int, int], float]:
    """Parse an OrthoFinder MCL matrix graph file.

    Each body line is ``node    neighbor:weight neighbor:weight ... $``.
    Node indices are global integer offsets (0..N-1). Mapping those offsets
    back to original IDs is done by the caller using the SequenceIDs.txt line
    order (TODO(remote): confirm line order equals global index order).
    """
    edges: dict[tuple[int, int], float] = {}
    in_body = False
    for line in _iter_lines(graph_path):
        stripped = line.strip()
        if stripped == "begin":
            in_body = True
            continue
        if not in_body:
            continue
        if stripped in ("", ")", "$"):
            continue
        fields = stripped.split()
        source = int(fields[0])
        for token in fields[1:]:
            if token == "$":
                continue
            target, weight = token.split(":", 1)
            target = int(target)
            left, right = sorted((source, target))
            edges[(left, right)] = float(weight)
    return edges


# ---------------------------------------------------------------------------
# Comparison primitives
# ---------------------------------------------------------------------------


def _remap_hits(
    hits: list[tuple[object, object, float]],
    remap: Callable[[object], str],
) -> Counter[tuple[str, str, int]]:
    """Round bitscore to a stable bucket for exact multiset comparison."""
    counter: Counter[tuple[str, str, int]] = Counter()
    for query, target, score in hits:
        counter[(remap(query), remap(target), round(score, 6))] += 1
    return counter


def compare_multiset(
    left: Counter, right: Counter, label_left: str, label_right: str, limit: int = 20
) -> dict:
    diffs: list[dict] = []
    keys = set(left) | set(right)
    for key in sorted(keys)[: limit]:
        if left[key] != right[key]:
            diffs.append(
                {
                    "key": list(key),
                    f"{label_left}_count": left[key],
                    f"{label_right}_count": right[key],
                }
            )
    return {
        "equal": left == right,
        "left_count": sum(left.values()),
        "right_count": sum(right.values()),
        "unique_left": len(left),
        "unique_right": len(right),
        "diff_count": len(keys),
        "diffs": diffs,
    }


def compare_edges(
    left: dict, right: dict, tolerance: float, limit: int = 20
) -> dict:
    left_keys = set(left)
    right_keys = set(right)
    only_left = sorted(left_keys - right_keys)[:limit]
    only_right = sorted(right_keys - left_keys)[:limit]
    weight_diffs: list[dict] = []
    for key in sorted(left_keys & right_keys):
        delta = abs(left[key] - right[key])
        scale = max(abs(left[key]), abs(right[key]), 1e-12)
        if delta > tolerance * scale:
            weight_diffs.append(
                {"key": list(key), "left_weight": left[key], "right_weight": right[key]}
            )
            if len(weight_diffs) >= limit:
                break
    equal = (
        left_keys == right_keys
        and not only_left
        and not only_right
        and not weight_diffs
    )
    return {
        "equal": equal,
        "left_count": len(left_keys),
        "right_count": len(right_keys),
        "only_left": only_left,
        "only_right": only_right,
        "weight_diff_count": len(weight_diffs),
        "weight_diffs": weight_diffs,
    }


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--v1-working-dir", required=True, type=Path)
    parser.add_argument("--v2-run", required=True, type=Path)
    parser.add_argument("--of-working-dir", required=True, type=Path)
    parser.add_argument("--of-graph", type=Path, default=None,
                        help="Override the OrthoFinder graph filename.")
    parser.add_argument("--out", required=True, type=Path)
    parser.add_argument("--tolerance", type=float, default=1e-6)
    args = parser.parse_args()

    v1 = args.v1_working_dir
    v2 = args.v2_run
    of = args.of_working_dir

    # --- ID maps -----------------------------------------------------------
    v1_remap = load_v1_recoded_to_original(v1 / "SequenceIDs.txt")
    v2_remap = load_v2_id_to_original(v2 / "input" / "proteins.parquet")
    of_remap, of_ordered = load_orthofinder_ids(of / "SequenceIDs.txt")

    # --- search hits -------------------------------------------------------
    v1_hits = parse_v1_blast(v1 / "BlastResults")
    v2_hits = parse_v2_hits(v2 / "search" / "hits.parquet")
    of_hits = parse_orthofinder_blast(of)

    v1_hits_counter = _remap_hits(v1_hits, lambda rid: v1_remap.get(rid, str(rid)))
    v2_hits_counter = _remap_hits(v2_hits, lambda pid: v2_remap.get(pid, str(pid)))
    of_hits_counter = _remap_hits(of_hits, lambda oid: of_remap.get(oid, str(oid)))

    # --- SSN edges ---------------------------------------------------------
    v1_edges_raw = parse_v1_ssn(v1 / "ssn.gml")
    v2_edges_raw = parse_v2_edges(v2 / "edges" / "retained_edges.parquet")
    of_graph_path = args.of_graph
    if of_graph_path is None:
        candidates = [of / "OrthoFinder_graph.txt"] + sorted(of.glob("graph_*.txt"))
        of_graph_path = next((p for p in candidates if p.is_file()), None)
    of_edges_raw = parse_orthofinder_graph(of_graph_path) if of_graph_path else {}

    def _remap_edges(
        edges: dict, remap: Callable[[object], str]
    ) -> dict[tuple[str, str], float]:
        return {
            (remap(u), remap(v)): weight for (u, v), weight in edges.items()
        }

    v1_edges = _remap_edges(v1_edges_raw, lambda rid: v1_remap.get(rid, str(rid)))
    v2_edges = _remap_edges(v2_edges_raw, lambda pid: v2_remap.get(pid, str(pid)))
    of_edges = _remap_edges(
        of_edges_raw,
        # OrthoFinder graph nodes are global ints; map via SequenceIDs order.
        lambda idx: of_ordered[int(idx)]
        if 0 <= int(idx) < len(of_ordered)
        else str(idx),
    )

    report = {
        "search_hits": {
            "v1_vs_v2": compare_multiset(v1_hits_counter, v2_hits_counter, "v1", "v2"),
            "v1_vs_orthofinder": compare_multiset(
                v1_hits_counter, of_hits_counter, "v1", "orthofinder"
            ),
            "v2_vs_orthofinder": compare_multiset(
                v2_hits_counter, of_hits_counter, "v2", "orthofinder"
            ),
        },
        "ssn_edges": {
            "v1_vs_v2": compare_edges(v1_edges, v2_edges, args.tolerance),
            "v1_vs_orthofinder": compare_edges(v1_edges, of_edges, args.tolerance),
            "v2_vs_orthofinder": compare_edges(v2_edges, of_edges, args.tolerance),
        },
        "counts": {
            "v1_hits": len(v1_hits),
            "v2_hits": len(v2_hits),
            "orthofinder_hits": len(of_hits),
            "v1_edges": len(v1_edges_raw),
            "v2_edges": len(v2_edges_raw),
            "orthofinder_edges": len(of_edges_raw),
        },
    }

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(
        json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps(report, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
