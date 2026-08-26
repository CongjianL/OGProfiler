"""Importer for OGProfiler 1 GML sequence-similarity networks."""

from __future__ import annotations

from pathlib import Path
from typing import Any

import igraph as ig

from ogprofiler.exceptions import InputError
from ogprofiler.graph.edges import EdgeTable


def _integer_value(value: Any) -> int | None:
    if isinstance(value, bool):
        return None
    if isinstance(value, int):
        return value
    if isinstance(value, float) and value.is_integer():
        return int(value)
    if isinstance(value, str):
        try:
            parsed = int(value)
        except ValueError:
            return None
        if str(parsed) == value.strip() or value.strip().lstrip("+").isdigit():
            return parsed
    return None


def _vertex_labels(graph: ig.Graph, requested_attribute: str | None) -> list[str]:
    attributes = set(graph.vs.attributes())
    candidates = [requested_attribute, "original_id", "name", "id"]
    attribute = next(
        (candidate for candidate in candidates if candidate and candidate in attributes), None
    )
    if attribute is None:
        return [str(index) for index in range(graph.vcount())]
    labels = [str(value) for value in graph.vs[attribute]]
    if len(set(labels)) != len(labels):
        raise InputError(f"Legacy SSN vertex attribute {attribute!r} is not unique")
    return labels


def _assign_integer_ids(labels: list[str]) -> tuple[list[int], dict[int, str]]:
    parsed = [_integer_value(label) for label in labels]
    if all(value is not None for value in parsed) and len(set(parsed)) == len(parsed):
        integer_ids = [int(value) for value in parsed if value is not None]
        return integer_ids, dict(zip(integer_ids, labels, strict=True))
    sorted_labels = sorted(labels, key=lambda value: value.encode("utf-8"))
    integer_by_label = {label: index for index, label in enumerate(sorted_labels)}
    integer_ids = [integer_by_label[label] for label in labels]
    return integer_ids, {identifier: label for label, identifier in integer_by_label.items()}


def import_legacy_ssn(
    path: Path,
    *,
    weight_attribute: str = "NBS",
    vertex_id_attribute: str | None = None,
) -> EdgeTable:
    if not path.is_file():
        raise InputError(f"Legacy SSN does not exist: {path}")
    try:
        graph = ig.Graph.Read_GML(str(path))
    except (OSError, ig.InternalError) as error:
        raise InputError(f"Failed to read legacy GML {path}: {error}") from error
    if graph.is_directed():
        raise InputError("Legacy SSN importer expects an undirected graph")
    if weight_attribute not in graph.es.attributes() and graph.ecount() > 0:
        available = ", ".join(sorted(graph.es.attributes())) or "none"
        raise InputError(
            f"Edge weight attribute {weight_attribute!r} is missing; available: {available}"
        )

    labels = _vertex_labels(graph, vertex_id_attribute)
    integer_ids, originals = _assign_integer_ids(labels)
    edges: list[tuple[int, int, float]] = []
    for edge in graph.es:
        source_index, target_index = edge.tuple
        weight = float(edge[weight_attribute])
        edges.append((integer_ids[source_index], integer_ids[target_index], weight))
    return EdgeTable.canonicalize(integer_ids, edges, originals)
