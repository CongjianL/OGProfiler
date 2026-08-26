"""Graph import, canonical edge tables, and connected components."""

from ogprofiler.graph.components import Component, extract_components
from ogprofiler.graph.edges import EdgeTable, WeightedEdge
from ogprofiler.graph.legacy import import_legacy_ssn
from ogprofiler.graph.partition import (
    load_component_edge_table,
    non_singleton_component_ids,
    read_component_edges,
)

__all__ = [
    "Component",
    "EdgeTable",
    "WeightedEdge",
    "extract_components",
    "import_legacy_ssn",
    "load_component_edge_table",
    "non_singleton_component_ids",
    "read_component_edges",
]
