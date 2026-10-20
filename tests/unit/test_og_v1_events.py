from __future__ import annotations

import json
from collections import deque
from pathlib import Path

import pytest

from benchmarks.og_extraction.reference_v1 import run_reference
from ogprofiler.exceptions import HierarchyError
from ogprofiler.orthogroups.legacy_events import annotate_v1_events, v1_event_for_node

FIXTURES = json.loads(
    (Path(__file__).parents[1] / "fixtures/og_extraction/v1_contract.json").read_text()
)["cases"]
pytestmark = pytest.mark.filterwarnings("ignore:The boolean vertex attribute.*:RuntimeWarning")


def component_inputs(spec):
    """Translate explicit test graph attributes, not production event labels."""
    by_name = {v["name"]: v for v in spec["vertices"]}
    ids = {name: i for i, name in enumerate(by_name)}
    neighbors = {name: [] for name in by_name}
    for a, b in spec["edges"]:
        neighbors[a].append(b)
        neighbors[b].append(a)
    species_ids = {
        s: i for i, s in enumerate(sorted({s for v in by_name.values() for s in v["species"]}))
    }
    protein_ids = {
        g: i for i, g in enumerate(sorted({g for v in by_name.values() for g in v["genes"]}))
    }
    species = {protein_ids[g]: species_ids[g.split("|")[0]] for g in protein_ids}
    seen = set()
    for root in by_name:
        if root in seen:
            continue
        parent = {root: None}
        depth = {root: 0}
        queue = deque([root])
        while queue:
            name = queue.popleft()
            seen.add(name)
            for other in neighbors[name]:
                if other not in parent:
                    parent[other] = name
                    depth[other] = depth[name] + 1
                    queue.append(other)
        nodes, members = [], []
        for name in by_name:
            if name not in parent:
                continue
            vertex = by_name[name]
            nodes.append(
                dict(
                    cluster_id=ids[name],
                    parent_id=ids[parent[name]] if parent[name] else None,
                    component_id=ids[root],
                    depth=depth[name],
                    n_genes=len(vertex["genes"]),
                    n_species=len(vertex["species"]),
                )
            )
            if not any(parent[child] == name for child in parent):
                members.extend(
                    dict(protein_id=protein_ids[g], terminal_cluster_id=ids[name])
                    for g in vertex["genes"]
                )
        yield nodes, members, species, ids


@pytest.mark.parametrize("spec", FIXTURES, ids=lambda spec: spec["name"])
@pytest.mark.parametrize("overlap_count", [0, 1, 10])
def test_component_event_labels_match_actual_v1(spec, overlap_count):
    reference = run_reference(spec, overlap_count)
    for nodes, members, species, ids in component_inputs(spec):
        annotations = annotate_v1_events(nodes, members, species, overlap_count)
        assert [a.reference_order for a in annotations] == list(range(len(nodes)))
        by_id = {a.cluster_id: a for a in annotations}
        for name in ids:
            if ids[name] in by_id:
                annotation = by_id[ids[name]]
                assert annotation.v1_event == reference["raw_events"][name]
                assert annotation.selection_event == (
                    "None"
                    if reference["raw_events"][name] is None
                    else reference["raw_events"][name]
                )


@pytest.mark.parametrize("threshold", [-1, 0.5, True, "0", None])
def test_overlap_is_an_integer_species_count(threshold):
    with pytest.raises(HierarchyError, match="nonnegative integer"):
        v1_event_for_node(3, [1, 2], 2, 2, threshold)


def test_identical_children_remain_II_even_with_large_overlap_tolerance():
    assert v1_event_for_node(3, [3, 3], 4, 2, 10) == "II"
    assert v1_event_for_node(7, [3, 5], 4, 2, 0) == "III-2"
    assert v1_event_for_node(7, [3, 5], 4, 2, 1) == "I"


def test_unary_internal_is_classified_not_silently_annotated():
    with pytest.raises(HierarchyError, match="unsupported legacy shape"):
        v1_event_for_node(3, [1], 3, 2)


def test_corrupt_members_and_node_statistics_are_rejected():
    nodes, members, species, _ = next(component_inputs(FIXTURES[0]))
    with pytest.raises(HierarchyError, match="duplicate terminal membership"):
        annotate_v1_events(nodes, members + [members[0]], species)
    with pytest.raises(HierarchyError, match="Member/species count mismatch"):
        annotate_v1_events(
            [{**row, "n_genes": 99} if row["parent_id"] is None else row for row in nodes],
            members,
            species,
        )
    with pytest.raises(HierarchyError, match="Missing/invalid species"):
        annotate_v1_events(nodes, members, {})
    with pytest.raises(HierarchyError, match="Invalid parent"):
        annotate_v1_events(
            [{**row, "parent_id": 999} if row["parent_id"] is not None else row for row in nodes],
            members,
            species,
        )


def test_deep_hierarchy_is_annotated_iteratively_without_descendant_lists():
    internal_count = 1500
    nodes, members, species = [], [], {}
    for i in range(internal_count):
        nodes.append(
            dict(
                cluster_id=i,
                parent_id=i - 1 if i else None,
                component_id=0,
                depth=i,
                n_genes=internal_count + 1 - i,
                n_species=2,
            )
        )
    for i in range(internal_count + 1):
        parent = min(i, internal_count - 1)
        nodes.append(
            dict(
                cluster_id=internal_count + i,
                parent_id=parent,
                component_id=0,
                depth=parent + 1,
                n_genes=1,
                n_species=1,
            )
        )
        members.append(dict(protein_id=i, terminal_cluster_id=internal_count + i))
        species[i] = i % 2
    assert len(annotate_v1_events(nodes, members, species)) == 3001
