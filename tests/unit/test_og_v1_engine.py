from __future__ import annotations

import copy
import json
from pathlib import Path

import pytest
from test_og_v1_events import component_inputs

from ogprofiler.orthogroups.engine import (
    OrthogroupConflictError,
    extract_component_orthogroups,
)

FIXTURES = json.loads(
    (Path(__file__).parents[1] / "fixtures/og_extraction/v1_contract.json").read_text()
)["cases"]


@pytest.mark.parametrize("spec", FIXTURES, ids=lambda spec: spec["name"])
def test_component_og_members_and_remaining_view_match_frozen_v1(spec):
    actual_groups, actual_remaining, unassigned, duplicates = [], [], [], []
    inputs_before = copy.deepcopy(spec)
    for nodes, members, species, ids in component_inputs(spec):
        names = {i: name for name, i in ids.items()}
        gene_ids = {
            g: i for i, g in enumerate(sorted({g for v in spec["vertices"] for g in v["genes"]}))
        }
        proteins = {i: g for g, i in gene_ids.items()}
        originals = {i: g.split("|", 1)[1] for i, g in proteins.items()}
        kwargs = dict(
            total_species=spec["n_species"],
            overlap_count=spec.get("overlap_count", 0),
            ssn_isolates=[gene_ids[g] for g in spec.get("isolates", []) if gene_ids[g] in species],
        )
        snapshot = copy.deepcopy((nodes, members, species, originals))
        if spec["expected"]["duplicate_members"]:
            with pytest.raises(OrthogroupConflictError) as raised:
                extract_component_orthogroups(nodes, members, species, originals, **kwargs)
            result = raised.value.result
            duplicates.extend(proteins[p] for p in raised.value.duplicate_members)
        else:
            result = extract_component_orthogroups(nodes, members, species, originals, **kwargs)
            assert result == extract_component_orthogroups(
                nodes, members, species, originals, **kwargs
            )
        assert (nodes, members, species, originals) == snapshot
        actual_groups.extend(
            (group.processing_level, tuple(proteins[p] for p in group.protein_ids))
            for group in result.groups
        )
        actual_remaining.extend(names[c] for c in result.remaining_cluster_ids)
        unassigned.extend(proteins[row.protein_id] for row in result.unassigned)
        assert all(
            group.n_species == len({species[p] for p in group.protein_ids})
            for group in result.groups
        )
    expected = spec["expected"]
    assert sorted(actual_groups) == sorted(
        (g["level"], tuple(g["members"])) for g in expected["groups"]
    )
    assert sorted(actual_remaining) == sorted(expected["remaining_nodes"])
    assert sorted(unassigned) == expected["unassigned"]
    assert sorted(duplicates) == expected["duplicate_members"]
    assert spec == inputs_before


def test_root_I_merges_multiple_terminal_families_and_traces_consumption():
    spec = FIXTURES[0]
    nodes, members, species, _ = next(component_inputs(spec))
    result = extract_component_orthogroups(
        nodes, members, species, {0: "a", 1: "b"}, total_species=3
    )
    assert len(result.groups) == 1
    assert result.groups[0].source_cluster_id == 0
    assert result.groups[0].selection_type == "EVENT_I"
    assert result.groups[0].protein_ids == (0, 1)
    assert result.remaining_cluster_ids == (0,)
    assert {row.cluster_id for row in result.trace if row.status == "DESCENDANT_CONSUMED"} == {1, 2}


def test_singleton_without_ssn_isolate_evidence_stays_unassigned():
    nodes = [dict(component_id=0, cluster_id=0, parent_id=None, n_genes=1, n_species=1)]
    members = [dict(protein_id=7, terminal_cluster_id=0)]
    result = extract_component_orthogroups(nodes, members, {7: 0}, {7: "same"}, total_species=2)
    assert result.groups == ()
    assert [u.protein_id for u in result.unassigned] == [7]


def test_hash_preserves_species_identity_and_is_independent_of_protein_id():
    nodes = [dict(component_id=0, cluster_id=0, parent_id=None, n_genes=2, n_species=2)]
    first = extract_component_orthogroups(
        nodes,
        [dict(protein_id=p, terminal_cluster_id=0) for p in (0, 1)],
        {0: 0, 1: 1},
        {0: "same", 1: "same"},
        total_species=2,
    )
    second = extract_component_orthogroups(
        nodes,
        [dict(protein_id=p, terminal_cluster_id=0) for p in (50, 9)],
        {50: 0, 9: 1},
        {50: "same", 9: "same"},
        total_species=2,
    )
    assert first.groups[0].membership_hash == second.groups[0].membership_hash


def test_nested_candidate_skip_records_the_consumption_source():
    spec = next(s for s in FIXTURES if s["name"] == "nested_same_coverage_largest_first")
    nodes, members, species, ids = next(component_inputs(spec))
    result = extract_component_orthogroups(
        nodes, members, species, {p: str(p) for p in species}, total_species=2, overlap_count=1
    )
    skipped = [row for row in result.trace if row.status == "SKIPPED_CONSUMED"]
    assert [(row.cluster_id, row.consumed_by) for row in skipped] == [(ids["a"], ids["r"])]


def test_deep_og_selection_uses_static_intervals_not_recursive_gene_caches():
    depth = 1500
    nodes, members = [], []
    for i in range(depth):
        nodes.append(
            dict(
                cluster_id=i,
                parent_id=i - 1 if i else None,
                component_id=0,
                n_genes=depth + 1 - i,
                n_species=1,
            )
        )
    for p in range(depth + 1):
        nodes.append(
            dict(
                cluster_id=depth + p,
                parent_id=min(p, depth - 1),
                component_id=0,
                n_genes=1,
                n_species=1,
            )
        )
        members.append(dict(protein_id=p, terminal_cluster_id=depth + p))
    result = extract_component_orthogroups(
        nodes,
        members,
        dict.fromkeys(range(depth + 1), 0),
        {p: str(p) for p in range(depth + 1)},
        total_species=2,
    )
    assert len(result.groups) == depth + 1
    assert all(
        group.n_genes == 1 and group.selection_type == "RESIDUAL_NONE" for group in result.groups
    )
    assert result.unassigned == ()
