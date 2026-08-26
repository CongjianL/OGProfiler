from __future__ import annotations

import json

from ogprofiler.evolution.network import annotate_network_events


def _annotate(child_species: list[set[int]], threshold: float = 0.75):
    all_species = set().union(*child_species)
    nodes = [
        {
            "cluster_id": 0,
            "parent_id": None,
            "component_id": 0,
            "depth": 0,
            "n_species": len(all_species),
        }
    ]
    members = []
    species = {}
    protein_id = 0
    for index, values in enumerate(child_species, start=1):
        nodes.append(
            {
                "cluster_id": index,
                "parent_id": 0,
                "component_id": 0,
                "depth": 1,
                "n_species": len(values),
            }
        )
        for species_id in sorted(values):
            members.append({"protein_id": protein_id, "terminal_cluster_id": index})
            species[protein_id] = species_id
            protein_id += 1
    return annotate_network_events(nodes, members, species, threshold)


def test_binary_disjoint_identical_and_partial_overlap() -> None:
    disjoint = _annotate([{0, 1}, {2, 3}])
    assert disjoint[0].network_event == "SPECIATION_LIKE"
    assert disjoint[0].overlap_count == 0
    assert disjoint[0].overlap_score == 0
    assert disjoint[0].legacy_event == "I"

    identical = _annotate([{0, 1}, {0, 1}])
    assert identical[0].network_event == "DUPLICATION_LIKE"
    assert identical[0].overlap_count == 2
    assert identical[0].overlap_score == 1
    assert identical[0].legacy_event == "II"

    partial = _annotate([{0, 1}, {1, 2}])
    assert partial[0].network_event == "MIXED"
    assert partial[0].overlap_count == 1
    assert partial[0].overlap_score == 0.5
    assert partial[0].legacy_event == "III-2"


def test_polytomy_zero_mixed_and_duplication_patterns() -> None:
    zero = _annotate([{0}, {1}, {2}])
    assert zero[0].network_event == "POLYTOMY"
    assert zero[0].legacy_event == "III-3"
    assert len(json.loads(zero[0].pairwise_overlap_summary)) == 3

    mixed = _annotate([{0, 1}, {1, 2}, {3}])
    assert mixed[0].network_event == "MIXED"

    duplication = _annotate([{0, 1}, {0, 1}, {0, 1}])
    assert duplication[0].network_event == "DUPLICATION_LIKE"
    assert duplication[0].confidence == 1


def test_one_species_terminal_family_is_species_specific() -> None:
    annotations = annotate_network_events(
        [{"cluster_id": 0, "parent_id": None, "component_id": 4, "depth": 0, "n_species": 1}],
        [{"protein_id": 9, "terminal_cluster_id": 0}],
        {9: 3},
        0.5,
    )
    assert annotations[0].network_event == "SPECIES_SPECIFIC"
    assert annotations[0].child_count == 0
    assert annotations[0].legacy_event is None
