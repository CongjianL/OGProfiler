from __future__ import annotations

import tracemalloc
from collections.abc import Iterator
from pathlib import Path

import pyarrow as pa
import pyarrow.parquet as pq

from ogprofiler.orthology.engine import OrthologCandidate, generate_component_candidates
from ogprofiler.orthology.stage import HEADER, run_orthology_stage, write_candidate_stream


def _node(cluster_id: int, parent_id: int | None) -> dict[str, int | None]:
    return {"cluster_id": cluster_id, "parent_id": parent_id}


def test_cross_child_pairs_filter_same_species_and_retain_coorthologs() -> None:
    nodes = [_node(0, None), _node(1, 0), _node(2, 0)]
    memberships = [
        {"protein_id": 0, "terminal_cluster_id": 1},
        {"protein_id": 1, "terminal_cluster_id": 1},
        {"protein_id": 2, "terminal_cluster_id": 2},
        {"protein_id": 3, "terminal_cluster_id": 2},
    ]
    events = [
        {"cluster_id": 0, "network_event": "SPECIATION_LIKE"},
        {"cluster_id": 1, "network_event": "AMBIGUOUS"},
        {"cluster_id": 2, "network_event": "AMBIGUOUS"},
    ]
    pairs = list(
        generate_component_candidates(
            component_id=7,
            nodes=nodes,
            terminal_memberships=memberships,
            events=events,
            species_by_protein={0: 0, 1: 0, 2: 1, 3: 0},
        )
    )
    assert [(pair.protein_a_id, pair.protein_b_id) for pair in pairs] == [(0, 2), (1, 2)]
    assert {pair.supporting_cluster_id for pair in pairs} == {0}
    assert {pair.relationship for pair in pairs} == {"CO_ORTHOLOG_CANDIDATE"}


def test_nested_speciation_nodes_emit_each_lca_pair_once() -> None:
    nodes = [
        _node(0, None),
        _node(1, 0),
        _node(2, 0),
        _node(3, 1),
        _node(4, 1),
    ]
    memberships = [
        {"protein_id": 0, "terminal_cluster_id": 3},
        {"protein_id": 1, "terminal_cluster_id": 4},
        {"protein_id": 2, "terminal_cluster_id": 2},
    ]
    events = [
        {"cluster_id": 0, "network_event": "SPECIATION_LIKE"},
        {"cluster_id": 1, "network_event": "SPECIATION_LIKE"},
        {"cluster_id": 2, "network_event": "AMBIGUOUS"},
        {"cluster_id": 3, "network_event": "SPECIES_SPECIFIC"},
        {"cluster_id": 4, "network_event": "SPECIES_SPECIFIC"},
    ]
    pairs = list(
        generate_component_candidates(
            component_id=0,
            nodes=nodes,
            terminal_memberships=memberships,
            events=events,
            species_by_protein={0: 0, 1: 1, 2: 2},
        )
    )
    assert [(p.protein_a_id, p.protein_b_id, p.supporting_cluster_id) for p in pairs] == [
        (0, 2, 0),
        (1, 2, 0),
        (0, 1, 1),
    ]


def test_disjoint_polytomy_matches_simple_v1_cross_species_pairs() -> None:
    nodes = [_node(0, None), *[_node(index, 0) for index in range(1, 5)]]
    memberships = [
        {"protein_id": index - 1, "terminal_cluster_id": index}
        for index in range(1, 5)
    ]
    events = [
        {"cluster_id": 0, "network_event": "POLYTOMY"},
        *[
            {"cluster_id": index, "network_event": "SPECIES_SPECIFIC"}
            for index in range(1, 5)
        ],
    ]
    pairs = list(
        generate_component_candidates(
            component_id=0,
            nodes=nodes,
            terminal_memberships=memberships,
            events=events,
            species_by_protein={index: index for index in range(4)},
        )
    )
    assert [(pair.protein_a_id, pair.protein_b_id) for pair in pairs] == [
        (0, 1),
        (0, 2),
        (0, 3),
        (1, 2),
        (1, 3),
        (2, 3),
    ]


def test_streaming_zstd_writer_has_chunk_bounded_memory(tmp_path: Path) -> None:
    pair_count = 100_000

    def candidates() -> Iterator[OrthologCandidate]:
        for index in range(pair_count):
            yield OrthologCandidate(index, index + pair_count, 0, 1, 0, 0)

    output = tmp_path / "pairs.tsv.zst"
    tracemalloc.start()
    count, chunks = write_candidate_stream(
        output,
        candidates(),
        chunk_size=1_000,
    )
    _, peak = tracemalloc.get_traced_memory()
    tracemalloc.stop()
    assert count == pair_count
    assert chunks == 100
    assert peak < 2_000_000
    with pa.input_stream(str(output), compression="zstd") as stream:
        text = stream.read().decode("utf-8")
    assert text.startswith(HEADER)
    assert text.count("\n") == pair_count + 1


def test_all_singleton_run_produces_header_only_output(tmp_path: Path) -> None:
    (tmp_path / "input").mkdir()
    (tmp_path / "components").mkdir()
    pq.write_table(
        pa.table({"protein_id": [0, 1], "species_id": [0, 1]}),
        tmp_path / "input/proteins.parquet",
    )
    pq.write_table(
        pa.table(
            {
                "protein_id": [0, 1],
                "component_id": [0, 1],
                "terminal_reason": ["SINGLETON", "SINGLETON"],
            }
        ),
        tmp_path / "components/singleton_terminal_families.parquet",
    )
    manifest, reused, pairs = run_orthology_stage(
        tmp_path, ["ogprofiler", "orthologs"], chunk_size=10
    )
    assert not reused
    assert pairs == 0
    assert manifest.is_file()
    with pa.input_stream(
        str(tmp_path / "results/ortholog_pairs.tsv.zst"), compression="zstd"
    ) as stream:
        assert stream.read().decode() == HEADER
