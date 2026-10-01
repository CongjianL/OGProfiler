from __future__ import annotations

import json
from pathlib import Path

import pytest

from benchmarks.og_extraction.reference_v1 import (
    graph_from_fixture,
    load_reference,
    run_reference,
)

FIXTURES = json.loads(
    (Path(__file__).parents[1] / "fixtures/og_extraction/v1_contract.json").read_text()
)["cases"]
pytestmark = pytest.mark.filterwarnings("ignore:The boolean vertex attribute.*:RuntimeWarning")


@pytest.mark.parametrize("spec", FIXTURES, ids=lambda spec: spec["name"])
def test_actual_v1_reference_matches_frozen_contract(spec):
    assert run_reference(spec) == spec["expected"]


def test_none_and_string_none_have_different_descendant_consumption():
    env = load_reference()
    graph, _ = graph_from_fixture(FIXTURES[0])
    members, deleted = env["GetGenesIDs"](graph.vs[0])
    assert sorted(members) == ["s0|a", "s1|b"]
    assert deleted == []  # Python None makes it an immediate member source
    graph.vs[0]["Event"] = "None"
    members, deleted = env["GetGenesIDs"](graph.vs[0])
    assert sorted(members) == ["s0|a", "s1|b"]
    assert [v["name"] for v in deleted] == ["a", "b"]


def test_missing_global_event_column_is_a_classified_v1_pipeline_error():
    spec = next(s for s in FIXTURES if s["name"] == "degree_zero_duplicate")
    with pytest.raises(KeyError):
        run_reference({**spec, "initialize_event": False})


def test_v1_unary_internal_degree_two_raises_not_a_binary_event():
    spec = dict(
        n_species=3,
        vertices=[
            dict(name="r", genes=["s0|a", "s1|b", "s2|c"], species=["s0", "s1", "s2"]),
            dict(name="u", genes=["s0|a", "s1|b"], species=["s0", "s1"]),
            dict(name="t", genes=["s0|a"], species=["s0"]),
        ],
        edges=[["r", "u"], ["u", "t"]],
    )
    with pytest.raises(IndexError):
        run_reference(spec)


def test_v1_forces_two_species_two_gene_singleton_children():
    env = load_reference()
    graph, ssn = graph_from_fixture(FIXTURES[0])
    # GetAttribution reads ssn names of the form species|gene and vertex index.
    ssn.vs["index"] = ssn.vs.indices
    lines = []
    children = env["HHN"](ssn, graph).RunCommunityDetection(
        lines, [("0", list(ssn.vs.indices))], "rber", "NBS", 1.0
    )
    assert children == [("0-0", [0]), ("0-1", [1])]
    assert lines[0].startswith("0,0-0 0-1,")
