import pytest

from benchmarks.qfo.fixed_tree import boundary_losses
from benchmarks.qfo.pure_lca_trace import classify


def test_pure_node_losses_preserve_exact_pair_accounting():
    nodes = [
        dict(cluster_id=0, parent_id=None),
        dict(cluster_id=1, parent_id=0),
        dict(cluster_id=2, parent_id=0),
    ]
    detail = {}
    result = boundary_losses(
        nodes, [(1, 1), (2, 2), (3, 2)], {1: "A", 2: "A", 3: "A"}, {1: "X", 2: "Y", 3: "Y"}, detail
    )
    assert detail == {
        0: dict(lost_pairs=2, reference_og="A", assigned_size=3, production_fragments=2)
    }
    assert (
        sum(r["lost_pairs"] for r in detail.values()) == result["pure_assigned_lca_tree_available"]
    )


@pytest.mark.parametrize(
    "event,species,status,expected",
    [
        ("II", 2, "DESCENDANT_CONSUMED", "excluded_by_event_qualification"),
        ("I", 1, None, "excluded_by_event_qualification"),
        ("I", 2, "SKIPPED_CONSUMED", "eligible_consumed_without_selection"),
        ("I", 2, "DESCENDANT_CONSUMED", "eligible_consumed_without_selection"),
        (None, 1, "SELECTED", "eligible_selected_but_reference_fragmented"),
        ("I", 2, None, "eligible_without_observed_selection"),
    ],
)
def test_observed_paths(event, species, status, expected):
    trace = [] if status is None else [dict(status=status)]
    assert classify(dict(n_species=species), event, trace) == expected


def test_selected_path_takes_precedence_over_later_consumption():
    assert (
        classify(
            dict(n_species=2), "I", [dict(status="SELECTED"), dict(status="DESCENDANT_CONSUMED")]
        )
        == "eligible_selected_but_reference_fragmented"
    )


def test_selected_ineligible_is_inconsistent():
    with pytest.raises(ValueError, match="Ineligible"):
        classify(dict(n_species=2), "III-1", [dict(status="SELECTED")])
