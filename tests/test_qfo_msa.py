import random

import pytest

from benchmarks.qfo.msa_qualification import msa_candidates


def tree():
    nodes = [
        dict(
            cluster_id=i,
            parent_id=None if i == 0 else (i - 1) // 2,
            n_genes=8 if i == 0 else 4 if i < 3 else 2,
            n_species=2,
        )
        for i in range(7)
    ]
    return nodes, [(p, 3 + p // 2) for p in range(8)], {i: "III-3" for i in range(7)}


def test_mutual_preference_and_root_exclusion():
    n, m, e = tree()
    rows = msa_candidates(n, m, [(0, 2, 4.0), (1, 3, 4.0), (0, 4, 0.1)], e)
    assert {r["cluster_id"] for r in rows if r["qualifies"]} == {1}
    assert all(r["cluster_id"] != 0 for r in rows)


def test_strict_tie_does_not_qualify():
    n, m, e = tree()
    # D(left,right)=1/4; D(left,external)=2/8=1/4.
    rows = msa_candidates(n, m, [(0, 2, 1.0), (0, 4, 2.0)], e)
    assert not rows[0]["qualifies"]


def test_sparse_weights_against_brute_subtree_density():
    n, m, events = tree()
    rng = random.Random(42)
    descendant = {c: {p for p, leaf in m if leaf == c} for c in range(7)}
    for c in reversed(range(3)):
        descendant[c] = descendant[2 * c + 1] | descendant[2 * c + 2]
    for _ in range(20):
        edges = [
            (u, v, rng.random()) for u in range(8) for v in range(u + 1, 8) if rng.random() < 0.5
        ]
        rows = msa_candidates(n, m, edges, events)

        def density(a, b, edges=edges):
            x, y = descendant[a], descendant[b]
            return sum(w for u, v, w in edges if (u in x and v in y) or (v in x and u in y)) / (
                len(x) * len(y)
            )

        for r in rows:
            c = r["cluster_id"]
            a, b = 2 * c + 1, 2 * c + 2
            external = 2 if c == 1 else 1
            assert r["sibling_density"] == pytest.approx(density(a, b))
            assert r["max_external_child_densities"] == pytest.approx(
                [density(a, external), density(b, external)]
            )
            assert r["qualifies"] == (
                density(a, b) > 0
                and density(a, b) > density(a, external)
                and density(a, b) > density(b, external)
            )


def test_weight_scaling_and_original_eligibility_unchanged():
    n, m, e = tree()
    edges = [(0, 2, 4.0), (0, 4, 0.1)]
    a = msa_candidates(n, m, edges, e)
    b = msa_candidates(n, m, [(u, v, w * 3) for u, v, w in edges], e)
    assert [r["qualifies"] for r in a] == [r["qualifies"] for r in b]
    assert msa_candidates(n, m, edges, {i: "I" for i in range(7)}) == []


def test_benchmark_adapter_restores_shared_engine_and_baseline():
    from benchmarks.qfo.msa_qualification import extract
    from ogprofiler.orthogroups import engine

    n, m, _ = tree()
    for node in n:
        node["component_id"] = 0
    membership = [dict(protein_id=p, terminal_cluster_id=c) for p, c in m]
    species = {p: p % 2 for p, c in m}
    originals = {p: str(p) for p, c in m}
    original_annotation = engine.annotate_v1_events
    baseline = engine.extract_component_orthogroups(
        n, membership, species, originals, total_species=2
    )
    replay = extract(n, membership, species, originals, 2)
    assert replay == baseline
    assert engine.annotate_v1_events is original_annotation
