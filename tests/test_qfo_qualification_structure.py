import pytest

from benchmarks.qfo.qualification_structure import structural_features


def fixture():
    nodes = [
        dict(cluster_id=i, parent_id=None if i == 0 else 0, n_genes=4 if i == 0 else 2)
        for i in range(3)
    ]
    return nodes, [(0, 1), (1, 1), (2, 2), (3, 2)], {0: 0, 1: 1, 2: 0, 3: 1}


def test_cross_relative_to_internal_uses_missing_pairs_as_zero():
    n, m, s = fixture()
    f = structural_features(n, m, s, [(0, 1, 4.0), (2, 3, 4.0), (0, 2, 2.0)])[0]
    assert f["within_density"] == 4
    assert f["cross_density"] == 0.5
    assert f["cross_to_within_density"] == 0.125
    assert f["species_intersection_fraction"] == 1
    assert f["external_weight_fraction"] == 0


def test_zero_weight_not_arbitrary_ratio():
    n, m, s = fixture()
    f = structural_features(n, m, s, [(0, 2, 2.0)])[0]
    assert f["cross_to_within_density"] is None
    assert f["density_status"] == "zero_within_weight"


def test_uniform_weight_rescaling_preserves_dimensionless_features():
    n, m, s = fixture()
    e = [(0, 1, 4.0), (2, 3, 4.0), (0, 2, 2.0)]
    a = structural_features(n, m, s, e)[0]
    b = structural_features(n, m, s, [(u, v, w * 7) for u, v, w in e])[0]
    for k in ("cross_to_within_density", "external_weight_fraction", "min_child_fraction"):
        assert a[k] == pytest.approx(b[k])


def test_duplicate_edges_fail():
    n, m, s = fixture()
    with pytest.raises(ValueError, match="canonical"):
        structural_features(n, m, s, [(0, 1, 1.0), (0, 1, 1.0)])
