from benchmarks.qfo.msa_path_trace import intervals, scope_profile


def test_scope_separates_eligibility_from_subset_and_expansion():
    assert scope_profile([1, 2], [2, 1])["scope"] == "complete_source"
    assert scope_profile([1, 2], [1])["omitted_from_source"] == 1
    assert scope_profile([1, 2], [1])["scope"] == "subset"
    assert scope_profile([1, 2], [1, 3])["outside_source"] == 1
    assert scope_profile([1, 2], [1, 3])["scope"] == "expanded_or_rearranged"


def test_intervals_represent_exact_original_subtrees():
    nodes = [
        dict(cluster_id=0, parent_id=None),
        dict(cluster_id=1, parent_id=0),
        dict(cluster_id=2, parent_id=0),
    ]
    ids, start, end = intervals(nodes, [(10, 1), (11, 1), (12, 2)])
    assert ids[start[0] : end[0]] == [10, 11, 12]
    assert ids[start[1] : end[1]] == [10, 11]
    assert ids[start[2] : end[2]] == [12]
