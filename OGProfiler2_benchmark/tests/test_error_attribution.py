from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).parents[1]/"scripts/diagnostics"))
from orthobench_pairwise_attribution import attribute, official_family_contribution


def test_perfect_refog_pairwise_combinatorics():
    tp,fp,fn,n=official_family_contribution({"a","b","c"},set(),{"a","b","c"})
    assert (tp,fp,fn,n)==(1.5,0.0,0.0,3)


def test_low_certainty_member_is_removed_from_reference_and_prediction():
    tp,fp,fn,n=official_family_contribution({"a","b","c"},{"c"},{"a","b","c"})
    assert (tp,fp,fn,n)==(1.0,0.0,0.0,2)


def test_disconnected_refog_has_family_attributable_fn():
    tp,fp,fn,n=official_family_contribution({"a","b","c","d"},set(),{"a","b"})
    assert (tp,fp,fn,n)==(1/3,0.0,2/3,2)


def test_hierarchy_split_is_same_pairwise_partition_shape():
    truth={"R":{"a","b","c","d"}};low={"R":set()}
    rows=attribute({"P1":{"a","b"},"P2":{"c","d"}},truth,low)
    assert sum(float(r["official_like_TP"]) for r in rows)==2/3
    assert sum(float(r["official_like_FN_if_applicable"]) for r in rows)==4/3


def test_two_refogs_catastrophically_fused():
    truth={"R1":{"a","b"},"R2":{"c","d"}};low={"R1":set(),"R2":set()}
    row=attribute({"P":{"a","b","c","d"}},truth,low)[0]
    assert row["official_like_TP"]==2.0
    assert row["official_like_FP"]==8.0
    assert row["n_distinct_refogs"]==2


def test_bridge_protein_fixture_connects_two_groups():
    edges={("a","x"),("x","c")}
    adjacency={v:set() for e in edges for v in e}
    for a,b in edges: adjacency[a].add(b);adjacency[b].add(a)
    assert len(adjacency["x"])==2 and len(adjacency["a"]|adjacency["c"])==1


def test_one_giant_family_dominates_fp():
    truth={"R":{"a","b"}};low={"R":set()}
    rows=attribute({"giant":{"a","b",*{f"x{i}" for i in range(100)}},"clean":{"x200"}},truth,low)
    assert rows[0]["family_id"]=="giant"
    assert rows[0]["fraction_of_all_FP"]==1.0
