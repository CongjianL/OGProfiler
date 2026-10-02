"""First separation is distinct from OG recovery, through the audit report seam."""

import json
from pathlib import Path

import pyarrow as pa
import pyarrow.parquet as pq
import pytest

from benchmarks.og_extraction.refog_split_audit import audit
from ogprofiler.core.manifest import sha256_file
from ogprofiler.orthogroups.legacy_events import v1_event_for_node


def test_single_copy_k_way_and_binary_have_different_v1_events():
    children = [1 << i for i in range(11)]
    parent = (1 << 11) - 1
    assert v1_event_for_node(parent, children, 11, 12) == "III-3"
    assert v1_event_for_node(parent, [children[0], parent ^ children[0]], 11, 3) == "I"


@pytest.mark.parametrize("rejoin", [False, True])
def test_audit_distinguishes_lost_pairs_and_rejoined_leaves(tmp_path: Path, rejoin: bool):
    run = tmp_path / "new-hierarchy"
    hierarchy = run / "hierarchy/components/component=00000000"
    og = run / "orthogroups/components/component=00000000"
    hierarchy.mkdir(parents=True)
    og.mkdir(parents=True)
    children = [1, 2] if rejoin else [1, 2, 3]
    memberships = [1, 2, 2] if rejoin else [1, 2, 3]
    nodes = [
        dict(
            cluster_id=0, parent_id=None, depth=0, n_genes=3, n_species=3, child_count=len(children)
        )
    ]
    nodes += [
        dict(
            cluster_id=c,
            parent_id=0,
            depth=1,
            n_genes=memberships.count(c),
            n_species=memberships.count(c),
            child_count=0,
        )
        for c in children
    ]

    def write(directory, name, data):
        pq.write_table(pa.Table.from_pylist(data), directory / name)

    write(hierarchy, "nodes.parquet", nodes)
    write(
        hierarchy,
        "members.parquet",
        [dict(protein_id=p, terminal_cluster_id=c) for p, c in enumerate(memberships)],
    )
    write(hierarchy, "candidates.parquet", [dict(cluster_id=0, selected=True, stability=1.0)])
    (hierarchy / "metrics.json").write_text("{}")
    write(
        og,
        "v1_events.parquet",
        [
            dict(
                cluster_id=n["cluster_id"],
                v1_event=("I" if rejoin else "III-3") if n["cluster_id"] == 0 else None,
            )
            for n in nodes
        ],
    )
    write(
        og,
        "selection_trace.parquet",
        [dict(cluster_id=c, status="SELECTED") for c in ([0] if rejoin else children)],
    )
    write(og, "groups.parquet", [dict(local_group_id=0)])
    write(
        og,
        "members.parquet",
        [dict(protein_id=p, local_group_id=0 if rejoin else p) for p in range(3)],
    )
    for directory, name in ((hierarchy, "hierarchy-manifest.json"), (og, "og-manifest.json")):
        checks = {p.name: sha256_file(p) for p in directory.iterdir()}
        (directory / name).write_text(json.dumps(dict(output_checksums=checks)))
    genes = ["A", "B", "C"]
    data = dict(
        protein={g: dict(protein_id=p) for p, g in enumerate(genes)},
        components={p: 0 for p in range(3)},
        selected_components=[0],
        truth={"RefOG001": genes},
        groups={g: "all" if rejoin else g for g in genes},
    )
    (tmp_path / "refog-audit-input.json").write_text(json.dumps(data))
    (tmp_path / "h5-report-compact.json").write_text(
        json.dumps(
            dict(report=dict(results=[{}, {}, {}, dict(official_recall=1.0 if rejoin else 0.0)]))
        )
    )
    report, refs, splits, examples = audit(tmp_path)
    assert report["counts"].get("ssn_lost_pairs", 0) == 0
    assert report["counts"].get("lost_pairs", 0) == (0 if rejoin else 3)
    if rejoin:
        assert report["counts"]["same_terminal_pairs"] == 1
        assert report["counts"]["hierarchy_separated_but_og_rejoined"] == 2
        assert not splits and not examples
    else:
        assert report["counts"]["NO_COMMON_SELECTABLE_ANCESTOR"] == 3
        assert refs[0]["first_split_nodes"] == [dict(component_id=0, cluster_id=0, lost_pairs=3)]
