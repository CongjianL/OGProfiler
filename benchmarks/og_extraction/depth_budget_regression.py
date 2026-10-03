"""ADR0005 H4-only control and immutable-prefix acceptance (not an OG decision)."""

from __future__ import annotations

import argparse
import copy
import json
from collections import defaultdict
from pathlib import Path

import pyarrow.parquet as pq
import yaml

from benchmarks.og_extraction.depth_path_audit import _tree
from ogprofiler.config import validate_config
from ogprofiler.core.manifest import sha256_file


def depth_config(before):
    after = copy.deepcopy(before)
    assert after["hierarchy"]["topology_policy"] == "soft_binary_24_v2"
    assert after["hierarchy"]["max_depth"] == 20
    after["hierarchy"]["max_depth"] = 42
    validate_config(after)
    return after


def compare_prefix(old_rows, old_members, old_candidates, rows, members, candidates, *, depth=20):
    old, _, old_desc, old_terminal = _tree(old_rows, old_members)
    new, _, new_desc, new_terminal = _tree(rows, members)
    same = {frozenset(ps): cid for cid, ps in new_desc.items()}
    old_children, new_children = defaultdict(list), defaultdict(list)
    for source, index in ((old_rows, old_children), (rows, new_children)):
        for n in source:
            if n["parent_id"] is not None:
                index[n["parent_id"]].append(n["cluster_id"])
    traces_old, traces_new = defaultdict(list), defaultdict(list)
    for source, index in ((old_candidates, traces_old), (candidates, traces_new)):
        for row in source:
            index[row["cluster_id"]].append({k: v for k, v in row.items() if k != "cluster_id"})
    mismatches = []
    count = 0
    for cid, n in old.items():
        if n["depth"] >= depth:
            continue
        count += 1
        nid = same.get(frozenset(old_desc[cid]))
        if nid is None:
            mismatches.append(dict(baseline_cluster_id=cid, reason="MISSING_CLADE"))
            continue
        ignore = {"cluster_id", "parent_id"}
        checks = dict(
            node_fields={k: v for k, v in n.items() if k not in ignore}
            == {k: v for k, v in new[nid].items() if k not in ignore},
            partition={frozenset(old_desc[c]) for c in old_children[cid]}
            == {frozenset(new_desc[c]) for c in new_children[nid]},
            candidates=traces_old[cid] == traces_new[nid],
        )
        if not all(checks.values()):
            mismatches.append(dict(baseline_cluster_id=cid, new_cluster_id=nid, checks=checks))
    return dict(
        passed=not mismatches and set(old_terminal) == set(new_terminal),
        checked_nodes=count,
        cutoff_depth=depth,
        mismatches=mismatches,
    )


def main():
    p = argparse.ArgumentParser()
    p.add_argument("command", choices=("configure", "validate"))
    p.add_argument("--baseline", type=Path, required=True)
    p.add_argument("--current", type=Path, required=True)
    p.add_argument("--out", type=Path, required=True)
    args = p.parse_args()
    before = yaml.safe_load((args.baseline / "new-hierarchy/run.yaml").read_text())
    after = depth_config(before)
    if args.command == "configure":
        (args.current / "new-hierarchy/run.yaml").write_text(yaml.safe_dump(after, sort_keys=True))
        args.out.write_text(
            json.dumps(
                dict(
                    baseline=str(args.baseline),
                    previous=before,
                    current=after,
                    changed_parameter="hierarchy.max_depth",
                    previous_value=20,
                    current_value=42,
                    scope="H4_ONLY",
                    h5_enabled=False,
                ),
                indent=2,
            )
        )
        return
    assert yaml.safe_load((args.current / "new-hierarchy/run.yaml").read_text()) == after

    def load(root):
        folder = root / "new-hierarchy/hierarchy/components/component=00000000"
        manifest = json.loads((folder / "hierarchy-manifest.json").read_text())
        for name, digest in manifest["output_checksums"].items():
            assert sha256_file(folder / name) == digest
        nodes = pq.ParquetFile(folder / "nodes.parquet").read().to_pylist()
        members = pq.ParquetFile(folder / "members.parquet").read().to_pylist()
        return (
            nodes,
            [(r["protein_id"], r["terminal_cluster_id"]) for r in members],
            pq.ParquetFile(folder / "candidates.parquet").read().to_pylist(),
        ), manifest

    old, old_manifest = load(args.baseline)
    current, manifest = load(args.current)
    old_effective = copy.deepcopy(old_manifest["parameters"]["hierarchy"])
    old_effective["max_depth"] = 42
    assert old_effective == manifest["parameters"]["hierarchy"]
    prefix = compare_prefix(*old, *current)
    old_metrics = old_manifest["metrics"]
    metrics = manifest["metrics"]
    unresolved = [n for n in current[0] if n["split_status"] == "UNRESOLVED"]
    extra_calls = metrics["leiden_calls"] - old_metrics["leiden_calls"]
    extra_nodes = len(current[0]) - len(old[0])
    checks = dict(
        prefix_unchanged=prefix["passed"],
        extra_calls_bound=0 <= extra_calls <= 30744,
        extra_nodes_bound=0 <= extra_nodes <= 854,
        no_depth_limit=not any(n.get("terminal_reason") == "DEPTH_LIMIT" for n in unresolved),
    )
    report = dict(
        passed=all(checks.values()),
        checks=checks,
        prefix=prefix,
        extra_calls=extra_calls,
        extra_nodes=extra_nodes,
        unresolved=unresolved,
        resolved=not unresolved,
        baseline_manifest=old_manifest,
        current_manifest=manifest,
    )
    args.out.write_text(json.dumps(report, indent=2))
    if not report["passed"]:
        raise RuntimeError("Depth control acceptance failed")


if __name__ == "__main__":
    main()
