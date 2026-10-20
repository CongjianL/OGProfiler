"""Replay frozen MSA decisions and diagnose source-clade versus active-view output."""

from __future__ import annotations

import argparse
import csv
import json
from collections import defaultdict
from pathlib import Path

import pyarrow.parquet as pq

from benchmarks.og_extraction.embleya_reference import read_of
from benchmarks.og_extraction.graph_merge_diagnosis import transitions
from benchmarks.og_extraction.tree_cut import topology
from benchmarks.qfo.msa_qualification import extract
from ogprofiler.core.manifest import sha256_file


def intervals(nodes, membership):
    by_id, children, order = topology(nodes)
    leaves = defaultdict(list)
    for p, c in membership:
        leaves[c].append(p)
    start, end, ids = {}, {}, []
    stack = [(order[0], False)]
    while stack:
        c, done = stack.pop()
        if done:
            end[c] = len(ids)
            continue
        start[c] = len(ids)
        stack.append((c, True))
        if children[c]:
            stack.extend((ch, False) for ch in reversed(children[c]))
        else:
            ids.extend(sorted(leaves[c]))
    return ids, start, end


def scope_profile(source_ids, selected_ids):
    source, selected = set(source_ids), set(selected_ids)
    return dict(
        source_full_size=len(source),
        output_full_size=len(selected),
        inside_source=len(source & selected),
        outside_source=len(selected - source),
        omitted_from_source=len(source - selected),
        scope="expanded_or_rearranged"
        if selected - source
        else "subset"
        if source - selected
        else "complete_source",
    )


def diagnose(run, reference_dir, experiment, out):
    hashes = json.loads((experiment / "input-hashes.json").read_text())

    def frozen(path):
        digest = sha256_file(path)
        if str(path) in hashes and hashes[str(path)] != digest:
            raise ValueError("Input identity changed")
        hashes[str(path)] = digest

    for name in ("summary.json", "members.tsv", "candidates.json", "input-hashes.json"):
        frozen(experiment / name)
    report = json.loads((experiment / "summary.json").read_text())
    if not report["baseline_replay_verified"] or report["failures"] or not report["score_emitted"]:
        raise ValueError("Experiment is not a validated full partition")
    candidates = json.loads((experiment / "candidates.json").read_text())
    decisions = defaultdict(dict)
    for c in candidates:
        decisions[c["component_id"]][c["cluster_id"]] = c

    def read(path):
        frozen(path)
        return pq.ParquetFile(path).read().to_pylist()

    def members(path):
        frozen(path)
        with path.open() as h:
            rows = list(csv.DictReader(h, delimiter="\t"))
        result = {int(r["protein_id"]): r["family_id"] for r in rows}
        if len(result) != len(rows):
            raise ValueError("Duplicate predicted member")
        return result

    actual = members(run / "results/members.tsv")
    prediction = members(experiment / "members.tsv")
    if set(actual) != set(prediction):
        raise ValueError("Partition universe mismatch")
    proteins = read(run / "input/proteins.parquet")
    species = {r["protein_id"]: r["species_id"] for r in proteins}
    original = {r["protein_id"]: r["original_id"] for r in proteins}
    keys = {(r["species_id"], r["original_id"]): r["protein_id"] for r in proteins}
    names = {r["species_name"]: r["species_id"] for r in read(run / "input/species.parquet")}
    frozen(reference_dir / "Orthogroups.tsv")
    reference = {
        keys[k]: v for k, v in read_of(reference_dir / "Orthogroups.tsv", names, keys).items()
    }
    component_ids = defaultdict(list)
    for r in read(run / "components/index.parquet"):
        component_ids[r["component_id"]].append(r["protein_id"])
    expected = defaultdict(set)
    for p, g in prediction.items():
        expected[g].add(p)
    rows, trace_outputs = [], []
    for cid, ids in sorted(component_ids.items()):
        if len(ids) == 1:
            continue
        folder = run / "hierarchy/components" / f"component={cid:08d}"
        nodes = read(folder / "nodes.parquet")
        membership = read(folder / "members.parquet")
        local = [(r["protein_id"], r["terminal_cluster_id"]) for r in membership]
        by_id, children, order = topology(nodes)
        all_ids, start, end = intervals(nodes, local)
        event_path = run / "orthogroups/components" / folder.name / "v1_events.parquet"
        events = {r["cluster_id"]: r["v1_event"] for r in read(event_path)}
        admitted = [c for c, d in decisions[cid].items() if d["qualifies"]]
        result = extract(nodes, membership, species, original, len(names), admitted)
        if {frozenset(g.protein_ids) for g in result.groups} != {
            frozenset(expected[prediction[p]]) for p in ids
        }:
            raise ValueError(f"MSA replay differs from saved partition in {cid}")
        # Also confirm IDs, because archived merge examples use component:local_group.
        if any(
            prediction[p] != f"{cid}:{g.local_group_id}"
            for g in result.groups
            for p in g.protein_ids
        ):
            raise ValueError("Experimental group ID replay differs")
        # No reference labels participate in replay or eligibility decisions.
        selected = {g.source_cluster_id for g in result.groups}
        trace_by_source = defaultdict(list)
        for i, t in enumerate(result.trace):
            record = dict(
                component_id=cid,
                trace_order=i,
                cluster_id=t.cluster_id,
                processing_level=t.processing_level,
                status=t.status,
                consumed_by=t.consumed_by,
                original_event=events.get(t.cluster_id),
                selection_event=t.selection_event,
            )
            if t.consumed_by is not None:
                trace_by_source[t.consumed_by].append(record)
            if t.cluster_id in selected:
                trace_by_source[t.cluster_id].append(record)
        relevant = set()
        for g in result.groups:
            c = g.source_cluster_id
            if c is None:
                continue
            source = all_ids[start[c] : end[c]]
            current = transitions(g.protein_ids, reference, actual)
            original_node = transitions(source, reference, actual)
            if not current["new_tp"] and not current["new_fp"]:
                continue
            relevant.add(c)
            rows.append(
                dict(
                    component_id=cid,
                    experimental_group=f"{cid}:{g.local_group_id}",
                    source_cluster_id=c,
                    original_event=events[c],
                    msa_admitted=c in admitted,
                    qualification_evidence=decisions[cid].get(c),
                    processing_level=g.processing_level,
                    output_transitions=current,
                    complete_source_transitions=original_node,
                    **scope_profile(source, g.protein_ids),
                    path_category=(
                        "admitted_complete_source"
                        if c in admitted
                        else "original_eligible_complete_source"
                    )
                    if set(source) == set(g.protein_ids)
                    else "active_view_scope_changed",
                    trace_records=len(trace_by_source[c]),
                    trace_preview=(trace_by_source[c][:5] + trace_by_source[c][5:][-5:]),
                    consumed_clades=len(
                        {
                            t["cluster_id"]
                            for t in trace_by_source[c]
                            if t["consumed_by"] == c and t["status"] == "DESCENDANT_CONSUMED"
                        }
                    ),
                )
            )
        trace_outputs.extend(
            dict(component_id=cid, source_cluster_id=c, trace=trace_by_source[c])
            for c in sorted(relevant)
        )
    for path, digest in hashes.items():
        if sha256_file(Path(path)) != digest:
            raise ValueError("Frozen input changed")
    bad = sorted(
        [r for r in rows if r["output_transitions"]["new_fp"]],
        key=lambda r: -r["output_transitions"]["new_fp"],
    )
    good = sorted(
        [
            r
            for r in rows
            if r["output_transitions"]["new_tp"] and not r["output_transitions"]["new_fp"]
        ],
        key=lambda r: -r["output_transitions"]["new_tp"],
    )
    strata = defaultdict(lambda: dict(groups=0, new_tp=0, new_fp=0))
    for r in rows:
        for key in (
            r["path_category"],
            f"event:{r['original_event']}",
            f"admitted:{r['msa_admitted']}|scope:{r['scope']}",
        ):
            strata[key]["groups"] += 1
            for k in ("new_tp", "new_fp"):
                strata[key][k] += r["output_transitions"][k]
    totals = {k: sum(r["output_transitions"][k] for r in rows) for k in ("new_tp", "new_fp")}
    if any(totals[k] != report["pair_transitions"][k] for k in totals):
        raise ValueError("Transition totals differ")
    out.mkdir(parents=True, exist_ok=False)
    (out / "groups.json").write_text(json.dumps(rows, indent=2) + "\n")
    (out / "source-traces.json").write_text(json.dumps(trace_outputs, indent=2) + "\n")
    (out / "input-hashes.json").write_text(json.dumps(hashes, indent=2) + "\n")
    (out / "summary.json").write_text(
        json.dumps(
            dict(
                diagnostic_only=True,
                production_changed=False,
                replay_verified=True,
                immutable_inputs_verified=True,
                totals=totals,
                strata=dict(strata),
                top_new_merges=bad[:30],
                top_clean_tp_gains=good[:30],
                interpretation=(
                    "Observed source/scope paths, not causal attribution; "
                    "same archived MSA decisions"
                ),
                experiment=str(experiment),
            ),
            indent=2,
        )
        + "\n"
    )


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for key in ("run", "reference-dir", "experiment", "out"):
        p.add_argument("--" + key, type=Path, required=True)
    a = p.parse_args()
    diagnose(a.run, a.reference_dir, a.experiment, a.out)


if __name__ == "__main__":
    main()
