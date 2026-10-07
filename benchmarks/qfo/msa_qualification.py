"""MSA-v1 benchmark-only qualification adapter; production code/defaults untouched."""

from __future__ import annotations

import argparse
import csv
import json
import math
from collections import Counter, defaultdict
from dataclasses import replace
from pathlib import Path
from unittest.mock import patch

import pyarrow.parquet as pq

from benchmarks.og_extraction.embleya_reference import choose2, compare, read_of
from benchmarks.og_extraction.tree_cut import topology
from ogprofiler.core.manifest import sha256_file
from ogprofiler.orthogroups import engine
from ogprofiler.orthogroups.legacy_events import annotate_v1_events


def msa_candidates(nodes, membership, edges, events):
    """Sparse subtree/ancestor-sibling weights in O(edges * tree height). No labels."""
    by_id, children, order = topology(nodes)
    leaf = dict(membership)
    if len(leaf) != len(membership) or any(c not in by_id or children[c] for c in leaf.values()):
        raise ValueError("Invalid membership")
    sizes, weights, depth = Counter(leaf.values()), defaultdict(float), {}
    for c in order:
        parent = by_id[c]["parent_id"]
        depth[c] = 0 if parent is None else depth[parent] + 1
    for c in reversed(order):
        sizes[c] += sum(sizes[ch] for ch in children[c])
        if sizes[c] != by_id[c]["n_genes"]:
            raise ValueError("Node size mismatch")
    seen = set()
    for u, v, w in edges:
        if (
            u not in leaf
            or v not in leaf
            or u >= v
            or (u, v) in seen
            or not math.isfinite(w)
            or w < 0
        ):
            raise ValueError("Invalid canonical edge")
        seen.add((u, v))
        a, b = leaf[u], leaf[v]
        pa, pb = [], []
        while depth[a] > depth[b]:
            pa.append(a)
            a = by_id[a]["parent_id"]
        while depth[b] > depth[a]:
            pb.append(b)
            b = by_id[b]["parent_id"]
        while a != b:
            pa.append(a)
            pb.append(b)
            a, b = by_id[a]["parent_id"], by_id[b]["parent_id"]
        if pa and pb:
            for x in pa:
                weights[x, pb[-1]] += w
            for x in pb:
                weights[x, pa[-1]] += w
    rows = []
    for c in order:
        if len(children[c]) != 2 or by_id[c]["n_species"] <= 1 or by_id[c]["parent_id"] is None:
            continue
        if events[c] not in ("II", "III-1", "III-2", "III-3"):
            continue
        a, b = children[c]
        d = weights[a, b] / (sizes[a] * sizes[b])
        outside = []
        branch = c
        while by_id[branch]["parent_id"] is not None:
            parent = by_id[branch]["parent_id"]
            outside.extend(x for x in children[parent] if x != branch)
            branch = parent
        maxima = [max(weights[ch, x] / (sizes[ch] * sizes[x]) for x in outside) for ch in (a, b)]
        rows.append(
            dict(
                cluster_id=c,
                original_event=events[c],
                sibling_density=d,
                max_external_child_densities=maxima,
                qualifies=d > 0 and all(d > m for m in maxima),
            )
        )
    return rows


def extract(nodes, membership, species, original, total_species, qualified=()):
    # Override only the in-memory annotation view inside the benchmark adapter.
    # Original event tables are not changed. Shared engine sorting/active-view is unchanged.
    admitted = set(qualified)

    def annotations(*args, **kwargs):
        return [
            replace(a, v1_event="I") if a.cluster_id in admitted else a
            for a in annotate_v1_events(*args, **kwargs)
        ]

    with patch.object(engine, "annotate_v1_events", annotations):
        result = engine.extract_component_orthogroups(
            nodes, membership, species, original, total_species=total_species
        )
    if result.unassigned:
        raise ValueError("Experiment extraction has unassigned proteins")
    return result


def diagnose(run, reference_dir, baseline, out):
    prior = json.loads(baseline.read_text())
    frozen = dict(prior["input_hashes"])
    frozen[str(baseline)] = sha256_file(baseline)
    boundary_path = baseline.parent / "boundary-summary.json"
    frozen[str(boundary_path)] = sha256_file(boundary_path)
    boundary = json.loads(boundary_path.read_text())
    if boundary["source_run"] != str(run):
        raise ValueError("Baseline source mismatch")
    frozen[str(run / "edges/retained_edges.parquet")] = boundary["fixed_graph_sha256"]

    def read(path):
        digest = sha256_file(path)
        if str(path) in frozen and frozen[str(path)] != digest:
            raise ValueError("Frozen input mismatch")
        frozen[str(path)] = digest
        return pq.ParquetFile(path).read().to_pylist()

    proteins = read(run / "input/proteins.parquet")
    species = {r["protein_id"]: r["species_id"] for r in proteins}
    original = {r["protein_id"]: r["original_id"] for r in proteins}
    names = {r["species_name"]: r["species_id"] for r in read(run / "input/species.parquet")}
    with (run / "results/members.tsv").open() as handle:
        actual = {
            int(r["protein_id"]): r["family_id"] for r in csv.DictReader(handle, delimiter="\t")
        }
    components = {
        r["protein_id"]: r["component_id"] for r in read(run / "components/index.parquet")
    }
    ids = defaultdict(list)
    for p, c in components.items():
        ids[c].append(p)
    baseline_partition = defaultdict(set)
    for p, g in actual.items():
        baseline_partition[g].add(p)
    predictions, evidence, failures = {}, [], []
    for cid, members in sorted(ids.items()):
        if len(members) == 1:
            if baseline_partition[actual[members[0]]] != {members[0]}:
                raise ValueError("Baseline isolate is not a singleton")
            predictions[members[0]] = f"{cid}:isolate"
            continue
        folder = run / "hierarchy/components" / f"component={cid:08d}"
        nodes = read(folder / "nodes.parquet")
        membership = read(folder / "members.parquet")
        local = [(r["protein_id"], r["terminal_cluster_id"]) for r in membership]
        replay = extract(nodes, membership, species, original, len(names))
        expected = {frozenset(baseline_partition[g]) for g in {actual[p] for p in members}}
        if {frozenset(g.protein_ids) for g in replay.groups} != expected:
            raise ValueError(f"Baseline replay differs in component {cid}")
        events = {a.cluster_id: a.v1_event for a in annotate_v1_events(nodes, membership, species)}
        event_path = run / "orthogroups/components" / folder.name / "v1_events.parquet"
        stored = {r["cluster_id"]: r["v1_event"] for r in read(event_path)}
        if events != stored:
            raise ValueError("Baseline event replay differs")
        paths = sorted((run / "components/edges" / folder.name).glob("*.parquet"))
        if not paths:
            raise ValueError("Missing graph fragments")

        def edges(paths=paths):
            for path in paths:
                for e in read(path):
                    yield e["u"], e["v"], e["weight"]

        candidates = msa_candidates(nodes, local, edges(), events)
        evidence.extend(dict(component_id=cid, **r) for r in candidates)
        try:
            result = extract(
                nodes,
                membership,
                species,
                original,
                len(names),
                [r["cluster_id"] for r in candidates if r["qualifies"]],
            )
        except (engine.OrthogroupConflictError, ValueError) as error:
            failures.append(dict(component_id=cid, error=str(error)))
            continue
        if {p for g in result.groups for p in g.protein_ids} != set(members):
            raise ValueError("Incomplete component partition")
        selected = {t.cluster_id for t in result.trace if t.status == "SELECTED"}
        for row in candidates:
            row["selected"] = row["cluster_id"] in selected
        for g in result.groups:
            for p in g.protein_ids:
                if p in predictions:
                    raise ValueError("Duplicate prediction")
                predictions[p] = f"{cid}:{g.local_group_id}"
    for path, digest in frozen.items():
        if sha256_file(Path(path)) != digest:
            raise ValueError("Frozen inputs changed")
    out.mkdir(parents=True, exist_ok=False)
    (out / "candidates.json").write_text(json.dumps(evidence, indent=2) + "\n")
    (out / "input-hashes.json").write_text(json.dumps(frozen, indent=2) + "\n")
    summary = dict(
        algorithm="mutual-sibling-affinity-v1",
        production_changed=False,
        eligibility_reference_labels_used=False,
        baseline_replay_verified=True,
        immutable_inputs_verified=True,
        failures=failures,
        score_emitted=not failures,
        candidates=len(evidence),
        admitted=sum(r["qualifies"] for r in evidence),
        admitted_selected=sum(r["qualifies"] and r.get("selected", False) for r in evidence),
        admitted_by_event=dict(Counter(r["original_event"] for r in evidence if r["qualifies"])),
    )
    if not failures:
        if set(predictions) != set(actual):
            raise ValueError("Incomplete full partition")
        with (out / "members.tsv").open("w") as handle:
            writer = csv.writer(handle, delimiter="\t")
            writer.writerow(["protein_id", "family_id"])
            writer.writerows(sorted(predictions.items()))
        # Labels are first loaded after candidate decisions and full extraction.
        keys = {(r["species_id"], r["original_id"]): r["protein_id"] for r in proteins}
        reference = {
            keys[k]: v for k, v in read_of(reference_dir / "Orthogroups.tsv", names, keys).items()
        }
        summary["actual"] = prior["actual"]
        summary["experimental"] = compare(reference, predictions, components)
        summary["delta"] = {
            k: summary["experimental"][k] - prior["actual"][k]
            for k in ("true_positive_pairs", "discordant_pairs", "pair_f1", "pair_recall")
        }
        joint = Counter((actual[p], predictions[p], reference[p]) for p in reference)
        shared_tp = sum(map(choose2, joint.values()))
        shared_pairs = sum(
            map(choose2, Counter((actual[p], predictions[p]) for p in reference).values())
        )
        shared_fp = shared_pairs - shared_tp
        summary["pair_transitions"] = dict(
            new_tp=summary["experimental"]["true_positive_pairs"] - shared_tp,
            removed_tp=prior["actual"]["true_positive_pairs"] - shared_tp,
            new_fp=summary["experimental"]["discordant_pairs"] - shared_fp,
            removed_fp=prior["actual"]["discordant_pairs"] - shared_fp,
        )
        for path, digest in frozen.items():
            if sha256_file(Path(path)) != digest:
                raise ValueError("Frozen inputs changed during scoring")
    (out / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for key in ("run", "reference-dir", "baseline", "out"):
        p.add_argument("--" + key, type=Path, required=True)
    a = p.parse_args()
    diagnose(a.run, a.reference_dir, a.baseline, a.out)


if __name__ == "__main__":
    main()
