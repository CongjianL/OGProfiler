"""Frozen OF-reference representation audit; oracle output is never production OG.

Three per-OG match levels plus a globally optimal complete-cut pair F1 oracle.
No Leiden, search, event relabeling, or production output files are written.
"""

from __future__ import annotations

import argparse
import csv
import json
from collections import Counter, defaultdict
from pathlib import Path

import pyarrow.parquet as pq

from benchmarks.og_extraction.embleya_reference import choose2, compare, evaluate, read_of
from benchmarks.og_extraction.tree_cut import pair_f1_oracle, topology
from ogprofiler.core.manifest import sha256_file
from ogprofiler.orthogroups.legacy_events import annotate_v1_events


def node_profiles(nodes, membership, reference, species):
    """Postorder counts only, no internal descendant protein-list duplication."""
    by_id, children, order = topology(nodes)
    counts = {c: Counter() for c in by_id}
    full = Counter()
    seen = set()
    for p, c in membership:
        if p in seen or c not in by_id or children[c]:
            raise ValueError("Duplicate protein or invalid terminal membership")
        seen.add(p)
        full[c] += 1
        if p in reference:
            counts[c][reference[p]] += 1
    for c in reversed(order):
        for u in children[c]:
            counts[c].update(counts[u])
            full[c] += full[u]
        if full[c] != by_id[c]["n_genes"]:
            raise ValueError("Node membership count mismatch")
    events = annotate_v1_events(
        nodes, [dict(protein_id=p, terminal_cluster_id=c) for p, c in membership], species
    )
    eligible = {
        e.cluster_id
        for e in events
        if e.selection_event == "None"
        or (e.selection_event == "I" and by_id[e.cluster_id]["n_species"] > 1)
    }
    if len(seen) == 1:  # explicit isolated component admission, not an event-rule change
        eligible.add(order[0])
    return counts, full, eligible, {e.cluster_id: e.v1_event for e in events}


def strata(reference, prediction, species):
    """Unordered within-species vs between-species co-membership counts."""
    counts = [Counter(), Counter(), Counter()]
    within = [Counter(), Counter(), Counter()]
    for p, r in reference.items():
        labels = (r, prediction[p], (r, prediction[p]))
        for i, label in enumerate(labels):
            counts[i][label] += 1
            within[i][(label, species[p])] += 1
    all_pairs = [sum(map(choose2, c.values())) for c in counts]
    same_pairs = [sum(map(choose2, c.values())) for c in within]
    result = {}
    for name, values in [
        ("same_species", same_pairs),
        ("cross_species", [a - b for a, b in zip(all_pairs, same_pairs, strict=True)]),
    ]:
        rp, pp, tp = values
        result[name] = dict(
            reference_pairs=rp,
            predicted_pairs=pp,
            true_positive_pairs=tp,
            precision=tp / pp if pp else None,
            recall=tp / rp if rp else None,
            pair_f1=2 * tp / (rp + pp) if rp + pp else None,
        )
    return result


def run_audit(run, reference_dir, out):
    out.mkdir(parents=True, exist_ok=False)
    baseline = evaluate(run, reference_dir, out / "actual-comparison.json")
    frozen = {Path(p): h for p, h in baseline["hashes"].items()}
    proteins = pq.read_table(run / "input/proteins.parquet").to_pylist()
    species = {r["protein_id"]: r["species_id"] for r in proteins}
    keys = {(r["species_id"], r["original_id"]): r["protein_id"] for r in proteins}
    names = {
        r["species_name"]: r["species_id"]
        for r in pq.read_table(run / "input/species.parquet").to_pylist()
    }
    ref = {keys[k]: v for k, v in read_of(reference_dir / "Orthogroups.tsv", names, keys).items()}
    sizes = Counter(ref.values())
    index = pq.read_table(run / "components/index.parquet").to_pylist()
    component = {r["protein_id"]: r["component_id"] for r in index}
    if len(component) != len(index) or set(component) != set(species):
        raise ValueError("Component mapping differs from protein universe")
    component_ids = defaultdict(list)
    for p, c in component.items():
        component_ids[c].append(p)
    actual = {
        int(r["protein_id"]): r["family_id"]
        for r in csv.DictReader((run / "results/members.tsv").open(), delimiter="\t")
    }
    if set(actual) != set(species):
        raise ValueError("Actual prediction must cover entire universe")
    best_all, best_eligible, best_actual = {}, {}, {}
    actual_counts = defaultdict(Counter)
    for p, r in ref.items():
        actual_counts[actual[p]][r] += 1

    def update(best, counts, n, identifier, **extra):
        for r, overlap in counts.items():
            value = 2 * overlap / (sizes[r] + n)
            item = dict(f1=value, overlap=overlap, size=n, identifier=identifier, **extra)
            if r not in best or (value, -n, identifier) > (
                best[r]["f1"],
                -best[r]["size"],
                best[r]["identifier"],
            ):
                best[r] = item

    for g, counts in actual_counts.items():
        update(best_actual, counts, sum(counts.values()), g)
    forest, locations = [], []
    for cid, ids in sorted(component_ids.items()):
        folder = run / "hierarchy/components" / f"component={cid:08d}"
        if len(ids) == 1 and not (folder / "nodes.parquet").exists():
            p = ids[0]
            nodes = [
                dict(
                    cluster_id=0,
                    parent_id=None,
                    component_id=cid,
                    depth=0,
                    child_count=0,
                    n_genes=1,
                    n_species=1,
                    split_status="TERMINAL",
                )
            ]
            membership = [(p, 0)]
        else:
            manifest_path = folder / "hierarchy-manifest.json"
            manifest = json.loads(manifest_path.read_text())
            if (
                manifest["hierarchy_status"] != "RESOLVED"
                or not manifest["structural_validation_passed"]
            ):
                raise ValueError(f"Unverified hierarchy {cid}")
            paths = {
                folder / name: digest
                for name, digest in manifest["output_checksums"].items()
                if name in ("nodes.parquet", "members.parquet")
            }
            if len(paths) != 2:
                raise ValueError("Missing required manifest hashes")
            for path, digest in paths.items():
                if sha256_file(path) != digest:
                    raise ValueError(f"Artifact mismatch {path}")
            for relative, digest in manifest["input_checksums"].items():
                path = run / relative
                if path not in frozen:
                    if sha256_file(path) != digest:
                        raise ValueError(f"Input mismatch {path}")
                    frozen[path] = digest
                elif frozen[path] != digest:
                    raise ValueError("Hierarchy input identity mismatch")
            frozen.update(paths)
            frozen[manifest_path] = sha256_file(manifest_path)
            nodes = pq.ParquetFile(folder / "nodes.parquet").read().to_pylist()
            membership = [
                (r["protein_id"], r["terminal_cluster_id"])
                for r in pq.ParquetFile(folder / "members.parquet").read().to_pylist()
            ]
        if {p for p, _ in membership} != set(ids):
            raise ValueError("Hierarchy differs from indexed component")
        counts, full, eligible, events = node_profiles(nodes, membership, ref, species)
        for c, values in counts.items():
            n = sum(values.values())
            extra = dict(
                component_id=cid,
                cluster_id=c,
                event=events[c],
                full_size=full[c],
                eligible=c in eligible,
            )
            update(best_all, values, n, f"{cid}:{c}", **extra)
            if c in eligible:
                update(best_eligible, values, n, f"{cid}:{c}", **extra)
        forest.append(
            dict(
                nodes=[dict(cluster_id=n["cluster_id"], parent_id=n["parent_id"]) for n in nodes],
                tp={c: sum(map(choose2, v.values())) for c, v in counts.items()},
                pp={c: choose2(sum(v.values())) for c, v in counts.items()},
            )
        )
        locations.append((cid, membership))
    rp = sum(map(choose2, sizes.values()))
    oracle = pair_f1_oracle(forest, rp)
    prediction = {}
    selected = []
    for tree, (cid, membership), cut in zip(forest, locations, oracle["cuts"], strict=True):
        by_id, children, order = topology(tree["nodes"])
        cut = set(cut)
        owner = {}
        for c in order:
            owner[c] = c if c in cut else owner.get(by_id[c]["parent_id"])
        for p, leaf in membership:
            if owner[leaf] is None or p in prediction:
                raise ValueError("Invalid non-overlapping complete cut")
            prediction[p] = (cid, owner[leaf])
        selected.extend(dict(component_id=cid, cluster_id=c) for c in sorted(cut))
    if set(prediction) != set(species):
        raise ValueError("Oracle cut lost proteins")
    score = compare(ref, prediction, component)
    if abs(score["pair_f1"] - oracle["pair_f1"]) > 1e-12:
        raise ValueError("Oracle scoring mismatch")
    rows = []
    ref_components = defaultdict(set)
    for p, r in ref.items():
        ref_components[r].add(component[p])
    for r in sorted(sizes):
        a = best_all[r]
        e = best_eligible.get(r, dict(f1=0))
        o = best_actual[r]
        rows.append(
            dict(
                reference_og=r,
                size=sizes[r],
                components=len(ref_components[r]),
                best_any_node=a,
                best_eligible_node=e,
                actual=o,
                eligibility_gap=a["f1"] - e["f1"],
                selection_gap=e["f1"] - o["f1"],
            )
        )
    for path, h in frozen.items():
        if sha256_file(path) != h:
            raise ValueError(f"Frozen input changed: {path}")
    (out / "per-reference-og.json").write_text(json.dumps(rows, indent=2))
    (out / "oracle-cut-DIAGNOSTIC-ONLY.json").write_text(json.dumps(selected, indent=2))

    def mean(field):
        return sum(r[field]["f1"] for r in rows) / len(rows)

    summary = dict(
        completed=True,
        reference_kind="method concordance, not biological truth",
        source_run=str(run),
        reference_proteins=len(ref),
        reference_groups=len(sizes),
        full_proteins=len(species),
        immutable_inputs_verified=True,
        independent_best_nodes_may_overlap=True,
        diagnostic_oracle_uses_reference_labels=True,
        oracle_objective="global pair F1 over non-overlapping complete fixed-tree cuts",
        oracle_is_not_macro_f1_upper_bound=True,
        whole_hierarchy_macro_best_f1=mean("best_any_node"),
        eligible_macro_best_f1=mean("best_eligible_node"),
        actual_macro_best_f1=mean("actual"),
        exact_any_node=sum(r["best_any_node"]["f1"] == 1 for r in rows),
        exact_eligible_node=sum(r["best_eligible_node"]["f1"] == 1 for r in rows),
        positive_eligibility_gap=sum(r["eligibility_gap"] > 1e-12 for r in rows),
        positive_selection_gap=sum(r["selection_gap"] > 1e-12 for r in rows),
        cross_component_reference_ogs=sum(r["components"] > 1 for r in rows),
        actual=baseline["primary_assigned_only"],
        oracle_cut=score,
        oracle_iterations=oracle["iterations"],
        actual_strata=strata(ref, actual, species),
        oracle_strata=strata(ref, prediction, species),
        input_hashes={str(p): h for p, h in frozen.items()},
        audit_source_sha256=sha256_file(Path(__file__)),
    )
    (out / "report.json").write_text(json.dumps(summary, indent=2))


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--run", type=Path, required=True)
    p.add_argument("--reference-dir", type=Path, required=True)
    p.add_argument("--out", type=Path, required=True)
    a = p.parse_args()
    run_audit(a.run, a.reference_dir, a.out)


if __name__ == "__main__":
    main()
