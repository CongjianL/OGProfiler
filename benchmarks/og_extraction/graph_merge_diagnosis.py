"""Read-only pair-transition and boundary evidence diagnosis; no score fitting."""

from __future__ import annotations

import argparse
import csv
import json
import math
from collections import Counter, defaultdict
from pathlib import Path

import pyarrow.parquet as pq

from benchmarks.og_extraction.embleya_reference import choose2, read_of
from benchmarks.og_extraction.tree_cut import topology
from ogprofiler.core.manifest import sha256_file


def transitions(ids, reference, baseline):
    """Exact unordered pair transitions inside one candidate group, assigned only."""
    ids = [p for p in ids if p in reference]
    ref = Counter(reference[p] for p in ids)
    old = Counter(baseline[p] for p in ids)
    joint = Counter((baseline[p], reference[p]) for p in ids)
    tp = sum(map(choose2, ref.values()))
    shared_tp = sum(map(choose2, joint.values()))
    fp = choose2(len(ids)) - tp
    shared_fp = sum(map(choose2, old.values())) - shared_tp
    return dict(
        assigned_size=len(ids),
        tp=tp,
        fp=fp,
        new_tp=tp - shared_tp,
        new_fp=fp - shared_fp,
        retained_fp=shared_fp,
    )


def relative_strength(branch, edges):
    """Within/cross densities include absent edges as zero, with no pseudocounts.

    Within weight counts each undirected edge once. Per-child incident strength
    counts its internal edges twice and its boundary edges once.
    """
    counts = Counter(branch.values())
    internal, outgoing, between = defaultdict(list), defaultdict(list), defaultdict(list)
    seen = set()
    for u, v, w in edges:
        if u not in branch or v not in branch or u >= v or (u, v) in seen:
            raise ValueError("Invalid or duplicate canonical edge")
        if not math.isfinite(w) or w < 0:
            raise ValueError("Invalid weight")
        seen.add((u, v))
        a, b = branch[u], branch[v]
        if a == b:
            internal[a].append(w)
        else:
            between[tuple(sorted((a, b)))].append(w)
            outgoing[a].append(w)
            outgoing[b].append(w)
    within = {c: math.fsum(internal[c]) for c in counts}
    cross = {c: math.fsum(outgoing[c]) for c in counts}
    wp = sum(choose2(n) for n in counts.values())
    cp = choose2(len(branch)) - wp
    ww = math.fsum(within.values())
    cw = math.fsum(math.fsum(ws) for ws in between.values())
    wd, cd = (ww / wp if wp else None), (cw / cp if cp else None)
    pair_count = choose2(len(counts))
    # Unobserved child pairs contribute zeros, without enumerating quadratic pairs.
    densities = [math.fsum(ws) / (counts[a] * counts[b]) for (a, b), ws in between.items()]
    mean = math.fsum(densities) / pair_count if pair_count else None
    variance = (
        (math.fsum((d - mean) ** 2 for d in densities) + (pair_count - len(densities)) * mean**2)
        / pair_count
        if pair_count
        else None
    )
    return dict(
        within_possible_pairs=wp,
        cross_possible_pairs=cp,
        within_weight=ww,
        cross_weight=cw,
        within_density=wd,
        cross_density=cd,
        cross_to_within_density=cd / wd if cd is not None and wd is not None and wd > 0 else None,
        within_status="no_possible_pairs" if not wp else "zero_weight" if not ww else "positive",
        cross_weight_fraction=cw / (cw + ww) if cw + ww else None,
        child_pair_density_mean=mean,
        child_pair_density_cv=math.sqrt(variance) / mean if mean is not None and mean > 0 else None,
        zero_weight_child_pairs=pair_count - sum(d > 0 for d in densities),
        child_strengths=[
            dict(
                child_id=c,
                size=counts[c],
                within_weight=within[c],
                outgoing_weight=cross[c],
                outgoing_density=cross[c] / (counts[c] * (len(branch) - counts[c]))
                if counts[c] < len(branch)
                else None,
                within_density=within[c] / choose2(counts[c]) if counts[c] > 1 else None,
                outgoing_strength_fraction=cross[c] / (cross[c] + 2 * within[c])
                if cross[c] + 2 * within[c]
                else None,
            )
            for c in sorted(counts)
        ],
    )


def boundary_evidence(branch, edges):
    """Label-free positive cross-child edges. No-edge statistics remain explicit."""
    cross = [(u, v, w) for u, v, w in edges if branch[u] != branch[v] and w > 0]
    weights = sorted((w for _, _, w in cross), reverse=True)
    total = math.fsum(weights)
    incident = defaultdict(list)
    touched = set()
    child_pairs = set()
    for u, v, w in cross:
        incident[u].append(w)
        incident[v].append(w)
        touched.update((u, v))
        child_pairs.add(tuple(sorted((branch[u], branch[v]))))
    counts = Counter(branch.values())
    covered = Counter(branch[p] for p in touched)
    endpoint_shares = [math.fsum(v) / (2 * total) for v in incident.values()] if total else []
    return dict(
        relative_strength=relative_strength(branch, edges),
        children=len(counts),
        boundary_edges=len(cross),
        boundary_weight=total,
        endpoint_coverage=len(touched) / len(branch) if branch else None,
        minimum_child_coverage=min((covered[c] / n for c, n in counts.items()), default=None),
        child_pair_coverage=len(child_pairs) / choose2(len(counts)) if len(counts) > 1 else None,
        largest_edge_share=weights[0] / total if total else None,
        top_decile_edge_share=math.fsum(weights[: max(1, math.ceil(len(weights) / 10))]) / total
        if total
        else None,
        endpoint_weight_hhi=math.fsum(s * s for s in endpoint_shares) if total else None,
        largest_endpoint_weight_share=max(endpoint_shares) if total else None,
    )


def diagnose(run, reference_dir, batch, out):
    out.mkdir(parents=True, exist_ok=False)
    hashes = {}

    def check(path, expected=None):
        value = sha256_file(path)
        if expected is not None and value != expected:
            raise ValueError(f"Checksum mismatch: {path}")
        hashes[path] = value

    check(batch / "report.json")
    report = json.loads((batch / "report.json").read_text())
    if not report["completed"] or report["protocol"]["source_run"] != str(run):
        raise ValueError("Wrong or incomplete batch")
    choices = [c for c in report["candidates"] if c["pair_penalty"] == 0.1]
    if len(choices) != 1:
        raise ValueError("Expected frozen penalty 0.1 candidate")
    candidate = batch / choices[0]["candidate"]
    for path, digest in report["input_hashes"].items():
        check(Path(path), digest)
    check(candidate / "evaluation.json")
    evaluation = json.loads((candidate / "evaluation.json").read_text())
    for path, digest in evaluation["hashes"].items():
        check(Path(path), digest)
    for name in ("Orthogroups.tsv", "Orthogroups_UnassignedGenes.tsv"):
        if reference_dir / name not in hashes:
            raise ValueError("Reference directory differs from frozen evaluation")
    # Prior baseline prediction is outside the graph scorer input hash set.
    check(run / "results/members.tsv")

    def read_members(path):
        rows = list(csv.DictReader(path.open(), delimiter="\t"))
        result = {int(r["protein_id"]): r["family_id"] for r in rows}
        if len(result) != len(rows):
            raise ValueError("Duplicate assignment")
        return result

    prediction = read_members(candidate / "members.tsv")
    baseline = read_members(run / "results/members.tsv")
    proteins = pq.read_table(run / "input/proteins.parquet").to_pylist()
    keys = {(r["species_id"], r["original_id"]): r["protein_id"] for r in proteins}
    names = {
        r["species_name"]: r["species_id"]
        for r in pq.read_table(run / "input/species.parquet").to_pylist()
    }
    reference = {
        keys[k]: v for k, v in read_of(reference_dir / "Orthogroups.tsv", names, keys).items()
    }
    if set(prediction) != set(baseline) or set(prediction) != {r["protein_id"] for r in proteins}:
        raise ValueError("Universe mismatch")
    groups = defaultdict(list)
    for p, g in prediction.items():
        groups[g].append(p)
    rows = {
        g: dict(group=g, full_size=len(ids), **transitions(ids, reference, baseline))
        for g, ids in groups.items()
    }
    selected = defaultdict(dict)
    # All error groups + all zero-FP groups gaining correct pairs: no top-k censoring.
    for g, r in rows.items():
        if r["new_fp"] or (r["fp"] == 0 and r["new_tp"]):
            c, n = g.split(":")
            selected[int(c[1:])][int(n[1:])] = g
    for cid, targets in sorted(selected.items()):
        folder = run / "hierarchy/components" / f"component={cid:08d}"
        if any(folder / name not in hashes for name in ("nodes.parquet", "members.parquet")):
            raise ValueError("Unrecorded hierarchy artifact")
        nodes = pq.ParquetFile(folder / "nodes.parquet").read().to_pylist()
        by_id, children, order = topology(nodes)
        membership = pq.ParquetFile(folder / "members.parquet").read().to_pylist()
        owner = {}
        branch = {}
        edge_groups = defaultdict(list)
        for c in order:
            parent = by_id[c]["parent_id"]
            owner[c] = c if c in targets else owner.get(parent)
            branch[c] = c if parent in targets else branch.get(parent)
        branches = defaultdict(dict)
        for m in membership:
            p, leaf = m["protein_id"], m["terminal_cluster_id"]
            if owner[leaf] is not None:
                g = targets[owner[leaf]]
                if prediction[p] != g:
                    raise ValueError("Prediction is not the indicated complete clade")
                branches[g][p] = branch[leaf] if branch[leaf] is not None else leaf
        for g, b in branches.items():
            if set(b) != set(groups[g]):
                raise ValueError("Incomplete clade")
        for path in sorted((run / "components/edges" / f"component={cid:08d}").glob("*.parquet")):
            if path not in hashes:
                raise ValueError("Unrecorded edge fragment")
            for e in pq.ParquetFile(path).read(columns=["u", "v", "weight"]).to_pylist():
                u, v, w = e["u"], e["v"], e["weight"]
                g = prediction[u]
                if g == prediction[v] and g in branches:
                    edge_groups[g].append((u, v, w))
        for node, g in targets.items():
            rows[g]["structural_leaf"] = not children[node]
            rows[g]["boundary"] = boundary_evidence(branches[g], edge_groups[g])
            # Reference-crossing edge weight is diagnostic only, never a model feature.
            internal = edge_groups[g]
            total = math.fsum(w for _, _, w in internal)
            bad = math.fsum(
                w
                for u, v, w in internal
                if u in reference and v in reference and reference[u] != reference[v]
            )
            rows[g]["reference_discordant_weight_fraction_all_internal"] = (
                bad / total if total else None
            )
    totals = {
        k: sum(r[k] for r in rows.values()) for k in ("tp", "fp", "new_tp", "new_fp", "retained_fp")
    }
    if (
        totals["fp"] != choices[0]["primary"]["discordant_pairs"]
        or totals["tp"] != choices[0]["primary"]["true_positive_pairs"]
    ):
        raise ValueError("Candidate pair totals mismatch")
    oldgroups = defaultdict(list)
    for p in reference:
        oldgroups[baseline[p]].append(p)
    oldtp = sum(
        sum(map(choose2, Counter(reference[p] for p in ids).values())) for ids in oldgroups.values()
    )
    oldfp = sum(choose2(len(ids)) for ids in oldgroups.values()) - oldtp
    totals.update(
        removed_fp=oldfp - totals["retained_fp"], lost_tp=oldtp - (totals["tp"] - totals["new_tp"])
    )
    for path, digest in hashes.items():
        if sha256_file(path) != digest:
            raise ValueError(f"Input changed: {path}")
    result = dict(
        completed=True,
        penalty=0.1,
        diagnostic_only=True,
        score_changed=False,
        totals=totals,
        error_groups=sum(r["new_fp"] > 0 for r in rows.values()),
        clean_gain_controls=sum(r["fp"] == 0 and r["new_tp"] > 0 for r in rows.values()),
        controls_are_not_size_matched=True,
        rows=sorted(rows.values(), key=lambda r: (-r["new_fp"], r["group"])),
        input_hashes={str(p): h for p, h in hashes.items()},
        source_sha256=sha256_file(Path(__file__)),
    )
    (out / "report.json").write_text(json.dumps(result, indent=2, allow_nan=False))
    return result


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for name in ("run", "reference-dir", "batch", "out"):
        p.add_argument("--" + name, type=Path, required=True)
    a = p.parse_args()
    diagnose(a.run, a.reference_dir, a.batch, a.out)


if __name__ == "__main__":
    main()
