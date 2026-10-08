"""Frozen-cut species and strength-null gain decomposition, diagnostic only."""

from __future__ import annotations

import argparse
import json
import math
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np
import pyarrow.parquet as pq

from ogprofiler.core.manifest import sha256_file


def decompose(prediction, species, edges):
    """No reference labels. Expected cut-edge W_ij=s_i*s_j/(2W)."""
    groups = sorted(set(prediction.values()))
    genes = Counter(prediction.values())
    composition = Counter((prediction[p], species[p]) for p in prediction)
    strength = Counter()
    species_strength = Counter()
    observed = Counter()
    pair_observed = Counter()
    seen = set()
    within_species_weight = Counter()
    within_cut_weight = Counter()
    total = 0.0
    for u, v, w in edges:
        if (
            u not in prediction
            or v not in prediction
            or u >= v
            or (u, v) in seen
            or not math.isfinite(w)
            or w < 0
        ):
            raise ValueError("Invalid edge")
        seen.add((u, v))
        total += w
        a, b = prediction[u], prediction[v]
        strength[a] += w
        strength[b] += w
        species_strength[a, species[u]] += w
        species_strength[b, species[v]] += w
        kind = "same_species" if species[u] == species[v] else "cross_species"
        within_species_weight[kind] += w
        if a == b:
            within_cut_weight[kind] += w
        else:
            observed[kind] += w
            pair_observed[min(a, b), max(a, b), kind] += w
    expected = Counter()
    pair_rows = []
    all_species = sorted({species[p] for p in prediction})
    for i, a in enumerate(groups):
        for b in groups[i + 1 :]:
            null = strength[a] * strength[b] / (2 * total) if total else 0.0
            same = (
                sum(species_strength[a, s] * species_strength[b, s] for s in all_species)
                / (2 * total)
                if total
                else 0.0
            )
            expected["same_species"] += same
            expected["cross_species"] += null - same
            shared = {s for s in all_species if composition[a, s] and composition[b, s]}
            pair_rows.append(
                dict(
                    a=a,
                    b=b,
                    expected_same_weight=same,
                    expected_cross_weight=null - same,
                    observed_same_weight=pair_observed[a, b, "same_species"],
                    observed_cross_weight=pair_observed[a, b, "cross_species"],
                    observed_density=(
                        pair_observed[a, b, "same_species"] + pair_observed[a, b, "cross_species"]
                    )
                    / (genes[a] * genes[b]),
                    expected_density=null / (genes[a] * genes[b]),
                    shared_species=len(shared),
                )
            )
    terms = {}
    for k in ("same_species", "cross_species"):
        terms[k] = dict(
            total_graph_weight=within_species_weight[k],
            inside_cut_groups_weight=within_cut_weight[k],
            observed_between_groups_weight=observed[k],
            expected_between_groups_weight=expected[k],
            gain_contribution=(expected[k] - observed[k]) / total if total else 0.0,
        )
    gain = sum(t["gain_contribution"] for t in terms.values())
    split_score = (
        sum(within_cut_weight.values()) / total
        - sum((strength[g] / (2 * total)) ** 2 for g in groups)
        if total
        else 0.0
    )
    if not math.isclose(gain, split_score, rel_tol=1e-9, abs_tol=1e-9):
        raise ValueError("Gain decomposition mismatch")
    composition_rows = [
        dict(
            group=g,
            n_genes=genes[g],
            n_species=sum(composition[g, s] > 0 for s in all_species),
            species_counts={str(s): composition[g, s] for s in all_species if composition[g, s]},
            strength=strength[g],
        )
        for g in groups
    ]
    return dict(
        total_weight=total,
        groups=len(groups),
        species=len(all_species),
        gain=gain,
        species_terms=terms,
        group_composition=composition_rows,
        group_species_pair_overlap=pair_rows,
    )


def diagnose(run, boundary_dir, out):
    hashes = json.loads((boundary_dir / "input-hashes.json").read_text())

    def freeze(path):
        h = sha256_file(path)
        if str(path) in hashes and hashes[str(path)] != h:
            raise ValueError("Input identity changed")
        hashes[str(path)] = h

    for name in ("summary.json", "cases.json", "cuts-DIAGNOSTIC-ONLY.json", "input-hashes.json"):
        freeze(boundary_dir / name)
    summary = json.loads((boundary_dir / "summary.json").read_text())
    if not summary["immutable_inputs_verified"] or not summary["coverage_verified"]:
        raise ValueError("Unverified cut source")
    cases = json.loads((boundary_dir / "cases.json").read_text())
    cuts = json.loads((boundary_dir / "cuts-DIAGNOSTIC-ONLY.json").read_text())

    def read(path):
        freeze(path)
        return pq.ParquetFile(path).read().to_pylist()

    species = {r["protein_id"]: r["species_id"] for r in read(run / "input/proteins.parquet")}
    index = {r["protein_id"]: r["component_id"] for r in read(run / "components/index.parquet")}
    grouped = defaultdict(list)
    for case, cut in zip(cases, cuts, strict=True):
        if (case["component_id"], case["source_cluster_id"]) != (
            cut["component_id"],
            cut["source_cluster_id"],
        ):
            raise ValueError("Case cut identity mismatch")
        cut["prediction"] = {int(p): g for p, g in cut["prediction"].items()}
        if any(index[p] != case["component_id"] for p in cut["prediction"]):
            raise ValueError("Component mismatch")
        grouped[case["component_id"]].append((case, cut))
    rows = []
    for cid, items in sorted(grouped.items()):
        owner = {}
        edges = defaultdict(list)
        for i, (_case, cut) in enumerate(items):
            for p in cut["prediction"]:
                if p in owner:
                    raise ValueError("Overlapping source cuts")
                owner[p] = i
        folder = run / "components/edges" / f"component={cid:08d}"
        paths = sorted(folder.glob("*.parquet"))
        if not paths:
            raise ValueError("Missing component edges")
        for path in paths:
            for e in read(path):
                i, j = owner.get(e["u"]), owner.get(e["v"])
                if i is not None and i == j:
                    edges[i].append((e["u"], e["v"], e["weight"]))
        for i, (case, cut) in enumerate(items):
            result = decompose(cut["prediction"], species, edges[i])
            if not math.isclose(
                result["gain"], case["objective"]["cut_gain"], rel_tol=1e-8, abs_tol=1e-8
            ):
                raise ValueError("Frozen objective gain differs")
            rows.append(
                dict(
                    component_id=cid,
                    source_cluster_id=case["source_cluster_id"],
                    cohort=case["reference_cohort"],
                    event=case["original_event"],
                    removed_merge_tp=case["removed_merge_tp"],
                    removed_merge_fp=case["removed_merge_fp"],
                    **result,
                )
            )
    strata = []
    groups = defaultdict(list)
    for r in rows:
        if r["groups"] > 1:
            groups[r["cohort"], r["event"]].append(r)
    for (cohort, event), rs in sorted(groups.items()):
        contributions = {
            k: [r["species_terms"][k]["gain_contribution"] for r in rs]
            for k in ("same_species", "cross_species")
        }
        strata.append(
            dict(
                cohort=cohort,
                event=event,
                split_nodes=len(rs),
                positive_same_species_terms=sum(v > 1e-10 for v in contributions["same_species"]),
                positive_cross_species_terms=sum(v > 1e-10 for v in contributions["cross_species"]),
                larger_same_species_term=sum(
                    a > b
                    for a, b in zip(
                        contributions["same_species"], contributions["cross_species"], strict=True
                    )
                ),
                term_q10_q50_q90={
                    k: np.quantile(v, [0.1, 0.5, 0.9]).tolist() for k, v in contributions.items()
                },
            )
        )
    for path, digest in hashes.items():
        if sha256_file(Path(path)) != digest:
            raise ValueError("Frozen inputs changed")
    out.mkdir(parents=True, exist_ok=False)
    (out / "cases.json").write_text(json.dumps(rows, indent=2) + "\n")
    (out / "input-hashes.json").write_text(json.dumps(hashes, indent=2) + "\n")
    ids = {(0, 10040), (0, 21394), (0, 22807), (0, 14737)}
    (out / "summary.json").write_text(
        json.dumps(
            dict(
                diagnostic_only=True,
                production_changed=False,
                fixed_cuts_verified=True,
                immutable_inputs_verified=True,
                reference_labels_used_for_features=False,
                gain_identity=(
                    "sum (expected_cut_weight-observed_cut_weight)/W, "
                    "decomposed by same/cross species"
                ),
                target_nodes=len(rows),
                strata=strata,
                key_cases=[r for r in rows if (r["component_id"], r["source_cluster_id"]) in ids],
                scope=(
                    "Frozen source cuts; signed contributions describe null deficit, "
                    "not causal biological mechanisms"
                ),
            ),
            indent=2,
        )
        + "\n"
    )


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for key in ("run", "boundary-dir", "out"):
        p.add_argument("--" + key, type=Path, required=True)
    a = p.parse_args()
    diagnose(a.run, a.boundary_dir, a.out)


if __name__ == "__main__":
    main()
