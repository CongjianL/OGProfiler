#!/usr/bin/env python3
"""Compare current V2 and an installed OF3 source snapshot on identical hits.

Diagnostic only: never change production SSNs or run search/clustering.
Both paths use the V2 protein-ID ordering, isolating algorithm semantics from
OrthoFinder's original FASTA/file ordering. SciPy is a diagnostic dependency.
"""

from __future__ import annotations

import argparse
import ast
import csv
import hashlib
import json
import os
import shutil
import sys
import time
import warnings
from collections import Counter, defaultdict
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pyarrow
import pyarrow.parquet as pq
import scipy
from scipy import sparse
from scipy.optimize import curve_fit
from scipy.sparse.csgraph import connected_components

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from ogprofiler.similarity.engine import (
    EdgeBuildConfig,
    _deduplicate,
    best_hit_keys,
    build_retained_edges,
    lrb_cutoffs,
)
from ogprofiler.similarity.io import iter_directional_hits
from ogprofiler.similarity.models import NormalizedHit
from ogprofiler.similarity.normalization import (
    _fit_parameters,
    legacy_nbs,
    max_bitscore_hits,
    retain_top_data,
)


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def load_reference(root):
    """Execute unmodified selected definitions; exclude module orchestration."""
    paths = {
        "waterfall": root / "tools/waterfall.py",
        "matrices": root / "utils/matrices.py",
        "blast": root / "utils/blast_file_processor.py",
        "gathering": root / "orthogroups/gathering.py",
    }
    env = dict(np=np, numeric=np, sparse=sparse, curve_fit=curve_fit)
    tree = ast.parse(paths["matrices"].read_text())
    nodes = [n for n in tree.body if isinstance(n, ast.FunctionDef) and n.name == "sparse_max_row"]
    if len(nodes) != 1:
        raise ValueError("Reference sparse_max_row not found")
    exec(compile(ast.Module(body=nodes, type_ignores=[]), str(paths["matrices"]), "exec"), env)
    env["matrices"] = SimpleNamespace(sparse_max_row=env["sparse_max_row"])
    tree = ast.parse(paths["waterfall"].read_text())
    nodes = [
        n
        for n in tree.body
        if isinstance(n, ast.ClassDef) and n.name in {"scnorm", "WaterfallMethod"}
    ]
    if len(nodes) != 2:
        raise ValueError("Reference scnorm/WaterfallMethod not found")
    exec(compile(ast.Module(body=nodes, type_ignores=[]), str(paths["waterfall"]), "exec"), env)
    return env["scnorm"], env["WaterfallMethod"], paths


class Audit:
    def __init__(self, rows, out, atol=1e-8, rtol=1e-6, sample_limit=20):
        self.rows = {int(r["protein_id"]): r for r in rows}
        if not rows or len(self.rows) != len(rows):
            raise ValueError("Protein metadata must be nonempty with unique protein_id")
        self.species = sorted({int(r["species_id"]) for r in rows})
        self.species_index = {s: i for i, s in enumerate(self.species)}
        self.ids = [
            sorted(p for p, r in self.rows.items() if int(r["species_id"]) == s)
            for s in self.species
        ]
        self.local = {p: i for group in self.ids for i, p in enumerate(group)}
        self.global_ids = [p for group in self.ids for p in group]
        self.global_local = {p: i for i, p in enumerate(self.global_ids)}
        self.lengths = {p: int(r["length"]) for p, r in self.rows.items()}
        if any(length <= 0 for length in self.lengths.values()):
            raise ValueError("Protein lengths must be positive")
        self.length_arrays = [
            np.array([self.lengths[p] for p in group], dtype=float) for group in self.ids
        ]
        self.info = SimpleNamespace(
            nSpecies=len(self.species),
            speciesToUse=list(range(len(self.species))),
            nSeqsPerSpecies={s: len(ids) for s, ids in enumerate(self.ids)},
        )
        self.atol, self.rtol, self.limit = atol, rtol, sample_limit
        self.stats = []
        self.handle = (out / "difference_samples.tsv").open("w", newline="")
        self.writer = csv.DictWriter(
            self.handle,
            delimiter="\t",
            fieldnames=[
                "stage",
                "query_id",
                "query_species",
                "query_original_id",
                "target_id",
                "target_species",
                "target_original_id",
                "v2",
                "orthofinder",
                "abs_error",
            ],
        )
        self.writer.writeheader()

    def blank(self):
        return [[sparse.lil_matrix((len(a), len(b))) for b in self.ids] for a in self.ids]

    def matrices(self, keyed, get_score):
        groups = defaultdict(list)
        for (q, t), value in keyed.items():
            s = self.species_index[int(self.rows[q]["species_id"])]
            j = self.species_index[int(self.rows[t]["species_id"])]
            groups[s, j].append((self.local[q], self.local[t], get_score(value)))
        result = self.blank()
        for (s, j), values in groups.items():
            rr, cc, vv = zip(*values, strict=True)
            result[s][j] = sparse.csr_matrix((vv, (rr, cc)), shape=result[s][j].shape).tolil()
        return result

    def compare(self, stage, left, right, query_ids, target_ids):
        a, b = left.tocsr(copy=True), right.tocsr(copy=True)
        a.eliminate_zeros()
        b.eliminate_zeros()
        aa, bb = a.astype(bool), b.astype(bool)
        shared = aa.multiply(bb).nnz
        delta = (a - b).tocoo()

        def gather(matrix):
            if delta.nnz == 0:
                return np.empty(0, dtype=float)
            values = matrix[delta.row, delta.col]
            return np.asarray(values.toarray() if sparse.issparse(values) else values).ravel()

        av, bv = gather(a), gather(b)
        errors = np.abs(av - bv)
        scales = np.maximum(np.abs(av), np.abs(bv))
        mask = errors > self.atol + self.rtol * scales
        relative = errors / np.maximum(scales, 1e-300)
        stat = dict(
            stage=stage,
            v2_nonzero=a.nnz,
            of_nonzero=b.nnz,
            shared_support=shared,
            only_v2=a.nnz - shared,
            only_of=b.nnz - shared,
            exact_difference_count=delta.nnz,
            outside_tolerance_count=int(mask.sum()),
            max_abs_error=float(errors.max()) if len(errors) else 0.0,
            max_relative_error=float(relative.max()) if len(relative) else 0.0,
        )
        self.stats.append(stat)
        for k in np.flatnonzero(mask)[: self.limit]:
            q = query_ids[int(delta.row[k])]
            t = target_ids[int(delta.col[k])] if target_ids is not None else None
            qr = self.rows[q]
            tr = self.rows[t] if t is not None else {}
            self.writer.writerow(
                dict(
                    stage=stage,
                    query_id=q,
                    query_species=qr["species_id"],
                    query_original_id=qr["original_id"],
                    target_id=t,
                    target_species=tr.get("species_id", ""),
                    target_original_id=tr.get("original_id", ""),
                    v2=float(av[k]),
                    orthofinder=float(bv[k]),
                    abs_error=float(errors[k]),
                )
            )
        self.handle.flush()
        return stat

    def compare_blocks(self, name, a, b):
        for s in range(len(self.ids)):
            for j in range(len(self.ids)):
                self.compare(
                    f"{name}/{self.species[s]}->{self.species[j]}",
                    a[s][j],
                    b[s][j],
                    self.ids[s],
                    self.ids[j],
                )


def positive_support(matrix):
    return (matrix > 0).astype(np.int8)


def graph_summary(matrix):
    graph = positive_support(matrix)
    n, labels = connected_components(graph, directed=False)
    sizes = np.bincount(labels)
    return dict(
        n_vertices=graph.shape[0],
        directed_nonzero=graph.nnz,
        undirected_edges=sparse.triu(graph.maximum(graph.T), k=1).nnz,
        components=n,
        singletons=int((sizes == 1).sum()),
        largest_component=int(sizes.max()) if len(sizes) else 0,
        top_component_sizes=sorted(map(int, sizes), reverse=True)[:20],
    )


def run(hits, rows, reference_root, out, atol=1e-8, rtol=1e-6, sample_limit=20):
    out.mkdir(parents=True, exist_ok=False)
    snapshot = out / "reference" / "orthofinder"
    for relative in (
        "tools/waterfall.py",
        "utils/matrices.py",
        "utils/blast_file_processor.py",
        "orthogroups/gathering.py",
    ):
        path = reference_root / relative
        target = snapshot / relative
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(path, target)
    sc, wf, paths = load_reference(snapshot)
    audit = Audit(rows, out, atol, rtol, sample_limit)
    config = EdgeBuildConfig()
    started = time.monotonic()
    report = dict(
        algorithm="fixed-hits-layer-audit-v1",
        identity_order="V2 protein_id ascending within species",
        mode="default full-run OF scoring, not assign/v2_scores",
        v2_config=repr(config),
        dependencies=dict(
            python=sys.version,
            numpy=np.__version__,
            scipy=scipy.__version__,
            pyarrow=pyarrow.__version__,
        ),
        reference_sha256={name: sha256(p) for name, p in paths.items()},
        tolerances=dict(atol=atol, rtol=rtol),
        stages_complete=[],
    )

    def checkpoint(stage):
        report["stages_complete"].append(stage)
        report["elapsed_seconds"] = time.monotonic() - started
        (out / "progress.json").write_text(json.dumps(report, indent=2))
        print(f"[fixed-hits-audit] {stage}: {report['elapsed_seconds']:.1f}s", flush=True)

    usable = []
    for h in hits:
        for p, s in ((h.query_id, h.query_species), (h.target_id, h.target_species)):
            if p not in audit.rows or int(audit.rows[p]["species_id"]) != s:
                raise ValueError(f"Protein metadata mismatch: {p}/{s}")
        if not np.isfinite(h.bitscore) or not np.isfinite(h.evalue):
            raise ValueError("Non-finite raw hit")
        if h.query_id != h.target_id and h.bitscore > 0:
            usable.append(h)
    maximum = {}
    for h in usable:
        key = (h.query_id, h.target_id)
        maximum[key] = max(maximum.get(key, 0.0), h.bitscore)
    raw = audit.matrices(maximum, float)
    report["hits"] = dict(
        raw=len(hits),
        usable=len(usable),
        unique_usable=len(maximum),
        duplicate_usable=len(usable) - len(maximum),
    )
    # Apply the actual V2 max operator to bitscores as an input diagnostic.
    # Production V2 still deduplicates AFTER fitting in the path below.
    vraw_keys = _deduplicate(NormalizedHit(h, h.bitscore) for h in usable)
    vraw = audit.matrices(vraw_keys, lambda h: h.normalized_score)
    audit.compare_blocks("max_bitscore_input_diagnostic", vraw, raw)
    del maximum, vraw_keys, vraw
    checkpoint("max_bitscore")

    raw_counts = Counter((h.query_species, h.target_species) for h in usable)
    grouped = defaultdict(list)
    for h in max_bitscore_hits(hits):
        grouped[h.query_species, h.target_species].append(h)
    pair_stats = []
    official = audit.blank()
    for s, species in enumerate(audit.species):
        for j, target_species in enumerate(audit.species):
            group = grouped[species, target_species]
            lv2 = [float(audit.lengths[h.query_id] * audit.lengths[h.target_id]) for h in group]
            sv2 = [h.bitscore for h in group]
            tv2, bv2 = retain_top_data(lv2, sv2)
            pv2 = _fit_parameters(lv2, sv2, "v1_zero") if group else None
            li, lj, scores = sc.GetLengthArraysForMatrix(
                raw[s][j], audit.length_arrays[s], audit.length_arrays[j]
            )
            to, bo = sc.GetTopPercentileOfScores(li * lj, scores, 95)
            with warnings.catch_warnings(record=True) as caught:
                warnings.simplefilter("always")
                po = sc.CalculateFittingParameters(to, bo) if len(bo) > 1 else None
                official[s][j] = (
                    sc.NormaliseScoresByLogLengthProduct(
                        raw[s][j], audit.length_arrays[s], audit.length_arrays[j], po
                    )
                    if po is not None
                    else sparse.lil_matrix(raw[s][j].shape)
                )
            cv2, co = Counter(zip(tv2, bv2, strict=True)), Counter(zip(to, bo, strict=True))
            pair_stats.append(
                dict(
                    query_species=species,
                    target_species=target_species,
                    raw_usable=raw_counts[species, target_species],
                    max_dedup_hits=raw[s][j].nnz,
                    v2_fit_samples=len(tv2),
                    of_fit_samples=len(to),
                    fit_sample_only_v2=sum((cv2 - co).values()),
                    fit_sample_only_of=sum((co - cv2).values()),
                    fit_sample_examples_only_v2=list((cv2 - co).items())[:sample_limit],
                    fit_sample_examples_only_of=list((co - cv2).items())[:sample_limit],
                    v2_parameters=list(pv2) if pv2 else None,
                    of_parameters=list(map(float, po)) if po is not None else None,
                    of_warnings=[str(w.message) for w in caught],
                )
            )
            if not np.all(np.isfinite(official[s][j].tocsr().data)):
                raise ValueError(f"Non-finite OF B: {species}->{target_species}")
    report["normalization_pairs"] = pair_stats
    del grouped, raw, usable, lv2, sv2, group, li, lj, scores, tv2, bv2, to, bo, cv2, co
    checkpoint("NBS_samples_and_parameters")

    normalized = legacy_nbs(hits, audit.lengths, "v1_zero")
    keyed = _deduplicate(normalized)
    v2 = audit.matrices(keyed, lambda h: h.normalized_score)
    audit.compare_blocks("B", v2, official)
    checkpoint("B")

    best = best_hit_keys(keyed, config.best_hit_tolerance)
    reciprocal = {
        k
        for k in best
        if (k[1], k[0]) in best and keyed[k].hit.query_species != keyed[k].hit.target_species
    }
    vbh = audit.matrices(dict.fromkeys(best, 1.0), float)
    vrbh = audit.matrices(dict.fromkeys(reciprocal, 1.0), float)
    obh = [wf.GetBH_s(official[s], audit.info, s) for s in range(len(audit.ids))]
    orbh = [
        [
            obh[s][j].multiply(obh[j][s].T) if s != j else sparse.csr_matrix(obh[s][j].shape)
            for j in range(len(audit.ids))
        ]
        for s in range(len(audit.ids))
    ]
    audit.compare_blocks("BH", vbh, obh)
    audit.compare_blocks("RBH_cross_species", vrbh, orbh)
    # Control: hold B fixed to V2's complete matrix. Distinguish intrinsic
    # BH/cutoff differences from differences propagated by normalization.
    shared_bh = [wf.GetBH_s(v2[s], audit.info, s) for s in range(len(audit.ids))]
    audit.compare_blocks("BH_sameB_control", vbh, shared_bh)
    shared_rbhs = [
        [
            shared_bh[s][j].multiply(shared_bh[j][s].T)
            if s != j
            else sparse.csr_matrix(shared_bh[s][j].shape)
            for j in range(len(audit.ids))
        ]
        for s in range(len(audit.ids))
    ]
    audit.compare_blocks("RBH_sameB_control", vrbh, shared_rbhs)
    del vbh, vrbh, best
    checkpoint("BH_RBH")

    cutoffs = dict.fromkeys(audit.rows, 1e-6)
    cutoffs.update(lrb_cutoffs(keyed, config.best_hit_tolerance))
    connections = []
    shared_connections = []
    for s in range(len(audit.ids)):
        oc = wf.GetMostDistant_s(orbh[s], official[s], audit.info, s)
        vc = np.array([cutoffs[p] for p in audit.ids[s]])
        audit.compare(
            f"cutoff/{audit.species[s]}",
            sparse.csr_matrix(vc[:, None]),
            sparse.csr_matrix(oc[:, None]),
            audit.ids[s],
            None,
        )
        connections.append(wf.ConnectAllBetterThanCutoff_s(official[s], oc, audit.info, s))
        scut = wf.GetMostDistant_s(shared_rbhs[s], v2[s], audit.info, s)
        audit.compare(
            f"cutoff_sameB_control/{audit.species[s]}",
            sparse.csr_matrix(vc[:, None]),
            sparse.csr_matrix(scut[:, None]),
            audit.ids[s],
            None,
        )
        shared_connections.append(wf.ConnectAllBetterThanCutoff_s(v2[s], scut, audit.info, s))
    del reciprocal, obh, orbh, shared_bh, shared_rbhs
    checkpoint("cutoff")

    selected, edges = build_retained_edges(normalized, config)
    selected_keys = {(h.hit.query_id, h.hit.target_id) for h in selected}
    expected = {key for key, h in keyed.items() if h.normalized_score >= cutoffs[key[0]]}
    if selected_keys != expected:
        raise AssertionError("Instrumented V2 cutoff disagrees with production selector")
    vconnect = audit.matrices(dict.fromkeys(selected_keys, 1.0), float)
    audit.compare_blocks("connect", vconnect, connections)
    audit.compare_blocks("connect_sameB_control", vconnect, shared_connections)
    del shared_connections
    del hits, normalized, keyed, selected, selected_keys, expected, cutoffs
    checkpoint("connect")

    # Same directional assembly on each path's B/connect isolates upstream
    # differences. This is NOT the production V2 undirected graph.
    ow = []
    vw = []
    for s in range(len(audit.ids)):
        ow.append(
            [
                (connections[s][j] + connections[j][s].T).multiply(official[s][j]).tocsr()
                for j in range(len(audit.ids))
            ]
        )
        vw.append(
            [
                (vconnect[s][j] + vconnect[j][s].T).multiply(v2[s][j]).tocsr()
                for j in range(len(audit.ids))
            ]
        )
    audit.compare_blocks("directional_W_same_assembly_DIAGNOSTIC", vw, ow)
    W = sparse.bmat(ow, format="csr")
    del vconnect, connections, official, v2, vw, ow
    checkpoint("directional_W")

    rr = []
    cc = []
    vv = []
    for e in edges:
        u, v = audit.global_local[e.u], audit.global_local[e.v]
        rr.extend((u, v))
        cc.extend((v, u))
        vv.extend((e.weight, e.weight))
    U = sparse.csr_matrix((vv, (rr, cc)), shape=W.shape)
    del edges, rr, cc, vv
    # U is an explicitly labelled symmetric lift of the production V2 SSN.
    audit.compare("OF_directional_W_vs_V2_symmetric_lift", U, W, audit.global_ids, audit.global_ids)
    rounded = W.copy()
    rounded.data = np.round(rounded.data, 3)
    rounded.eliminate_zeros()
    projections = {"mean": (W + W.T) * 0.5, "max": W.maximum(W.T)}
    for name, P in projections.items():
        audit.compare(
            f"OF_{name}_projection_vs_V2_undirected",
            sparse.triu(U, k=1),
            sparse.triu(P, k=1),
            audit.global_ids,
            audit.global_ids,
        )
    report["graph_objects"] = dict(
        V2_undirected=graph_summary(U),
        OF_directional_full_precision=graph_summary(W),
        OF_directional_written_3dp=graph_summary(rounded),
    )
    report["projection_contract"] = dict(
        production_v2=(
            "one canonical undirected edge; mean of complete OF-assembled directional W; "
            "igraph directed=False"
        ),
        orthofinder=(
            "directional W=(C+C.T) elementwise B; writer formats %.3f; MCL consumes the matrix"
        ),
        diagnostic_directional_v2=(
            "same OF assembly applied to V2 B/connect; not V2 production output"
        ),
        diagnostic_symmetric_lift="V2 Uuv placed in both directions solely for comparisons",
        candidate_projections={"mean": "(W+W.T)/2", "max": "max(W,W.T)"},
        selected_projection=config.symmetrization,
        note=(
            "Candidate projections are diagnostics, "
            "not exact replacements for OF's directed MCL input."
        ),
    )
    report["stage_comparisons"] = audit.stats
    totals = {}
    for row in audit.stats:
        name = row["stage"].split("/")[0]
        total = totals.setdefault(name, {"stage": name})
        for key, value in row.items():
            if key == "stage":
                continue
            total[key] = (
                max(total.get(key, 0), value)
                if key.startswith("max_")
                else total.get(key, 0) + value
            )
    report["stage_totals"] = list(totals.values())
    checkpoint("graph_objects_and_projections")
    audit.handle.close()
    with (out / "stage_summary.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(audit.stats[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(audit.stats)
    with (out / "stage_totals.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(audit.stats[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(totals.values())
    with (out / "normalization_pairs.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(pair_stats[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(pair_stats)
    report["completed"] = True
    (out / "report.json").write_text(json.dumps(report, indent=2))
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run", required=True, type=Path)
    parser.add_argument("--of-source", required=True, type=Path)
    parser.add_argument("--out", required=True, type=Path)
    parser.add_argument("--sample-limit", type=int, default=20)
    args = parser.parse_args()
    if args.sample_limit < 0:
        parser.error("sample-limit must be nonnegative")
    hits_path = args.run / "search/hits.parquet"
    proteins_path = args.run / "input/proteins.parquet"
    for path in (hits_path, proteins_path):
        if not path.is_file():
            parser.error(f"Missing fixed input: {path}")
    if args.out.exists():
        parser.error("Output must be a fresh directory")
    checksums = {str(p): sha256(p) for p in (hits_path, proteins_path)}
    rows = pq.read_table(
        proteins_path, columns=["protein_id", "species_id", "original_id", "length"]
    ).to_pylist()
    hits = list(iter_directional_hits(hits_path))
    report = run(hits, rows, args.of_source, args.out, sample_limit=args.sample_limit)
    report["input_checksums"] = checksums
    report["input_run"] = str(args.run)
    report["scheduler_provenance"] = {
        k: os.environ.get(k)
        for k in (
            "DEV_RUN_ID",
            "DEV_GIT_COMMIT",
            "DEV_GIT_DIRTY",
            "DEV_SOURCE_HASH",
            "SLURM_JOB_ID",
        )
    }
    report["completed"] = False
    report["input_integrity_verified"] = False
    (args.out / "report.json").write_text(json.dumps(report, indent=2))
    for path in (hits_path, proteins_path):
        if sha256(path) != checksums[str(path)]:
            raise ValueError(f"Input changed during audit: {path}")
    report["input_integrity_verified"] = True
    report["completed"] = True
    (args.out / "report.json").write_text(json.dumps(report, indent=2))
    print(json.dumps(report["graph_objects"], indent=2), flush=True)


if __name__ == "__main__":
    main()
