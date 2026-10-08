#!/usr/bin/env python3
"""Add verified default V2 and soft42 to all Fig2 panels, preserving prior rows."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
import statistics
from collections import Counter, defaultdict
from pathlib import Path

import yaml

TABLES = (
    "official_precision_recall.tsv",
    "refog_F1_long.tsv",
    "split_contamination.tsv",
    "fp_concentration_long.tsv",
)
CONDITIONS = {
    "OGProfiler2_default": ("1410777", "kway_v1", 20),
    "OGProfiler2_soft42": ("1411546", "soft_binary_24_v2", 42),
}


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read(path):
    with path.open() as h:
        return list(csv.DictReader(h, delimiter="\t"))


def write(path, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as h:
        w = csv.DictWriter(h, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)


def predictions(path):
    groups = {}
    all_ids = set()
    for line in path.read_text().splitlines():
        gid, members = line.split(":", 1)
        members = members.split()
        if gid in groups or not members or len(members) != len(set(members)):
            raise ValueError("Invalid prediction group")
        if all_ids.intersection(members):
            raise ValueError("Duplicate protein assignment")
        all_ids.update(members)
        groups[gid] = set(members)
    return groups


def diagnostics(groups, truth, low):
    """Raw best-group metrics match existing Fig2; FP matches official exclusion."""
    owner = {p: gid for gid, genes in groups.items() for p in genes}
    refs, burden = [], defaultdict(float)
    tp_total, fp_total, fn_total = 0.0, 0.0, 0.0
    for name, raw in sorted(truth.items()):
        overlap = Counter(owner[p] for p in raw if p in owner)
        choices = [(2 * n / (len(raw) + len(groups[g])), g) for g, n in overlap.items()]
        f1, best = max(choices) if choices else (0.0, "")
        genes = groups.get(best, set())
        n = len(raw & genes)
        refs.append(
            dict(
                refog=name,
                n_true=len(raw),
                best_predicted_group=best,
                intersection=n,
                precision=n / len(genes) if genes else 0,
                recall=n / len(raw),
                F1=f1,
                exact=str(genes == raw).lower(),
                split_count=max(0, len(overlap) - 1),
                contamination=len(genes - raw) / len(genes) if genes else 0,
                missing_fraction=len(raw - genes) / len(raw),
            )
        )
        uncertain = low.get(name, set())
        ref = raw - uncertain
        denom = len(ref) - 1
        if denom <= 0:
            raise ValueError("Unsupported singleton confident RefOG")
        confident = Counter(owner[p] for p in ref if p in owner)
        for gid, n in confident.items():
            pred_size = len(groups[gid] - uncertain)
            tp_total += n * (n - 1) / 2 / denom
            fp = n * (pred_size - n) / denom
            fp_total += fp
            burden[gid] += fp
            fn_total += n * (len(ref) - n) / 2 / denom
        if ref - owner.keys():
            raise ValueError("Missing reference genes in full-coverage figure inputs")
    ranked = sorted(burden.items(), key=lambda x: (-x[1], x[0]))
    concentration = {
        n: sum(v for _, v in ranked[:n]) / fp_total if fp_total else 0 for n in (1, 5, 10)
    }
    p = tp_total / (tp_total + fp_total) if tp_total + fp_total else 0
    r = tp_total / (tp_total + fn_total) if tp_total + fn_total else 0
    return (
        refs,
        concentration,
        {"precision": p, "recall": r, "F_score": 2 * p * r / (p + r) if p + r else 0},
    )


def build(base, runs, out, metrics_out):
    tables = {name: read(base / name) for name in TABLES}
    for rows in tables.values():
        if any(r["method"] in CONDITIONS for r in rows):
            raise ValueError("Use the preserved pre-update base, not an already updated table")
    raw_root = runs["OGProfiler2_default"] / "benchmark/RefOGs"
    truth = {p.stem: set(p.read_text().splitlines()) for p in sorted(raw_root.glob("RefOG*.txt"))}
    low = {
        p.stem: set(p.read_text().splitlines())
        for p in (raw_root / "low_certainty_assignments").glob("RefOG*.txt")
    }
    assert len(truth) == 70
    provenance = {
        "base_input_sha256": {n: sha(base / n) for n in TABLES},
        "reference_sha256": {
            str(p.relative_to(raw_root)): sha(p) for p in raw_root.rglob("*") if p.is_file()
        },
        "conditions": {},
        "existing_rows_preserved": True,
        "refog_diagnostics": "raw RefOG best-group metrics, unchanged Fig2 convention",
        "fp_concentration": (
            "official-like normalized per-RefOG contributions with low-certainty exclusions"
        ),
        "runtime_comparison": "not added: fixed-SSN runs are not end-to-end runtime controls",
    }
    metrics_out.mkdir(parents=True, exist_ok=True)
    for method, run in runs.items():
        job, policy, depth = CONDITIONS[method]
        report = json.loads((run / "h5-metrics/report.json").read_text())
        for k in (
            "evaluation_completed",
            "strategy_parity_passed",
            "fixed_ssn_unchanged",
            "benchmark_unchanged",
        ):
            assert report[k]
        config = yaml.safe_load((run / "new-hierarchy/run.yaml").read_text())
        assert config["hierarchy"].get("topology_policy", "kway_v1") == policy
        assert (
            config["hierarchy"]["max_depth"] == depth
            and config["hierarchy"]["leiden_iterations"] == 10
        )
        assert config["orthogroups"]["strategy"] == "v1_compatible"
        folder = run / "h5-metrics/new_bounded_hierarchy_v1_compatible"
        metrics = json.loads((folder / "metrics.json").read_text())
        prediction = folder / "prediction.txt"
        assert sha(prediction) == metrics["prediction_sha256"]
        assert metrics["mapping_validated"] and metrics["official_reader_exact"]
        assert metrics["coverage"] == metrics["refog_raw_coverage"] == 1
        groups = predictions(prediction)
        assert sum(map(len, groups.values())) == metrics["input_proteins"] == 251378
        refs, fp, pairwise = diagnostics(groups, truth, low)
        for key in ("precision", "recall"):
            assert abs(pairwise[key] - metrics["official_" + key]) < 1e-10
        assert abs(pairwise["F_score"] - metrics["official_f1"]) < 1e-10
        exact = re.search(
            r"(\d+)\s+orthogroups? exactly correct", (folder / "official.stdout").read_text(), re.I
        )
        if not exact:
            raise ValueError("Official exact count missing")
        version = f"job{job}; {policy}; depth{depth}; fixed mean SSN"
        tables[TABLES[0]].append(
            dict(
                method=method,
                version=version,
                precision=metrics["official_precision"] * 100,
                recall=metrics["official_recall"] * 100,
                F_score=metrics["official_f1"] * 100,
                official_exact=int(exact[1]),
                status="PASS",
            )
        )
        tables[TABLES[1]].extend(dict(RefOG=r["refog"], method=method, F1=r["F1"]) for r in refs)
        tables[TABLES[2]].append(
            dict(
                method=method,
                median_split=statistics.median(r["split_count"] for r in refs),
                mean_split=statistics.mean(r["split_count"] for r in refs),
                median_contamination=statistics.median(r["contamination"] for r in refs),
                mean_contamination=statistics.mean(r["contamination"] for r in refs),
            )
        )
        tables[TABLES[3]].extend(
            dict(method=method, top_n=n, cumulative_FP_fraction=fp[n]) for n in (1, 5, 10)
        )
        write(metrics_out / method / "refog_metrics.tsv", refs)
        write(
            metrics_out / method / "fp_concentration.tsv",
            [dict(top_n=n, cumulative_FP_fraction=fp[n]) for n in (1, 5, 10)],
        )
        provenance["conditions"][method] = dict(
            job_id=job,
            run_id={
                "OGProfiler2_default": "20261002T060705Z_cb0c0ba98049_1cb4fcd7_14702",
                "OGProfiler2_soft42": "20261005T232747Z_8bcccc9d3302_512d821e_28231",
            }[method],
            prediction_sha256=sha(prediction),
            run_config_sha256=sha(run / "new-hierarchy/run.yaml"),
            official_metrics=metrics,
            calculated_pairwise=pairwise,
        )
    for name, rows in tables.items():
        write(out / name, rows)
    provenance["output_sha256"] = {n: sha(out / n) for n in TABLES}
    (out / "V2_CONFIG_UPDATE_PROVENANCE.json").write_text(json.dumps(provenance, indent=2) + "\n")
    write(metrics_out / "official_summary.tsv", tables[TABLES[0]])
    write(metrics_out / "split_contamination.tsv", tables[TABLES[2]])


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for name in ("base", "default-run", "soft42-run", "out", "metrics-out"):
        p.add_argument("--" + name, type=Path, required=True)
    a = p.parse_args()
    build(
        a.base,
        {"OGProfiler2_default": a.default_run, "OGProfiler2_soft42": a.soft42_run},
        a.out,
        a.metrics_out,
    )


if __name__ == "__main__":
    main()
