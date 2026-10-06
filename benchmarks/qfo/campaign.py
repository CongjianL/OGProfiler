"""Stage bacteria-only method inputs with a lossless historical V1 ID adapter."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import shutil
from itertools import islice
from pathlib import Path

from ogprofiler.core.manifest import sha256_file
from ogprofiler.input.fasta import parse_fasta


def stage(origin: Path, out: Path, expected_digest: str, *, mode: str):
    if mode not in {"smoke", "full"}:
        raise ValueError("Unknown campaign mode")
    report = json.loads((origin / "preflight/qfo-input-audit.json").read_text())
    completion = json.loads((origin / "qfo-preflight-completion.json").read_text())
    if not completion["preflight_completed"] or not report["input_ready"]:
        raise ValueError("Origin preflight incomplete")
    if report["collection"] != "bacteria" or report["dataset_sha256"] != expected_digest:
        raise ValueError("Unexpected collection/dataset digest")
    frozen = origin / "preflight/input/bacteria"
    full_files = sorted(frozen.glob("*.fasta"))
    full_hashes = {p.name: sha256_file(p) for p in full_files}
    full_digest = hashlib.sha256(json.dumps(full_hashes, sort_keys=True).encode()).hexdigest()
    if full_digest != expected_digest:
        raise ValueError("Frozen input changed")
    source = frozen if mode == "full" else origin / "preflight/smoke-input"
    files = sorted(source.glob("*.fasta"))
    if not files:
        raise ValueError("No method inputs")
    if mode == "smoke":
        if [p.name for p in files] != [p.name for p in full_files[:3]]:
            raise ValueError("Smoke species selection changed")
        for p in files:
            expected = list(islice(parse_fasta(frozen / p.name), 20))
            observed = list(parse_fasta(p))
            if [(r.identifier, r.sequence) for r in expected] != [
                (r.identifier, r.sequence) for r in observed
            ]:
                raise ValueError("Smoke records changed")
    hashes = {p.name: sha256_file(p) for p in files}
    digest = hashlib.sha256(json.dumps(hashes, sort_keys=True).encode()).hexdigest()
    if mode == "full" and digest != expected_digest:
        raise ValueError("Frozen input changed")
    out.mkdir(parents=True, exist_ok=False)
    common, v1 = out / "input", out / "v1-input"
    common.mkdir()
    v1.mkdir()
    seen = set()
    with (out / "v1-id-map.tsv").open("w") as h:
        w = csv.writer(h, delimiter="\t", lineterminator="\n")
        w.writerow(["adapted_id", "original_id"])
        for species, p in enumerate(files, 1):
            destination = common / p.name
            shutil.copy2(p, destination)
            if sha256_file(destination) != hashes[p.name] or sha256_file(p) != hashes[p.name]:
                raise ValueError("Input changed during staging")
            with (v1 / p.name).open("w") as f:
                for index, record in enumerate(parse_fasta(destination), 1):
                    if record.identifier in seen:
                        raise ValueError("Duplicate original ID")
                    seen.add(record.identifier)
                    adapted = f"OGPV1_{species:02d}|QFO_{species:02d}_{index:08d}"
                    f.write(f">{adapted}\n{record.sequence}\n")
                    w.writerow([adapted, record.identifier])
    if mode == "full" and (len(files) != report["n_species"] or len(seen) != report["n_proteins"]):
        raise ValueError("Frozen input counts changed")
    shutil.copy2(origin / "preflight/soft42-default.yaml", out / "soft42-default.yaml")
    summary = dict(
        origin=str(origin),
        mode=mode,
        collection="bacteria",
        n_species=len(files),
        n_proteins=len(seen),
        dataset_sha256=digest,
        full_dataset_sha256=expected_digest,
        official_qfo_scores_produced=False,
        v1_adapter="bijective species|record ID; sequences and species unchanged",
        v1_seed="unavailable in historical CLI",
        input_sha256=hashes,
    )
    (out / "campaign-input.json").write_text(json.dumps(summary, indent=2) + "\n")
    return summary


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--origin", type=Path, required=True)
    p.add_argument("--out", type=Path, required=True)
    p.add_argument("--expected-digest", required=True)
    p.add_argument("--mode", choices=("smoke", "full"), required=True)
    a = p.parse_args()
    print(json.dumps(stage(a.origin, a.out, a.expected_digest, mode=a.mode), indent=2))


if __name__ == "__main__":
    main()
