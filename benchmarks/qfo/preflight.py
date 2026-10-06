"""Freeze and validate supplied QFO FASTA inputs; no orthology accuracy claims."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import shutil
from collections import Counter
from pathlib import Path

from ogprofiler.config import load_config
from ogprofiler.core.manifest import sha256_file
from ogprofiler.input.fasta import parse_fasta


def inspect_file(path: Path):
    n, residues = 0, 0
    accessions, header_types = set(), Counter()
    samples = []
    for record in parse_fasta(path):
        parts = record.identifier.split("|")
        if len(parts) != 3 or parts[0] not in ("sp", "tr") or not parts[1]:
            raise ValueError(f"Unexpected UniProt ID: {record.identifier} in {path}")
        accession = parts[1]
        if accession in accessions:
            raise ValueError(f"Duplicate UniProt accession {accession} in {path}")
        accessions.add(accession)
        header_types[parts[0]] += 1
        n += 1
        residues += len(record.sequence)
        if len(samples) < 20:
            samples.append(record)
    return dict(
        n_proteins=n,
        n_residues=residues,
        sp=header_types["sp"],
        tr=header_types["tr"],
        sha256=sha256_file(path),
        size_bytes=path.stat().st_size,
    ), samples


def audit(source: Path, out: Path):
    if out.exists():
        raise ValueError(f"Output already exists: {out}")
    out.mkdir(parents=True)
    frozen = out / "input/all"
    frozen.mkdir(parents=True)
    paths = sorted((source / "all").glob("*.fasta"))
    if not paths:
        raise ValueError("Missing all/*.fasta")
    manifests, rows, owners, accession_owners = {}, [], {}, {}
    smoke = out / "smoke-input"
    smoke.mkdir()
    map_path = out / "id-map.tsv"
    try:
        with map_path.open("w") as h:
            writer = csv.writer(h, delimiter="\t", lineterminator="\n")
            writer.writerow(["species", "taxid_from_filename", "original_id", "uniprot_accession"])
            for i, p in enumerate(paths):
                destination = frozen / p.name
                before = sha256_file(p)
                shutil.copy2(p, destination)
                metadata, samples = inspect_file(destination)
                if metadata["sha256"] != before or sha256_file(p) != before:
                    raise ValueError(f"Input changed during snapshot: {p}")
                manifests[p.name] = metadata
                rows.append(dict(file=p.name, species=p.stem, **metadata))
                for record in parse_fasta(destination):
                    accession = record.identifier.split("|")[1]
                    if record.identifier in owners or accession in accession_owners:
                        raise ValueError(
                            f"Cross-species duplicate original ID/accession: {record.identifier}"
                        )
                    owners[record.identifier] = p.stem
                    accession_owners[accession] = p.stem
                    writer.writerow(
                        [p.stem, p.stem.rsplit("_", 1)[-1], record.identifier, accession]
                    )
                if i < 3:
                    (smoke / p.name).write_text(
                        "".join(f">{r.identifier}\n{r.sequence}\n" for r in samples)
                    )
        subsets = {}
        for subset in ("bacteria", "eukaryota"):
            subset_paths = sorted((source / subset).glob("*.fasta"))
            mismatches = [
                p.name
                for p in subset_paths
                if p.name not in manifests or sha256_file(p) != manifests[p.name]["sha256"]
            ]
            if mismatches:
                raise ValueError(f"Subset differs from all: {subset}: {mismatches}")
            subsets[subset] = dict(
                files=[p.name for p in subset_paths],
                n_species=len(subset_paths),
                n_proteins=sum(manifests[p.name]["n_proteins"] for p in subset_paths),
            )
        if set(subsets["bacteria"]["files"]) & set(subsets["eukaryota"]["files"]):
            raise ValueError("Overlapping supplied domain subsets")
        assigned = set(subsets["bacteria"]["files"]) | set(subsets["eukaryota"]["files"])
        digest = hashlib.sha256(
            json.dumps({k: v["sha256"] for k, v in manifests.items()}, sort_keys=True).encode()
        ).hexdigest()
        report = dict(
            input_ready=True,
            source_root=str(source),
            release="UNIDENTIFIED",
            release_scope="supplied local QFO FASTA collection; release not inferred",
            all_species=len(paths),
            all_proteins=len(owners),
            total_residues=sum(r["n_residues"] for r in rows),
            dataset_sha256=digest,
            original_ids_unique=True,
            accessions_unique=True,
            subset_checks_passed=True,
            subsets=subsets,
            other_species=sorted(set(manifests) - assigned),
            frozen_inputs=str(frozen),
            official_qfo_scoring_ready=False,
            official_scoring_scope=(
                "requires reference/scoring package and validated orthology representation"
            ),
        )
        (out / "qfo-input-audit.json").write_text(json.dumps(report, indent=2) + "\n")
        (out / "input-sha256.json").write_text(
            json.dumps({k: v["sha256"] for k, v in manifests.items()}, indent=2) + "\n"
        )
        with (out / "input-manifest.tsv").open("w") as h:
            w = csv.DictWriter(h, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
            w.writeheader()
            w.writerows(rows)
        config = load_config()
        import yaml

        (out / "soft42-default.yaml").write_text(yaml.safe_dump(config, sort_keys=True))
        return report
    except Exception as error:
        (out / "preflight-failure.json").write_text(
            json.dumps(dict(input_ready=False, error=str(error)), indent=2) + "\n"
        )
        raise


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--source", type=Path, required=True)
    p.add_argument("--out", type=Path, required=True)
    a = p.parse_args()
    report = audit(a.source, a.out)
    print(
        json.dumps(
            {k: v for k, v in report.items() if k not in ("subsets", "other_species")}, indent=2
        )
    )


if __name__ == "__main__":
    main()
