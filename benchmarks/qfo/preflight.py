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


def sequence_map(path: Path):
    """Accession-keyed sequence identity, independent of FASTA wrapping/description/order."""
    result = {}
    for record in parse_fasta(path):
        parts = record.identifier.split("|")
        if len(parts) != 3 or parts[0] not in ("sp", "tr") or not parts[1]:
            raise ValueError(f"Unexpected UniProt ID: {record.identifier} in {path}")
        accession = parts[1]
        if accession in result:
            raise ValueError(f"Duplicate UniProt accession {accession} in {path}")
        result[accession] = (
            record.identifier,
            hashlib.sha256(record.sequence.encode()).hexdigest(),
        )
    return result


def compare_subset_file(path: Path, main: Path | None):
    before = sha256_file(path)
    if main is None:
        result = dict(
            file=path.name,
            classification="NOT_IN_ALL",
            sha256=before,
            main_sha256=None,
            usable_as_identical_subcollection=False,
        )
    elif before == sha256_file(main):
        result = dict(
            file=path.name,
            classification="BYTE_IDENTICAL",
            sha256=before,
            main_sha256=before,
            usable_as_identical_subcollection=True,
        )
    else:
        other, base = sequence_map(path), sequence_map(main)
        common = other.keys() & base.keys()
        added, missing = sorted(other.keys() - base.keys()), sorted(base.keys() - other.keys())
        changed = sorted(a for a in common if other[a][1] != base[a][1])
        id_changes = sorted(a for a in common if other[a][0] != base[a][0])
        sequence_identical = not (added or missing or changed)
        identical = sequence_identical and not id_changes
        result = dict(
            file=path.name,
            classification="FORMAT_OR_DESCRIPTION_ONLY"
            if identical
            else "ID_DIFFERENCE"
            if sequence_identical
            else "SEQUENCE_SET_DIFFERENCE",
            sha256=before,
            main_sha256=sha256_file(main),
            usable_as_identical_subcollection=identical,
            subset_proteins=len(other),
            main_proteins=len(base),
            added_accessions=len(added),
            missing_accessions=len(missing),
            changed_sequences=len(changed),
            changed_original_ids=len(id_changes),
            added_examples=added[:10],
            missing_examples=missing[:10],
            changed_sequence_examples=changed[:10],
            changed_id_examples=id_changes[:10],
        )
    if sha256_file(path) != before:
        raise ValueError(f"Subset input changed during audit: {path}")
    return result


def audit(source: Path, out: Path, *, collection: str = "all"):
    if collection not in {"all", "bacteria"}:
        raise ValueError("Unknown QFO collection")
    if out.exists():
        raise ValueError(f"Output already exists: {out}")
    out.mkdir(parents=True)
    frozen = out / "input" / collection
    frozen.mkdir(parents=True)
    primary = source / "all" if collection == "all" else source
    paths = sorted(primary.glob("*.fasta"))
    if not paths:
        raise ValueError(f"Missing FASTA in selected collection: {primary}")
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
        subsets, subset_diagnostics = {}, []
        for subset in ("bacteria", "eukaryota") if collection == "all" else ():
            subset_paths = sorted((source / subset).glob("*.fasta"))
            comparisons = []
            for p in subset_paths:
                main = frozen / p.name if p.name in manifests else None
                comparisons.append(compare_subset_file(p, main))
            subset_diagnostics.extend(dict(subset=subset, **r) for r in comparisons)
            subsets[subset] = dict(
                files=[p.name for p in subset_paths],
                n_species=len(subset_paths),
                n_proteins_in_all=sum(
                    manifests[p.name]["n_proteins"] for p in subset_paths if p.name in manifests
                ),
                missing_species_in_all=[p.name for p in subset_paths if p.name not in manifests],
                identical_subcollection=all(
                    r["usable_as_identical_subcollection"] for r in comparisons
                ),
                comparison_classifications=dict(Counter(r["classification"] for r in comparisons)),
                source_role="directory species list only; sequence input remains frozen all/",
            )
        (out / "subset-comparison.json").write_text(json.dumps(subset_diagnostics, indent=2) + "\n")
        if collection == "all" and (
            set(subsets["bacteria"]["files"]) & set(subsets["eukaryota"]["files"])
        ):
            raise ValueError("Overlapping supplied domain subsets")
        assigned = (
            set(subsets["bacteria"]["files"]) | set(subsets["eukaryota"]["files"])
            if collection == "all"
            else set(manifests)
        )
        digest = hashlib.sha256(
            json.dumps({k: v["sha256"] for k, v in manifests.items()}, sort_keys=True).encode()
        ).hexdigest()
        report = dict(
            input_ready=True,
            source_root=str(source),
            release="UNIDENTIFIED",
            release_scope="supplied local QFO FASTA collection; release not inferred",
            collection=collection,
            primary_input=str(primary),
            n_species=len(paths),
            n_proteins=len(owners),
            total_residues=sum(r["n_residues"] for r in rows),
            dataset_sha256=digest,
            original_ids_unique=True,
            accessions_unique=True,
            subset_audit_completed=collection == "all",
            subset_checks_passed=(
                all(r["identical_subcollection"] for r in subsets.values())
                if collection == "all"
                else None
            ),
            domain_stratification_scope=(
                "species lists from supplied directories, using all sequences only"
                if collection == "all"
                else "single bacterial collection; no cross-directory inputs"
            ),
            subsets=subsets,
            other_species=sorted(set(manifests) - assigned),
            frozen_inputs=str(frozen),
            official_qfo_scoring_ready=False,
            official_scoring_scope=(
                "requires reference/scoring package and validated orthology representation"
            ),
        )
        if collection == "all":
            report.update(all_species=len(paths), all_proteins=len(owners))
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
    p.add_argument("--collection", choices=("all", "bacteria"), default="all")
    a = p.parse_args()
    report = audit(a.source, a.out, collection=a.collection)
    print(
        json.dumps(
            {k: v for k, v in report.items() if k not in ("subsets", "other_species")}, indent=2
        )
    )


if __name__ == "__main__":
    main()
