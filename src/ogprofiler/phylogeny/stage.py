"""Selected-family alignment, gene-tree, rooting, and reconciliation stage."""

from __future__ import annotations

import csv
import json
import os
import shutil
import uuid
from pathlib import Path
from typing import Any, cast

import pyarrow.parquet as pq

from ogprofiler.core.manifest import sha256_file, sha256_json, write_json
from ogprofiler.exceptions import PhylogenyError
from ogprofiler.input.fasta import parse_fasta
from ogprofiler.phylogeny.backends import AlignmentBackend, TreeBackend
from ogprofiler.phylogeny.newick import (
    leaf_names,
    midpoint_root,
    outgroup_root,
    parse_newick,
    prune_tree,
    to_newick,
)
from ogprofiler.phylogeny.reconciliation import (
    ReconciliationBackend,
    species_tree_aware_root,
)
from ogprofiler.phylogeny.selection import (
    RefinementFamily,
    select_refinement_families,
    terminal_result_paths,
)

PHYLOGENY_ALGORITHM_VERSION = "selected-terminal-family-phylogeny-v2"


def _parquet_rows(path: Path) -> list[dict[str, Any]]:
    try:
        return cast(list[dict[str, Any]], pq.read_table(path).to_pylist())
    except Exception as error:
        raise PhylogenyError(f"Failed to read phylogeny input {path}: {error}") from error


def _tsv(path: Path) -> list[dict[str, str]]:
    if not path.is_file():
        raise PhylogenyError(f"Missing phylogeny input: {path}")
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def _write_tsv(path: Path, fields: list[str], rows: list[dict[str, Any]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def _prepared_sequences(path: Path) -> dict[int, str]:
    sequences: dict[int, str] = {}
    for record in parse_fasta(path, "error"):
        if not record.identifier.startswith("OGP2P"):
            raise PhylogenyError(f"Unexpected prepared sequence ID: {record.identifier}")
        sequences[int(record.identifier[5:])] = record.sequence
    return sequences


def _family_members(run_root: Path) -> dict[str, list[int]]:
    result: dict[str, list[int]] = {}
    for row in _tsv(terminal_result_paths(run_root)[1]):
        result.setdefault(row["family_id"], []).append(int(row["protein_id"]))
    for values in result.values():
        values.sort()
    return result


def classify_event_conflict(network_event: str, phylo_event: str) -> str:
    network_class = {
        "SPECIATION_LIKE": "SPECIATION",
        "POLYTOMY": "SPECIATION",
        "DUPLICATION_LIKE": "DUPLICATION",
    }.get(network_event)
    if network_class is None or phylo_event == "UNRESOLVED":
        return "UNRESOLVED"
    return "CONCORDANT" if network_class == phylo_event else "CONFLICT"


def _verified(directory: Path, identity: str) -> bool:
    manifest_path = directory / "phylogeny-manifest.json"
    if not manifest_path.is_file():
        return False
    try:
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        if manifest["identity_sha256"] != identity:
            return False
        return all(
            (directory / name).is_file() and sha256_file(directory / name) == digest
            for name, digest in manifest["output_checksums"].items()
        )
    except (OSError, KeyError, TypeError, ValueError, json.JSONDecodeError):
        return False


def _write_family_fasta(
    path: Path,
    protein_ids: list[int],
    sequences: dict[int, str],
) -> None:
    with path.open("w", encoding="utf-8", newline="\n") as handle:
        for protein_id in protein_ids:
            try:
                sequence = sequences[protein_id]
            except KeyError as error:
                raise PhylogenyError(
                    f"Missing prepared sequence for protein {protein_id}"
                ) from error
            handle.write(f">P{protein_id:012d}\n")
            for start in range(0, len(sequence), 80):
                handle.write(sequence[start : start + 80] + "\n")


def _refine_family(
    *,
    run_root: Path,
    family: RefinementFamily,
    members: list[int],
    sequences: dict[int, str],
    species_name_by_id: dict[int, str],
    species_by_protein: dict[int, int],
    species_tree_text: str | None,
    rooting: str,
    outgroup: str | None,
    threads: int,
    alignment_backend: AlignmentBackend,
    tree_backend: TreeBackend,
    reconciliation_backend: ReconciliationBackend,
    command: list[str],
) -> tuple[dict[str, Any], bool]:
    output = run_root / "evolution" / "phylogenetic" / f"family={family.family_id}"
    parameters = {
        "algorithm_version": PHYLOGENY_ALGORITHM_VERSION,
        "grouping_source": "terminal_family",
        "component_id": family.component_id,
        "cluster_id": family.cluster_id,
        "network_event": family.network_event,
        "family_id": family.family_id,
        "members": members,
        "rooting": rooting,
        "outgroup": outgroup,
        "threads": threads,
        "alignment_backend": alignment_backend.name,
        "tree_backend": tree_backend.name,
        "reconciliation_backend": reconciliation_backend.name,
        "species_tree_sha256": sha256_json(species_tree_text),
        "member_sequence_sha256": sha256_json(
            [[protein_id, sequences[protein_id]] for protein_id in members]
        ),
        "leaf_species": [
            [protein_id, species_name_by_id[species_by_protein[protein_id]]]
            for protein_id in members
        ],
        "alignment_backend_version": alignment_backend.version(),
        "tree_backend_version": tree_backend.version(),
    }
    identity = sha256_json(parameters)
    if _verified(output, identity):
        manifest = json.loads((output / "phylogeny-manifest.json").read_text(encoding="utf-8"))
        return dict(manifest["summary"]), True

    if len(members) < 2:
        summary = {
            "family_id": family.family_id,
            "component_id": family.component_id,
            "cluster_id": family.cluster_id,
            "network_event": family.network_event,
            "phylo_event": "UNRESOLVED",
            "event_confidence": 0.0,
            "supporting_node": "",
            "conflict_status": "UNRESOLVED",
            "selection_reasons": ",".join(family.selection_reasons),
        }
        if output.exists():
            shutil.rmtree(output)
        output.mkdir(parents=True)
        write_json(
            output / "phylogeny-manifest.json",
            {
                "algorithm_version": PHYLOGENY_ALGORITHM_VERSION,
                "identity_sha256": identity,
                "parameters": parameters,
                "command": command,
                "output_checksums": {},
                "summary": summary,
                "status": "TOO_FEW_SEQUENCES",
            },
        )
        return summary, False

    build = output.with_name(f".{output.name}.{uuid.uuid4().hex}.building")
    if build.exists():
        shutil.rmtree(build)
    build.mkdir(parents=True)
    try:
        input_fasta = build / "input.faa"
        alignment = build / "alignment.faa"
        unrooted_path = build / "gene_tree.unrooted.nwk"
        rooted_path = build / "gene_tree.rooted.nwk"
        _write_family_fasta(input_fasta, members, sequences)
        alignment_run = alignment_backend.align(input_fasta, alignment, threads)
        tree_run = tree_backend.infer(alignment, unrooted_path)
        gene_tree = parse_newick(unrooted_path.read_text(encoding="utf-8"))
        expected_leaves = {f"P{protein_id:012d}" for protein_id in members}
        if set(leaf_names(gene_tree)) != expected_leaves:
            raise PhylogenyError(f"Gene-tree leaves differ from family {family.family_id}")
        leaf_species = {
            f"P{protein_id:012d}": species_name_by_id[species_by_protein[protein_id]]
            for protein_id in members
        }
        species_tree = None
        if species_tree_text is not None:
            parsed_species_tree = parse_newick(species_tree_text)
            missing = sorted(set(leaf_species.values()) - set(leaf_names(parsed_species_tree)))
            if missing:
                raise PhylogenyError(f"Species tree is missing selected species: {missing[0]}")
            species_tree = prune_tree(parsed_species_tree, set(leaf_species.values()))
            (build / "species_tree.pruned.nwk").write_text(
                to_newick(species_tree) + "\n", encoding="utf-8"
            )
        if rooting == "midpoint":
            rooted = midpoint_root(gene_tree)
        elif rooting == "outgroup":
            if outgroup is None:
                raise PhylogenyError("Outgroup rooting requires --outgroup")
            rooted = outgroup_root(gene_tree, outgroup)
        else:
            if species_tree is None:
                raise PhylogenyError("Species-tree-aware rooting requires --species-tree")
            rooted = species_tree_aware_root(
                gene_tree, species_tree, leaf_species, reconciliation_backend
            )
        rooted_path.write_text(to_newick(rooted) + "\n", encoding="utf-8")
        reconciliation = reconciliation_backend.annotate(rooted, species_tree, leaf_species)
        reconciliation_rows = [
            {
                "supporting_node": event.supporting_node,
                "phylo_event": event.phylo_event,
                "confidence": event.confidence,
                "species_node": event.species_node,
            }
            for event in reconciliation.events
        ]
        _write_tsv(
            build / "reconciliation.tsv",
            ["supporting_node", "phylo_event", "confidence", "species_node"],
            reconciliation_rows,
        )
        supporting_node = reconciliation.events[-1].supporting_node if reconciliation.events else ""
        summary = {
            "family_id": family.family_id,
            "component_id": family.component_id,
            "cluster_id": family.cluster_id,
            "network_event": family.network_event,
            "phylo_event": reconciliation.root_event,
            "event_confidence": reconciliation.confidence,
            "supporting_node": supporting_node,
            "conflict_status": classify_event_conflict(
                family.network_event, reconciliation.root_event
            ),
            "selection_reasons": ",".join(family.selection_reasons),
        }
        outputs = sorted(path for path in build.iterdir() if path.is_file())
        write_json(
            build / "phylogeny-manifest.json",
            {
                "algorithm_version": PHYLOGENY_ALGORITHM_VERSION,
                "identity_sha256": identity,
                "parameters": parameters,
                "command": command,
                "backend_versions": {
                    "alignment": alignment_run.version,
                    "tree": tree_run.version,
                    "reconciliation": reconciliation_backend.name,
                },
                "backend_commands": {
                    "alignment": list(alignment_run.command),
                    "tree": list(tree_run.command),
                },
                "output_checksums": {path.name: sha256_file(path) for path in outputs},
                "summary": summary,
                "status": "DONE",
            },
        )
        if output.exists():
            shutil.rmtree(output)
        os.replace(build, output)
        return summary, False
    except Exception:
        shutil.rmtree(build, ignore_errors=True)
        raise


def run_phylogenetic_refinement_stage(
    run_root: Path,
    command: list[str],
    *,
    explicit_family_ids: tuple[str, ...],
    selection_events: set[str],
    large_family_size: int,
    max_families: int,
    rooting: str,
    outgroup: str | None,
    species_tree_path: Path | None,
    threads: int,
    alignment_backend: AlignmentBackend,
    tree_backend: TreeBackend,
    reconciliation_backend: ReconciliationBackend,
) -> tuple[Path, int, int]:
    selected = select_refinement_families(
        run_root,
        explicit_family_ids=explicit_family_ids,
        selection_events=selection_events,
        large_family_size=large_family_size,
        max_families=max_families,
    )
    members_by_family = _family_members(run_root)
    sequences = _prepared_sequences(run_root / "input" / "proteins.faa")
    proteins = _parquet_rows(run_root / "input" / "proteins.parquet")
    species_rows = _parquet_rows(run_root / "input" / "species.parquet")
    species_by_protein = {int(row["protein_id"]): int(row["species_id"]) for row in proteins}
    species_name_by_id = {int(row["species_id"]): str(row["species_name"]) for row in species_rows}
    try:
        species_tree_text = (
            species_tree_path.read_text(encoding="utf-8") if species_tree_path is not None else None
        )
    except OSError as error:
        raise PhylogenyError(f"Failed to read species tree {species_tree_path}: {error}") from error
    summaries: list[dict[str, Any]] = []
    reused = 0
    for family in selected:
        summary, was_reused = _refine_family(
            run_root=run_root,
            family=family,
            members=members_by_family[family.family_id],
            sequences=sequences,
            species_name_by_id=species_name_by_id,
            species_by_protein=species_by_protein,
            species_tree_text=species_tree_text,
            rooting=rooting,
            outgroup=outgroup,
            threads=threads,
            alignment_backend=alignment_backend,
            tree_backend=tree_backend,
            reconciliation_backend=reconciliation_backend,
            command=command,
        )
        summaries.append(summary)
        reused += int(was_reused)
    output_root = run_root / "evolution" / "phylogenetic"
    output_root.mkdir(parents=True, exist_ok=True)
    fields = [
        "family_id",
        "component_id",
        "cluster_id",
        "network_event",
        "phylo_event",
        "event_confidence",
        "supporting_node",
        "conflict_status",
        "selection_reasons",
    ]
    events_path = output_root / "phylogenetic-events.tsv"
    temporary = events_path.with_name(f".{events_path.name}.{uuid.uuid4().hex}.tmp")
    _write_tsv(temporary, fields, summaries)
    os.replace(temporary, events_path)
    manifest_path = output_root / "phylogenetic-manifest.json"
    write_json(
        manifest_path,
        {
            "algorithm_version": PHYLOGENY_ALGORITHM_VERSION,
            "command": command,
            "selection": {
                "explicit_family_ids": list(explicit_family_ids),
                "selection_events": sorted(selection_events),
                "large_family_size": large_family_size,
                "max_families": max_families,
            },
            "rooting": rooting,
            "species_tree": str(species_tree_path) if species_tree_path else None,
            "counts": {"selected": len(selected), "completed": len(summaries), "reused": reused},
            "output_checksums": {events_path.name: sha256_file(events_path)},
        },
    )
    return manifest_path, len(selected), reused
