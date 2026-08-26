#!/usr/bin/env python3
"""Generate deterministic Phase 0 datasets A-E and their ground truth."""

from __future__ import annotations

import argparse
import hashlib
import json
import random
import shutil
import tempfile
from dataclasses import dataclass
from pathlib import Path

AMINO_ACIDS = "ACDEFGHIKLMNPQRSTVWY"
SEED = 20260825
ROOT = Path(__file__).resolve().parent / "datasets"
SPECIES = ("species_01", "species_02", "species_03", "species_04")


@dataclass(frozen=True)
class Entry:
    species: str
    protein_id: str
    family: str
    role: str
    sequence: str


def external_id(entry: Entry) -> str:
    """Return a globally unique identifier compatible with frozen V1 output code."""

    return f"{entry.species}|{entry.protein_id}"


def sequence_for(label: str, length: int = 180) -> str:
    rng = random.Random(f"{SEED}:{label}")
    return "".join(rng.choice(AMINO_ACIDS) for _ in range(length))


def mutate(sequence: str, label: str, fraction: float) -> str:
    rng = random.Random(f"{SEED}:mutation:{label}")
    result = list(sequence)
    count = max(1, round(len(result) * fraction))
    for position in rng.sample(range(len(result)), min(count, len(result))):
        choices = AMINO_ACIDS.replace(result[position], "")
        result[position] = rng.choice(choices)
    return "".join(result)


def background_entries(dataset: str, family_count: int) -> list[Entry]:
    entries: list[Entry] = []
    for family_index in range(family_count):
        family = f"{dataset}_BG{family_index:03d}"
        ancestor = sequence_for(family)
        for species_index, species in enumerate(SPECIES):
            protein_id = f"{family.lower()}_{species_index + 1:02d}"
            entries.append(
                Entry(
                    species,
                    protein_id,
                    family,
                    "one_to_one",
                    mutate(ancestor, protein_id, 0.04 + species_index * 0.005),
                )
            )
    return entries


def dataset_a() -> tuple[list[Entry], dict[str, object]]:
    entries = background_entries("A", 10)
    for species_index, species in enumerate(SPECIES):
        for unique_index in range(2):
            protein_id = f"a_unique_{species_index + 1:02d}_{unique_index + 1:02d}"
            entries.append(
                Entry(
                    species,
                    protein_id,
                    f"A_UNIQUE_{protein_id}",
                    "species_specific",
                    sequence_for(protein_id),
                )
            )
    return entries, {
        "name": "A_small_sanity",
        "purpose": (
            "Complete-pipeline smoke test with ten one-to-one families and "
            "species-specific proteins."
        ),
        "expected": {"one_to_one_families": 10, "proteins_per_species": 12},
    }


def dataset_b() -> tuple[list[Entry], dict[str, object]]:
    entries = background_entries("B", 8)
    ancestor = sequence_for("B_PARALOG_ANCESTOR", 210)
    copies = (3, 2, 1, 1)
    for species_index, (species, copy_count) in enumerate(zip(SPECIES, copies, strict=True)):
        species_parent = mutate(ancestor, f"B_species_{species_index}", 0.05)
        for copy_index in range(copy_count):
            protein_id = f"b_paralog_s{species_index + 1:02d}_c{copy_index + 1:02d}"
            entries.append(
                Entry(
                    species,
                    protein_id,
                    "B_PARALOG",
                    "recent_duplication" if copy_count > 1 else "single_copy",
                    mutate(species_parent, protein_id, 0.015 * copy_index + 0.005),
                )
            )
        while sum(1 for entry in entries if entry.species == species) < 12:
            unique_index = sum(1 for entry in entries if entry.species == species)
            protein_id = f"b_unique_s{species_index + 1:02d}_{unique_index:02d}"
            entries.append(
                Entry(
                    species,
                    protein_id,
                    f"B_UNIQUE_{protein_id}",
                    "species_specific",
                    sequence_for(protein_id),
                )
            )
    return entries, {
        "name": "B_paralog",
        "purpose": "RBH/non-RBH/LRB behavior around recent species-specific duplications.",
        "expected": {"expanded_family": "B_PARALOG", "copy_counts": list(copies)},
    }


def dataset_c() -> tuple[list[Entry], dict[str, object]]:
    entries = background_entries("C", 5)
    root = sequence_for("C_EXPANSION_ROOT", 220)
    clades = [mutate(root, f"C_clade_{index}", 0.12) for index in range(3)]
    for species_index, species in enumerate(SPECIES):
        for clade_index, clade in enumerate(clades):
            for copy_index in range(5):
                protein_id = (
                    f"c_expand_s{species_index + 1:02d}_"
                    f"clade{clade_index + 1:02d}_c{copy_index + 1:02d}"
                )
                entries.append(
                    Entry(
                        species,
                        protein_id,
                        "C_EXPANSION",
                        f"clade_{clade_index + 1}",
                        mutate(clade, protein_id, 0.025 + copy_index * 0.005),
                    )
                )
    return entries, {
        "name": "C_gene_family_expansion",
        "purpose": "A 60-member expanded family with three planted sequence clades.",
        "expected": {
            "expanded_family_size": 60,
            "planted_subclades": 3,
            "proteins_per_species": 20,
        },
    }


def dataset_d() -> tuple[list[Entry], dict[str, object]]:
    entries = background_entries("D", 8)
    domain_a = sequence_for("D_DOMAIN_A", 110)
    domain_b = sequence_for("D_DOMAIN_B", 105)
    linker = "GGGGSGGGGS"
    for species_index, species in enumerate(SPECIES):
        a_id = f"d_domain_a_s{species_index + 1:02d}"
        b_id = f"d_domain_b_s{species_index + 1:02d}"
        fusion_id = f"d_fusion_s{species_index + 1:02d}"
        entries.extend(
            [
                Entry(species, a_id, "D_DOMAIN_A", "domain_only", mutate(domain_a, a_id, 0.04)),
                Entry(species, b_id, "D_DOMAIN_B", "domain_only", mutate(domain_b, b_id, 0.04)),
                Entry(
                    species,
                    fusion_id,
                    "D_FUSION",
                    "fusion_bridge",
                    mutate(domain_a, fusion_id + "a", 0.03)
                    + linker
                    + mutate(domain_b, fusion_id + "b", 0.03),
                ),
                Entry(
                    species,
                    f"d_unique_s{species_index + 1:02d}",
                    f"D_UNIQUE_{species_index + 1:02d}",
                    "species_specific",
                    sequence_for(f"D_unique_{species_index}"),
                ),
            ]
        )
    return entries, {
        "name": "D_fusion_multidomain",
        "purpose": "Coverage filtering and bridge sensitivity for domain-only and fusion hits.",
        "expected": {"fusion_proteins": 4, "domain_families": 2, "proteins_per_species": 12},
    }


def dataset_e() -> tuple[list[Entry], dict[str, object]]:
    entries: list[Entry] = []
    root = sequence_for("E_LARGE_CC_ROOT", 240)
    subfamily_roots = [mutate(root, f"E_subfamily_{index}", 0.08) for index in range(6)]
    for species_index, species in enumerate(SPECIES):
        for copy_index in range(60):
            subfamily = copy_index % len(subfamily_roots)
            protein_id = f"e_large_s{species_index + 1:02d}_p{copy_index + 1:03d}"
            entries.append(
                Entry(
                    species,
                    protein_id,
                    "E_LARGE_CC",
                    f"subfamily_{subfamily + 1}",
                    mutate(subfamily_roots[subfamily], protein_id, 0.04 + (copy_index % 5) * 0.01),
                )
            )
    return entries, {
        "name": "E_large_connected_component",
        "purpose": "Runtime and peak-RSS benchmark with a planted 240-protein connected component.",
        "expected": {"largest_component_size": 240, "planted_subfamilies": 6},
    }


DATASETS = {
    "A_small_sanity": dataset_a,
    "B_paralog": dataset_b,
    "C_gene_family_expansion": dataset_c,
    "D_fusion_multidomain": dataset_d,
    "E_large_connected_component": dataset_e,
}


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def reference_edges(entries: list[Entry]) -> list[tuple[int, int, float]]:
    index_by_id = {entry.protein_id: index for index, entry in enumerate(entries)}
    by_family: dict[str, list[Entry]] = {}
    for entry in entries:
        by_family.setdefault(entry.family, []).append(entry)
    edge_weights: dict[tuple[int, int], float] = {}

    def add(left: Entry, right: Entry, weight: float) -> None:
        pair = tuple(sorted((index_by_id[left.protein_id], index_by_id[right.protein_id])))
        if pair[0] != pair[1]:
            edge_weights[pair] = max(edge_weights.get(pair, 0.0), weight)

    for family, members in by_family.items():
        by_role: dict[str, list[Entry]] = {}
        for member in members:
            by_role.setdefault(member.role, []).append(member)
        structured = family in {"C_EXPANSION", "E_LARGE_CC"}
        if structured:
            for role_members in by_role.values():
                for left_index, left in enumerate(role_members):
                    for right in role_members[left_index + 1 :]:
                        add(left, right, 1.0)
            representatives = [
                sorted(group, key=lambda item: item.protein_id)[0]
                for group in by_role.values()
            ]
            representatives.sort(key=lambda item: item.role)
            for left, right in zip(representatives, representatives[1:], strict=False):
                add(left, right, 0.05)
        else:
            for left_index, left in enumerate(members):
                for right in members[left_index + 1 :]:
                    add(left, right, 1.0)

    if any(entry.family == "D_FUSION" for entry in entries):
        by_species = {
            species: [entry for entry in entries if entry.species == species]
            for species in SPECIES
        }
        for species_entries in by_species.values():
            fusion = next(entry for entry in species_entries if entry.family == "D_FUSION")
            domain_a = next(entry for entry in species_entries if entry.family == "D_DOMAIN_A")
            domain_b = next(entry for entry in species_entries if entry.family == "D_DOMAIN_B")
            add(fusion, domain_a, 0.35)
            add(fusion, domain_b, 0.35)
    return [(left, right, edge_weights[(left, right)]) for left, right in sorted(edge_weights)]


def write_reference_gml(path: Path, entries: list[Entry]) -> None:
    edges = reference_edges(entries)
    with path.open("w", encoding="utf-8", newline="\n") as handle:
        handle.write('Creator "OGProfiler Phase 0 deterministic generator"\n')
        handle.write("Version 1\n")
        handle.write("graph\n[\n  directed 0\n")
        for index, entry in enumerate(entries):
            handle.write("  node\n  [\n")
            handle.write(f"    id {index}\n")
            handle.write(f'    name "{external_id(entry)}"\n')
            handle.write(f'    species "{entry.species}"\n')
            handle.write(f'    family "{entry.family}"\n')
            handle.write(f'    role "{entry.role}"\n')
            handle.write("  ]\n")
        for source, target, weight in edges:
            handle.write("  edge\n  [\n")
            handle.write(f"    source {source}\n    target {target}\n    NBS {weight:.8g}\n")
            handle.write("  ]\n")
        handle.write("]\n")


def write_dataset(root: Path, directory_name: str) -> None:
    entries, description = DATASETS[directory_name]()
    target = root / directory_name
    proteomes = target / "proteomes"
    proteomes.mkdir(parents=True)
    for species in SPECIES:
        selected = sorted(
            (entry for entry in entries if entry.species == species),
            key=lambda entry: entry.protein_id,
        )
        with (proteomes / f"{species}.faa").open("w", encoding="utf-8", newline="\n") as handle:
            for entry in selected:
                handle.write(
                    f">{external_id(entry)} family={entry.family} role={entry.role}\n"
                )
                for start in range(0, len(entry.sequence), 80):
                    handle.write(entry.sequence[start : start + 80] + "\n")

    with (target / "ground_truth.tsv").open("w", encoding="utf-8", newline="\n") as handle:
        handle.write("species\tprotein_id\tfamily\trole\tlength\n")
        for entry in sorted(entries, key=lambda item: (item.species, item.protein_id)):
            handle.write(
                f"{entry.species}\t{external_id(entry)}\t{entry.family}\t{entry.role}\t"
                f"{len(entry.sequence)}\n"
            )
    write_reference_gml(target / "reference_ssn.gml", entries)

    files = sorted(path for path in target.rglob("*") if path.is_file())
    manifest = {
        **description,
        "generator_seed": SEED,
        "n_species": len(SPECIES),
        "n_proteins": len(entries),
        "files": {path.relative_to(target).as_posix(): sha256(path) for path in files},
    }
    (target / "manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )


def generate(root: Path) -> None:
    if root.exists():
        shutil.rmtree(root)
    root.mkdir(parents=True)
    for directory_name in DATASETS:
        write_dataset(root, directory_name)


def tree_digest(root: Path) -> str:
    digest = hashlib.sha256()
    for path in sorted(item for item in root.rglob("*") if item.is_file()):
        digest.update(path.relative_to(root).as_posix().encode())
        digest.update(b"\0")
        digest.update(path.read_bytes())
        digest.update(b"\0")
    return digest.hexdigest()


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    if args.check:
        with tempfile.TemporaryDirectory() as temporary:
            generated = Path(temporary) / "datasets"
            generate(generated)
            if not ROOT.exists() or tree_digest(ROOT) != tree_digest(generated):
                raise SystemExit("Phase 0 datasets are stale; regenerate without --check")
        print(f"Phase 0 datasets verified: {tree_digest(ROOT)}")
        return 0
    generate(ROOT)
    print(f"Generated Phase 0 datasets: {tree_digest(ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
