"""Deterministic protein-family evolution fixtures with explicit ground truth."""

from __future__ import annotations

import csv
import random
from collections import defaultdict
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any

from ogprofiler.core.manifest import checksum_manifest, sha256_json, write_json
from ogprofiler.exceptions import InputError

AMINO_ACIDS = "ACDEFGHIKLMNPQRSTVWY"
SPECIES = ("species_00", "species_01", "species_02", "species_03")


@dataclass(frozen=True, slots=True)
class SyntheticScenario:
    divergence: float = 0.10
    duplication_rate: float = 0.10
    loss_rate: float = 0.05
    expansion: int = 1
    fusion_rate: float = 0.0
    replicate: int = 0
    seed: int = 42

    def __post_init__(self) -> None:
        for name in ("divergence", "duplication_rate", "loss_rate", "fusion_rate"):
            value = float(getattr(self, name))
            if not 0.0 <= value <= 1.0:
                raise InputError(f"{name} must be between 0 and 1")
        if self.expansion < 1:
            raise InputError("expansion must be at least 1")
        if self.replicate < 0:
            raise InputError("replicate must be non-negative")

    @property
    def scenario_id(self) -> str:
        identity = sha256_json(asdict(self))[:12]
        return f"synthetic-{identity}"


DEFAULT_SCENARIO_AXES: dict[str, tuple[float | int, ...]] = {
    "divergence": (0.05, 0.20, 0.40),
    "duplication_rate": (0.0, 0.25, 0.50),
    "loss_rate": (0.0, 0.20, 0.40),
    "expansion": (1, 2, 4),
    "fusion_rate": (0.0, 0.10, 0.25),
}


def generate_scenario_matrix(
    *,
    replicates: int = 1,
    seed: int = 42,
    axes: dict[str, tuple[float | int, ...]] | None = None,
) -> tuple[SyntheticScenario, ...]:
    """Create a deterministic full-factorial scenario matrix."""
    if replicates < 1:
        raise InputError("replicates must be at least 1")
    selected = axes or DEFAULT_SCENARIO_AXES
    required = tuple(DEFAULT_SCENARIO_AXES)
    if tuple(selected) != required:
        raise InputError(f"Synthetic axes must be ordered as {required}")
    scenarios: list[SyntheticScenario] = []
    for divergence in selected["divergence"]:
        for duplication_rate in selected["duplication_rate"]:
            for loss_rate in selected["loss_rate"]:
                for expansion in selected["expansion"]:
                    for fusion_rate in selected["fusion_rate"]:
                        for replicate in range(replicates):
                            scenarios.append(
                                SyntheticScenario(
                                    divergence=float(divergence),
                                    duplication_rate=float(duplication_rate),
                                    loss_rate=float(loss_rate),
                                    expansion=int(expansion),
                                    fusion_rate=float(fusion_rate),
                                    replicate=replicate,
                                    seed=seed,
                                )
                            )
    return tuple(scenarios)


def write_scenario_matrix(
    output: Path, scenarios: tuple[SyntheticScenario, ...]
) -> tuple[Path, Path]:
    output.mkdir(parents=True, exist_ok=True)
    json_path = output / "scenario-matrix.json"
    tsv_path = output / "scenario-matrix.tsv"
    write_json(
        json_path,
        {
            "schema_version": "synthetic-scenario-matrix-v1",
            "design": "full-factorial",
            "scenario_count": len(scenarios),
            "scenarios": [
                {"scenario_id": item.scenario_id, **asdict(item)} for item in scenarios
            ],
        },
    )
    fields = ["scenario_id", *asdict(scenarios[0])] if scenarios else ["scenario_id"]
    with tsv_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for item in scenarios:
            writer.writerow({"scenario_id": item.scenario_id, **asdict(item)})
    return json_path, tsv_path


@dataclass(frozen=True, slots=True)
class _Gene:
    species: str
    protein_id: str
    family: str
    ancestral_family: str
    lineage: int
    sequence: str
    role: str


def _random_sequence(rng: random.Random, length: int) -> str:
    return "".join(rng.choice(AMINO_ACIDS) for _ in range(length))


def _mutate(sequence: str, probability: float, rng: random.Random) -> str:
    residues: list[str] = []
    for residue in sequence:
        if rng.random() >= probability:
            residues.append(residue)
            continue
        options = AMINO_ACIDS.replace(residue, "")
        residues.append(rng.choice(options))
    return "".join(residues)


def _node(
    rows: list[dict[str, Any]],
    *,
    family: str,
    node_id: str,
    parent_id: str,
    event: str,
    species: str = "",
    protein_id: str = "",
    lineage: int,
) -> None:
    rows.append(
        {
            "family": family,
            "node_id": node_id,
            "parent_id": parent_id,
            "event": event,
            "species": species,
            "protein_id": protein_id,
            "lineage": lineage,
        }
    )


def _write_tsv(path: Path, rows: list[dict[str, Any]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def _write_fasta(path: Path, genes: list[_Gene]) -> None:
    with path.open("w", encoding="utf-8", newline="\n") as handle:
        for gene in sorted(genes, key=lambda item: item.protein_id):
            handle.write(f">{gene.protein_id}\n")
            for start in range(0, len(gene.sequence), 80):
                handle.write(gene.sequence[start : start + 80] + "\n")


def generate_synthetic_dataset(
    output: Path,
    scenario: SyntheticScenario,
    *,
    ancestral_families: int = 12,
    sequence_length: int = 120,
) -> Path:
    """Generate FASTA plus family, genealogy, event, orthology, and domain truth."""
    if ancestral_families < 1:
        raise InputError("ancestral_families must be at least 1")
    if sequence_length < 40:
        raise InputError("sequence_length must be at least 40")
    output.mkdir(parents=True, exist_ok=True)
    proteomes = output / "proteomes"
    proteomes.mkdir(exist_ok=True)
    scenario_seed = int(sha256_json(asdict(scenario))[:16], 16) ^ scenario.seed
    rng = random.Random(scenario_seed)
    genes: list[_Gene] = []
    genealogy: list[dict[str, Any]] = []

    for family_index in range(ancestral_families):
        ancestral_family = f"AF{family_index:05d}"
        ancestor = _random_sequence(rng, sequence_length)
        root_id = f"{ancestral_family}.root"
        _node(
            genealogy,
            family=ancestral_family,
            node_id=root_id,
            parent_id="",
            event="ANCESTRAL_FAMILY",
            lineage=-1,
        )
        for lineage in range(scenario.expansion):
            true_family = f"{ancestral_family}.L{lineage:02d}"
            lineage_root = f"{true_family}.lineage"
            _node(
                genealogy,
                family=true_family,
                node_id=lineage_root,
                parent_id=root_id,
                event="DUPLICATION" if scenario.expansion > 1 else "LINEAGE",
                lineage=lineage,
            )
            lineage_sequence = _mutate(
                ancestor, scenario.divergence * 0.35 * lineage, rng
            )
            clade_sequences = {
                "AB": _mutate(lineage_sequence, scenario.divergence * 0.5, rng),
                "CD": _mutate(lineage_sequence, scenario.divergence * 0.5, rng),
            }
            for clade, species_indices in (("AB", (0, 1)), ("CD", (2, 3))):
                clade_node = f"{true_family}.{clade}"
                _node(
                    genealogy,
                    family=true_family,
                    node_id=clade_node,
                    parent_id=lineage_root,
                    event="SPECIATION",
                    lineage=lineage,
                )
                for species_index in species_indices:
                    species = SPECIES[species_index]
                    if rng.random() < scenario.loss_rate:
                        loss_node = f"{true_family}.{species}.loss"
                        _node(
                            genealogy,
                            family=true_family,
                            node_id=loss_node,
                            parent_id=clade_node,
                            event="LOSS",
                            species=species,
                            lineage=lineage,
                        )
                        continue
                    species_sequence = _mutate(
                        clade_sequences[clade], scenario.divergence * 0.5, rng
                    )
                    copies = 2 if rng.random() < scenario.duplication_rate else 1
                    parent_id = clade_node
                    if copies == 2:
                        parent_id = f"{true_family}.{species}.dup"
                        _node(
                            genealogy,
                            family=true_family,
                            node_id=parent_id,
                            parent_id=clade_node,
                            event="DUPLICATION",
                            species=species,
                            lineage=lineage,
                        )
                    for copy in range(copies):
                        protein_id = (
                            f"{species}_{ancestral_family}_L{lineage:02d}_C{copy:02d}"
                        )
                        sequence = _mutate(
                            species_sequence, scenario.divergence * 0.15 * copy, rng
                        )
                        genes.append(
                            _Gene(
                                species,
                                protein_id,
                                true_family,
                                ancestral_family,
                                lineage,
                                sequence,
                                "terminal_duplicate" if copies == 2 else "single_copy",
                            )
                        )
                        _node(
                            genealogy,
                            family=true_family,
                            node_id=f"{true_family}.{species}.gene{copy:02d}",
                            parent_id=parent_id,
                            event="GENE",
                            species=species,
                            protein_id=protein_id,
                            lineage=lineage,
                        )

    domains: list[dict[str, Any]] = []
    genes_by_family: dict[str, list[_Gene]] = defaultdict(list)
    for gene in genes:
        genes_by_family[gene.family].append(gene)
    mutable = list(genes)
    candidates = list(range(len(mutable)))
    rng.shuffle(candidates)
    fusion_count = round(len(mutable) * scenario.fusion_rate)
    family_names = sorted(genes_by_family)
    for index in candidates[:fusion_count]:
        gene = mutable[index]
        choices = [name for name in family_names if name != gene.family]
        if not choices:
            break
        secondary_family = rng.choice(choices)
        donor = rng.choice(genes_by_family[secondary_family])
        breakpoint = len(gene.sequence)
        fused_sequence = gene.sequence + donor.sequence[len(donor.sequence) // 2 :]
        mutable[index] = _Gene(
            gene.species,
            gene.protein_id,
            gene.family,
            gene.ancestral_family,
            gene.lineage,
            fused_sequence,
            "fusion",
        )
        domains.append(
            {
                "protein_id": gene.protein_id,
                "primary_family": gene.family,
                "secondary_family": secondary_family,
                "fusion_breakpoint": breakpoint,
            }
        )
    genes = mutable

    for species in SPECIES:
        _write_fasta(proteomes / f"{species}.faa", [g for g in genes if g.species == species])
    truth = [
        {
            "species": gene.species,
            "protein_id": gene.protein_id,
            "family": gene.family,
            "role": gene.role,
            "length": len(gene.sequence),
        }
        for gene in sorted(genes, key=lambda item: (item.species, item.protein_id))
    ]
    _write_tsv(
        output / "ground_truth.tsv",
        truth,
        ["species", "protein_id", "family", "role", "length"],
    )
    _write_tsv(
        output / "genealogy.tsv",
        genealogy,
        ["family", "node_id", "parent_id", "event", "species", "protein_id", "lineage"],
    )
    event_rows = [
        row for row in genealogy if row["event"] in {"SPECIATION", "DUPLICATION", "LOSS"}
    ]
    _write_tsv(
        output / "true_events.tsv",
        event_rows,
        ["family", "node_id", "parent_id", "event", "species", "protein_id", "lineage"],
    )
    orthologs: list[dict[str, Any]] = []
    by_family: dict[str, list[_Gene]] = defaultdict(list)
    for gene in genes:
        by_family[gene.family].append(gene)
    for family, members in sorted(by_family.items()):
        for left_index, left in enumerate(members):
            for right in members[left_index + 1 :]:
                if left.species == right.species:
                    continue
                first, second = sorted((left.protein_id, right.protein_id))
                orthologs.append(
                    {"protein_a": first, "protein_b": second, "family": family}
                )
    orthologs.sort(key=lambda row: (row["protein_a"], row["protein_b"]))
    _write_tsv(
        output / "true_orthologs.tsv",
        orthologs,
        ["protein_a", "protein_b", "family"],
    )
    _write_tsv(
        output / "domain_architecture.tsv",
        domains,
        ["protein_id", "primary_family", "secondary_family", "fusion_breakpoint"],
    )
    (output / "species_tree.nwk").write_text(
        "((species_00,species_01),(species_02,species_03));\n", encoding="utf-8"
    )
    artifacts = sorted(
        [
            output / "ground_truth.tsv",
            output / "genealogy.tsv",
            output / "true_events.tsv",
            output / "true_orthologs.tsv",
            output / "domain_architecture.tsv",
            output / "species_tree.nwk",
            *proteomes.glob("*.faa"),
        ]
    )
    manifest_path = output / "manifest.json"
    write_json(
        manifest_path,
        {
            "schema_version": "synthetic-evolution-dataset-v1",
            "scenario_id": scenario.scenario_id,
            "scenario": asdict(scenario),
            "species_tree": "((species_00,species_01),(species_02,species_03));",
            "ancestral_families": ancestral_families,
            "sequence_length": sequence_length,
            "protein_count": len(genes),
            "terminal_family_count": len(by_family),
            "genealogy_node_count": len(genealogy),
            "true_event_count": len(event_rows),
            "true_ortholog_pair_count": len(orthologs),
            "fusion_count": len(domains),
            "artifact_checksums": checksum_manifest(artifacts, output),
        },
    )
    return manifest_path
