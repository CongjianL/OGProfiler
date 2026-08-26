from __future__ import annotations

import json
import shutil
from pathlib import Path

import pyarrow as pa
import pyarrow.parquet as pq
import pytest

from ogprofiler.exceptions import PhylogenyError
from ogprofiler.phylogeny.backends import BackendRun, FastTreeBackend, MafftBackend
from ogprofiler.phylogeny.newick import (
    leaf_names,
    midpoint_root,
    outgroup_root,
    parse_newick,
    prune_tree,
    to_newick,
)
from ogprofiler.phylogeny.reconciliation import LcaReconciliationBackend
from ogprofiler.phylogeny.selection import select_refinement_families
from ogprofiler.phylogeny.stage import (
    classify_event_conflict,
    run_phylogenetic_refinement_stage,
)


class FakeAlignment:
    name = "fake-alignment"

    def version(self) -> str:
        return "fake-align 1"

    def align(self, input_fasta: Path, output_fasta: Path, threads: int) -> BackendRun:
        shutil.copyfile(input_fasta, output_fasta)
        return BackendRun(("fake-align", str(threads)), self.version(), "", "")


class FakeTree:
    name = "fake-tree"

    def version(self) -> str:
        return "fake-tree 1"

    def infer(self, alignment_fasta: Path, output_newick: Path) -> BackendRun:
        tree = "(P000000000000:1,(P000000000001:1,P000000000002:1):1);\n"
        output_newick.write_text(tree, encoding="utf-8")
        return BackendRun(("fake-tree", str(alignment_fasta)), self.version(), tree, "")


def _write_tsv(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")


def _fixture(root: Path) -> None:
    (root / "input").mkdir(parents=True)
    pq.write_table(
        pa.table(
            {
                "protein_id": [0, 1, 2],
                "species_id": [0, 1, 2],
                "original_id": ["a", "b", "c"],
                "length": [4, 4, 4],
            }
        ),
        root / "input/proteins.parquet",
    )
    pq.write_table(
        pa.table(
            {
                "species_id": [0, 1, 2],
                "species_name": ["species_0", "species_1", "species_2"],
                "source_file": ["0.faa", "1.faa", "2.faa"],
            }
        ),
        root / "input/species.parquet",
    )
    (root / "input/proteins.faa").write_text(
        ">OGP2P000000000000\nAAAA\n"
        ">OGP2P000000000001\nAAAA\n"
        ">OGP2P000000000002\nAAAA\n",
        encoding="utf-8",
    )
    _write_tsv(
        root / "results/families.tsv",
        "family_id\tcomponent_id\tcluster_id\tn_genes\tn_species\tterminal_reason\t"
        "network_event\nOG000000000\t0\t2\t3\t3\tNO_SPLIT\tAMBIGUOUS\n",
    )
    _write_tsv(
        root / "results/members.tsv",
        "family_id\tprotein_id\tspecies_id\toriginal_id\n"
        "OG000000000\t0\t0\ta\nOG000000000\t1\t1\tb\nOG000000000\t2\t2\tc\n",
    )
    _write_tsv(
        root / "results/hierarchy.tsv",
        "cluster_id\tparent_id\tcomponent_id\tdepth\tn_genes\tn_species\tresolution\t"
        "quality\tchild_count\tterminal_reason\n"
        "0\t\t0\t0\t3\t3\t1\t1\t2\t\n"
        "1\t0\t0\t1\t1\t1\t\t\t0\tONE_SPECIES\n"
        "2\t0\t0\t1\t3\t3\t\t\t0\tNO_SPLIT\n",
    )
    _write_tsv(
        root / "results/events.tsv",
        "component_id\tcluster_id\tnetwork_event\toverlap_score\tconfidence\n"
        "0\t0\tDUPLICATION_LIKE\t1\t1\n0\t1\tSPECIES_SPECIFIC\t0\t1\n"
        "0\t2\tAMBIGUOUS\t0\t0\n",
    )


def test_newick_parse_prune_and_rooting_preserve_leaf_identity() -> None:
    tree = parse_newick("(('a one':1,b:2):3,c:4);")
    assert set(leaf_names(tree)) == {"a one", "b", "c"}
    assert set(leaf_names(parse_newick(to_newick(tree)))) == {"a one", "b", "c"}
    assert set(leaf_names(prune_tree(tree, {"a one", "c"}))) == {"a one", "c"}
    assert set(leaf_names(midpoint_root(tree))) == {"a one", "b", "c"}
    assert set(leaf_names(outgroup_root(tree, "c"))) == {"a one", "b", "c"}


def test_lca_reconciliation_separates_duplication_and_speciation() -> None:
    gene = parse_newick("((a1:1,a2:1):1,b1:1);")
    species = parse_newick("(A:1,B:1);")
    result = LcaReconciliationBackend().annotate(
        gene, species, {"a1": "A", "a2": "A", "b1": "B"}
    )
    assert result.root_event == "SPECIATION"
    assert result.duplication_count == 1
    assert [event.phylo_event for event in result.events] == ["DUPLICATION", "SPECIATION"]
    assert classify_event_conflict("DUPLICATION_LIKE", "SPECIATION") == "CONFLICT"
    assert classify_event_conflict("SPECIATION_LIKE", "SPECIATION") == "CONCORDANT"
    with pytest.raises(PhylogenyError, match="missing selected species"):
        LcaReconciliationBackend().annotate(
            gene, parse_newick("(A:1,C:1);"), {"a1": "A", "a2": "A", "b1": "B"}
        )


def test_external_backend_discovery_versions_and_standardized_subprocess(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    binary_root = tmp_path / "bin"
    binary_root.mkdir()
    mafft = binary_root / "mafft-test"
    mafft.write_text(
        "#!/bin/bash\n"
        "if [[ $1 == --version ]]; then echo 'MAFFT v-test' >&2; exit 0; fi\n"
        "/bin/cat \"${!#}\"\n",
        encoding="utf-8",
    )
    fasttree = binary_root / "fasttree-test"
    fasttree.write_text(
        "#!/bin/bash\n"
        "if [[ $# == 0 ]]; then echo 'FastTree Version test' >&2; exit 0; fi\n"
        "echo '(P000000000000:1,P000000000001:1);'\n",
        encoding="utf-8",
    )
    mafft.chmod(0o755)
    fasttree.chmod(0o755)
    monkeypatch.setenv("PATH", str(binary_root))
    input_fasta = tmp_path / "input.faa"
    input_fasta.write_text(">P000000000000\nAA\n>P000000000001\nAA\n", encoding="utf-8")
    alignment = tmp_path / "alignment.faa"
    alignment_run = MafftBackend(mafft.name).align(input_fasta, alignment, 2)
    assert alignment_run.version == "MAFFT v-test"
    assert alignment.read_text() == input_fasta.read_text()
    tree_path = tmp_path / "tree.nwk"
    tree_run = FastTreeBackend(fasttree.name).infer(alignment, tree_path)
    assert tree_run.version == "FastTree Version test"
    assert tree_path.read_text().endswith(";\n")
    with pytest.raises(PhylogenyError, match="Executable not found"):
        MafftBackend("missing-mafft").version()


def test_selection_and_refinement_preserve_network_and_phylo_evidence(tmp_path: Path) -> None:
    _fixture(tmp_path)
    selected = select_refinement_families(
        tmp_path,
        explicit_family_ids=(),
        selection_events={"AMBIGUOUS", "MIXED"},
        large_family_size=100,
        max_families=10,
    )
    assert len(selected) == 1
    assert selected[0].selection_reasons == ("AMBIGUOUS", "DUPLICATION_RICH")
    species_tree = tmp_path / "species.nwk"
    species_tree.write_text("((species_0,species_1),species_2,extra);\n", encoding="utf-8")
    kwargs = {
        "explicit_family_ids": (),
        "selection_events": {"AMBIGUOUS", "MIXED"},
        "large_family_size": 100,
        "max_families": 10,
        "rooting": "species-tree-aware",
        "outgroup": None,
        "species_tree_path": species_tree,
        "threads": 1,
        "alignment_backend": FakeAlignment(),
        "tree_backend": FakeTree(),
        "reconciliation_backend": LcaReconciliationBackend(),
    }
    manifest, count, reused = run_phylogenetic_refinement_stage(
        tmp_path, ["ogprofiler", "annotate"], **kwargs
    )
    assert (count, reused) == (1, 0)
    assert manifest.is_file()
    event_text = (
        tmp_path / "evolution/phylogenetic/phylogenetic-events.tsv"
    ).read_text(encoding="utf-8")
    assert "network_event\tphylo_event" in event_text
    assert "AMBIGUOUS\tSPECIATION" in event_text
    family_root = tmp_path / "evolution/phylogenetic/family=OG000000000"
    family_manifest = json.loads(
        (family_root / "phylogeny-manifest.json").read_text(encoding="utf-8")
    )
    assert family_manifest["backend_versions"]["alignment"] == "fake-align 1"
    assert "extra" not in (family_root / "species_tree.pruned.nwk").read_text()
    _, _, reused = run_phylogenetic_refinement_stage(
        tmp_path, ["ogprofiler", "annotate"], **kwargs
    )
    assert reused == 1
    (family_root / "gene_tree.rooted.nwk").write_bytes(b"corrupt")
    _, _, reused = run_phylogenetic_refinement_stage(
        tmp_path, ["ogprofiler", "annotate"], **kwargs
    )
    assert reused == 0
