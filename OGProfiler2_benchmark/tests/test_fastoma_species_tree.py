from __future__ import annotations

import re
from pathlib import Path


TREE = Path(__file__).parents[1] / "02_configs/fastoma/open_orthobench_ncbi_taxonomy_tree.nwk"


def test_fastoma_species_tree_is_balanced_and_has_frozen_leaves():
    text = TREE.read_text(encoding="utf-8").strip()
    assert text.endswith(";")
    depth = 0
    for char in text[:-1]:
        if char == "(":
            depth += 1
        elif char == ")":
            depth -= 1
            assert depth >= 0
    assert depth == 0
    observed = set(re.findall(r"(?<=[(,])([^(),;]+?)(?=[,)])", text))
    expected = {
        "Caenorhabditis_elegans.WBcel235.pep.all",
        "Canis_familiaris.CanFam3.1.pep.all",
        "Ciona_intestinalis.KH.pep.all",
        "Danio_rerio.GRCz11.pep.all",
        "Drosophila_melanogaster.BDGP6.28.pep.all",
        "Gallus_gallus.GRCg6a.pep.all",
        "Homo_sapiens.GRCh38.pep.all",
        "Monodelphis_domestica.ASM229v1.pep.all",
        "Mus_musculus.GRCm38.pep.all",
        "Pan_troglodytes.Pan_tro_3.0.pep.all",
        "Rattus_norvegicus.Rnor_6.0.pep.all",
        "Tetraodon_nigroviridis.TETRAODON8.pep.all",
    }
    assert observed == expected
