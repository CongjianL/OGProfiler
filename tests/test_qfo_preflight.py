from pathlib import Path

import pytest

from benchmarks.qfo.preflight import audit


def fixture(root: Path):
    for d in ("all", "bacteria", "eukaryota"):
        (root / d).mkdir(parents=True)
    a = ">sp|A001|ONE\nAAUU\n>tr|A002|TWO\nACDX\n"
    b = ">sp|B001|THREE\nACDE\n"
    for d in ("all", "bacteria"):
        (root / d / "UP000000001_1.fasta").write_text(a)
    for d in ("all", "eukaryota"):
        (root / d / "UP000000002_2.fasta").write_text(b)


def test_audit_freezes_valid_input_and_subset_identity(tmp_path):
    source = tmp_path / "source"
    fixture(source)
    result = audit(source, tmp_path / "out")
    assert result["input_ready"]
    assert result["all_species"] == 2 and result["all_proteins"] == 3
    assert result["subsets"]["bacteria"]["n_proteins"] == 2
    assert result["subsets"]["eukaryota"]["n_proteins"] == 1
    assert not result["official_qfo_scoring_ready"]
    assert result["release"] == "UNIDENTIFIED"
    assert (tmp_path / "out/input/all/UP000000001_1.fasta").read_bytes() == (
        source / "all/UP000000001_1.fasta"
    ).read_bytes()
    assert len((tmp_path / "out/id-map.tsv").read_text().splitlines()) == 4
    with pytest.raises(ValueError, match="already exists"):
        audit(source, tmp_path / "out")


def test_audit_rejects_cross_species_duplicate_accession(tmp_path):
    source = tmp_path / "source"
    fixture(source)
    (source / "all/UP000000002_2.fasta").write_text(">tr|A001|OTHER\nACDE\n")
    with pytest.raises(ValueError, match="Cross-species duplicate"):
        audit(source, tmp_path / "out")
    assert (tmp_path / "out/preflight-failure.json").is_file()


def test_audit_rejects_subsets_differing_from_all(tmp_path):
    source = tmp_path / "source"
    fixture(source)
    (source / "bacteria/UP000000001_1.fasta").write_text(">sp|CHANGED|ONE\nACDE\n")
    with pytest.raises(ValueError, match="Subset differs"):
        audit(source, tmp_path / "out")
