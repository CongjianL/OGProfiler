from pathlib import Path

import pytest

from benchmarks.qfo.preflight import audit, compare_subset_file


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
    assert result["subsets"]["bacteria"]["n_proteins_in_all"] == 2
    assert result["subsets"]["eukaryota"]["n_proteins_in_all"] == 1
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


def test_audit_records_different_subset_without_mixing_main_input(tmp_path):
    source = tmp_path / "source"
    fixture(source)
    main_before = (source / "all/UP000000001_1.fasta").read_bytes()
    (source / "bacteria/UP000000001_1.fasta").write_text(">sp|CHANGED|ONE\nACDE\n")
    result = audit(source, tmp_path / "out")
    assert result["input_ready"] and result["subset_audit_completed"]
    assert not result["subset_checks_passed"]
    assert result["all_proteins"] == 3
    assert result["subsets"]["bacteria"]["comparison_classifications"] == {
        "SEQUENCE_SET_DIFFERENCE": 1
    }
    assert (tmp_path / "out/input/all/UP000000001_1.fasta").read_bytes() == main_before


def test_format_and_description_differences_preserve_sequence_identity(tmp_path):
    main, subset = tmp_path / "main.fa", tmp_path / "subset.fa"
    main.write_text(">sp|A001|ONE first description\nAAUU\n>tr|A002|TWO\nACDX\n")
    subset.write_text(">tr|A002|TWO changed description\nAC\nDX\n>sp|A001|ONE\nAAUU\n")
    result = compare_subset_file(subset, main)
    assert result["classification"] == "FORMAT_OR_DESCRIPTION_ONLY"
    assert result["usable_as_identical_subcollection"]
    assert result["sha256"] != result["main_sha256"]


def test_changed_sequence_and_accessions_are_counted(tmp_path):
    main, subset = tmp_path / "main.fa", tmp_path / "subset.fa"
    main.write_text(">sp|A001|ONE\nAAUU\n>tr|A002|TWO\nACDX\n")
    subset.write_text(">sp|A001|ONE\nAAAA\n>tr|A003|THREE\nACDE\n")
    result = compare_subset_file(subset, main)
    assert result["classification"] == "SEQUENCE_SET_DIFFERENCE"
    assert result["changed_sequences"] == 1
    assert result["added_accessions"] == result["missing_accessions"] == 1
    assert not result["usable_as_identical_subcollection"]


def test_accession_matching_does_not_hide_original_id_changes(tmp_path):
    main, subset = tmp_path / "main.fa", tmp_path / "subset.fa"
    main.write_text(">sp|A001|ONE\nAAUU\n")
    subset.write_text(">tr|A001|RENAMED\nAAUU\n")
    result = compare_subset_file(subset, main)
    assert result["classification"] == "ID_DIFFERENCE"
    assert result["changed_original_ids"] == 1
    assert not result["usable_as_identical_subcollection"]


def test_bacteria_is_independent_collection_without_reading_all(tmp_path, monkeypatch):
    source = tmp_path / "bacteria"
    source.mkdir()
    (source / "UP000000001_1.fasta").write_text(">sp|A001|ONE\nACDE\n")
    # Neither all nor eukaryota directories exist. Cross-directory audit is forbidden.
    monkeypatch.setattr(
        "benchmarks.qfo.preflight.compare_subset_file",
        lambda *args: pytest.fail("Unexpected cross-directory read"),
    )
    result = audit(source, tmp_path / "out", collection="bacteria")
    assert result["input_ready"]
    assert result["collection"] == "bacteria"
    assert result["n_species"] == result["n_proteins"] == 1
    assert not result["subset_audit_completed"]
    assert result["subset_checks_passed"] is None
    assert result["subsets"] == {} and result["other_species"] == []
    assert (tmp_path / "out/input/bacteria/UP000000001_1.fasta").is_file()
    assert not (tmp_path / "out/input/all").exists()
    assert "all_proteins" not in result
