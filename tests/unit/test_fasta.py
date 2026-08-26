from __future__ import annotations

from pathlib import Path

import pytest

from ogprofiler.exceptions import InputError
from ogprofiler.input.fasta import parse_fasta


def test_parse_and_normalize_fasta(tmp_path: Path) -> None:
    path = tmp_path / "species.faa"
    path.write_text(">b description\nac d\n>a\nMZX*\n", encoding="utf-8")
    records = list(parse_fasta(path))
    assert [(record.identifier, record.sequence) for record in records] == [
        ("b", "ACD"),
        ("a", "MZX*"),
    ]


def test_duplicate_protein_id_is_rejected(tmp_path: Path) -> None:
    path = tmp_path / "species.faa"
    path.write_text(">same\nAAA\n>same duplicate\nCCC\n", encoding="utf-8")
    with pytest.raises(InputError, match="Duplicate protein ID"):
        list(parse_fasta(path))


def test_empty_sequence_is_rejected(tmp_path: Path) -> None:
    path = tmp_path / "species.faa"
    path.write_text(">empty\n>next\nAAA\n", encoding="utf-8")
    with pytest.raises(InputError, match="Empty sequence"):
        list(parse_fasta(path))


def test_empty_proteome_is_rejected(tmp_path: Path) -> None:
    path = tmp_path / "species.faa"
    path.write_text("\n", encoding="utf-8")
    with pytest.raises(InputError, match="Empty proteome"):
        list(parse_fasta(path))


def test_illegal_residue_policy(tmp_path: Path) -> None:
    path = tmp_path / "species.faa"
    path.write_text(">p\nAA?C\n", encoding="utf-8")
    with pytest.raises(InputError, match="Illegal sequence"):
        list(parse_fasta(path))
    assert list(parse_fasta(path, "replace_with_x"))[0].sequence == "AAXC"
