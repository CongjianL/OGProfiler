from __future__ import annotations

import json
from pathlib import Path

import pyarrow.parquet as pq

from ogprofiler.cli import main


def _write_inputs(directory: Path, reverse_records: bool = False) -> None:
    directory.mkdir()
    first = [("z", "AAAA"), ("a", "CCCCC")]
    if reverse_records:
        first.reverse()
    (directory / "zeta.faa").write_text(
        "".join(f">{name}\n{sequence}\n" for name, sequence in first),
        encoding="utf-8",
    )
    (directory / "alpha.faa").write_text(">p2\nMMMM\n>p1\nGGG\n", encoding="utf-8")


def _metadata(run: Path) -> tuple[list[dict[str, object]], list[dict[str, object]], str, dict]:
    species = pq.read_table(run / "input/species.parquet").to_pylist()
    proteins = pq.read_table(run / "input/proteins.parquet").to_pylist()
    fasta = (run / "input/proteins.faa").read_text(encoding="utf-8")
    manifest = json.loads((run / "manifest.json").read_text(encoding="utf-8"))
    manifest.pop("command")
    return species, proteins, fasta, manifest


def test_prepare_is_deterministic_across_record_order(tmp_path: Path) -> None:
    inputs_a = tmp_path / "inputs-a"
    inputs_b = tmp_path / "inputs-b"
    _write_inputs(inputs_a)
    _write_inputs(inputs_b, reverse_records=True)
    run_a = tmp_path / "run-a"
    run_b = tmp_path / "run-b"

    assert main(["prepare", "--proteomes", str(inputs_a), "--out", str(run_a)]) == 0
    assert main(["prepare", "--proteomes", str(inputs_b), "--out", str(run_b)]) == 0

    metadata_a = _metadata(run_a)
    metadata_b = _metadata(run_b)
    assert metadata_a[:3] == metadata_b[:3]
    assert metadata_a[3]["dataset"]["n_species"] == 2
    assert metadata_a[3]["dataset"]["n_proteins"] == 4
    assert metadata_a[3]["dataset"]["id_algorithm"] == "sorted-path+sorted-original-id-v1"
    assert (run_a / "run.db").is_file()
    assert (run_a / "run.yaml").is_file()
