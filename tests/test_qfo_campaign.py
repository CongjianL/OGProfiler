import csv
import json
import subprocess
import sys
from pathlib import Path

import pytest

from benchmarks.qfo.campaign import stage
from benchmarks.qfo.preflight import audit
from ogprofiler.input.fasta import parse_fasta


def origin_fixture(tmp_path):
    source = tmp_path / "bacteria"
    source.mkdir()
    for i in range(2):
        (source / f"species{i}.fasta").write_text(
            f">sp|P{i}|NAME{i} description\nACDE\n>tr|Q{i}|OTHER{i}\nAAAA\n"
        )
    origin = tmp_path / "origin"
    report = audit(source, origin / "preflight", collection="bacteria")
    (origin / "qfo-preflight-completion.json").write_text('{"preflight_completed":true}')
    return origin, report["dataset_sha256"]


def test_campaign_v1_adapter_roundtrip_preserves_sequences_and_species(tmp_path):
    origin, digest = origin_fixture(tmp_path)
    out = tmp_path / "campaign"
    summary = stage(origin, out, digest, mode="full")
    assert summary["n_species"] == 2 and summary["n_proteins"] == 4
    assert summary["dataset_sha256"] == digest
    with (out / "v1-id-map.tsv").open() as h:
        mapping = {r["adapted_id"]: r["original_id"] for r in csv.DictReader(h, delimiter="\t")}
    assert len(set(mapping.values())) == len(mapping) == 4
    for p in (out / "input").glob("*.fasta"):
        adapted = list(parse_fasta(out / "v1-input" / p.name))
        original = list(parse_fasta(p))
        assert [(mapping[r.identifier], r.sequence) for r in adapted] == [
            (r.identifier, r.sequence) for r in original
        ]
        assert all(len(r.identifier.split("|")) == 2 for r in adapted)
    primary = out / "fake-v1.tsv"
    primary.write_text("OG0\t2\t4\t" + " ".join(mapping) + "\n")
    converter = (
        Path(__file__).resolve().parents[1]
        / "OGProfiler2_benchmark/scripts/converters/convert_competitor_groups.py"
    )
    subprocess.run(
        [
            sys.executable,
            str(converter),
            "--tool",
            "ogprofiler_v1",
            "--input",
            str(primary),
            "--id-map",
            str(out / "v1-id-map.tsv"),
            "--fasta",
            str(out / "input"),
            "--out",
            str(out / "groups.tsv"),
        ],
        check=True,
    )
    with (out / "groups.tsv").open() as h:
        assert {r["protein_id"] for r in csv.DictReader(h, delimiter="\t")} == set(mapping.values())
    with pytest.raises(FileExistsError):
        stage(origin, out, digest, mode="full")


def test_campaign_rejects_modified_frozen_input(tmp_path):
    origin, digest = origin_fixture(tmp_path)
    p = origin / "preflight/input/bacteria/species0.fasta"
    p.write_text(p.read_text().replace("ACDE", "ACDF"))
    with pytest.raises(ValueError, match="Frozen input changed"):
        stage(origin, tmp_path / "campaign", digest, mode="full")


def test_campaign_rejects_all_origin(tmp_path):
    origin, digest = origin_fixture(tmp_path)
    p = origin / "preflight/qfo-input-audit.json"
    report = json.loads(p.read_text())
    report["collection"] = "all"
    p.write_text(json.dumps(report))
    with pytest.raises(ValueError, match="Unexpected collection"):
        stage(origin, tmp_path / "campaign", digest, mode="full")


def test_campaign_smoke_is_labeled_separately(tmp_path):
    origin, digest = origin_fixture(tmp_path)
    result = stage(origin, tmp_path / "smoke", digest, mode="smoke")
    assert result["mode"] == "smoke"
    assert not result["official_qfo_scores_produced"]
    assert result["full_dataset_sha256"] == digest


def test_campaign_rejects_modified_smoke_records(tmp_path):
    origin, digest = origin_fixture(tmp_path)
    p = origin / "preflight/smoke-input/species0.fasta"
    p.write_text(p.read_text().replace("ACDE", "ACDF"))
    with pytest.raises(ValueError, match="Smoke records changed"):
        stage(origin, tmp_path / "campaign", digest, mode="smoke")


def test_proteome_smoke_uses_first_three_full_species_and_own_digest(tmp_path):
    origin, _ = origin_fixture(tmp_path)
    source = origin / "preflight/input/bacteria"
    for i in (2, 3):
        (source / f"species{i}.fasta").write_text(f">sp|NEW{i}|NAME\nACDE\n")
    import hashlib

    from ogprofiler.core.manifest import sha256_file

    digest = hashlib.sha256(
        json.dumps(
            {p.name: sha256_file(p) for p in sorted(source.glob("*.fasta"))}, sort_keys=True
        ).encode()
    ).hexdigest()
    report_path = origin / "preflight/qfo-input-audit.json"
    report = json.loads(report_path.read_text())
    report.update(dataset_sha256=digest, n_species=4, n_proteins=6)
    report_path.write_text(json.dumps(report))
    result = stage(origin, tmp_path / "smoke", digest, mode="smoke-proteomes")
    assert result["n_species"] == 3 and result["n_proteins"] == 5
    assert result["dataset_sha256"] != digest
    assert result["full_dataset_sha256"] == digest
    assert set(result["input_sha256"]) == {"species0.fasta", "species1.fasta", "species2.fasta"}
    for p in (tmp_path / "smoke/input").glob("*.fasta"):
        assert p.read_bytes() == (source / p.name).read_bytes()
