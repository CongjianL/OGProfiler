from __future__ import annotations

from pathlib import Path

import pytest
import yaml

from ogprofiler.config import load_config
from ogprofiler.exceptions import InputError


def test_default_leiden_weights_use_mean_projection() -> None:
    assert load_config()["edges"]["symmetrization"] == "mean"


def test_yaml_and_cli_override(tmp_path: Path) -> None:
    config_path = tmp_path / "config.yaml"
    config_path.write_text("search:\n  threads: 12\nhierarchy:\n  seed: 9\n", encoding="utf-8")
    config = load_config(str(config_path), ["search.threads=3", "output.compression=gzip"])
    assert config["search"]["threads"] == 3
    assert config["hierarchy"]["seed"] == 9
    assert config["output"]["compression"] == "gzip"


def test_unknown_config_key_is_rejected(tmp_path: Path) -> None:
    config_path = tmp_path / "config.yaml"
    config_path.write_text("search:\n  mystery: true\n", encoding="utf-8")
    with pytest.raises(InputError, match="mystery"):
        load_config(str(config_path))


def test_invalid_coverage_is_rejected() -> None:
    with pytest.raises(InputError, match="coverage"):
        load_config(overrides=["edges.min_query_coverage=101"])


def test_section_type_mismatch_is_rejected(tmp_path: Path) -> None:
    config_path = tmp_path / "config.yaml"
    config_path.write_text("search: disabled\n", encoding="utf-8")
    with pytest.raises(InputError, match="wrong type"):
        load_config(str(config_path))


def test_search_configuration_is_explicit_and_validated() -> None:
    config = load_config(
        overrides=["search.sensitivity=very-sensitive", "search.max_target_seqs=0"]
    )
    assert config["search"]["sensitivity"] == "very-sensitive"
    assert config["search"]["max_target_seqs"] == 0
    with pytest.raises(InputError, match="max_target_seqs"):
        load_config(overrides=["search.max_target_seqs=-1"])


def test_edge_thresholds_are_validated() -> None:
    config = load_config(
        overrides=["edges.min_bidirectional_coverage=25", "edges.best_hit_tolerance=0.01"]
    )
    assert config["edges"]["min_bidirectional_coverage"] == 25
    with pytest.raises(InputError, match="best_hit_tolerance"):
        load_config(overrides=["edges.best_hit_tolerance=-0.1"])


def test_component_batch_configuration_is_validated() -> None:
    with pytest.raises(InputError, match="components.edge_batch_size"):
        load_config(overrides=["components.edge_batch_size=0"])


def test_default_configuration_exactly_matches_frozen_soft42_preset() -> None:
    preset = Path(__file__).resolve().parents[2] / "presets/embleya-soft42.yaml"
    assert load_config() == load_config(str(preset))
    assert load_config() == yaml.safe_load(preset.read_text(encoding="utf-8"))


def test_previous_kway_depth20_is_available_as_explicit_override() -> None:
    config = load_config(overrides=["hierarchy.topology_policy=kway_v1", "hierarchy.max_depth=20"])
    assert config["hierarchy"]["topology_policy"] == "kway_v1"
    assert config["hierarchy"]["max_depth"] == 20


def test_partial_configuration_inherits_soft42(tmp_path: Path) -> None:
    path = tmp_path / "partial.yaml"
    path.write_text("runtime:\n  workers: 3\n", encoding="utf-8")
    config = load_config(str(path))
    assert config["hierarchy"]["topology_policy"] == "soft_binary_24_v2"
    assert config["hierarchy"]["max_depth"] == 42
    assert config["runtime"]["workers"] == 3
