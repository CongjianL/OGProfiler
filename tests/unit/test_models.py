from __future__ import annotations

from ogprofiler.core.models import (
    ComponentSummary,
    DatasetManifest,
    HierarchyNode,
    Protein,
    RunManifest,
    Species,
    SplitResult,
)


def test_core_models_round_trip() -> None:
    models = [
        Protein(1, 2, "protein-a", 123),
        Species(2, "species-a", "species-a.faa"),
        ComponentSummary(3, 10, 12, 4),
        HierarchyNode(1, None, 3, 0, 10, 4),
        SplitResult((0, 0, 1), 0.5, 0.91, 2, 1, 2 / 3, True, "accepted"),
    ]
    for model in models:
        assert type(model).from_dict(model.to_dict()) == model


def test_run_manifest_nested_round_trip() -> None:
    dataset = DatasetManifest("abc", {"a.faa": "def"}, 1, 2)
    manifest = RunManifest("2.0", "prepare-v1", ["ogprofiler", "prepare"], 42, "cfg", dataset)
    assert RunManifest.from_dict(manifest.to_dict()) == manifest
