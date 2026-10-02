"""Validated YAML configuration with deterministic CLI overrides."""

from __future__ import annotations

from copy import deepcopy
from typing import TYPE_CHECKING, Any

if TYPE_CHECKING:
    from ogprofiler.hierarchy.engine import HierarchyConfig

import yaml

from ogprofiler.exceptions import InputError

DEFAULT_CONFIG: dict[str, dict[str, Any]] = {
    "input": {
        "extensions": [".faa", ".fa", ".fasta", ".fna"],
        "illegal_character_policy": "error",
    },
    "search": {
        "backend": "diamond",
        "executable": "diamond",
        "evalue": 1e-3,
        "threads": 8,
        "sensitivity": "more-sensitive",
        "max_target_seqs": 0,
        "max_hsps": 0,
        "mmseqs_sensitivity": "sensitive",
    },
    "similarity": {"normalization": "legacy_nbs", "nbs_fallback": "v1_zero"},
    "edges": {
        "method": "lrb",
        "apply_coverage_filter": False,
        "min_query_coverage": 0.0,
        "min_target_coverage": 0.0,
        "min_bidirectional_coverage": 0.0,
        "best_hit_tolerance": 1e-3,
        "symmetrization": "mean",
    },
    "components": {"edge_batch_size": 65_536, "max_open_files": 64},
    "hierarchy": {
        "method": "rber",
        "seed": 42,
        "min_family_size": None,
        "admission_policy": "nonempty_children_v1",
        "topology_policy": "kway_v1",
        "recursion_stop_size": 1,
        "max_candidate_evaluations": 24,
        "max_coarse_candidates": 10,
        "rescue_grid_points": 8,
        "component_leiden_call_budget": None,
        "leiden_iterations": 10,
        "max_depth": 20,
        "resolution_strategy": "bounded_adaptive_v2",
        "gamma_min": 0.01,
        "gamma_max": 10.0,
        "gamma_growth": 2.0,
        "local_grid_points": 5,
        "max_child_fraction": 0.95,
        "tiny_fragment_size": 2,
        "max_tiny_fragment_fraction": 1.0,
        "stability_threshold": 0.9,
        "publication_seeds": 5,
        "min_split_quality": None,
        "stability_mode": "robust",
        "subtree_workers": 1,
        "subtree_release_size": 50_000,
    },
    "evolution": {
        "network_overlap_threshold": 0.0,
        "phylogenetic_refinement": False,
    },
    "orthogroups": {"strategy": "v1_compatible", "species_overlap_count": 0, "refinement": False},
    "phylogeny": {
        "alignment_backend": "mafft",
        "alignment_executable": "mafft",
        "tree_backend": "fasttree",
        "tree_executable": "FastTree",
        "threads": 1,
        "rooting": "midpoint",
        "selection_events": ["AMBIGUOUS", "MIXED"],
        "large_family_size": 100,
        "max_families": 100,
    },
    "output": {
        "compression": "zstd",
        "emit_pairwise_orthologs": False,
        "ortholog_pair_chunk_size": 100_000,
    },
    "runtime": {"workers": 1, "component_retries": 1, "log_level": "INFO"},
}

_ALLOWED_BACKENDS = {"diamond", "mmseqs", "blastp"}
_ALLOWED_NORMALIZATION = {"legacy_nbs", "raw_bitscore", "length_scaled_bitscore"}
_ALLOWED_NBS_FALLBACK = {"v1_zero", "v2_max"}
_ALLOWED_DIAMOND_SENSITIVITY = {
    "fast",
    "mid-sensitive",
    "sensitive",
    "more-sensitive",
    "very-sensitive",
    "ultra-sensitive",
}
_ALLOWED_EDGE_METHODS = {"lrb", "rbh", "ar", "arb"}
_ALLOWED_SYMMETRIZATION = {"forward", "max", "min", "mean", "geometric_mean"}
_ALLOWED_HIERARCHY_METHODS = {"rber", "rbcv", "cpm", "modularity"}
_ALLOWED_RESOLUTION_STRATEGIES = {"adaptive", "log_grid", "bounded_adaptive_v2"}
_ALLOWED_STABILITY = {"fast", "robust", "publication"}
_ALLOWED_CHARACTER_POLICIES = {"error", "replace_with_x"}
_ALLOWED_COMPRESSION = {"none", "gzip", "snappy", "zstd"}
_ALLOWED_LOG_LEVELS = {"DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"}
_ALLOWED_ROOTING = {"midpoint", "species-tree-aware", "outgroup"}


def deep_merge(base: dict[str, Any], updates: dict[str, Any]) -> dict[str, Any]:
    result = deepcopy(base)
    for key, value in updates.items():
        if key not in result:
            raise InputError(f"Unknown configuration section or key: {key}")
        current_is_mapping = isinstance(result[key], dict)
        update_is_mapping = isinstance(value, dict)
        if current_is_mapping != update_is_mapping:
            raise InputError(f"Configuration value has the wrong type: {key}")
        if current_is_mapping:
            result[key] = deep_merge(result[key], value)
        else:
            result[key] = value
    return result


def _set_dotted(config: dict[str, Any], dotted_key: str, value: Any) -> None:
    parts = dotted_key.split(".")
    if len(parts) != 2 or parts[0] not in config or parts[1] not in config[parts[0]]:
        raise InputError(f"Unknown configuration override: {dotted_key}")
    config[parts[0]][parts[1]] = value


def parse_override(expression: str) -> tuple[str, Any]:
    if "=" not in expression:
        raise InputError(f"Configuration override must be KEY=VALUE: {expression}")
    key, raw_value = expression.split("=", 1)
    key = key.strip()
    if not key:
        raise InputError("Configuration override key is empty")
    return key, yaml.safe_load(raw_value)


def load_config(path: str | None = None, overrides: list[str] | None = None) -> dict[str, Any]:
    config = deepcopy(DEFAULT_CONFIG)
    iterations_explicit = False
    if path is not None:
        try:
            with open(path, encoding="utf-8") as handle:
                loaded = yaml.safe_load(handle) or {}
        except (OSError, yaml.YAMLError) as error:
            raise InputError(f"Failed to read configuration {path}: {error}") from error
        if not isinstance(loaded, dict):
            raise InputError("Configuration root must be a mapping")
        config = deep_merge(config, loaded)
        iterations_explicit = "leiden_iterations" in loaded.get("hierarchy", {})
    for expression in overrides or []:
        key, value = parse_override(expression)
        _set_dotted(config, key, value)
        iterations_explicit |= key == "hierarchy.leiden_iterations"
    # Preserve historical replay defaults, without rewriting explicit run budgets.
    if not iterations_explicit and config["hierarchy"]["admission_policy"] == "legacy_strict":
        config["hierarchy"]["leiden_iterations"] = 2
    validate_config(config)
    return config


def validate_config(config: dict[str, Any]) -> None:
    from ogprofiler.orthogroups.models import OrthogroupConfig

    OrthogroupConfig(**config["orthogroups"])

    def choice(section: str, key: str, allowed: set[str]) -> None:
        value = config[section][key]
        if value not in allowed:
            expected = ", ".join(sorted(allowed))
            raise InputError(f"{section}.{key} must be one of: {expected}")

    choice("input", "illegal_character_policy", _ALLOWED_CHARACTER_POLICIES)
    choice("search", "backend", _ALLOWED_BACKENDS)
    choice("search", "sensitivity", _ALLOWED_DIAMOND_SENSITIVITY)
    choice("search", "mmseqs_sensitivity", _ALLOWED_DIAMOND_SENSITIVITY)
    choice("similarity", "normalization", _ALLOWED_NORMALIZATION)
    choice("similarity", "nbs_fallback", _ALLOWED_NBS_FALLBACK)
    choice("edges", "method", _ALLOWED_EDGE_METHODS)
    choice("edges", "symmetrization", _ALLOWED_SYMMETRIZATION)
    choice("hierarchy", "topology_policy", {"kway_v1", "soft_binary_24_v2"})
    choice("hierarchy", "method", _ALLOWED_HIERARCHY_METHODS)
    choice("hierarchy", "resolution_strategy", _ALLOWED_RESOLUTION_STRATEGIES)
    choice("hierarchy", "stability_mode", _ALLOWED_STABILITY)
    choice("output", "compression", _ALLOWED_COMPRESSION)
    choice("phylogeny", "rooting", _ALLOWED_ROOTING)

    if config["phylogeny"]["alignment_backend"] != "mafft":
        raise InputError("phylogeny.alignment_backend must be mafft")
    if config["phylogeny"]["tree_backend"] != "fasttree":
        raise InputError("phylogeny.tree_backend must be fasttree")
    for key in ("alignment_executable", "tree_executable"):
        if not isinstance(config["phylogeny"][key], str) or not config["phylogeny"][key]:
            raise InputError(f"phylogeny.{key} must be a non-empty string")
    for key in ("threads", "large_family_size", "max_families"):
        value = config["phylogeny"][key]
        if not isinstance(value, int) or isinstance(value, bool) or value < 1:
            raise InputError(f"phylogeny.{key} must be a positive integer")
    selection_events = config["phylogeny"]["selection_events"]
    if not isinstance(selection_events, list) or not all(
        isinstance(value, str) for value in selection_events
    ):
        raise InputError("phylogeny.selection_events must be a list of event names")

    if not isinstance(config["output"]["emit_pairwise_orthologs"], bool):
        raise InputError("output.emit_pairwise_orthologs must be boolean")
    chunk_size = config["output"]["ortholog_pair_chunk_size"]
    if not isinstance(chunk_size, int) or isinstance(chunk_size, bool) or chunk_size < 1:
        raise InputError("output.ortholog_pair_chunk_size must be a positive integer")

    extensions = config["input"]["extensions"]
    if (
        not isinstance(extensions, list)
        or not extensions
        or not all(isinstance(item, str) and item.startswith(".") for item in extensions)
    ):
        raise InputError("input.extensions must be a non-empty list of dot-prefixed suffixes")
    if not isinstance(config["search"]["threads"], int) or config["search"]["threads"] < 1:
        raise InputError("search.threads must be a positive integer")
    if (
        not isinstance(config["search"]["executable"], str)
        or not config["search"]["executable"].strip()
    ):
        raise InputError("search.executable must be a non-empty string")
    if (
        not isinstance(config["search"]["max_target_seqs"], int)
        or config["search"]["max_target_seqs"] < 0
    ):
        raise InputError("search.max_target_seqs must be a non-negative integer")
    if (
        not isinstance(config["search"]["max_hsps"], int)
        or isinstance(config["search"]["max_hsps"], bool)
        or config["search"]["max_hsps"] < 0
    ):
        raise InputError("search.max_hsps must be a non-negative integer")
    if not isinstance(config["runtime"]["workers"], int) or config["runtime"]["workers"] < 1:
        raise InputError("runtime.workers must be a positive integer")
    if (
        not isinstance(config["runtime"]["component_retries"], int)
        or isinstance(config["runtime"]["component_retries"], bool)
        or config["runtime"]["component_retries"] < 0
    ):
        raise InputError("runtime.component_retries must be a non-negative integer")
    for key in ("edge_batch_size", "max_open_files"):
        value = config["components"][key]
        if not isinstance(value, int) or isinstance(value, bool) or value < 1:
            raise InputError(f"components.{key} must be a positive integer")
    if not isinstance(config["hierarchy"]["seed"], int):
        raise InputError("hierarchy.seed must be an integer")
    hierarchy_config(config["hierarchy"])
    if (
        not isinstance(config["hierarchy"]["max_depth"], int)
        or config["hierarchy"]["max_depth"] < 1
    ):
        raise InputError("hierarchy.max_depth must be a positive integer")
    for key in ("subtree_workers", "subtree_release_size"):
        value = config["hierarchy"][key]
        if not isinstance(value, int) or isinstance(value, bool) or value < 1:
            raise InputError(f"hierarchy.{key} must be a positive integer")
    for key in ("gamma_min", "gamma_max", "gamma_growth", "max_child_fraction"):
        value = config["hierarchy"][key]
        if not isinstance(value, (int, float)) or isinstance(value, bool):
            raise InputError(f"hierarchy.{key} must be numeric")
    if config["hierarchy"]["gamma_min"] <= 0:
        raise InputError("hierarchy.gamma_min must be positive")
    if config["hierarchy"]["gamma_max"] < config["hierarchy"]["gamma_min"]:
        raise InputError("hierarchy.gamma_max must be at least gamma_min")
    if config["hierarchy"]["gamma_growth"] <= 1:
        raise InputError("hierarchy.gamma_growth must be greater than 1")
    if not 0 < config["hierarchy"]["max_child_fraction"] < 1:
        raise InputError("hierarchy.max_child_fraction must be between 0 and 1")
    if (
        not isinstance(config["hierarchy"]["local_grid_points"], int)
        or config["hierarchy"]["local_grid_points"] < 2
    ):
        raise InputError("hierarchy.local_grid_points must be at least 2")
    if (
        not isinstance(config["hierarchy"]["tiny_fragment_size"], int)
        or config["hierarchy"]["tiny_fragment_size"] < 1
    ):
        raise InputError("hierarchy.tiny_fragment_size must be positive")
    if (
        not isinstance(config["hierarchy"]["publication_seeds"], int)
        or config["hierarchy"]["publication_seeds"] < 5
    ):
        raise InputError("hierarchy.publication_seeds must be at least 5")
    for key in ("max_tiny_fragment_fraction", "stability_threshold"):
        value = config["hierarchy"][key]
        if not isinstance(value, (int, float)) or isinstance(value, bool) or not 0 <= value <= 1:
            raise InputError(f"hierarchy.{key} must be between 0 and 1")
    minimum_quality = config["hierarchy"]["min_split_quality"]
    if minimum_quality is not None and (
        not isinstance(minimum_quality, (int, float)) or isinstance(minimum_quality, bool)
    ):
        raise InputError("hierarchy.min_split_quality must be numeric or null")
    evalue = config["search"]["evalue"]
    if not isinstance(evalue, (int, float)) or isinstance(evalue, bool) or evalue <= 0:
        raise InputError("search.evalue must be positive")
    for key in (
        "min_query_coverage",
        "min_target_coverage",
        "min_bidirectional_coverage",
    ):
        value = config["edges"][key]
        if not isinstance(value, (int, float)) or not 0 <= value <= 100:
            raise InputError(f"edges.{key} must be between 0 and 100")
    tolerance = config["edges"]["best_hit_tolerance"]
    if not isinstance(tolerance, (int, float)) or isinstance(tolerance, bool) or tolerance < 0:
        raise InputError("edges.best_hit_tolerance must be non-negative")
    if not isinstance(config["edges"]["apply_coverage_filter"], bool):
        raise InputError("edges.apply_coverage_filter must be boolean")
    overlap = config["evolution"]["network_overlap_threshold"]
    if not isinstance(overlap, (int, float)) or not 0 <= overlap <= 1:
        raise InputError("evolution.network_overlap_threshold must be between 0 and 1")
    log_level = config["runtime"]["log_level"]
    if not isinstance(log_level, str) or log_level.upper() not in _ALLOWED_LOG_LEVELS:
        expected = ", ".join(sorted(_ALLOWED_LOG_LEVELS))
        raise InputError(f"runtime.log_level must be one of: {expected}")
    config["runtime"]["log_level"] = log_level.upper()


def dump_config(config: dict[str, Any]) -> str:
    return str(yaml.safe_dump(config, sort_keys=True, allow_unicode=True))


def hierarchy_config(options: dict[str, Any]) -> HierarchyConfig:
    """One validated config seam for production, prototype and publication gates."""
    from ogprofiler.exceptions import HierarchyError
    from ogprofiler.hierarchy.engine import HierarchyConfig
    from ogprofiler.hierarchy.resolution import ResolutionSearchConfig

    policy = options["admission_policy"]
    legacy = policy == "legacy_strict"
    if not legacy and options["min_family_size"] is not None:
        raise InputError(
            "min_family_size is legacy-only; explicitly migrate to recursion_stop_size"
        )
    if legacy and options["resolution_strategy"] == "bounded_adaptive_v2":
        raise InputError("legacy_strict requires adaptive or log_grid")
    if not legacy and options["resolution_strategy"] != "bounded_adaptive_v2":
        raise InputError("Historical searches require explicit admission_policy=legacy_strict")
    config = HierarchyConfig(
        method=options["method"],
        seed=options["seed"],
        max_depth=options["max_depth"],
        stability_mode=options["stability_mode"],
        subtree_workers=options["subtree_workers"],
        subtree_release_size=options["subtree_release_size"],
        recursion_stop_size=options["recursion_stop_size"],
        component_leiden_call_budget=options["component_leiden_call_budget"],
        leiden_iterations=options["leiden_iterations"],
        resolution=ResolutionSearchConfig(
            strategy=options["resolution_strategy"],
            admission_policy=policy,
            topology_policy=options.get("topology_policy", "kway_v1"),
            gamma_min=options["gamma_min"],
            gamma_max=options["gamma_max"],
            growth_factor=options["gamma_growth"],
            local_grid_points=options["local_grid_points"],
            min_child_size=options["min_family_size"]
            if options["min_family_size"] is not None
            else 2,
            max_child_fraction=options["max_child_fraction"],
            tiny_fragment_size=options["tiny_fragment_size"],
            max_tiny_fragment_fraction=options["max_tiny_fragment_fraction"],
            stability_threshold=options["stability_threshold"],
            publication_seeds=options["publication_seeds"],
            min_quality=options["min_split_quality"],
            max_candidate_evaluations=options["max_candidate_evaluations"],
            max_coarse_candidates=options["max_coarse_candidates"],
            rescue_grid_points=options["rescue_grid_points"],
        ),
    )
    try:
        config.resolution.validate()
    except HierarchyError as error:
        raise InputError(str(error)) from error
    for key in ("recursion_stop_size", "leiden_iterations"):
        if type(options[key]) is not int or options[key] < 1:
            raise InputError(f"hierarchy.{key} must be a positive integer")
    cap = options["component_leiden_call_budget"]
    if cap is not None and (type(cap) is not int or cap < 1):
        raise InputError("component_leiden_call_budget must be a positive integer or null")
    if cap is not None and legacy:
        raise InputError("component budget requires bounded_adaptive_v2")
    return config
