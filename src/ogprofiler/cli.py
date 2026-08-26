"""Command-line interface for explicit OGProfiler 2 pipeline stages."""

from __future__ import annotations

import argparse
import sys
from collections.abc import Sequence
from dataclasses import asdict
from pathlib import Path

from ogprofiler import __version__
from ogprofiler.config import dump_config, load_config
from ogprofiler.core.manifest import sha256_file, sha256_json, write_json
from ogprofiler.core.models import RunManifest
from ogprofiler.core.workspace import Workspace
from ogprofiler.exceptions import OGProfilerError
from ogprofiler.input.proteomes import prepare_proteomes
from ogprofiler.logging import configure_logging

PREPARE_ALGORITHM_VERSION = "prepare-v1"
HIERARCHY_PROTOTYPE_VERSION = "hierarchy-prototype-v1"


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="ogprofiler")
    parser.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    subparsers = parser.add_subparsers(dest="command", required=True)

    run = subparsers.add_parser("run", help="execute or resume the standard end-to-end pipeline")
    run.add_argument("--proteomes", type=Path)
    run.add_argument("--out", required=True, type=Path)
    run.add_argument("--config", type=Path)
    run.add_argument("--set", dest="overrides", action="append", default=[])
    run.add_argument(
        "--from-stage",
        choices=(
            "prepare",
            "search",
            "edges",
            "components",
            "hierarchy",
            "annotate-network",
            "export",
        ),
        default="prepare",
    )
    run.add_argument(
        "--until-stage",
        choices=(
            "prepare",
            "search",
            "edges",
            "components",
            "hierarchy",
            "annotate-network",
            "export",
        ),
        default="export",
    )
    run.add_argument("--dry-run", action="store_true")

    status = subparsers.add_parser("status", help="summarize pipeline and checkpoint state")
    status.add_argument("--run", required=True, type=Path)
    status.add_argument("--json", action="store_true")

    inspect = subparsers.add_parser("inspect", help="inspect dataset and hierarchy summaries")
    inspect.add_argument("--run", required=True, type=Path)
    inspect.add_argument("--component", type=int)
    inspect.add_argument("--json", action="store_true")

    prepare = subparsers.add_parser(
        "prepare",
        help="validate proteomes and write deterministic input metadata",
    )
    prepare.add_argument("--proteomes", required=True, type=Path)
    prepare.add_argument("--out", required=True, type=Path)
    prepare.add_argument("--config", type=Path)
    prepare.add_argument(
        "--set",
        dest="overrides",
        action="append",
        default=[],
        metavar="KEY=VALUE",
        help="override a config value, for example hierarchy.seed=7",
    )

    prototype = subparsers.add_parser(
        "prototype-hierarchy",
        help="import a legacy GML SSN and run component-centric hierarchical Leiden",
    )
    prototype.add_argument("--ssn", required=True, type=Path)
    prototype.add_argument("--out", required=True, type=Path)
    prototype.add_argument("--weight-attribute", default="NBS")
    prototype.add_argument("--vertex-id-attribute")
    prototype.add_argument("--config", type=Path)
    prototype.add_argument(
        "--set",
        dest="overrides",
        action="append",
        default=[],
        metavar="KEY=VALUE",
    )

    search = subparsers.add_parser(
        "search",
        help="run a directional all-vs-all homolog search on prepared proteins",
    )
    search.add_argument("--run", required=True, type=Path)
    search.add_argument("--backend", choices=("diamond", "mmseqs", "blastp"))
    search.add_argument("--config", type=Path)
    search.add_argument(
        "--set",
        dest="overrides",
        action="append",
        default=[],
        metavar="KEY=VALUE",
    )

    edges = subparsers.add_parser(
        "edges",
        help="normalize directional hits and build canonical retained edges",
    )
    edges.add_argument("--run", required=True, type=Path)
    edges.add_argument("--method", choices=("lrb", "rbh", "ar", "arb"))
    edges.add_argument("--config", type=Path)
    edges.add_argument(
        "--set",
        dest="overrides",
        action="append",
        default=[],
        metavar="KEY=VALUE",
    )

    components = subparsers.add_parser(
        "components",
        help="build a Union-Find component index and partition retained edges",
    )
    components.add_argument("--run", required=True, type=Path)
    components.add_argument("--config", type=Path)
    components.add_argument(
        "--set",
        dest="overrides",
        action="append",
        default=[],
        metavar="KEY=VALUE",
    )
    hierarchy = subparsers.add_parser(
        "hierarchy",
        help="run production hierarchical Leiden for one partitioned component",
    )
    hierarchy.add_argument("--run", required=True, type=Path)
    hierarchy.add_argument("--component-id", required=True, type=int)
    hierarchy.add_argument("--config", type=Path)
    hierarchy.add_argument("--set", dest="overrides", action="append", default=[])
    hierarchy_all = subparsers.add_parser(
        "hierarchy-all",
        help="schedule every non-singleton component with checkpointed resume",
    )
    hierarchy_all.add_argument("--run", required=True, type=Path)
    hierarchy_all.add_argument("--config", type=Path)
    hierarchy_all.add_argument("--failed-only", action="store_true")
    hierarchy_all.add_argument("--set", dest="overrides", action="append", default=[])
    annotate = subparsers.add_parser(
        "annotate-network",
        help="annotate hierarchy nodes using child species-overlap heuristics",
    )
    annotate.add_argument("--run", required=True, type=Path)
    annotate.add_argument("--config", type=Path)
    annotate.add_argument("--set", dest="overrides", action="append", default=[])
    export = subparsers.add_parser(
        "export",
        help="write stable terminal-family and hierarchy exchange tables",
    )
    export.add_argument("kind", nargs="?", choices=("results", "graph"), default="results")
    export.add_argument("--run", required=True, type=Path)
    export.add_argument("--component", type=int, help="component ID for graph export")
    export.add_argument("--format", choices=("graphml",), default="graphml")
    export.add_argument(
        "--family-fasta",
        dest="fasta_families",
        action="append",
        default=[],
        metavar="FAMILY_ID",
        help="export FASTA for one selected family; may be repeated",
    )
    export.add_argument(
        "--all-family-fasta",
        action="store_true",
        help="explicitly export one FASTA file for every terminal family",
    )
    orthologs = subparsers.add_parser(
        "orthologs",
        help="stream hierarchy-derived pairwise ortholog candidates",
    )
    orthologs.add_argument("--run", required=True, type=Path)
    orthologs.add_argument("--config", type=Path)
    orthologs.add_argument("--set", dest="overrides", action="append", default=[])
    orthologs.add_argument(
        "--emit-pairwise-orthologs",
        action=argparse.BooleanOptionalAction,
        default=None,
        help="explicitly enable or disable the potentially quadratic pair output",
    )
    phylogeny = subparsers.add_parser(
        "annotate",
        help="add optional phylogenetic refinement evidence for selected families",
    )
    phylogeny.add_argument("--run", required=True, type=Path)
    phylogeny.add_argument("--phylogenetic-refinement", action="store_true")
    phylogeny.add_argument("--family", dest="families", action="append", default=[])
    phylogeny.add_argument("--species-tree", type=Path)
    phylogeny.add_argument(
        "--rooting", choices=("midpoint", "species-tree-aware", "outgroup")
    )
    phylogeny.add_argument("--outgroup")
    phylogeny.add_argument("--config", type=Path)
    phylogeny.add_argument("--set", dest="overrides", action="append", default=[])
    benchmark = subparsers.add_parser(
        "benchmark",
        help="plan or evaluate scientific parameter and method benchmarks",
    )
    benchmark.add_argument(
        "kind",
        choices=(
            "plan",
            "metrics",
            "compare",
            "synthetic-plan",
            "synthetic-generate",
            "synthetic-metrics",
            "synthetic-map",
        ),
    )
    benchmark.add_argument("--out", required=True, type=Path)
    benchmark.add_argument("--run", type=Path)
    benchmark.add_argument("--ground-truth", type=Path)
    benchmark.add_argument("--method", default="ogprofiler2")
    benchmark.add_argument("--dataset")
    benchmark.add_argument("--metrics", type=Path, action="append", default=[])
    benchmark.add_argument("--legacy-normalized", type=Path)
    benchmark.add_argument("--large-family-size", type=int, default=20)
    benchmark.add_argument("--dataset-root", type=Path)
    benchmark.add_argument("--replicates", type=int, default=1)
    benchmark.add_argument("--index", type=int, default=0)
    benchmark.add_argument("--ancestral-families", type=int, default=12)
    benchmark.add_argument("--sequence-length", type=int, default=120)
    return parser


def _prepare(args: argparse.Namespace, command: list[str]) -> int:
    config = load_config(str(args.config) if args.config else None, args.overrides)
    workspace = Workspace.create(args.out)
    run_yaml = dump_config(config)
    (workspace.root / "run.yaml").write_text(run_yaml, encoding="utf-8", newline="\n")

    logger = configure_logging(
        workspace.root / "ogprofiler.log",
        level=config["runtime"]["log_level"],
    ).bind(stage="prepare", task="dataset")
    logger.info("Preparing proteomes from %s", args.proteomes)
    result = prepare_proteomes(args.proteomes, workspace.root / "input", config)
    manifest = RunManifest(
        ogprofiler_version=__version__,
        algorithm_version=PREPARE_ALGORITHM_VERSION,
        command=command,
        random_seed=config["hierarchy"]["seed"],
        resolved_config_sha256=sha256_json(config),
        dataset=result.dataset_manifest,
    )
    write_json(workspace.root / "manifest.json", manifest.to_dict())
    logger.info(
        "Prepared %d proteins from %d species; dataset_sha256=%s",
        result.dataset_manifest.n_proteins,
        result.dataset_manifest.n_species,
        result.dataset_manifest.dataset_sha256,
    )
    return 0


def _run_pipeline(args: argparse.Namespace, command: list[str]) -> int:
    from ogprofiler.ux import pipeline_commands, write_run_provenance

    commands = pipeline_commands(
        args.out,
        args.proteomes,
        config=args.config,
        overrides=tuple(args.overrides),
        from_stage=str(args.from_stage),
        until_stage=str(args.until_stage),
    )
    if args.dry_run:
        for stage_command in commands:
            print("ogprofiler " + " ".join(stage_command))
        return 0
    for stage_command in commands:
        status = main(stage_command)
        if status != 0:
            return status
    config_path = args.config
    if config_path is None and (args.out / "run.yaml").is_file():
        config_path = args.out / "run.yaml"
    resolved = load_config(str(config_path) if config_path else None, args.overrides)
    provenance = write_run_provenance(args.out, resolved, command)
    print(f"pipeline complete: run={args.out} provenance={provenance}")
    return 0


def _print_report(report: dict[str, object], *, as_json: bool) -> None:
    import json

    if as_json:
        print(json.dumps(report, indent=2, sort_keys=True))
        return
    print(f"run: {report['run_root']}")
    print(f"overall: {report['overall']}")
    stages = report.get("stages", [])
    if isinstance(stages, list):
        for row in stages:
            if isinstance(row, dict):
                print(f"{row['stage']:<18} {row['status']:<11} {row['artifact']}")
    counts = report.get("counts")
    if isinstance(counts, dict) and counts:
        print("counts: " + ", ".join(f"{key}={value}" for key, value in counts.items()))
    component = report.get("component")
    if isinstance(component, dict):
        print("component: " + ", ".join(f"{key}={value}" for key, value in component.items()))


def _status(args: argparse.Namespace) -> int:
    from ogprofiler.ux import run_status

    _print_report(run_status(args.run), as_json=bool(args.json))
    return 0


def _inspect(args: argparse.Namespace) -> int:
    from ogprofiler.ux import inspect_run

    _print_report(inspect_run(args.run, args.component), as_json=bool(args.json))
    return 0


def _prototype_hierarchy(args: argparse.Namespace, command: list[str]) -> int:
    from ogprofiler.graph.components import extract_components
    from ogprofiler.graph.legacy import import_legacy_ssn
    from ogprofiler.hierarchy.engine import HierarchyConfig, infer_component_hierarchy
    from ogprofiler.hierarchy.resolution import ResolutionSearchConfig
    from ogprofiler.hierarchy.validation import validate_hierarchy
    from ogprofiler.storage.hierarchy import write_hierarchy_result

    config = load_config(str(args.config) if args.config else None, args.overrides)
    workspace = Workspace.create(args.out)
    (workspace.root / "run.yaml").write_text(dump_config(config), encoding="utf-8", newline="\n")
    logger = configure_logging(
        workspace.root / "ogprofiler.log", level=config["runtime"]["log_level"]
    ).bind(stage="prototype-hierarchy", task="legacy-ssn")
    logger.info("Importing legacy SSN from %s", args.ssn)
    edge_table = import_legacy_ssn(
        args.ssn,
        weight_attribute=args.weight_attribute,
        vertex_id_attribute=args.vertex_id_attribute,
    )
    edge_table.write_parquet(workspace.root / "edges" / "legacy_ssn_edges.parquet")
    components = extract_components(edge_table)
    hierarchy_options = config["hierarchy"]
    hierarchy_config = HierarchyConfig(
        method=hierarchy_options["method"],
        seed=hierarchy_options["seed"],
        max_depth=hierarchy_options["max_depth"],
        stability_mode=hierarchy_options["stability_mode"],
        resolution=ResolutionSearchConfig(
            strategy=hierarchy_options["resolution_strategy"],
            gamma_min=hierarchy_options["gamma_min"],
            gamma_max=hierarchy_options["gamma_max"],
            growth_factor=hierarchy_options["gamma_growth"],
            local_grid_points=hierarchy_options["local_grid_points"],
            min_child_size=hierarchy_options["min_family_size"],
            max_child_fraction=hierarchy_options["max_child_fraction"],
            tiny_fragment_size=hierarchy_options["tiny_fragment_size"],
            max_tiny_fragment_fraction=hierarchy_options["max_tiny_fragment_fraction"],
            stability_threshold=hierarchy_options["stability_threshold"],
            publication_seeds=hierarchy_options["publication_seeds"],
            min_quality=hierarchy_options["min_split_quality"],
        ),
    )
    summaries: list[dict[str, object]] = []
    for component in components:
        component.write(workspace.root / "components")
        component_logger = logger.bind(component=component.component_id)
        component_logger.info(
            "Inferring component with %d vertices and %d edges",
            component.n_vertices,
            component.n_edges,
        )
        result = infer_component_hierarchy(component, hierarchy_config)
        validate_hierarchy(component, result)
        result_dir = workspace.root / "hierarchy" / "components" / f"{component.component_id:08d}"
        write_hierarchy_result(result_dir, result)
        summaries.append(
            {
                "component_id": component.component_id,
                "n_vertices": component.n_vertices,
                "n_edges": component.n_edges,
                **asdict(result.metrics),
            }
        )
    write_json(
        workspace.root / "manifest.json",
        {
            "ogprofiler_version": __version__,
            "algorithm_version": HIERARCHY_PROTOTYPE_VERSION,
            "command": command,
            "ssn_sha256": sha256_file(args.ssn),
            "resolved_config_sha256": sha256_json(config),
            "weight_attribute": args.weight_attribute,
            "component_count": len(components),
            "components": summaries,
        },
    )
    logger.info("Completed %d connected components", len(components))
    return 0


def _search(args: argparse.Namespace, command: list[str]) -> int:
    from ogprofiler.search.base import SearchParameters
    from ogprofiler.search.stage import create_backend, run_search_stage

    default_config = args.config
    if default_config is None and (args.run / "run.yaml").is_file():
        default_config = args.run / "run.yaml"
    overrides = list(args.overrides)
    if args.backend is not None:
        overrides.append(f"search.backend={args.backend}")
    config = load_config(str(default_config) if default_config else None, overrides)
    search_options = config["search"]
    backend = create_backend(config)
    parameters = SearchParameters(
        threads=int(search_options["threads"]),
        evalue=float(search_options["evalue"]),
        sensitivity=str(search_options["sensitivity"]),
        max_target_seqs=int(search_options["max_target_seqs"]),
    )
    logger = configure_logging(
        args.run / "ogprofiler.log", level=config["runtime"]["log_level"]
    ).bind(stage="search", task=backend.name)
    logger.info("Starting %s all-vs-all search", backend.name)
    result = run_search_stage(args.run, backend, parameters, command)
    logger.info(
        "%s search; hits=%s manifest=%s",
        "Reused verified" if result.reused else "Completed",
        result.hits_path,
        result.manifest_path,
    )
    return 0


def _edges(args: argparse.Namespace, command: list[str]) -> int:
    from ogprofiler.similarity.engine import EdgeBuildConfig
    from ogprofiler.similarity.stage import run_edge_stage

    default_config = args.config
    if default_config is None and (args.run / "run.yaml").is_file():
        default_config = args.run / "run.yaml"
    overrides = list(args.overrides)
    if args.method is not None:
        overrides.append(f"edges.method={args.method}")
    config = load_config(str(default_config) if default_config else None, overrides)
    options = config["edges"]
    edge_config = EdgeBuildConfig(
        method=str(options["method"]),
        normalization=str(config["similarity"]["normalization"]),
        min_query_coverage=float(options["min_query_coverage"]),
        min_target_coverage=float(options["min_target_coverage"]),
        min_bidirectional_coverage=float(options["min_bidirectional_coverage"]),
        best_hit_tolerance=float(options["best_hit_tolerance"]),
        symmetrization=str(options["symmetrization"]),
    )
    logger = configure_logging(
        args.run / "ogprofiler.log", level=config["runtime"]["log_level"]
    ).bind(stage="edges", task=edge_config.method)
    logger.info("Building canonical retained edges")
    path, reused = run_edge_stage(args.run, edge_config, command)
    logger.info("%s edge artifact: %s", "Reused verified" if reused else "Completed", path)
    return 0


def _components(args: argparse.Namespace, command: list[str]) -> int:
    from ogprofiler.graph.stage import ComponentStageConfig, run_component_stage

    default_config = args.config
    if default_config is None and (args.run / "run.yaml").is_file():
        default_config = args.run / "run.yaml"
    config = load_config(str(default_config) if default_config else None, args.overrides)
    options = config["components"]
    component_config = ComponentStageConfig(
        edge_batch_size=int(options["edge_batch_size"]),
        max_open_files=int(options["max_open_files"]),
    )
    logger = configure_logging(
        args.run / "ogprofiler.log", level=config["runtime"]["log_level"]
    ).bind(stage="components", task="union-find")
    logger.info("Building component index and partitioned edges")
    path, reused = run_component_stage(args.run, component_config, command)
    logger.info("%s component artifacts: %s", "Reused verified" if reused else "Completed", path)
    return 0


def _hierarchy(args: argparse.Namespace, command: list[str]) -> int:
    from ogprofiler.hierarchy.engine import HierarchyConfig
    from ogprofiler.hierarchy.resolution import ResolutionSearchConfig
    from ogprofiler.hierarchy.stage import run_hierarchy_component_stage

    default_config = args.config
    if default_config is None and (args.run / "run.yaml").is_file():
        default_config = args.run / "run.yaml"
    config = load_config(str(default_config) if default_config else None, args.overrides)
    options = config["hierarchy"]
    hierarchy_config = HierarchyConfig(
        method=options["method"],
        seed=options["seed"],
        max_depth=options["max_depth"],
        stability_mode=options["stability_mode"],
        subtree_workers=options["subtree_workers"],
        subtree_release_size=options["subtree_release_size"],
        resolution=ResolutionSearchConfig(
            strategy=options["resolution_strategy"],
            gamma_min=options["gamma_min"],
            gamma_max=options["gamma_max"],
            growth_factor=options["gamma_growth"],
            local_grid_points=options["local_grid_points"],
            min_child_size=options["min_family_size"],
            max_child_fraction=options["max_child_fraction"],
            tiny_fragment_size=options["tiny_fragment_size"],
            max_tiny_fragment_fraction=options["max_tiny_fragment_fraction"],
            stability_threshold=options["stability_threshold"],
            publication_seeds=options["publication_seeds"],
            min_quality=options["min_split_quality"],
        ),
    )
    logger = configure_logging(
        args.run / "ogprofiler.log", level=config["runtime"]["log_level"]
    ).bind(stage="hierarchy", task=str(args.component_id), component=args.component_id)
    logger.info("Inferring production component hierarchy")
    path, reused = run_hierarchy_component_stage(
        args.run, args.component_id, hierarchy_config, command
    )
    logger.info("%s hierarchy artifacts: %s", "Reused verified" if reused else "Completed", path)
    return 0


def _hierarchy_all(args: argparse.Namespace, command: list[str]) -> int:
    from ogprofiler.exceptions import HierarchyError
    from ogprofiler.hierarchy.engine import HierarchyConfig
    from ogprofiler.hierarchy.resolution import ResolutionSearchConfig
    from ogprofiler.hierarchy.scheduler import run_hierarchy_scheduler

    default_config = args.config
    if default_config is None and (args.run / "run.yaml").is_file():
        default_config = args.run / "run.yaml"
    config = load_config(str(default_config) if default_config else None, args.overrides)
    options = config["hierarchy"]
    hierarchy_config = HierarchyConfig(
        method=options["method"],
        seed=options["seed"],
        max_depth=options["max_depth"],
        stability_mode=options["stability_mode"],
        subtree_workers=options["subtree_workers"],
        subtree_release_size=options["subtree_release_size"],
        resolution=ResolutionSearchConfig(
            strategy=options["resolution_strategy"],
            gamma_min=options["gamma_min"],
            gamma_max=options["gamma_max"],
            growth_factor=options["gamma_growth"],
            local_grid_points=options["local_grid_points"],
            min_child_size=options["min_family_size"],
            max_child_fraction=options["max_child_fraction"],
            tiny_fragment_size=options["tiny_fragment_size"],
            max_tiny_fragment_fraction=options["max_tiny_fragment_fraction"],
            stability_threshold=options["stability_threshold"],
            publication_seeds=options["publication_seeds"],
            min_quality=options["min_split_quality"],
        ),
    )
    logger = configure_logging(
        args.run / "ogprofiler.log", level=config["runtime"]["log_level"]
    ).bind(stage="hierarchy", task="scheduler")
    logger.info("Scheduling component-local hierarchy workers")
    result = run_hierarchy_scheduler(
        args.run,
        hierarchy_config,
        command,
        workers=int(config["runtime"]["workers"]),
        retries=int(config["runtime"]["component_retries"]),
        failed_only=bool(args.failed_only),
    )
    logger.info(
        "Hierarchy schedule complete: completed=%d skipped=%d failed=%d singletons=%d",
        result.completed,
        result.skipped,
        result.failed,
        result.singleton_components,
    )
    if result.failed:
        raise HierarchyError(f"Hierarchy scheduler reported {result.failed} failed components")
    return 0


def _annotate_network(args: argparse.Namespace, command: list[str]) -> int:
    from ogprofiler.evolution.stage import run_network_annotation_stage

    default_config = args.config
    if default_config is None and (args.run / "run.yaml").is_file():
        default_config = args.run / "run.yaml"
    config = load_config(str(default_config) if default_config else None, args.overrides)
    threshold = float(config["evolution"]["network_overlap_threshold"])
    logger = configure_logging(
        args.run / "ogprofiler.log", level=config["runtime"]["log_level"]
    ).bind(stage="evolution", task="network-overlap")
    logger.info("Annotating hierarchy nodes with network events")
    path, components, reused = run_network_annotation_stage(args.run, threshold, command)
    logger.info(
        "Completed network annotations: components=%d reused=%d manifest=%s",
        components,
        reused,
        path,
    )
    return 0


def _export(args: argparse.Namespace, command: list[str]) -> int:
    if args.kind == "graph":
        from ogprofiler.output.graph import export_component_graphml

        if args.component is None:
            raise OGProfilerError("export graph requires --component")
        path = export_component_graphml(args.run, int(args.component), command)
        logger = configure_logging(args.run / "ogprofiler.log").bind(
            stage="output", task="graph", component=args.component
        )
        logger.info("Completed component GraphML export: %s", path)
        return 0

    from ogprofiler.output.stage import run_export_stage

    logger = configure_logging(args.run / "ogprofiler.log").bind(
        stage="output", task="terminal-families"
    )
    logger.info("Exporting stable terminal-family result tables")
    path, reused, families = run_export_stage(
        args.run,
        command,
        fasta_families=tuple(args.fasta_families),
        all_family_fasta=bool(args.all_family_fasta),
    )
    logger.info(
        "%s final result export: families=%d manifest=%s",
        "Reused verified" if reused else "Completed",
        families,
        path,
    )
    return 0


def _orthologs(args: argparse.Namespace, command: list[str]) -> int:
    from ogprofiler.orthology.stage import run_orthology_stage

    default_config = args.config
    if default_config is None and (args.run / "run.yaml").is_file():
        default_config = args.run / "run.yaml"
    overrides = list(args.overrides)
    if args.emit_pairwise_orthologs is not None:
        value = "true" if args.emit_pairwise_orthologs else "false"
        overrides.append(f"output.emit_pairwise_orthologs={value}")
    config = load_config(str(default_config) if default_config else None, overrides)
    logger = configure_logging(
        args.run / "ogprofiler.log", level=config["runtime"]["log_level"]
    ).bind(stage="orthology", task="pairwise-candidates")
    if not bool(config["output"]["emit_pairwise_orthologs"]):
        logger.info("Pairwise ortholog output is disabled by configuration")
        return 0

    def progress(pairs: int, chunks: int) -> None:
        logger.info("Streaming ortholog candidates: pairs=%d chunks=%d", pairs, chunks)

    logger.info("Generating cross-child ortholog candidates")
    path, reused, pairs = run_orthology_stage(
        args.run,
        command,
        chunk_size=int(config["output"]["ortholog_pair_chunk_size"]),
        progress=progress,
    )
    logger.info(
        "%s pairwise ortholog output: pairs=%d manifest=%s",
        "Reused verified" if reused else "Completed",
        pairs,
        path,
    )
    return 0


def _phylogenetic_refinement(args: argparse.Namespace, command: list[str]) -> int:
    from ogprofiler.phylogeny.backends import FastTreeBackend, MafftBackend
    from ogprofiler.phylogeny.reconciliation import LcaReconciliationBackend
    from ogprofiler.phylogeny.stage import run_phylogenetic_refinement_stage

    if not args.phylogenetic_refinement:
        raise OGProfilerError("annotate requires --phylogenetic-refinement")
    default_config = args.config
    if default_config is None and (args.run / "run.yaml").is_file():
        default_config = args.run / "run.yaml"
    overrides = list(args.overrides)
    if args.rooting is not None:
        overrides.append(f"phylogeny.rooting={args.rooting}")
    config = load_config(str(default_config) if default_config else None, overrides)
    options = config["phylogeny"]
    logger = configure_logging(
        args.run / "ogprofiler.log", level=config["runtime"]["log_level"]
    ).bind(stage="phylogeny", task="selected-families")
    logger.info("Starting optional phylogenetic refinement")
    path, selected, reused = run_phylogenetic_refinement_stage(
        args.run,
        command,
        explicit_family_ids=tuple(args.families),
        selection_events=set(options["selection_events"]),
        large_family_size=int(options["large_family_size"]),
        max_families=int(options["max_families"]),
        rooting=str(options["rooting"]),
        outgroup=args.outgroup,
        species_tree_path=args.species_tree,
        threads=int(options["threads"]),
        alignment_backend=MafftBackend(str(options["alignment_executable"])),
        tree_backend=FastTreeBackend(str(options["tree_executable"])),
        reconciliation_backend=LcaReconciliationBackend(),
    )
    logger.info(
        "Completed phylogenetic refinement: selected=%d reused=%d manifest=%s",
        selected,
        reused,
        path,
    )
    return 0


def _benchmark(args: argparse.Namespace) -> int:
    if args.kind in {"synthetic-plan", "synthetic-generate"}:
        from ogprofiler.benchmark.synthetic import (
            generate_scenario_matrix,
            generate_synthetic_dataset,
            write_scenario_matrix,
        )

        scenarios = generate_scenario_matrix(replicates=int(args.replicates))
        if args.kind == "synthetic-plan":
            json_path, _ = write_scenario_matrix(args.out, scenarios)
            print(f"synthetic scenario matrix: runs={len(scenarios)} manifest={json_path}")
            return 0
        if not 0 <= args.index < len(scenarios):
            raise OGProfilerError(
                f"synthetic-generate --index must be between 0 and {len(scenarios) - 1}"
            )
        manifest = generate_synthetic_dataset(
            args.out,
            scenarios[int(args.index)],
            ancestral_families=int(args.ancestral_families),
            sequence_length=int(args.sequence_length),
        )
        print(f"synthetic dataset: {manifest}")
        return 0
    if args.kind == "synthetic-metrics":
        from ogprofiler.benchmark.synthetic_metrics import (
            evaluate_synthetic_run,
            write_synthetic_evaluation,
        )

        if args.run is None or args.dataset_root is None:
            raise OGProfilerError(
                "benchmark synthetic-metrics requires --run and --dataset-root"
            )
        result = evaluate_synthetic_run(args.run, args.dataset_root, method=str(args.method))
        output = args.out / "synthetic-metrics.json" if args.out.suffix == "" else args.out
        write_synthetic_evaluation(output, result)
        print(f"synthetic recovery metrics: {output}")
        return 0
    if args.kind == "synthetic-map":
        from ogprofiler.benchmark.synthetic_metrics import aggregate_applicability

        if not args.metrics:
            raise OGProfilerError("benchmark synthetic-map requires --metrics files")
        json_path, _ = aggregate_applicability(args.metrics, args.out)
        print(f"Leiden applicability map: {json_path}")
        return 0
    if args.kind == "plan":
        from ogprofiler.benchmark.matrix import generate_ofat_matrix, write_matrix

        runs = generate_ofat_matrix()
        json_path, _ = write_matrix(args.out, runs)
        print(f"parameter matrix: runs={len(runs)} manifest={json_path}")
        return 0
    if args.kind == "metrics":
        from ogprofiler.benchmark.metrics import (
            evaluate_legacy_membership,
            evaluate_run,
            write_evaluation,
        )

        if args.ground_truth is None or args.dataset is None:
            raise OGProfilerError(
                "benchmark metrics requires --ground-truth and --dataset"
            )
        if args.legacy_normalized is not None:
            result = evaluate_legacy_membership(
                args.legacy_normalized,
                args.ground_truth,
                method=str(args.method),
                dataset=str(args.dataset),
                large_family_size=int(args.large_family_size),
            )
        else:
            if args.run is None:
                raise OGProfilerError(
                    "benchmark metrics requires --run or --legacy-normalized"
                )
            result = evaluate_run(
                args.run,
                args.ground_truth,
                method=str(args.method),
                dataset=str(args.dataset),
                large_family_size=int(args.large_family_size),
            )
        output = args.out / "scientific-metrics.json" if args.out.suffix == "" else args.out
        write_evaluation(output, result)
        print(f"scientific metrics: {output}")
        return 0
    from ogprofiler.benchmark.comparison import compare_methods

    if not args.metrics:
        raise OGProfilerError("benchmark compare requires at least one --metrics file")
    json_path, _ = compare_methods(args.metrics, args.out)
    print(f"method comparison: {json_path}")
    return 0


def main(argv: Sequence[str] | None = None) -> int:
    parser = build_parser()
    supplied = list(argv) if argv is not None else sys.argv[1:]
    try:
        args = parser.parse_args(supplied)
        if args.command == "run":
            return _run_pipeline(args, ["ogprofiler", *supplied])
        if args.command == "status":
            return _status(args)
        if args.command == "inspect":
            return _inspect(args)
        if args.command == "prepare":
            return _prepare(args, ["ogprofiler", *supplied])
        if args.command == "prototype-hierarchy":
            return _prototype_hierarchy(args, ["ogprofiler", *supplied])
        if args.command == "search":
            return _search(args, ["ogprofiler", *supplied])
        if args.command == "edges":
            return _edges(args, ["ogprofiler", *supplied])
        if args.command == "components":
            return _components(args, ["ogprofiler", *supplied])
        if args.command == "hierarchy":
            return _hierarchy(args, ["ogprofiler", *supplied])
        if args.command == "hierarchy-all":
            return _hierarchy_all(args, ["ogprofiler", *supplied])
        if args.command == "annotate-network":
            return _annotate_network(args, ["ogprofiler", *supplied])
        if args.command == "export":
            return _export(args, ["ogprofiler", *supplied])
        if args.command == "orthologs":
            return _orthologs(args, ["ogprofiler", *supplied])
        if args.command == "annotate":
            return _phylogenetic_refinement(args, ["ogprofiler", *supplied])
        if args.command == "benchmark":
            return _benchmark(args)
        parser.error(f"Unknown command: {args.command}")
    except OGProfilerError as error:
        print(f"error: {error}", file=sys.stderr)
        return 2
    return 2


if __name__ == "__main__":
    raise SystemExit(main())
