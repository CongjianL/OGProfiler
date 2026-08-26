"""Opt-in deterministic component GraphML export."""

from __future__ import annotations

import os
import uuid
from pathlib import Path
from xml.sax.saxutils import escape, quoteattr

import pyarrow.parquet as pq

from ogprofiler.core.manifest import sha256_file, write_json
from ogprofiler.exceptions import ExportError
from ogprofiler.graph.partition import read_component_edges


def export_component_graphml(run_root: Path, component_id: int, command: list[str]) -> Path:
    index_path = run_root / "components" / "index.parquet"
    proteins_path = run_root / "input" / "proteins.parquet"
    if not index_path.is_file() or not proteins_path.is_file():
        raise ExportError("Component index and protein metadata are required for graph export")
    index_rows = pq.read_table(
        index_path, filters=[("component_id", "=", component_id)]
    ).to_pylist()
    protein_ids = sorted(int(row["protein_id"]) for row in index_rows)
    if not protein_ids:
        raise ExportError(f"Unknown component ID: {component_id}")
    protein_rows = pq.read_table(
        proteins_path,
        filters=[("protein_id", "in", protein_ids)],
        columns=["protein_id", "species_id", "original_id"],
    ).to_pylist()
    metadata = {int(row["protein_id"]): row for row in protein_rows}
    edges = read_component_edges(run_root / "components", component_id).to_pylist()
    output_root = run_root / "results" / "graphs"
    output_root.mkdir(parents=True, exist_ok=True)
    output = output_root / f"component={component_id:08d}.graphml"
    temporary = output.with_name(f".{output.name}.{uuid.uuid4().hex}.tmp")
    with temporary.open("w", encoding="utf-8", newline="\n") as handle:
        handle.write('<?xml version="1.0" encoding="UTF-8"?>\n')
        handle.write('<graphml xmlns="http://graphml.graphdrawing.org/xmlns">\n')
        handle.write('<key id="species_id" for="node" attr.name="species_id" attr.type="int"/>\n')
        handle.write(
            '<key id="original_id" for="node" attr.name="original_id" '
            'attr.type="string"/>\n'
        )
        handle.write('<key id="weight" for="edge" attr.name="weight" attr.type="double"/>\n')
        handle.write('<graph id="G" edgedefault="undirected">\n')
        for protein_id in protein_ids:
            row = metadata[protein_id]
            handle.write(f'<node id={quoteattr(str(protein_id))}>')
            handle.write(f'<data key="species_id">{int(row["species_id"])}</data>')
            handle.write(
                f'<data key="original_id">{escape(str(row["original_id"]))}</data></node>\n'
            )
        for edge_id, row in enumerate(edges):
            handle.write(
                f'<edge id="e{edge_id}" source={quoteattr(str(int(row["u"])))} '
                f'target={quoteattr(str(int(row["v"])))}>'
                f'<data key="weight">{float(row["weight"]):.17g}</data></edge>\n'
            )
        handle.write("</graph>\n</graphml>\n")
    os.replace(temporary, output)
    write_json(
        output.with_suffix(".manifest.json"),
        {
            "algorithm_version": "component-graphml-v1",
            "command": command,
            "component_id": component_id,
            "counts": {"nodes": len(protein_ids), "edges": len(edges)},
            "output_sha256": sha256_file(output),
        },
    )
    return output
