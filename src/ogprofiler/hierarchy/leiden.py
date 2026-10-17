"""Small deterministic wrapper around leidenalg partition methods."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import igraph as ig
import leidenalg as la

from ogprofiler.exceptions import HierarchyError


@dataclass(slots=True)
class LeidenCallCounter:
    count: int = 0
    n_iterations: int = 2


@dataclass(frozen=True, slots=True)
class LeidenResult:
    membership: tuple[int, ...]
    quality: float


_PARTITIONS: dict[str, Any] = {
    "rber": la.RBERVertexPartition,
    "rbcv": la.RBConfigurationVertexPartition,
    "cpm": la.CPMVertexPartition,
    "modularity": la.ModularityVertexPartition,
}


def run_leiden(
    graph: ig.Graph,
    gamma: float,
    method: str,
    weights: str | list[float] | None,
    seed: int,
    counter: LeidenCallCounter | None = None,
) -> LeidenResult:
    """Run Leiden once and return only membership and quality."""

    if method not in _PARTITIONS:
        allowed = ", ".join(sorted(_PARTITIONS))
        raise HierarchyError(f"Unknown Leiden method {method!r}; expected: {allowed}")
    if graph.vcount() == 0:
        raise HierarchyError("Leiden requires at least one vertex")
    kwargs: dict[str, Any] = {
        "weights": weights,
        "seed": seed,
        "n_iterations": counter.n_iterations if counter else 2,
    }
    if method != "modularity":
        if gamma <= 0:
            raise HierarchyError("Leiden resolution gamma must be positive")
        kwargs["resolution_parameter"] = gamma
    if counter is not None:
        counter.count += 1
    try:
        partition = la.find_partition(graph, _PARTITIONS[method], **kwargs)
    except (ValueError, ig.InternalError) as error:
        raise HierarchyError(f"Leiden failed at gamma={gamma}: {error}") from error
    return LeidenResult(
        tuple(int(value) for value in partition.membership), float(partition.quality())
    )
