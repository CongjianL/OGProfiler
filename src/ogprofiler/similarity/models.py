"""Typed records for the numeric similarity pipeline."""

from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True, slots=True)
class DirectionalHit:
    query_id: int
    target_id: int
    query_species: int
    target_species: int
    bitscore: float
    identity: float
    query_coverage: float
    target_coverage: float
    evalue: float


@dataclass(frozen=True, slots=True)
class NormalizedHit:
    hit: DirectionalHit
    normalized_score: float


@dataclass(frozen=True, slots=True)
class RetainedEdge:
    u: int
    v: int
    u_species: int
    v_species: int
    score_uv: float
    score_vu: float
    weight: float
    coverage: float
    edge_type: str
