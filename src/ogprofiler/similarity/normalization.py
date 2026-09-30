"""OrthoFinder full-run normalized bit-score model."""

from __future__ import annotations

import math
from collections import defaultdict

import numpy as np
from numpy.typing import NDArray
from scipy.optimize import curve_fit

from ogprofiler.similarity.models import DirectionalHit, NormalizedHit


def retain_top_data(
    length_products: list[float], bit_scores: list[float]
) -> tuple[list[float], list[float]]:
    """OF bins: full, disjoint bins only; preserve matrix order for length ties."""

    ordered = sorted(zip(length_products, bit_scores, strict=True), key=lambda item: item[0])
    hit_count = len(ordered)
    if hit_count < 100:
        return [item[0] for item in ordered], [item[1] for item in ordered]
    scale = 1000 if hit_count > 5000 else (200 if hit_count > 1000 else 20)
    top_lengths: list[float] = []
    top_scores: list[float] = []
    for start in range(0, hit_count - scale + 1, scale):
        chunk = ordered[start : start + scale]
        cutoff = float(np.percentile([item[1] for item in chunk], 95))
        for length_product, bit_score in chunk:
            if bit_score >= cutoff:
                top_lengths.append(length_product)
                top_scores.append(bit_score)
    return top_lengths, top_scores


def _fit_parameters(
    length_products: list[float],
    bit_scores: list[float],
    nbs_fallback: str,
) -> tuple[float, float] | None:
    top_lengths, top_scores = retain_top_data(length_products, bit_scores)
    v2_max = (0.0, math.log10(max(bit_scores)))
    if len(top_lengths) < 2:
        return v2_max if nbs_fallback == "v2_max" else None
    # Use OF's solver even for equal length products. A max-score fallback in
    # that case solves a different problem. Solver failures must surface.
    parameters, _ = curve_fit(_loglinear, top_lengths, np.log10(top_scores))
    return float(parameters[0]), float(parameters[1])


def _loglinear(x: NDArray[np.float64], a: float, b: float) -> NDArray[np.float64]:
    return np.asarray(a * np.log10(x) + b, dtype=np.float64)


def max_bitscore_hits(hits: list[DirectionalHit]) -> list[DirectionalHit]:
    """Max HSP per direction, then OF sparse row/column order for fitting."""
    keyed: dict[tuple[int, int], DirectionalHit] = {}
    for hit in hits:
        if hit.query_id == hit.target_id or hit.bitscore <= 0:
            continue
        key = (hit.query_id, hit.target_id)
        previous = keyed.get(key)
        if previous is None or (hit.bitscore, -hit.evalue) > (previous.bitscore, -previous.evalue):
            keyed[key] = hit
    return [keyed[key] for key in sorted(keyed)]


def legacy_nbs(
    hits: list[DirectionalHit],
    protein_lengths: dict[int, int],
    nbs_fallback: str = "v1_zero",
) -> list[NormalizedHit]:
    """Normalize max-HSP matrices using OF3 default full-run semantics."""

    usable = max_bitscore_hits(hits)
    groups: dict[tuple[int, int], list[DirectionalHit]] = defaultdict(list)
    for hit in usable:
        groups[(hit.query_species, hit.target_species)].append(hit)
    parameters: dict[tuple[int, int], tuple[float, float] | None] = {}
    for key, group in groups.items():
        products = [
            float(protein_lengths[item.query_id] * protein_lengths[item.target_id])
            for item in group
        ]
        parameters[key] = _fit_parameters(products, [item.bitscore for item in group], nbs_fallback)

    result: list[NormalizedHit] = []
    for hit in usable:
        params = parameters[(hit.query_species, hit.target_species)]
        if params is None:
            # Too few fit points: OF emits an all-zero matrix.
            continue
        a, b = params
        score = (
            (10 ** (-b) * float(protein_lengths[hit.query_id]) ** (-a))
            * hit.bitscore
            * float(protein_lengths[hit.target_id]) ** (-a)
        )
        result.append(NormalizedHit(hit, score if math.isfinite(score) else 0.0))
    return result


def normalize_hits(
    hits: list[DirectionalHit],
    protein_lengths: dict[int, int],
    method: str,
    nbs_fallback: str = "v1_zero",
) -> list[NormalizedHit]:
    if method == "legacy_nbs":
        return legacy_nbs(hits, protein_lengths, nbs_fallback)
    usable = [hit for hit in hits if hit.query_id != hit.target_id and hit.bitscore > 0]
    if method == "raw_bitscore":
        return [NormalizedHit(hit, hit.bitscore) for hit in usable]
    if method == "length_scaled_bitscore":
        return [
            NormalizedHit(
                hit,
                hit.bitscore / min(protein_lengths[hit.query_id], protein_lengths[hit.target_id]),
            )
            for hit in usable
        ]
    raise ValueError(f"Unknown normalization method: {method}")
