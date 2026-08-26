"""Frozen-V1-compatible normalized bit-score model."""

from __future__ import annotations

import math
from collections import defaultdict

from ogprofiler.similarity.models import DirectionalHit, NormalizedHit


def retain_top_data(
    length_products: list[float], bit_scores: list[float]
) -> tuple[list[float], list[float]]:
    """Reproduce V1 length sorting, bin widths, overlap, and 95th percentile."""

    ordered = sorted(zip(length_products, bit_scores, strict=True), key=lambda item: item[0])
    hit_count = len(ordered)
    if hit_count < 100:
        return [item[0] for item in ordered], [item[1] for item in ordered]
    scale = 1000 if hit_count > 5000 else (200 if hit_count > 1000 else 20)
    top_lengths: list[float] = []
    top_scores: list[float] = []
    for start in range(0, hit_count, scale):
        chunk = ordered[start : start + scale + 1]
        scores = sorted(item[1] for item in chunk)
        position = 0.95 * (len(scores) - 1)
        lower = math.floor(position)
        upper = math.ceil(position)
        cutoff = scores[lower] + (position - lower) * (scores[upper] - scores[lower])
        for length_product, bit_score in chunk:
            if bit_score >= cutoff:
                top_lengths.append(length_product)
                top_scores.append(bit_score)
    return top_lengths, top_scores


def _fit_parameters(length_products: list[float], bit_scores: list[float]) -> tuple[float, float]:
    top_lengths, top_scores = retain_top_data(length_products, bit_scores)
    fallback = (0.0, math.log10(max(bit_scores)))
    if len(top_lengths) < 2 or len(set(float(value) for value in top_lengths)) < 2:
        return fallback
    x_values = [math.log10(value) for value in top_lengths]
    y_values = [math.log10(value) for value in top_scores]
    x_mean = sum(x_values) / len(x_values)
    y_mean = sum(y_values) / len(y_values)
    denominator = sum((value - x_mean) ** 2 for value in x_values)
    if denominator <= 0:
        return fallback
    a = (
        sum(
            (x_value - x_mean) * (y_value - y_mean)
            for x_value, y_value in zip(x_values, y_values, strict=True)
        )
        / denominator
    )
    b = y_mean - a * x_mean
    return (a, b) if math.isfinite(a) and math.isfinite(b) else fallback


def legacy_nbs(hits: list[DirectionalHit], protein_lengths: dict[int, int]) -> list[NormalizedHit]:
    """Normalize each directed species-pair group using frozen V1 semantics."""

    usable = [hit for hit in hits if hit.query_id != hit.target_id and hit.bitscore > 0]
    groups: dict[tuple[int, int], list[DirectionalHit]] = defaultdict(list)
    for hit in usable:
        groups[(hit.query_species, hit.target_species)].append(hit)
    parameters: dict[tuple[int, int], tuple[float, float]] = {}
    for key, group in groups.items():
        products = [
            float(protein_lengths[item.query_id] * protein_lengths[item.target_id])
            for item in group
        ]
        parameters[key] = _fit_parameters(products, [item.bitscore for item in group])

    result: list[NormalizedHit] = []
    for hit in usable:
        a, b = parameters[(hit.query_species, hit.target_species)]
        product = protein_lengths[hit.query_id] * protein_lengths[hit.target_id]
        denominator = (10**b) * (product**a)
        score = hit.bitscore / denominator if denominator > 0 else 0.0
        result.append(NormalizedHit(hit, score if math.isfinite(score) else 0.0))
    return result


def normalize_hits(
    hits: list[DirectionalHit], protein_lengths: dict[int, int], method: str
) -> list[NormalizedHit]:
    if method == "legacy_nbs":
        return legacy_nbs(hits, protein_lengths)
    usable = [hit for hit in hits if hit.query_id != hit.target_id and hit.bitscore > 0]
    if method == "raw_bitscore":
        return [NormalizedHit(hit, hit.bitscore) for hit in usable]
    if method == "length_scaled_bitscore":
        return [
            NormalizedHit(
                hit,
                hit.bitscore
                / min(protein_lengths[hit.query_id], protein_lengths[hit.target_id]),
            )
            for hit in usable
        ]
    raise ValueError(f"Unknown normalization method: {method}")
