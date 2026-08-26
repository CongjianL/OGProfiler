"""Replaceable homolog-search backends and stage orchestration."""

from ogprofiler.search.base import SearchBackend, SearchParameters, SearchStageResult
from ogprofiler.search.diamond import DiamondBackend
from ogprofiler.search.stage import run_search_stage

__all__ = [
    "DiamondBackend",
    "SearchBackend",
    "SearchParameters",
    "SearchStageResult",
    "run_search_stage",
]
