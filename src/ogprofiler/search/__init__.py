"""Replaceable homolog-search backends and stage orchestration."""

from ogprofiler.search.base import SearchBackend, SearchParameters, SearchStageResult
from ogprofiler.search.blast import BlastBackend
from ogprofiler.search.diamond import DiamondBackend
from ogprofiler.search.mmseqs import MmseqsBackend
from ogprofiler.search.stage import run_search_stage

__all__ = [
    "BlastBackend",
    "DiamondBackend",
    "MmseqsBackend",
    "SearchBackend",
    "SearchParameters",
    "SearchStageResult",
    "run_search_stage",
]
