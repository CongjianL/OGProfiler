"""Application exception hierarchy."""


class OGProfilerError(Exception):
    """Base class for expected OGProfiler failures."""


class InputError(OGProfilerError):
    """Input data is missing, malformed, or inconsistent."""


class SearchError(OGProfilerError):
    """Homology search failed."""


class EdgeConstructionError(OGProfilerError):
    """Similarity edge construction failed."""


class ComponentError(OGProfilerError):
    """Connected-component indexing or partitioning failed."""


class HierarchyError(OGProfilerError):
    """Hierarchy inference failed or produced invalid structure."""


class ExportError(OGProfilerError):
    """Final result export failed validation."""


class PhylogenyError(OGProfilerError):
    """Phylogenetic refinement failed."""


class CheckpointError(OGProfilerError):
    """Checkpoint state is invalid or inconsistent."""
