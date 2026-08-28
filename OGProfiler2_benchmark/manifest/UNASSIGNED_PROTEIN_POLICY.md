# B1 unassigned-protein normalization policy

The supplied `benchmark.py` checks the expected protein universe, and the standardized
partition validator requires every input protein exactly once. Some competitor outputs
omit singleton or unassigned proteins. For every method, after parsing the frozen primary
group output, each input protein absent from all reported groups is emitted as its own
group named `TOOL_UNASSIGNED_SINGLETON_<protein_id>`.

This operation does not insert proteins into an inferred family and does not create an
orthology claim. It represents the method's lack of assignment as a singleton, preserves
the complete benchmark universe, and is applied before both official and extended scoring.
Unknown IDs and duplicate assignments remain hard failures. OGProfiler already reports a
complete partition and is unchanged.
