# Protocol

Raw input paths are never modified. Raw command outputs belong in `03_runs`; interoperable tables belong in `04_standardized`; metrics belong in `05_metrics`. Every execution goes through `workflows/run_timed.sh`; failed processes retain provenance. QfO pairwise scoring is prohibited until a validated orthology/OrthoXML representation exists.

Open Orthobench official score is produced only by its supplied `BENCHMARKS/benchmark.py`, invoked with one argument: a one-orthogroup-per-line prediction file. The local extended metrics do not replace it. ARI is deferred because a full protein-universe truth partition is not supplied by RefOGs; VI here is explicitly RefOG-universe restricted.
