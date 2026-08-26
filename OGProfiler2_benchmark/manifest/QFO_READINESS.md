# QFO readiness

**NOT_READY_FOR_OFFICIAL_QFO**

A pairwise stream exists only as `CO_ORTHOLOG_CANDIDATE` (`src/ogprofiler/orthology/engine.py:OrthologCandidate`), generated across children of `SPECIATION_LIKE`/`POLYTOMY` nodes. Those labels are network species-overlap heuristics (`src/ogprofiler/evolution/network.py:annotate_network_events`); they are not reconciliation-derived speciation/duplication calls. The README explicitly states that the network hierarchy is not a gene tree and that `network_event` must remain separate from `phylo_event`.

Consequently, current output does not provide validated ortholog/paralog relations, internal-node speciation/duplication semantics, or an OrthoXML exporter. Terminal family membership will not be pair-expanded for QfO.
