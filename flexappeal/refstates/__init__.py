"""Reference-state annotation: what state a trajectory's motifs are drifting towards.

Built to `FLEXAPPEAL_REFSTATES_BUILD_SPEC.md`. Motif annotation itself comes from
BoltzMaker's `sse_comparison`, vendored at a pinned SHA and wrapped rather than
reimplemented; this package adds trajectory iteration, reference ensembles, the
state-similarity score and the membrane analysis that `compare-sse` does not have.
"""

from __future__ import annotations

from .sources import (Catalogue, Ensemble, SourceError, kinase_id,
                      kinase_references, receptor_references, receptor_residues)

__all__ = ["Catalogue", "Ensemble", "SourceError", "kinase_id",
           "kinase_references", "receptor_references", "receptor_residues"]
