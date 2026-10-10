"""The state-similarity score: how active-like one motif is in one frame.

Spec section 5.6. Two independent readings of the same motif are kept side by
side and never merged:

    the geometric one   RMSD to the nearest active reference against RMSD to the
                        nearest inactive one, as the signed score `s`
    the rule-based one  the DFG call, the salt bridge, the ionic lock: a
                        published criterion with a published threshold

When they agree the label is `active` or `inactive`. When they disagree the label
is `conflict`, which is the point of measuring both. A frame whose DFG geometry
reads out but whose A-loop has drifted to an inactive reference is a real and
interesting thing, and a scorer that averaged the two would report an
intermediate and hide it.

`s` is a normalised difference rather than a raw RMSD gap so that motifs of very
different sizes land on one axis. A three-residue motif and a twenty-residue loop
both saturate towards +/-1, and neither dominates the whole-protein mean through
having larger absolute RMSDs.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

# Section 5.6. Inside this band of `s` the geometry is not claiming either state,
# so the label is `intermediate` whatever the rule says.
DECISIVE_S = 0.2

# Display smoothing only; the raw per-frame values are what gets stored.
SMOOTH_WINDOW_NS = 2.0


@dataclass
class MotifState:
    """One motif in one frame, read both ways.

    `d_act_a` and `d_inact_a` travel with the score because `s` alone cannot
    distinguish a motif sitting between two well-separated references from one
    where both references are equally far away and the whole comparison is
    meaningless.
    """

    motif: str
    d_act_a: float | None = None
    d_inact_a: float | None = None
    rule_label: str | None = None
    rule_metric: str | None = None
    n_residues: int = 0
    nearest_active: str | None = None
    nearest_inactive: str | None = None

    @property
    def s(self) -> float | None:
        """+1 active-like, -1 inactive-like, 0 equidistant."""
        if self.d_act_a is None or self.d_inact_a is None:
            return None
        total = self.d_inact_a + self.d_act_a
        if total <= 0.0:
            return None
        return (self.d_inact_a - self.d_act_a) / total

    @property
    def geometric_label(self) -> str | None:
        value = self.s
        if value is None:
            return None
        if value > DECISIVE_S:
            return "active"
        if value < -DECISIVE_S:
            return "inactive"
        return "intermediate"

    @property
    def label(self) -> str | None:
        """The combined call, with disagreement shown rather than resolved."""
        geometric = self.geometric_label
        if geometric is None:
            return None
        rule = rule_state(self.rule_metric or "", self.rule_label)
        if rule not in ("active", "inactive"):
            # Nothing to agree or disagree with: no rule, an unresolved one, or
            # one whose own reading is intermediate. The geometry stands alone
            # and `rule_label` shows the caller it was not corroborated.
            return geometric
        if geometric == "intermediate":
            return "intermediate"
        return geometric if geometric == rule else "conflict"

    def as_row(self) -> dict[str, Any]:
        """The stored form. Keys whose value is unknown are omitted, not zeroed."""
        row: dict[str, Any] = {"motif": self.motif, "n_residues": self.n_residues}
        for key, value in (("d_act_a", self.d_act_a),
                           ("d_inact_a", self.d_inact_a),
                           ("s", self.s),
                           ("rule_label", self.rule_label),
                           ("rule_metric", self.rule_metric),
                           ("rule_state", rule_state(self.rule_metric or "",
                                                     self.rule_label)),
                           ("geometric_label", self.geometric_label),
                           ("label", self.label),
                           ("nearest_active", self.nearest_active),
                           ("nearest_inactive", self.nearest_inactive)):
            if value is not None:
                row[key] = round(value, 4) if isinstance(value, float) else value
        return row


# How each rule-based vocabulary maps onto the two states, keyed by the metric
# as well as the word. A single shared table was wrong and quietly so: `broken`
# means inactive for the R-spine, whose four side chains stack in the active
# kinase, and active for the GPCR ionic lock, which holds the inactive receptor
# shut. One dict mapping the word alone would have labelled every broken lock
# inactive and then agreed with the geometry for the wrong reason.
#
# So direction belongs to the metric, not to the word, and each entry here is
# named for the metric that produces it.
RULE_VOCABULARY = {
    "dfg": {"in": "active", "out": "inactive", "other": "intermediate"},
    "alphac": {"in": "active", "out": "inactive"},
    "r_spine": {"assembled": "active", "broken": "inactive"},
    "ionic_lock": {"intact": "inactive", "broken": "active"},
    "cleft": {"open": "active", "closed": "inactive"},
}


def rule_state(metric: str, label: str | None) -> str | None:
    """The state a rule label implies, read in its own metric's vocabulary.

    An unknown metric or label returns None rather than a guess, which makes the
    motif fall back to geometry alone and say so, instead of agreeing or
    conflicting on a mapping nobody defined.
    """
    if label is None or label == "unresolved":
        return None
    if label in ("active", "inactive", "intermediate"):
        return label
    return RULE_VOCABULARY.get(metric, {}).get(label)


@dataclass
class FrameState:
    """Every motif in one frame, plus the residue-weighted whole-protein score."""

    time_ns: float
    motifs: list[MotifState] = field(default_factory=list)

    @property
    def whole_s(self) -> float | None:
        """Residue-weighted mean of the motif scores.

        Weighted by residue count because a twenty-residue activation loop is
        more of the protein than a single toggle tryptophan, and an unweighted
        mean would let the smallest motifs set the headline number.
        """
        pairs = [(m.s, m.n_residues) for m in self.motifs
                 if m.s is not None and m.n_residues > 0]
        if not pairs:
            return None
        weight = sum(n for _, n in pairs)
        return sum(s * n for s, n in pairs) / weight

    @property
    def whole_label(self) -> str | None:
        value = self.whole_s
        if value is None:
            return None
        if value > DECISIVE_S:
            return "active"
        if value < -DECISIVE_S:
            return "inactive"
        return "intermediate"

    @property
    def n_conflicts(self) -> int:
        return sum(1 for m in self.motifs if m.label == "conflict")

    def as_row(self) -> dict[str, Any]:
        row: dict[str, Any] = {"time_ns": round(self.time_ns, 3),
                               "motifs": [m.as_row() for m in self.motifs],
                               "n_conflicts": self.n_conflicts}
        if self.whole_s is not None:
            row["whole_s"] = round(self.whole_s, 4)
            row["whole_label"] = self.whole_label
        return row


def nearest(distances: dict[str, float]) -> tuple[str | None, float | None]:
    """The closest reference and its distance, or (None, None) if there are none."""
    if not distances:
        return None, None
    name = min(distances, key=lambda key: distances[key])
    return name, distances[name]


def score_motif(motif: str, *, active: dict[str, float],
                inactive: dict[str, float], n_residues: int,
                rule_label: str | None = None,
                rule_metric: str | None = None) -> MotifState:
    """One motif's state, from its RMSD to each reference and its rule call.

    `active` and `inactive` map a reference's PDB code to this frame's motif RMSD
    against it. Keyed by code rather than passed as bare lists so the reported
    nearest reference is traceable back to a structure a reader can open.
    """
    act_name, act_value = nearest(active)
    inact_name, inact_value = nearest(inactive)
    return MotifState(motif=motif, d_act_a=act_value, d_inact_a=inact_value,
                      rule_label=rule_label, rule_metric=rule_metric,
                      n_residues=n_residues,
                      nearest_active=act_name, nearest_inactive=inact_name)


def smooth(values: list[float | None], frame_ps: float,
           window_ns: float = SMOOTH_WINDOW_NS) -> list[float | None]:
    """A centred moving mean, for display only (section 5.6).

    Gaps stay gaps: a frame whose own value is None comes back None rather than
    being filled from its neighbours, because an interpolated point on a state
    trace reads as a measurement that was never taken.
    """
    half = max(1, int(round((window_ns * 1000.0 / frame_ps) / 2)))
    out: list[float | None] = []
    for index, value in enumerate(values):
        if value is None:
            out.append(None)
            continue
        lo, hi = max(0, index - half), min(len(values), index + half + 1)
        window = [v for v in values[lo:hi] if v is not None]
        out.append(sum(window) / len(window))
    return out


def transitions(labels: list[str | None], times_ns: list[float],
                *, min_dwell_ns: float = 1.0) -> list[dict[str, Any]]:
    """Label changes that persist, which is what a state transition means.

    A bare label-change detector on per-frame data reports every excursion that
    one frame of thermal noise produced. A change only counts here if the new
    label holds for `min_dwell_ns`, so what comes back is the handful of real
    events rather than hundreds of flickers.
    """
    runs: list[tuple[str | None, float, float]] = []
    for label, time_ns in zip(labels, times_ns):
        if runs and runs[-1][0] == label:
            runs[-1] = (label, runs[-1][1], time_ns)
        else:
            runs.append((label, time_ns, time_ns))

    kept = [r for r in runs if r[2] - r[1] >= min_dwell_ns or r is runs[0]]
    events = []
    for before, after in zip(kept, kept[1:]):
        if before[0] == after[0]:
            continue
        events.append({"at_ns": round(after[1], 3), "from": before[0],
                       "to": after[0],
                       "dwell_before_ns": round(before[2] - before[1], 3)})
    return events
