"""The read-only SQLite the explorer serves from. Spec section 7.

The schema is section 7's, with one change: its `references` table is renamed
`reference_structures`. `REFERENCES` is a SQL keyword and
`CREATE TABLE references (...)` is a syntax error in SQLite, so the spec's DDL
does not execute as written. Quoting it everywhere would work and would mean
every query carrying quotes around one table name for the rest of the project.

Written once by the build and only read at serve time, which is what makes
SQLite the right choice here: no write concurrency to arrange, one file to copy
to the droplet, and the per-frame table indexed on the two columns the explorer
filters by.

Per-residue DSSP does not go in here, which is the point section 7 makes: it is
a residues-by-frames matrix and would be a few hundred thousand rows per run in
a table. Section 7 names Parquet for it; this build keeps it in the run's own
analysis JSON instead, as one string per frame of one character per residue.

That is the same saving by the same argument. b2ar_apo's matrix is 499 frames
by 322 residues, 160,678 cells, which is 162 KiB on disk and 10 KiB over the
wire gzipped, against the 160,678 SQLite rows the spec is avoiding. Parquet
would hold it in a comparable space, and pyarrow is already installed, so the
reason for not using it is not cost: it is that the analysis JSON is already
the one artefact written per run and already read by the page that draws the
ribbon. A second file would be another thing to write, to rsync and to keep in
step with the first, for a matrix this size. Recorded as a deviation in
FINDINGS section 4 rather than left as a difference nobody wrote down.
"""

from __future__ import annotations

import json
import statistics
import sqlite3

from .state import DECISIVE_S
from contextlib import contextmanager
from pathlib import Path
from typing import Any, Iterator

SCHEMA = """
CREATE TABLE IF NOT EXISTS systems (
    system_id TEXT PRIMARY KEY, family TEXT, uniprot TEXT, name TEXT,
    environment TEXT, n_atoms INTEGER, forcefield TEXT, notes TEXT);

CREATE TABLE IF NOT EXISTS reference_structures (
    ref_id TEXT PRIMARY KEY, system_id TEXT REFERENCES systems,
    pdb TEXT, chain TEXT, source TEXT, source_state TEXT,
    our_state TEXT, resolution REAL, metrics_json TEXT);

CREATE TABLE IF NOT EXISTS motifs (
    motif_id TEXT PRIMARY KEY, system_id TEXT REFERENCES systems,
    name TEXT, token TEXT, residues_json TEXT, generic_labels_json TEXT,
    source TEXT);

CREATE TABLE IF NOT EXISTS runs (
    run_id TEXT PRIMARY KEY, system_id TEXT REFERENCES systems,
    start_state TEXT, start_model TEXT, start_model_rung INTEGER,
    arm TEXT, replica INTEGER, tier INTEGER, length_ns REAL,
    platform TEXT, ns_per_day REAL, status TEXT);

CREATE TABLE IF NOT EXISTS motif_frames (
    run_id TEXT REFERENCES runs, frame INTEGER, time_ns REAL,
    motif_id TEXT REFERENCES motifs, metric TEXT, value REAL,
    s_score REAL, label TEXT,
    PRIMARY KEY (run_id, frame, motif_id, metric));

CREATE TABLE IF NOT EXISTS membrane (
    run_id TEXT REFERENCES runs, time_ns REAL, apl REAL,
    thickness REAL, tilt REAL, core_waters INTEGER);

CREATE INDEX IF NOT EXISTS ix_mf_run ON motif_frames(run_id);
CREATE INDEX IF NOT EXISTS ix_mf_motif ON motif_frames(motif_id);
CREATE INDEX IF NOT EXISTS ix_ref_sys ON reference_structures(system_id);
CREATE INDEX IF NOT EXISTS ix_runs_sys ON runs(system_id);
"""

# Section 6.1's motif colour tokens. Shared by every view, so they live with the
# data rather than in a stylesheet: the explorer reads a motif's token from the
# same row it reads the motif's residues from, and a motif cannot then be one
# colour on the timeline and another on the viewer.
# Section 6.1 fixes these in a table, and this is that table wherever it can
# be. It cannot be everywhere, because the motifs the build ended up scoring
# are not the ones the table names: the spec gives tokens to the R-spine, the
# ionic lock, the PIF connector and the toggle switch, all four of which this
# build measures but as pocket states, micro-switch metrics and corroborating
# rules rather than as motifs with reference ensembles of their own, so none of
# them has a row here or anything to colour. Meanwhile `c.l`, `TM5` and `TM7`
# are scored as motifs and the table does not name them.
#
# So every assignment in section 6.1 that applies to a motif this build renders
# is kept, and only the tokens left free by a motif the build does not render
# are reused. Kinase: 1 to 5 are the spec's, and `c.l` takes 6 from the R-spine.
# Receptor: 1, 3 and 6 are the spec's, `TM5` takes 2 from the ionic lock and
# `TM7` takes 4 from the PIF connector, leaving 5 unused because nothing
# replaces the toggle switch.
#
# A token for a motif that does not exist is harmless; a motif with no token is
# not, because it draws in whatever colour the chart library offers next, which
# is the opposite of a shared palette.
MOTIF_TOKENS = {
    "kinase": {"alphaC": "--motif-1", "xDFG": "--motif-2", "A-loop": "--motif-3",
               "g.l": "--motif-4", "hinge": "--motif-5", "c.l": "--motif-6"},
    "gpcr": {"TM6": "--motif-1", "TM5": "--motif-2", "NPxxY": "--motif-3",
             "TM7": "--motif-4", "H8": "--motif-6"},
}


@contextmanager
def connect(path: Path, *, read_only: bool = False) -> Iterator[sqlite3.Connection]:
    """A connection with foreign keys on and rows addressable by name.

    Foreign keys are off by default in SQLite, which makes the REFERENCES
    clauses in the schema decoration rather than constraint. Turned on here so a
    motif frame cannot name a run that does not exist.
    """
    if read_only:
        connection = sqlite3.connect(f"file:{path}?mode=ro", uri=True)
    else:
        Path(path).parent.mkdir(parents=True, exist_ok=True)
        connection = sqlite3.connect(str(path))
    connection.row_factory = sqlite3.Row
    connection.execute("PRAGMA foreign_keys = ON")
    try:
        yield connection
        if not read_only:
            connection.commit()
    finally:
        connection.close()


def create(path: Path) -> None:
    with connect(path) as connection:
        connection.executescript(SCHEMA)


def _upsert(connection, table: str, keys: tuple[str, ...], row: dict[str, Any]) -> None:
    """Insert or replace, so re-ingesting a run is idempotent.

    The build re-runs. An ingest that appended would silently double every
    per-frame row and halve every mean computed from it.
    """
    columns = sorted(row)
    placeholders = ", ".join("?" for _ in columns)
    updates = ", ".join(f"{c} = excluded.{c}" for c in columns if c not in keys)
    sql = (f"INSERT INTO {table} ({', '.join(columns)}) VALUES ({placeholders}) "
           f"ON CONFLICT({', '.join(keys)}) DO UPDATE SET {updates}"
           if updates else
           f"INSERT OR REPLACE INTO {table} ({', '.join(columns)}) "
           f"VALUES ({placeholders})")
    connection.execute(sql, [row[c] for c in columns])


def put_system(connection, system_id: str, **fields: Any) -> None:
    _upsert(connection, "systems", ("system_id",),
            {"system_id": system_id, **fields})


def put_run(connection, run_id: str, system_id: str, **fields: Any) -> None:
    _upsert(connection, "runs", ("run_id",),
            {"run_id": run_id, "system_id": system_id, **fields})


def put_motif(connection, motif_id: str, system_id: str, name: str,
              residues: list[int], *, token: str = "",
              generic_labels: dict | None = None, source: str = "") -> None:
    _upsert(connection, "motifs", ("motif_id",), {
        "motif_id": motif_id, "system_id": system_id, "name": name,
        "token": token, "residues_json": json.dumps(residues),
        "generic_labels_json": json.dumps(generic_labels or {}),
        "source": source})


def put_reference(connection, ref_id: str, system_id: str, **fields: Any) -> None:
    metrics = fields.pop("metrics", None)
    if metrics is not None:
        fields["metrics_json"] = json.dumps(metrics)
    _upsert(connection, "reference_structures", ("ref_id",),
            {"ref_id": ref_id, "system_id": system_id, **fields})


def put_frames(connection, run_id: str, rows: list[dict[str, Any]]) -> int:
    """Per-frame motif values, in one transaction.

    `executemany` rather than a loop of `_upsert`: a 50 ns run at one frame per
    100 ps with six motifs is 3,000 rows per metric, and per-statement commits
    on that are the difference between a second and a minute.
    """
    prepared = [(run_id, int(r["frame"]), float(r["time_ns"]), r["motif_id"],
                 r["metric"], r.get("value"), r.get("s_score"), r.get("label"))
                for r in rows]
    connection.executemany(
        "INSERT INTO motif_frames (run_id, frame, time_ns, motif_id, metric, "
        "value, s_score, label) VALUES (?, ?, ?, ?, ?, ?, ?, ?) "
        "ON CONFLICT(run_id, frame, motif_id, metric) DO UPDATE SET "
        "value = excluded.value, s_score = excluded.s_score, "
        "label = excluded.label, time_ns = excluded.time_ns", prepared)
    return len(prepared)


def put_membrane(connection, run_id: str, chunks: list[dict[str, Any]]) -> int:
    """Membrane chunks. Cleared for the run first, because a chunk has no key.

    Section 7 gives this table no primary key, so a re-ingest cannot upsert and
    would append instead. Deleting the run's rows first keeps a re-run
    idempotent in the same way every other table here is.
    """
    connection.execute("DELETE FROM membrane WHERE run_id = ?", (run_id,))
    prepared = [(run_id, c.get("start_ns"), c.get("apl_a2"),
                 c.get("thickness_a"), c.get("tilt_deg"),
                 c.get("core_waters")) for c in chunks]
    connection.executemany(
        "INSERT INTO membrane (run_id, time_ns, apl, thickness, tilt, "
        "core_waters) VALUES (?, ?, ?, ?, ?, ?)", prepared)
    return len(prepared)


# ===========================================================================
#  Ingest
# ===========================================================================


def ingest_analysis(path: Path, analysis: dict[str, Any], *,
                    system_id: str, run_id: str,
                    family: str, environment: str,
                    run_fields: dict[str, Any] | None = None,
                    system_fields: dict[str, Any] | None = None) -> dict[str, int]:
    """Load one `tools/analyse_run.py` payload into the store.

    Returns the row counts, which is what the caller reports: an ingest that
    wrote nothing and an ingest that wrote everything both finish silently
    otherwise.
    """
    tokens = MOTIF_TOKENS.get(family, {})
    counts = {"motifs": 0, "frames": 0, "references": 0}

    with connect(path) as connection:
        connection.executescript(SCHEMA)
        # `system_fields` carries section 7's `uniprot`, `n_atoms` and
        # `forcefield`, which the schema has had from the start and which
        # nothing had ever written: all three were NULL for both systems, and
        # that is why section 6.2's card could not show an atom count. Only
        # keys the caller could actually produce are passed, so a system whose
        # build manifest is missing leaves the column null rather than carrying
        # a zero that reads like a measurement.
        put_system(connection, system_id, family=family,
                   name=analysis.get("name", ""),
                   environment=environment,
                   **{k: v for k, v in (system_fields or {}).items()
                      if v is not None},
                   notes=json.dumps({"core": analysis.get("core", {}),
                                     "notes": analysis.get("notes", [])}))
        # The caller's fields win, so a run still in flight keeps its own
        # status rather than being overwritten with "complete" by the default.
        fields = {"status": "complete",
                  "length_ns": (analysis.get("frames") or [{}])[-1].get("time_ns")}
        fields.update(run_fields or {})
        put_run(connection, run_id, system_id, **fields)

        for name, residues in (analysis.get("motif_residues") or {}).items():
            put_motif(connection, f"{system_id}:{name}", system_id, name,
                      residues, token=tokens.get(name, ""),
                      source=analysis.get("core", {}).get("source", ""))
            counts["motifs"] += 1

        for state, ensemble in (analysis.get("ensembles") or {}).items():
            for code in ensemble.get("codes", []):
                put_reference(connection, f"{system_id}:{code}", system_id,
                              pdb=code, our_state=state,
                              source="KLIFS" if family == "kinase" else "GPCRdb")
                counts["references"] += 1

        rows = []
        for index, frame in enumerate(analysis.get("frames") or []):
            for motif in frame.get("motifs", []):
                name = motif["motif"]
                base = {"frame": index, "time_ns": frame["time_ns"],
                        "motif_id": f"{system_id}:{name}",
                        "s_score": motif.get("s"), "label": motif.get("label")}
                # One row per metric, which is what the schema's primary key
                # expects: the s score repeats across a motif's metrics and the
                # value is the metric's own.
                for metric, value in (("d_act_a", motif.get("d_act_a")),
                                      ("d_inact_a", motif.get("d_inact_a"))):
                    if value is not None:
                        rows.append({**base, "metric": metric, "value": value})
        counts["frames"] = put_frames(connection, run_id, rows) if rows else 0
    return counts


def summarise(path: Path, run_id: str) -> dict[str, Any]:
    """A run's headline numbers, for the Examples card in section 6.2.

    The one-line finding on that card has to be generated from the data rather
    than written by hand, so the numbers it needs come from here: each motif's
    first and last values, and its final-tenth mean.
    """
    with connect(path, read_only=True) as connection:
        run = connection.execute(
            "SELECT * FROM runs WHERE run_id = ?", (run_id,)).fetchone()
        if run is None:
            return {}
        motifs = connection.execute(
            "SELECT DISTINCT motif_id FROM motif_frames WHERE run_id = ?",
            (run_id,)).fetchall()
        span = connection.execute(
            "SELECT MIN(time_ns) AS lo, MAX(time_ns) AS hi FROM motif_frames "
            "WHERE run_id = ?", (run_id,)).fetchone()

        out: dict[str, Any] = {"run": dict(run), "motifs": {},
                               "span_ns": [span["lo"], span["hi"]]}
        if span["hi"] is None:
            return out
        # The final tenth, which is the window section 6.2's motif table uses.
        cut = span["hi"] - 0.1 * (span["hi"] - (span["lo"] or 0.0))
        for row in motifs:
            motif_id = row["motif_id"]
            first = connection.execute(
                "SELECT s_score, label FROM motif_frames WHERE run_id = ? AND "
                "motif_id = ? ORDER BY frame LIMIT 1", (run_id, motif_id)).fetchone()
            last = connection.execute(
                "SELECT s_score, label FROM motif_frames WHERE run_id = ? AND "
                "motif_id = ? ORDER BY frame DESC LIMIT 1",
                (run_id, motif_id)).fetchone()
            tail = connection.execute(
                "SELECT AVG(s_score) AS mean FROM motif_frames WHERE run_id = ? "
                "AND motif_id = ? AND time_ns >= ?",
                (run_id, motif_id, cut)).fetchone()
            # The halves and the scatter, which is what separates a motif
            # that moved from one that is noisy. A first-half against
            # second-half comparison on its own calls a step change a drift,
            # so the standard deviation goes beside it: a change smaller than
            # the scatter is not a trend, and the reader needs both numbers to
            # see that. Computed here rather than in the report, because the
            # page reads this too and the two must not derive it differently.
            scores = [row["s_score"] for row in connection.execute(
                "SELECT DISTINCT frame, s_score FROM motif_frames "
                "WHERE run_id = ? AND motif_id = ? AND s_score IS NOT NULL "
                "ORDER BY frame", (run_id, motif_id))]
            halves: dict[str, Any] = {}
            if len(scores) >= 4:
                # Not `cut`: that name already holds the final-tenth time
                # boundary, computed once before this loop and used by the
                # `tail` query above. Shadowing it here replaced 13.14 ns with
                # a frame count of 73, so every motif after the first asked for
                # frames past 73 ns of a 14.6 ns run and got NULL back. The
                # first motif was computed before the overwrite and looked
                # fine, which is why the page still rendered a number.
                midpoint = len(scores) // 2
                halves = {
                    "first_half_mean_s": round(
                        sum(scores[:midpoint]) / midpoint, 4),
                    "second_half_mean_s": round(
                        sum(scores[midpoint:]) / len(scores[midpoint:]), 4),
                    "sd_s": round(statistics.pstdev(scores), 4),
                    "n_scored_frames": len(scores),
                }

            out["motifs"][motif_id] = {
                "start_s": first["s_score"] if first else None,
                "start_label": first["label"] if first else None,
                "final_s": last["s_score"] if last else None,
                "final_label": last["label"] if last else None,
                "final_tenth_mean_s": tail["mean"] if tail else None,
                **halves,
            }
        return out


def whole_score(path: Path, run_id: str) -> dict[str, Any]:
    """The residue-weighted whole-protein score at the start and the end.

    Recomputed from the stored per-motif scores rather than stored itself, so
    there is one definition of the weighting and it cannot drift between what
    the analysis wrote and what a page shows. The weights are the motifs' own
    residue counts, which is what makes a twenty-residue loop count for more
    than a single toggle tryptophan.

    Reported with the number of motifs behind it, because a mean over two
    motifs and a mean over six are not the same claim and the figure alone
    does not say which it is.
    """
    import json as _json

    summary = summarise(path, run_id)
    if not summary.get("motifs"):
        return {}

    with connect(path, read_only=True) as connection:
        weights = {}
        for row in connection.execute(
                "SELECT motif_id, residues_json FROM motifs"):
            try:
                weights[row["motif_id"]] = len(
                    _json.loads(row["residues_json"] or "[]"))
            except (TypeError, _json.JSONDecodeError):
                weights[row["motif_id"]] = 0

    def weighted(key: str) -> float | None:
        pairs = [(values[key], weights.get(motif_id, 0))
                 for motif_id, values in summary["motifs"].items()
                 if values.get(key) is not None and weights.get(motif_id, 0) > 0]
        if not pairs:
            return None
        total = sum(n for _, n in pairs)
        return sum(v * n for v, n in pairs) / total

    start = weighted("start_s")
    final = weighted("final_tenth_mean_s")
    out: dict[str, Any] = {
        "n_motifs": sum(1 for v in summary["motifs"].values()
                        if v.get("final_tenth_mean_s") is not None),
        "span_ns": summary.get("span_ns"),
    }
    if start is not None:
        out["start_s"] = start
    if final is not None:
        out["final_s"] = final
        # `state.DECISIVE_S`, not a literal. This band is the difference
        # between "active" and "intermediate" on every page, and it was written
        # out three times: here, in `web.held`, and as the shaded region in the
        # run template. Three copies of one threshold is three chances for the
        # summary label to disagree with the per-frame label it summarises,
        # silently, with nothing on the page to say which one moved.
        out["label"] = ("active" if final > DECISIVE_S else
                        "inactive" if final < -DECISIVE_S else "intermediate")
    return out
