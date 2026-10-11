# FlexAppeal: reference-state build

Section 12.2 of `FLEXAPPEAL_REFSTATES_BUILD_SPEC.md` gives this file's contents.
Every line the spec writes is here. Two are **changed**, both because the user
overrode them during the build, and both say so where they appear. Some are
added: the units and omitted-key rules, and the security section, which the
user gave during the build and which the spec does not contain.

- `FLEXAPPEAL_REFSTATES_BUILD_SPEC.md` is authoritative. Read the clause before
  changing what it describes: several hours of this build were spent correcting
  decisions that a table in the spec had already made.
- Never stop except at Gates G1 to G6; follow the gate protocol.
- Log to `BUILD_LOG.md` and `DECISIONS.md`. Every stage writes
  `data/manifests/<stage>.jsonl` and is resumable.
  **Changed:** the spec also lists `PROGRESS.md`. The user waived that on
  10 October ("Ignore rule 5") and `tools/progress.py` was deleted. Do not
  reinstate it without being asked.
- Commit at the end of every phase.
  **Changed:** the spec says "never push". The user authorised pushing to
  `github.com/bellcheddar/FlexAppeal` on 10 October and it has been pushed
  since. The repository is public, so the security rules below are not
  optional.
- Motif annotation comes from `vendor/boltzmaker/sse_comparison` at the pinned
  SHA. Wrap it, never edit it; a needed change goes in a FlexAppeal subclass or
  adapter and is logged in `DECISIONS.md` as a candidate upstream change.
- No hand-typed residue numbers. KLIFS pocket positions and GPCRdb generic
  numbers only, resolved at run time.
- Hardware settings come from `config/tuning.toml`, written by
  `tools/hwprobe.py`. The OpenMM platform is chosen by `tools/run_md.py` against
  the real system rather than read from a table.
- The UI uses FlexAppeal's existing design system in `flexappeal/static/brand.css`.
  Motif and state colours come from the `--motif-*` and `--state-*` tokens,
  which section 6.1 fixes; do not pick colours in a chart.
- House style: British English, no em dashes, and none of "groundbreaking",
  "revolutionary", "paradigm-shifting", "game-changing", "cutting-edge",
  "unprecedented", "seamless", "leverage" as a verb, or "delve".
- Units are angstrom and nanoseconds. Strain is the excess over frame 0.
- Omit a key you cannot produce. Never write `0`, `-1` or `null` as a stand-in.

## Security, because the repository is public

- Test fixtures use public structures only: lysozyme, trypsin, glycophorin, and
  the public PDB entries this build already cites. Never a real compound name
  or coordinates from the GLP1R/GIPR work, which is why DeepFlex is private.
- Generate fixtures in the test; do not commit them.
- No structure strings in any member or manifest field. The roots to avoid in a
  field name are `smiles`, `inchi`, `molblock`, `molfile`, `tripos`, `sdf`.
- No `.sdf`, `.mol`, `.mol2` or `.pdbqt` members, no force-field libraries, no
  trajectory files.
