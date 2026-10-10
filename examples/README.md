# examples

Complete, reproducible runs. Each one is a real simulation whose artefacts are
committed, not a written-up description of one.

| Example | System | Length | Hardware | Wall time |
|---|---|---|---|---|
| [`lysozyme_benzene_100ns`](lysozyme_benzene_100ns/) | T4 lysozyme L99A with benzene in its buried cavity, PDB 181L, 25,463 atoms with water | 100 ns | Apple M2 Ultra, OpenCL | 8h 08m |

## What is in each directory

| File | What it is |
|---|---|
| `config.json` | The complete option set, validated against the registry |
| `flexappeal_lysozyme_benzene_100ns.command` | The self-contained bundle, exactly as the Prepare tab emits it |
| `run/` | The bundle extracted and run in place. Not committed: it holds the full trajectory, and `rebuild.sh` recreates it from the bundle |
| `terminal.log.gz` | Raw `script(1)` capture of the run, escape sequences and all (78 MB → 416 KB) |
| `analyse.log.gz` | The same for the analysis pass |
| `screenshots/*.html` | Terminal captures for the web page, named by phase (`run-*`, `analyse-*`) |
| `output/*.fxa` | The results file, which is what the Example tab renders |
| `rebuild.sh` | Regenerates all of the above |

## Why it is committed rather than generated on demand

The Example tab is served from `output/*.fxa` through the same code path as an
upload: `plots.build_all()` and `plots.build_binding()` build the panels, Mol\*
loads the topology and trajectory out of the same file. So the page cannot drift
from the product. If a change breaks the Analysis tab, this page breaks too, and
it breaks in CI rather than in front of someone.

The example is a complex rather than an apo protein for the same reason. An apo
run draws the whole-system panels and not one of the eleven binding-site panels,
so an apo example would leave most of the Analysis tab untested by the page whose
job is to prove it works.

The alternative -- a hand-written page with static images -- would keep looking
correct long after it had stopped being true.

## Regenerating

```bash
bash examples/lysozyme_benzene_100ns/rebuild.sh
```

Roughly eight and a half hours on an M2 Ultra, eight of them the simulation at
298 ns/day on OpenCL. It needs the development environment, through `pixi run`
rather than by calling `.pixi/envs/default/bin/python` directly: the ligand needs
AM1-BCC charges, which means antechamber and sqm on PATH and `AMBERHOME` set, and
only pixi's activation does that.

One thing to know before committing a new example. The run script prints each
ligand's structure to the console beside its atom count and formal charge,
because the charge is a real check, and that output is what `terminal.log`
captures and this directory commits. So a committed example publishes its
ligand's structure, and only a public molecule belongs here.
