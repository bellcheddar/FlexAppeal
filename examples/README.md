# examples

Complete, reproducible runs. Each one is a real simulation whose artefacts are
committed, not a written-up description of one.

| Example | System | Length | Hardware | Wall time |
|---|---|---|---|---|
| [`lysozyme_10ns`](lysozyme_10ns/) | Hen egg-white lysozyme, PDB 1AKI, 19,433 atoms with water | 10 ns | Apple M1 Max, OpenCL | ~75 min |

## What is in each directory

| File | What it is |
|---|---|
| `config.json` | The complete option set, validated against the registry |
| `flexappeal_lysozyme_10ns.command` | The self-contained bundle, exactly as the Prepare tab emits it |
| `run/` | The bundle extracted and run in place. Not committed: it holds the full trajectory, and `rebuild.sh` recreates it from the bundle |
| `terminal.log.gz` | Raw `script(1)` capture of the run, escape sequences and all (3 MB → 80 KB) |
| `analyse.log.gz` | The same for the analysis pass |
| `screenshots/*.html` | Terminal captures for the web page, named by phase (`run-*`, `analyse-*`) |
| `output/*.fxa` | The results file, which is what the Example tab renders |
| `media/lysozyme.webp` | The rotating cartoon on the Example tab: animated WebP with a real alpha channel, keyed from a flat-background source |
| `rebuild.sh` | Regenerates all of the above |

## Why it is committed rather than generated on demand

The Example tab is served from `output/*.fxa` through the same code path as an
upload: `plots.build_all()` builds the panels, Mol\* loads the topology and
trajectory out of the same file. So the page cannot drift from the product. If
a change breaks the Analysis tab, this page breaks too, and it breaks in CI
rather than in front of someone.

The alternative -- a hand-written page with static images -- would keep looking
correct long after it had stopped being true.

## Regenerating

```bash
bash examples/lysozyme_10ns/rebuild.sh
```

The animation is only re-keyed if its source video is to hand (set
`ANIMATION_SOURCE`, default `~/Desktop/lysozyme.mp4`); otherwise the committed
`media/` is left alone, since the keyed WebP is the deliverable and the source
is not in the repository.

Roughly ninety minutes on an M1 Max, most of it the simulation. It needs the
development environment (`pixi install`), not the bundle's own: the bundle would
solve a second multi-gigabyte environment to do identical work.
