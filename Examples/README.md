# Examples

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
| `run/` | The bundle extracted: `run.py`, `analyse.py`, `pixi.toml`, `input.pdb` |
| `terminal.log.gz` | Raw `script(1)` capture of the run, escape sequences and all (3 MB → 80 KB) |
| `analyse.log.gz` | The same for the analysis pass |
| `screenshots/*.html` | Terminal captures rendered for the web page, from that log |
| `output/*.fxa` | The results file, which is what the Example tab renders |
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
bash Examples/lysozyme_10ns/rebuild.sh
```

Roughly ninety minutes on an M1 Max, most of it the simulation. It needs the
development environment (`pixi install`), not the bundle's own: the bundle would
solve a second multi-gigabyte environment to do identical work.
