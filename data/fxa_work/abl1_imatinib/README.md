# abl1_imatinib

ABL1 kinase domain in its inactive DFG-out conformation with imatinib bound, in explicit water. The positive control for an inactive start staying inactive, and for a drug holding a state.

A FlexAppeal run bundle, generated 2026-10-14 02:56 UTC by FlexAppeal 0.1.0.
Prepared by Marc C. Deller.
## Running it

```bash
chmod +x flexappeal_abl1_imatinib.command
./flexappeal_abl1_imatinib.command
```

The `chmod` is required: a file downloaded from a browser arrives without the
executable bit, so double-clicking will not work until you have run it.

The first run installs pixi (if you do not already have it) and solves this
bundle's environment, which is several gigabytes and takes 5 to 15 minutes.
Later runs skip straight to the simulation.

Re-running resumes from the last checkpoint. To start over, delete
`./flexappeal_output/`.

## What it will do

| Stage | Setting |
|---|---|
| Structure | 2HYY.pdb |
| Force field | amber14-all.xml, amber14/tip3p.xml |
| Solvent | tip3p in a dodecahedron, 1.2 nm padding, 0.15 M Na+Cl- |
| Minimisation | to 10.0 kJ/mol/nm |
| Heating | 50.0 → 300.0 K over 100.0 ps, 5 stages |
| Equilibration | 1000.0 ps, restraints released 1000.0 → 500.0 → 100.0 → 10.0 → 0.0 kJ/mol/nm² |
| Production | 50.0 ns at 300.0 K, 1.0 bar |
| Integrator | LangevinMiddle, 4.0 fs, hydrogen mass repartitioned to 1.5 amu |
| Constraints | HBonds |
| Trajectory | XTC every 10.0 ps, protein_ligand — about 5000 frames, 43.9 KB |
| Platform | benchmarked at start-up |

Rough estimate: **load a structure to estimate**. The bundle benchmarks your hardware
before it starts and prints a real figure.

## What you get

Everything lands in `./flexappeal_output/`:

| File | What it is |
|---|---|
| `abl1_imatinib.fxa` | **The results file.** Upload this to the Analysis tab. |
| `prepared.pdb` | after repair, before solvation |
| `solvated.pdb` | the complete system as simulated |
| `topology.pdb` | first frame, matched to the trajectory's atom subset |
| `trajectory.xtc` | the full trajectory — stays here, never uploaded |
| `state_data.csv` | energies, temperature, pressure, density over time |
| `state.chk` | checkpoint, for resuming |
| `system.xml`, `integrator.xml`, `final_state.xml` | full serialisation, for exact reproduction or extension |
| `run_manifest.json` | versions, platform, achieved ns/day, wall time |

The `.fxa` carries the metrics computed against your **full** trajectory plus a
decimated copy for viewing (standard tier). The full trajectory
never leaves your machine.

## Changing something

`run.py` is a plain, readable OpenMM script with every value written out
literally — edit it and re-run if you want to change something small. For
anything structural, go back to https://flexappeal.mdeller.com and build a new bundle.

`config.json` holds the same settings in machine-readable form and is copied
into the results manifest, so a `.fxa` always knows exactly how it was made.

## If something goes wrong

| Symptom | Likely cause |
|---|---|
| `pixi: command not found` after install | Run `source ~/.zshrc`, then re-run the bundle. |
| Environment fails to solve | Network, usually. Try again; the solve is cached. |
| `No template found for residue X` | A heteroatom was kept that has no parameters. Rebuild the bundle without it, or supply parameters. |
| Energy rises steadily through the run | The timestep is too long for the constraints in use. |
| Everything is much slower than the benchmark | Something else is competing for CPU, or the machine is swapping. |
| `Particle coordinate is NaN` | The system blew up: almost always a clash minimisation did not resolve, or too long a timestep. |

---

Built by [Marc C. Deller, D.Phil.](https://marcdeller.com) · [marc@marcdeller.com](mailto:marc@marcdeller.com)
