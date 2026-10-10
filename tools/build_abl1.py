#!/usr/bin/env python3
"""Build the two ABL1 systems: apo from 2GQG and imatinib-bound from 2HYY.

Spec section 5.1 for the protocol and 4.2 for the construct. The construct rule
is the part worth reading twice: both systems are trimmed to the residues present
in *both* start structures plus up to five flanking, so the two trajectories are
of the same molecule and a motif RMSD between them means something. Trimming each
to its own resolved span would give two different proteins and a comparison that
quietly is not one.

No residue number is written down here. The span comes from what the two
structures actually resolve, and the ligand comes from the RCSB chemical
definition by its three-letter code.
"""

from __future__ import annotations

import argparse
import gzip
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "vendor" / "boltzmaker"))

# Section 4.2: the seeds, and the chains KLIFS lists for them.
SEEDS = {
    "apo":     {"pdb": "2GQG", "chain": "A", "ligand": None,  "state": "DFG-in"},
    "imatinib": {"pdb": "2HYY", "chain": "C", "ligand": "STI", "state": "DFG-out"},
}
FLANK = 5


def _is_protein(residue) -> bool:
    """gemmi puts the amino-acid test on the tabulated entry, not the residue."""
    import gemmi

    return gemmi.find_tabulated_residue(residue.name) is not None and \
        gemmi.find_tabulated_residue(residue.name).is_amino_acid()


def resolved_span(cif_gz: Path, chain_id: str) -> tuple[int, int, set[int]]:
    """The first and last modelled residue of one chain, and every number in it."""
    import gemmi

    text = gzip.open(cif_gz, "rt").read()
    structure = gemmi.read_structure_string(text, format=gemmi.CoorFormat.Mmcif)
    structure.setup_entities()
    structure.remove_alternative_conformations()
    structure.remove_hydrogens()
    for chain in structure[0]:
        if chain.name != chain_id:
            continue
        numbers = {r.seqid.num for r in chain if _is_protein(r)}
        if numbers:
            return min(numbers), max(numbers), numbers
    raise SystemExit(f"{cif_gz.name}: chain {chain_id} has no amino acids")


def write_trimmed(cif_gz: Path, chain_id: str, keep: set[int], ligand: str | None,
                  out_pdb: Path) -> dict:
    """One chain, trimmed to `keep`, waters and other chains gone.

    The ligand is kept only when asked for, and then only the copy in this chain:
    a crystal form with four ligand copies would otherwise put all four into a
    system that has one binding site.
    """
    import gemmi

    text = gzip.open(cif_gz, "rt").read()
    structure = gemmi.read_structure_string(text, format=gemmi.CoorFormat.Mmcif)
    structure.setup_entities()
    structure.remove_alternative_conformations()
    structure.remove_hydrogens()
    structure.remove_waters()

    out = gemmi.Structure()
    out.spacegroup_hm = "P 1"
    model = gemmi.Model("1")
    protein = gemmi.Chain("A")
    het = gemmi.Chain("B")
    n_lig = 0
    for chain in structure[0]:
        if chain.name != chain_id:
            continue
        for residue in chain:
            if _is_protein(residue) and residue.seqid.num in keep:
                protein.add_residue(residue)
            elif ligand and residue.name == ligand and n_lig == 0:
                het.add_residue(residue)
                n_lig += 1
    model.add_chain(protein)
    if len(het):
        model.add_chain(het)
    out.add_model(model)
    out_pdb.parent.mkdir(parents=True, exist_ok=True)
    out.write_pdb(str(out_pdb))
    return {"residues": len(protein), "ligand_copies": n_lig}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--which", choices=sorted(SEEDS), required=True)
    parser.add_argument("--out", type=Path, default=ROOT / "data" / "systems")
    args = parser.parse_args()

    from flexappeal.refstates.refs import fetch_mmcif

    cif_dir = ROOT / "data" / "refs" / "cif"
    spans = {}
    for name, seed in SEEDS.items():
        path = fetch_mmcif(seed["pdb"], cif_dir)
        lo, hi, numbers = resolved_span(path, seed["chain"])
        spans[name] = (lo, hi, numbers)
        print(f"  {seed['pdb']}/{seed['chain']}: residues {lo}-{hi}, {len(numbers)} modelled")

    # The shared construct: present in both, plus flanking where both have it.
    common = spans["apo"][2] & spans["imatinib"][2]
    lo, hi = min(common), max(common)
    keep = {n for n in common if lo <= n <= hi}
    print(f"  shared span {lo}-{hi}: {len(keep)} residues in both start structures")

    seed = SEEDS[args.which]
    out_pdb = Path(args.out) / f"abl1_{args.which}_start.pdb"
    info = write_trimmed(fetch_mmcif(seed["pdb"], cif_dir), seed["chain"], keep,
                         seed["ligand"], out_pdb)
    print(f"  wrote {out_pdb.name}: {info['residues']} residues, "
          f"{info['ligand_copies']} ligand copies")
    (Path(args.out) / f"abl1_{args.which}_start.json").write_text(json.dumps({
        "which": args.which, "pdb": seed["pdb"], "chain": seed["chain"],
        "ligand": seed["ligand"], "expected_state": seed["state"],
        "span": [lo, hi], "n_residues": info["residues"],
        "shared_construct": True,
    }, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
