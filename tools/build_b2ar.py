#!/usr/bin/env python3
"""Build the two β2AR start structures from GPCRdb refined models.

Spec section 4.3. The refined models are the top rung of its fallback ladder:
fusions removed, thermostabilising mutations reverted, missing segments modelled.
Not every entry has one, and the spec says so: 2RH1 and 3SN6 both report none,
which this build confirmed. Discovered at run time rather than assumed.

    active, apo       4LDE refined  2.79 A, the best-resolution active entry
                                    that has a refined model, stripped of its
                                    agonist
    inactive, ligand  5D5A refined  2.48 A, and it carries (S)-carazolol, an
                                    inverse agonist, already in place

The second of those was a piece of luck worth recording. The plan in D007 was to
transplant carazolol onto a refined model by superposing 2RH1, because refined
models were assumed to carry no ligand. They do carry their ligands, and 5D5A's
is the exact inverse agonist wanted, so the transplant and its clash risk are not
needed at all.

The download is a zip, not a PDB, and the URL pattern is
`/structure/refined/<PDB>_full/download_pdb`. `_noloops` exists beside it and
drops ICL3 entirely; `_full` is taken because section 4.3 says to keep whatever
the refined model contains.
"""

from __future__ import annotations

import argparse
import io
import json
import sys
import zipfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

REFINED_URL = "https://gpcrdb.org/structure/refined/{pdb}_full/download_pdb"

STARTS = {
    "apo":       {"pdb": "4LDE", "state": "active",   "keep_ligand": None},
    "carazolol": {"pdb": "5D5A", "state": "inactive", "keep_ligand": "CAU"},
}

# Everything in these refined models that is not the receptor or its ligand:
# crystallisation additives, lipids and ions that belong to the experiment rather
# than to the system being simulated. The bilayer is built fresh in section 5.1,
# so a few crystallographic cholesterols would be in the way rather than useful.
ADDITIVES = {"BU1", "ACM", "CLR", "PLM", "12P", "SO4", "NA", "CL", "HOH",
             "OLA", "OLB", "EDO", "PEG", "GOL", "MES", "TRS"}


OPM_URL = "https://opm-assets.storage.googleapis.com/pdb/{pdb}.pdb"
# The canonical OPM-oriented beta2AR. Its own refined model does not exist and
# OPM has no entry for 5D5A or 4LDE, so the orientation is borrowed from the
# receptor rather than from the entry, which is what section 5.1 allows: the
# bilayer frame belongs to the protein fold, not to one crystal form.
OPM_REFERENCE = "2rh1"


def orient_to_opm(start_pdb: Path, numbering) -> dict:
    """Rotate the start model into OPM's membrane frame, in place.

    `addMembrane` needs the bilayer normal along z and the membrane centred at
    z = 0, and a refined model carries whatever frame its deposition had. These
    two arrived 50 A apart and tilted differently from each other.

    Superposed on TM Calpha only, matched by residue number. Loops are excluded
    because ICL3 is fully modelled here and is the most mobile part of the
    receptor: fitting on it would tilt the helices to satisfy a loop.
    """
    import gemmi
    import numpy as np
    import requests

    response = requests.get(OPM_URL.format(pdb=OPM_REFERENCE), timeout=90)
    response.raise_for_status()
    reference = gemmi.read_pdb_string(response.text)
    reference.setup_entities()

    # OPM marks the bilayer with dummy atoms at the two headgroup planes.
    half_thickness = max(
        abs(atom.pos.z)
        for chain in reference[0] for residue in chain if residue.name == "DUM"
        for atom in residue)

    tm_residues = {n for n, seg in numbering.segment.items() if seg.startswith("TM")}

    def tm_ca(structure):
        found = {}
        for chain in structure[0]:
            for residue in chain:
                if residue.seqid.num not in tm_residues:
                    continue
                atom = residue.find_atom("CA", "*")
                if atom is not None:
                    found.setdefault(residue.seqid.num, atom.pos)
        return found

    mobile = gemmi.read_structure(str(start_pdb))
    mobile.setup_entities()
    target_ca, mobile_ca = tm_ca(reference), tm_ca(mobile)
    shared = sorted(set(target_ca) & set(mobile_ca))
    if len(shared) < 50:
        raise SystemExit(f"only {len(shared)} shared TM Ca; cannot orient reliably")

    target = gemmi.Position_list() if hasattr(gemmi, "Position_list") else []
    sup = gemmi.superpose_positions([target_ca[n] for n in shared],
                                    [mobile_ca[n] for n in shared])
    for model in mobile:
        for chain in model:
            for residue in chain:
                for atom in residue:
                    atom.pos = gemmi.Position(sup.transform.apply(atom.pos))
    mobile.write_pdb(str(start_pdb))

    # Section 5.1's validation, on the oriented coordinates.
    z = {}
    for chain in mobile[0]:
        for residue in chain:
            atom = residue.find_atom("CA", "*")
            if atom is not None:
                z[residue.seqid.num] = atom.pos.z
    belt = [abs(z[n]) for n, letter in numbering.residue_letter.items()
            if letter in ("W", "Y") and n in z
            and numbering.segment.get(n, "").startswith("TM")]
    basics = [z[n] for n, letter in numbering.residue_letter.items()
              if letter in ("R", "K") and n in z
              and numbering.segment.get(n, "") in ("ICL1", "ICL2", "ICL3", "H8")]
    belt_median = float(np.median(belt)) if belt else 0.0
    basic_median = float(np.median(basics)) if basics else 0.0
    return {
        "rmsd": sup.rmsd, "n_atoms": len(shared),
        "half_thickness": round(half_thickness, 1),
        "belt_median": belt_median, "basic_median": basic_median,
        "passes": 10.0 <= belt_median <= 17.0 and basic_median < -13.0,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--which", choices=sorted(STARTS), required=True)
    parser.add_argument("--out", type=Path, default=ROOT / "data" / "systems")
    args = parser.parse_args()

    import gemmi
    import requests

    from flexappeal.refstates import Catalogue, receptor_residues
    from flexappeal.refstates.gpcr import build_numbering

    numbering = build_numbering(receptor_residues(Catalogue(ROOT / "data" / "cache" / "gpcrdb")))
    lo = min(numbering.segment_residues("TM1")) - 3
    hi = max(numbering.segment_residues("H8")) + 3
    icl3 = set(numbering.segment_residues("ICL3"))
    print(f"  construct {lo}-{hi} (TM1 start minus 3 to H8 end plus 3)")

    spec = STARTS[args.which]
    response = requests.get(REFINED_URL.format(pdb=spec["pdb"]), timeout=90)
    response.raise_for_status()
    archive = zipfile.ZipFile(io.BytesIO(response.content))
    inner = [n for n in archive.namelist() if n.lower().endswith(".pdb")][0]
    print(f"  {spec['pdb']} refined: {inner}")

    out_dir = Path(args.out)
    out_dir.mkdir(parents=True, exist_ok=True)
    raw = out_dir / f"b2ar_{args.which}_refined.pdb"
    raw.write_bytes(archive.read(inner))

    structure = gemmi.read_structure(str(raw))
    structure.setup_entities()
    structure.remove_alternative_conformations()
    structure.remove_hydrogens()

    out = gemmi.Structure()
    out.spacegroup_hm = "P 1"
    model = gemmi.Model("1")
    protein = gemmi.Chain("A")
    het = gemmi.Chain("B")
    kept_icl3 = 0
    for chain in structure[0]:
        for residue in chain:
            info = gemmi.find_tabulated_residue(residue.name)
            if info and info.is_amino_acid():
                if lo <= residue.seqid.num <= hi:
                    protein.add_residue(residue)
                    if residue.seqid.num in icl3:
                        kept_icl3 += 1
            elif spec["keep_ligand"] and residue.name == spec["keep_ligand"]:
                het.add_residue(residue)
            elif residue.name not in ADDITIVES:
                print(f"    note: dropping unrecognised heteroatom {residue.name}")
    model.add_chain(protein)
    if len(het):
        model.add_chain(het)
    out.add_model(model)
    start_pdb = out_dir / f"b2ar_{args.which}_start.pdb"
    out.write_pdb(str(start_pdb))

    orientation = orient_to_opm(start_pdb, numbering)
    print(f"  oriented onto OPM 2RH1: rmsd {orientation['rmsd']:.2f} A on "
          f"{orientation['n_atoms']} TM Ca, bilayer at z = +/-{orientation['half_thickness']}")
    print(f"    aromatic belt |z| median {orientation['belt_median']:.1f} A "
          f"(want 10 to 17), intracellular basics median z "
          f"{orientation['basic_median']:.1f} A (want below -13)  "
          f"{'PASS' if orientation['passes'] else 'FAIL'}")

    # Section 4.3: more than 20 modelled ICL3 residues means keep them and flag
    # the loop, rather than quietly letting it dominate every RMSF and RMSD.
    icl3_flag = "modelled, flexible" if kept_icl3 > 20 else "as deposited"
    meta = {
        "which": args.which, "pdb": spec["pdb"], "state": spec["state"],
        "source": "GPCRdb refined model (ladder rung 1)",
        "url": REFINED_URL.format(pdb=spec["pdb"]),
        "ligand": spec["keep_ligand"], "span": [lo, hi],
        "n_residues": len(protein), "n_ligand": len(het),
        "icl3_modelled": kept_icl3, "icl3_flag": icl3_flag,
        "icl3_excluded_from_rmsd": kept_icl3 > 20,
    }
    (out_dir / f"b2ar_{args.which}_start.json").write_text(json.dumps(meta, indent=2))
    print(f"  wrote {start_pdb.name}: {len(protein)} residues, {len(het)} ligand, "
          f"ICL3 {kept_icl3}/24 ({icl3_flag})")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
