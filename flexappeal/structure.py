"""Light structure introspection for the Prepare tab.

Parses a PDB or mmCIF file with gemmi and reports what the form needs to offer
intelligent choices: which chains exist, what heteroatoms are present, where the
chain breaks are, which cysteines are bonded, and roughly how big the solvated
system will be.

Deliberately light. This runs on a 3.8 GB droplet shared with four other
applications, so it does no repair, no protonation and no parameterisation --
gemmi reads the file, we count things, and the answer is a few kilobytes of JSON.
All the real preparation work happens inside the bundle, on the user's machine,
where PDBFixer and OpenMM actually live.

(Named `structure` rather than `inspect` on purpose: a module called `inspect`
inside a package sits one careless relative import away from shadowing the
standard library one.)
"""

from __future__ import annotations

import math
import re
from dataclasses import asdict, dataclass, field
from typing import Any

import gemmi

# Residue names that are solvent or common crystallisation additives rather than
# anything a user means by "ligand". Offering to keep PEG fragments and cryo-
# protectant as parameterisable ligands is noise, so they are flagged rather than
# hidden -- a user who genuinely wants glycerol in their system can still say so.
_CRYSTALLISATION_ADDITIVES = {
    "GOL", "EDO", "PEG", "PG4", "PGE", "1PE", "MPD", "DMS", "TRS", "MES",
    "EPE", "ACT", "ACY", "FMT", "SO4", "PO4", "CIT", "TLA", "IMD", "BME",
}

_COMMON_IONS = {
    "NA", "K", "MG", "CA", "ZN", "MN", "FE", "FE2", "CU", "CU1", "CO", "NI",
    "CD", "HG", "CL", "BR", "IOD", "F", "CS", "RB", "LI", "SR", "BA",
}

# Cofactors worth calling out by name, because keeping one has real consequences
# for parameterisation and a user should be told rather than left to notice.
_KNOWN_COFACTORS = {
    "HEM": "haem b", "HEC": "haem c", "NAD": "NAD⁺", "NAI": "NADH",
    "NAP": "NADP⁺", "NDP": "NADPH", "FAD": "FAD", "FMN": "FMN",
    "ATP": "ATP", "ADP": "ADP", "AMP": "AMP", "GTP": "GTP", "GDP": "GDP",
    "SAM": "S-adenosylmethionine", "SAH": "S-adenosylhomocysteine",
    "COA": "coenzyme A", "PLP": "pyridoxal phosphate", "BTN": "biotin",
    "TPP": "thiamine pyrophosphate", "MTX": "methotrexate",
}

# Van der Waals volume per protein atom, in nm3. Used only for estimating how
# much water a box will need; the value is fitted to the observed ratio of
# protein volume to atom count across typical globular proteins, not derived.
_NM3_PER_SOLUTE_ATOM = 0.0145

# Water molecules per nm3 at 1 g/mL and 310 K.
_WATERS_PER_NM3 = 33.4

# Volume of the periodic cell relative to a cube enclosing the same minimum-image
# distance. These are the standard geometric factors -- the reason a dodecahedron
# is worth choosing.
_BOX_SHAPE_FACTOR = {"cube": 1.0, "dodecahedron": 0.7071, "octahedron": 0.7698}

_ATOMS_PER_WATER = {
    "tip3p": 3, "tip3pfb": 3, "spce": 3, "opc3": 3, "charmm_tip3p": 3,
    "tip4pew": 4, "tip4pfb": 4, "opc": 4,
}


class StructureError(ValueError):
    """The uploaded file could not be understood as a structure."""


def _is_polymer_residue(residue: gemmi.Residue) -> bool:
    """Whether a residue is part of the protein/nucleic chain rather than a hetero group.

    gemmi puts this on the tabulated residue info rather than on Residue itself,
    so it needs the lookup -- an unknown residue name returns a record whose
    `found` is False and which answers False to everything, which is the correct
    behaviour for a modified or non-standard residue we cannot classify.
    """
    info = gemmi.find_tabulated_residue(residue.name)
    return bool(info and (info.is_amino_acid() or info.is_nucleic_acid()))


@dataclass
class ChainInfo:
    id: str
    kind: str  # protein | nucleic | water | other
    residues: int
    observed_residues: int
    sequence: str = ""
    gaps: list[dict[str, Any]] = field(default_factory=list)
    atoms: int = 0

    @property
    def missing_residues(self) -> int:
        return max(0, self.residues - self.observed_residues)


# Elements no fixed-charge small-molecule force field can parameterise. A ligand
# containing one of these cannot go through GAFF, OpenFF or espaloma at all --
# they are organic force fields, and a transition metal centre needs either
# bespoke bonded parameters or a QM treatment. Detected so the app can refuse
# clearly at build time rather than failing inside the user's run.
_UNPARAMETERISABLE_METALS = {
    "FE", "CU", "ZN", "MN", "CO", "NI", "MO", "W", "V", "CR", "CD", "HG",
    "PT", "PD", "RU", "RH", "IR", "AU", "AG", "TI", "SN", "PB", "AS", "SB",
}


@dataclass
class HeteroInfo:
    name: str
    count: int
    chains: list[str]
    atoms: int
    category: str  # cofactor | ligand | ion | water | additive
    description: str = ""
    formula: str = ""
    elements: list[str] = field(default_factory=list)

    @property
    def metals(self) -> list[str]:
        return sorted({e for e in self.elements if e.upper() in _UNPARAMETERISABLE_METALS})

    @property
    def parameterisable(self) -> bool:
        """Whether a small-molecule force field could handle this at all."""
        return self.category in ("ligand", "cofactor", "additive") and not self.metals


@dataclass
class DisulfideInfo:
    chain_a: str
    resid_a: int
    chain_b: str
    resid_b: int
    distance: float


@dataclass
class StructureReport:
    """Everything the Prepare form learns from the uploaded file."""

    name: str
    format: str
    models: int
    chains: list[ChainInfo]
    heteroatoms: list[HeteroInfo]
    disulfides: list[DisulfideInfo]
    has_hydrogens: bool
    has_altlocs: bool
    altlocs: list[str]
    resolution: float | None
    method: str
    title: str
    warnings: list[str]

    # Geometry, used for the size and time estimates on the form.
    extent_nm: tuple[float, float, float]
    max_extent_nm: float
    solute_atoms: int

    def to_dict(self) -> dict[str, Any]:
        d = asdict(self)
        d["extent_nm"] = list(self.extent_nm)
        return d

    @property
    def protein_chains(self) -> list[ChainInfo]:
        return [c for c in self.chains if c.kind == "protein"]

    @property
    def total_residues(self) -> int:
        return sum(c.observed_residues for c in self.protein_chains)


# ===========================================================================
#  Parsing
# ===========================================================================


def _detect_format(data: bytes, filename: str) -> str:
    lowered = filename.lower()
    if lowered.endswith((".cif", ".mmcif", ".bcif")):
        return "mmcif"
    if lowered.endswith((".pdb", ".ent")):
        return "pdb"
    head = data[:4096].decode("utf-8", errors="replace")
    if "data_" in head or "_atom_site." in head:
        return "mmcif"
    return "pdb"


def read_structure(data: bytes, filename: str = "structure.pdb") -> gemmi.Structure:
    """Parse bytes into a gemmi Structure, with errors a user can act on."""
    fmt = _detect_format(data, filename)
    text = data.decode("utf-8", errors="replace")
    try:
        if fmt == "mmcif":
            doc = gemmi.cif.read_string(text)
            st = gemmi.make_structure_from_block(doc.sole_block())
        else:
            st = gemmi.read_pdb_string(text)
    except Exception as exc:  # noqa: BLE001 -- gemmi raises several unrelated types
        raise StructureError(
            f"this does not parse as a {fmt.upper()} file: {exc}. Check that the "
            f"upload is a coordinate file rather than a header, a map or an archive."
        ) from None

    if not len(st):
        raise StructureError("the file contains no models.")
    st.setup_entities()
    return st


def _merge_chains_by_id(chains: list[ChainInfo]) -> list[ChainInfo]:
    """Collapse the per-entity fragments gemmi produces back into author chains.

    Kind resolution is ordered rather than last-wins: a chain that has any
    polymer content is a polymer chain, whatever its hetero fragment says.
    """
    order = {"protein": 0, "nucleic": 1, "other": 2, "water": 3}
    merged: dict[str, ChainInfo] = {}
    for chain in chains:
        existing = merged.get(chain.id)
        if existing is None:
            merged[chain.id] = chain
            continue
        if order[chain.kind] < order[existing.kind]:
            existing.kind = chain.kind
        existing.atoms += chain.atoms
        existing.residues += chain.residues
        existing.observed_residues += chain.observed_residues
        existing.gaps.extend(chain.gaps)
        if len(chain.sequence) > len(existing.sequence):
            existing.sequence = chain.sequence
    return list(merged.values())


def _classify_chain(chain: gemmi.Chain) -> str:
    polymer = chain.get_polymer()
    if len(polymer):
        kind = polymer.check_polymer_type()
        if kind in (gemmi.PolymerType.PeptideL, gemmi.PolymerType.PeptideD):
            return "protein"
        if kind in (gemmi.PolymerType.Dna, gemmi.PolymerType.Rna, gemmi.PolymerType.DnaRnaHybrid):
            return "nucleic"
    names = {r.name for r in chain}
    if names and names <= {"HOH", "WAT", "DOD"}:
        return "water"
    return "other"


def _categorise_hetero(name: str) -> tuple[str, str]:
    upper = name.upper()
    if upper in ("HOH", "WAT", "DOD"):
        return "water", "water"
    if upper in _COMMON_IONS:
        return "ion", "ion"
    if upper in _KNOWN_COFACTORS:
        return "cofactor", _KNOWN_COFACTORS[upper]
    if upper in _CRYSTALLISATION_ADDITIVES:
        return "additive", "crystallisation additive"
    return "ligand", ""


def _chain_gaps(chain: gemmi.Chain, entity: gemmi.Entity | None) -> list[dict[str, Any]]:
    """Find breaks in the observed sequence.

    Two independent signals, because neither alone is reliable: a jump in
    residue numbering (catches most cases, misses gaps where numbering was
    renumbered contiguously) and SEQRES length versus observed count (catches
    the rest, but says nothing about where).
    """
    gaps: list[dict[str, Any]] = []
    residues = [r for r in chain if _is_polymer_residue(r)]
    for prev, curr in zip(residues, residues[1:]):
        delta = curr.seqid.num - prev.seqid.num
        if delta > 1:
            gaps.append({
                "after": prev.seqid.num,
                "before": curr.seqid.num,
                "length": delta - 1,
                "source": "numbering",
            })
    if entity is not None and entity.full_sequence:
        declared = len(entity.full_sequence)
        observed = len(residues)
        accounted = sum(g["length"] for g in gaps)
        unexplained = declared - observed - accounted
        if unexplained > 0:
            gaps.append({
                "after": None,
                "before": None,
                "length": unexplained,
                "source": "seqres",
            })
    return gaps


def _find_disulfides(st: gemmi.Structure, cutoff: float = 2.5) -> list[DisulfideInfo]:
    """SG-SG pairs closer than the cutoff, in ångström."""
    sg_atoms: list[tuple[str, int, gemmi.Position]] = []
    for chain in st[0]:
        for residue in chain:
            if residue.name != "CYS":
                continue
            for atom in residue:
                if atom.name == "SG" and atom.altloc in ("\0", "", "A"):
                    sg_atoms.append((chain.name, residue.seqid.num, atom.pos))

    found: list[DisulfideInfo] = []
    for i, (chain_a, resid_a, pos_a) in enumerate(sg_atoms):
        for chain_b, resid_b, pos_b in sg_atoms[i + 1:]:
            distance = pos_a.dist(pos_b)
            if distance <= cutoff:
                found.append(DisulfideInfo(chain_a, resid_a, chain_b, resid_b, round(distance, 2)))
    return found


def analyse(data: bytes, filename: str = "structure.pdb",
            disulfide_cutoff: float = 2.5) -> StructureReport:
    """Produce the full report the Prepare form renders from."""
    st = read_structure(data, filename)
    model = st[0]
    warnings: list[str] = []
    # st.info is a gemmi map rather than a dict; copy it so .get() is available
    # and the values are plain strings by the time they reach the template.
    info = {k: str(v) for k, v in dict(st.info).items()}

    # Entities are keyed by subchain id ('Axp', 'Axw'), not by chain name, so the
    # SEQRES length for a chain has to be reached through its subchains.
    entity_by_subchain = {sub: e for e in st.entities for sub in e.subchains}

    chains: list[ChainInfo] = []
    hetero_counts: dict[str, dict[str, Any]] = {}
    altlocs: set[str] = set()
    has_hydrogens = False
    solute_atoms = 0

    positions: list[gemmi.Position] = []

    for chain in model:
        kind = _classify_chain(chain)
        polymer = chain.get_polymer()
        entity = None
        for sub in chain.subchains():
            entity = entity_by_subchain.get(sub.subchain_id())
            if entity is not None and entity.full_sequence:
                break
        declared = len(entity.full_sequence) if entity and entity.full_sequence else 0

        observed = 0
        atom_count = 0
        for residue in chain:
            is_polymer_residue = _is_polymer_residue(residue)
            if is_polymer_residue:
                observed += 1
            for atom in residue:
                atom_count += 1
                if atom.element == gemmi.Element("H"):
                    has_hydrogens = True
                if atom.altloc not in ("\0", ""):
                    altlocs.add(atom.altloc)
                if not residue.is_water():
                    solute_atoms += 1
                    positions.append(atom.pos)

            if not is_polymer_residue:
                category, description = _categorise_hetero(residue.name)
                entry = hetero_counts.setdefault(residue.name, {
                    "count": 0, "chains": set(), "atoms": 0,
                    "category": category, "description": description,
                    "elements": set(),
                })
                entry["count"] += 1
                entry["chains"].add(chain.name)
                entry["atoms"] += len(residue)
                entry["elements"].update(a.element.name.upper() for a in residue)

        sequence = polymer.make_one_letter_sequence() if len(polymer) else ""

        chains.append(ChainInfo(
            id=chain.name,
            kind=kind,
            residues=declared or observed,
            observed_residues=observed,
            sequence=sequence,
            gaps=_chain_gaps(chain, entity),
            atoms=atom_count,
        ))

    # setup_entities() splits one author chain into separate gemmi Chain objects
    # per entity -- so chain A of a haem protein arrives as "A (polymer)" plus
    # "A (non-polymer, 0 residues)". A user thinks of that as one chain, and
    # showing it twice in the chain picker is both confusing and would let them
    # select the same chain inconsistently.
    chains = _merge_chains_by_id(chains)

    heteroatoms = [
        HeteroInfo(
            name=name,
            count=info["count"],
            chains=sorted(info["chains"]),
            atoms=info["atoms"],
            category=info["category"],
            description=info["description"],
            elements=sorted(info["elements"]),
        )
        for name, info in sorted(hetero_counts.items())
    ]

    disulfides = _find_disulfides(st, disulfide_cutoff)

    # --- geometry, in nm ---
    if positions:
        xs = [p.x for p in positions]
        ys = [p.y for p in positions]
        zs = [p.z for p in positions]
        extent = (
            (max(xs) - min(xs)) / 10.0,
            (max(ys) - min(ys)) / 10.0,
            (max(zs) - min(zs)) / 10.0,
        )
    else:
        extent = (0.0, 0.0, 0.0)

    # --- warnings the user should see before they commit to a run ---
    if not chains:
        raise StructureError("the file contains no chains.")
    if not any(c.kind == "protein" for c in chains):
        warnings.append(
            "no protein chain was found. FlexAppeal is built for protein molecular "
            "dynamics; nucleic-acid and ligand-only systems are not supported."
        )
    for c in chains:
        missing = c.missing_residues
        if missing > 30:
            warnings.append(
                f"chain {c.id} is missing {missing} residues. Rebuilt loops that long "
                f"are guesses and will dominate the RMSD; consider simulating a "
                f"single well-resolved domain instead."
            )
        elif missing > 0:
            warnings.append(
                f"chain {c.id} is missing {missing} residue(s), which PDBFixer will rebuild."
            )
    if len(st) > 1:
        warnings.append(
            f"the file contains {len(st)} models. Model 1 is used unless you choose "
            f"another in the advanced options."
        )
    if altlocs:
        warnings.append(
            f"alternate conformations are present ({', '.join(sorted(altlocs))}). "
            f"Only one can be simulated; the selection rule is in the advanced options."
        )
    ligands = [h for h in heteroatoms if h.category in ("ligand", "cofactor")]
    if ligands:
        warnings.append(
            f"{len(ligands)} ligand/cofactor type(s) found ({', '.join(h.name for h in ligands)}). "
            f"Each one you keep needs small-molecule parameters, which adds a few "
            f"minutes to the bundle's first run."
        )
    metal_ligands = [h for h in heteroatoms
                     if h.metals and h.category in ("ligand", "cofactor")]
    if metal_ligands:
        warnings.append(
            "metal-containing cofactor(s) found ("
            + ", ".join(f"{h.name}: {', '.join(h.metals)}" for h in metal_ligands)
            + "). No fixed-charge small-molecule force field can parameterise a "
            "transition-metal centre, so these cannot be kept. Simulating the apo "
            "protein is the supported option; a metalloprotein needs bespoke bonded "
            "parameters or a QM/MM treatment."
        )

    metal_ions = [h for h in heteroatoms if h.metals and h.category == "ion"]
    if metal_ions:
        warnings.append(
            "structural metal ion(s) present ("
            + ", ".join(h.name for h in metal_ions)
            + "). Standard force fields treat these as non-bonded point charges, "
            "which handles a loosely bound ion adequately and a catalytic metal "
            "centre badly."
        )

    return StructureReport(
        name=st.name or filename,
        format=_detect_format(data, filename),
        models=len(st),
        chains=chains,
        heteroatoms=heteroatoms,
        disulfides=disulfides,
        has_hydrogens=has_hydrogens,
        has_altlocs=bool(altlocs),
        altlocs=sorted(altlocs),
        resolution=st.resolution or None,
        method=info.get("_exptl.method", ""),
        title=info.get("_struct.title", ""),
        warnings=warnings,
        extent_nm=extent,
        max_extent_nm=max(extent) if extent else 0.0,
        solute_atoms=solute_atoms,
    )


# ===========================================================================
#  System size estimation
# ===========================================================================


def estimate_system_size(report: StructureReport, cfg: dict[str, Any]) -> dict[str, Any]:
    """Estimate the solvated atom count for the live readouts on the form.

    Approximate by construction, and labelled as such in the UI. The exact count
    is not knowable until Modeller.addSolvent has actually run, which happens on
    the user's machine -- but an estimate that is right to within about 10% is
    what makes the size and wall-time readouts useful while someone is still
    choosing options.
    """
    solute = report.solute_atoms
    mode = cfg.get("solvent_mode", "explicit")

    if mode in ("implicit", "vacuum") or not solute:
        return {
            "solute_atoms": solute,
            "water_molecules": 0,
            "ion_count": 0,
            "total_atoms": solute,
            "box_nm": None,
            "basis": "no explicit solvent",
        }

    if cfg.get("use_membrane"):
        return _estimate_membrane_size(report, cfg)

    padding = float(cfg.get("padding") or 1.2)
    shape = str(cfg.get("box_shape") or "dodecahedron")
    water_model = str(cfg.get("water_model") or "tip3p")

    edge = report.max_extent_nm + 2.0 * padding
    cube_volume = edge ** 3
    volume = cube_volume * _BOX_SHAPE_FACTOR.get(shape, 1.0)

    solute_volume = solute * _NM3_PER_SOLUTE_ATOM
    water_volume = max(volume - solute_volume, 0.0)
    waters = int(water_volume * _WATERS_PER_NM3)

    # Ions needed to reach the target ionic strength, plus whatever neutralises
    # the net charge. The charge is not known here (it depends on protonation,
    # which happens in the bundle), so this covers the ionic strength only.
    ionic_strength = float(cfg.get("ionic_strength") or 0.0)
    ions = int(round(2 * ionic_strength * water_volume * 0.6022))

    atoms_per_water = _ATOMS_PER_WATER.get(water_model, 3)
    total = solute + waters * atoms_per_water + ions

    return {
        "solute_atoms": solute,
        "water_molecules": waters,
        "ion_count": ions,
        "total_atoms": total,
        "box_nm": round(edge, 2),
        "box_volume_nm3": round(volume, 1),
        "basis": f"{shape} box, {padding:g} nm padding",
    }


def _estimate_membrane_size(report: StructureReport, cfg: dict[str, Any]) -> dict[str, Any]:
    """Membrane systems are dominated by lipid and by the water slabs above and below."""
    solute = report.solute_atoms
    padding = float(cfg.get("membrane_padding") or 1.0)
    water_model = str(cfg.get("water_model") or "charmm_tip3p")

    # addMembrane sizes the patch to the solute's xy footprint plus padding.
    xy = max(report.extent_nm[0], report.extent_nm[1]) + 2.0 * padding
    area = xy ** 2

    # ~0.65 nm2 per lipid per leaflet for POPC, ~130 atoms per lipid.
    lipids = int((area / 0.65) * 2)
    lipid_atoms = lipids * 130

    # Bilayer is ~4 nm thick; water slabs fill the rest of z.
    z = report.extent_nm[2] + 2.0 * padding
    water_thickness = max(z - 4.0, 2.0 * padding)
    water_volume = area * water_thickness
    waters = int(water_volume * _WATERS_PER_NM3)

    ionic_strength = float(cfg.get("ionic_strength") or 0.15)
    ions = int(round(2 * ionic_strength * water_volume * 0.6022))

    atoms_per_water = _ATOMS_PER_WATER.get(water_model, 3)
    total = solute + lipid_atoms + waters * atoms_per_water + ions

    return {
        "solute_atoms": solute,
        "lipid_count": lipids,
        "lipid_atoms": lipid_atoms,
        "water_molecules": waters,
        "ion_count": ions,
        "total_atoms": total,
        "box_nm": round(xy, 2),
        "basis": f"{cfg.get('lipid_type', 'POPC')} bilayer, {xy:.1f} × {xy:.1f} nm patch",
    }


def dynamic_choices(report: StructureReport) -> dict[str, list[dict[str, str]]]:
    """Fill in the choices the registry cannot know until a structure is loaded."""
    chain_choices = [
        {
            "value": c.id,
            "label": f"Chain {c.id}",
            "help": (
                f"{c.kind}, {c.observed_residues} residues"
                + (f", {c.missing_residues} missing" if c.missing_residues else "")
            ),
        }
        for c in report.chains
        if c.kind != "water"
    ]

    hetero_choices = [
        {
            "value": h.name,
            "label": f"{h.name}{' × ' + str(h.count) if h.count > 1 else ''}",
            "help": (
                (h.description or h.category)
                + f", {h.atoms // max(h.count, 1)} atoms"
                + (
                    " -- needs small-molecule parameters"
                    if h.category in ("ligand", "cofactor") else ""
                )
            ),
        }
        for h in report.heteroatoms
        if h.category != "water"
    ]

    return {"chains": chain_choices, "keep_heteroatoms": hetero_choices}
