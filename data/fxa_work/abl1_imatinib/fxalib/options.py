"""The FlexAppeal option registry -- the single source of truth for every
OpenMM parameter the app exposes.

Everything else reads from here and nothing else declares an option:

  * ``flexappeal.webapp``      renders the Prepare form from ``GROUPS``/``OPTIONS``
  * ``flexappeal.schema``      coerces and validates a POST against these declarations
  * ``flexappeal.runtime``     templates the generated ``run.py`` from the validated config
  * ``docs/options.md``        is generated from here (``FlexAppeal.py docs``)

Adding an OpenMM option must mean editing exactly one file. That constraint is
the whole reason this module exists -- the alternative (a hand-written HTML form,
a hand-written validator and a hand-written script template) drifts within a week.

Design notes
------------
``requires`` is a predicate string over other option ids, evaluated by
``schema.evaluate_predicate`` against the current config. It drives three things
at once: which fields the browser shows, which fields the validator demands, and
which lines the script template emits. Keeping that in one place is what stops
the form offering (say) a surface tension field for a non-membrane system.

``openmm`` records the API symbol an option maps to. ``tests/test_options.py``
resolves every one of them against the installed OpenMM, so a name that changes
upstream fails a test rather than a user's run three steps later.

Help text is a full sentence explaining *why* the default is what it is, not a
restatement of the label. A user who has never run MD should be able to accept
every default and get a physically sensible trajectory; a user who knows what
they are doing should be able to find out why we disagree with them.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

FLEXAPPEAL_VERSION = "0.1.0"

# The OpenMM release the generated scripts and the pinned bundle environment are
# written against. Bumping this means re-running tests/test_openmm_symbols.py,
# which resolves every `openmm=` symbol below against the installed package.
OPENMM_VERSION = "8.5"


# Water model -> (force-field XML, the `model=` argument to Modeller.addSolvent).
#
# These are two different things and conflating them is a silent-corruption bug.
# addSolvent's `model` argument only accepts 'tip3p', 'spce', 'tip4pew', 'tip5p'
# and 'swm4ndp' -- it uses it purely to choose the water *geometry* to pack into
# the box. The actual parameters come from whichever water XML the ForceField
# loaded, matched by residue name.
#
# So OPC (a four-site model) is placed with tip4pew geometry and parameterised
# from amber14/opc.xml, and OPC3 (three-site) is placed with tip3p geometry.
# Passing model='opc' does not raise -- it is simply not one of the recognised
# names -- which is exactly why this mapping is explicit and tested by actually
# running addSolvent for every entry rather than merely asserting the files exist.
WATER_MODEL_XML: dict[str, tuple[str, str]] = {
    "tip3p":        ("amber14/tip3p.xml", "tip3p"),
    "tip3pfb":      ("amber14/tip3pfb.xml", "tip3p"),
    "spce":         ("amber14/spce.xml", "spce"),
    "tip4pew":      ("amber14/tip4pew.xml", "tip4pew"),
    "tip4pfb":      ("amber14/tip4pfb.xml", "tip4pew"),
    "opc":          ("amber14/opc.xml", "tip4pew"),
    "opc3":         ("amber14/opc3.xml", "tip3p"),
    "charmm_tip3p": ("charmm36/water.xml", "tip3p"),
}


# ===========================================================================
#  Data model
# ===========================================================================


@dataclass(frozen=True)
class Choice:
    """One option in a select/multiselect."""

    value: str
    label: str
    help: str = ""
    # A choice can itself be conditional -- e.g. the CHARMM-modified water model
    # is only meaningful once the CHARMM36 protein force field is selected.
    requires: str | None = None
    experimental: bool = False


@dataclass(frozen=True)
class Option:
    """One user-facing parameter."""

    id: str
    group: str
    label: str
    widget: str  # text | number | int | select | multiselect | checkbox | textarea | file | pdbid
    default: Any
    help: str
    units: str | None = None
    choices: tuple[Choice, ...] = ()
    minimum: float | None = None
    maximum: float | None = None
    step: float | None = None
    placeholder: str | None = None
    # The OpenMM (or PDBFixer/MDTraj) symbol this maps to, for the resolver test.
    openmm: str | None = None
    # Predicate over other option ids; the field is hidden and skipped when false.
    requires: str | None = None
    # Collapsed behind "Advanced" in the UI. Sensible defaults; rarely touched.
    advanced: bool = False
    # Shown with a warning badge. Works, but the user should know what they are
    # signing up for (third-party plugin, deprecated back end, expensive method).
    experimental: bool = False
    # Choices are not known until a structure is parsed (chains, heteroatoms).
    # inspect.py fills them in per-request; the registry cannot.
    dynamic: bool = False

    def choice_values(self) -> tuple[str, ...]:
        return tuple(c.value for c in self.choices)


@dataclass(frozen=True)
class Group:
    """A section of the Prepare form."""

    id: str
    title: str
    blurb: str
    icon: str
    advanced: bool = False


# ===========================================================================
#  Groups -- the accordion sections of the Prepare tab, in form order
# ===========================================================================

GROUPS: tuple[Group, ...] = (
    Group("job", "Job", "What this run is called and where its output lands.", "🏷️"),
    Group("input", "Input structure", "The starting coordinates and what to keep from them.", "🧬"),
    Group("repair", "Structure repair", "PDBFixer: gaps, missing atoms, protonation.", "🩹"),
    Group("forcefield", "Force field", "The physics. Protein, water, lipid and ligand parameters.", "⚛️"),
    Group("solvation", "Solvation and box", "Explicit water, counter-ions and periodic box.", "💧"),
    Group("membrane", "Membrane", "Lipid bilayer construction for membrane proteins.", "🧱"),
    Group("system", "System", "Cutoffs, constraints and long-range electrostatics.", "🔧"),
    Group("integrator", "Integrator", "How the equations of motion are advanced.", "⏱️"),
    Group("barostat", "Pressure control", "Constant-pressure coupling, if any.", "🎈"),
    Group("restraints", "Restraints", "Holding parts of the system in place.", "📌"),
    Group("protocol", "Protocol", "Minimise, heat, equilibrate, produce.", "📈"),
    Group("output", "Output and reporting", "What gets written, how often, and how big it will be.", "💾"),
    Group("platform", "Platform and performance", "Which hardware back end runs the maths.", "🚀"),
    Group("analysis", "Analysis", "What the local run computes before packing the results file.", "📊"),
)

GROUP_ORDER: tuple[str, ...] = tuple(g.id for g in GROUPS)


# ===========================================================================
#  Shared choice sets
# ===========================================================================

_POSITIVE_IONS = (
    Choice("Na+", "Na⁺ (sodium)", "The default counter-ion for almost all biomolecular work."),
    Choice("K+", "K⁺ (potassium)", "Closer to intracellular conditions than sodium."),
    Choice("Li+", "Li⁺ (lithium)"),
    Choice("Rb+", "Rb⁺ (rubidium)"),
    Choice("Cs+", "Cs⁺ (caesium)"),
)

_NEGATIVE_IONS = (
    Choice("Cl-", "Cl⁻ (chloride)", "The default counter-ion for almost all biomolecular work."),
    Choice("Br-", "Br⁻ (bromide)"),
    Choice("F-", "F⁻ (fluoride)"),
    Choice("I-", "I⁻ (iodide)"),
)

_ATOM_SELECTIONS = (
    Choice("protein", "Protein", "All protein atoms including hydrogens."),
    Choice("protein_heavy", "Protein heavy atoms", "Protein excluding hydrogens: half the size, and hydrogens are rarely analysed."),
    Choice("backbone", "Backbone", "N, CA, C, O only."),
    Choice("ca", "Cα only", "One atom per residue: the smallest useful representation."),
    Choice("protein_ligand", "Protein and ligands", "Everything except water, ions and lipids."),
    Choice("all", "Everything", "Including water, ions and lipids. Large."),
    Choice("custom", "Custom selection", "An MDTraj selection expression, entered below."),
)


# ===========================================================================
#  The options
# ===========================================================================

OPTIONS: tuple[Option, ...] = (
    # -----------------------------------------------------------------------
    #  1. Job metadata
    # -----------------------------------------------------------------------
    Option(
        "job_name", "job", "Job name", "text", "flexappeal_run",
        "Used for the bundle filename, the output directory and the title on every "
        "analysis plot. Letters, digits, hyphens and underscores only, because it "
        "becomes a path on the user's machine.",
        placeholder="lysozyme_apo_100ns",
    ),
    Option(
        "job_description", "job", "Description", "textarea", "",
        "Free text carried verbatim into the results manifest. Six months from now "
        "this is the only thing that will tell you why you ran this.",
        placeholder="Apo lysozyme, 100 ns, ff14SB/TIP3P, checking loop flexibility",
    ),
    Option(
        "job_author", "job", "Author", "text", "",
        "Recorded in the manifest for provenance. Optional.",
        placeholder="M. Deller",
    ),
    Option(
        "output_dir", "job", "Output directory", "text", "./flexappeal_output",
        "Where the run writes on your machine, relative to wherever you put the "
        "bundle. Everything the run produces stays inside this one directory so it "
        "is trivial to archive or delete.",
        advanced=True,
    ),

    # -----------------------------------------------------------------------
    #  2. Input structure
    # -----------------------------------------------------------------------
    Option(
        "input_source", "input", "Structure source", "select", "upload",
        "Where the starting coordinates come from. Fetching by ID is reproducible "
        "and records the exact accession in the manifest; uploading is right for "
        "anything unreleased, modelled or hand-edited.",
        choices=(
            Choice("upload", "Upload a file", "PDB or PDBx/mmCIF, up to 25 MB."),
            Choice("rcsb", "Fetch from the RCSB PDB", "By four-character PDB ID."),
            Choice("opm", "Fetch from OPM (pre-oriented)",
                   "The Orientations of Proteins in Membranes database returns the structure "
                   "already positioned in a bilayer frame, which is exactly what addMembrane needs."),
            Choice("alphafold", "Fetch from the AlphaFold DB", "By UniProt accession."),
        ),
    ),
    Option(
        "input_file", "input", "Structure file", "file", None,
        "A PDB or PDBx/mmCIF file. It is parsed server-side only to populate the "
        "options below; the copy embedded in your bundle is the sanitised version "
        "produced by the choices you make here.",
        requires="input_source == 'upload'",
    ),
    Option(
        "pdb_id", "input", "PDB ID", "pdbid", "",
        "Four characters, for example 1AKI. Fetched over HTTPS when the bundle is "
        "built, not at run time, so your run does not depend on the RCSB being up.",
        requires="input_source in ('rcsb', 'opm')",
        placeholder="1AKI",
    ),
    Option(
        "uniprot_id", "input", "UniProt accession", "text", "",
        "For example P00698. Retrieves the AlphaFold model; remember that predicted "
        "structures have no waters, no ligands and unreliable loop conformations.",
        requires="input_source == 'alphafold'",
        placeholder="P00698",
    ),
    Option(
        "assembly", "input", "Assembly", "select", "asymmetric",
        "Crystal structures are deposited as an asymmetric unit, which is often not "
        "the biological molecule. If your protein is a physiological dimer, simulating "
        "one chain is simulating the wrong thing.",
        choices=(
            Choice("asymmetric", "Asymmetric unit", "Exactly what is in the deposited file."),
            Choice("biological", "First biological assembly", "Built from the deposited transformation matrices."),
        ),
        requires="input_source in ('rcsb', 'opm')",
    ),
    Option(
        "model_index", "input", "Model", "int", 1,
        "NMR ensembles and some predictions contain many models. Model 1 is "
        "conventionally the representative one.",
        minimum=1, advanced=True,
    ),
    Option(
        "chains", "input", "Chains to keep", "multiselect", ["*"],
        "Populated from the structure once it is loaded. Dropping chains you do not "
        "need is the single cheapest way to make a simulation faster.",
        dynamic=True,
    ),
    Option(
        "altloc", "input", "Alternate locations", "select", "occupancy",
        "Crystal structures often model two conformations for the same side chain. "
        "OpenMM cannot simulate both, so exactly one has to be chosen.",
        choices=(
            Choice("occupancy", "Highest occupancy", "Whichever the crystallographer refined as dominant."),
            Choice("first", "First listed", "Usually altloc A."),
            Choice("a", "Always altloc A"),
            Choice("b", "Always altloc B"),
        ),
        advanced=True,
    ),
    Option(
        "keep_waters", "input", "Keep crystallographic waters", "checkbox", True,
        "Buried and bridging waters are frequently structurally important, and "
        "deleting them can collapse an active site during equilibration. Bulk waters "
        "are replaced by the solvation step regardless.",
    ),
    Option(
        "water_shell", "input", "Water shell radius", "number", 5.0,
        "Only crystallographic waters within this distance of the protein are kept. "
        "Anything further out is bulk solvent that addSolvent will regenerate anyway.",
        units="Å", minimum=0.0, maximum=20.0, step=0.5,
        requires="keep_waters", advanced=True,
    ),
    Option(
        "keep_heteroatoms", "input", "Heteroatoms to keep", "multiselect", [],
        "Cofactors, ligands, metals and modified residues found in the file. Anything "
        "kept here must have force-field parameters, which is what the ligand section "
        "below is for. Anything not listed is deleted.",
        dynamic=True,
    ),
    Option(
        "strip_hydrogens", "input", "Discard input hydrogens", "checkbox", True,
        "Deposited hydrogen positions are usually either absent or placed by a "
        "different program with different conventions. Re-adding them with PDBFixer "
        "at a known pH is more consistent than trusting what arrived.",
        advanced=True,
    ),
    Option(
        "disulfides", "input", "Disulfide bonds", "select", "auto",
        "Getting these wrong changes the fold. Automatic detection finds SG–SG pairs "
        "closer than the cutoff below, which is reliable for well-resolved structures.",
        choices=(
            Choice("auto", "Detect automatically", "By SG–SG distance."),
            Choice("manual", "Specify manually", "Enter the residue pairs yourself."),
            Choice("none", "None", "Treat every cysteine as free thiol."),
        ),
    ),
    Option(
        "disulfide_cutoff", "input", "Disulfide detection cutoff", "number", 2.5,
        "An S–S bond is about 2.05 Å; 2.5 Å catches slightly strained or "
        "lower-resolution geometry without picking up unbonded neighbours.",
        units="Å", minimum=1.5, maximum=4.0, step=0.05,
        requires="disulfides == 'auto'", advanced=True,
    ),
    Option(
        "disulfide_pairs", "input", "Disulfide pairs", "textarea", "",
        "One pair per line as CHAIN:RESID-CHAIN:RESID, for example A:6-A:127.",
        requires="disulfides == 'manual'",
        placeholder="A:6-A:127\nA:30-A:115",
    ),
    Option(
        "termini_caps", "input", "Terminal caps", "select", "charged",
        "Charged termini are correct for a complete protein. Neutral ACE/NME caps are "
        "correct for a fragment or a single domain excised from something larger, "
        "where an artificial charge would distort the local electrostatics.",
        choices=(
            Choice("charged", "Charged (NH₃⁺ / COO⁻)", "The right choice for an intact protein."),
            Choice("capped", "Capped (ACE / NME)", "The right choice for a fragment."),
        ),
        advanced=True,
    ),

    # -----------------------------------------------------------------------
    #  3. Structure repair (PDBFixer)
    # -----------------------------------------------------------------------
    Option(
        "fix_missing_residues", "repair", "Build missing residues", "checkbox", True,
        "Disordered loops are absent from most crystal structures. Leaving a gap "
        "means simulating a chain break, which is almost never what you want; "
        "PDBFixer rebuilds them, though the conformation it invents is a guess.",
        openmm="pdbfixer.PDBFixer.findMissingResidues",
    ),
    Option(
        "max_gap_length", "repair", "Longest gap to build", "int", 25,
        "A backstop against modelling something unusable, not an endorsement of "
        "everything below it. A rebuilt loop is placed by geometry rather than "
        "physics: up to about eight residues it is usually reasonable, by fifteen "
        "it is a guess, and by twenty-five it is fiction that will dominate your "
        "RMSD and is the likeliest thing to destabilise the run. The default "
        "admits structures like oncostatin M (a 21-residue gap) because refusing "
        "them outright helps nobody, but the run warns loudly and you should treat "
        "that region as decoration.",
        units="residues", minimum=0, maximum=100,
        requires="fix_missing_residues", advanced=True,
    ),
    Option(
        "build_terminal_residues", "repair", "Build missing terminal residues", "checkbox", False,
        "Unresolved termini are usually genuinely disordered rather than merely "
        "unmodelled. Adding a flapping tail costs atoms and adds noise to global "
        "metrics without adding information.",
        requires="fix_missing_residues", advanced=True,
    ),
    Option(
        "replace_nonstandard", "repair", "Replace non-standard residues", "checkbox", True,
        "Selenomethionine, phosphoserine and similar have no parameters in the "
        "standard force fields. PDBFixer swaps them for their closest standard "
        "equivalent, which for MSE→MET is exact and for a phospho-residue is not.",
        openmm="pdbfixer.PDBFixer.replaceNonstandardResidues",
    ),
    Option(
        "add_missing_atoms", "repair", "Add missing heavy atoms", "checkbox", True,
        "Partially resolved side chains are common at the surface. Without this the "
        "force field will refuse to build a template for them.",
        openmm="pdbfixer.PDBFixer.addMissingAtoms",
    ),
    Option(
        "ph", "repair", "pH", "number", 7.4,
        "Sets the protonation of ionisable side chains. 7.4 is physiological; use the "
        "pH of your actual experiment if you are comparing to one, since a histidine "
        "flipping protonation state can change a binding site entirely.",
        minimum=0.0, maximum=14.0, step=0.1,
        openmm="pdbfixer.PDBFixer.addMissingHydrogens",
    ),
    Option(
        "protonation_overrides", "repair", "Protonation overrides", "textarea", "",
        "One per line as CHAIN:RESID=STATE, for example A:64=HIP. Overrides the pH "
        "rule for individual residues, which matters most for catalytic histidines "
        "and buried carboxylates whose real pKa is nowhere near the model value. "
        "Valid states: HID, HIE, HIP, ASH, GLH, LYN, CYM, CYX.",
        advanced=True,
        placeholder="A:64=HIP\nA:35=GLH",
    ),
    Option(
        "mutations", "repair", "Mutations", "textarea", "",
        "One per line as CHAIN:WT-RESID-MUT using three-letter codes, for example "
        "A:ALA-57-GLY. Applied before any repair, so a rebuilt loop reflects the mutant.",
        openmm="pdbfixer.PDBFixer.applyMutations", advanced=True,
        placeholder="A:ALA-57-GLY",
    ),

    # -----------------------------------------------------------------------
    #  4. Force field
    # -----------------------------------------------------------------------
    Option(
        "protein_ff", "forcefield", "Protein force field", "select", "amber14-all.xml",
        "ff14SB (inside amber14-all) is the most widely validated modern protein force "
        "field and the safe default. CHARMM36m is its main rival and is the better "
        "choice for disordered regions and for membrane work, where the lipid "
        "parameters are natively matched.",
        choices=(
            Choice("amber14-all.xml", "AMBER ff14SB", "The default. Well validated for folded globular proteins."),
            Choice("amber/protein.ff19SB.xml", "AMBER ff19SB",
                   "Improved amino-acid-specific backbone profiles. Designed for and effectively "
                   "requires OPC water, which is selected automatically below. Ships with the "
                   "openmmforcefields package rather than OpenMM itself, so choosing it adds that "
                   "dependency to your bundle even for an apo run."),
            Choice("charmm36.xml", "CHARMM36m",
                   "Better balanced for intrinsically disordered proteins, and the native partner "
                   "for the CHARMM36 lipid parameters used by the membrane builder."),
            Choice("amber99sbildn.xml", "AMBER ff99SB-ILDN", "The previous generation. Use for continuity with older work."),
            Choice("amber03.xml", "AMBER ff03"),
            Choice("amber10.xml", "AMBER ff10"),
            Choice("amber96.xml", "AMBER ff96"),
            Choice("amoeba2018.xml", "AMOEBA 2018",
                   "Polarisable. Far more accurate in principle and roughly two orders of magnitude "
                   "slower in practice; only viable for very small systems.",
                   experimental=True),
        ),
        openmm="openmm.app.ForceField",
    ),
    Option(
        "water_model", "forcefield", "Water model", "select", "tip3p",
        "Must be compatible with the protein force field, and the app enforces that. "
        "TIP3P is the AMBER default and by far the cheapest. OPC and TIP4P-Ew "
        "reproduce bulk water properties considerably better at roughly 1.3× the cost, "
        "and matter most when you care about solvation or disordered states.",
        choices=(
            Choice("tip3p", "TIP3P", "Three sites. The default partner for the AMBER force fields.",
                   requires="protein_ff != 'charmm36.xml'"),
            Choice("tip3pfb", "TIP3P-FB", "Force-balance reparameterised TIP3P: better density and dielectric at the same cost.",
                   requires="protein_ff != 'charmm36.xml'"),
            Choice("spce", "SPC/E", "Three sites, better bulk properties than plain TIP3P.",
                   requires="protein_ff != 'charmm36.xml'"),
            Choice("tip4pew", "TIP4P-Ew", "Four sites. Excellent bulk water; the extra site costs about 25%.",
                   requires="protein_ff != 'charmm36.xml'"),
            Choice("tip4pfb", "TIP4P-FB", "Force-balance four-site model.",
                   requires="protein_ff != 'charmm36.xml'"),
            Choice("opc", "OPC", "Four sites. The best-performing rigid model for biomolecular work, and required by ff19SB.",
                   requires="protein_ff != 'charmm36.xml'"),
            Choice("opc3", "OPC3", "Three-site cousin of OPC: most of the accuracy at TIP3P cost.",
                   requires="protein_ff != 'charmm36.xml'"),
            # TIP5P is deliberately absent: OpenMM ships no amber14/tip5p.xml, so it
            # cannot be paired with the AMBER force fields offered here without
            # mixing force-field generations. It is also rarely worth its cost for
            # protein work. Confirmed by tests/test_openmm_symbols.py, which is how
            # the omission stays honest rather than becoming a stale comment.
            Choice("charmm_tip3p", "CHARMM-modified TIP3P",
                   "TIP3P with Lennard-Jones terms on the hydrogens. The only correct partner for CHARMM36.",
                   requires="protein_ff == 'charmm36.xml'"),
        ),
        requires="solvent_mode == 'explicit'",
    ),
    Option(
        "solvent_mode", "forcefield", "Solvent treatment", "select", "explicit",
        "Explicit water is the physically correct choice and what almost all published "
        "protein MD uses. Implicit solvent is perhaps ten times faster and useful for "
        "long-timescale conformational sampling, but it misrepresents solvent-mediated "
        "interactions and salt effects.",
        choices=(
            Choice("explicit", "Explicit solvent", "A real box of water molecules and ions."),
            Choice("implicit", "Implicit solvent (GB)", "A continuum dielectric. Much faster, much less accurate."),
            Choice("vacuum", "Vacuum", "No solvent at all. For testing and for gas-phase questions only.",
                   experimental=True),
        ),
    ),
    Option(
        "implicit_model", "forcefield", "Implicit solvent model", "select", "implicit/gbn2.xml",
        "GBn2 is the most accurate of the generalised Born variants OpenMM ships and "
        "the one to use unless you are reproducing older work.",
        choices=(
            Choice("implicit/gbn2.xml", "GBn2", "The recommended generalised Born model."),
            Choice("implicit/obc2.xml", "OBC2", "Onufriev-Bashford-Case II. The long-standing default elsewhere."),
            Choice("implicit/obc1.xml", "OBC1"),
            Choice("implicit/gbn.xml", "GBn"),
            Choice("implicit/hct.xml", "HCT", "The oldest and least accurate."),
        ),
        requires="solvent_mode == 'implicit'",
    ),
    Option(
        "implicit_salt", "forcefield", "Implicit salt concentration", "number", 0.15,
        "Screens electrostatics via the Debye-Hückel term. 0.15 M approximates "
        "physiological ionic strength.",
        units="M", minimum=0.0, maximum=2.0, step=0.01,
        requires="solvent_mode == 'implicit'",
    ),
    Option(
        "solute_dielectric", "forcefield", "Solute dielectric", "number", 1.0,
        "The interior dielectric constant. 1.0 is the standard choice; raising it is a "
        "crude way to mimic electronic polarisation.",
        minimum=1.0, maximum=20.0, step=0.5,
        requires="solvent_mode == 'implicit'", advanced=True,
    ),
    Option(
        "solvent_dielectric", "forcefield", "Solvent dielectric", "number", 78.5,
        "The exterior dielectric constant. 78.5 is water at room temperature.",
        minimum=1.0, maximum=100.0, step=0.5,
        requires="solvent_mode == 'implicit'", advanced=True,
    ),
    Option(
        "ligand_ff", "forcefield", "Small-molecule force field", "select", "openff-2.2.1",
        "Applies to every ligand and cofactor kept above. OpenFF Sage is the current "
        "best general-purpose small-molecule force field and is actively maintained. "
        "GAFF2 is the classical AMBER partner and the safer choice if you need to match "
        "older AMBER work.",
        choices=(
            Choice("openff-2.2.1", "OpenFF 2.2.1 (Sage)", "The recommended modern choice."),
            Choice("openff-2.0.0", "OpenFF 2.0.0 (Sage)", "For continuity with work published against 2.0."),
            Choice("gaff-2.11", "GAFF 2.11", "The classical AMBER general force field."),
            Choice("espaloma-0.3.2", "Espaloma 0.3.2",
                   "A graph-neural-network force field. Promising accuracy, much less field experience.",
                   experimental=True),
        ),
        openmm="openmmforcefields.generators.SystemGenerator",
        requires="has_ligands",
    ),
    Option(
        "ligand_charge_method", "forcefield", "Ligand partial charges", "select", "am1bcc",
        "AM1-BCC is the standard and what both GAFF2 and OpenFF were parameterised "
        "against; it takes a minute or two per ligand. Gasteiger is instant and "
        "noticeably worse, and should only be used to get a pipeline working.",
        choices=(
            Choice("am1bcc", "AM1-BCC", "The correct choice. Requires AmberTools, which the bundle installs."),
            Choice("am1bccelf10", "AM1-BCC ELF10", "Conformer-averaged AM1-BCC. Slower, more reproducible."),
            Choice("gasteiger", "Gasteiger", "Fast and crude. For smoke-testing only."),
            Choice("from_file", "Take from the input file", "Use charges already present in the SDF/MOL2."),
        ),
        requires="has_ligands", advanced=True,
    ),
    Option(
        "has_ligands", "forcefield", "System contains ligands or cofactors", "checkbox", False,
        "Set automatically when you keep a heteroatom above. Turning it on reveals the "
        "small-molecule parameterisation options.",
        advanced=True,
    ),
    Option(
        "extra_ff_xml", "forcefield", "Additional force-field XML", "file", None,
        "An OpenMM force-field XML appended after the ones selected above, for custom "
        "residue templates or parameters you have derived yourself.",
        advanced=True,
    ),

    # -----------------------------------------------------------------------
    #  5. Solvation and box
    # -----------------------------------------------------------------------
    Option(
        "box_shape", "solvation", "Box shape", "select", "dodecahedron",
        "A rhombic dodecahedron encloses the same minimum-image distance in about 71% "
        "of the volume of a cube, so roughly 30% fewer water molecules for identical "
        "physics. The only reason to prefer a cube is a downstream tool that cannot "
        "read triclinic boxes.",
        choices=(
            Choice("dodecahedron", "Rhombic dodecahedron", "The efficient default: ~29% fewer waters than a cube."),
            Choice("octahedron", "Truncated octahedron", "Almost as efficient; the traditional AMBER choice."),
            Choice("cube", "Cube", "Simplest, largest, slowest. Choose it only if something downstream demands it."),
        ),
        openmm="openmm.app.Modeller.addSolvent",
        requires="solvent_mode == 'explicit' and not use_membrane",
    ),
    Option(
        "box_sizing", "solvation", "Box sizing", "select", "padding",
        "Padding scales with the molecule and is what you want in almost every case. "
        "Explicit vectors are for reproducing someone else's box exactly.",
        choices=(
            Choice("padding", "Solvent padding", "A minimum solvent shell around the solute."),
            Choice("vectors", "Explicit box vectors", "Give the three box lengths yourself."),
            Choice("num_waters", "Fixed water count", "Add exactly this many water molecules."),
        ),
        requires="solvent_mode == 'explicit' and not use_membrane", advanced=True,
    ),
    Option(
        "padding", "solvation", "Solvent padding", "number", 1.2,
        "The minimum distance from any solute atom to the box edge. It must exceed the "
        "non-bonded cutoff, or a molecule can interact with its own periodic image -- "
        "a physically meaningless artefact that is easy to miss in the output. The "
        "default leaves 0.2 nm of headroom above the default 1.0 nm cutoff, because "
        "the box shrinks under the barostat during equilibration.",
        units="nm", minimum=0.5, maximum=5.0, step=0.1,
        requires="solvent_mode == 'explicit' and box_sizing == 'padding' and not use_membrane",
    ),
    Option(
        "box_vectors", "solvation", "Box vectors", "text", "",
        "Three lengths in nanometres, comma separated, for example 6.0, 6.0, 6.0.",
        units="nm",
        requires="solvent_mode == 'explicit' and box_sizing == 'vectors'",
        placeholder="6.0, 6.0, 6.0",
    ),
    Option(
        "num_waters", "solvation", "Number of waters", "int", 10000,
        "Exact water count. The box is sized to accommodate them.",
        minimum=100, maximum=1000000,
        requires="solvent_mode == 'explicit' and box_sizing == 'num_waters'",
    ),
    Option(
        "positive_ion", "solvation", "Positive ion", "select", "Na+",
        "The cation used for neutralisation and for reaching the target ionic strength.",
        choices=_POSITIVE_IONS,
        requires="solvent_mode == 'explicit'",
    ),
    Option(
        "negative_ion", "solvation", "Negative ion", "select", "Cl-",
        "The anion used for neutralisation and for reaching the target ionic strength.",
        choices=_NEGATIVE_IONS,
        requires="solvent_mode == 'explicit'",
    ),
    Option(
        "ionic_strength", "solvation", "Ionic strength", "number", 0.15,
        "Salt added on top of neutralisation. 0.15 M is roughly physiological. Zero "
        "gives you only the counter-ions needed to neutralise the net charge, which "
        "under-screens surface electrostatics.",
        units="M", minimum=0.0, maximum=2.0, step=0.01,
        requires="solvent_mode == 'explicit'",
    ),
    Option(
        "neutralize", "solvation", "Neutralise the system", "checkbox", True,
        "A non-zero net charge under particle-mesh Ewald is handled by an implicit "
        "uniform background charge, which introduces real artefacts near the solute. "
        "There is essentially never a good reason to turn this off.",
        requires="solvent_mode == 'explicit'", advanced=True,
    ),

    # -----------------------------------------------------------------------
    #  6. Membrane
    # -----------------------------------------------------------------------
    Option(
        "use_membrane", "membrane", "Embed in a lipid bilayer", "checkbox", False,
        "Builds a bilayer around the protein using OpenMM's own addMembrane. The "
        "protein must already be oriented with the membrane normal along z, which is "
        "exactly what fetching from OPM gives you.",
        openmm="openmm.app.Modeller.addMembrane",
    ),
    Option(
        "lipid_type", "membrane", "Lipid", "select", "POPC",
        "POPC is the standard general-purpose model membrane and the right default "
        "for a plasma-membrane protein. Mixed compositions are not supported by "
        "addMembrane; build those elsewhere and upload the result.",
        choices=(
            Choice("POPC", "POPC", "Palmitoyl-oleoyl phosphatidylcholine. The usual choice."),
            Choice("POPE", "POPE", "Phosphatidylethanolamine. Common for bacterial inner membranes."),
            Choice("DOPC", "DOPC", "Dioleoyl phosphatidylcholine. More fluid than POPC."),
            Choice("DPPC", "DPPC", "Dipalmitoyl phosphatidylcholine. Gel phase at room temperature."),
            Choice("DMPC", "DMPC", "Dimyristoyl phosphatidylcholine. Short tails, thin bilayer."),
            Choice("DLPC", "DLPC", "Dilauroyl phosphatidylcholine. Thinner still."),
            Choice("DLPE", "DLPE", "Dilauroyl phosphatidylethanolamine."),
        ),
        requires="use_membrane",
    ),
    Option(
        "membrane_orientation", "membrane", "Orientation", "select", "opm",
        "The bilayer is built in the xy plane at a fixed z, so a protein that is not "
        "already aligned to the membrane normal will be built into the membrane "
        "sideways. This is the most common way a membrane setup goes silently wrong.",
        choices=(
            Choice("opm", "Fetch pre-oriented from OPM", "The reliable route. Uses the OPM database's own alignment."),
            Choice("preoriented", "Already oriented", "Trust the input file's frame as given."),
            Choice("principal", "Align principal axis to z",
                   "A geometric guess. Works for a clean single bundle of transmembrane helices and not much else.",
                   experimental=True),
        ),
        requires="use_membrane",
    ),
    Option(
        "membrane_center_z", "membrane", "Membrane centre (z)", "number", 0.0,
        "Where the bilayer midplane sits along z, relative to the structure's own "
        "frame. Zero is correct for anything fetched from OPM.",
        units="nm", minimum=-10.0, maximum=10.0, step=0.1,
        requires="use_membrane", advanced=True,
    ),
    Option(
        "membrane_padding", "membrane", "Minimum padding", "number", 1.0,
        "Water thickness above and below the bilayer. Too little and the protein's "
        "extramembrane domains see their own periodic image through the solvent.",
        units="nm", minimum=0.5, maximum=5.0, step=0.1,
        requires="use_membrane",
    ),

    # -----------------------------------------------------------------------
    #  7. System
    # -----------------------------------------------------------------------
    Option(
        "nonbonded_method", "system", "Non-bonded method", "select", "PME",
        "Particle-mesh Ewald is the only correct treatment of long-range "
        "electrostatics in a periodic box and is what essentially all modern "
        "biomolecular MD uses. The cutoff methods are appropriate only for implicit "
        "solvent or vacuum.",
        choices=(
            Choice("PME", "PME", "Particle-mesh Ewald. The default for any periodic system.",
                   requires="solvent_mode == 'explicit'"),
            Choice("LJPME", "LJPME", "PME applied to dispersion as well as electrostatics. More accurate, ~15% slower.",
                   requires="solvent_mode == 'explicit'"),
            Choice("Ewald", "Ewald", "Exact but scales badly. For small reference calculations only.",
                   requires="solvent_mode == 'explicit'"),
            Choice("CutoffPeriodic", "Cutoff (periodic)", "Reaction-field style. Not appropriate for charged solutes."),
            Choice("CutoffNonPeriodic", "Cutoff (non-periodic)", "For implicit solvent.",
                   requires="solvent_mode != 'explicit'"),
            Choice("NoCutoff", "No cutoff", "All pairs. Only tractable for small implicit-solvent or vacuum systems.",
                   requires="solvent_mode != 'explicit'"),
        ),
        openmm="openmm.app.PME",
    ),
    Option(
        "nonbonded_cutoff", "system", "Non-bonded cutoff", "number", 1.0,
        "Where direct-space non-bonded interactions stop. 1.0 nm is standard for the "
        "AMBER force fields; CHARMM36 was parameterised at 1.2 nm with a switching "
        "function and should use that. Must be less than half the shortest box vector.",
        units="nm", minimum=0.6, maximum=2.0, step=0.05,
        openmm="openmm.app.ForceField.createSystem",
    ),
    Option(
        "use_switching", "system", "Use a switching function", "checkbox", False,
        "Smoothly tapers Lennard-Jones to zero before the cutoff instead of truncating "
        "it. Required for CHARMM36 as parameterised; unnecessary for AMBER, which uses "
        "a long-range dispersion correction instead.",
        advanced=True,
    ),
    Option(
        "switch_distance", "system", "Switching distance", "number", 0.9,
        "Where the taper begins. Conventionally 0.1 nm inside the cutoff.",
        units="nm", minimum=0.4, maximum=1.9, step=0.05,
        requires="use_switching", advanced=True,
    ),
    Option(
        "ewald_error_tolerance", "system", "Ewald error tolerance", "number", 0.0005,
        "Controls the PME grid density. The default is accurate enough that reducing "
        "it further costs time without changing results.",
        minimum=1e-6, maximum=1e-2, step=1e-5,
        requires="nonbonded_method in ('PME', 'LJPME', 'Ewald')", advanced=True,
    ),
    Option(
        "constraints", "system", "Bond constraints", "select", "HBonds",
        "Constraining bonds to hydrogen removes the fastest vibrations in the system, "
        "which is what allows a 2 fs timestep instead of 0.5 fs. It is standard "
        "practice and has a negligible effect on everything people actually measure.",
        choices=(
            Choice("HBonds", "Bonds to hydrogen", "The standard choice. Permits a 2 fs timestep."),
            Choice("AllBonds", "All bonds", "Permits slightly longer steps; distorts some flexibility."),
            Choice("HAngles", "All bonds and H-X-H angles", "Aggressive. Used with very long timesteps."),
            Choice("None", "None", "Requires a timestep of 0.5 fs or less."),
        ),
        openmm="openmm.app.HBonds",
    ),
    Option(
        "rigid_water", "system", "Rigid water", "checkbox", True,
        "Water models are parameterised as rigid bodies. Making them flexible is both "
        "slower and wrong relative to how the model was fitted.",
        requires="solvent_mode == 'explicit'", advanced=True,
    ),
    Option(
        "use_hmr", "system", "Hydrogen mass repartitioning", "checkbox", True,
        "Moves mass from heavy atoms onto the hydrogens bonded to them, slowing the "
        "fastest remaining motions so the timestep can double to 4 fs. Total mass and "
        "equilibrium properties are unchanged. This is close to a free 2× speedup and "
        "is the single most effective performance option on this page.",
        openmm="openmm.app.ForceField.createSystem",
    ),
    Option(
        "hydrogen_mass", "system", "Hydrogen mass", "number", 1.5,
        "The repartitioned hydrogen mass. 1.5 amu is the conservative, well-tested "
        "value that supports a 4 fs step; 4.0 amu permits longer steps but has been "
        "shown to perturb some kinetics.",
        units="amu", minimum=1.0, maximum=4.0, step=0.5,
        requires="use_hmr", advanced=True,
    ),
    Option(
        "remove_cm_motion", "system", "Remove centre-of-mass motion", "checkbox", True,
        "Without this the whole system slowly acquires a net drift velocity, which "
        "steals kinetic energy from the thermostat and corrupts diffusion measurements.",
        openmm="openmm.CMMotionRemover", advanced=True,
    ),
    Option(
        "dispersion_correction", "system", "Long-range dispersion correction", "checkbox", True,
        "An analytic estimate of the Lennard-Jones energy beyond the cutoff. Needed "
        "for correct densities in constant-pressure simulations. Turn it off only when "
        "using a switching function, or you are counting the same energy twice.",
        openmm="openmm.NonbondedForce.setUseDispersionCorrection", advanced=True,
    ),

    # -----------------------------------------------------------------------
    #  8. Integrator
    # -----------------------------------------------------------------------
    Option(
        "integrator", "integrator", "Integrator", "select", "LangevinMiddle",
        "LangevinMiddle is OpenMM's recommended integrator: it thermostats and "
        "integrates in one step, and its middle-scheme discretisation gives noticeably "
        "better configurational sampling at large timesteps than the older Langevin "
        "integrator.",
        choices=(
            Choice("LangevinMiddle", "Langevin (middle scheme)", "The recommended default for constant-temperature MD."),
            Choice("Langevin", "Langevin (leapfrog)", "The older scheme. Use for continuity with previous work."),
            Choice("NoseHoover", "Nosé-Hoover", "Deterministic thermostat. Gives a true canonical ensemble without stochastic forces."),
            Choice("Verlet", "Verlet", "Constant energy. Needs a separate thermostat, or none for NVE."),
            Choice("Brownian", "Brownian", "Overdamped. For coarse-grained or implicit-solvent sampling.", experimental=True),
            Choice("VariableLangevin", "Variable-timestep Langevin", "Adapts the step to the error tolerance.", experimental=True),
            Choice("VariableVerlet", "Variable-timestep Verlet", experimental=True),
            Choice("MTSLangevin", "Multiple-timestep Langevin",
                   "Evaluates slow forces less often. Real speedups, and real care needed with the force grouping.",
                   experimental=True),
        ),
        openmm="openmm.LangevinMiddleIntegrator",
    ),
    Option(
        "temperature", "integrator", "Temperature", "number", 310.0,
        "310 K is physiological (37 °C). Use 298 K if you are comparing against "
        "room-temperature biophysical measurements.",
        units="K", minimum=1.0, maximum=1000.0, step=1.0,
    ),
    Option(
        "friction", "integrator", "Friction coefficient", "number", 1.0,
        "How strongly the thermostat couples to the bath. 1 ps⁻¹ is the standard "
        "compromise: enough to control temperature, weak enough not to over-damp the "
        "dynamics. Values above about 5 ps⁻¹ noticeably slow real conformational motion.",
        units="ps⁻¹", minimum=0.01, maximum=100.0, step=0.1,
        requires="integrator in ('LangevinMiddle', 'Langevin', 'Brownian', 'VariableLangevin', 'MTSLangevin')",
    ),
    Option(
        "timestep", "integrator", "Timestep", "number", 4.0,
        "4 fs is safe with hydrogen mass repartitioning and H-bond constraints, and is "
        "the default here because both are on. Without HMR use 2 fs; without any "
        "constraints use 0.5 fs. Too large a step shows up as a slow, steady rise in "
        "total energy rather than an immediate crash.",
        units="fs", minimum=0.1, maximum=10.0, step=0.5,
    ),
    Option(
        "error_tolerance", "integrator", "Error tolerance", "number", 0.0001,
        "The per-step error target for the variable-timestep integrators.",
        minimum=1e-8, maximum=1e-2,
        requires="integrator in ('VariableLangevin', 'VariableVerlet')", advanced=True,
    ),
    Option(
        "use_thermostat", "integrator", "Add an Andersen thermostat", "checkbox", False,
        "Only relevant with the Verlet integrator, which has no temperature control of "
        "its own. Langevin integrators already thermostat.",
        openmm="openmm.AndersenThermostat",
        requires="integrator == 'Verlet'", advanced=True,
    ),
    Option(
        "collision_frequency", "integrator", "Collision frequency", "number", 1.0,
        "How often the Andersen thermostat randomises velocities.",
        units="ps⁻¹", minimum=0.01, maximum=100.0, step=0.1,
        requires="use_thermostat", advanced=True,
    ),
    Option(
        "random_seed", "integrator", "Random seed", "int", 0,
        "Zero means OpenMM picks a fresh seed each run. Set a specific value to make a "
        "trajectory reproducible, though bitwise reproducibility also needs a fixed "
        "platform and deterministic forces.",
        minimum=0, advanced=True,
    ),

    # -----------------------------------------------------------------------
    #  9. Barostat
    # -----------------------------------------------------------------------
    Option(
        "barostat", "barostat", "Pressure control", "select", "MonteCarlo",
        "Production runs are normally NPT, so the density equilibrates to the correct "
        "value rather than being fixed by whatever the solvation step happened to "
        "produce. Membrane systems need the membrane barostat, which couples the "
        "bilayer plane and the normal separately.",
        choices=(
            Choice("MonteCarlo", "Monte Carlo (isotropic)", "The standard choice for a soluble protein in water.",
                   requires="not use_membrane"),
            Choice("MonteCarloMembrane", "Monte Carlo (membrane)",
                   "Scales the membrane plane and normal independently, with a surface-tension term.",
                   requires="use_membrane"),
            Choice("MonteCarloAnisotropic", "Monte Carlo (anisotropic)",
                   "Independent pressure control per axis. For anisotropic systems that are not membranes."),
            Choice("MonteCarloFlexible", "Monte Carlo (flexible)", "Allows the box angles to change. For crystal-like systems.",
                   experimental=True),
            Choice("none", "None (NVT)", "Constant volume. Correct once the density has already equilibrated."),
        ),
        openmm="openmm.MonteCarloBarostat",
    ),
    Option(
        "pressure", "barostat", "Pressure", "number", 1.0,
        "1 bar is atmospheric and what essentially every biomolecular simulation uses.",
        units="bar", minimum=0.0, maximum=10000.0, step=0.5,
        requires="barostat != 'none'",
    ),
    Option(
        "barostat_frequency", "barostat", "Barostat interval", "int", 25,
        "How many steps between volume-change attempts. 25 is OpenMM's default and is "
        "frequent enough to equilibrate the density without the Monte Carlo moves "
        "costing meaningful time.",
        units="steps", minimum=1, maximum=1000,
        requires="barostat != 'none'", advanced=True,
    ),
    Option(
        "surface_tension", "barostat", "Surface tension", "number", 0.0,
        "Zero is the right choice for a normal bilayer: a tensionless membrane is what "
        "the CHARMM36 lipid parameters were designed to reproduce. Non-zero values are "
        "for deliberately stretching or compressing the membrane.",
        units="bar·nm", minimum=-100.0, maximum=100.0, step=1.0,
        requires="barostat == 'MonteCarloMembrane'",
    ),
    Option(
        "membrane_xymode", "barostat", "XY coupling", "select", "XYIsotropic",
        "Whether the two in-plane dimensions scale together. Isotropic is correct for "
        "a homogeneous bilayer.",
        choices=(
            Choice("XYIsotropic", "Isotropic in xy", "The two membrane-plane axes scale together."),
            Choice("XYAnisotropic", "Anisotropic in xy", "The axes scale independently."),
        ),
        requires="barostat == 'MonteCarloMembrane'", advanced=True,
    ),
    Option(
        "membrane_zmode", "barostat", "Z coupling", "select", "ZFree",
        "Whether the membrane normal is allowed to fluctuate. Free is standard.",
        choices=(
            Choice("ZFree", "Free", "The z dimension fluctuates independently."),
            Choice("ZFixed", "Fixed", "The z dimension is held constant."),
            Choice("ConstantVolume", "Constant volume", "z compensates for xy so total volume is conserved."),
        ),
        requires="barostat == 'MonteCarloMembrane'", advanced=True,
    ),
    Option(
        "anisotropic_pressure", "barostat", "Per-axis pressure", "text", "1.0, 1.0, 1.0",
        "Three pressures in bar for x, y and z. Equal values behave like the isotropic "
        "barostat but scale each axis independently, which lets a box become "
        "non-cubic in response to an anisotropic solute.",
        units="bar",
        requires="barostat == 'MonteCarloAnisotropic'", advanced=True,
    ),

    # -----------------------------------------------------------------------
    #  10. Restraints
    # -----------------------------------------------------------------------
    Option(
        "use_positional_restraints", "restraints", "Restrain the solute during equilibration", "checkbox", True,
        "Holds the solute near its starting coordinates while the freshly added water "
        "relaxes around it. Without this, a structure can distort during the first few "
        "picoseconds simply because the solvent has not yet settled. The restraints are "
        "released in stages before production.",
        openmm="openmm.CustomExternalForce",
    ),
    Option(
        "restraint_selection", "restraints", "Restrained atoms", "select", "protein_heavy",
        "Heavy atoms is the usual choice: it holds the fold without fighting hydrogen "
        "placement. Backbone-only allows side chains to relax, which is what you want "
        "if the side-chain packing came from a prediction.",
        choices=_ATOM_SELECTIONS[:5],
        requires="use_positional_restraints",
    ),
    Option(
        "restraint_force", "restraints", "Restraint force constant", "number", 1000.0,
        "1000 kJ/mol/nm² holds atoms within roughly 0.05 nm at 310 K: firm but not "
        "rigid. Much stiffer and minimisation struggles; much weaker and it is not "
        "really doing anything.",
        units="kJ/mol/nm²", minimum=0.0, maximum=100000.0, step=100.0,
        requires="use_positional_restraints",
    ),
    Option(
        "restraint_schedule", "restraints", "Release schedule", "text", "1000, 500, 100, 10, 0",
        "Force constants applied in equal slices of the equilibration, from the first "
        "value down to the last. A stepped release lets the structure adapt gradually "
        "rather than being let go all at once.",
        units="kJ/mol/nm²",
        requires="use_positional_restraints", advanced=True,
    ),
    Option(
        "custom_restraints", "restraints", "Distance restraints", "textarea", "",
        "One per line as SEL1 | SEL2 | distance_nm | force_constant, using MDTraj "
        "selection syntax. Applied throughout production, not just equilibration.",
        advanced=True,
        placeholder="resid 45 and name CA | resid 88 and name CA | 0.8 | 1000",
    ),

    # -----------------------------------------------------------------------
    #  11. Protocol
    # -----------------------------------------------------------------------
    Option(
        "minimize", "protocol", "Energy minimisation", "checkbox", True,
        "Removes the steric clashes that solvation and hydrogen placement inevitably "
        "create. Skipping it on a solvated system will usually blow the run up in the "
        "first few steps.",
        openmm="openmm.app.Simulation.minimizeEnergy",
    ),
    Option(
        "minimize_tolerance", "protocol", "Minimisation tolerance", "number", 10.0,
        "The force convergence criterion. 10 kJ/mol/nm is loose enough to be quick and "
        "tight enough that dynamics starts cleanly; there is no benefit in minimising "
        "hard, because the next step is to add thermal energy anyway.",
        units="kJ/mol/nm", minimum=0.1, maximum=1000.0, step=1.0,
        requires="minimize", advanced=True,
    ),
    Option(
        "minimize_max_iterations", "protocol", "Maximum minimisation steps", "int", 0,
        "Zero means run until the tolerance above is met, however long that takes. Set "
        "a ceiling only if a pathological starting structure is making minimisation "
        "run away rather than converge.",
        minimum=0, maximum=100000,
        requires="minimize", advanced=True,
    ),
    Option(
        "heat_duration", "protocol", "Heating", "number", 100.0,
        "Raises the temperature from near zero to the target in stages. A gradual ramp "
        "avoids the local hot spots that a sudden velocity assignment produces.",
        units="ps", minimum=0.0, maximum=10000.0, step=10.0,
    ),
    Option(
        "heat_start_temp", "protocol", "Starting temperature", "number", 50.0,
        "Where the ramp begins. Starting from a low but non-zero temperature is gentler "
        "than starting from absolute zero.",
        units="K", minimum=0.0, maximum=500.0, step=10.0,
        requires="heat_duration > 0", advanced=True,
    ),
    Option(
        "heat_stages", "protocol", "Heating stages", "int", 5,
        "How many discrete temperature steps the ramp uses.",
        minimum=1, maximum=100,
        requires="heat_duration > 0", advanced=True,
    ),
    Option(
        "equilibration_duration", "protocol", "Equilibration", "number", 1000.0,
        "Constant-pressure equilibration with the restraints releasing. 1 ns is enough "
        "for the density and the solvation shell to settle for a typical globular "
        "protein; membrane systems usually need considerably longer, because the lipids "
        "have to reorganise around the protein.",
        units="ps", minimum=0.0, maximum=100000.0, step=100.0,
    ),
    Option(
        "production_duration", "protocol", "Production", "number", 100.0,
        "The trajectory you will actually analyse. 100 ns is a reasonable first look at "
        "a globular protein's local dynamics; it is nowhere near enough to see a folding "
        "event or a large domain rearrangement.",
        units="ns", minimum=0.001, maximum=100000.0, step=1.0,
    ),
    Option(
        "replicates", "protocol", "Replicates", "int", 1,
        "Independent runs from the same starting structure with different velocity "
        "seeds. Three short replicates tell you far more about whether a result is real "
        "than one run three times as long, because they show you the variance.",
        minimum=1, maximum=20,
    ),
    Option(
        "seed_strategy", "protocol", "Replicate seeds", "select", "sequential",
        "How each replicate's random seed is chosen.",
        choices=(
            Choice("sequential", "Sequential from the base seed", "Reproducible: seed, seed+1, seed+2."),
            Choice("random", "Random per replicate", "Recorded in each replicate's manifest."),
        ),
        requires="replicates > 1", advanced=True,
    ),

    # -----------------------------------------------------------------------
    #  12. Output and reporting
    # -----------------------------------------------------------------------
    Option(
        "traj_format", "output", "Trajectory format", "select", "xtc",
        "XTC is compressed with a controlled precision loss and is typically three to "
        "four times smaller than DCD for coordinates that will only ever be analysed "
        "and viewed. DCD is uncompressed and universally readable.",
        choices=(
            Choice("xtc", "XTC", "Compressed. The default: much smaller, lossy only in the third decimal place."),
            Choice("dcd", "DCD", "Uncompressed, universally supported."),
            Choice("hdf5", "HDF5 (MDTraj)", "Carries topology and unit cell inside the file."),
        ),
        openmm="openmm.app.XTCReporter",
    ),
    Option(
        "traj_interval", "output", "Trajectory write interval", "number", 10.0,
        "How often a frame is saved. This is the main determinant of output size. "
        "10 ps gives 10,000 frames from a 100 ns run, which is ample for every metric "
        "on the analysis tab; sub-picosecond intervals are only needed for vibrational "
        "or fast-kinetics questions.",
        units="ps", minimum=0.001, maximum=10000.0, step=1.0,
    ),
    Option(
        "traj_selection", "output", "Atoms to save", "select", "protein_ligand",
        "Excluding water typically cuts the trajectory by 80-90%. Keep water only if "
        "you specifically need solvent structure, water-mediated contacts or diffusion.",
        choices=_ATOM_SELECTIONS,
    ),
    Option(
        "traj_custom_selection", "output", "Custom atom selection", "text", "",
        "An MDTraj selection expression, for example 'protein or resname LIG'.",
        requires="traj_selection == 'custom'",
        placeholder="protein or resname LIG",
    ),
    Option(
        "enforce_periodic_box", "output", "Wrap coordinates into the box", "checkbox", True,
        "Keeps every molecule inside the periodic box in the saved frames. Convenient "
        "for viewing, but it splits molecules that straddle a boundary, so the analysis "
        "step re-images the solute before computing anything geometric.",
        advanced=True,
    ),
    Option(
        "state_interval", "output", "Energy log interval", "number", 10.0,
        "How often energies, temperature, pressure and density are written to the CSV. "
        "These are tiny, so log them often: they are your only diagnostic if a run goes "
        "wrong.",
        units="ps", minimum=0.001, maximum=10000.0, step=1.0,
    ),
    Option(
        "state_fields", "output", "Energy log fields", "multiselect",
        ["step", "time", "potentialEnergy", "kineticEnergy", "totalEnergy",
         "temperature", "volume", "density", "speed", "progress", "remainingTime"],
        "What goes into the CSV. The defaults cover every convergence check the "
        "analysis tab plots.",
        choices=(
            Choice("step", "Step number"),
            Choice("time", "Simulated time"),
            Choice("potentialEnergy", "Potential energy"),
            Choice("kineticEnergy", "Kinetic energy"),
            Choice("totalEnergy", "Total energy", "The drift diagnostic: a steady rise means the timestep is too large."),
            Choice("temperature", "Temperature"),
            Choice("volume", "Box volume"),
            Choice("density", "Density", "Should plateau near 1.0 g/mL for a solvated system."),
            Choice("speed", "Speed (ns/day)"),
            Choice("progress", "Progress"),
            Choice("remainingTime", "Estimated time remaining"),
            Choice("elapsedTime", "Elapsed wall time"),
        ),
        openmm="openmm.app.StateDataReporter", advanced=True,
    ),
    Option(
        "checkpoint_interval", "output", "Checkpoint interval", "number", 100.0,
        "How often a restart file is written. A long run will be interrupted at some "
        "point, and a checkpoint is the difference between resuming and starting again.",
        units="ps", minimum=1.0, maximum=100000.0, step=10.0,
        openmm="openmm.app.CheckpointReporter",
    ),
    Option(
        "save_system_xml", "output", "Save serialised system", "checkbox", True,
        "Writes the System, Integrator and final State as XML. This is what makes a run "
        "genuinely reproducible and lets you extend it later on different hardware "
        "without rebuilding anything.",
        advanced=True,
    ),
    Option(
        "save_solvated_pdb", "output", "Save the solvated structure", "checkbox", True,
        "The complete prepared system before dynamics. Essential for diagnosing a "
        "preparation problem after the fact.",
        advanced=True,
    ),

    # -----------------------------------------------------------------------
    #  13. Platform and performance
    # -----------------------------------------------------------------------
    Option(
        "platform", "platform", "Compute platform", "select", "auto",
        "OpenMM has no official Metal back end. On Apple Silicon the realistic choices "
        "are the multithreaded CPU platform and OpenCL, which does reach the GPU but is "
        "deprecated by Apple. Automatic benchmarking runs a few hundred steps on each "
        "available platform and picks the fastest, which takes under a minute and is "
        "more reliable than guessing.",
        choices=(
            Choice("auto", "Benchmark and choose", "The recommended default."),
            Choice("CPU", "CPU", "Multithreaded. Consistently fast on Apple Silicon's wide cores."),
            Choice("OpenCL", "OpenCL", "Reaches the integrated GPU. Deprecated by Apple but functional on M-series."),
            Choice("Metal", "Metal", "The third-party openmm-metal plugin. Not installed by default and not officially supported.",
                   experimental=True),
            Choice("CUDA", "CUDA", "NVIDIA only. Present for completeness if you run the bundle on a Linux box.",
                   experimental=True),
            Choice("Reference", "Reference", "The unoptimised reference implementation. For debugging only; roughly 100× slower.",
                   experimental=True),
        ),
        openmm="openmm.Platform.getPlatformByName",
    ),
    Option(
        "precision", "platform", "Precision", "select", "mixed",
        "Mixed precision computes forces in single precision and accumulates in double. "
        "It is as accurate as full double precision for everything that matters here "
        "and considerably faster. Full double is only needed for energy-conservation "
        "tests.",
        choices=(
            Choice("mixed", "Mixed", "The right choice. Single-precision forces, double-precision accumulation."),
            Choice("single", "Single", "Fastest, and prone to visible energy drift in long runs."),
            Choice("double", "Double", "Slowest. For NVE energy-conservation testing."),
        ),
        requires="platform in ('auto', 'OpenCL', 'CUDA', 'Metal')",
    ),
    Option(
        "cpu_threads", "platform", "CPU threads", "int", 0,
        "Zero lets OpenMM use every core. On an Apple Silicon machine with efficiency "
        "cores it is often faster to set this to the number of performance cores only, "
        "because the slowest thread paces the whole step.",
        minimum=0, maximum=256,
        requires="platform in ('auto', 'CPU')", advanced=True,
    ),
    Option(
        "deterministic_forces", "platform", "Deterministic forces", "checkbox", False,
        "Forces identical results from identical input at a small performance cost. "
        "Needed only when you require bitwise reproducibility.",
        advanced=True,
    ),
    Option(
        "device_index", "platform", "Device index", "text", "",
        "Which GPU to use, when there is more than one. Blank means the first.",
        requires="platform in ('OpenCL', 'CUDA', 'Metal')", advanced=True,
    ),
    Option(
        "memory_guard", "platform", "Refuse to start if memory is short", "checkbox", True,
        "Estimates the system's memory footprint and stops before starting if it would "
        "push the machine into swap. Swapping during MD does not merely slow the run "
        "down, it can extend a job from hours to days.",
        advanced=True,
    ),

    # -----------------------------------------------------------------------
    #  14. Analysis
    # -----------------------------------------------------------------------
    Option(
        "analysis_metrics", "analysis", "Metrics to compute", "multiselect",
        ["rmsd", "rmsf", "rgyr", "sasa", "dssp", "hbonds", "contacts", "pca", "clusters"],
        "Computed on your machine against the full trajectory, then packed into the "
        "results file. Doing it here rather than server-side is what keeps the upload "
        "small and means nothing is lost to decimation.",
        choices=(
            Choice("rmsd", "RMSD", "Backbone deviation from the reference over time. The first thing to look at."),
            Choice("rmsf", "RMSF", "Per-residue fluctuation. Shows you which loops move."),
            Choice("rgyr", "Radius of gyration", "Global compactness. Detects unfolding or collapse."),
            Choice("sasa", "Solvent-accessible surface area", "Shrake-Rupley. Slow on long trajectories."),
            Choice("dssp", "Secondary structure", "Per-residue DSSP assignment over time."),
            Choice("hbonds", "Hydrogen bonds", "Baker-Hubbard occupancy analysis."),
            Choice("saltbridges", "Salt bridges", "Charged-pair contact occupancy."),
            Choice("contacts", "Contact map", "Residue-residue minimum distances."),
            Choice("nativecontacts", "Native contacts (Q)", "Fraction of the starting structure's contacts retained."),
            Choice("pca", "Principal component analysis", "Essential dynamics: the dominant collective motions."),
            Choice("clusters", "Conformational clustering", "RMSD clustering with representative frames."),
            Choice("ligand_rmsd", "Ligand RMSD", "Ligand pose stability after aligning on the protein.",
                   requires="has_ligands"),
            Choice("ligand_contacts", "Ligand contacts", "Per-residue contact occupancy with the ligand.",
                   requires="has_ligands"),
            Choice("membrane_scd", "Lipid order parameters", "Deuterium order parameters per tail carbon.",
                   requires="use_membrane"),
            Choice("membrane_apl", "Area per lipid", "The standard bilayer equilibration diagnostic.",
                   requires="use_membrane"),
            Choice("membrane_thickness", "Bilayer thickness", "Phosphate-to-phosphate distance.",
                   requires="use_membrane"),
        ),
    ),
    Option(
        "analysis_reference", "analysis", "Reference structure", "select", "first",
        "What RMSD and native contacts are measured against. The minimised structure is "
        "usually the more meaningful reference, because the first production frame has "
        "already drifted through equilibration.",
        choices=(
            Choice("first", "First production frame"),
            Choice("minimized", "Minimised structure", "The recommended reference."),
            Choice("input", "Input structure", "The original deposited coordinates, before any repair."),
            Choice("average", "Trajectory average", "Useful for RMSF, misleading for RMSD."),
        ),
        advanced=True,
    ),
    Option(
        "analysis_selection", "analysis", "Alignment selection", "select", "backbone",
        "Which atoms are superposed before geometric metrics are computed. Backbone is "
        "standard; align on a single domain if the molecule has two that move relative "
        "to each other, or the RMSD will just report the hinge.",
        choices=_ATOM_SELECTIONS[:5],
        advanced=True,
    ),
    Option(
        "custom_trackers", "analysis", "Custom measurements", "textarea", "",
        "One per line as NAME | TYPE | selections, where TYPE is distance, angle or "
        "dihedral. Selections use MDTraj syntax and are separated by semicolons.",
        advanced=True,
        placeholder="gate_width | distance | resid 45 and name CA; resid 88 and name CA",
    ),
    Option(
        "payload_tier", "analysis", "Results file size", "select", "standard",
        "How much trajectory travels back with the metrics. Standard is right for "
        "almost everyone. Full is only needed if you expect to ask the web app for "
        "extra analyses later; the full trajectory always stays on your machine "
        "regardless of what you choose here.",
        choices=(
            Choice("light", "Light (~5 MB)", "Metrics plus 100 frames. Enough to view and plot."),
            Choice("standard", "Standard (~25 MB)", "Metrics plus 500 protein-only frames. The recommended default."),
            Choice("full", "Full (~150 MB)", "Metrics plus 2000 frames. Enables server-side re-analysis."),
        ),
    ),
)


# ===========================================================================
#  Derived lookups
# ===========================================================================

BY_ID: dict[str, Option] = {o.id: o for o in OPTIONS}

BY_GROUP: dict[str, tuple[Option, ...]] = {
    g.id: tuple(o for o in OPTIONS if o.group == g.id) for g in GROUPS
}

GROUP_BY_ID: dict[str, Group] = {g.id: g for g in GROUPS}


def defaults() -> dict[str, Any]:
    """A complete config with every option at its default value."""
    return {o.id: (list(o.default) if isinstance(o.default, list) else o.default) for o in OPTIONS}


def option(option_id: str) -> Option:
    """Look up one option, with a useful error rather than a KeyError."""
    try:
        return BY_ID[option_id]
    except KeyError:
        raise KeyError(
            f"unknown option '{option_id}' -- options are declared only in "
            f"flexappeal/options.py; did you mean one of "
            f"{sorted(i for i in BY_ID if i.startswith(option_id[:4]))}?"
        ) from None
