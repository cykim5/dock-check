"""Is the site actually open?

A control that fails tells you the screen is not trustworthy. It does
not tell you which part was wrong, and the part people reach for last
is the receptor, because a structure feels like a fact rather than a
choice. It is a choice. A structure solved without the ligand bound has
had nothing holding its site open, and side chains relax into empty
space. Dock into that and every ligand is competing with the protein
for the same volume.

This measures that directly: put the ligand where the crystal says it
goes, and count how much of the receptor is already in the way.

The distinction that matters is backbone against side chain. A side
chain in the way is a bad starting conformation, and flexible side
chain docking can sometimes recover it. Backbone in the way is not
recoverable by any docking setting, because nothing moves a backbone
except a different structure.
"""

import math

# van der Waals radii in angstroms, Bondi 1964 with the usual Rowland
# and Taylor revision for hydrogen. Only the elements a ligand or a
# protein actually presents are listed; anything else falls back to
# carbon, which is the most common and middle of the range.
VDW = {
    "C": 1.70, "N": 1.55, "O": 1.52, "S": 1.80, "P": 1.80,
    "F": 1.47, "CL": 1.75, "BR": 1.85, "I": 1.98, "SE": 1.90,
    "B": 1.92, "SI": 2.10,
}
VDW_DEFAULT = 1.70

# Two atoms are allowed to sit closer than their radii suggest, because
# a contact is not a clash. 0.5 A is the usual slack for ordinary
# contacts. Nitrogen and oxygen get another 0.6 A on top, because a
# hydrogen bond legitimately pulls two polar atoms to around 2.7 A and
# flagging every hydrogen bond as a clash would make the check useless
# on exactly the sites it is meant for.
TOLERANCE = 0.5
HBOND_EXTRA = 0.6

# Nothing further apart than this can overlap, so it is the cutoff for
# bothering to measure at all.
NEAR = 5.0
POLAR = ("N", "O")

BACKBONE_ATOMS = ("N", "CA", "C", "O", "OXT")

# Chemically modified amino acids are written as HETATM but they are
# part of the chain, not something bound to it. Selenomethionine is the
# common one, because it is how the phase problem gets solved. Treat
# them as receptor: leave them out and the protein has holes in it
# exactly where a clash would have been found.
MODIFIED_RESIDUES = {
    "MSE", "SEP", "TPO", "PTR", "CSO", "CME", "CSD", "OCS", "KCX",
    "LLP", "PCA", "HIC", "MLY", "M3L", "ALY", "SAC", "FME", "CAS",
    "CSS", "CSX", "NEP", "HYP", "DAL", "SAH", "ABA", "CGU", "TYS",
}

# HETATM records that are not the thing anyone means by "the ligand":
# water, the salts and buffers that come with a crystallisation
# condition, and the cryoprotectants.
NOT_A_LIGAND = {
    "HOH", "DOD", "WAT", "SO4", "PO4", "CL", "NA", "K", "MG", "CA",
    "ZN", "MN", "FE", "CU", "NI", "CD", "HG", "IOD", "BR", "ACT",
    "EDO", "GOL", "PEG", "PG4", "PGE", "MPD", "DMS", "TRS", "EPE",
    "MES", "IMD", "FMT", "ACY", "NO3", "CO3", "AZI", "BME", "UNX",
}


def _element(line):
    """The element for a PDB ATOM or HETATM line.

    Columns 77 to 78 hold it when the file is well formed. Plenty of
    files in the wild leave that blank, so the atom name is the
    fallback: its first character is the element for everything except
    the two letter metals and halogens, which the explicit column
    almost always carries anyway.
    """
    el = line[76:78].strip().upper()
    if el == "":
        name = line[12:16].strip()
        el = name[0:1].upper()
    return el


def read_atoms(path, want_hetatm=False, resname=None, chain=None, resseq=None):
    """Heavy atoms from a PDB file.

    Hydrogens are dropped. Most crystal structures do not have them,
    the ones that do were given them by a program rather than by the
    data, and a clash is decided by the heavy atoms either way.

    Only the first model is read, so an NMR ensemble gives its first
    conformer rather than all of them piled on top of each other.
    """
    atoms = []
    in_model = True
    for line in open(path):
        if line.startswith("MODEL"):
            # A second MODEL means the first one has ended.
            in_model = len(atoms) == 0
            continue
        if line.startswith("ENDMDL"):
            in_model = False
            continue
        if in_model == False:
            continue
        is_atom = line.startswith("ATOM")
        is_het = line.startswith("HETATM")
        if is_atom == False and is_het == False:
            continue
        residue = line[17:20].strip().upper()
        # A modified amino acid is chain, so it answers to the receptor
        # read and never to the ligand read, whichever record type it
        # happens to be written as. This has to be decided before the
        # record type is used to filter, not after.
        if residue in MODIFIED_RESIDUES:
            is_het = False
        if want_hetatm == True and is_het == False:
            continue
        if want_hetatm == False and is_het == True:
            continue
        if resname is not None and residue != resname.upper():
            continue
        this_chain = line[21:22].strip()
        if chain is not None and this_chain != chain:
            continue
        this_resseq = line[22:27].strip()
        if resseq is not None and this_resseq != str(resseq).strip():
            continue
        element = _element(line)
        if element == "H" or element == "D":
            continue
        # An alternate location flag of anything but blank or A is a
        # second copy of the same atom, and counting both would
        # double every contact it makes.
        altloc = line[16:17].strip()
        if altloc != "" and altloc != "A":
            continue
        atoms.append({
            "name": line[12:16].strip(),
            "element": element,
            "residue": residue,
            "chain": this_chain,
            "resseq": this_resseq,
            "xyz": (float(line[30:38]), float(line[38:46]), float(line[46:54])),
        })
    return atoms


def bound_ligands(path):
    """What is bound in this structure, ignoring solvent and buffer.

    An empty list is the warning this whole module exists for: a
    structure with nothing in it was solved with nothing holding its
    sites open.
    """
    found = {}
    for atom in read_atoms(path, want_hetatm=True):
        if atom["residue"] in NOT_A_LIGAND:
            continue
        if atom["residue"] in MODIFIED_RESIDUES:
            continue
        key = (atom["residue"], atom["chain"], atom["resseq"])
        found[key] = found.get(key, 0) + 1
    ligands = []
    for (residue, chain, resseq), count in found.items():
        ligands.append({"residue": residue, "chain": chain,
                        "resseq": resseq, "heavy_atoms": count})
    ligands.sort(key=lambda l: -l["heavy_atoms"])
    return ligands


def _allowance(element_a, element_b):
    """How far inside their radii two atoms may sit before it counts."""
    slack = TOLERANCE
    if element_a in POLAR and element_b in POLAR:
        slack = slack + HBOND_EXTRA
    return slack


def occlusion(receptor_atoms, ligand_atoms):
    """How much of the receptor is standing in the ligand's volume.

    Every ligand atom is compared against every receptor atom, which is
    the honest way to do it and fast enough: a site ligand is tens of
    atoms and a chain is a few thousand, so this is a few hundred
    thousand distances, well under a second. Growing a neighbour grid
    for it would be speed nobody needs at the cost of a bug nobody
    sees.
    """
    clashes = []
    worst_contact = None
    for lig in ligand_atoms:
        lx, ly, lz = lig["xyz"]
        lig_radius = VDW.get(lig["element"], VDW_DEFAULT)
        worst_here = None
        for rec in receptor_atoms:
            rx, ry, rz = rec["xyz"]
            # A cheap box test first. Most receptor atoms are nowhere
            # near, and skipping the square root for those is most of
            # the work avoided.
            if abs(lx - rx) > NEAR or abs(ly - ry) > NEAR or abs(lz - rz) > NEAR:
                continue
            distance = math.dist((lx, ly, lz), (rx, ry, rz))
            if worst_contact is None or distance < worst_contact:
                worst_contact = distance
            rec_radius = VDW.get(rec["element"], VDW_DEFAULT)
            limit = lig_radius + rec_radius - _allowance(lig["element"], rec["element"])
            overlap = limit - distance
            if overlap <= 0:
                continue
            is_backbone = rec["name"].upper() in BACKBONE_ATOMS
            record = {
                "ligand_atom": lig["name"],
                "receptor_atom": rec["name"],
                "residue": rec["residue"],
                "chain": rec["chain"],
                "resseq": rec["resseq"],
                "distance": round(distance, 3),
                "overlap": round(overlap, 3),
                "backbone": is_backbone,
            }
            if worst_here is None or overlap > worst_here["overlap"]:
                worst_here = record
        if worst_here is not None:
            clashes.append(worst_here)

    backbone = [c for c in clashes if c["backbone"] == True]
    result = {
        "ligand_atoms": len(ligand_atoms),
        "clashing_atoms": len(clashes),
        "backbone_clashes": len(backbone),
        "worst_contact": round(worst_contact, 3) if worst_contact is not None else None,
        "worst_overlap": max([c["overlap"] for c in clashes], default=0.0),
        "worst_backbone_overlap": max([c["overlap"] for c in backbone], default=0.0),
        "residues": sorted({"%s %s%s" % (c["residue"], c["chain"], c["resseq"])
                            for c in clashes}),
        "clashes": sorted(clashes, key=lambda c: -c["overlap"]),
    }
    result["verdict"] = verdict(result)
    return result


def verdict(result):
    """One of three words, because there are three different actions.

    clear: dock here.
    side-chain: the site is occupied by something that can move, so
        flexible side chain docking is worth a try before giving up on
        the structure.
    backbone: nothing in a docking program moves a backbone. Find a
        structure solved with something in the site.
    """
    if result["clashing_atoms"] == 0:
        return "clear"
    if result["backbone_clashes"] > 0:
        return "backbone"
    return "side-chain"


def summary(result):
    """The verdict as a sentence, for the report and the terminal."""
    clashing = result["clashing_atoms"]
    total = result["ligand_atoms"]
    if result["verdict"] == "clear":
        if result["worst_contact"] is None:
            return ("the site is open: no receptor atom comes within %.0f A of "
                    "any of the %d ligand atoms, which is far enough apart to "
                    "be worth checking the two files share a coordinate frame"
                    % (NEAR, total))
        return ("the site is open: none of the %d ligand atoms clash with the "
                "receptor, closest contact %.2f A" % (total, result["worst_contact"]))
    where = ", ".join(result["residues"][:4])
    if len(result["residues"]) > 4:
        where = where + " and %d more" % (len(result["residues"]) - 4)
    if result["verdict"] == "backbone":
        return ("the site is occluded: %d of %d ligand atoms clash with the "
                "receptor, worst by %.2f A, and %d of those clashes are with "
                "backbone (worst %.2f A), which no flexible side chain setting "
                "can move. Dock against a structure solved with something in "
                "this site. Clashing residues: %s"
                % (clashing, total, result["worst_overlap"],
                   result["backbone_clashes"], result["worst_backbone_overlap"], where))
    return ("the site is partly occupied: %d of %d ligand atoms clash with "
            "side chains, worst by %.2f A, none with backbone. Flexible side "
            "chain docking may recover this. Clashing residues: %s"
            % (clashing, total, result["worst_overlap"], where))


def check(receptor_path, ligand_path, ligand_resname=None, ligand_chain=None,
          ligand_resseq=None, receptor_chain=None):
    """Run the check and return everything worth reporting.

    The ligand may be its own file or a residue inside the receptor
    file. The second case is the common one and the easier one to get
    right, because nothing has been superposed.
    """
    same_file = receptor_path == ligand_path
    ligand_atoms = read_atoms(ligand_path, want_hetatm=True,
                              resname=ligand_resname, chain=ligand_chain,
                              resseq=ligand_resseq)
    if len(ligand_atoms) == 0:
        # A ligand given as its own file is often written as ATOM
        # records rather than HETATM, so fall back before giving up.
        ligand_atoms = read_atoms(ligand_path, want_hetatm=False,
                                  resname=ligand_resname, chain=ligand_chain,
                                  resseq=ligand_resseq)
    if len(ligand_atoms) == 0:
        raise ValueError("no ligand atoms found in %s" % ligand_path)

    receptor_atoms = read_atoms(receptor_path, want_hetatm=False,
                                chain=receptor_chain)
    if len(receptor_atoms) == 0:
        raise ValueError("no receptor atoms found in %s" % receptor_path)

    result = occlusion(receptor_atoms, ligand_atoms)
    result["receptor"] = receptor_path
    result["ligand"] = ligand_path
    result["ligand_from_receptor_file"] = same_file
    result["bound_ligands"] = bound_ligands(receptor_path)
    result["summary"] = summary(result)
    return result
