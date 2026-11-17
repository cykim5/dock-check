"""Turn SMILES into docking-ready 3D structures.

Three steps people skip, each of which quietly ruins a decoy set:

  largest fragment   ChEMBL stores salts, so "drug" often means drug
                     plus a sodium. Docking the sodium is not useful.
  explicit hydrogens Geometry and hydrogen bonding both need them.
  one conformer      Enough for a decoy. A native ligand you care about
                     deserves more, which is what n_confs is for.
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.rdBase import DisableLog

from meeko import MoleculePreparation, PDBQTWriterLegacy

DisableLog("rdApp.*")


def largest_fragment(mol):
    frags = Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=True)
    if len(frags) == 1:
        return mol
    return max(frags, key=lambda m: m.GetNumHeavyAtoms())


def embed(smiles, seed=42, n_confs=1):
    """SMILES to a 3D molecule, or (None, reason)."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None, "RDKit could not parse the SMILES"
    mol = largest_fragment(mol)
    mol = Chem.AddHs(mol)
    params = AllChem.ETKDGv3()
    params.randomSeed = seed
    if n_confs == 1:
        if AllChem.EmbedMolecule(mol, params) != 0:
            params.useRandomCoords = True
            if AllChem.EmbedMolecule(mol, params) != 0:
                return None, "could not generate a 3D conformer"
    else:
        if len(AllChem.EmbedMultipleConfs(mol, n_confs, params)) == 0:
            return None, "could not generate 3D conformers"
    try:
        if AllChem.MMFFHasAllMoleculeParams(mol):
            AllChem.MMFFOptimizeMoleculeConfs(mol, maxIters=500)
        else:
            AllChem.UFFOptimizeMoleculeConfs(mol, maxIters=500)
    except Exception:
        pass
    return mol, ""


def to_pdbqt(mol, path):
    """Write a Meeko PDBQT. Returns (ok, reason)."""
    try:
        setups = MoleculePreparation().prepare(mol)
    except Exception as err:
        return False, f"Meeko failed: {str(err)[:120]}"
    if len(setups) == 0:
        return False, "Meeko produced no setup"
    text, ok, err = PDBQTWriterLegacy.write_string(setups[0])
    if ok == False:
        return False, f"Meeko could not write: {str(err)[:120]}"
    with open(path, "w") as f:
        f.write(text)
    return True, ""


def heavy_atoms(mol):
    return mol.GetNumHeavyAtoms()
