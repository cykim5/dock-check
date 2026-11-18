"""One definition of what the molecule actually is.

A drug library is a list of formulated products, not a list of
molecules. ChEMBL stores paroxetine hydrochloride and paroxetine
mesylate under different identifiers, and both are the same molecule
once the counterion comes off. Two things go wrong if you take the
stored SMILES at face value.

  the set is smaller than it says   Two salts of one parent become two
                                    decoys. The percentile is then
                                    computed over a set with repeats in
                                    it, and a repeated parent counts
                                    twice wherever it happens to score.

  the size is the wrong molecule's  The stored molecular weight belongs
                                    to the salt, counterion included.
                                    Filter on it and you are selecting
                                    on a weight the docking never sees,
                                    because the counterion is stripped
                                    before anything is docked.

The counterion is not an inert label either. It carries the charge and
the polar groups that a scoring function reads as chemistry, so a
mesylate can make its parent look like it owns a sulfonate. Judge the
parent and the question goes away.

One limit worth stating. An entry that is nothing but an inorganic salt
has no drug in it to find, so the parent of zinc sulfate comes back as
sulfuric acid and the parent of sodium chloride as hydrogen chloride.
Those are not useful answers, but they are small, and the heavy atom
window a decoy set is matched on removes them without being asked to.
"""
from rdkit import Chem
from rdkit.Chem import Descriptors, inchi
from rdkit.Chem.MolStandardize import rdMolStandardize
from rdkit.rdBase import DisableLog

DisableLog("rdApp.*")


def parent_mol(smiles):
    """The molecule left after the counterions come off, or None.

    FragmentParent hands back a molecule with no ring information, so
    anything that asks about rings afterwards raises until it has been
    sanitised again. Sanitising here means callers never meet that.
    """
    if smiles is None or smiles == "":
        return None
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    try:
        mol = rdMolStandardize.FragmentParent(mol)
        Chem.SanitizeMol(mol)
    except Exception:
        return None
    if mol.GetNumHeavyAtoms() == 0:
        return None
    return mol


def parent_key(smiles):
    """An identity for the parent, so salt forms collapse together.

    The InChIKey is used rather than canonical SMILES because it is
    insensitive to how the input was written.
    """
    mol = parent_mol(smiles)
    if mol is None:
        return None
    try:
        return inchi.MolToInchiKey(mol)
    except Exception:
        return None


def parent_weight(smiles):
    """Molecular weight of the parent, not of the formulated salt."""
    mol = parent_mol(smiles)
    if mol is None:
        return None
    return Descriptors.MolWt(mol)


def stripped_fragments(smiles):
    """The pieces that came off, largest first, as SMILES.

    Reported rather than discarded, so a run can say what it removed
    instead of asking anyone to take it on trust.
    """
    if smiles is None or smiles == "":
        return []
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return []
    frags = Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=False)
    if len(frags) < 2:
        return []
    biggest = max(range(len(frags)),
                  key=lambda i: frags[i].GetNumHeavyAtoms())
    rest = [f for i, f in enumerate(frags) if i != biggest]
    rest.sort(key=lambda m: m.GetNumHeavyAtoms(), reverse=True)
    out = []
    for f in rest:
        try:
            out.append(Chem.MolToSmiles(f))
        except Exception:
            continue
    return out


def is_salt(smiles):
    """True when the stored structure is more than its parent."""
    return len(stripped_fragments(smiles)) > 0


def collapse_to_parents(rows, smiles_key="smiles", id_key="id"):
    """Drop rows whose parent already appeared.

    Returns (kept, collapsed, unparsable).

      kept        one row per distinct parent, in the order given, each
                  with parent_key, parent_mw and stripped added
      collapsed   the dropped rows, each with the id it duplicates
      unparsable  rows whose parent could not be worked out at all

    The first row to carry a parent wins, so shuffling the pool before
    this is called is what decides which salt form represents it.
    """
    kept, collapsed, unparsable = [], [], []
    first_seen = {}
    for r in rows:
        smiles = r.get(smiles_key)
        key = parent_key(smiles)
        if key is None:
            unparsable.append(dict(r))
            continue
        if key in first_seen:
            dup = dict(r)
            dup["duplicate_of"] = first_seen[key]
            dup["parent_key"] = key
            collapsed.append(dup)
            continue
        first_seen[key] = r.get(id_key)
        row = dict(r)
        row["parent_key"] = key
        row["parent_mw"] = parent_weight(smiles)
        row["stripped"] = stripped_fragments(smiles)
        kept.append(row)
    return kept, collapsed, unparsable
