"""Build a decoy set matched to the native ligand's size.

Size matching is the part people get wrong. Empirical scoring functions
correlate strongly with how many heavy atoms a molecule has, so a decoy
set of small molecules makes any large native ligand look like a
winner, and the test passes without meaning anything. Decoys here are
drawn from the same heavy atom window as the native ligand, so size
cannot be what separates them.

The pool is approved small molecules from ChEMBL: real, synthesizable,
and not chosen to be easy.
"""
import csv
import json
import time
import urllib.parse
import urllib.request

API = "https://www.ebi.ac.uk/chembl/api/data/molecule.json"
PAGE = 1000

# ChEMBL reports full_mwt for the formulated salt, counterion included,
# while docking only ever sees the parent. So the fetch asks for a
# window this much wider at the top and the real window is applied to
# the parent weight afterwards. Without the allowance a 380 Da drug
# sold as a 476 Da mesylate is dropped for being too heavy, when the
# molecule that would be docked is well inside the window.
SALT_ALLOWANCE = 200.0


def _get(url):
    req = urllib.request.Request(
        url, headers={"Accept": "application/json",
                      "User-Agent": "dock-check"})
    for attempt in range(4):
        try:
            with urllib.request.urlopen(req, timeout=120) as r:
                return json.loads(r.read().decode())
        except Exception as err:
            if attempt == 3:
                raise
            print(f"  retry {attempt + 1} after {err}", flush=True)
            time.sleep(5)


def fetch_pool(mw_min=150.0, mw_max=600.0, max_phase=4, limit=None):
    """Every approved small molecule in a molecular weight window."""
    rows, seen, offset = [], set(), 0
    while True:
        params = {"max_phase": max_phase, "molecule_type": "Small molecule",
                  "limit": PAGE, "offset": offset}
        data = _get(API + "?" + urllib.parse.urlencode(params))
        page = data.get("molecules", [])
        if len(page) == 0:
            break
        for m in page:
            cid = m.get("molecule_chembl_id")
            struct = m.get("molecule_structures") or {}
            smiles = struct.get("canonical_smiles")
            props = m.get("molecule_properties") or {}
            mw = props.get("full_mwt")
            if cid is None or smiles is None or mw is None:
                continue
            if cid in seen:
                continue
            mw = float(mw)
            if mw < mw_min or mw > mw_max:
                continue
            seen.add(cid)
            rows.append({"id": cid, "smiles": smiles, "mw": mw,
                         "name": m.get("pref_name") or ""})
        offset = offset + PAGE
        print(f"  fetched {len(rows)} so far", flush=True)
        if limit is not None and len(rows) >= limit:
            break
        if data.get("page_meta", {}).get("next") is None:
            break
    return rows


def heavy_atom_count(sdf_path):
    """Heavy atoms in a V2000 SDF, without needing RDKit."""
    lines = open(sdf_path).read().splitlines()
    n = int(lines[3][0:3])
    return sum(1 for l in lines[4:4 + n] if l[31:34].strip() != "H")


HYDROGEN_TYPES = ("H", "HD", "HS")


def heavy_atom_count_pdbqt(path):
    """Heavy atoms in a PDBQT, by Vina atom type rather than element.

    The type is the last field on the line. Columns 78 to 79 are the
    standard position, but files written by different tools drift, so
    the last field is read first and the fixed columns are the fallback.
    """
    n = 0
    for l in open(path):
        if l[:6].strip() not in ("ATOM", "HETATM"):
            continue
        fields = l.split()
        t = fields[-1].upper() if len(fields) > 0 else ""
        if len(t) > 2 or t == "":
            t = l[77:79].strip().upper()
        if t in HYDROGEN_TYPES:
            continue
        n = n + 1
    return n


def heavy_atoms_of(path):
    """Heavy atom count from whichever format this file is."""
    if path.endswith(".pdbqt") == True:
        return heavy_atom_count_pdbqt(path)
    return heavy_atom_count(path)


def size_window(native_heavy, tolerance=6):
    """The heavy atom window decoys must fall in."""
    return (native_heavy - tolerance, native_heavy + tolerance)


def write_pool(rows, path):
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["id", "smiles", "mw", "name"])
        w.writeheader()
        w.writerows(rows)
    return len(rows)
