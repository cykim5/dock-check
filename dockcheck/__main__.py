"""Command line entry point.

  python -m dockcheck decoys --native ligand.sdf --out decoys/ --n 200
  python -m dockcheck run --receptor r.pdbqt --native ligand.sdf \
      --decoys decoys/ --center 62.3 16.4 34.2 --size 28 --out report/
  python -m dockcheck receptor --receptor 2b69.pdb --resname UDP --resseq 900
"""
import argparse
import csv
import json
import os
import random
import sys

from . import decoys as dk
from . import receptor as rc
from .control import run_control
from .report import PASS_PERCENTILE, write_report


def cmd_decoys(a):
    heavy = dk.heavy_atom_count(a.native) if a.heavy is None else a.heavy
    lo, hi = dk.size_window(heavy, a.tolerance)
    print(f"native has {heavy} heavy atoms; decoys must have {lo} to {hi}")
    os.makedirs(a.out, exist_ok=True)
    pool_csv = os.path.join(a.out, "pool.csv")
    if os.path.exists(pool_csv) == True and a.refetch == False:
        rows = list(csv.DictReader(open(pool_csv)))
        print(f"reusing {len(rows)} molecules from {pool_csv}")
    else:
        print("fetching approved small molecules from ChEMBL ...")
        rows = dk.fetch_pool(a.mw_min, a.mw_max + dk.SALT_ALLOWANCE)
        dk.write_pool(rows, pool_csv)
        print(f"wrote {len(rows)} to {pool_csv}")

    from .prepare import embed, heavy_atoms, to_pdbqt
    from .salts import collapse_to_parents
    random.Random(a.seed).shuffle(rows)

    # One row per molecule, not one per formulated product. Done after
    # the shuffle, so which salt form represents a parent is decided by
    # the seed rather than by ChEMBL's ordering.
    rows, collapsed, unparsable = collapse_to_parents(rows)
    print(f"{len(rows)} distinct molecules; dropped {len(collapsed)} salt "
          f"forms of a parent already in the set")
    if len(unparsable) > 0:
        print(f"  {len(unparsable)} could not be reduced to a parent")

    # The window the user asked for, applied to the molecule that will
    # actually be docked rather than to the salt's weight.
    before = len(rows)
    rows = [r for r in rows if r["parent_mw"] is not None
            and a.mw_min <= r["parent_mw"] <= a.mw_max]
    if before != len(rows):
        print(f"  {before - len(rows)} fell outside {a.mw_min:.0f} to "
              f"{a.mw_max:.0f} Da once the counterion came off")

    # Counted on the final set, so these numbers describe the same
    # molecules the decoys are drawn from.
    salt_rows = [r for r in rows if len(r["stripped"]) > 0]

    lig_dir = os.path.join(a.out, "ligands")
    os.makedirs(lig_dir, exist_ok=True)
    kept, skipped = [], 0
    for r in rows:
        if len(kept) >= a.n:
            break
        path = os.path.join(lig_dir, r["id"] + ".pdbqt")
        if os.path.exists(path) == True:
            kept.append({"id": r["id"], "pdbqt": path, "heavy": None})
            continue
        mol, why = embed(r["smiles"], seed=a.seed)
        if mol is None:
            skipped = skipped + 1
            continue
        h = heavy_atoms(mol)
        if h < lo or h > hi:
            continue
        ok, why = to_pdbqt(mol, path)
        if ok == False:
            skipped = skipped + 1
            continue
        kept.append({"id": r["id"], "pdbqt": path, "heavy": h})
        if len(kept) % 25 == 0:
            print(f"  prepared {len(kept)}/{a.n}", flush=True)
    with open(os.path.join(a.out, "decoys.json"), "w") as f:
        json.dump({"native_heavy": heavy, "window": [lo, hi],
                   "decoys": kept,
                   "salts": {
                       "distinct_parents": len(rows),
                       "collapsed": len(collapsed),
                       "unparsable": len(unparsable),
                       "drawn_from_a_salt": len(salt_rows),
                       "examples": [
                           {"id": c.get("id"),
                            "duplicate_of": c.get("duplicate_of")}
                           for c in collapsed[:12]],
                   }}, f, indent=2)
    print(f"prepared {len(kept)} decoys, skipped {skipped}")


def cmd_receptor(a):
    try:
        result = rc.check(a.receptor, a.ligand or a.receptor,
                          ligand_resname=a.resname, ligand_chain=a.ligand_chain,
                          ligand_resseq=a.resseq, receptor_chain=a.chain)
    except ValueError as why:
        sys.exit(str(why))

    print(f"# dock-check receptor: {result['verdict'].upper()}")
    print()
    print(result["summary"])
    print()
    if len(result["bound_ligands"]) == 0:
        print("nothing is bound in this structure apart from solvent, which is "
              "the condition this check exists for: a site with nothing in it "
              "had nothing holding it open.")
    else:
        bound = ", ".join("%s %s%s (%d atoms)"
                          % (l["residue"], l["chain"], l["resseq"], l["heavy_atoms"])
                          for l in result["bound_ligands"][:6])
        print("bound in this structure:", bound)
    if a.out is not None:
        os.makedirs(os.path.dirname(a.out) or ".", exist_ok=True)
        with open(a.out, "w") as f:
            json.dump(result, f, indent=2)
        print()
        print("written to", a.out)
    # A blocked site is the one failure worth a non zero exit, so a
    # script can stop before spending a day docking into a wall.
    if result["verdict"] == "backbone":
        sys.exit(2)


def cmd_run(a):
    meta = json.load(open(os.path.join(a.decoys, "decoys.json")))
    decoy_list = [d for d in meta["decoys"] if d["heavy"] is not None]
    natives = []
    for path in a.native:
        pdbqt = path
        if path.endswith(".sdf") == True:
            from .prepare import to_pdbqt
            from rdkit import Chem
            mol = Chem.AddHs(Chem.MolFromMolFile(path, removeHs=False),
                             addCoords=True)
            pdbqt = os.path.splitext(path)[0] + ".pdbqt"
            ok, why = to_pdbqt(mol, pdbqt)
            if ok == False:
                sys.exit(f"could not prepare {path}: {why}")
        natives.append({"id": os.path.basename(path).rsplit(".", 1)[0],
                        "pdbqt": pdbqt,
                        "heavy": dk.heavy_atoms_of(path)})
    site = None
    if a.site_ligand is not None:
        try:
            site = rc.check(a.site_receptor or a.receptor, a.site_ligand,
                            ligand_resname=a.site_resname,
                            ligand_resseq=a.site_resseq)
        except ValueError as why:
            sys.exit(str(why))
        print(f"receptor: {site['verdict']}, {site['summary']}")
        if site["verdict"] == "backbone" and a.force == False:
            sys.exit("\nStopping. Docking into this receptor measures the "
                     "protein, not the compounds. Pass --force to run anyway.")
        print()

    print(f"{len(natives)} known binders against {len(decoy_list)} decoys")
    result = run_control(natives, decoy_list, a.receptor, a.center,
                         [a.size] * 3, a.exhaustiveness, a.seed,
                         a.scoring, progress=25)
    if meta.get("salts") is not None:
        result["salts"] = meta["salts"]
    if site is not None:
        result["receptor_check"] = site
    v, md = write_report(result, a.out, a.bar)
    print()
    print(md)


def main(argv=None):
    p = argparse.ArgumentParser(prog="dockcheck")
    sub = p.add_subparsers(dest="cmd", required=True)

    d = sub.add_parser("decoys", help="build a size matched decoy set")
    d.add_argument("--native", required=True, help="native ligand SDF")
    d.add_argument("--out", required=True)
    d.add_argument("--n", type=int, default=200)
    d.add_argument("--tolerance", type=int, default=6)
    d.add_argument("--heavy", type=int, default=None)
    d.add_argument("--mw-min", type=float, default=150.0)
    d.add_argument("--mw-max", type=float, default=600.0)
    d.add_argument("--seed", type=int, default=42)
    d.add_argument("--refetch", action="store_true")
    d.set_defaults(func=cmd_decoys)

    rec = sub.add_parser("receptor",
                         help="check the site is open before docking into it")
    rec.add_argument("--receptor", required=True, help="receptor PDB")
    rec.add_argument("--ligand", default=None,
                     help="PDB holding the bound ligand; defaults to the receptor file")
    rec.add_argument("--resname", default=None,
                     help="residue name of the ligand, e.g. UDP")
    rec.add_argument("--resseq", default=None,
                     help="residue number, when the ligand appears more than once")
    rec.add_argument("--ligand-chain", default=None)
    rec.add_argument("--chain", default=None, help="receptor chain to check")
    rec.add_argument("--out", default=None, help="write the full result as JSON")
    rec.set_defaults(func=cmd_receptor)

    r = sub.add_parser("run", help="run the control and report a verdict")
    r.add_argument("--receptor", required=True)
    r.add_argument("--native", required=True, action="append",
                   help="repeat for each known binder")
    r.add_argument("--decoys", required=True)
    r.add_argument("--center", required=True, type=float, nargs=3)
    r.add_argument("--size", type=float, default=24.0)
    r.add_argument("--out", required=True)
    r.add_argument("--exhaustiveness", type=int, default=8)
    r.add_argument("--seed", type=int, default=42)
    r.add_argument("--scoring", default="vina", choices=["vina", "vinardo"])
    r.add_argument("--bar", type=float, default=PASS_PERCENTILE)
    r.add_argument("--site-ligand", default=None,
                   help="PDB of the ligand that marks the site, to check the "
                        "receptor is open before docking into it")
    r.add_argument("--site-receptor", default=None,
                   help="receptor PDB for that check; defaults to --receptor")
    r.add_argument("--site-resname", default=None)
    r.add_argument("--site-resseq", default=None)
    r.add_argument("--force", action="store_true",
                   help="dock even when the site is blocked by backbone")
    r.set_defaults(func=cmd_run)

    a = p.parse_args(argv)
    a.func(a)


if __name__ == "__main__":
    main()
