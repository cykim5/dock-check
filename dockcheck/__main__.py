"""Command line entry point.

  python -m dockcheck decoys --native ligand.sdf --out decoys/ --n 200
  python -m dockcheck run --receptor r.pdbqt --native ligand.sdf \
      --decoys decoys/ --center 62.3 16.4 34.2 --size 28 --out report/
"""
import argparse
import csv
import json
import os
import random
import sys

from . import decoys as dk
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
    print(f"{len(natives)} known binders against {len(decoy_list)} decoys")
    result = run_control(natives, decoy_list, a.receptor, a.center,
                         [a.size] * 3, a.exhaustiveness, a.seed,
                         a.scoring, progress=25)
    if meta.get("salts") is not None:
        result["salts"] = meta["salts"]
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
    r.set_defaults(func=cmd_run)

    a = p.parse_args(argv)
    a.func(a)


if __name__ == "__main__":
    main()
