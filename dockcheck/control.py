"""Run the control: dock the natives and the decoys, same box, same settings.

Everything that could differ between the two groups is held identical
on purpose. Same receptor, same box, same exhaustiveness, same seed.
The only difference is which molecule is in the box, so the percentile
means what it says.
"""
import time
from statistics import median

from scipy.stats import spearmanr


def dock_one(ligand_pdbqt, receptor_pdbqt, center, size,
             exhaustiveness=8, seed=42, scoring="vina"):
    """Best score for one ligand. Imported lazily so tests need no Vina."""
    from vina import Vina
    v = Vina(sf_name=scoring, verbosity=0, seed=seed)
    v.set_receptor(receptor_pdbqt)
    v.set_ligand_from_file(ligand_pdbqt)
    v.compute_vina_maps(center=list(center), box_size=list(size))
    v.dock(exhaustiveness=exhaustiveness, n_poses=9)
    return float(v.energies(n_poses=1)[0][0])


def percentile_of(score, decoy_scores):
    """What share of decoys this score beats. More negative is better."""
    n = len(decoy_scores)
    if n == 0:
        return None
    beaten = sum(1 for d in decoy_scores if d > score)
    return 100.0 * beaten / n


def size_bias(scores, heavy_counts):
    """Spearman correlation between score and heavy atom count.

    A strongly negative value means the function is largely reporting
    molecular size. Reported because it is the most common reason a
    screen looks like it worked.
    """
    if len(scores) < 5:
        return None
    rho, p = spearmanr(heavy_counts, scores)
    return {"rho": float(rho), "p": float(p), "n": len(scores)}


def run_control(natives, decoys, receptor, center, size,
                exhaustiveness=8, seed=42, scoring="vina", progress=None):
    """natives and decoys are lists of {id, pdbqt, heavy}.

    Returns a dict the report module turns into a verdict.
    """
    started = time.time()
    decoy_scores, decoy_heavy, decoy_rows = [], [], []
    for i, d in enumerate(decoys):
        s = dock_one(d["pdbqt"], receptor, center, size,
                     exhaustiveness, seed, scoring)
        decoy_scores.append(s)
        decoy_heavy.append(d["heavy"])
        decoy_rows.append({"id": d["id"], "score": s, "heavy": d["heavy"]})
        if progress is not None and (i + 1) % progress == 0:
            print(f"  {i + 1}/{len(decoys)} decoys, "
                  f"{(time.time() - started) / 60:.1f} min", flush=True)

    native_rows = []
    for n in natives:
        s = dock_one(n["pdbqt"], receptor, center, size,
                     exhaustiveness, seed, scoring)
        native_rows.append({"id": n["id"], "score": s, "heavy": n["heavy"],
                            "percentile": percentile_of(s, decoy_scores),
                            "beaten_by": sum(1 for d in decoy_scores if d < s)})

    return {
        "scoring": scoring, "exhaustiveness": exhaustiveness, "seed": seed,
        "receptor": receptor, "center": list(center), "size": list(size),
        "natives": native_rows, "decoys": decoy_rows,
        "n_decoys": len(decoy_rows),
        "decoy_median": median(decoy_scores) if decoy_scores else None,
        "size_bias": size_bias(decoy_scores, decoy_heavy),
        "minutes": (time.time() - started) / 60.0,
    }
