# dock-check

**Before you trust a virtual screen, check that docking works on your target.**

A screen ranks compounds and hands you a list. Nothing in that output tells you whether the ranking means anything. The cheapest way to find out is to dock a ligand your target is already known to bind, against decoys matched to it in size, and see where it lands.

If the known binder cannot beat most random approved drugs, then nothing in your screen has been demonstrated, however good the top scores look.

This runs that check and gives you a verdict.

```
# dock-check: FAIL

the weakest known binder is at the 13.6th percentile, so a compound that
outscores it has shown nothing: 86% of random drugs do too.
```

## Why this is not already routine

Most people know in principle that a control is a good idea. Almost no published virtual screen reports one, and the reason is that doing it properly takes a day of fiddly work: fetching a compound library, stripping salts, generating 3D conformers, matching decoys on size rather than on molecular weight, and keeping every docking parameter identical between the two groups.

Each of those has a way to go quietly wrong, and the most common one is size. Empirical scoring functions correlate strongly with how many heavy atoms a molecule has. Pick decoys without matching on size and a large native ligand sails past them, the control passes, and it has told you nothing.

The second one is subtler, and it is the reason a decoy set is usually smaller than it says. A drug library is a list of formulated products, not a list of molecules. ChEMBL holds paroxetine hydrochloride and paroxetine mesylate under different identifiers, and they are the same molecule once the counterion comes off. Take the stored entries at face value and you draw both, so a decoy set of 200 contains fewer than 200 structures and a repeated parent counts twice wherever it happens to score.

The counterion is not an inert label while it is still attached. It carries the charge and the polar groups that a scoring function reads as chemistry, so the mesylate makes paroxetine look like it owns a sulfonate, and the mirror mistake throws out fluorescein for having a sodium. The stored molecular weight belongs to the salt too, which means a weight filter is selecting on a molecule that never reaches the docking.

dock-check reduces every entry to its parent before anything is docked, keeps one decoy per distinct parent, applies your weight window to the parent rather than to the product, and says in the report how many salt forms it collapsed. That last number is the useful one: it tells you how much smaller your library really is.

## Install

The docking itself needs AutoDock Vina, which is easiest through conda:

```bash
conda create -n dockcheck -c conda-forge python=3.12 vina=1.2.7 numpy scipy
conda activate dockcheck
pip install rdkit meeko pytest
```

## Use

Two steps. Build a decoy set matched to your native ligand, then run the control.

```bash
# 1. decoys the same size as your known binder
python -m dockcheck decoys \
    --native ligand.sdf \
    --out decoys/ \
    --n 200

# 2. the control
python -m dockcheck run \
    --receptor receptor.pdbqt \
    --native ligand.sdf \
    --decoys decoys/ \
    --center 62.3 16.4 34.2 \
    --size 28 \
    --out report/
```

You get `report/dock-check.md` and `report/dock-check.json`.

Pass `--native` more than once if your target has several known binders. The verdict is set by the worst of them, because a protocol that finds one binder and misses another is not one you can point at an unknown compound.

## Reading the verdict

| Verdict | Meaning |
|---|---|
| **PASS** | every known binder beats 95% of size matched decoys |
| **WEAK** | the worst known binder is between the 80th and 95th percentile |
| **FAIL** | below the 80th percentile: the ranking is not usable |

The bar is 95% by default and set with `--bar`. Set it before you run, not after you see the number.

Two extra warnings appear when they apply:

- **Not size matched.** A known binder whose heavy atom count falls outside the decoy range. Its percentile compares molecules of different sizes and should not be read as a result.
- **Size bias.** The correlation between score and heavy atom count across the decoys. Strongly negative means the ranking is substantially reporting molecular size.

## A worked example

[`examples/uxs1/`](examples/uxs1/) holds a real control from a real project: a screen of every approved drug against UXS1, an enzyme with a solved structure and no known inhibitor.

The screen produced a confident looking ranked list with a plausible chemical series near the top. The control says the list is worthless:

| Ligand | Heavy atoms | Score | Percentile |
| --- | --- | --- | --- |
| UDP-glucuronate, the enzyme's own substrate | 37 | -7.938 | **58.3** |
| UDP-galactose | 36 | -7.866 | 54.8 |
| UDP-xylose, the product | 34 | -7.296 | 30.2 |

83 of 199 random approved drugs outscore the molecule this enzyme evolved to bind. The chemical series near the top of the screen did not survive a corrected search box either.

That result is why this exists.

## What a PASS does and does not mean

**Does:** the protocol found a binder it was not told about, on this target, in this box, against decoys it cannot beat on size alone. Your ranking is worth looking at.

**Does not:** that any particular compound in your screen binds. That the result transfers to a different box, a different protonation state, or a different target. That the top of your list is real rather than the known false positive classes, promiscuous binders and aggregators.

A FAIL is the more useful outcome, because it is close to conclusive. A known binder that cannot beat random drugs means the scoring function is outside its competence on this target, and no amount of screening will fix that.

## When docking is expected to struggle

Worth knowing before you spend compute. Empirical scoring functions were fitted on neutral, drug-like molecules in enclosed greasy pockets, because that is what the structural data is made of. They do poorly on:

- charged ligands, phosphates and carboxylates in particular, where the cost of stripping water off the charge is large and not modelled
- polar, water-filled pockets
- large, floppy molecules
- ranking chemically unrelated molecules against each other, where the errors do not cancel

If your target hits several of those, run the control first. That is the case this tool was built from.

## Tests

```bash
python -m pytest tests/ -q
```

The suite needs neither Vina nor a network. The last test replays the real UXS1 numbers and asserts the tool calls that screen a failure, which is the regression test that matters.

## Limits

- One scoring function at a time. Vina and Vinardo are supported; smina and GNINA would be the useful next backends.
- The decoy pool is approved drugs from ChEMBL. That is a deliberately hard, realistic set, and it is not the right pool for every target.
- A single conformer per decoy. Fine for a decoy, not enough for a compound you care about.
- The control is about the protocol, not about any compound in it.
