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

Three steps. Check the site is open, build a decoy set matched to your native ligand, then run the control.

```bash
# 0. is the site even open?
python -m dockcheck receptor \
    --receptor structure.pdb \
    --resname UDP --resseq 900

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

Step 0 takes a second and is the one people skip. It is described below. You can fold it into step 2 with `--site-ligand`, which refuses to dock into a blocked site unless you pass `--force`.

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

## Check the receptor before you blame the compounds

A failing control tells you the screen is untrustworthy. It does not tell you which part was wrong, and people work through the possibilities in a predictable order: the box, the scoring function, the protonation states, the conformer generation. The receptor comes last, because a structure feels like a fact rather than a choice.

It is a choice. A structure solved **without** the ligand bound has had nothing holding its site open, and side chains relax into the empty space. Dock into that and every compound is competing with the protein for the same volume. The ranking that comes out is mostly a measure of which molecules can pile up bulk nearby, which is to say a measure of size.

So put the ligand where the crystal says it goes, and count how much of the receptor is already standing there:

```
# dock-check receptor: BACKBONE

the site is occluded: 15 of 25 ligand atoms clash with the receptor, worst by
2.23 A, and 4 of those clashes are with backbone (worst 1.48 A), which no
flexible side chain setting can move. Dock against a structure solved with
something in this site. Clashing residues: ARG A272, ASN A260, ASN A276, GLY A271

nothing is bound in this structure apart from solvent, which is the condition
this check exists for: a site with nothing in it had nothing holding it open.
```

Three verdicts, because there are three different things to do:

| Verdict | What it means | What to do |
| --- | --- | --- |
| `clear` | nothing is in the way | dock |
| `side-chain` | something is in the way, and it can move | try flexible side chain docking before giving up on the structure |
| `backbone` | something is in the way, and it cannot move | get a structure solved with something in this site |

`backbone` exits with status 2, so a pipeline can stop before spending a day on it.

Two details that decide whether the check is any use. Nitrogen and oxygen are allowed to sit 0.6 A closer than their radii suggest, because a hydrogen bond legitimately pulls two polar atoms to about 2.7 A, and a check that flagged every hydrogen bond would condemn every site it was pointed at. And chemically modified residues such as selenomethionine are read as protein rather than as bound ligands, even though the file writes them as `HETATM`: 2B69 has 80 heavy atoms of selenomethionine, and treating those as cargo rather than chain leaves ten holes in the receptor exactly where a clash might be found.

**What it cannot do.** The check needs a bound ligand somewhere to say where the site is. If all you have is an apo structure, which is the case that most needs this warning, there is nothing to measure against and you are back to pocket detection. Finding a holo structure of a close homologue and superposing is usually the answer, and that is what the example below does.

## A worked example

[`examples/uxs1/`](examples/uxs1/) holds the same control run twice against the same enzyme, once into each of two structures of it. UXS1 has a solved structure and no known inhibitor, and the project around it was a screen of every approved drug.

The screen produced a confident looking ranked list with a plausible chemical series near the top. The control said the list was worthless. The enzyme's own substrate sat at the 58th percentile, so 83 of 199 random approved drugs outscored the molecule this enzyme evolved to bind.

Then came four days of fixing the wrong thing. Each repair was reasonable, each one helped, and none of them worked:

| What was changed | Substrate percentile |
| --- | --- |
| nothing, as first run | 31.9 |
| different scoring function | 44.2 |
| search box, which was 8.6 A off the crystallographic ligand | 58.3 |
| ligands ionized, since the substrate carries about three negative charges | 67.3 |
| scoring function and ionization together | 78.9 |

Still short of the bar, which had been fixed at the 95th percentile before the screen ran.

The cause was none of those. The structure was **4GLL, solved with the cofactor bound and the substrate site empty**, and ARG272 had relaxed into that empty site. Re-running against 2B69, the same enzyme at 1.21 A solved with the substrate actually in place, changed everything that the four repairs could not:

| Ligand | Heavy atoms | 4GLL | 2B69 |
| --- | --- | --- | --- |
| UDP-xylose | 34 | 30.2 | **99.5** |
| UDP-galactose | 36 | 54.8 | **99.0** |
| UDP-glucuronate, the substrate | 37 | 58.3 | **88.4** |

And the tell, which is the part worth remembering. Correlation between score and molecular size, over the same decoys:

| Receptor | Spearman rho |
| --- | --- |
| 4GLL, blocked site | **-0.67** |
| 2B69, open site | **+0.11** |

A blocked pocket cannot form real contacts, so the only thing separating one compound from another is how much bulk it can pile up nearby. The ranking was substantially a ranking by size. Open the pocket and that bias does not shrink, it disappears.

The receptor check takes about a second and reports 15 of the substrate's 25 atoms clashing, 4 of them with backbone. Running it first would have replaced all four repairs with one instruction: use a different structure.

That is why this exists.

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
