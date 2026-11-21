# dock-check

**Before you trust a virtual screen, check that docking works on your target.**

A screen ranks compounds and hands you a list. Nothing in that list tells you whether the ranking means anything.

There is a cheap way to find out. Dock a ligand your target is already known to bind. Put it up against decoys matched to it in size. See where it lands. If the known binder cannot beat most random drugs then your screen has demonstrated nothing, however good the top scores look.

This runs that check and gives you a verdict.

```
# dock-check: FAIL

the weakest known binder is at the 13.6th percentile, so a compound that
outscores it has shown nothing: 86% of random drugs do too.
```

## Why this is not already routine

Everyone agrees a control is a good idea. Almost no published screen reports one. The reason is that doing it properly takes a day of fiddly work: fetch a library, strip salts, generate 3D conformers, match decoys on size, and hold every docking parameter identical between the two groups.

Each step has a way to go quietly wrong. Two matter most.

**Size.** Empirical scoring functions track heavy atom count. Pick decoys without matching on size and a large native ligand sails past them. The control passes and tells you nothing.

**Counterions.** A drug library is a list of formulated products, not molecules. ChEMBL stores paroxetine hydrochloride and paroxetine mesylate under different IDs. They are one molecule once the salt comes off. Take the entries at face value and you draw both, so a 200 decoy set holds fewer than 200 structures and a repeated parent counts twice wherever it scores.

The counterion is not an inert label while it is still attached. It carries charge and polar groups that a scoring function reads as the drug's own chemistry. The mesylate makes paroxetine look like it owns a sulfonate. The mirror mistake throws out fluorescein for carrying a sodium. The stored molecular weight belongs to the salt too, so a weight filter selects on a molecule that never reaches the docking.

dock-check reduces every entry to its parent before anything is docked. It keeps one decoy per distinct parent and applies your weight window to the parent. The report says how many salt forms it collapsed, which tells you how much smaller your library really is.

## Install

Docking needs AutoDock Vina. Conda is the easiest route:

```bash
conda create -n dockcheck -c conda-forge python=3.12 vina=1.2.7 numpy scipy
conda activate dockcheck
pip install rdkit meeko pytest
```

## Use

Three steps. Check the site is open. Build a size matched decoy set. Run the control.

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

Step 0 takes a second. It is also the step people skip, so it has its own section below. Fold it into step 2 with `--site-ligand` and the run will refuse a blocked site unless you pass `--force`.

Pass `--native` more than once if your target has several known binders. The worst of them sets the verdict. A protocol that finds one binder and misses another is not one you can point at an unknown compound.

## Reading the verdict

| Verdict | Meaning |
|---|---|
| **PASS** | every known binder beats 95% of size matched decoys |
| **WEAK** | the worst known binder lands between the 80th and 95th percentile |
| **FAIL** | below the 80th percentile. The ranking is not usable |

The bar defaults to 95% and moves with `--bar`. Set it before you run, not after you see the number.

Two warnings appear when they apply.

- **Not size matched.** A known binder whose heavy atom count falls outside the decoy range. Its percentile compares molecules of different sizes. Do not read it as a result.
- **Size bias.** The correlation between score and heavy atom count across the decoys. Strongly negative means your ranking is substantially reporting molecular size.

## Check the receptor before you blame the compounds

A failing control tells you the screen is untrustworthy. It does not tell you which part was wrong.

People work the possibilities in a predictable order: the box, the scoring function, the protonation states, the conformers. The receptor comes last because a structure feels like a fact rather than a choice.

It is a choice. A structure solved **without** the ligand bound had nothing holding its site open, so side chains relaxed into the empty space. Dock into that and every compound competes with the protein for the same volume. What comes out is mostly a measure of which molecules can pile up bulk nearby. That is a measure of size.

So put the ligand where the crystal says it goes. Count how much of the receptor is already standing there.

```
# dock-check receptor: BACKBONE

the site is occluded: 15 of 25 ligand atoms clash with the receptor, worst by
2.23 A, and 4 of those clashes are with backbone (worst 1.48 A), which no
flexible side chain setting can move. Dock against a structure solved with
something in this site. Clashing residues: ARG A272, ASN A260, ASN A276, GLY A271

nothing is bound in this structure apart from solvent, which is the condition
this check exists for: a site with nothing in it had nothing holding it open.
```

Three verdicts, because there are three different things to do.

| Verdict | What it means | What to do |
| --- | --- | --- |
| `clear` | nothing is in the way | dock |
| `side-chain` | something is in the way and it can move | try flexible side chain docking before giving up on the structure |
| `backbone` | something is in the way and it cannot move | get a structure solved with something in this site |

`backbone` exits with status 2 so a pipeline can stop before spending a day on it.

Two details decide whether the check is any use.

**Hydrogen bonds are not clashes.** Nitrogen and oxygen may sit 0.6 A closer than their radii suggest. A hydrogen bond pulls two polar atoms to about 2.7 A, and a check that flagged every one of those would condemn every site worth pointing it at.

**Modified residues are protein.** Selenomethionine and its relatives are written as `HETATM` but they are part of the chain. 2B69 carries 80 heavy atoms of selenomethionine. Read those as cargo rather than chain and the receptor has ten holes in it, exactly where a clash might be found.

**What it cannot do.** The check needs a bound ligand to say where the site is. An apo structure is the case that most needs this warning and the one case it cannot handle, which puts you back on pocket detection. The usual answer is to superpose a holo structure of a close homologue. That is what the example below does.

## A worked example

[`examples/uxs1/`](examples/uxs1/) holds the same control run twice against the same enzyme, once into each of two structures of it. UXS1 has a solved structure and no known inhibitor. The project around it screened every approved drug.

The screen produced a confident looking list with a plausible chemical series near the top. The control said the list was worthless. The enzyme's own substrate sat at the 58th percentile, beaten by 83 of 199 random approved drugs.

Then came four days of fixing the wrong thing. Each repair was reasonable. Each one helped. None of them worked.

| What was changed | Substrate percentile |
| --- | --- |
| nothing, as first run | 31.9 |
| different scoring function | 44.2 |
| search box, which sat 8.6 A off the crystallographic ligand | 58.3 |
| ligands ionized, since the substrate carries about three negative charges | 67.3 |
| scoring function and ionization together | 78.9 |

Still short of the bar, fixed at the 95th percentile before the screen ran.

The cause was none of those. The structure was **4GLL, solved with the cofactor bound and the substrate site empty**. ARG272 had relaxed into that empty site. Re-running against 2B69, the same enzyme at 1.21 A with the substrate in place, changed everything the four repairs could not.

| Ligand | Heavy atoms | 4GLL | 2B69 |
| --- | --- | --- | --- |
| UDP-xylose | 34 | 30.2 | **99.5** |
| UDP-galactose | 36 | 54.8 | **99.0** |
| UDP-glucuronate, the substrate | 37 | 58.3 | **88.4** |

Here is the tell. Correlation between score and molecular size over the same decoys:

| Receptor | Spearman rho |
| --- | --- |
| 4GLL, blocked site | **-0.67** |
| 2B69, open site | **+0.11** |

A blocked pocket cannot form real contacts. The only thing left to separate one compound from another is bulk piled up nearby, so the ranking was substantially a ranking by size. Open the pocket and that bias does not shrink. It disappears.

The receptor check takes about a second. It reports 15 of the substrate's 25 atoms clashing and 4 of them against backbone. Running it first would have replaced all four repairs with one instruction: use a different structure.

That is why this exists.

## What a PASS does and does not mean

**It does mean** the protocol found a binder it was not told about, on this target, in this box, against decoys it cannot beat on size alone. Your ranking is worth looking at.

**It does not mean** any particular compound in your screen binds. It does not transfer to a different box, a different protonation state or a different target. It says nothing about whether your top hits are real or the usual false positives, promiscuous binders and aggregators.

A FAIL is the more useful outcome because it is close to conclusive. A known binder that cannot beat random drugs means the scoring function is outside its competence on this target. More screening will not fix that.

## When docking is expected to struggle

Worth knowing before you spend compute. Empirical scoring functions were fitted on neutral drug-like molecules in enclosed greasy pockets, because that is what the structural data is made of. They do poorly on:

- charged ligands, especially phosphates and carboxylates, where the cost of stripping water off the charge is large and not modelled
- polar, water-filled pockets
- large floppy molecules
- ranking chemically unrelated molecules against each other, where the errors do not cancel

If your target hits several of those, run the control first. That is the case this tool was built from.

## Tests

```bash
python -m pytest tests/ -q
```

The suite needs neither Vina nor a network. The last test replays the real UXS1 numbers and asserts the tool calls that screen a failure. That is the regression test that matters.

## Limits

- One scoring function at a time. Vina and Vinardo are supported. smina and GNINA would be the useful next backends.
- The decoy pool is approved drugs from ChEMBL. That is a deliberately hard and realistic set. It is not the right pool for every target.
- One conformer per decoy. Fine for a decoy. Not enough for a compound you care about.
- The control is about the protocol, not about any compound in it.
