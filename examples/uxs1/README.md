# UXS1: the same control, two structures

The same enzyme, the same five known ligands, the same decoys from the
same library, the same size window, the same exhaustiveness and the
same seed. The only difference is which crystal structure the ligands
were docked into.

- [`4gll-blocked/`](4gll-blocked/) is PDB 4GLL, solved with the
  cofactor bound and the substrate site empty. ARG272 had relaxed into
  that empty site.
- [`2b69-open/`](2b69-open/) is PDB 2B69, the same enzyme at 1.21 A
  solved with the substrate actually in place.

| Ligand | Heavy atoms | 4GLL | 2B69 |
| --- | --- | --- | --- |
| UDP-xylose | 34 | 30.2 | 99.5 |
| UDP-galactose | 36 | 54.8 | 99.0 |
| UDP-glucuronate, the substrate | 37 | 58.3 | 88.4 |

Four separate repairs were tried on the 4GLL run first: a corrected
search box, a second scoring function, ionized ligands, and the last
two combined. The best of them reached the 78.9th percentile against a
bar of 95. None of them was the problem.

`dock-check receptor` reports the 4GLL site as occluded in about a
second: 15 of the substrate's 25 atoms clash with the receptor, 4 of
them with backbone, which no docking setting can move.
