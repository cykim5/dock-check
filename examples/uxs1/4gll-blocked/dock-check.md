# dock-check: FAIL

**the weakest known binder is at the 13.6th percentile, so a compound that outscores it has shown nothing: 86% of random drugs do too.**

> **Not size matched:** these known binders fall outside the decoy size window, so their percentiles compare molecules of different sizes and should not be read as a result: UDP (25 heavy atoms), UMP (21 heavy atoms).

> Among the size matched binders the worst is at the 30.2th percentile.

> Also: score correlates with heavy atom count at rho=-0.67 across 2738 decoys, so the ranking is substantially reporting molecular size.

## Known binders against size matched decoys

| Ligand | Heavy atoms | Score | Decoys beating it | Percentile |
| --- | --- | --- | --- | --- |
| UDP | 25 | -6.960 | 172 of 199 | 13.6 |
| UMP | 21 | -6.984 | 169 of 199 | 15.1 |
| UDP-xylose | 34 | -7.296 | 139 of 199 | 30.2 |
| UDP-galactose | 36 | -7.866 | 90 of 199 | 54.8 |
| UDP-glucuronate | 37 | -7.938 | 83 of 199 | 58.3 |

## How this was run

- scoring function: `vina`
- exhaustiveness 8, seed 42
- box centre [62.3, 16.4, 34.2], size [28, 28, 28]
- 199 decoys, median score -7.724
- 78.5 minutes

- score against heavy atom count: rho -0.67 (p=1.0e-30, n=2738)

## What this does and does not tell you

Passing means the protocol can find a binder it was not told
about, on this target, in this box. It does not mean any
particular compound in a screen binds, and it does not
transfer to a different box or a different target.

Failing is the stronger result. If a known binder cannot beat
random approved drugs here, then nothing in a ranked screen
on this target has been demonstrated, however good the
scores look.

