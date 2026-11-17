"""Turn the control into a verdict someone can act on.

The bar is set before the run and stated in the output, so the result
cannot be reinterpreted afterwards to suit what was found.
"""
import json

PASS_PERCENTILE = 95.0      # a native ligand should beat 95% of its decoys
WARN_PERCENTILE = 80.0
SIZE_BIAS_CONCERN = -0.4    # Spearman rho below this means size dominates


def _unmatched(result):
    """Known binders whose size falls outside the decoy size range.

    A percentile is only meaningful against decoys of comparable size,
    because the score tracks size. Comparing a 21 atom ligand against
    decoys of 31 to 43 atoms measures size, not binding.
    """
    heavy = [d["heavy"] for d in result.get("decoys", [])
             if d.get("heavy") is not None]
    if len(heavy) == 0:
        return []
    lo, hi = min(heavy), max(heavy)
    return [n for n in result["natives"]
            if n.get("heavy") is not None and (n["heavy"] < lo or n["heavy"] > hi)]


def verdict(result, bar=PASS_PERCENTILE):
    """PASS, WEAK or FAIL, plus the reason, from the worst native."""
    pcts = [n["percentile"] for n in result["natives"]
            if n["percentile"] is not None]
    if len(pcts) == 0:
        return {"verdict": "FAIL", "bar": bar, "worst": None,
                "reason": "no native ligand was scored"}
    worst = min(pcts)
    if worst >= bar:
        v = "PASS"
        reason = (f"every known binder beats at least {bar:.0f}% of "
                  f"size matched decoys")
    elif worst >= WARN_PERCENTILE:
        v = "WEAK"
        reason = (f"the weakest known binder is at the {worst:.1f}th "
                  f"percentile, short of the {bar:.0f}th percentile bar")
    else:
        v = "FAIL"
        reason = (f"the weakest known binder is at the {worst:.1f}th "
                  f"percentile, so a compound that outscores it has shown "
                  f"nothing: {100 - worst:.0f}% of random drugs do too")
    out = {"verdict": v, "bar": bar, "worst_percentile": worst,
           "reason": reason}
    unmatched = _unmatched(result)
    if len(unmatched) > 0:
        out["size_match_warning"] = (
            "these known binders fall outside the decoy size window, so "
            "their percentiles compare molecules of different sizes and "
            "should not be read as a result: "
            + ", ".join(f"{u['id']} ({u['heavy']} heavy atoms)"
                        for u in unmatched))
        matched = [n["percentile"] for n in result["natives"]
                   if n["percentile"] is not None
                   and n["id"] not in {u["id"] for u in unmatched}]
        if len(matched) > 0:
            out["worst_percentile_matched"] = min(matched)
    sb = result.get("size_bias")
    if sb is not None and sb["rho"] <= SIZE_BIAS_CONCERN:
        out["size_bias_warning"] = (
            f"score correlates with heavy atom count at rho={sb['rho']:.2f} "
            f"across {sb['n']} decoys, so the ranking is substantially "
            f"reporting molecular size")
    return out


def to_markdown(result, v):
    L = [f"# dock-check: {v['verdict']}", "",
         f"**{v['reason']}.**", ""]
    if "size_match_warning" in v:
        L += [f"> **Not size matched:** {v['size_match_warning']}.", ""]
        if "worst_percentile_matched" in v:
            L += [f"> Among the size matched binders the worst is at the "
                  f"{v['worst_percentile_matched']:.1f}th percentile.", ""]
    if "size_bias_warning" in v:
        L += [f"> Also: {v['size_bias_warning']}.", ""]
    L += ["## Known binders against size matched decoys", "",
          "| Ligand | Heavy atoms | Score | Decoys beating it | Percentile |",
          "| --- | --- | --- | --- | --- |"]
    for n in sorted(result["natives"], key=lambda r: r["percentile"] or 0):
        L.append(f"| {n['id']} | {n['heavy']} | {n['score']:.3f} | "
                 f"{n['beaten_by']} of {result['n_decoys']} | "
                 f"{n['percentile']:.1f} |")
    L += ["", "## How this was run", "",
          f"- scoring function: `{result['scoring']}`",
          f"- exhaustiveness {result['exhaustiveness']}, seed {result['seed']}",
          f"- box centre {result['center']}, size {result['size']}",
          f"- {result['n_decoys']} decoys, median score "
          f"{result['decoy_median']:.3f}",
          f"- {result['minutes']:.1f} minutes", ""]
    sb = result.get("size_bias")
    if sb is not None:
        L.append(f"- score against heavy atom count: rho {sb['rho']:.2f} "
                 f"(p={sb['p']:.1e}, n={sb['n']})")
    L += ["", "## What this does and does not tell you", "",
          "Passing means the protocol can find a binder it was not told",
          "about, on this target, in this box. It does not mean any",
          "particular compound in a screen binds, and it does not",
          "transfer to a different box or a different target.", "",
          "Failing is the stronger result. If a known binder cannot beat",
          "random approved drugs here, then nothing in a ranked screen",
          "on this target has been demonstrated, however good the",
          "scores look.", ""]
    return "\n".join(L)


def write_report(result, out_dir, bar=PASS_PERCENTILE):
    import os
    os.makedirs(out_dir, exist_ok=True)
    v = verdict(result, bar)
    with open(os.path.join(out_dir, "dock-check.json"), "w") as f:
        json.dump({"verdict": v, "result": result}, f, indent=2)
    md = to_markdown(result, v)
    with open(os.path.join(out_dir, "dock-check.md"), "w") as f:
        f.write(md + "\n")
    return v, md
