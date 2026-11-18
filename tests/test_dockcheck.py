"""Tests. None of these need Vina or a network.

The last one is the important one: it replays a real control from a
real project and asserts the tool calls it the way a human did.
"""
import json
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from dockcheck.control import percentile_of, size_bias
from dockcheck.decoys import size_window
from dockcheck.report import to_markdown, verdict


def result_from(native_scores, decoy_scores, heavy=37):
    decoys = [{"id": f"D{i}", "score": s, "heavy": heavy}
              for i, s in enumerate(decoy_scores)]
    natives = []
    for name, s in native_scores.items():
        natives.append({"id": name, "score": s, "heavy": heavy,
                        "percentile": percentile_of(s, decoy_scores),
                        "beaten_by": sum(1 for d in decoy_scores if d < s)})
    return {"scoring": "vina", "exhaustiveness": 8, "seed": 42,
            "receptor": "r.pdbqt", "center": [0, 0, 0], "size": [24, 24, 24],
            "natives": natives, "decoys": decoys, "n_decoys": len(decoys),
            "decoy_median": sorted(decoy_scores)[len(decoy_scores) // 2],
            "size_bias": None, "minutes": 1.0}


# ----------------------------------------------------------- percentile
def test_percentile_best_score_is_100():
    assert percentile_of(-10.0, [-5.0, -6.0, -7.0]) == 100.0


def test_percentile_worst_score_is_zero():
    assert percentile_of(-1.0, [-5.0, -6.0, -7.0]) == 0.0


def test_percentile_counts_only_decoys_it_beats():
    # more negative is better, so -6 beats -5 only
    assert percentile_of(-6.0, [-5.0, -7.0, -8.0, -4.0]) == 50.0


def test_percentile_of_nothing_is_none():
    assert percentile_of(-6.0, []) is None


# -------------------------------------------------------- size matching
def test_size_window_is_symmetric():
    assert size_window(37, 6) == (31, 43)


def test_size_bias_detects_a_pure_size_effect():
    heavy = list(range(20, 60))
    scores = [-0.2 * h for h in heavy]      # score is exactly size
    sb = size_bias(scores, heavy)
    assert sb["rho"] < -0.99


def test_size_bias_is_none_on_too_few_points():
    assert size_bias([-1.0], [30]) is None


# -------------------------------------------------------------- verdict
def test_pass_when_every_native_beats_the_bar():
    r = result_from({"known": -12.0}, [-5.0] * 99 + [-13.0])
    v = verdict(r)
    assert v["verdict"] == "PASS"


def test_fail_when_the_native_is_mid_pack():
    r = result_from({"known": -6.0}, [-5.0] * 50 + [-7.0] * 50)
    v = verdict(r)
    assert v["verdict"] == "FAIL"
    assert "shown nothing" in v["reason"]


def test_weak_sits_between_pass_and_fail():
    decoys = [-5.0 - 0.01 * i for i in range(100)]   # -5.00 to -5.99
    r = result_from({"known": -5.85}, decoys)
    v = verdict(r)
    assert v["verdict"] == "WEAK"


def test_the_worst_native_decides_the_verdict():
    r = result_from({"good": -12.0, "bad": -6.0}, [-5.0] * 50 + [-7.0] * 50)
    assert verdict(r)["verdict"] == "FAIL"


def test_size_bias_warning_appears_in_the_verdict():
    r = result_from({"known": -12.0}, [-5.0] * 99 + [-13.0])
    r["size_bias"] = {"rho": -0.73, "p": 1e-20, "n": 100}
    v = verdict(r)
    assert "size_bias_warning" in v
    assert "molecular size" in v["size_bias_warning"]


# --------------------------------------------------------------- report
def test_markdown_leads_with_the_verdict():
    r = result_from({"known": -6.0}, [-5.0] * 50 + [-7.0] * 50)
    md = to_markdown(r, verdict(r))
    assert md.startswith("# dock-check: FAIL")
    assert "does not mean any" in md


# --------------------------------------------- the real historical case
UXS1 = Path(__file__).resolve().parents[1] / "examples" / "uxs1" / "scores.json"


@pytest.mark.skipif(UXS1.exists() == False, reason="example data not present")
def test_replays_the_uxs1_screen_and_fails_it():
    """The real numbers from a real screen that produced no candidate.

    UDP-glucuronate is the enzyme's own substrate. If the tool does not
    call this a failure, the tool is broken.
    """
    data = json.loads(UXS1.read_text())
    r = result_from(data["natives"], data["decoys"], heavy=37)
    v = verdict(r)
    assert v["verdict"] == "FAIL"
    worst = min(n["percentile"] for n in r["natives"])
    assert worst < 70.0


# ------------------------------------------------- size matching caveat
def test_flags_a_native_outside_the_decoy_size_range():
    r = result_from({"small": -6.0}, [-5.0] * 50 + [-7.0] * 50, heavy=37)
    r["natives"][0]["heavy"] = 21          # far below the decoys
    v = verdict(r)
    assert "size_match_warning" in v
    assert "21 heavy atoms" in v["size_match_warning"]


def test_no_warning_when_every_native_is_in_range():
    r = result_from({"ok": -6.0}, [-5.0] * 50 + [-7.0] * 50, heavy=37)
    assert "size_match_warning" not in verdict(r)


def test_the_unmatched_binder_is_named_but_not_judged_on():
    r = result_from({"small": -9.0, "big": -6.0},
                    [-5.0] * 50 + [-7.0] * 50, heavy=37)
    r["natives"][0]["heavy"] = 21
    v = verdict(r)
    assert v["not_judged_on"] == ["small"]
    assert v["judged_on"] == ["big"]
    assert v["worst_percentile"] == pytest.approx(50.0)


# ------------------------------------- heavy atoms, whatever the format
def test_heavy_atoms_from_pdbqt_skips_polar_hydrogens(tmp_path):
    from dockcheck.decoys import heavy_atoms_of
    p = tmp_path / "x.pdbqt"
    # exactly the layout Meeko writes
    p.write_text(
        "ATOM      1  O   UNL     1      -0.884  -2.540  -0.698  1.00  0.00    -0.185 OA\n"
        "ATOM      2  C   UNL     1       0.000   0.000   0.000  1.00  0.00     0.100 C \n"
        "ATOM      3  H   UNL     1       2.000   0.000   0.000  1.00  0.00     0.200 HD\n")
    assert heavy_atoms_of(str(p)) == 2


# ------------------------------------------- one bad ligand is survivable
def test_a_broken_decoy_is_skipped_not_fatal(monkeypatch):
    import dockcheck.control as C

    def fake(path, *a, **k):
        if "BORON" in path:
            raise TypeError("PDBQT parsing error: Atom type B is not a "
                            "valid AutoDock type")
        return -7.0 if "NATIVE" in path else -6.0

    monkeypatch.setattr(C, "dock_one", fake)
    out = C.run_control(
        [{"id": "n", "pdbqt": "NATIVE.pdbqt", "heavy": 37}],
        [{"id": "good", "pdbqt": "A.pdbqt", "heavy": 37},
         {"id": "bad", "pdbqt": "BORON.pdbqt", "heavy": 37}],
        "r.pdbqt", [0, 0, 0], [24, 24, 24])
    assert out["n_decoys"] == 1
    assert len(out["skipped"]) == 1
    assert "Atom type B" in out["skipped"][0]["why"]


def test_a_broken_known_binder_is_fatal(monkeypatch):
    import dockcheck.control as C

    def fake(path, *a, **k):
        if "NATIVE" in path:
            raise TypeError("PDBQT parsing error: Atom type B is not valid")
        return -6.0

    monkeypatch.setattr(C, "dock_one", fake)
    with pytest.raises(RuntimeError, match="control cannot be run"):
        C.run_control([{"id": "n", "pdbqt": "NATIVE.pdbqt", "heavy": 37}],
                      [{"id": "g", "pdbqt": "A.pdbqt", "heavy": 37}],
                      "r.pdbqt", [0, 0, 0], [24, 24, 24])


# --------------------- the verdict must not rest on an unreadable number
def test_an_unmatched_native_does_not_decide_the_verdict():
    """A ligand the report calls unreadable cannot set the outcome."""
    decoys = [-9.0 - 0.01 * i for i in range(100)]      # -9.00 to -9.99
    r = result_from({"matched": -9.98, "tiny": -9.40}, decoys, heavy=37)
    r["natives"][1]["heavy"] = 21                       # outside the window
    v = verdict(r)
    assert v["judged_on"] == ["matched"]
    assert v["not_judged_on"] == ["tiny"]
    assert v["verdict"] == "PASS"


def test_with_no_matched_natives_it_falls_back_to_all_of_them():
    r = result_from({"tiny": -6.0}, [-5.0] * 50 + [-7.0] * 50, heavy=37)
    r["natives"][0]["heavy"] = 21
    v = verdict(r)
    assert v["judged_on"] == ["tiny"]
