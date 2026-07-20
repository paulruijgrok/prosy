"""Unit tests for prosy.core. Run: python -m pytest tests/ (or python tests/test_core.py)."""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from prosy.core import enzymes
from prosy.core.cloning import Flanks, pad_to_length
from prosy.core.codon import HighestFrequencyBackend, get_backend
from prosy.core.constraints import ConstraintSet
from prosy.core.optimize import optimize_cds
from prosy.core.plate import PLATE_96
from prosy.core.sequence import (
    PointMutation, apply_mutation, reverse_complement, translate,
)

PROTEIN = "QVQLVESGGGLVQAGGSLRLSCKTSGFTFEDYVIGWFRQAPGKQREWVATFTSSGDANYADSVKGRFTISRDNAKSTVYLQMNSLKPEDTAVYYCNADVYGWGYTSYSDYWGQGTQVTVSS"


def test_genetic_code_translation():
    assert translate("ATGGCC") == "MA"
    assert reverse_complement("GGTCTC") == "GAGACC"


def test_point_mutation_apply():
    m = PointMutation.parse("Q1E")
    assert str(m) == "Q1E"
    assert apply_mutation("QVQ", m) == "EVQ"


def test_mutation_wt_check():
    try:
        apply_mutation("QVQ", PointMutation.parse("A1E"))
    except Exception as e:
        assert "expected A" in str(e)
    else:
        raise AssertionError("wt mismatch not caught")


def test_enzyme_registry():
    bsai = enzymes.get("BsaI")
    assert bsai.site == "GGTCTC"
    assert set(enzymes.patterns_for("BsaI")) == {"GGTCTC", "GAGACC"}


def test_fallback_optimize_bsai_clean():
    cons = ConstraintSet(avoid_enzymes=["BsaI"])
    res = optimize_cds(PROTEIN, constraints=cons, backend=HighestFrequencyBackend(seed=3))
    assert res.verify()
    assert "GGTCTC" not in res.dna and "GAGACC" not in res.dna


def test_optimize_translation_roundtrip():
    res = optimize_cds("MAAKLEER", backend="highest_frequency")
    assert translate(res.dna) == "MAAKLEER"


def test_windowed_gc_cap_with_flank_context():
    # Windowed GC is enforced by DNAChisel's global solver; skip if unavailable.
    try:
        import dnachisel  # noqa: F401
    except Exception:
        return
    from prosy.core.codon import DnaChiselBackend

    left, right = "TGTATCGGTCTCGAGGA", "GGTTCCGGAGACCTCTAGT"
    cap, window = 0.75, 50
    cons = ConstraintSet(avoid_enzymes=["BsaI"], gc_bounds=(0.0, cap), gc_window=window)
    res = optimize_cds(PROTEIN, constraints=cons, backend=DnaChiselBackend(),
                       left_context=left, right_context=right)
    assert res.verify()
    assert "GGTCTC" not in res.dna and "GAGACC" not in res.dna  # coding stays clean
    final = (left + res.dna + right).upper()
    worst = max(
        (final[i:i + window].count("G") + final[i:i + window].count("C")) / window
        for i in range(len(final) - window + 1)
    )
    assert worst <= cap + 1e-9  # every window, junction included, respects the cap


def test_padding_reaches_min_length_and_preserves_core():
    cons = ConstraintSet(avoid_enzymes=["BsaI"], max_homopolymer=4,
                         forbid_low_complexity=True, gc_bounds=(0.30, 0.70))
    core = Flanks("TGTATCGGTCTCAGGA", "TATTCCGAGACCTCTAGT").apply("ATG" * 10)
    padded = pad_to_length(core, 300, constraints=cons, seed=5)
    assert len(padded) >= 300
    assert core in padded


def test_plate_layout_column_major():
    wells = PLATE_96.wells(order="column")
    assert wells[0] == "A1" and wells[1] == "B1" and wells[7] == "H1"
    assert wells[8] == "A2"
    assert len(wells) == 96


def test_get_backend_auto():
    b = get_backend(None, seed=0)
    assert b.name in {"dnachisel", "highest_frequency"}


def test_platemap_svg_render():
    from prosy.core import platemap
    cells = {
        "A1": platemap.Cell(top="Nb01", bottom="(parent)", category="ADL1", emphasize=True),
        "B1": platemap.Cell(top="Nb73", bottom="N96A", category="ADL1"),
        "A2": platemap.Cell(top="Nb50", bottom="(parent)", category="KRU1", emphasize=True),
    }
    svg = platemap.render_svg(cells, plate=96, title="test")
    assert svg.startswith("<svg") and svg.rstrip().endswith("</svg>")
    assert "Nb73" in svg and "N96A" in svg
    assert "ADL1" in svg and "KRU1" in svg  # legend categories


def _run_all():
    fns = [v for k, v in sorted(globals().items()) if k.startswith("test_") and callable(v)]
    passed = 0
    for fn in fns:
        fn()
        passed += 1
        print(f"  ok  {fn.__name__}")
    print(f"\n{passed}/{len(fns)} tests passed")


if __name__ == "__main__":
    _run_all()
