"""Unit tests for scripts/nanobodies/rare_codon_scan.py.
Run: python -m pytest tests/ (or python tests/test_rare_codon_scan.py).
"""

import sys
from pathlib import Path

_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(_ROOT))
sys.path.insert(0, str(_ROOT / "scripts" / "nanobodies"))

import rare_codon_scan as rcs  # noqa: E402


def test_codons_split_drops_partial():
    assert rcs.codons_of("ATGAGA") == ["ATG", "AGA"]
    assert rcs.codons_of("ATGAG") == ["ATG"]  # trailing 2 nt dropped


def test_rare_set_membership():
    # classic tRNA-limited rares are flagged; common codons are not
    for c in ("AGA", "AGG", "CGA", "CGG", "ATA", "CTA", "CCC", "GGA"):
        assert c in rcs.RARE
    for c in ("GGC", "GGT", "CTG", "GCA", "ATG"):
        assert c not in rcs.RARE


def test_tandem_run_detected():
    # M | AGA AGG (tandem run, res 2-3) | GGT | GGA (lone rare, res 5) | GAA
    codons = rcs.codons_of("ATG" "AGA" "AGG" "GGT" "GGA" "GAA")
    runs = rcs.rare_runs(codons)
    assert runs == [(2, ["AGA", "AGG"])]          # lone GGA at res 5 excluded


def test_longer_run_and_positions():
    codons = rcs.codons_of("ATG" "CGA" "CGG" "ATA" "GGT")  # run of 3 at res 2-4
    runs = rcs.rare_runs(codons)
    assert runs == [(2, ["CGA", "CGG", "ATA"])]


def test_no_rare_returns_empty():
    codons = rcs.codons_of("ATG" "GGC" "GGT" "CTG" "GCA")
    assert rcs.rare_runs(codons) == []


def test_scan_aggregates():
    rows = [
        {"NbID_parent": "Nb01", "Well": "A1",
         "Coding_DNA": "ATG" "AGA" "AGG" "GGT" "GGA" "GAA"},   # run len2 + 1 lone
        {"NbID_parent": "Nb02", "Well": "B1",
         "Coding_DNA": "ATG" "GGC" "GGT" "CTG"},               # clean
    ]
    s = rcs.scan(rows, "Coding_DNA")
    assert s["n_variants"] == 2
    assert s["total_rare"] == 3                                # AGA, AGG, GGA
    assert s["rare_by_codon"] == {"AGA": 1, "AGG": 1, "GGA": 1}
    assert s["total_runs"] == 1 and s["variants_with_runs"] == 1
    assert s["longest_run"] == 2
    assert s["total_watch"] == 0


def test_scan_counts_ggg_as_watch_not_rare():
    rows = [{"NbID": "x", "Coding_DNA": "ATG" "GGG" "GGG"}]
    s = rcs.scan(rows, "Coding_DNA")
    assert s["total_rare"] == 0        # GGG is 'watch', not counted rare
    assert s["total_watch"] == 2
    assert s["total_runs"] == 0        # watch codons don't form rare runs


def _run_all():
    fns = [v for k, v in sorted(globals().items()) if k.startswith("test_") and callable(v)]
    for fn in fns:
        fn()
        print(f"  ok  {fn.__name__}")
    print(f"\n{len(fns)} tests passed")


if __name__ == "__main__":
    _run_all()
