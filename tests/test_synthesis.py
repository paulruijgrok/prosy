"""Tests for prosy.core.synthesis: the manufacturability checks a vendor runs.

``VENDOR_REJECTED`` is a **synthetic** 1035 bp fragment that reproduces the
failure mode of a real rejected order: a 333 aa design with de-novo-like
composition, codon-optimized with plain ``use_best_codon`` and no repeat
constraint. It breaches the same three rules the manufacturer quoted, at
almost the same magnitudes (70.2% repeat coverage in both cases).

The numbers asserted below are therefore measurements of *this* sequence, not
quotes from the vendor report. What they pin is the measurement code: if these
drift, the pipeline has stopped measuring what the vendor measures and the gate
is worthless. The implementation was originally validated by reproducing the
vendor's own three figures exactly on the real fragment.
"""

from __future__ import annotations

import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from prosy.core import synthesis  # noqa: E402

# Synthetic stand-in for a rejected fragment. Measured:
#   - repeated 8-mers cover 70.2% of the sequence (limit 40%)
#   - 262 is the start of a 90 base window that is >88% repeat
#   - a 20 base window at position 884 is 100% GC (limit 90%)
VENDOR_REJECTED = (
    "TGTATCGGTCTCGAGGACTGGAAAACGAAAAAGTGGCGGGCCCGGCGGATAGCTGGGGCGAAGTGGGCC"
    "TGCTGAACGTGGTGCCGAGCCGCAACATTGAAGAAGTGCCGGCGGTGTTTCCGAAATTTGAACTGATTA"
    "ACCCGAGCACCAACGTGCTGACCAGCAGCTGCCTGCTGGCGGGCCTGGCGGCGCTGATTTGGAGCAGCC"
    "GCGAACAGGCGGATGGCGAAGCGCTGGTGGCGGAAGTGCTGAAAATTAGCAAACATGGCCGCGCGATTA"
    "ACGCGGGCGGCAACCTGCTGGGCATTGTGGCGGCGGATGCGGATGGCCGCGTGGCGCTGGAAGTGTATC"
    "CGGCGAAAATTCTGGAAGCGCTGAGCGGCGGCGCGCTGAGCAAATGCGCGAGCCAGGAACTGGCGCTGG"
    "CGGGCAAAAAATTTCCGGAACGCAGCCTGGTGGAATGGGCGCTGGAAGCGTTTGGCGCGTTTGTGCGCC"
    "CGGTGAAACGCGTGCTGAAAGATGTGGGCGAAGAACGCTTTCCGAACGTGAGCGTGGAAGCGGCGGAAG"
    "CGGATCCGATTGAACTGGAAGGCGGCGCGGTGGTGGGCAACGCGGGCCTGATTGCGAAAGAAAACCCGC"
    "TGACCCGCGCGTGGCCGCTGAAAAGCGGCCCGGGCCTGCTGCTGCTGGCGGAAGGCGGCGGCGGCGGCA"
    "AAGATGCGGTGCCGAGCGGCGTGATTAAAGCGGAAAACCGCGGCAGCCCGGCGGAATATCTGGCGCCGG"
    "GCGCGGTGCATCCGGGCGTGCTGGGCTTTGTGAACTGGGATATTCCGTGGCGCAAAGTGCGCGGCGGCC"
    "TGTATAACCATCCGGCGCTGATGCGCGTGGCGAACAAATTTAAAGAAGAAGTGATGGCGGGCGCGGGCG"
    "GCGCGGCGATTTATTTTGTGGCGAGCGATCCGGATGTGGTGGAATATCCGGCGACCATTATTAAAGAAG"
    "TGCAGCGCAGCGAACAGCGCGGCCCGCTGACCCTGGGCATGCGCCCGCATGGTTCCGGAGACCTCTAGT"
)


def test_repeat_fraction_of_the_rejected_fragment():
    assert len(VENDOR_REJECTED) == 1035
    frac = synthesis.repeat_fraction(VENDOR_REJECTED, k=8)
    assert round(frac * 100, 1) == 70.2


def test_repeat_fraction_must_count_both_strands():
    # Single-stranded counting scores lower (65.1% here, 59.2% on the real
    # rejected fragment against its 70.2%) and can pass the 40% limit by
    # accident; the vendor counts a k-mer and its reverse complement as one
    # repeat, so this must too.
    one_strand = synthesis.repeat_fraction(
        VENDOR_REJECTED, k=8, include_reverse_complement=False)
    assert round(one_strand * 100, 1) == 65.1


def test_repeat_window_and_gc_hotspot_are_located():
    report = synthesis.check(VENDOR_REJECTED)
    assert report.repeat_window_start == 262
    assert report.repeat_window_fraction > 0.88
    # ...and the GC hotspot, elsewhere in the sequence.
    _, _, hi_pos, hi = report.gc_window_extremes[20]
    assert hi_pos + 1 == 884
    assert round(hi * 100) == 100


def test_the_rejected_sequence_is_reported_as_rejected():
    report = synthesis.check(VENDOR_REJECTED)
    assert not report.ok
    assert any("repeated 8-mers cover" in p for p in report.problems)
    assert any("100% GC" in p for p in report.problems)


def test_a_clean_sequence_passes():
    import random

    rng = random.Random(0)
    seq = "".join(rng.choice("ACGT") for _ in range(1035))
    report = synthesis.check(seq)
    assert report.repeat_fraction < 0.40
    assert report.ok, report.problems


# --------------------------------------------------------------------------- #
# Metric primitives                                                            #
# --------------------------------------------------------------------------- #


def test_repeat_coverage_marks_exactly_the_repeated_kmer():
    # "ACGTACGT" repeated back to back: every base is inside a repeated 8-mer.
    seq = "ACGTACGT" * 4
    assert sum(synthesis.repeat_coverage(seq, 8)) == len(seq)


def test_no_repeat_means_zero_coverage():
    seq = "AAAAAAAACCCCCCCCGGGGGGGGTTTTTTTT"
    mask = synthesis.repeat_coverage(seq, 8, include_reverse_complement=False)
    # Each homopolymer block is its own 8-mer, all distinct.
    assert sum(mask) == 0
    # Counting both strands, A8 pairs with T8 and C8 with G8, so every block
    # becomes a repeat: this is exactly the effect one-strand counting misses.
    assert sum(synthesis.repeat_coverage(seq, 8)) == 32


def test_sequence_shorter_than_k_has_no_repeats():
    assert synthesis.repeat_fraction("ACGT", k=8) == 0.0


def test_densest_repeat_window_finds_the_right_place():
    mask = bytearray([0] * 100 + [1] * 90 + [0] * 100)
    start, frac = synthesis.densest_repeat_window(mask, 90)
    assert (start, frac) == (100, 1.0)


def test_extreme_gc_window_matches_a_naive_scan():
    seq = "ATATATATATGCGCGCGCGCATATATATAT"
    w = 10
    naive = [sum(1 for b in seq[i:i + w] if b in "GC") / w
             for i in range(len(seq) - w + 1)]
    lo_i, lo, hi_i, hi = synthesis.extreme_gc_window(seq, w)
    assert (lo, hi) == pytest.approx((min(naive), max(naive)))
    assert naive[lo_i] == pytest.approx(lo) and naive[hi_i] == pytest.approx(hi)


def test_longest_homopolymers():
    assert synthesis.longest_homopolymers("AAGGGGGTTTC") == {
        "A": 2, "G": 5, "T": 3, "C": 1}


def test_spec_thresholds_are_parameters():
    seq = "ACGTACGT" * 40                      # 100% repeat coverage
    assert not synthesis.check(seq).ok
    permissive = synthesis.SynthesisSpec(
        max_repeat_fraction=1.0, max_repeat_window_fraction=1.0,
        gc_windows=(), gc_bounds=(0.0, 1.0), max_homopolymer={})
    assert synthesis.check(seq, permissive).ok


def test_report_summary_is_readable():
    s = synthesis.check(VENDOR_REJECTED).summary()
    assert "1035 bp" in s and "repeat8 70.2%" in s


if __name__ == "__main__":
    raise SystemExit(pytest.main([__file__, "-q"]))
