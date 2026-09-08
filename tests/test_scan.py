"""Tests for prosy.core.scan: position selection and the three scan types."""

from __future__ import annotations

import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from prosy.core.scan import (  # noqa: E402
    REDUCED_ALPHABETS,
    STANDARD_ALPHABET,
    alanine_scan,
    degenerate_residues,
    degenerate_saturation_positions,
    expand_degenerate,
    parse_positions,
    reduced_alphabet_scan,
    resolve_positions,
    run_scan,
    saturation_scan,
)
from prosy.core.sequence import SequenceError  # noqa: E402

PROTEIN = "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQ"


# --------------------------------------------------------------------------- #
# Position selection                                                           #
# --------------------------------------------------------------------------- #


def test_parse_positions_ranges_and_singletons():
    assert parse_positions("1-3,10,7-8") == [1, 2, 3, 7, 8, 10]
    assert parse_positions("5,5,5") == [5]
    assert parse_positions("all") == []
    assert parse_positions("") == []


def test_parse_positions_rejects_garbage():
    with pytest.raises(SequenceError):
        parse_positions("3-1")
    with pytest.raises(SequenceError):
        parse_positions("x")


def test_resolve_positions_defaults_to_whole_sequence():
    assert resolve_positions(PROTEIN) == list(range(1, len(PROTEIN) + 1))
    assert resolve_positions(PROTEIN, "2-4") == [2, 3, 4]
    assert resolve_positions(PROTEIN, [5, 1, 5]) == [1, 5]
    assert resolve_positions(PROTEIN, "1-5", exclude="2,4") == [1, 3, 5]


def test_resolve_positions_bounds_checked():
    with pytest.raises(SequenceError):
        resolve_positions(PROTEIN, [len(PROTEIN) + 1])
    with pytest.raises(SequenceError):
        resolve_positions(PROTEIN, [0])


# --------------------------------------------------------------------------- #
# 1. Alanine scan                                                              #
# --------------------------------------------------------------------------- #


def test_alanine_scan_whole_sequence_by_default():
    variants = alanine_scan(PROTEIN, parent_name="P")
    assert variants[0].is_parent and variants[0].protein == PROTEIN
    mutants = variants[1:]
    # One per position; native Ala positions become Gly rather than being skipped.
    assert len(mutants) == len(PROTEIN)
    for v in mutants:
        pos = v.position
        expected = "G" if PROTEIN[pos - 1] == "A" else "A"
        assert v.protein[pos - 1] == expected
        assert v.protein[: pos - 1] == PROTEIN[: pos - 1]
        assert v.protein[pos:] == PROTEIN[pos:]
        assert v.scan == "alanine"


def test_alanine_scan_selected_positions_and_skipping_native():
    variants = alanine_scan(PROTEIN, "1-5", include_parent=False, native_substitute=None)
    positions = [v.position for v in variants]
    # Position 4 is Ala in PROTEIN and is skipped when native_substitute is None.
    assert positions == [1, 2, 3, 5]
    assert all(str(v.mutations[0]).endswith("A") for v in variants)


def test_alanine_scan_names_encode_the_mutation():
    v = alanine_scan(PROTEIN, [2], include_parent=False, parent_name="Nb01")[0]
    assert v.name == "Nb01_K2A"
    assert v.mutation_label == "K2A"


# --------------------------------------------------------------------------- #
# 2. Site-saturation                                                           #
# --------------------------------------------------------------------------- #


def test_saturation_scan_gives_19_variants_per_position():
    variants = saturation_scan(PROTEIN, "3-5", include_parent=False)
    assert len(variants) == 3 * 19
    for pos in (3, 4, 5):
        at_pos = [v for v in variants if v.position == pos]
        muts = {v.mutations[0].mut for v in at_pos}
        assert muts == set(STANDARD_ALPHABET) - {PROTEIN[pos - 1]}
        assert all(v.protein[pos - 1] == v.mutations[0].mut for v in at_pos)


def test_saturation_scan_accepts_a_custom_alphabet():
    variants = saturation_scan(PROTEIN, [1], include_parent=False, alphabet="MKW")
    # Position 1 is M, so only K and W remain.
    assert {v.mutations[0].mut for v in variants} == {"K", "W"}


def test_degenerate_codon_expansion():
    assert len(expand_degenerate("NNK")) == 32
    residues, has_stop = degenerate_residues("NNK")
    assert len(residues) == 20 and has_stop is True
    residues, has_stop = degenerate_residues("NDT")
    assert has_stop is False and len(residues) == 12


def test_degenerate_saturation_positions_reports_wild_type():
    rows = degenerate_saturation_positions(PROTEIN, "2-3", codon="NNK")
    assert [(p, wt) for p, wt, _, _ in rows] == [(2, "K"), (3, "T")]


# --------------------------------------------------------------------------- #
# 3. Reduced alphabet                                                          #
# --------------------------------------------------------------------------- #


def test_reduced_alphabet_scan_uses_the_preset():
    variants = reduced_alphabet_scan(PROTEIN, "1-2", include_parent=False)
    alphabet = set(REDUCED_ALPHABETS["adklw"])
    assert {v.mutations[0].mut for v in variants} <= alphabet
    # Position 1 is M (not in the alphabet) -> 5 variants; position 2 is K -> 4.
    assert len([v for v in variants if v.position == 1]) == 5
    assert len([v for v in variants if v.position == 2]) == 4


def test_reduced_alphabet_scan_labels_property_categories():
    variants = reduced_alphabet_scan(PROTEIN, [1], include_parent=False, alphabet="ADKLW")
    categories = {v.mutations[0].mut: v.category for v in variants}
    assert categories == {"A": "small", "D": "negative", "K": "positive",
                          "L": "hydrophobic", "W": "aromatic"}


def test_reduced_alphabet_scan_accepts_explicit_residues():
    variants = reduced_alphabet_scan(PROTEIN, [1], include_parent=False, alphabet="GDKVY")
    assert {v.mutations[0].mut for v in variants} == set("GDKVY")


# --------------------------------------------------------------------------- #
# Dispatch                                                                     #
# --------------------------------------------------------------------------- #


def test_run_scan_dispatches_by_name():
    assert run_scan("alanine", PROTEIN, positions=[1], include_parent=False)[0].scan == "alanine"
    with pytest.raises(SequenceError):
        run_scan("nonsense", PROTEIN)


if __name__ == "__main__":
    raise SystemExit(pytest.main([__file__, "-q"]))
