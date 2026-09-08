"""Tests for prosy.core.antibody CDR annotation.

The expected CDRs are the ones recorded by hand in the lab's nanobody
spreadsheet ("CDR 1/2/3" columns of the Nanobodies sheet), so this doubles as a
regression test against the annotations already in use.
"""

from __future__ import annotations

import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from prosy.core.antibody import (  # noqa: E402
    AnnotationError,
    annotate,
    cdr_positions,
)

# name -> (sequence, IMGT CDR1, CDR2, CDR3) as recorded in the lab spreadsheet.
NANOBODIES = {
    "Nb01": (
        "QVQLVESGGGLVQAGGSLRLSCKTSGFTFEDYVIGWFRQAPGKQREWVATFTSSGDANYADSVKGRFTISRDNAKSTVYLQMNSLKPEDTAVYYCNADVYGWGYTSYSDYWGQGTQVTVSS",
        "GFTFEDYV", "FTSSGDA", "NADVYGWGYTSYSDY",
    ),
    "Nb02": (
        "QVQLVESGGGLVQAGGSLTLSCAASGLTFSSYAMAWYRQAPGKEREFVAAISGRGISTAYADSVKGRFTISRDNAKNTMYLQMNSLNAEDTAVYYCNALEDDGGLLRMHWGQGTQVTVSS",
        "GLTFSSYA", "ISGRGIST", "NALEDDGGLLRMH",
    ),
    "Nb50": (
        "QVQLQESGGGLVQAGGSLRLSCAASGTIFDDRIMGWYRQAPGKERELVAAIDRGGNTNYADSVKGRFTISRDNAKNTVYLQMNSLKPEDTAVYYCAVSWSADVQDPYTHVYWGQGTQVTVSS",
        "GTIFDDRI", "IDRGGNT", "AVSWSADVQDPYTHVY",
    ),
    "Nb58": (
        "QVQLQESGGGLVQAGGSLRLSCAASGNIFSYTYMGWYRQAPGKERELVASIDNGTNTYYADSVKGRFTISRDNAKNTVYLQMNSLKPEDTAVYYCAAETWYAGTEYVYLDYWGQGTQVTVSS",
        "GNIFSYTY", "IDNGTNT", "AAETWYAGTEYVYLDY",
    ),
}


@pytest.mark.parametrize("name", sorted(NANOBODIES))
def test_imgt_cdrs_match_the_lab_annotation(name):
    protein, cdr1, cdr2, cdr3 = NANOBODIES[name]
    ann = annotate(protein, scheme="imgt")
    assert ann.region("CDR1").sequence(protein) == cdr1
    assert ann.region("CDR2").sequence(protein) == cdr2
    assert ann.region("CDR3").sequence(protein) == cdr3


def test_regions_tile_the_sequence_without_gaps_or_overlaps():
    protein = NANOBODIES["Nb01"][0]
    ann = annotate(protein)
    assert ann.regions[0].start == 1
    assert ann.regions[-1].end == len(protein)
    for left, right in zip(ann.regions, ann.regions[1:]):
        assert right.start == left.end + 1
    assert "".join(r.sequence(protein) for r in ann.regions) == protein


def test_kabat_scheme_gives_the_canonical_h1_and_h2():
    protein = NANOBODIES["Nb01"][0]
    ann = annotate(protein, scheme="kabat")
    assert (ann.region("CDR1").start, ann.region("CDR1").end) == (31, 35)
    assert (ann.region("CDR2").start, ann.region("CDR2").end) == (50, 65)
    assert ann.region("CDR1").sequence(protein) == "DYVIG"


def test_extended_scheme_includes_the_h1_stem():
    protein = NANOBODIES["Nb01"][0]
    assert annotate(protein, scheme="extended").region("CDR1").sequence(protein) == "GFTFEDYVIG"


def test_cdr_positions_helper():
    protein = NANOBODIES["Nb01"][0]
    all_cdrs = cdr_positions(protein)
    cdr3 = cdr_positions(protein, which="CDR3")
    assert cdr3 == list(range(96, 111))
    assert set(cdr3) < set(all_cdrs)
    assert len(all_cdrs) == 8 + 7 + 15


def test_framework_and_cdr_positions_partition_the_domain():
    protein = NANOBODIES["Nb50"][0]
    ann = annotate(protein)
    assert set(ann.cdr_positions()) | set(ann.framework_positions()) == set(
        range(1, len(protein) + 1))
    assert not set(ann.cdr_positions()) & set(ann.framework_positions())


def test_non_antibody_sequence_raises():
    with pytest.raises(AnnotationError):
        annotate("MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQ" * 4)


def test_unknown_scheme_raises():
    with pytest.raises(AnnotationError):
        annotate(NANOBODIES["Nb01"][0], scheme="chothia")


if __name__ == "__main__":
    raise SystemExit(pytest.main([__file__, "-q"]))
