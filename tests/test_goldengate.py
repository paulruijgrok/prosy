"""Tests for prosy.core.goldengate: digestion, adapter design and assembly.

Uses the real FP01 destination vector (``data/plasmids/FP01.fa``) so the
overhangs and the resulting fusion protein are checked against the actual
plasmid rather than a toy construct.
"""

from __future__ import annotations

import sys
from pathlib import Path

import pytest

_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(_ROOT))

from prosy.core.goldengate import (  # noqa: E402
    CloningError,
    assemble,
    check_fragment,
    count_sites,
    design_insert_flanks,
    destination_ends,
    digest,
    find_cuts,
)
from prosy.core.library import naive_cds, orf_protein  # noqa: E402
from prosy.core.sequence import reverse_complement  # noqa: E402

NB01 = (
    "QVQLVESGGGLVQAGGSLRLSCKTSGFTFEDYVIGWFRQAPGKQREWVATFTSSGDANYADSVKGRFTIS"
    "RDNAKSTVYLQMNSLKPEDTAVYYCNADVYGWGYTSYSDYWGQGTQVTVSS"
)
# The vector's payload downstream of the insert, as expressed.
LINKER_AND_TAGS_PREFIX = "GSGSGSGSGSAEIGTGFPFDPHYVEVLGERMHY"


@pytest.fixture(scope="module")
def fp01() -> str:
    text = (_ROOT / "data" / "plasmids" / "FP01.fa").read_text(encoding="utf-8")
    return "".join(ln.strip() for ln in text.splitlines()
                   if ln.strip() and not ln.startswith(">")).upper()


# --------------------------------------------------------------------------- #
# Cut geometry                                                                 #
# --------------------------------------------------------------------------- #


def test_top_strand_cut_leaves_the_overhang_on_the_downstream_fragment():
    # BsaI: GGTCTC(N1) then a 4-nt overhang.
    seq = "AAAAGGTCTCTACGTTTTTTTTTT"
    (cut,) = find_cuts(seq, "BsaI")
    assert cut.strand == 1
    assert cut.overhang == "ACGT"
    assert seq[cut.position : cut.position + 4] == "ACGT"


def test_bottom_strand_site_cuts_to_the_left():
    core = "ACGT"
    seq = "TTTTTTTTTT" + core + "T" + reverse_complement("GGTCTC") + "AAAA"
    (cut,) = find_cuts(seq, "BsaI")
    assert cut.strand == -1
    assert cut.overhang == core


def test_linear_digest_fragment_count_and_overhang_ownership():
    seq = "AAAAGGTCTCTACGTCCCCCCCCCC"
    left, right = digest(seq, "BsaI")
    assert left.left_overhang is None and left.right_overhang == "ACGT"
    assert right.left_overhang == "ACGT" and right.right_overhang is None
    assert right.top.startswith("ACGT")     # downstream fragment owns the overhang
    assert not left.top.endswith("ACGT")    # upstream fragment does not
    assert left.top + right.top == seq


def test_enzyme_without_type_iis_geometry_is_rejected():
    with pytest.raises(CloningError):
        find_cuts("GAATTCAAAA", "EcoRI")


# --------------------------------------------------------------------------- #
# The FP01 destination                                                         #
# --------------------------------------------------------------------------- #


def test_fp01_has_exactly_two_outward_facing_bsai_sites(fp01):
    assert count_sites(fp01, "BsaI", circular=True) == 2
    frags = digest(fp01, "BsaI", circular=True)
    assert len(frags) == 2
    assert sum(1 for f in frags if not f.has_site) == 1


def test_fp01_overhangs_place_the_start_codon_and_open_the_linker(fp01):
    ends = destination_ends(fp01, "BsaI")
    assert (ends.five_overhang, ends.three_overhang) == ("TATG", "GGTT")
    # The ATG of the 5' overhang completes the vector's CATATG start context.
    assert ends.upstream.endswith("CA")
    # The 3' junction runs straight into the GS linker.
    assert ends.downstream.startswith("GGTTCCGGTTCTGGT")
    assert len(ends.backbone) + len(ends.dropout) == len(fp01)


def test_designed_adapters_carry_one_site_each_and_the_right_overhangs(fp01):
    ends = destination_ends(fp01, "BsaI")
    flanks = design_insert_flanks(ends, "BsaI")
    assert flanks.five_prime.endswith("TATG")
    assert flanks.three_prime.startswith("GGTT")
    assert flanks.five_prime.count("GGTCTC") == 1
    assert flanks.three_prime.count("GAGACC") == 1
    assert check_fragment(flanks.apply(naive_cds(NB01)), ends) == []


def test_check_fragment_reports_a_wrong_overhang(fp01):
    ends = destination_ends(fp01, "BsaI")
    bad = "TGTATCGGTCTCAAAAA" + naive_cds(NB01) + "GGTTAGAGACCTCTAGT"
    problems = check_fragment(bad, ends)
    assert any("5' overhang" in p for p in problems)


# --------------------------------------------------------------------------- #
# Assembly                                                                     #
# --------------------------------------------------------------------------- #


def test_assembly_reconstitutes_the_expected_fusion_protein(fp01):
    ends = destination_ends(fp01, "BsaI")
    flanks = design_insert_flanks(ends, "BsaI")
    product = assemble(fp01, [flanks.apply(naive_cds(NB01))], "BsaI")

    # The stuffer is gone and the insert is in.
    assert len(product) == len(ends.backbone) + 3 * len(NB01) + 4
    assert count_sites(product, "BsaI", circular=True) == 0

    protein = orf_protein(product)
    assert protein.startswith("M" + NB01)
    assert protein[1 + len(NB01):].startswith(LINKER_AND_TAGS_PREFIX)
    assert "DYKDDDDK" in protein          # FLAG tag downstream of HaloTag


def test_assembly_fails_loudly_on_mismatched_overhangs(fp01):
    ends = destination_ends(fp01, "BsaI")
    flanks = design_insert_flanks(ends, "BsaI")
    wrong = flanks.five_prime.replace("TATG", "TACG") + naive_cds(NB01) + flanks.three_prime
    with pytest.raises(CloningError):
        assemble(fp01, [wrong], "BsaI")


def test_assembly_rejects_an_insert_with_an_internal_site(fp01):
    ends = destination_ends(fp01, "BsaI")
    flanks = design_insert_flanks(ends, "BsaI")
    coding = naive_cds(NB01)
    contaminated = coding[:30] + "GGTCTC" + coding[36:]
    with pytest.raises(CloningError):
        assemble(fp01, [flanks.apply(contaminated)], "BsaI")


if __name__ == "__main__":
    raise SystemExit(pytest.main([__file__, "-q"]))
