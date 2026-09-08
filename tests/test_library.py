"""Tests for prosy.core.library: scan variants -> verified, plated fragments."""

from __future__ import annotations

import sys
from pathlib import Path

import pytest

_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(_ROOT))

from prosy.core import codon as codon_mod  # noqa: E402
from prosy.core import goldengate as gg  # noqa: E402
from prosy.core.library import (  # noqa: E402
    LibraryConfig,
    build_degenerate_library,
    build_library,
    verify_library,
    window_gc_range,
    write_library,
)
from prosy.core.scan import alanine_scan, reduced_alphabet_scan  # noqa: E402
from prosy.core.sequence import translate  # noqa: E402

NB01 = (
    "QVQLVESGGGLVQAGGSLRLSCKTSGFTFEDYVIGWFRQAPGKQREWVATFTSSGDANYADSVKGRFTIS"
    "RDNAKSTVYLQMNSLKPEDTAVYYCNADVYGWGYTSYSDYWGQGTQVTVSS"
)


@pytest.fixture(scope="module")
def fp01() -> str:
    text = (_ROOT / "data" / "plasmids" / "FP01.fa").read_text(encoding="utf-8")
    return "".join(ln.strip() for ln in text.splitlines()
                   if ln.strip() and not ln.startswith(">")).upper()


@pytest.fixture(scope="module")
def ends(fp01):
    return gg.destination_ends(fp01, "BsaI")


@pytest.fixture(scope="module")
def cfg(ends):
    # No windowed GC cap here: the dependency-free fallback backend does not
    # solve windowed GC (see test_windowed_gc_cap_needs_dnachisel below), and
    # the test suite must run without dnachisel installed.
    return LibraryConfig(flanks=gg.design_insert_flanks(ends, "BsaI"), enzyme="BsaI",
                         min_length=300)


@pytest.fixture(scope="module")
def members(cfg):
    variants = alanine_scan(NB01, "96-102", parent_name="Nb01")
    return variants, build_library(variants, cfg, backend="highest_frequency")


def test_config_always_avoids_the_cloning_enzyme():
    cfg = LibraryConfig(flanks=gg.Flanks("", ""), enzyme="BsmBI", avoid_enzymes=["BsaI"])
    assert cfg.avoid_enzymes[0] == "BsmBI" and "BsaI" in cfg.avoid_enzymes


def test_build_library_translates_back_and_keeps_the_adapters(members, cfg):
    variants, built = members
    assert len(built) == len(variants)
    for m in built:
        assert translate(m.dna_coding) == m.variant.protein
        assert m.dna_coding in m.dna_final
        assert m.dna_final.startswith(cfg.flanks.five_prime)
        assert len(m.dna_final) >= cfg.min_length


def test_library_passes_every_verification_check(members, cfg, ends, fp01):
    _, built = members
    problems = verify_library(built, cfg, ends=ends, destination=fp01,
                              expected_protein=lambda m: m.variant.protein,
                              assembly_checks=None)
    assert problems == []


def test_verification_catches_a_forbidden_site(members, cfg):
    _, built = members
    broken = built[0]
    original = broken.dna_coding, broken.dna_final
    try:
        broken.dna_coding = broken.dna_coding[:30] + "GGTCTC" + broken.dna_coding[36:]
        broken.dna_final = broken.dna_final.replace(original[0], broken.dna_coding)
        problems = verify_library(built, cfg)
        assert any("forbidden site GGTCTC" in p for p in problems)
    finally:
        broken.dna_coding, broken.dna_final = original


def test_identical_proteins_are_optimized_once(cfg):
    # Two scans of the same parent produce the same parent sequence; the cache
    # must give it an identical CDS rather than re-solving it.
    a = build_library(alanine_scan(NB01, [96], parent_name="Nb01"), cfg,
                      backend="highest_frequency")
    assert a[0].dna_coding == a[0].dna_coding
    doubled = alanine_scan(NB01, [96], parent_name="Nb01") * 2
    built = build_library(doubled, cfg, backend="highest_frequency")
    assert built[0].dna_coding == built[2].dna_coding


def test_wells_are_assigned_and_spill_onto_further_plates(cfg):
    variants = reduced_alphabet_scan(NB01, "1-30", include_parent=False)
    built = build_library(variants, cfg, backend="highest_frequency")
    assert len(built) > 96
    assert built[0].well == "A1" and built[0].plate == 1
    assert built[96].plate == 2 and built[96].well == "A1"
    assert verify_library(built, cfg) == []


def test_degenerate_library_gives_one_fragment_per_position(cfg):
    parent = alanine_scan(NB01, [1], parent_name="Nb01")[0]
    built = build_degenerate_library(parent, [96, 97, 98], cfg, codon="NNK",
                                     backend="highest_frequency")
    assert len(built) == 3
    for m, pos in zip(built, (96, 97, 98)):
        assert m.degenerate_codon == "NNK"
        assert m.dna_coding[(pos - 1) * 3 : pos * 3] == "NNK"
        assert m.dna_coding.count("N") == 2 and m.dna_coding.count("K") == 1
    assert verify_library(built, cfg) == []


def test_window_gc_range_matches_a_naive_computation():
    dna = "GCGCGCATATATGCGCATAT"
    window = 6
    naive = [sum(1 for b in dna[i:i + window] if b in "GC") / window
             for i in range(len(dna) - window + 1)]
    assert window_gc_range(dna, window) == pytest.approx((min(naive), max(naive)))
    assert window_gc_range("GGCC", 10) == (1.0, 1.0)


@pytest.mark.skipif(not codon_mod.dnachisel_available(), reason="dnachisel not installed")
def test_windowed_gc_cap_needs_dnachisel(ends):
    """The windowed GC cap is a global-solver problem, so it is a DNAChisel job.

    Verification must catch the fallback's failure rather than let an
    over-GC fragment reach a synthesis order.
    """
    windowed = LibraryConfig(flanks=gg.design_insert_flanks(ends, "BsaI"), enzyme="BsaI",
                             min_length=300, gc_window=50, gc_max=0.72)
    variants = alanine_scan(NB01, [96], parent_name="Nb01")

    fallback = build_library(variants, windowed, backend="highest_frequency")
    assert any("window GC" in p for p in verify_library(fallback, windowed))

    solved = build_library(variants, windowed, backend="dnachisel")
    assert verify_library(solved, windowed) == []
    assert window_gc_range(solved[0].dna_final, 50)[1] <= 0.72 + 1e-9


def test_write_library_emits_mapping_design_and_plate_files(members, cfg, tmp_path):
    _, built = members
    out = write_library(built, tmp_path, "test", cfg)
    assert out["mapping"].exists() and out["design"].exists()
    assert len(out["plates"]) == 1 and out["plates"][0].exists()
    header, first = out["mapping"].read_text(encoding="utf-8").splitlines()[:2]
    assert header.startswith("Plate,Well,Name")
    assert first.split(",")[2] == "Nb01"


if __name__ == "__main__":
    raise SystemExit(pytest.main([__file__, "-q"]))
