"""Tests for synthesis profiles, the constraint ladder, and reproducibility.

The reproducibility tests matter more than they look: DNAChisel draws on the
*global* random state, so before the backend seeded it, two identical runs
produced entirely different sequences and an order could not be regenerated.
"""

from __future__ import annotations

import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from prosy.core import codon as codon_mod  # noqa: E402
from prosy.core import synthesis  # noqa: E402
from prosy.core.constraints import ConstraintSet  # noqa: E402
from prosy.core.optimize import build_fragment_with_ladder  # noqa: E402
from prosy.core.sequence import translate  # noqa: E402

PROTEIN = (
    "REDACTED_SEQUENCE"
)
needs_dnachisel = pytest.mark.skipif(
    not codon_mod.dnachisel_available(), reason="dnachisel not installed")


# --------------------------------------------------------------------------- #
# Profiles                                                                     #
# --------------------------------------------------------------------------- #


def test_vendor_standard_constrains_what_the_vendor_measures():
    p = synthesis.get_profile("vendor-standard")
    assert p.unique_kmer_size == 8            # the repeat fix
    assert p.min_codon_frequency == 0.10      # ...without reintroducing rare codons
    assert p.spec.max_repeat_fraction == 0.40
    # The design cap sits below the acceptance limit on purpose: a solver
    # satisfies a cap by sitting against it, so targeting the limit leaves no margin.
    assert p.gc_global[1] < p.spec.gc_bounds[1]


def test_unconstrained_profile_adds_no_manufacturability_constraints():
    p = synthesis.get_profile("none")
    assert p.unique_kmer_size is None and p.min_codon_frequency is None
    (constraints, label), = p.constraint_ladder(["BsaI"])
    assert label == "none"
    assert constraints.unique_kmer_size is None
    assert constraints.soft_unique_kmer_size is None


def test_unknown_profile_raises():
    with pytest.raises(KeyError, match="Unknown synthesis profile"):
        synthesis.get_profile("nonesuch")


# --------------------------------------------------------------------------- #
# The relaxation ladder                                                        #
# --------------------------------------------------------------------------- #


def test_ladder_degrades_hard_constraint_to_soft_objective():
    steps = synthesis.get_profile("vendor-standard").constraint_ladder(["BsaI"])
    labels = [label for _, label in steps]
    assert labels == ["k=8", "k=9,soft_k=8", "k=10,soft_k=8", "soft_k=8"]

    first = steps[0][0]
    assert first.unique_kmer_size == 8 and first.soft_unique_kmer_size is None
    middle = steps[1][0]
    assert middle.unique_kmer_size == 9 and middle.soft_unique_kmer_size == 8
    # The last step must never be infeasible: no hard uniqueness constraint.
    last = steps[-1][0]
    assert last.unique_kmer_size is None and last.soft_unique_kmer_size == 8


def test_every_ladder_step_keeps_the_enzyme_and_gc_constraints():
    steps = synthesis.get_profile("vendor-standard").constraint_ladder(["BsaI", "BsmBI"])
    for constraints, _ in steps:
        assert constraints.avoid_enzymes == ["BsaI", "BsmBI"]
        assert "GGTCTC" in constraints.patterns()
        windows = {w for w, _, _ in constraints.gc_bands()}
        assert {20, 50, None} <= windows
        assert constraints.min_codon_frequency == 0.10


def test_build_fragment_with_ladder_takes_the_first_workable_step():
    impossible = ConstraintSet(extra_patterns=["A", "C", "G", "T"])  # unsatisfiable
    workable = ConstraintSet(avoid_enzymes=["BsaI"])
    fragment, label, relaxed = build_fragment_with_ladder(
        PROTEIN, [(impossible, "strict"), (workable, "loose")],
        backend="highest_frequency")
    assert (label, relaxed) == ("loose", True)
    assert translate(fragment.coding) == PROTEIN


def test_build_fragment_with_ladder_reports_no_step_working():
    impossible = ConstraintSet(extra_patterns=["A", "C", "G", "T"])
    with pytest.raises(RuntimeError, match="No ladder step"):
        build_fragment_with_ladder(PROTEIN, [(impossible, "strict")],
                                   backend="highest_frequency")


def test_build_fragment_with_ladder_rejects_an_empty_ladder():
    with pytest.raises(ValueError):
        build_fragment_with_ladder(PROTEIN, [], backend="highest_frequency")


# --------------------------------------------------------------------------- #
# Reproducibility                                                              #
# --------------------------------------------------------------------------- #


@needs_dnachisel
def test_same_inputs_give_the_same_sequence():
    constraints = synthesis.get_profile("vendor-standard").constraint_ladder(["BsaI"])[0][0]
    a = codon_mod.DnaChiselBackend(seed=0)
    b = codon_mod.DnaChiselBackend(seed=0)
    kwargs = dict(species="e_coli", avoid_patterns=constraints.patterns(),
                  gc_bands=constraints.gc_bands(),
                  unique_kmer_size=constraints.unique_kmer_size,
                  min_codon_frequency=constraints.min_codon_frequency)
    assert a.optimize(PROTEIN, **kwargs) == b.optimize(PROTEIN, **kwargs)
    # ...and repeated calls on the same instance are stable too.
    assert a.optimize(PROTEIN, **kwargs) == a.optimize(PROTEIN, **kwargs)


@needs_dnachisel
def test_any_seed_yields_a_valid_sequence():
    # A different seed may or may not change the result - for a well-determined
    # problem the solver converges to the same answer either way. What must
    # hold for every seed is that the output is correct.
    kwargs = dict(species="e_coli", avoid_patterns=["GGTCTC", "GAGACC"],
                  unique_kmer_size=8)
    for seed in (0, 99, None):
        dna = codon_mod.DnaChiselBackend(seed=seed).optimize(PROTEIN, **kwargs)
        assert translate(dna) == PROTEIN
        assert "GGTCTC" not in dna and "GAGACC" not in dna


@needs_dnachisel
def test_optimization_does_not_disturb_the_callers_random_state():
    import random

    random.seed(1234)
    expected = [random.random() for _ in range(3)]
    random.seed(1234)
    codon_mod.DnaChiselBackend(seed=0).optimize(
        PROTEIN, avoid_patterns=["GGTCTC"], unique_kmer_size=8)
    assert [random.random() for _ in range(3)] == expected


@needs_dnachisel
def test_a_sequence_does_not_depend_on_what_else_is_in_the_batch():
    # Seeding per-input rather than per-call means adding a design to a set
    # leaves every other design's sequence untouched.
    backend = codon_mod.DnaChiselBackend(seed=0)
    kwargs = dict(avoid_patterns=["GGTCTC", "GAGACC"], unique_kmer_size=8)
    alone = backend.optimize(PROTEIN, **kwargs)
    backend.optimize(PROTEIN[:90], **kwargs)          # another job in between
    assert backend.optimize(PROTEIN, **kwargs) == alone


@needs_dnachisel
def test_the_profile_actually_produces_an_orderable_sequence():
    profile = synthesis.get_profile("vendor-standard")
    constraints, _ = profile.constraint_ladder(["BsaI"])[0]
    dna = codon_mod.DnaChiselBackend(seed=0).optimize(
        PROTEIN, avoid_patterns=constraints.patterns(),
        gc_bands=constraints.gc_bands(),
        unique_kmer_size=constraints.unique_kmer_size,
        min_codon_frequency=constraints.min_codon_frequency)
    report = synthesis.check(dna, profile.spec)
    assert report.ok, report.problems
    assert report.repeat_fraction < 0.10


if __name__ == "__main__":
    raise SystemExit(pytest.main([__file__, "-q"]))
