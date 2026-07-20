"""High-level coding-sequence optimization: protein -> synthesis-ready CDS.

This is the general entry point scripts call. It mirrors the original
``codon_optimize.py`` semantics (reverse-translate, enforce translation, avoid
restriction sites, codon-optimize for a species) but takes the enzyme list,
species and GC bounds as parameters instead of hard-coding BsaI.
"""

from __future__ import annotations

from dataclasses import dataclass

from prosy.core import codon
from prosy.core.constraints import ConstraintSet
from prosy.core.sequence import translate, validate_protein


@dataclass
class OptimizationResult:
    protein: str
    dna: str
    backend: str

    def verify(self) -> bool:
        return translate(self.dna) == self.protein


def optimize_cds(
    protein: str,
    *,
    species: str = "e_coli",
    constraints: ConstraintSet | None = None,
    backend: codon.CodonBackend | str | None = None,
    seed: int | None = 0,
    left_context: str = "",
    right_context: str = "",
) -> OptimizationResult:
    """Optimize a single protein into a coding DNA sequence.

    Parameters
    ----------
    protein:
        Single-letter amino-acid sequence.
    species:
        Codon-usage species (e.g. ``"e_coli"``, ``"h_sapiens"`` for DNAChisel).
    constraints:
        A :class:`ConstraintSet`; defaults to avoiding BsaI only (matching the
        original script's default). If ``constraints.gc_window`` is set, its
        ``gc_bounds`` are enforced over every sliding window of that width.
    backend:
        A backend instance, a name (``"dnachisel"``/``"highest_frequency"``),
        or None to auto-select.
    left_context, right_context:
        Fixed flanking sequences (e.g. cloning adapters) held immutable during
        optimization so a windowed GC bound spans the flank/coding junctions.
    """
    protein = validate_protein(protein)
    if constraints is None:
        constraints = ConstraintSet(avoid_enzymes=["BsaI"])
    if backend is None or isinstance(backend, str):
        backend = codon.get_backend(backend, seed=seed)

    dna = backend.optimize(
        protein,
        species=species,
        avoid_patterns=constraints.patterns(),
        gc_bounds=constraints.gc_bounds,
        gc_window=constraints.gc_window,
        left_context=left_context,
        right_context=right_context,
    )
    result = OptimizationResult(protein=protein, dna=dna, backend=backend.name)
    if not result.verify():
        raise RuntimeError("Optimized DNA does not translate back to the input protein.")
    return result


@dataclass
class Fragment:
    """A synthesis-ready fragment: optimized coding region plus flanks/padding."""

    protein: str
    coding: str          # codon-optimized CDS (constraint-clean)
    final: str           # coding wrapped in flanks and padded to min length
    backend: str


def build_fragment(
    protein: str,
    *,
    species: str = "e_coli",
    constraints: ConstraintSet | None = None,
    flanks=None,
    min_length: int = 0,
    pad_constraints: ConstraintSet | None = None,
    backend: codon.CodonBackend | str | None = None,
    seed: int = 0,
) -> Fragment:
    """End-to-end single-protein pipeline: optimize CDS -> add flanks -> pad.

    This is the reusable composition of :func:`optimize_cds`,
    :func:`prosy.core.cloning.add_flanks` and
    :func:`prosy.core.cloning.pad_to_length` that task scripts repeat.
    ``flanks`` is an optional :class:`prosy.core.cloning.Flanks`.
    """
    from prosy.core.cloning import Flanks, add_flanks, pad_to_length

    if backend is None or isinstance(backend, str):
        backend = codon.get_backend(backend, seed=seed)

    flanks_obj = None
    if flanks is not None:
        flanks_obj = flanks if isinstance(flanks, Flanks) else Flanks(*flanks)

    # Optimize the coding region inside its flank context so any windowed GC
    # bound is satisfied across the flank/coding junctions of the final construct.
    result = optimize_cds(
        protein, species=species, constraints=constraints, backend=backend,
        left_context=flanks_obj.five_prime if flanks_obj else "",
        right_context=flanks_obj.three_prime if flanks_obj else "",
    )
    seq = result.dna
    if flanks_obj is not None:
        seq = add_flanks(seq, flanks_obj)
    if min_length and len(seq) < min_length:
        seq = pad_to_length(seq, min_length, constraints=pad_constraints,
                            species=species, seed=seed)
    return Fragment(protein=result.protein, coding=result.dna, final=seq, backend=backend.name)
