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
        original script's default).
    backend:
        A backend instance, a name (``"dnachisel"``/``"highest_frequency"``),
        or None to auto-select.
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
    )
    result = OptimizationResult(protein=protein, dna=dna, backend=backend.name)
    if not result.verify():
        raise RuntimeError("Optimized DNA does not translate back to the input protein.")
    return result
