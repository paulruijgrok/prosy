"""Codon usage tables and codon-optimization backends.

Two backends share one interface:

* :class:`DnaChiselBackend` (default when ``dnachisel`` is importable) reproduces
  the behaviour of the original ``codon_optimize.py``: reverse-translate, enforce
  translation, avoid forbidden patterns, optimize codon usage for a species.
* :class:`HighestFrequencyBackend` is a dependency-free fallback that picks the
  most frequent synonymous codon per residue and greedily swaps codons to remove
  forbidden patterns. It lets the whole pipeline run (and be tested) without
  ``dnachisel`` installed; it is not a substitute for DNAChisel's global solver.

``get_backend()`` auto-selects DNAChisel if available, else the fallback.
"""

from __future__ import annotations

import random
from typing import Protocol

from prosy.core.sequence import (
    SequenceError,
    translate,
    validate_protein,
)

# --------------------------------------------------------------------------- #
# Codon usage tables                                                          #
# --------------------------------------------------------------------------- #

# Relative synonymous-codon frequencies for E. coli K-12 (fraction within each
# amino acid; values approximate Kazusa/HEG usage). Used to rank codons.
ECOLI_K12: dict[str, dict[str, float]] = {
    "F": {"TTT": 0.58, "TTC": 0.42},
    "L": {"CTG": 0.47, "TTA": 0.14, "TTG": 0.13, "CTT": 0.12, "CTC": 0.10, "CTA": 0.04},
    "I": {"ATT": 0.51, "ATC": 0.39, "ATA": 0.11},
    "M": {"ATG": 1.0},
    "V": {"GTG": 0.35, "GTT": 0.28, "GTC": 0.20, "GTA": 0.17},
    "S": {"AGC": 0.25, "TCT": 0.17, "TCC": 0.15, "AGT": 0.16, "TCA": 0.14, "TCG": 0.14},
    "P": {"CCG": 0.49, "CCA": 0.20, "CCT": 0.18, "CCC": 0.13},
    "T": {"ACC": 0.40, "ACG": 0.25, "ACT": 0.19, "ACA": 0.17},
    "A": {"GCG": 0.34, "GCC": 0.26, "GCA": 0.23, "GCT": 0.18},
    "Y": {"TAT": 0.59, "TAC": 0.41},
    "H": {"CAT": 0.57, "CAC": 0.43},
    "Q": {"CAG": 0.66, "CAA": 0.34},
    "N": {"AAC": 0.51, "AAT": 0.49},
    "K": {"AAA": 0.74, "AAG": 0.26},
    "D": {"GAT": 0.63, "GAC": 0.37},
    "E": {"GAA": 0.68, "GAG": 0.32},
    "C": {"TGC": 0.54, "TGT": 0.46},
    "W": {"TGG": 1.0},
    "R": {"CGC": 0.36, "CGT": 0.36, "CGG": 0.11, "CGA": 0.07, "AGA": 0.07, "AGG": 0.04},
    "G": {"GGC": 0.37, "GGT": 0.35, "GGG": 0.15, "GGA": 0.13},
    "*": {"TAA": 0.61, "TGA": 0.30, "TAG": 0.09},
}

# Registry of named tables for the fallback backend. DNAChisel ships its own
# species tables; these mirror the species names used by the original script.
USAGE_TABLES: dict[str, dict[str, dict[str, float]]] = {
    "e_coli": ECOLI_K12,
}


class CodonUsageTable:
    """Ranks synonymous codons for an amino acid by usage frequency (desc)."""

    def __init__(self, table: dict[str, dict[str, float]]):
        self._weights = table
        self._ranked = {
            aa: sorted(codons, key=lambda c: codons[c], reverse=True)
            for aa, codons in table.items()
        }

    def ranked(self, aa: str) -> list[str]:
        if aa not in self._ranked:
            raise SequenceError(f"No codons for residue {aa!r}")
        return self._ranked[aa]

    def best(self, aa: str) -> str:
        return self.ranked(aa)[0]

    def weight(self, codon: str, aa: str) -> float:
        return self._weights[aa].get(codon, 0.0)


def usage_table(species: str) -> CodonUsageTable:
    try:
        return CodonUsageTable(USAGE_TABLES[species])
    except KeyError as exc:
        raise KeyError(
            f"No built-in fallback codon table for {species!r}. "
            f"Available: {sorted(USAGE_TABLES)}. "
            "Install dnachisel to use its species tables, or register a table."
        ) from exc


# --------------------------------------------------------------------------- #
# Backend interface                                                           #
# --------------------------------------------------------------------------- #


class CodonBackend(Protocol):
    name: str

    def optimize(
        self,
        protein: str,
        *,
        species: str = "e_coli",
        avoid_patterns: list[str] | None = None,
        gc_bounds: tuple[float, float] | None = None,
        gc_window: int | None = None,
        left_context: str = "",
        right_context: str = "",
    ) -> str:
        """Return a coding DNA sequence for ``protein`` meeting the constraints.

        ``gc_bounds`` are applied over the whole coding region; when
        ``gc_window`` is given the GC bounds instead apply to every sliding
        window of that width (bp). ``left_context``/``right_context`` are fixed
        flanking sequences (e.g. cloning adapters) held immutable during
        optimization so that windowed GC spans the flank/coding junctions; only
        the coding region is returned.
        """
        ...


def _aa_index_for_position(pos: int) -> int:
    return pos // 3


class HighestFrequencyBackend:
    """Dependency-free optimizer: best-codon choice + greedy pattern removal.

    Good enough to produce valid, constraint-satisfying sequences for testing
    and for environments without dnachisel. Determinism is controlled by
    ``seed``.
    """

    name = "highest_frequency"

    def __init__(self, seed: int | None = 0):
        self._rng = random.Random(seed)

    def optimize(
        self,
        protein: str,
        *,
        species: str = "e_coli",
        avoid_patterns: list[str] | None = None,
        gc_bounds: tuple[float, float] | None = None,
        gc_window: int | None = None,
        left_context: str = "",
        right_context: str = "",
    ) -> str:
        # This dependency-free fallback enforces GC over the whole coding region
        # only; it does not solve the windowed-GC / fixed-context problem (that
        # needs DNAChisel's global solver). gc_window and the flank contexts are
        # accepted for interface parity but not enforced here.
        protein = validate_protein(protein)
        table = usage_table(species)
        avoid = [p.upper() for p in (avoid_patterns or [])]

        codons = [table.best(aa) for aa in protein]
        codons = self._remove_patterns(codons, protein, table, avoid)
        if gc_bounds is not None:
            codons = self._nudge_gc(codons, protein, table, avoid, gc_bounds)

        dna = "".join(codons)
        # Safety: translation must be preserved.
        if translate(dna) != protein:
            raise SequenceError("Fallback optimizer changed the translation.")
        return dna

    # -- internal helpers ---------------------------------------------------- #

    def _find_pattern(self, dna: str, avoid: list[str]) -> int:
        for pat in avoid:
            idx = dna.find(pat)
            if idx != -1:
                return idx
        return -1

    def _remove_patterns(self, codons, protein, table, avoid, max_iter=10000):
        if not avoid:
            return codons
        for _ in range(max_iter):
            dna = "".join(codons)
            idx = self._find_pattern(dna, avoid)
            if idx == -1:
                return codons
            # Codons spanning [idx, idx+patlen) are candidates to change.
            patlen = max(len(p) for p in avoid)
            first = idx // 3
            last = min(len(codons) - 1, (idx + patlen) // 3)
            changed = False
            order = list(range(first, last + 1))
            self._rng.shuffle(order)
            for ci in order:
                aa = protein[ci]
                alts = [c for c in table.ranked(aa) if c != codons[ci]]
                self._rng.shuffle(alts)
                for alt in alts:
                    trial = codons[:]
                    trial[ci] = alt
                    # Accept a synonymous swap if it clears the occurrence at idx.
                    if alt_removes(trial, avoid, idx):
                        codons = trial
                        changed = True
                        break
                if changed:
                    break
            if not changed:
                raise SequenceError(
                    f"Fallback could not remove forbidden pattern near position {idx}. "
                    "Use the DNAChisel backend for hard cases."
                )
        raise SequenceError("Fallback pattern removal did not converge.")

    def _nudge_gc(self, codons, protein, table, avoid, gc_bounds, max_iter=20000):
        from prosy.core.sequence import gc_content

        lo, hi = gc_bounds
        for _ in range(max_iter):
            dna = "".join(codons)
            gc = gc_content(dna)
            if lo <= gc <= hi:
                return codons
            want_more_gc = gc < lo
            order = list(range(len(codons)))
            self._rng.shuffle(order)
            moved = False
            for ci in order:
                aa = protein[ci]
                cur = codons[ci]
                cur_gc = cur.count("G") + cur.count("C")
                for alt in table.ranked(aa):
                    if alt == cur:
                        continue
                    alt_gc = alt.count("G") + alt.count("C")
                    if (want_more_gc and alt_gc > cur_gc) or (not want_more_gc and alt_gc < cur_gc):
                        trial = codons[:]
                        trial[ci] = alt
                        if self._find_pattern("".join(trial), avoid) == -1:
                            codons = trial
                            moved = True
                            break
                if moved:
                    break
            if not moved:
                return codons  # best effort
        return codons


def alt_removes(codons: list[str], avoid: list[str], idx: int) -> bool:
    """True if no forbidden pattern overlaps original position ``idx`` anymore."""
    dna = "".join(codons)
    for pat in avoid:
        start = 0
        while True:
            j = dna.find(pat, start)
            if j == -1:
                break
            if j <= idx < j + len(pat):
                return False
            start = j + 1
    return True


class DnaChiselBackend:
    """Wraps DNAChisel to reproduce the original codon_optimize.py behaviour."""

    name = "dnachisel"

    def __init__(self):
        import dnachisel  # noqa: F401  (import error surfaces to caller)

    def optimize(
        self,
        protein: str,
        *,
        species: str = "e_coli",
        avoid_patterns: list[str] | None = None,
        gc_bounds: tuple[float, float] | None = None,
        gc_window: int | None = None,
        left_context: str = "",
        right_context: str = "",
    ) -> str:
        from dnachisel import (
            AvoidChanges,
            AvoidPattern,
            CodonOptimize,
            DnaOptimizationProblem,
            EnforceGCContent,
            EnforceTranslation,
            reverse_translate,
        )

        protein = validate_protein(protein)
        coding = reverse_translate(protein)
        left = (left_context or "").upper()
        right = (right_context or "").upper()

        # Optimize the coding region embedded in its fixed flank context, so a
        # windowed GC constraint sees the real neighbourhood at each junction.
        full = left + coding + right
        cstart, cend = len(left), len(left) + len(coding)
        coding_loc = (cstart, cend)

        constraints = [EnforceTranslation(location=coding_loc)]
        if left:
            constraints.append(AvoidChanges(location=(0, cstart)))
        if right:
            constraints.append(AvoidChanges(location=(cend, len(full))))
        # Forbidden patterns apply to the coding region only; the flanks are
        # fixed and may legitimately contain e.g. the Type IIS site themselves.
        for pat in avoid_patterns or []:
            constraints.append(AvoidPattern(pat, location=coding_loc))
        if gc_bounds is not None or gc_window is not None:
            mini, maxi = gc_bounds if gc_bounds is not None else (0.0, 1.0)
            gc_kwargs = {"mini": mini, "maxi": maxi}
            if gc_window is not None:
                gc_kwargs["window"] = gc_window
            constraints.append(EnforceGCContent(**gc_kwargs))

        problem = DnaOptimizationProblem(
            sequence=full,
            constraints=constraints,
            objectives=[CodonOptimize(species=species, location=coding_loc)],
            logger=None,
        )
        problem.resolve_constraints()
        problem.optimize()
        return problem.sequence[cstart:cend]


def dnachisel_available() -> bool:
    try:
        import dnachisel  # noqa: F401

        return True
    except Exception:
        return False


def get_backend(prefer: str | None = None, *, seed: int | None = 0) -> CodonBackend:
    """Return a codon backend.

    ``prefer`` may be ``"dnachisel"``, ``"highest_frequency"`` or ``None``
    (auto: DNAChisel if importable, else the fallback).
    """
    if prefer == "highest_frequency":
        return HighestFrequencyBackend(seed=seed)
    if prefer == "dnachisel":
        return DnaChiselBackend()
    if prefer is None:
        return DnaChiselBackend() if dnachisel_available() else HighestFrequencyBackend(seed=seed)
    raise ValueError(f"Unknown backend {prefer!r}")
