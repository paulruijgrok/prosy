"""Helpers to assemble the set of forbidden patterns / GC bounds that get
handed to a codon backend. Keeps constraint vocabulary in one place so scripts
read declaratively (``avoid_enzymes=["BsaI"]``) instead of listing raw motifs.
"""

from __future__ import annotations

from dataclasses import dataclass, field

from prosy.core import enzymes


# Low-complexity dinucleotide repeats that synthesis vendors often flag.
DEFAULT_LOW_COMPLEXITY = [
    "ATATAT", "TATATA", "GCGCGC", "CGCGCG",
    "ACACAC", "CACACA", "AGAGAG", "GAGAGA",
]


@dataclass
class ConstraintSet:
    """A declarative bundle of synthesis constraints.

    ``avoid_patterns`` is the flat list a codon backend consumes; build it from
    higher-level fields with :meth:`patterns`.
    """

    avoid_enzymes: list[str] = field(default_factory=list)
    extra_patterns: list[str] = field(default_factory=list)
    max_homopolymer: int | None = None  # e.g. 4 -> forbid 5+ identical bases
    max_homopolymer_by_base: dict[str, int] | None = None  # e.g. {"G": 5, "A": 8}
    forbid_low_complexity: bool = False
    gc_bounds: tuple[float, float] | None = None
    gc_window: int | None = None  # bp; if set, gc_bounds apply per sliding window

    # Additional GC bands, as (window_bp or None, min, max). These stack on top
    # of gc_bounds/gc_window, so a sequence can be held under 72% over 50 bp
    # *and* under 85% over 20 bp - a vendor's local-GC rule is usually a
    # narrower window than the one used to smooth synthesis difficulty.
    gc_windows: list[tuple[int | None, float, float]] = field(default_factory=list)

    # Forbid any repeated k-mer of this size. This is the constraint that keeps
    # a codon-optimized CDS out of a vendor's repeat-complexity penalty; see
    # prosy.core.synthesis. Only the DNAChisel backend can enforce it.
    unique_kmer_size: int | None = None
    unique_kmers_include_rc: bool = True

    # The same idea as an *objective* rather than a hard constraint: push
    # repeats down as far as possible without failing when they cannot be
    # eliminated. Some proteins carry a repeated peptide motif that forces a
    # repeated k-mer whatever the codons, which makes the hard constraint
    # infeasible - the soft form still gets those sequences from ~50% repeat
    # coverage to a few percent.
    soft_unique_kmer_size: int | None = None
    soft_unique_kmer_boost: float = 20.0

    # Minimum relative frequency for any codon used (a "no rare codons" floor).
    min_codon_frequency: float | None = None

    def gc_bands(self) -> list[tuple[int | None, float, float]]:
        """All GC constraints as a flat ``(window, min, max)`` list."""
        bands: list[tuple[int | None, float, float]] = []
        if self.gc_bounds is not None:
            bands.append((self.gc_window, self.gc_bounds[0], self.gc_bounds[1]))
        bands.extend(self.gc_windows)
        return bands

    def patterns(self) -> list[str]:
        out: list[str] = []

        def add(p: str) -> None:
            p = p.upper()
            if p not in out:
                out.append(p)

        for site in enzymes.patterns_for(self.avoid_enzymes):
            add(site)
        if self.max_homopolymer is not None:
            for base in "ACGT":
                add(base * (self.max_homopolymer + 1))
        for base, limit in (self.max_homopolymer_by_base or {}).items():
            add(base.upper() * (limit + 1))
        if self.forbid_low_complexity:
            for p in DEFAULT_LOW_COMPLEXITY:
                add(p)
        for p in self.extra_patterns:
            add(p)
        return out
