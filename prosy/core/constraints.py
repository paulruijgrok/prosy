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
    forbid_low_complexity: bool = False
    gc_bounds: tuple[float, float] | None = None
    gc_window: int | None = None  # bp; if set, gc_bounds apply per sliding window

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
        if self.forbid_low_complexity:
            for p in DEFAULT_LOW_COMPLEXITY:
                add(p)
        for p in self.extra_patterns:
            add(p)
        return out
