"""Framework / CDR annotation for VHH (nanobody) and VH domains.

Dependency-free, motif-anchored annotation. It locates the conserved framework
landmarks that are essentially invariant in VH/VHH domains and reads the CDRs
off as the segments between them:

    FR1 [ C1 ] -- CDR1 -- [ W2 ] FR2 -- CDR2 -- [ FR3 ] FR3 [ C2 ] CDR3 [ WGxG ] FR4

===========  =============================================================
landmark     how it is found
===========  =============================================================
``C1``       last Cys before the FR2 tryptophan (canonical Kabat C22)
``W2``       Trp of the ``W-x-R-Q`` FR2 motif (canonical Kabat W36)
``FR3``      the ``R-[FVLIA]-[TASV]-[ILVM]-[ST]`` motif (canonical R66)
``C2``       last Cys before the FR4 motif (canonical Kabat C92)
``W103``     Trp of the ``W-G-x-G`` FR4 motif (canonical Kabat W103)
===========  =============================================================

Three boundary conventions are available; CDR3 is the same in all of them
(``C2 + 1`` to the residue before the FR4 Trp):

=============  ==========================  ==========================
scheme         CDR1                        CDR2
=============  ==========================  ==========================
``imgt``       ``C1+4`` .. ``W2-3``        ``W2+15`` .. ``FR3-9``
``kabat``      ``C1+9`` .. ``W2-1``        ``W2+14`` .. ``FR3-1``
``extended``   ``C1+4`` .. ``W2-1``        ``W2+14`` .. ``FR3-1``
=============  ==========================  ==========================

``imgt`` is the default because it reproduces the CDR calls already recorded in
the lab's nanobody spreadsheet.

This is a pragmatic annotator for designing mutagenesis libraries, not a
numbering engine. For publication-grade numbering use ANARCI. When annotation
fails, :class:`AnnotationError` is raised so a caller can fall back to explicit
positions rather than silently scanning the wrong residues.
"""

from __future__ import annotations

import re
from dataclasses import dataclass

from prosy.core.sequence import validate_protein

#: Offsets defining each convention:
#: ``(cdr1 start from C1, cdr1 end from W2, cdr2 start from W2, cdr2 end from FR3)``
SCHEME_OFFSETS: dict[str, tuple[int, int, int, int]] = {
    "imgt": (4, -3, 15, -9),
    "kabat": (9, -1, 14, -1),
    "extended": (4, -1, 14, -1),
}
SCHEMES = tuple(SCHEME_OFFSETS)
DEFAULT_SCHEME = "imgt"

_FR2_MOTIF = re.compile(r"W[A-Z]RQ")
_FR2_MOTIF_LOOSE = re.compile(r"W[A-Z]{2}Q")
_FR3_MOTIF = re.compile(r"R[FVLIAY][TASVN][ILVMF][ST]")
_FR4_MOTIF = re.compile(r"WG[A-Z]G")

# Search windows (1-based, inclusive) that keep the motifs from matching inside
# a CDR. Generous enough for the length variation seen across VHH domains.
_FR2_WINDOW = (30, 50)
_FR4_START = 85


class AnnotationError(ValueError):
    """Raised when the conserved framework landmarks cannot be located."""


@dataclass(frozen=True)
class Region:
    """A named, inclusive 1-based residue span."""

    name: str
    start: int
    end: int

    @property
    def positions(self) -> list[int]:
        return list(range(self.start, self.end + 1))

    def __len__(self) -> int:  # noqa: D105
        return self.end - self.start + 1

    def sequence(self, protein: str) -> str:
        return protein[self.start - 1 : self.end]


@dataclass(frozen=True)
class DomainAnnotation:
    """CDR/framework regions of one VH or VHH domain."""

    protein: str
    scheme: str
    regions: tuple[Region, ...]

    def region(self, name: str) -> Region:
        for r in self.regions:
            if r.name.upper() == name.upper():
                return r
        raise KeyError(f"No region named {name!r}; have {[r.name for r in self.regions]}")

    @property
    def cdrs(self) -> tuple[Region, ...]:
        return tuple(r for r in self.regions if r.name.startswith("CDR"))

    def cdr_positions(self, which: str | None = None) -> list[int]:
        """All CDR positions, or those of one CDR (``"CDR3"``)."""
        regions = self.cdrs if which is None else (self.region(which),)
        return sorted(p for r in regions for p in r.positions)

    def framework_positions(self) -> list[int]:
        return sorted(p for r in self.regions
                      if r.name.startswith("FR") for p in r.positions)

    def describe(self) -> str:
        return "  ".join(
            f"{r.name} {r.start}-{r.end} {r.sequence(self.protein)}" for r in self.cdrs
        )


def _search(pattern: re.Pattern, protein: str, window: tuple[int, int]) -> list[int]:
    """1-based start positions of ``pattern`` within an inclusive window."""
    lo, hi = window
    return [m.start() + 1 for m in pattern.finditer(protein) if lo <= m.start() + 1 <= hi]


def annotate(protein: str, *, scheme: str = DEFAULT_SCHEME) -> DomainAnnotation:
    """Annotate the CDRs and frameworks of a VHH/VH domain.

    Raises :class:`AnnotationError` if the conserved landmarks are missing or
    inconsistent (e.g. a truncated domain, or a non-antibody sequence).
    """
    protein = validate_protein(protein)
    try:
        cdr1_from_c1, cdr1_to_w2, cdr2_from_w2, cdr2_to_fr3 = SCHEME_OFFSETS[scheme]
    except KeyError as exc:
        raise AnnotationError(f"Unknown scheme {scheme!r}; expected one of {SCHEMES}.") from exc
    n = len(protein)
    if n < 90:
        raise AnnotationError(f"Sequence is only {n} residues; too short for a VHH domain.")

    # --- FR2 tryptophan (Kabat W36) ---------------------------------------- #
    hits = (_search(_FR2_MOTIF, protein, _FR2_WINDOW)
            or _search(_FR2_MOTIF_LOOSE, protein, _FR2_WINDOW))
    if not hits:
        raise AnnotationError(
            "Could not locate the FR2 'WxRQ' motif between residues "
            f"{_FR2_WINDOW[0]} and {_FR2_WINDOW[1]}.")
    w2 = hits[0]

    # --- first conserved cysteine (Kabat C22) ------------------------------ #
    cys1 = protein.rfind("C", 0, w2 - 5) + 1
    if cys1 == 0:
        raise AnnotationError("Could not locate the FR1 cysteine upstream of the FR2 Trp.")

    # --- FR4 'WGxG' motif (Kabat W103) ------------------------------------- #
    fr4_hits = _search(_FR4_MOTIF, protein, (min(_FR4_START, n), n))
    if not fr4_hits:
        raise AnnotationError("Could not locate the FR4 'WGxG' motif.")
    w103 = fr4_hits[-1]

    # --- second conserved cysteine (Kabat C92) ----------------------------- #
    cys2 = protein.rfind("C", w2, w103 - 1) + 1
    if cys2 == 0:
        raise AnnotationError("Could not locate the FR3 cysteine upstream of 'WGxG'.")

    # --- FR3 anchor (Kabat R66) -------------------------------------------- #
    fr3_hits = _search(_FR3_MOTIF, protein, (w2 + 15, cys2 - 5))
    if not fr3_hits:
        raise AnnotationError(
            "Could not locate the FR3 'R[FVLIA][TASV][ILVM][ST]' motif between "
            f"residues {w2 + 15} and {cys2 - 5}.")
    fr3 = fr3_hits[0]

    cdr1 = Region("CDR1", cys1 + cdr1_from_c1, w2 + cdr1_to_w2)
    cdr2 = Region("CDR2", w2 + cdr2_from_w2, fr3 + cdr2_to_fr3)
    cdr3 = Region("CDR3", cys2 + 1, w103 - 1)
    regions = (
        Region("FR1", 1, cdr1.start - 1), cdr1,
        Region("FR2", cdr1.end + 1, cdr2.start - 1), cdr2,
        Region("FR3", cdr2.end + 1, cys2), cdr3,
        Region("FR4", w103, n),
    )
    for r in regions:
        if r.start > r.end:
            raise AnnotationError(
                f"Inconsistent annotation: {r.name} spans {r.start}-{r.end}. The "
                "sequence may be truncated or not a VH/VHH domain.")
    return DomainAnnotation(protein=protein, scheme=scheme, regions=regions)


def cdr_positions(protein: str, *, scheme: str = DEFAULT_SCHEME,
                  which: str | None = None) -> list[int]:
    """Convenience wrapper: 1-based CDR positions of a VHH/VH domain."""
    return annotate(protein, scheme=scheme).cdr_positions(which)
