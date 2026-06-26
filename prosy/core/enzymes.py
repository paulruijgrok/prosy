"""Restriction-enzyme registry.

The point of this module is to keep enzyme knowledge in one place so callers
say ``avoid("BsaI")`` instead of hard-coding ``GGTCTC``/``GAGACC`` everywhere.
Type IIS enzymes (used for Golden Gate) cut at a defined offset outside their
recognition site, captured here for adapter design.
"""

from __future__ import annotations

from dataclasses import dataclass, field

from prosy.core.sequence import reverse_complement


@dataclass(frozen=True)
class Enzyme:
    name: str
    site: str  # recognition site, 5'->3'
    # Type IIS cut offset relative to the END of the recognition site, on the
    # top strand (None for non-IIS or when not needed). e.g. BsaI = 1.
    cut_offset_top: int | None = None
    overhang_length: int | None = None  # length of the single-stranded overhang
    aliases: tuple[str, ...] = field(default_factory=tuple)

    @property
    def site_rc(self) -> str:
        return reverse_complement(self.site)

    def patterns(self) -> list[str]:
        """Both strands' recognition sequences (deduplicated, upper-case)."""
        fwd = self.site.upper()
        rev = self.site_rc.upper()
        return [fwd] if fwd == rev else [fwd, rev]


# Common enzymes. Extend freely; this is intentionally a small, practical set.
_ENZYMES: dict[str, Enzyme] = {}


def register(enzyme: Enzyme) -> None:
    _ENZYMES[enzyme.name.upper()] = enzyme
    for alias in enzyme.aliases:
        _ENZYMES[alias.upper()] = enzyme


for _enz in [
    # Type IIS (Golden Gate workhorses)
    Enzyme("BsaI", "GGTCTC", cut_offset_top=1, overhang_length=4, aliases=("Eco31I",)),
    Enzyme("BsmBI", "CGTCTC", cut_offset_top=1, overhang_length=4, aliases=("Esp3I",)),
    Enzyme("BbsI", "GAAGAC", cut_offset_top=2, overhang_length=4, aliases=("BpiI",)),
    Enzyme("SapI", "GCTCTTC", cut_offset_top=1, overhang_length=3, aliases=("LguI",)),
    Enzyme("AarI", "CACCTGC", cut_offset_top=4, overhang_length=4),
    # Common Type IIP (palindromic) cloning enzymes
    Enzyme("EcoRI", "GAATTC"),
    Enzyme("BamHI", "GGATCC"),
    Enzyme("XhoI", "CTCGAG"),
    Enzyme("NotI", "GCGGCCGC"),
    Enzyme("NdeI", "CATATG"),
    Enzyme("HindIII", "AAGCTT"),
    Enzyme("NcoI", "CCATGG"),
    Enzyme("XbaI", "TCTAGA"),
]:
    register(_enz)


def get(name: str) -> Enzyme:
    try:
        return _ENZYMES[name.upper()]
    except KeyError as exc:
        raise KeyError(
            f"Unknown enzyme {name!r}. Known: {sorted({e.name for e in _ENZYMES.values()})}"
        ) from exc


def patterns_for(names: str | Enzyme | list) -> list[str]:
    """Flatten one or more enzyme names/objects into a list of patterns to avoid."""
    if isinstance(names, (str, Enzyme)):
        names = [names]
    out: list[str] = []
    for n in names:
        enz = n if isinstance(n, Enzyme) else get(n)
        for p in enz.patterns():
            if p not in out:
                out.append(p)
    return out
