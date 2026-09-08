"""Mutational scans: turn one protein into a set of single-mutant proteins.

The module is deliberately generic. Everything that differs between a scan
types is a parameter:

* *which* positions are scanned -> :func:`resolve_positions`
* *what* each position is mutated to -> the ``alphabet`` argument of :func:`scan`

The three named scans in this module are thin wrappers around :func:`scan`:

===========================  ====================================================
:func:`alanine_scan`         one variant per position, wild-type -> Ala
:func:`saturation_scan`      site-saturation: 19 variants per position
:func:`reduced_alphabet_scan` a small property-spanning alphabet per position
===========================  ====================================================

Position selection is explicit by design: pass ``positions`` to scan a subset,
or leave it as ``None`` to scan the whole sequence. Domain-specific defaults
(e.g. "CDRs only" for antibodies) belong in the caller - see
:mod:`prosy.core.antibody` and ``scripts/nanobodies/nanobody_scan.py``.
"""

from __future__ import annotations

import itertools
from collections.abc import Callable, Iterable, Sequence
from dataclasses import dataclass, field

from prosy.core.sequence import (
    AMINO_ACIDS,
    PointMutation,
    SequenceError,
    apply_mutation,
    validate_protein,
)

# The 20 proteinogenic residues in a stable, conventional order.
STANDARD_ALPHABET = "ACDEFGHIKLMNPQRSTVWY"

# --------------------------------------------------------------------------- #
# Reduced alphabets                                                            #
# --------------------------------------------------------------------------- #

# Property class of each residue, used to label reduced-alphabet variants.
RESIDUE_CLASS: dict[str, str] = {
    "A": "small", "G": "small", "S": "polar", "T": "polar", "C": "polar",
    "N": "polar", "Q": "polar", "P": "special",
    "D": "negative", "E": "negative",
    "K": "positive", "R": "positive", "H": "positive",
    "V": "hydrophobic", "L": "hydrophobic", "I": "hydrophobic", "M": "hydrophobic",
    "F": "aromatic", "W": "aromatic", "Y": "aromatic",
}

#: Named 5-residue alphabets that span charge / hydrophobicity / size.
REDUCED_ALPHABETS: dict[str, str] = {
    # small-neutral, negative, positive, aliphatic-hydrophobic, aromatic-bulky
    "adklw": "ADKLW",
    # small, negative, positive (Arg), aromatic, polar-neutral
    "adrfs": "ADRFS",
    # tiny/flexible, negative, positive, beta-branched hydrophobic, aromatic-polar
    "gdkvy": "GDKVY",
}
DEFAULT_REDUCED_ALPHABET = "adklw"


# --------------------------------------------------------------------------- #
# Degenerate codons (for library-scale saturation at the DNA level)            #
# --------------------------------------------------------------------------- #

IUPAC: dict[str, str] = {
    "A": "A", "C": "C", "G": "G", "T": "T",
    "R": "AG", "Y": "CT", "S": "CG", "W": "AT", "K": "GT", "M": "AC",
    "B": "CGT", "D": "AGT", "H": "ACT", "V": "ACG", "N": "ACGT",
}

#: Common degenerate codons for saturation libraries, with what they encode.
DEGENERATE_CODONS = ("NNK", "NNS", "NNN", "NDT", "DBK", "VHG")


def expand_degenerate(codon: str) -> list[str]:
    """Expand an IUPAC-degenerate codon into its concrete DNA codons."""
    codon = codon.strip().upper()
    if len(codon) != 3:
        raise SequenceError(f"Degenerate codon must be 3 characters: {codon!r}")
    try:
        pools = [IUPAC[c] for c in codon]
    except KeyError as exc:
        raise SequenceError(f"Unknown IUPAC symbol in {codon!r}: {exc}") from exc
    return ["".join(t) for t in itertools.product(*pools)]


def degenerate_residues(codon: str) -> tuple[set[str], bool]:
    """Residues encoded by a degenerate codon, and whether it includes a stop."""
    from prosy.core.sequence import CODON_TABLE

    aas = {CODON_TABLE[c] for c in expand_degenerate(codon)}
    return ({a for a in aas if a != "*"}, "*" in aas)


# --------------------------------------------------------------------------- #
# Variant model                                                                #
# --------------------------------------------------------------------------- #


@dataclass(frozen=True)
class ScanVariant:
    """One member of a scan: a protein sequence plus how it was derived."""

    name: str                       # e.g. "Nb01_T50A"
    parent: str                     # parent name, e.g. "Nb01"
    protein: str                    # full mutated protein sequence
    mutations: tuple[PointMutation, ...] = ()
    scan: str = "parent"            # "parent" | "alanine" | "saturation" | "reduced"
    category: str | None = None     # property class of the substituted residue
    meta: dict = field(default_factory=dict)

    @property
    def is_parent(self) -> bool:
        return not self.mutations

    @property
    def mutation_label(self) -> str:
        if self.mutations:
            return "+".join(str(m) for m in self.mutations)
        if self.scan == "degenerate":
            return f"{self.meta['wild_type']}{self.meta['position']}{self.meta['codon']}"
        return "WT"

    @property
    def position(self) -> int | None:
        """Scanned position, for single mutants."""
        return self.mutations[0].position if len(self.mutations) == 1 else None


# --------------------------------------------------------------------------- #
# Position selection                                                           #
# --------------------------------------------------------------------------- #


def parse_positions(spec: str) -> list[int]:
    """Parse ``"1-5,10,20-25"`` into a sorted list of unique 1-based positions.

    ``"all"`` / ``""`` return an empty list, meaning "caller decides"
    (:func:`resolve_positions` then falls back to the whole sequence).
    """
    spec = (spec or "").strip().lower()
    if spec in {"", "all", "*"}:
        return []
    out: set[int] = set()
    for chunk in spec.replace(";", ",").split(","):
        chunk = chunk.strip()
        if not chunk:
            continue
        if "-" in chunk:
            lo_s, _, hi_s = chunk.partition("-")
            try:
                lo, hi = int(lo_s), int(hi_s)
            except ValueError as exc:
                raise SequenceError(f"Bad position range {chunk!r}") from exc
            if lo > hi:
                raise SequenceError(f"Inverted position range {chunk!r}")
            out.update(range(lo, hi + 1))
        else:
            try:
                out.add(int(chunk))
            except ValueError as exc:
                raise SequenceError(f"Bad position {chunk!r}") from exc
    return sorted(out)


def resolve_positions(
    protein: str,
    positions: Iterable[int] | str | None = None,
    *,
    exclude: Iterable[int] | str | None = None,
) -> list[int]:
    """Return the sorted 1-based positions to scan.

    ``positions`` may be an iterable of ints, a spec string (see
    :func:`parse_positions`), or ``None``/empty - in which case **every**
    residue of ``protein`` is scanned. ``exclude`` is subtracted afterwards.
    All positions are bounds-checked against ``protein``.
    """
    protein = validate_protein(protein)
    n = len(protein)

    if isinstance(positions, str):
        positions = parse_positions(positions)
    chosen = sorted(set(positions)) if positions else list(range(1, n + 1))

    if isinstance(exclude, str):
        exclude = parse_positions(exclude)
    if exclude:
        chosen = [p for p in chosen if p not in set(exclude)]

    bad = [p for p in chosen if not 1 <= p <= n]
    if bad:
        raise SequenceError(
            f"Position(s) {bad} out of range for a {n}-residue protein."
        )
    return chosen


# --------------------------------------------------------------------------- #
# The scan engine                                                              #
# --------------------------------------------------------------------------- #

# An alphabet is either a fixed set of residues or a function of (position, wt).
AlphabetSpec = str | Sequence[str] | Callable[[int, str], Sequence[str]]


def _alphabet_at(alphabet: AlphabetSpec, position: int, wt: str) -> list[str]:
    residues = alphabet(position, wt) if callable(alphabet) else alphabet
    out = [r.upper() for r in residues]
    bad = sorted(set(out) - AMINO_ACIDS)
    if bad:
        raise SequenceError(f"Alphabet contains non-standard residues: {bad}")
    return out


def scan(
    protein: str,
    *,
    alphabet: AlphabetSpec,
    positions: Iterable[int] | str | None = None,
    exclude: Iterable[int] | str | None = None,
    parent_name: str = "parent",
    include_parent: bool = True,
    skip_wild_type: bool = True,
    scan_name: str = "scan",
) -> list[ScanVariant]:
    """Generic single-substitution scan.

    Emits one :class:`ScanVariant` per (position, substitution) pair, grouped by
    position in ascending order, optionally preceded by the unmutated parent.

    Parameters
    ----------
    alphabet:
        Residues to substitute in. Either a fixed alphabet (``"ADKLW"``) or a
        callable ``(position, wt) -> residues`` for position-dependent sets.
    positions:
        Positions to scan (ints or a ``"1-5,10"`` spec). ``None`` = whole sequence.
    skip_wild_type:
        Drop substitutions that would re-create the wild-type residue.
    """
    protein = validate_protein(protein)
    targets = resolve_positions(protein, positions, exclude=exclude)

    variants: list[ScanVariant] = []
    if include_parent:
        variants.append(
            ScanVariant(name=parent_name, parent=parent_name, protein=protein,
                        mutations=(), scan="parent")
        )

    for pos in targets:
        wt = protein[pos - 1]
        for mut_aa in _alphabet_at(alphabet, pos, wt):
            if skip_wild_type and mut_aa == wt:
                continue
            mutation = PointMutation(wt=wt, position=pos, mut=mut_aa)
            variants.append(ScanVariant(
                name=f"{parent_name}_{mutation}",
                parent=parent_name,
                protein=apply_mutation(protein, mutation, check_wt=True),
                mutations=(mutation,),
                scan=scan_name,
                category=RESIDUE_CLASS.get(mut_aa),
            ))
    return variants


# --------------------------------------------------------------------------- #
# 1. Alanine scan                                                              #
# --------------------------------------------------------------------------- #


def alanine_scan(
    protein: str,
    positions: Iterable[int] | str | None = None,
    *,
    exclude: Iterable[int] | str | None = None,
    parent_name: str = "parent",
    include_parent: bool = True,
    substitute: str = "A",
    native_substitute: str | None = "G",
) -> list[ScanVariant]:
    """Alanine scan of ``protein`` (whole sequence unless ``positions`` given).

    Positions whose wild-type residue is already the ``substitute`` residue
    cannot be scanned with it; ``native_substitute`` (Gly by convention) is used
    there instead. Set it to ``None`` to skip those positions entirely.
    """
    substitute = substitute.upper()

    def alphabet(position: int, wt: str) -> list[str]:
        if wt != substitute:
            return [substitute]
        return [native_substitute.upper()] if native_substitute else []

    return scan(
        protein, alphabet=alphabet, positions=positions, exclude=exclude,
        parent_name=parent_name, include_parent=include_parent,
        scan_name="alanine",
    )


# --------------------------------------------------------------------------- #
# 2. Site-saturation mutagenesis                                               #
# --------------------------------------------------------------------------- #


def saturation_scan(
    protein: str,
    positions: Iterable[int] | str | None = None,
    *,
    exclude: Iterable[int] | str | None = None,
    parent_name: str = "parent",
    include_parent: bool = True,
    alphabet: str = STANDARD_ALPHABET,
) -> list[ScanVariant]:
    """Site-saturation mutagenesis: every non-wild-type residue at each position.

    With the default alphabet this yields 19 variants per position, each an
    explicitly defined protein sequence (as opposed to a degenerate-codon
    library - see :func:`degenerate_saturation_positions`).
    """
    return scan(
        protein, alphabet=alphabet, positions=positions, exclude=exclude,
        parent_name=parent_name, include_parent=include_parent,
        scan_name="saturation",
    )


def degenerate_saturation_positions(
    protein: str,
    positions: Iterable[int] | str | None = None,
    *,
    exclude: Iterable[int] | str | None = None,
    codon: str = "NNK",
) -> list[tuple[int, str, set[str], bool]]:
    """Positions for a degenerate-codon saturation library.

    Returns ``(position, wild_type, encoded_residues, encodes_stop)`` per
    position. One DNA fragment per position is built by
    :func:`prosy.core.library.build_degenerate_library`; the protein-level
    identity of each library member is only defined after sequencing, which is
    why this returns positions rather than :class:`ScanVariant` objects.
    """
    protein = validate_protein(protein)
    targets = resolve_positions(protein, positions, exclude=exclude)
    residues, has_stop = degenerate_residues(codon)
    return [(p, protein[p - 1], residues, has_stop) for p in targets]


# --------------------------------------------------------------------------- #
# 3. Reduced-alphabet (property) scan                                          #
# --------------------------------------------------------------------------- #


def reduced_alphabet_scan(
    protein: str,
    positions: Iterable[int] | str | None = None,
    *,
    exclude: Iterable[int] | str | None = None,
    parent_name: str = "parent",
    include_parent: bool = True,
    alphabet: str = DEFAULT_REDUCED_ALPHABET,
) -> list[ScanVariant]:
    """Saturation restricted to a small alphabet spanning charge/hydrophobicity.

    ``alphabet`` is either a preset name from :data:`REDUCED_ALPHABETS`
    (default ``"adklw"`` = A/D/K/L/W: small, negative, positive, aliphatic,
    aromatic) or an explicit residue string such as ``"ADKLW"``.

    At most ``len(alphabet)`` variants per position; a position whose wild-type
    residue is in the alphabet yields one fewer.
    """
    residues = REDUCED_ALPHABETS.get(alphabet.lower(), alphabet).upper()
    if len(set(residues)) != len(residues):
        raise SequenceError(f"Reduced alphabet has duplicates: {residues!r}")
    return scan(
        protein, alphabet=residues, positions=positions, exclude=exclude,
        parent_name=parent_name, include_parent=include_parent,
        scan_name="reduced",
    )


# --------------------------------------------------------------------------- #
# Dispatch helper                                                              #
# --------------------------------------------------------------------------- #

SCAN_TYPES = ("alanine", "saturation", "reduced")


def run_scan(scan_type: str, protein: str, **kwargs) -> list[ScanVariant]:
    """Dispatch to :func:`alanine_scan`, :func:`saturation_scan` or
    :func:`reduced_alphabet_scan` by name."""
    funcs = {
        "alanine": alanine_scan,
        "saturation": saturation_scan,
        "reduced": reduced_alphabet_scan,
    }
    try:
        func = funcs[scan_type.lower()]
    except KeyError as exc:
        raise SequenceError(
            f"Unknown scan type {scan_type!r}; expected one of {SCAN_TYPES}."
        ) from exc
    return func(protein, **kwargs)
