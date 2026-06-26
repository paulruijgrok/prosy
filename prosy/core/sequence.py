"""Amino-acid and DNA sequence primitives: validation, translation and the
genetic code. Kept dependency-free so it works in any environment.
"""

from __future__ import annotations

import re
from dataclasses import dataclass

# Standard genetic code (DNA codon -> single-letter amino acid, * = stop).
CODON_TABLE: dict[str, str] = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
    "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}

# Reverse map: amino acid -> list of synonymous codons.
SYNONYMOUS_CODONS: dict[str, list[str]] = {}
for _codon, _aa in CODON_TABLE.items():
    SYNONYMOUS_CODONS.setdefault(_aa, []).append(_codon)

AMINO_ACIDS = set("ACDEFGHIKLMNPQRSTVWY")
_DNA_RE = re.compile(r"^[ACGT]+$", re.IGNORECASE)
_AA_RE = re.compile(r"^[ACDEFGHIKLMNPQRSTVWY]+$", re.IGNORECASE)

_COMPLEMENT = str.maketrans("ACGTacgt", "TGCAtgca")


class SequenceError(ValueError):
    """Raised for malformed or inconsistent sequences."""


@dataclass(frozen=True)
class PointMutation:
    """A single residue substitution in 1-based protein coordinates.

    e.g. ``T27F`` -> wt='T', position=27, mut='F'.
    """

    wt: str
    position: int
    mut: str

    @classmethod
    def parse(cls, token: str) -> "PointMutation":
        m = re.fullmatch(r"([A-Z])(\d+)([A-Z])", token.strip().upper())
        if not m:
            raise SequenceError(f"Cannot parse mutation token: {token!r}")
        return cls(wt=m.group(1), position=int(m.group(2)), mut=m.group(3))

    def __str__(self) -> str:  # noqa: D105
        return f"{self.wt}{self.position}{self.mut}"


def is_dna(seq: str) -> bool:
    return bool(seq) and bool(_DNA_RE.match(seq))


def is_protein(seq: str) -> bool:
    return bool(seq) and bool(_AA_RE.match(seq))


def validate_protein(seq: str) -> str:
    """Return the upper-cased protein sequence or raise ``SequenceError``."""
    seq = seq.strip().upper()
    if not is_protein(seq):
        bad = sorted(set(seq) - AMINO_ACIDS)
        raise SequenceError(f"Invalid amino-acid sequence; unexpected: {bad}")
    return seq


def reverse_complement(dna: str) -> str:
    return dna.translate(_COMPLEMENT)[::-1]


def translate(dna: str, to_stop: bool = False) -> str:
    """Translate a coding DNA sequence (length must be a multiple of 3)."""
    dna = dna.upper()
    if len(dna) % 3 != 0:
        raise SequenceError("DNA length is not a multiple of 3; cannot translate.")
    out = []
    for i in range(0, len(dna), 3):
        aa = CODON_TABLE.get(dna[i : i + 3])
        if aa is None:
            raise SequenceError(f"Unknown codon {dna[i:i + 3]!r} at position {i}.")
        if aa == "*" and to_stop:
            break
        out.append(aa)
    return "".join(out)


def gc_content(dna: str) -> float:
    dna = dna.upper()
    if not dna:
        return 0.0
    return (dna.count("G") + dna.count("C")) / len(dna)


def apply_mutation(protein: str, mutation: PointMutation, *, check_wt: bool = True) -> str:
    """Return ``protein`` with ``mutation`` applied (1-based position).

    Raises ``SequenceError`` if the position is out of range or, when
    ``check_wt`` is True, the wild-type residue does not match.
    """
    protein = validate_protein(protein)
    idx = mutation.position - 1
    if not 0 <= idx < len(protein):
        raise SequenceError(
            f"Mutation {mutation} position out of range for length {len(protein)}."
        )
    if check_wt and protein[idx] != mutation.wt:
        raise SequenceError(
            f"Mutation {mutation}: expected {mutation.wt} at position "
            f"{mutation.position} but found {protein[idx]}."
        )
    return protein[:idx] + mutation.mut + protein[idx + 1 :]
