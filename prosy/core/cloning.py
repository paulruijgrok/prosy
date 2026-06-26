"""Cloning helpers: append flanking adapters and pad fragments to a minimum
length. Mirrors the flank/padding behaviour of the original script, but the
user's sequence (coding region + flanks) is always preserved verbatim - only
generated padding is constraint-checked.
"""

from __future__ import annotations

import random
from dataclasses import dataclass

from prosy.core import codon as codon_mod
from prosy.core.constraints import ConstraintSet
from prosy.core.sequence import SequenceError


@dataclass(frozen=True)
class Flanks:
    """5' and 3' adapter sequences added around a coding region."""

    five_prime: str = ""
    three_prime: str = ""

    def apply(self, coding: str) -> str:
        return self.five_prime.upper() + coding.upper() + self.three_prime.upper()


def add_flanks(coding: str, flanks: Flanks) -> str:
    return flanks.apply(coding)


def generate_padding(
    length: int,
    *,
    constraints: ConstraintSet | None = None,
    species: str = "e_coli",
    seed: int | None = None,
) -> str:
    """Generate a benign padding sequence of ``length`` bp satisfying constraints.

    Padding is meant to sit outside the cloning adapters purely to reach a
    vendor minimum length, so it avoids restriction sites, homopolymers,
    low-complexity repeats and extreme GC by default.
    """
    if length <= 0:
        return ""
    if constraints is None:
        constraints = ConstraintSet(
            avoid_enzymes=["BsaI"],
            max_homopolymer=4,
            forbid_low_complexity=True,
            gc_bounds=(0.30, 0.70),
        )
    avoid = [p.upper() for p in constraints.patterns()]
    rng = random.Random(seed)

    bases = "ATCG"
    out: list[str] = []
    max_pat = max((len(p) for p in avoid), default=1)
    attempts = 0
    max_attempts = length * 200 + 1000
    while len(out) < length:
        attempts += 1
        if attempts > max_attempts:
            raise SequenceError("Could not generate constraint-satisfying padding.")
        ch = rng.choice(bases)
        tail = "".join(out[-(max_pat - 1):]) + ch
        if any(p in tail for p in avoid):
            continue
        # GC guard near the end so we don't paint into a corner.
        out.append(ch)
    seq = "".join(out)
    if constraints.gc_bounds is not None:
        seq = _fix_padding_gc(seq, avoid, constraints.gc_bounds, rng)
    return seq


def _fix_padding_gc(seq, avoid, gc_bounds, rng, max_iter=20000):
    from prosy.core.sequence import gc_content

    lo, hi = gc_bounds
    s = list(seq)
    for _ in range(max_iter):
        gc = gc_content("".join(s))
        if lo <= gc <= hi:
            return "".join(s)
        want_gc = gc < lo
        i = rng.randrange(len(s))
        cur = s[i]
        cur_is_gc = cur in "GC"
        if want_gc == cur_is_gc:
            continue
        for repl in ("G", "C") if want_gc else ("A", "T"):
            trial = s[:]
            trial[i] = repl
            window = "".join(trial[max(0, i - 6): i + 7])
            if not any(p in window for p in avoid):
                s = trial
                break
    return "".join(s)


def pad_to_length(
    sequence: str,
    min_length: int,
    *,
    constraints: ConstraintSet | None = None,
    species: str = "e_coli",
    seed: int | None = None,
) -> str:
    """Pad ``sequence`` symmetrically to reach ``min_length`` bp.

    The input ``sequence`` (already containing coding region + flanks) is never
    modified; padding is generated independently for each end.
    """
    current = len(sequence)
    if min_length <= 0 or current >= min_length:
        return sequence
    needed = min_length - current
    left = needed // 2
    right = needed - left
    rng = random.Random(seed)
    pad5 = generate_padding(left, constraints=constraints, species=species,
                            seed=rng.randint(0, 2**31 - 1)) if left else ""
    pad3 = generate_padding(right, constraints=constraints, species=species,
                            seed=rng.randint(0, 2**31 - 1)) if right else ""
    return pad5 + sequence + pad3
