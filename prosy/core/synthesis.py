"""Synthesis manufacturability: the checks a DNA vendor runs on an order.

Vendors reject or surcharge sequences on a "complexity score" built from a few
measurable properties. The two that dominate for codon-optimized designs are

* **repeat content** - the fraction of bases belonging to any k-mer (k = 8 by
  convention) that occurs more than once in the sequence, *counting both
  strands*. A maximally codon-optimized CDS is pathological here: always
  picking the single most frequent synonymous codon means every Ala is ``GCG``
  and every Leu ``CTG``, so 8-mers recur constantly. Coverage of 70%+ is normal
  for ``use_best_codon`` output and is typically rejected above 40%.
* **local repeat density** - the same measure taken over a sliding window
  (90 bp by convention), which catches a locally repetitive patch inside an
  otherwise acceptable sequence.

plus windowed GC extremes, homopolymers, and overall GC.

:data:`DEFAULT_SPEC` encodes a common vendor rule set. Everything is a
parameter, so a different vendor's thresholds are a different
:class:`SynthesisSpec` rather than a different code path.

These functions *measure*; they do not design. Use them to gate an order:
build with :mod:`prosy.core.codon` constraints, then prove the result here.
"""

from __future__ import annotations

import re
from collections import Counter
from dataclasses import dataclass, field

from prosy.core.sequence import reverse_complement


@dataclass(frozen=True)
class SynthesisSpec:
    """A vendor's acceptance thresholds."""

    repeat_k: int = 8                     # k-mer size defining a "repeat"
    max_repeat_fraction: float = 0.40     # of the whole sequence
    repeat_window: int = 90               # bp, for local repeat density
    max_repeat_window_fraction: float = 0.88
    gc_windows: tuple[tuple[int, float, float], ...] = ((20, 0.00, 0.90),)
    gc_bounds: tuple[float, float] = (0.25, 0.68)
    max_homopolymer: dict[str, int] = field(
        default_factory=lambda: {"A": 8, "T": 8, "G": 5, "C": 5})
    include_reverse_complement: bool = True


#: A common vendor rule set (the one the 260917 order was scored against).
DEFAULT_SPEC = SynthesisSpec()


# --------------------------------------------------------------------------- #
# Repeat measurement                                                           #
# --------------------------------------------------------------------------- #


def repeat_coverage(
    seq: str, k: int = 8, *, include_reverse_complement: bool = True
) -> bytearray:
    """Per-base mask: 1 where the base belongs to a k-mer seen more than once.

    With ``include_reverse_complement`` (the vendor default) a k-mer and its
    reverse complement count as the same repeat, because they anneal.
    """
    seq = seq.upper()
    n = len(seq)
    mask = bytearray(n)
    if n < k:
        return mask
    counts: Counter[str] = Counter()
    for i in range(n - k + 1):
        kmer = seq[i : i + k]
        counts[kmer] += 1
        if include_reverse_complement:
            counts[reverse_complement(kmer)] += 1
    for i in range(n - k + 1):
        if counts[seq[i : i + k]] > 1:
            for j in range(i, i + k):
                mask[j] = 1
    return mask


def repeat_fraction(seq: str, k: int = 8, *, include_reverse_complement: bool = True) -> float:
    """Fraction of the sequence that sits inside a repeated k-mer."""
    if not seq:
        return 0.0
    mask = repeat_coverage(seq, k, include_reverse_complement=include_reverse_complement)
    return sum(mask) / len(mask)


def densest_repeat_window(mask: bytearray, window: int) -> tuple[int, float]:
    """``(0-based start, fraction)`` of the window with most repeat bases."""
    n = len(mask)
    if n == 0:
        return 0, 0.0
    if n <= window:
        return 0, sum(mask) / n
    running = sum(mask[:window])
    best_i, best = 0, running / window
    for i in range(window, n):
        running += mask[i] - mask[i - window]
        if running / window > best:
            best_i, best = i - window + 1, running / window
    return best_i, best


def extreme_gc_window(seq: str, window: int) -> tuple[int, float, int, float]:
    """``(argmin, min, argmax, max)`` GC fraction over sliding windows."""
    seq = seq.upper()
    n = len(seq)
    if n == 0:
        return 0, 0.0, 0, 0.0
    if n <= window:
        gc = sum(1 for b in seq if b in "GC") / n
        return 0, gc, 0, gc
    running = sum(1 for b in seq[:window] if b in "GC")
    lo_i = hi_i = 0
    lo = hi = running / window
    for i in range(window, n):
        running += (seq[i] in "GC") - (seq[i - window] in "GC")
        frac = running / window
        if frac < lo:
            lo_i, lo = i - window + 1, frac
        if frac > hi:
            hi_i, hi = i - window + 1, frac
    return lo_i, lo, hi_i, hi


def longest_homopolymers(seq: str) -> dict[str, int]:
    """Longest single-base run per base present in ``seq``."""
    out: dict[str, int] = {}
    for run in re.finditer(r"(A+|C+|G+|T+)", seq.upper()):
        base = run.group()[0]
        out[base] = max(out.get(base, 0), len(run.group()))
    return out


# --------------------------------------------------------------------------- #
# Report                                                                       #
# --------------------------------------------------------------------------- #


@dataclass
class SynthesisReport:
    """Measured manufacturability of one sequence against a spec."""

    length: int
    gc: float
    repeat_fraction: float
    repeat_window_start: int          # 1-based, as vendors report positions
    repeat_window_fraction: float
    gc_window_extremes: dict          # window -> (min_pos, min, max_pos, max)
    homopolymers: dict[str, int]
    problems: list[str]

    @property
    def ok(self) -> bool:
        return not self.problems

    def summary(self) -> str:
        parts = [f"{self.length} bp", f"GC {self.gc * 100:.1f}%",
                 f"repeat8 {self.repeat_fraction * 100:.1f}%",
                 f"dens {self.repeat_window_fraction * 100:.1f}%"]
        for w, (_, _, _, hi) in sorted(self.gc_window_extremes.items()):
            parts.append(f"GC{w}max {hi * 100:.0f}%")
        return "  ".join(parts)


def check(seq: str, spec: SynthesisSpec = DEFAULT_SPEC) -> SynthesisReport:
    """Measure ``seq`` against ``spec`` and list every threshold it breaches."""
    seq = seq.upper()
    n = len(seq)
    tol = 1e-9
    mask = repeat_coverage(seq, spec.repeat_k,
                           include_reverse_complement=spec.include_reverse_complement)
    frac = sum(mask) / n if n else 0.0
    win_start, win_frac = densest_repeat_window(mask, spec.repeat_window)
    gc = (seq.count("G") + seq.count("C")) / n if n else 0.0
    extremes = {w: extreme_gc_window(seq, w) for w, _, _ in spec.gc_windows}
    homo = longest_homopolymers(seq)

    problems: list[str] = []
    if frac > spec.max_repeat_fraction + tol:
        problems.append(
            f"repeated {spec.repeat_k}-mers cover {frac * 100:.1f}% of the sequence "
            f"(limit {spec.max_repeat_fraction * 100:.0f}%)")
    if win_frac > spec.max_repeat_window_fraction + tol:
        problems.append(
            f"{win_frac * 100:.1f}% of a {spec.repeat_window} bp window at position "
            f"{win_start + 1} is repeat (limit {spec.max_repeat_window_fraction * 100:.0f}%)")
    for w, lo_lim, hi_lim in spec.gc_windows:
        lo_i, lo, hi_i, hi = extremes[w]
        if hi > hi_lim + tol:
            problems.append(
                f"{w} bp window at position {hi_i + 1} is {hi * 100:.0f}% GC "
                f"(limit {hi_lim * 100:.0f}%)")
        if lo < lo_lim - tol:
            problems.append(
                f"{w} bp window at position {lo_i + 1} is {lo * 100:.0f}% GC "
                f"(minimum {lo_lim * 100:.0f}%)")
    lo_gc, hi_gc = spec.gc_bounds
    if not lo_gc - tol <= gc <= hi_gc + tol:
        problems.append(f"overall GC {gc * 100:.1f}% outside "
                        f"{lo_gc * 100:.0f}-{hi_gc * 100:.0f}%")
    for base, limit in spec.max_homopolymer.items():
        if homo.get(base, 0) > limit:
            problems.append(f"{homo[base]}x{base} homopolymer (limit {limit})")

    return SynthesisReport(
        length=n, gc=gc, repeat_fraction=frac, repeat_window_start=win_start + 1,
        repeat_window_fraction=win_frac, gc_window_extremes=extremes,
        homopolymers=homo, problems=problems,
    )
