"""Type IIS (Golden Gate) digestion, adapter design and assembly simulation.

Everything here works in **top-strand coordinates**, which makes sticky ends
easy to reason about:

    a Type IIS cut at top-strand index ``p`` severs the top strand before ``p``
    and the bottom strand before ``p + k`` (``k`` = overhang length).

Consequently the fragment lying to the *right* of ``p`` owns the ``k``-base
overhang ``seq[p : p + k]`` in its top strand, and the fragment to the *left*
of ``p`` does not contain those bases at all. Ligation is therefore plain
concatenation of top strands, with the requirement that the right-hand sticky
end of one fragment equals the left-hand sticky end of the next. That single
rule drives :func:`digest`, :func:`assemble` and :func:`design_insert_flanks`.

Typical use - work out what a destination plasmid needs, build matching flanks,
then prove the design by simulating the reaction:

>>> ends = destination_ends(FP01, "BsaI")          # doctest: +SKIP
>>> flanks = design_insert_flanks(ends, "BsaI")    # doctest: +SKIP
>>> product = assemble(FP01, [fragment], "BsaI")   # doctest: +SKIP
"""

from __future__ import annotations

from dataclasses import dataclass

from prosy.core import enzymes as enzymes_mod
from prosy.core.cloning import Flanks
from prosy.core.enzymes import Enzyme
from prosy.core.sequence import SequenceError, is_dna, reverse_complement


class CloningError(ValueError):
    """Raised for impossible or ambiguous cloning designs."""


# --------------------------------------------------------------------------- #
# Sites and cuts                                                               #
# --------------------------------------------------------------------------- #


@dataclass(frozen=True)
class Cut:
    """A Type IIS cut, in top-strand coordinates."""

    position: int      # top strand is cut immediately before this index
    overhang: str      # the k-base overhang, written 5'->3' on the top strand
    site_start: int    # 0-based start of the recognition site on the top strand
    strand: int        # +1 if the site reads on the top strand, -1 on the bottom

    @property
    def bottom_position(self) -> int:
        return self.position + len(self.overhang)


def _resolve(enzyme: str | Enzyme) -> Enzyme:
    enz = enzyme if isinstance(enzyme, Enzyme) else enzymes_mod.get(enzyme)
    if enz.cut_offset_top is None or enz.overhang_length is None:
        raise CloningError(
            f"{enz.name} has no Type IIS cut geometry registered; "
            "it cannot be used for Golden Gate."
        )
    return enz


def find_cuts(seq: str, enzyme: str | Enzyme, *, circular: bool = False) -> list[Cut]:
    """All Type IIS cuts ``enzyme`` makes in ``seq``, sorted by position.

    Sites whose cut would fall off the end of a linear sequence are ignored
    (the enzyme cannot cut what is not there). On a circular sequence cut
    positions wrap around the origin.
    """
    seq = _clean(seq)
    enz = _resolve(enzyme)
    n = len(seq)
    site, site_rc = enz.site.upper(), enz.site_rc.upper()
    k, offset = enz.overhang_length, enz.cut_offset_top

    search = seq + seq[: len(site) - 1] if circular else seq
    cuts: list[Cut] = []
    for strand, pattern in ((+1, site), (-1, site_rc)):
        start = 0
        while True:
            i = search.find(pattern, start)
            if i == -1 or i >= n:
                break
            start = i + 1
            if strand == +1:
                p = i + len(pattern) + offset
            else:
                p = i - offset - k
            if circular:
                p %= n
            elif not 0 <= p <= n - k:
                continue
            cuts.append(Cut(position=p, overhang=_slice(seq, p, k, circular),
                            site_start=i, strand=strand))
    cuts.sort(key=lambda c: c.position)
    return cuts


def _clean(seq: str) -> str:
    seq = "".join(seq.split()).upper()
    if not is_dna(seq):
        raise SequenceError("Sequence must be non-empty DNA (ACGT only).")
    return seq


def _slice(seq: str, start: int, length: int, circular: bool) -> str:
    if start + length <= len(seq):
        return seq[start : start + length]
    if not circular:
        raise CloningError("Cut runs off the end of a linear sequence.")
    return (seq + seq)[start : start + length]


# --------------------------------------------------------------------------- #
# Digestion                                                                    #
# --------------------------------------------------------------------------- #


@dataclass(frozen=True)
class DigestFragment:
    """One product of a Type IIS digest.

    ``top`` is the fragment's own top strand. It *starts with*
    ``left_overhang`` (the overhang it owns) and *stops before*
    ``right_overhang`` (which the neighbouring fragment owns). ``None`` marks an
    original end of a linear substrate rather than a cut end.
    """

    top: str
    left_overhang: str | None
    right_overhang: str | None
    has_site: bool

    def __len__(self) -> int:  # noqa: D105
        return len(self.top)


def digest(seq: str, enzyme: str | Enzyme, *, circular: bool = False) -> list[DigestFragment]:
    """Digest ``seq`` with a Type IIS ``enzyme``.

    A circular substrate with *n* cuts yields *n* fragments; a linear substrate
    yields *n + 1* (the two outermost keep one uncut, ``None``-overhang end).
    Returns ``[whole sequence]`` when the enzyme does not cut.
    """
    seq = _clean(seq)
    enz = _resolve(enzyme)
    cuts = find_cuts(seq, enz, circular=circular)
    if not cuts:
        return [DigestFragment(top=seq, left_overhang=None, right_overhang=None,
                               has_site=False)]

    positions = [c.position for c in cuts]
    overhangs = {c.position: c.overhang for c in cuts}
    frags: list[DigestFragment] = []

    if circular:
        for idx, p in enumerate(positions):
            q = positions[(idx + 1) % len(positions)]
            top = seq[p:q] if q > p else seq[p:] + seq[:q]
            frags.append(DigestFragment(
                top=top, left_overhang=overhangs[p], right_overhang=overhangs[q],
                has_site=_has_site(top, enz)))
    else:
        bounds = [0, *positions, len(seq)]
        for idx in range(len(bounds) - 1):
            p, q = bounds[idx], bounds[idx + 1]
            top = seq[p:q]
            if not top:
                continue
            frags.append(DigestFragment(
                top=top,
                left_overhang=overhangs.get(p),
                right_overhang=overhangs.get(q),
                has_site=_has_site(top, enz)))
    return frags


def _has_site(seq: str, enz: Enzyme) -> bool:
    return any(p in seq.upper() for p in enz.patterns())


def count_sites(seq: str, enzyme: str | Enzyme, *, circular: bool = False) -> int:
    """Number of recognition sites for ``enzyme`` (both strands)."""
    seq = _clean(seq)
    enz = enzyme if isinstance(enzyme, Enzyme) else enzymes_mod.get(enzyme)
    search = seq + seq[: len(enz.site) - 1] if circular else seq
    total = 0
    for pattern in enz.patterns():
        start = 0
        while True:
            i = search.find(pattern, start)
            if i == -1 or i >= len(seq):
                break
            total += 1
            start = i + 1
    return total


# --------------------------------------------------------------------------- #
# Destination analysis                                                         #
# --------------------------------------------------------------------------- #


@dataclass(frozen=True)
class DestinationEnds:
    """What a linearised destination vector demands of an insert.

    ``five_overhang`` is the overhang the insert must carry at its 5' end and
    ``three_overhang`` the one at its 3' end, both written 5'->3' on the
    insert's top strand. ``upstream``/``downstream`` are the vector sequences
    immediately flanking the insertion point, which let a caller check reading
    frame and fusion junctions.
    """

    enzyme: str
    backbone: str            # backbone top strand, starting at its own overhang
    five_overhang: str
    three_overhang: str
    upstream: str            # last 60 bp of vector before the insert's 5' overhang
    downstream: str          # first 60 bp of vector after the insert's 3' overhang
    dropout: str             # the released stuffer (top strand), '' if none


def destination_ends(
    plasmid: str,
    enzyme: str | Enzyme,
    *,
    circular: bool = True,
    context: int = 60,
) -> DestinationEnds:
    """Analyse a Golden Gate destination vector.

    The backbone is the digest fragment that no longer carries a recognition
    site (the sites leave with the drop-out stuffer, which is what makes the
    reaction directional and self-selecting). Fails loudly if the plasmid does
    not have exactly one such fragment.
    """
    enz = _resolve(enzyme)
    frags = digest(plasmid, enz, circular=circular)
    if len(frags) < 2:
        raise CloningError(
            f"{enz.name} makes {len(frags) - 1} cut(s) in the destination; "
            "a Golden Gate destination needs two outward-facing sites."
        )
    backbones = [f for f in frags if not f.has_site]
    if len(backbones) != 1:
        raise CloningError(
            f"Expected exactly one {enz.name}-free backbone fragment, found "
            f"{len(backbones)} (fragment sizes: {[len(f) for f in frags]}). "
            "Check that the destination's sites face outward from the stuffer."
        )
    backbone = backbones[0]
    if backbone.left_overhang is None or backbone.right_overhang is None:
        raise CloningError("Backbone fragment has an uncut end; is the plasmid circular?")
    dropout = "".join(f.top for f in frags if f is not backbone)

    # The insert's 5' end ligates to the backbone's right-hand sticky end.
    return DestinationEnds(
        enzyme=enz.name,
        backbone=backbone.top,
        five_overhang=backbone.right_overhang,
        three_overhang=backbone.left_overhang,
        upstream=backbone.top[-context:],
        downstream=backbone.top[:context],
        dropout=dropout,
    )


# --------------------------------------------------------------------------- #
# Adapter design                                                               #
# --------------------------------------------------------------------------- #

# Bases tried as the single spacer nucleotide between site and overhang.
_SPACER_CHOICES = ("A", "T", "C", "G")
# Short outer buffer so the enzyme has room to bind at a fragment end.
DEFAULT_OUTER_5 = "TGTATC"
DEFAULT_OUTER_3 = "TCTAGT"


def design_insert_flanks(
    ends: DestinationEnds,
    enzyme: str | Enzyme | None = None,
    *,
    outer_5: str = DEFAULT_OUTER_5,
    outer_3: str = DEFAULT_OUTER_3,
    spacer: str | None = None,
) -> Flanks:
    """Build the 5'/3' adapters an insert needs for ``ends``.

    Layout of the ordered fragment::

        outer_5 + SITE + N + five_overhang + <coding> + three_overhang + N + rc(SITE) + outer_3

    The overhangs are part of the adapters, so ``<coding>`` is exactly the
    sequence that sits between the two vector junctions. The spacer base is
    chosen so the adapters introduce no additional recognition site; pass
    ``spacer`` to pin it.
    """
    enz = _resolve(enzyme or ends.enzyme)
    site = enz.site.upper()
    site_rc = reverse_complement(site)

    for base in ([spacer.upper()] if spacer else _SPACER_CHOICES):
        five = f"{outer_5}{site}{base}{ends.five_overhang}".upper()
        three = f"{ends.three_overhang}{base}{site_rc}{outer_3}".upper()
        if _count_patterns(five, enz) == 1 and _count_patterns(three, enz) == 1:
            return Flanks(five_prime=five, three_prime=three)
    raise CloningError(
        f"No spacer base leaves exactly one {enz.name} site per adapter; "
        "supply outer buffers that do not complete a recognition site."
    )


def _count_patterns(seq: str, enz: Enzyme) -> int:
    return sum(seq.upper().count(p) for p in enz.patterns())


def check_fragment(
    fragment: str,
    ends: DestinationEnds,
    enzyme: str | Enzyme | None = None,
) -> list[str]:
    """Return the reasons ``fragment`` is not a valid insert (empty = valid)."""
    enz = _resolve(enzyme or ends.enzyme)
    problems: list[str] = []
    frags = digest(fragment, enz, circular=False)
    usable = [f for f in frags
              if f.left_overhang is not None and f.right_overhang is not None]
    if len(usable) != 1:
        problems.append(
            f"Fragment yields {len(usable)} excisable piece(s) with {enz.name}; "
            "expected exactly one (two inward-facing sites, none in the coding region)."
        )
        return problems
    core = usable[0]
    if core.has_site:
        problems.append(f"Excised insert still contains a {enz.name} site.")
    if core.left_overhang != ends.five_overhang:
        problems.append(
            f"5' overhang {core.left_overhang} != required {ends.five_overhang}.")
    if core.right_overhang != ends.three_overhang:
        problems.append(
            f"3' overhang {core.right_overhang} != required {ends.three_overhang}.")
    return problems


# --------------------------------------------------------------------------- #
# Assembly simulation                                                          #
# --------------------------------------------------------------------------- #


def assemble(
    destination: str,
    inserts: list[str],
    enzyme: str | Enzyme,
    *,
    circular_destination: bool = True,
) -> str:
    """Simulate the Golden Gate reaction and return the circular product.

    The product is returned as a top strand starting at the backbone's own
    sticky end. Raises :class:`CloningError` if the overhangs do not chain into
    a single closed circle - which is exactly the failure mode a design review
    is meant to catch.
    """
    enz = _resolve(enzyme)
    dest_frags = [f for f in digest(destination, enz, circular=circular_destination)
                  if not f.has_site]
    if len(dest_frags) != 1:
        raise CloningError(
            f"Destination gave {len(dest_frags)} site-free fragments; expected 1.")

    pieces: list[DigestFragment] = list(dest_frags)
    for idx, ins in enumerate(inserts):
        usable = [f for f in digest(ins, enz, circular=False)
                  if f.left_overhang is not None and f.right_overhang is not None
                  and not f.has_site]
        if len(usable) != 1:
            raise CloningError(
                f"Insert {idx + 1} gave {len(usable)} usable fragments; expected 1.")
        pieces.append(usable[0])

    chain = [pieces[0]]
    remaining = pieces[1:]
    while remaining:
        want = chain[-1].right_overhang
        matches = [f for f in remaining if f.left_overhang == want]
        if not matches:
            raise CloningError(
                f"No fragment starts with overhang {want}; available: "
                f"{[f.left_overhang for f in remaining]}.")
        if len(matches) > 1:
            raise CloningError(f"Overhang {want} is ambiguous ({len(matches)} matches).")
        chain.append(matches[0])
        remaining.remove(matches[0])

    if chain[-1].right_overhang != chain[0].left_overhang:
        raise CloningError(
            f"Assembly does not circularise: {chain[-1].right_overhang} != "
            f"{chain[0].left_overhang}.")
    product = "".join(f.top for f in chain)
    if _has_site(product, enz):
        raise CloningError(f"Assembled product still contains a {enz.name} site.")
    return product
