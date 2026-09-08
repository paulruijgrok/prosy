"""Protein variants -> ordered, cloning-ready DNA fragments.

This is the pipeline that sits between :mod:`prosy.core.scan` (which produces
protein sequences) and a synthesis order. For every variant it

1. codon-optimizes the coding region under a :class:`~prosy.core.constraints.ConstraintSet`
   that forbids the cloning enzyme's recognition site,
2. adds the Golden Gate adapters derived from the destination vector,
3. pads the fragment to the vendor's minimum length,
4. verifies it: translation, absence of the enzyme's site in the coding region,
   correct excised overhangs, sliding-window GC, and - once per distinct
   protein - a full simulated assembly whose ORF is translated and compared to
   the expected fusion protein,
5. lays it out across as many plates as it needs and writes the order files.

Nothing here is nanobody- or vector-specific: the destination plasmid, enzyme,
species and plate format are all parameters.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

from prosy.core import codon as codon_mod
from prosy.core import goldengate as gg
from prosy.core import io
from prosy.core import plate as plate_mod
from prosy.core.cloning import Flanks
from prosy.core.constraints import ConstraintSet
from prosy.core.optimize import build_fragment
from prosy.core.scan import ScanVariant
from prosy.core.sequence import SYNONYMOUS_CODONS, SequenceError, translate

# --------------------------------------------------------------------------- #
# Configuration                                                                #
# --------------------------------------------------------------------------- #


@dataclass
class LibraryConfig:
    """Every knob of a fragment-building run."""

    flanks: Flanks
    enzyme: str = "BsaI"
    species: str = "e_coli"
    avoid_enzymes: list[str] = field(default_factory=list)
    min_length: int = 300
    gc_window: int | None = None      # bp; None disables the windowed GC cap
    gc_max: float = 0.72
    gc_min: float = 0.0
    plate_size: int = 96
    plate_order: str = "column"       # "column" or "row"
    seed: int = 0

    def __post_init__(self) -> None:
        if not self.avoid_enzymes:
            self.avoid_enzymes = [self.enzyme]
        elif self.enzyme not in self.avoid_enzymes:
            self.avoid_enzymes = [self.enzyme, *self.avoid_enzymes]

    def coding_constraints(self) -> ConstraintSet:
        return ConstraintSet(
            avoid_enzymes=self.avoid_enzymes,
            gc_bounds=(self.gc_min, self.gc_max) if self.gc_window else None,
            gc_window=self.gc_window,
        )

    def padding_constraints(self) -> ConstraintSet:
        return ConstraintSet(
            avoid_enzymes=self.avoid_enzymes, max_homopolymer=4,
            forbid_low_complexity=True, gc_bounds=(0.30, 0.70),
        )


# --------------------------------------------------------------------------- #
# Library members                                                              #
# --------------------------------------------------------------------------- #


@dataclass
class LibraryMember:
    """One well of the order: a variant and the DNA that will be synthesised."""

    variant: ScanVariant
    dna_coding: str = ""     # codon-optimized coding region only
    dna_final: str = ""      # coding + adapters + padding (this is what is ordered)
    plate: int = 1
    well: str = ""
    degenerate_codon: str | None = None

    @property
    def name(self) -> str:
        return self.variant.name

    @property
    def plate_well(self) -> str:
        return f"P{self.plate}:{self.well}" if self.plate > 1 else self.well


# --------------------------------------------------------------------------- #
# Building                                                                     #
# --------------------------------------------------------------------------- #


def build_library(
    variants: list[ScanVariant],
    cfg: LibraryConfig,
    *,
    backend=None,
    progress: bool = False,
) -> list[LibraryMember]:
    """Codon-optimize, flank and pad every variant.

    Variants that share a protein sequence share one optimization, so an
    alanine scan of *N* positions costs *N + 1* optimizations rather than
    re-solving identical sequences.
    """
    if backend is None or isinstance(backend, str):
        backend = codon_mod.get_backend(backend, seed=cfg.seed)

    coding_constraints = cfg.coding_constraints()
    padding_constraints = cfg.padding_constraints()

    cache: dict[str, tuple[str, str]] = {}
    members: list[LibraryMember] = []
    for i, variant in enumerate(variants):
        cached = cache.get(variant.protein)
        if cached is None:
            frag = build_fragment(
                variant.protein, species=cfg.species, constraints=coding_constraints,
                flanks=cfg.flanks, min_length=cfg.min_length,
                pad_constraints=padding_constraints, backend=backend,
                seed=cfg.seed + i,
            )
            cached = (frag.coding, frag.final)
            cache[variant.protein] = cached
        members.append(LibraryMember(variant=variant, dna_coding=cached[0],
                                     dna_final=cached[1]))
        if progress and (i + 1) % 25 == 0:
            print(f"  optimized {i + 1}/{len(variants)}", flush=True)
    assign_wells(members, cfg)
    return members


def build_degenerate_library(
    parent: ScanVariant,
    positions: list[int],
    cfg: LibraryConfig,
    *,
    codon: str = "NNK",
    backend=None,
) -> list[LibraryMember]:
    """One fragment per position, carrying a degenerate codon at that position.

    The parent CDS is optimized once; each fragment is that CDS with the target
    codon replaced by ``codon``. Every concrete codon the degenerate symbol can
    take is checked against the forbidden patterns, so no library member can
    contain the cloning enzyme's site.
    """
    from prosy.core.cloning import pad_to_length
    from prosy.core.scan import expand_degenerate

    if backend is None or isinstance(backend, str):
        backend = codon_mod.get_backend(backend, seed=cfg.seed)

    base = build_fragment(
        parent.protein, species=cfg.species, constraints=cfg.coding_constraints(),
        flanks=cfg.flanks, min_length=0, backend=backend, seed=cfg.seed,
    )
    forbidden = cfg.coding_constraints().patterns()
    concrete = expand_degenerate(codon)

    members: list[LibraryMember] = []
    for pos in positions:
        start = (pos - 1) * 3
        cds = base.coding
        conflicts = _degenerate_conflicts(cds, start, concrete, forbidden)
        if conflicts:
            repaired = _repair_degenerate_context(
                cds, start, concrete, forbidden, parent.protein)
            if repaired is None:
                raise gg.CloningError(
                    f"Degenerate codon {codon} at position {pos} can create "
                    f"forbidden pattern(s) {sorted(conflicts)}, and no synonymous "
                    "change to the neighbouring codons removes the risk. Exclude "
                    "this position or use a narrower degenerate codon."
                )
            cds = repaired
        coding = cds[:start] + codon.upper() + cds[start + 3 :]
        final = pad_to_length(cfg.flanks.apply(coding), cfg.min_length,
                              constraints=cfg.padding_constraints(),
                              species=cfg.species, seed=cfg.seed + pos)
        wt = parent.protein[pos - 1]
        members.append(LibraryMember(
            variant=ScanVariant(
                name=f"{parent.name}_{wt}{pos}{codon.upper()}",
                parent=parent.name, protein=parent.protein,
                mutations=(), scan="degenerate",
                meta={"position": pos, "wild_type": wt, "codon": codon.upper()},
            ),
            dna_coding=coding, dna_final=final, degenerate_codon=codon.upper(),
        ))
    assign_wells(members, cfg)
    return members


def _degenerate_conflicts(
    cds: str, start: int, concrete: list[str], forbidden: list[str]
) -> set[str]:
    """Forbidden patterns any concrete form of the degenerate codon could create.

    Only the neighbourhood of the codon is inspected: a pattern the degenerate
    codon can complete must overlap it, and no registered pattern is longer
    than the +/-10 bp window used here.
    """
    bad: set[str] = set()
    for alt in concrete:
        window = (cds[:start] + alt + cds[start + 3 :])[max(0, start - 10) : start + 13]
        bad |= {p for p in forbidden if p in window}
    return bad


def _repair_degenerate_context(
    cds: str, start: int, concrete: list[str], forbidden: list[str], protein: str
) -> str | None:
    """Silence a degenerate-codon conflict with one synonymous neighbour swap.

    Tries the immediate neighbours outwards; returns the repaired CDS, or
    ``None`` if no single synonymous change clears every concrete expansion
    without introducing a forbidden pattern elsewhere.
    """
    codon_index = start // 3
    for delta in (-1, 1, -2, 2, -3, 3):
        ci = codon_index + delta
        if not 0 <= ci < len(protein):
            continue
        offset = ci * 3
        current = cds[offset : offset + 3]
        for alt in SYNONYMOUS_CODONS[protein[ci]]:
            if alt == current:
                continue
            trial = cds[:offset] + alt + cds[offset + 3 :]
            if any(p in trial for p in forbidden):
                continue
            if not _degenerate_conflicts(trial, start, concrete, forbidden):
                return trial
    return None


def assign_wells(members: list[LibraryMember], cfg: LibraryConfig) -> None:
    """Fill plates in order, spilling onto plate 2, 3, ... as needed."""
    fmt = plate_mod.resolve_plate_format(cfg.plate_size)
    wells = fmt.wells(cfg.plate_order)
    for i, member in enumerate(members):
        member.plate = i // len(wells) + 1
        member.well = wells[i % len(wells)]


# --------------------------------------------------------------------------- #
# Verification                                                                 #
# --------------------------------------------------------------------------- #


def window_gc_range(dna: str, window: int) -> tuple[float, float]:
    """(min, max) GC fraction over sliding windows; whole-sequence if shorter."""
    dna = dna.upper()
    n = len(dna)
    if n == 0:
        return 0.0, 0.0
    if n <= window:
        gc = (dna.count("G") + dna.count("C")) / n
        return gc, gc
    lo, hi = 1.0, 0.0
    gc = sum(1 for b in dna[:window] if b in "GC")
    lo = hi = gc / window
    for i in range(window, n):
        gc += (dna[i] in "GC") - (dna[i - window] in "GC")
        frac = gc / window
        lo, hi = min(lo, frac), max(hi, frac)
    return lo, hi


def verify_library(
    members: list[LibraryMember],
    cfg: LibraryConfig,
    *,
    ends: gg.DestinationEnds | None = None,
    destination: str | None = None,
    expected_protein=None,
    assembly_checks: int | None = 3,
) -> list[str]:
    """Return a list of problems; an empty list means every check passed.

    ``expected_protein`` is an optional ``(member) -> str`` callable giving the
    full fusion protein the assembled plasmid should express. When it and
    ``destination`` are supplied, up to ``assembly_checks`` members (``None`` =
    all) are assembled in silico and their ORF translated and compared.
    """
    from prosy.core import enzymes as enzymes_mod

    problems: list[str] = []
    avoided = enzymes_mod.patterns_for(cfg.avoid_enzymes)
    gc_tol = 1e-9
    seen: set[tuple[int, str]] = set()

    for m in members:
        tag = m.name
        coding = m.dna_coding.upper()
        final = m.dna_final.upper()

        if m.degenerate_codon is None:
            if translate(coding) != m.variant.protein:
                problems.append(f"{tag}: coding DNA does not translate to the variant protein.")
        else:
            pos = int(m.variant.meta["position"])
            cut = (pos - 1) * 3
            if coding[cut : cut + 3] != m.degenerate_codon:
                problems.append(f"{tag}: degenerate codon not at position {pos}.")
            if (translate(coding[:cut]) != m.variant.protein[: pos - 1]
                    or translate(coding[cut + 3 :]) != m.variant.protein[pos:]):
                problems.append(
                    f"{tag}: sequence flanking the degenerate codon no longer "
                    "translates to the parent protein.")
        for pattern in avoided:
            if pattern in coding:
                problems.append(f"{tag}: forbidden site {pattern} inside the coding region.")
        if coding not in final:
            problems.append(f"{tag}: coding region not preserved in the final fragment.")
        if cfg.flanks.five_prime.upper() not in final or cfg.flanks.three_prime.upper() not in final:
            problems.append(f"{tag}: adapters not intact in the final fragment.")
        if len(final) < cfg.min_length:
            problems.append(f"{tag}: length {len(final)} < minimum {cfg.min_length}.")
        if cfg.gc_window:
            lo, hi = window_gc_range(final, cfg.gc_window)
            if hi > cfg.gc_max + gc_tol:
                problems.append(
                    f"{tag}: max {cfg.gc_window}bp-window GC {hi * 100:.1f}% > "
                    f"cap {cfg.gc_max * 100:.1f}%.")
            if cfg.gc_min > 0 and lo < cfg.gc_min - gc_tol:
                problems.append(
                    f"{tag}: min {cfg.gc_window}bp-window GC {lo * 100:.1f}% < "
                    f"floor {cfg.gc_min * 100:.1f}%.")
        if ends is not None and m.degenerate_codon is None:
            problems.extend(f"{tag}: {p}" for p in gg.check_fragment(final, ends, cfg.enzyme))
        key = (m.plate, m.well)
        if key in seen:
            problems.append(f"{tag}: duplicate well {m.plate_well}.")
        seen.add(key)

    if destination and expected_protein is not None:
        subset = members if assembly_checks is None else members[:assembly_checks]
        for m in subset:
            if m.degenerate_codon is not None:
                continue
            try:
                product = gg.assemble(destination, [m.dna_final], cfg.enzyme)
            except gg.CloningError as exc:
                problems.append(f"{m.name}: simulated assembly failed - {exc}")
                continue
            want = expected_protein(m)
            got = orf_protein(product)
            if got is None or want not in got:
                problems.append(
                    f"{m.name}: assembled ORF does not contain the expected protein.")
    return problems


def orf_protein(product: str, *, start_motif: str = "ATG") -> str | None:
    """Translate the longest ORF of a circular product starting at ``start_motif``.

    The product is a circular top strand written from an arbitrary origin, so
    every frame is tried on the doubled sequence and the longest stop-terminated
    peptide is returned.
    """
    seq = product.upper()
    doubled = seq + seq
    best: str | None = None
    start = 0
    while True:
        i = doubled.find(start_motif, start)
        if i == -1 or i >= len(seq):
            break
        start = i + 1
        chunk = doubled[i : i + 3 * (len(seq) // 3)]
        chunk = chunk[: len(chunk) // 3 * 3]
        try:
            peptide = translate(chunk, to_stop=True)
        except SequenceError:
            # A frame that runs off the doubled sequence with a partial codon;
            # not a real ORF start, so move on.
            continue
        if best is None or len(peptide) > len(best):
            best = peptide
    return best


# --------------------------------------------------------------------------- #
# Outputs                                                                      #
# --------------------------------------------------------------------------- #

MAPPING_COLUMNS = [
    "Plate", "Well", "Name", "Parent", "Scan", "Mutation", "Position",
    "WT", "Mut", "Category", "Protein_len", "Final_len", "GC",
    "Protein", "Coding_DNA", "Final_DNA",
]


def write_library(
    members: list[LibraryMember],
    out_dir: str | Path,
    stamp: str,
    cfg: LibraryConfig,
    *,
    write_plates: bool = True,
) -> dict:
    """Write the mapping CSV, the design CSV and one vendor sheet per plate."""
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    result: dict = {"plates": []}

    rows = []
    for m in members:
        v = m.variant
        mut = v.mutations[0] if v.mutations else None
        gc = (m.dna_final.upper().count("G") + m.dna_final.upper().count("C"))
        rows.append({
            "Plate": m.plate, "Well": m.well, "Name": m.name, "Parent": v.parent,
            "Scan": v.scan, "Mutation": v.mutation_label,
            "Position": mut.position if mut else v.meta.get("position", ""),
            "WT": mut.wt if mut else v.meta.get("wild_type", ""),
            "Mut": mut.mut if mut else (m.degenerate_codon or ""),
            "Category": v.category or "", "Protein_len": len(v.protein),
            "Final_len": len(m.dna_final),
            "GC": f"{gc / len(m.dna_final):.3f}" if m.dna_final else "",
            "Protein": v.protein, "Coding_DNA": m.dna_coding,
            "Final_DNA": m.dna_final,
        })
    mapping_path = out_dir / f"{stamp}_library_mapping.csv"
    io.write_csv(rows, mapping_path, MAPPING_COLUMNS)
    result["mapping"] = mapping_path

    design_rows = [{
        "Name": m.name, "Sequence": m.variant.protein,
        "5' nucleotides": cfg.flanks.five_prime,
        "3' nucleotides": cfg.flanks.three_prime,
    } for m in members]
    design_path = out_dir / f"{stamp}_design_seqs.csv"
    io.write_csv(design_rows, design_path, io.DESIGN_COLUMNS)
    result["design"] = design_path

    if write_plates:
        by_plate: dict[int, list[LibraryMember]] = {}
        for m in members:
            by_plate.setdefault(m.plate, []).append(m)
        for plate_no, group in sorted(by_plate.items()):
            suffix = "" if len(by_plate) == 1 else f"_p{plate_no}"
            path = out_dir / f"{stamp}_plate{suffix}.xlsx"
            io.write_plate_xlsx(
                [{"well": m.well, "name": m.name, "sequence": m.dna_final} for m in group],
                path,
            )
            result["plates"].append(path)
    return result


# --------------------------------------------------------------------------- #
# Small helpers                                                                #
# --------------------------------------------------------------------------- #


def naive_cds(protein: str) -> str:
    """Most-frequent-codon CDS with no constraint solving (tests, previews)."""
    return "".join(SYNONYMOUS_CODONS[aa][0] for aa in protein.upper())
