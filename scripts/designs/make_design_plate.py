#!/usr/bin/env python3
"""Turn a CSV of designed protein sequences into a plate of orderable fragments.

The design-set analogue of ``scripts/nanobodies/make_nanobody_plate.py``: read a
sequence table, lay the designs out **one design group per plate row**,
codon-optimize each one, add the Golden Gate adapters derived from a
destination plasmid, verify everything, and write the same family of outputs
(plate upload sheet, design CSV, mapping CSV, plate map).

Layout differs from the nanobody script in one way that matters here: design
groups are not all the same size, so each group gets its own row via
``layout.layout_lines()`` and a group larger than 12 spills onto as many
further rows as it needs, balanced (14 -> 7 + 7) unless ``--fill-rows`` is
given (-> 12 + 2).

Examples
--------
# The 260917 designs design set into gg002 (BsaI):
python make_design_plate.py \\
    --csv "Working folder/260917_designs_enzyme/sequences.csv" \\
    --destination "Working folder/260917_designs_enzyme/gg002.txt" \\
    --out-dir "Working folder/260917_designs_enzyme" --stamp 260918_designs

# A different column layout, 384-well, no windowed GC cap:
python make_design_plate.py --csv designs.csv --destination gg002 \\
    --id-column name --group-column family --sequence-column aa_seq \\
    --plate-size 384 --gc-window 0
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

_HERE = Path(__file__).resolve()
sys.path.insert(0, str(_HERE.parents[2]))

from prosy.core import codon as codon_mod  # noqa: E402
from prosy.core import goldengate as gg  # noqa: E402
from prosy.core import io  # noqa: E402
from prosy.core import layout as layout_mod  # noqa: E402
from prosy.core import platemap as platemap_mod  # noqa: E402
from prosy.core import synthesis  # noqa: E402
from prosy.core.cloning import Flanks  # noqa: E402
from prosy.core.constraints import ConstraintSet  # noqa: E402
from prosy.core.library import orf_protein, window_gc_range  # noqa: E402
from prosy.core.optimize import build_fragment_with_ladder  # noqa: E402
from prosy.core.sequence import gc_content, translate, validate_protein  # noqa: E402

MAPPING_COLUMNS = [
    "Well", "Row", "Order", "Name", "Design_group", "Group_index", "Protein_len",
    "Coding_len", "Final_len", "GC", "GC_win_min", "GC_win_max",
    "Repeat8_pct", "Repeat_dens90_pct", "GC20_max_pct", "Ladder_step",
    "Product_len", "Product_protein", "Protein", "Coding_DNA", "Final_DNA",
]


def read_sequence_file(path: Path) -> str:
    lines = [ln.strip() for ln in Path(path).read_text(encoding="utf-8").splitlines()
             if ln.strip() and not ln.startswith(">")]
    return "".join(lines).upper()


def load_designs(args: argparse.Namespace) -> list[dict]:
    """Read the CSV into ``[{name, group, protein}]``, preserving file order."""
    rows = io.read_delimited(args.csv)
    if not rows:
        raise SystemExit(f"No rows in {args.csv}")
    cols = io.find_columns(list(rows[0].keys()), {
        "name": args.id_column, "group": args.group_column,
        "protein": args.sequence_column,
    })
    missing = [c for c in args.carry_columns if c not in rows[0]]
    if missing:
        raise SystemExit(f"--carry-columns not in {args.csv}: {missing}")
    designs = []
    for i, r in enumerate(rows):
        designs.append({
            "order": i + 1,          # CSV order, preserved end to end
            "name": str(r[cols["name"]]).strip(),
            "group": str(r[cols["group"]]).strip(),
            "protein": validate_protein(str(r[cols["protein"]]).strip()),
            "carried": {c: r[c] for c in args.carry_columns},
        })
    names = [d["name"] for d in designs]
    if len(set(names)) != len(names):
        dupes = sorted({n for n in names if names.count(n) > 1})
        raise SystemExit(f"Duplicate IDs in {args.id_column}: {dupes}")
    return designs


def assign_wells(designs: list[dict], args: argparse.Namespace) -> list[str]:
    """Place every design in a well; return the group names in plate order.

    ``--layout group-rows`` gives each design group its own row (spilling onto
    further rows when it is bigger than one). ``--layout sequential`` ignores
    the grouping and fills the plate in CSV order, which is what you want when
    the file is already ordered by something meaningful (a pLDDT ranking, say).
    """
    order: list[str] = []
    for d in designs:
        if d["group"] not in order:
            order.append(d["group"])

    if args.layout == "sequential":
        plate = layout_mod.resolve_plate(args.plate_size)
        wells = plate.wells(order=args.orientation)
        if len(designs) > len(wells):
            raise SystemExit(
                f"{len(designs)} designs exceed the {plate.size}-well plate. "
                "Use --plate-size 384 or split the input.")
        seen: dict[str, int] = {}
        for d, well in zip(designs, wells):
            d["well"] = well
            seen[d["group"]] = seen.get(d["group"], 0) + 1
            d["group_index"] = seen[d["group"]]
        return order

    by_group = {g: [d for d in designs if d["group"] == g] for g in order}
    groups = layout_mod.layout_lines(
        [len(by_group[g]) for g in order], plate=args.plate_size,
        orientation=args.orientation, balance_overflow=not args.fill_rows,
    )
    for g_name, placement in zip(order, groups):
        for i, (d, well) in enumerate(zip(by_group[g_name], placement.wells)):
            d["well"] = well
            d["group_index"] = i + 1
    return order


def main() -> int:
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    src = p.add_argument_group("input")
    src.add_argument("--csv", required=True, type=Path, help="Sequence table.")
    src.add_argument("--id-column", default="sample_id")
    src.add_argument("--group-column", default="design_group")
    src.add_argument("--sequence-column", default="protein_sequence_full")
    src.add_argument("--trim-start-met", action="store_true",
                     help="Drop a leading Met from each design (the vector supplies one).")
    src.add_argument("--carry-columns", nargs="*", default=[], metavar="COL",
                     help="Extra CSV columns to copy through into the mapping file "
                          "(e.g. selection_bucket complex_plddt).")

    lay = p.add_argument_group("layout")
    lay.add_argument("--layout", choices=["group-rows", "sequential"],
                     default="group-rows",
                     help="'group-rows' (default) gives each design group its own "
                          "row; 'sequential' fills the plate in CSV order, ignoring "
                          "the grouping (use when the file is already ranked).")
    lay.add_argument("--plate-size", type=int, choices=[96, 384], default=96)
    lay.add_argument("--orientation", choices=layout_mod.ORIENTATIONS, default="row",
                     help="Fill direction: 'row' = A1..A12,B1.. (default).")
    lay.add_argument("--fill-rows", action="store_true",
                     help="group-rows only: an oversized group fills each row to "
                          "capacity (12+2) instead of spreading evenly (7+7).")
    lay.add_argument("--plate-map", choices=["none", "svg", "png", "both"], default="both")

    dna = p.add_argument_group("fragments")
    dna.add_argument("--destination", required=True,
                     help="Destination plasmid: a path, or raw DNA.")
    dna.add_argument("--enzyme", default="BsaI")
    dna.add_argument("--flank-5", default=None, help="Override the derived 5' adapter.")
    dna.add_argument("--flank-3", default=None, help="Override the derived 3' adapter.")
    dna.add_argument("--species", default="e_coli")
    dna.add_argument("--codon-method", default="use_best_codon",
                     choices=["use_best_codon", "match_codon_usage", "harmonize_rca"],
                     help="DNAChisel codon-optimization objective.")
    dna.add_argument("--avoid", nargs="*", default=[], metavar="ENZYME")
    dna.add_argument("--min-length", type=int, default=300)
    dna.add_argument("--gc-window", type=int, default=None,
                     help="Sliding-window width (bp) for the GC cap. Overrides the "
                          "profile's band of the same width; omit to use the profile.")
    dna.add_argument("--gc-max", type=float, default=0.72)
    dna.add_argument("--gc-min", type=float, default=0.0)
    dna.add_argument("--backend", default=None)
    dna.add_argument("--seed", type=int, default=0)

    syn = p.add_argument_group("synthesis manufacturability")
    syn.add_argument("--synthesis-profile", default="vendor-standard",
                     choices=sorted(synthesis.PROFILES),
                     help="Manufacturability defaults and the acceptance spec the "
                          "run is gated on. Every flag below overrides it. Plain "
                          "codon optimization ('none') leaves ~70%% repeated-8-mer "
                          "coverage, which DNA vendors reject.")
    syn.add_argument("--no-check-synthesis", action="store_true",
                     help="Build under the profile but do not fail on its "
                          "acceptance spec. For inspection, not for ordering.")
    syn.add_argument("--unique-kmer", type=int, default=None, metavar="K",
                     help="Forbid any repeated K-mer (both strands). This is the "
                          "cure for a vendor's repeat-complexity rejection: plain "
                          "codon optimization reuses the same best codon everywhere "
                          "and lands around 70%% repeat coverage. Try 8.")
    syn.add_argument("--kmer-ladder", type=int, nargs="*", default=None, metavar="K",
                     help="K values to try in order when --unique-kmer is infeasible "
                          "for a given protein (default: K, K+1, K+2).")
    syn.add_argument("--soft-kmer-boost", type=float, default=20.0,
                     help="Weight of the soft k-mer-uniqueness objective used on "
                          "ladder steps where the hard constraint is infeasible.")
    syn.add_argument("--min-codon-frequency", type=float, default=None, metavar="F",
                     help="Floor on synonymous-codon usage, so breaking repeats does "
                          "not reintroduce rare codons. 0.10 works well for E. coli.")
    syn.add_argument("--gc-window-2", type=int, default=None, metavar="BP",
                     help="A second, narrower GC window (e.g. 20) with its own band.")
    syn.add_argument("--gc-window-2-max", type=float, default=0.85)
    syn.add_argument("--gc-window-2-min", type=float, default=0.0)
    syn.add_argument("--gc-global-max", type=float, default=None, metavar="F",
                     help="Cap on the GC content of the whole fragment, adapters "
                          "included. Codon optimizers drift GC-rich (E. coli's "
                          "preferred codons are GC-rich), so a vendor's overall-GC "
                          "rule usually needs this set explicitly. Try 0.56.")
    syn.add_argument("--gc-global-min", type=float, default=0.30, metavar="F",
                     help="Floor for whole-fragment GC (only applied with "
                          "--gc-global-max).")
    syn.add_argument("--max-homopolymer-a", type=int, default=None)
    syn.add_argument("--max-homopolymer-g", type=int, default=None)
    syn.add_argument("--check-synthesis", action="store_true",
                     help="Deprecated: the gate is on by default. Kept so existing "
                          "command lines keep working; use --no-check-synthesis to "
                          "turn it off.")
    syn.add_argument("--max-repeat-fraction", type=float, default=0.40)
    syn.add_argument("--max-gc-20", type=float, default=0.90,
                     help="Vendor limit for a 20 bp GC window (the check, not the "
                          "design constraint).")
    syn.add_argument("--check-gc-min", type=float, default=0.25,
                     help="Vendor limit: minimum overall GC (the check).")
    syn.add_argument("--check-gc-max", type=float, default=0.68,
                     help="Vendor limit: maximum overall GC (the check). Set 0.58 "
                          "to gate on the rule that flagged the 260919 order.")
    dna.add_argument("--assembly-checks", type=int, default=-1,
                     help="How many fragments to assemble in silico (-1 = all).")
    dna.add_argument("--expect-prefix", default=None, metavar="AA",
                     help="Residues the vector adds before the design (e.g. MSG). "
                          "With --expect-suffix, the assembled ORF must equal "
                          "prefix+design+suffix exactly - this is what catches a "
                          "frame slip that merely shifts a downstream tag.")
    dna.add_argument("--expect-suffix", default=None, metavar="AA",
                     help="Residues the vector adds after the design (e.g. GSHHHHHH).")

    out = p.add_argument_group("output")
    out.add_argument("--out-dir", type=Path, default=Path("."))
    out.add_argument("--stamp", default="designs")
    args = p.parse_args()

    # ---- load & lay out -------------------------------------------------- #
    designs = load_designs(args)
    if args.trim_start_met:
        for d in designs:
            d["protein"] = d["protein"][1:] if d["protein"].startswith("M") else d["protein"]
    group_order = assign_wells(designs, args)
    print(f"{len(designs)} designs in {len(group_order)} group(s), "
          f"{args.layout} layout ({args.orientation}-major):")
    if args.layout == "sequential":
        print(f"  CSV order preserved: {designs[0]['name']} -> {designs[0]['well']} ... "
              f"{designs[-1]['name']} -> {designs[-1]['well']}")
    # Report whichever axis the groups actually occupy: rows for a row-major
    # fill, columns for a column-major one.
    unit = "row" if args.orientation == "row" else "col"
    key = (lambda w: w[0]) if args.orientation == "row" else (lambda w: int(w[1:]))
    for g in group_order:
        wells = [d["well"] for d in designs if d["group"] == g]
        spans = sorted({key(w) for w in wells}, key=str if unit == "row" else int)
        print(f"  {g:<40} {len(wells):>3} -> {unit}(s) "
              f"{','.join(str(s) for s in spans)} ({wells[0]}-{wells[-1]})")

    # ---- destination ------------------------------------------------------ #
    dest_path = Path(args.destination)
    destination = (read_sequence_file(dest_path) if dest_path.exists()
                   else "".join(args.destination.split()).upper())
    ends = gg.destination_ends(destination, args.enzyme)
    flanks = gg.design_insert_flanks(ends, args.enzyme)
    if args.flank_5 is not None:
        flanks = Flanks(args.flank_5.upper(), flanks.three_prime)
    if args.flank_3 is not None:
        flanks = Flanks(flanks.five_prime, args.flank_3.upper())
    print(f"\n{dest_path.name if dest_path.exists() else 'destination'}: "
          f"{len(destination)} bp, {args.enzyme} overhangs "
          f"{ends.five_overhang}/{ends.three_overhang}, backbone {len(ends.backbone)} bp, "
          f"drop-out {len(ends.dropout)} bp")
    print(f"adapters: 5' {flanks.five_prime}   3' {flanks.three_prime}")

    # ---- optimize --------------------------------------------------------- #
    avoid = [args.enzyme, *[e for e in args.avoid if e != args.enzyme]]

    # The profile supplies the manufacturability defaults; any explicit flag
    # overrides it. Same interface as the other plate/scan scripts, so no run
    # silently falls back to unconstrained codon optimization.
    prof = synthesis.get_profile(args.synthesis_profile)
    unique_kmer = (args.unique_kmer if args.unique_kmer is not None
                   else prof.unique_kmer_size)
    min_freq = (args.min_codon_frequency if args.min_codon_frequency is not None
                else prof.min_codon_frequency)

    # An explicitly given window replaces the profile's band of that width
    # rather than stacking on top of it.
    overridden = {args.gc_window_2, args.gc_window or None}
    extra_gc: list[tuple[int | None, float, float]] = [
        band for band in prof.gc_bands if band[0] not in overridden]
    if args.gc_window_2:
        extra_gc.append((args.gc_window_2, args.gc_window_2_min, args.gc_window_2_max))
    if args.gc_global_max:
        extra_gc.append((None, args.gc_global_min, args.gc_global_max))
    elif prof.gc_global:
        extra_gc.append((None, prof.gc_global[0], prof.gc_global[1]))

    explicit_homo = {b: v for b, v in
                     (("A", args.max_homopolymer_a), ("T", args.max_homopolymer_a),
                      ("G", args.max_homopolymer_g), ("C", args.max_homopolymer_g)) if v}
    homo = explicit_homo or prof.max_homopolymer_by_base or None
    print(f"synthesis profile: {prof.name}, gate "
          + ("off (NOT order-ready)" if args.no_check_synthesis else "on")
          + (f", unique {unique_kmer}-mers, codon freq >= {min_freq}"
             if unique_kmer else ""))

    def constraints_for(kmer: int | None, soft: int | None,
                        min_freq: float | None) -> ConstraintSet:
        return ConstraintSet(
            avoid_enzymes=avoid,
            gc_bounds=(args.gc_min, args.gc_max) if args.gc_window else None,
            gc_window=args.gc_window or None,
            gc_windows=extra_gc,
            max_homopolymer_by_base=homo,
            unique_kmer_size=kmer,
            soft_unique_kmer_size=soft,
            soft_unique_kmer_boost=args.soft_kmer_boost,
            min_codon_frequency=min_freq,
        )

    # Relaxation ladder. K-mer uniqueness is not always satisfiable with
    # synonymous changes alone - a protein carrying a repeated peptide motif
    # forces a repeated k-mer whatever the codons - so fall back to a larger k
    # before giving up, and record which step each design needed.
    ladder = build_ladder(unique_kmer, min_freq, args.kmer_ladder)
    pad_constraints = ConstraintSet(avoid_enzymes=avoid, max_homopolymer=4,
                                    forbid_low_complexity=True, gc_bounds=(0.30, 0.70))
    backend = codon_mod.get_backend(args.backend, seed=args.seed)
    print(f"\noptimizing {len(designs)} sequences with the {backend.name} backend"
          + (f", ladder {ladder}" if len(ladder) > 1 else "") + "...")
    attempts = [(constraints_for(k, s, f), describe_step(k, s, f)) for k, s, f in ladder]
    for i, d in enumerate(designs):
        try:
            frag, step, relaxed = build_fragment_with_ladder(
                d["protein"], attempts, species=args.species,
                flanks=flanks, min_length=args.min_length,
                pad_constraints=pad_constraints, backend=backend,
                seed=args.seed + i, codon_method=args.codon_method,
            )
        except RuntimeError as exc:
            raise SystemExit(f"{d['name']}: {exc}") from exc
        d["coding"], d["final"] = frag.coding, frag.final
        d["ladder_step"], d["relaxed"] = step, relaxed
        if (i + 1) % 10 == 0:
            print(f"  {i + 1}/{len(designs)}", flush=True)

    relaxed = [d for d in designs if d.get("relaxed")]
    if relaxed:
        print(f"  {len(relaxed)} design(s) needed a relaxed step:")
        for d in relaxed:
            print(f"    {d['name']:<44} {d['ladder_step']}")

    # ---- verify ----------------------------------------------------------- #
    problems = verify(designs, args, ends, destination, flanks, avoid)

    # ---- write ------------------------------------------------------------ #
    outputs = write_outputs(designs, group_order, args, flanks)

    print()
    for label, path in outputs.items():
        if path:
            print(f"  wrote {path}  [{label}]")
    if problems:
        print(f"\nFAILED {len(problems)} check(s):")
        for msg in problems[:40]:
            print(f"  - {msg}")
        if len(problems) > 40:
            print(f"  ... and {len(problems) - 40} more")
        return 1
    print("\nAll checks passed.")
    return 0


def verify(designs, args, ends, destination, flanks, avoid) -> list[str]:
    """Translation, forbidden sites, overhangs, GC, wells, and real assembly."""
    from prosy.core import enzymes as enzymes_mod

    problems: list[str] = []
    patterns = enzymes_mod.patterns_for(avoid)
    tol = 1e-9
    seen: set[str] = set()
    spec = synthesis.SynthesisSpec(
        max_repeat_fraction=args.max_repeat_fraction,
        gc_windows=((20, 0.0, args.max_gc_20),),
        gc_bounds=(args.check_gc_min, args.check_gc_max),
    )
    if not args.no_check_synthesis:
        print(f"checking manufacturability: repeats <{args.max_repeat_fraction * 100:.0f}%, "
              f"20 bp GC <{args.max_gc_20 * 100:.0f}%, overall GC "
              f"{args.check_gc_min * 100:.0f}-{args.check_gc_max * 100:.0f}%...")

    n_assemble = len(designs) if args.assembly_checks < 0 else args.assembly_checks
    if n_assemble:
        print(f"assembling {min(n_assemble, len(designs))} construct(s) in silico...")

    for i, d in enumerate(designs):
        tag = d["name"]
        coding, final = d["coding"].upper(), d["final"].upper()
        if translate(coding) != d["protein"]:
            problems.append(f"{tag}: coding DNA does not translate to the design.")
        for pat in patterns:
            if pat in coding:
                problems.append(f"{tag}: forbidden site {pat} in the coding region.")
        if coding not in final:
            problems.append(f"{tag}: coding region not preserved in the final fragment.")
        if not final.startswith(flanks.five_prime.upper()) or \
                not final.endswith(flanks.three_prime.upper()):
            problems.append(f"{tag}: adapters not intact at the fragment ends.")
        if len(final) < args.min_length:
            problems.append(f"{tag}: length {len(final)} < minimum {args.min_length}.")
        problems.extend(f"{tag}: {m}" for m in gg.check_fragment(final, ends, args.enzyme))
        if args.gc_window:
            lo, hi = window_gc_range(final, args.gc_window)
            d["gc_lo"], d["gc_hi"] = lo, hi
            if hi > args.gc_max + tol:
                problems.append(f"{tag}: max {args.gc_window}bp-window GC "
                                f"{hi * 100:.1f}% > cap {args.gc_max * 100:.1f}%.")
            if args.gc_min > 0 and lo < args.gc_min - tol:
                problems.append(f"{tag}: min {args.gc_window}bp-window GC "
                                f"{lo * 100:.1f}% < floor {args.gc_min * 100:.1f}%.")
        else:
            d["gc_lo"], d["gc_hi"] = None, None
        if d["well"] in seen:
            problems.append(f"{tag}: duplicate well {d['well']}.")
        seen.add(d["well"])

        if not args.no_check_synthesis:
            report = synthesis.check(final, spec)
            d["synthesis"] = report
            problems.extend(f"{tag}: {m}" for m in report.problems)

        if i < n_assemble:
            try:
                product = gg.assemble(destination, [final], args.enzyme)
            except gg.CloningError as exc:
                problems.append(f"{tag}: simulated assembly failed - {exc}")
                continue
            d["product_len"] = len(product)
            protein = orf_protein(product)
            if protein is None or d["protein"] not in protein:
                problems.append(f"{tag}: assembled ORF does not contain the design.")
                continue
            d["product_protein"] = protein
            if args.expect_prefix is not None and args.expect_suffix is not None:
                want = args.expect_prefix.upper() + d["protein"] + args.expect_suffix.upper()
                if protein != want:
                    problems.append(
                        f"{tag}: assembled ORF is {protein[:len(args.expect_prefix) + 4]}..."
                        f"{protein[-(len(args.expect_suffix) + 4):]} "
                        f"({len(protein)} aa), expected ...{args.expect_suffix.upper()} "
                        f"({len(want)} aa) - check the reading frame at the 3' junction.")
    return problems


def write_outputs(designs, group_order, args, flanks) -> dict:
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    stamp = args.stamp
    result: dict = {}

    xlsx = out_dir / f"{stamp}_plate.xlsx"
    io.write_plate_xlsx(
        [{"well": d["well"], "name": d["name"], "sequence": d["final"]} for d in designs],
        xlsx)
    result["plate upload sheet"] = xlsx
    result["plate upload sheet (.xls)"] = io.convert_to_xls(xlsx)

    design_path = out_dir / f"{stamp}_design_seqs.csv"
    io.write_csv([{
        "Name": d["name"], "Sequence": d["protein"],
        "5' nucleotides": flanks.five_prime, "3' nucleotides": flanks.three_prime,
    } for d in designs], design_path, io.DESIGN_COLUMNS)
    result["design sheet"] = design_path

    mapping_path = out_dir / f"{stamp}_mapping.csv"
    carried = list(designs[0]["carried"])
    io.write_csv([{
        **d["carried"],
        "Well": d["well"], "Row": d["well"][0], "Order": d["order"],
        "Name": d["name"],
        "Design_group": d["group"], "Group_index": d["group_index"],
        "Protein_len": len(d["protein"]), "Coding_len": len(d["coding"]),
        "Final_len": len(d["final"]),
        "GC": f"{gc_content(d['final']):.3f}",
        "GC_win_min": "" if d["gc_lo"] is None else f"{d['gc_lo']:.3f}",
        "GC_win_max": "" if d["gc_hi"] is None else f"{d['gc_hi']:.3f}",
        "Repeat8_pct": (f"{d['synthesis'].repeat_fraction * 100:.1f}"
                        if d.get("synthesis") else ""),
        "Repeat_dens90_pct": (f"{d['synthesis'].repeat_window_fraction * 100:.1f}"
                              if d.get("synthesis") else ""),
        "GC20_max_pct": (f"{d['synthesis'].gc_window_extremes[20][3] * 100:.0f}"
                         if d.get("synthesis") else ""),
        "Ladder_step": d.get("ladder_step", ""),
        "Product_len": d.get("product_len", ""),
        "Product_protein": d.get("product_protein", ""),
        "Protein": d["protein"], "Coding_DNA": d["coding"], "Final_DNA": d["final"],
    } for d in designs], mapping_path,
        MAPPING_COLUMNS[:6] + carried + MAPPING_COLUMNS[6:])
    result["mapping"] = mapping_path

    fasta = out_dir / f"{stamp}_fragments.fasta"
    with open(fasta, "w", encoding="utf-8") as fh:
        for d in designs:
            fh.write(f">{d['name']} well={d['well']} group={d['group']} "
                     f"len={len(d['final'])}\n")
            for j in range(0, len(d["final"]), 70):
                fh.write(d["final"][j : j + 70] + "\n")
    result["fragments (FASTA)"] = fasta

    if args.plate_map != "none":
        sequential = args.layout == "sequential"
        sizes = {g: sum(1 for x in designs if x["group"] == g) for g in group_order}
        cells = {
            d["well"]: platemap_mod.Cell(
                top=short_label(d["name"], d["group"]),
                bottom=(f"#{d['order']}" if sequential
                        else f"{d['group_index']}/{sizes[d['group']]}"),
                category=d["group"],
                emphasize=(not sequential) and d["group_index"] == 1,
            )
            for d in designs
        }
        shape = ("filled in input order" if sequential
                 else f"one group per {args.orientation}")
        title = (f"{stamp} - {len(designs)} designs, {len(group_order)} groups "
                 f"({shape}, {args.plate_size}-well)")
        if args.plate_map in ("svg", "both"):
            result["plate map (SVG)"] = platemap_mod.write_svg(
                cells, out_dir / f"{stamp}_platemap.svg",
                plate=args.plate_size, title=title)
        if args.plate_map in ("png", "both"):
            result["plate map (PNG)"] = platemap_mod.write_png(
                cells, out_dir / f"{stamp}_platemap.png",
                plate=args.plate_size, title=title)
    return result


def describe_step(kmer: int | None, soft: int | None, min_freq: float | None) -> str:
    bits = []
    if kmer:
        bits.append(f"k={kmer}")
    if soft:
        bits.append(f"soft_k={soft}")
    if min_freq:
        bits.append(f"minfreq={min_freq}")
    return ",".join(bits) or "none"


def build_ladder(k: int | None, freq: float | None,
                 ks: list[int] | None = None
                 ) -> list[tuple[int | None, int | None, float | None]]:
    """``(hard k, soft k, min_codon_frequency)`` steps to try, in order.

    Hard k-mer uniqueness is preferred because it guarantees zero repeats, but
    it is not always satisfiable with synonymous changes alone. Later steps
    raise the hard k (a weaker guarantee) while keeping the *original* k as a
    soft objective, which is what actually drives repeat coverage down; the
    final step drops the hard constraint entirely. Whichever step a design
    lands on, prosy.core.synthesis.check() is still the gate on the result.
    """
    if not k:
        return [(None, None, freq)]
    ks = ks or [k, k + 1, k + 2]
    ladder: list[tuple[int | None, int | None, float | None]] = [(ks[0], None, freq)]
    ladder += [(bigger, k, freq) for bigger in ks[1:]]
    ladder.append((None, k, freq))          # soft-only: never infeasible
    return ladder


def short_label(name: str, group: str) -> str:
    """Trim the group prefix off a sample ID so plate cells stay readable."""
    label = name[len(group):].lstrip("_-") if name.startswith(group) else name
    return label or name


if __name__ == "__main__":
    raise SystemExit(main())
