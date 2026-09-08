"""Shared plumbing for the mutational-scan command-line tools.

Holds the argument groups, input loading, run orchestration and reporting that
``mutational_scan.py`` and ``scripts/nanobodies/nanobody_scan.py`` have in
common, so the nanobody entry point differs from the generic one only in how it
picks its default positions and where it reads the parent sequence from.
"""

from __future__ import annotations

import argparse
import sys
from dataclasses import dataclass
from pathlib import Path

_REPO_ROOT = Path(__file__).resolve().parents[2]
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))

from prosy.core import codon as codon_mod  # noqa: E402
from prosy.core import goldengate as gg  # noqa: E402
from prosy.core import io  # noqa: E402
from prosy.core import library as lib  # noqa: E402
from prosy.core import scan as scan_mod  # noqa: E402
from prosy.core.cloning import Flanks  # noqa: E402
from prosy.core.sequence import validate_protein  # noqa: E402

DEFAULT_PLASMID_DIR = _REPO_ROOT / "data" / "plasmids"


# --------------------------------------------------------------------------- #
# Input loading                                                                #
# --------------------------------------------------------------------------- #


def read_sequence_file(path: Path) -> str:
    """Read a plain-text or FASTA sequence file into one uppercase string."""
    text = Path(path).read_text(encoding="utf-8")
    lines = [ln.strip() for ln in text.splitlines() if ln.strip() and not ln.startswith(">")]
    return "".join(lines).upper()


def load_protein(args: argparse.Namespace) -> tuple[str, str]:
    """Return ``(name, protein)`` from ``--protein`` or ``--protein-file``."""
    if args.protein:
        protein = validate_protein(args.protein)
        return args.name or "parent", protein
    if args.protein_file:
        protein = validate_protein(read_sequence_file(Path(args.protein_file)))
        return args.name or Path(args.protein_file).stem, protein
    raise SystemExit("Provide --protein or --protein-file.")


def load_plasmid(spec: str) -> str:
    """Resolve a destination plasmid given a path, a bare name, or raw DNA."""
    candidate = Path(spec)
    if candidate.exists():
        return read_sequence_file(candidate)
    for suffix in (".fa", ".fasta", ".txt", ".seq"):
        alt = DEFAULT_PLASMID_DIR / f"{spec}{suffix}"
        if alt.exists():
            return read_sequence_file(alt)
    stripped = "".join(spec.split()).upper()
    if stripped and set(stripped) <= set("ACGT"):
        return stripped
    raise SystemExit(
        f"Could not resolve destination {spec!r}: not a file, not a plasmid in "
        f"{DEFAULT_PLASMID_DIR}, and not raw DNA."
    )


# --------------------------------------------------------------------------- #
# Argument groups                                                              #
# --------------------------------------------------------------------------- #


def add_scan_args(p: argparse.ArgumentParser, *, default_positions_help: str) -> None:
    g = p.add_argument_group("scan")
    g.add_argument("--scan", choices=scan_mod.SCAN_TYPES, default="alanine",
                   help="alanine = one Ala variant per position; saturation = all 19 "
                        "substitutions; reduced = a 5-residue property alphabet.")
    g.add_argument("--positions", default=None, metavar="SPEC",
                   help=f"Positions to scan, e.g. '31-35,50-65,96-110'. {default_positions_help}")
    g.add_argument("--exclude", default=None, metavar="SPEC",
                   help="Positions to drop from the scan (same syntax as --positions).")
    g.add_argument("--alphabet", default=None,
                   help="For --scan reduced: a preset name "
                        f"({', '.join(sorted(scan_mod.REDUCED_ALPHABETS))}) or an explicit "
                        "residue string like ADKLW. For --scan saturation: the residues "
                        "to substitute in (default all 20).")
    g.add_argument("--alanine-to", default="A", metavar="AA",
                   help="Residue an alanine scan substitutes in (default A).")
    g.add_argument("--native-alanine-to", default="G", metavar="AA",
                   help="Substitute used where the wild type already is the alanine-scan "
                        "residue; 'none' to skip those positions instead.")
    g.add_argument("--no-parent", action="store_true",
                   help="Omit the unmutated parent from the output set.")
    g.add_argument("--degenerate", default=None, metavar="CODON",
                   help="Build a degenerate-codon library (e.g. NNK) instead of explicit "
                        "variants: one fragment per position. Requires --fragments.")


def add_library_args(p: argparse.ArgumentParser) -> None:
    g = p.add_argument_group("fragments (DNA)")
    g.add_argument("--fragments", action="store_true",
                   help="Also build synthesis-ready DNA fragments for the scan.")
    g.add_argument("--destination", default=None, metavar="PLASMID",
                   help="Destination plasmid (path, name under data/plasmids/, or raw DNA). "
                        "Adapters and overhangs are derived from it.")
    g.add_argument("--enzyme", default="BsaI", help="Type IIS cloning enzyme (default BsaI).")
    g.add_argument("--flank-5", default=None, help="Override the derived 5' adapter.")
    g.add_argument("--flank-3", default=None, help="Override the derived 3' adapter.")
    g.add_argument("--species", default="e_coli", help="Codon-usage species.")
    g.add_argument("--avoid", nargs="*", default=[], metavar="ENZYME",
                   help="Extra enzymes whose sites the coding region must avoid.")
    g.add_argument("--min-length", type=int, default=300,
                   help="Vendor minimum fragment length; shorter fragments are padded.")
    g.add_argument("--gc-window", type=int, default=None,
                   help="Sliding-window width (bp) for the GC cap; omit to disable.")
    g.add_argument("--gc-max", type=float, default=0.72)
    g.add_argument("--gc-min", type=float, default=0.0)
    g.add_argument("--plate-size", type=int, choices=[96, 384], default=96)
    g.add_argument("--plate-order", choices=["column", "row"], default="column")
    g.add_argument("--backend", default=None,
                   help="Codon backend: 'dnachisel', 'highest_frequency', or omit to auto-select.")
    g.add_argument("--seed", type=int, default=0)
    g.add_argument("--assembly-checks", type=int, default=3,
                   help="How many fragments to assemble in silico and translate "
                        "(-1 = all, 0 = none).")


def add_output_args(p: argparse.ArgumentParser) -> None:
    g = p.add_argument_group("output")
    g.add_argument("--out-dir", type=Path, default=Path("."),
                   help="Directory for the output files.")
    g.add_argument("--stamp", default=None,
                   help="Output filename prefix (default: <name>_<scan>).")


# --------------------------------------------------------------------------- #
# Running                                                                      #
# --------------------------------------------------------------------------- #


@dataclass
class ScanRun:
    name: str
    protein: str
    variants: list
    members: list
    problems: list[str]
    outputs: dict
    ends: gg.DestinationEnds | None


def build_variants(args: argparse.Namespace, name: str, protein: str,
                   positions) -> list:
    """Apply the scan requested on the command line."""
    kwargs = {
        "positions": positions, "exclude": args.exclude, "parent_name": name,
        "include_parent": not args.no_parent,
    }
    if args.scan == "alanine":
        native = None if str(args.native_alanine_to).lower() == "none" else args.native_alanine_to
        return scan_mod.alanine_scan(protein, substitute=args.alanine_to,
                                     native_substitute=native, **kwargs)
    if args.scan == "saturation":
        if args.alphabet:
            kwargs["alphabet"] = args.alphabet
        return scan_mod.saturation_scan(protein, **kwargs)
    if args.alphabet:
        kwargs["alphabet"] = args.alphabet
    return scan_mod.reduced_alphabet_scan(protein, **kwargs)


def write_variants_csv(variants: list, path: Path) -> Path:
    rows = [{
        "Name": v.name, "Parent": v.parent, "Scan": v.scan,
        "Mutation": v.mutation_label,
        "Position": v.position if v.position is not None else "",
        "WT": v.mutations[0].wt if v.mutations else "",
        "Mut": v.mutations[0].mut if v.mutations else "",
        "Category": v.category or "", "Length": len(v.protein), "Protein": v.protein,
    } for v in variants]
    io.write_csv(rows, path, ["Name", "Parent", "Scan", "Mutation", "Position",
                              "WT", "Mut", "Category", "Length", "Protein"])
    return path


def run(args: argparse.Namespace, name: str, protein: str, positions) -> ScanRun:
    """Scan -> (optionally) fragments -> verify -> write. Shared by both CLIs."""
    stamp = args.stamp or f"{name}_{args.scan}"
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    if args.degenerate:
        if not args.fragments:
            raise SystemExit("--degenerate builds DNA fragments; add --fragments.")
        # A degenerate library has no enumerable protein set: its members are
        # only defined after sequencing. Carry the parent through instead.
        variants = [scan_mod.ScanVariant(name=name, parent=name, protein=protein)]
        outputs: dict = {}
    else:
        variants = build_variants(args, name, protein, positions)
        outputs = {"variants": write_variants_csv(variants,
                                                  out_dir / f"{stamp}_variants.csv")}

    if not args.fragments:
        return ScanRun(name, protein, variants, [], [], outputs, None)

    destination = load_plasmid(args.destination) if args.destination else None
    ends = None
    if destination:
        ends = gg.destination_ends(destination, args.enzyme)
        flanks = gg.design_insert_flanks(ends, args.enzyme)
    elif args.flank_5 is not None and args.flank_3 is not None:
        flanks = Flanks(args.flank_5, args.flank_3)
    else:
        raise SystemExit("--fragments needs either --destination or both --flank-5/--flank-3.")
    if args.flank_5 is not None:
        flanks = Flanks(args.flank_5, flanks.three_prime)
    if args.flank_3 is not None:
        flanks = Flanks(flanks.five_prime, args.flank_3)

    cfg = lib.LibraryConfig(
        flanks=flanks, enzyme=args.enzyme, species=args.species,
        avoid_enzymes=list(args.avoid), min_length=args.min_length,
        gc_window=args.gc_window, gc_max=args.gc_max, gc_min=args.gc_min,
        plate_size=args.plate_size, plate_order=args.plate_order, seed=args.seed,
    )
    backend = codon_mod.get_backend(args.backend, seed=args.seed)

    if args.degenerate:
        targets = [p for p, *_ in scan_mod.degenerate_saturation_positions(
            protein, positions, exclude=args.exclude, codon=args.degenerate)]
        members = lib.build_degenerate_library(variants[0], targets, cfg,
                                               codon=args.degenerate, backend=backend)
    else:
        members = lib.build_library(variants, cfg, backend=backend, progress=True)

    checks = None if args.assembly_checks < 0 else args.assembly_checks
    problems = lib.verify_library(
        members, cfg, ends=ends, destination=destination if checks else None,
        expected_protein=(lambda m: m.variant.protein) if checks else None,
        assembly_checks=checks,
    )
    outputs |= lib.write_library(members, out_dir, stamp, cfg)
    return ScanRun(name, protein, variants, members, problems, outputs, ends)


def report(run_result: ScanRun, args: argparse.Namespace) -> int:
    """Print a human-readable summary; return the process exit code."""
    r = run_result
    if args.degenerate:
        residues, has_stop = scan_mod.degenerate_residues(args.degenerate)
        print(f"\n{r.name}: {len(r.protein)} residues, {args.degenerate} library at "
              f"{len(r.members)} position(s); {args.degenerate} encodes "
              f"{len(residues)} residue(s)"
              f"{' plus a stop codon' if has_stop else ' and no stop codon'}")
    else:
        print(f"\n{r.name}: {len(r.protein)} residues, {len(r.variants)} variant(s) "
              f"({args.scan} scan)")
    if r.ends is not None:
        print(f"  destination {args.enzyme}: 5' overhang {r.ends.five_overhang}, "
              f"3' overhang {r.ends.three_overhang}, backbone {len(r.ends.backbone)} bp, "
              f"drop-out {len(r.ends.dropout)} bp")
    if r.members:
        lengths = [len(m.dna_final) for m in r.members]
        plates = max(m.plate for m in r.members)
        print(f"  {len(r.members)} fragment(s), {min(lengths)}-{max(lengths)} bp, "
              f"{plates} plate(s)")
    for value in r.outputs.values():
        if isinstance(value, list):
            for v in value:
                print(f"  wrote {v}")
        elif value:
            print(f"  wrote {value}")
    if r.problems:
        print(f"\nFAILED {len(r.problems)} check(s):")
        for p in r.problems[:40]:
            print(f"  - {p}")
        if len(r.problems) > 40:
            print(f"  ... and {len(r.problems) - 40} more")
        return 1
    if r.members:
        print("\nAll checks passed.")
    return 0
