#!/usr/bin/env python3
"""Mutational scans of a nanobody, CDRs by default, into FP01-ready fragments.

Same three scans as ``scripts/scan/mutational_scan.py`` (alanine, saturation,
reduced alphabet) and the same fragment pipeline, with two nanobody-specific
conveniences:

* **CDR-only by default.** The parent is annotated with
  :mod:`prosy.core.antibody` and only CDR1/2/3 are scanned. Override with
  ``--region`` (``cdr`` / ``cdr3`` / ``all`` / ``framework``) or with an
  explicit ``--positions`` spec, which always wins.
* **Parents by ID.** ``--parent Nb01`` looks the sequence up in the lab
  spreadsheet; ``--protein`` still takes a raw sequence.

The destination defaults to **FP01**, whose two outward-facing BsaI sites drop
out a chloramphenicol stuffer and leave a ``TATG`` / ``GGTT`` pair of overhangs:
the ``ATG`` of the vector's start codon sits in the 5' overhang and the 3'
overhang opens the GS linker, so the insert encodes the nanobody alone and the
product expresses ``M-<nanobody>-GS-linker-HaloTag-FLAG-SNAP``. That fusion is
checked by assembling the plasmid in silico and translating the ORF.

Examples
--------
# Alanine scan of Nb01's CDRs, through to a plate of orderable fragments:
python nanobody_scan.py --parent Nb01 --scan alanine --fragments

# Site-saturation of CDR3 only (19 x |CDR3| fragments):
python nanobody_scan.py --parent Nb01 --scan saturation --region cdr3 --fragments

# Five-residue property scan across all three CDRs:
python nanobody_scan.py --parent Nb01 --scan reduced --alphabet adklw --fragments

# NNK library, one fragment per CDR3 position:
python nanobody_scan.py --parent Nb01 --region cdr3 --degenerate NNK --fragments
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

_HERE = Path(__file__).resolve()
sys.path.insert(0, str(_HERE.parents[2] / "scripts" / "scan"))
sys.path.insert(0, str(_HERE.parents[2]))

import scan_common as sc  # noqa: E402

from prosy.core import antibody  # noqa: E402
from prosy.core import io  # noqa: E402
from prosy.core.sequence import validate_protein  # noqa: E402

DEFAULT_XLSX = _HERE.parents[2] / "Working folder" / "260825" / "Malaria DX plasmids.xlsx"
DEFAULT_SHEET = "Nanobodies"
DEFAULT_DESTINATION = "FP01"
REGIONS = ("cdr", "cdr1", "cdr2", "cdr3", "framework", "all")


def load_parent(args: argparse.Namespace) -> tuple[str, str]:
    """Return ``(name, protein)`` from ``--parent`` (lab sheet) or ``--protein``."""
    if args.protein:
        return args.name or "parent", validate_protein(args.protein)
    if not args.parent:
        raise SystemExit("Provide --parent <NbID> or --protein <sequence>.")
    xlsx = Path(args.xlsx)
    if not xlsx.exists():
        raise SystemExit(f"Nanobody spreadsheet not found: {xlsx}")
    for row in io.read_xlsx_sheet(xlsx, sheet=args.sheet):
        if str(row.get(args.id_column) or "").strip() == args.parent:
            return args.parent, validate_protein(str(row[args.seq_column]).strip())
    raise SystemExit(f"{args.parent!r} not found in {xlsx.name} sheet {args.sheet!r}.")


def resolve_region_positions(protein: str, args: argparse.Namespace) -> str | None:
    """Positions implied by ``--region``; ``--positions`` overrides it."""
    if args.positions:
        return args.positions
    if args.region == "all":
        return None
    try:
        annotation = antibody.annotate(protein, scheme=args.cdr_scheme)
    except antibody.AnnotationError as exc:
        raise SystemExit(
            f"CDR annotation failed ({exc}). Pass --positions explicitly or "
            "--region all to scan the whole domain."
        ) from exc
    print(f"  CDRs ({args.cdr_scheme}): {annotation.describe()}")
    if args.region == "framework":
        positions = annotation.framework_positions()
    elif args.region == "cdr":
        positions = annotation.cdr_positions()
    else:
        positions = annotation.cdr_positions(args.region.upper())
    return ",".join(str(p) for p in positions)


def main() -> int:
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    src = p.add_argument_group("parent")
    src.add_argument("--parent", default=None, metavar="NbID",
                     help="Nanobody ID to look up in the lab spreadsheet.")
    src.add_argument("--protein", default=None, help="Raw protein sequence instead of an ID.")
    src.add_argument("--name", default=None, help="Name used in variant IDs.")
    src.add_argument("--xlsx", default=str(DEFAULT_XLSX), help="Lab spreadsheet path.")
    src.add_argument("--sheet", default=DEFAULT_SHEET)
    src.add_argument("--id-column", default="Nanobody ID")
    src.add_argument("--seq-column", default="AA Sequence")
    src.add_argument("--region", choices=REGIONS, default="cdr",
                     help="Which residues to scan when --positions is not given "
                          "(default: all three CDRs).")
    src.add_argument("--cdr-scheme", choices=antibody.SCHEMES,
                     default=antibody.DEFAULT_SCHEME,
                     help="CDR boundary convention (default imgt, matching the "
                          "CDR columns in the lab spreadsheet).")

    sc.add_scan_args(p, default_positions_help="Overrides --region. Default: the CDRs.")
    sc.add_library_args(p)
    sc.add_output_args(p)
    p.set_defaults(destination=DEFAULT_DESTINATION)
    args = p.parse_args()

    name, protein = load_parent(args)
    print(f"{name}: {len(protein)} residues")
    positions = resolve_region_positions(protein, args)
    return sc.report(sc.run(args, name, protein, positions), args)


if __name__ == "__main__":
    raise SystemExit(main())
