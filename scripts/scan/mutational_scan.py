#!/usr/bin/env python3
"""Mutational scans of any protein, optionally through to ordered DNA.

Three scan types, all over an explicitly chosen set of positions (or the whole
sequence when none is given):

* ``--scan alanine``    - one variant per position, wild type -> Ala.
* ``--scan saturation`` - site-saturation: all 19 non-wild-type residues.
* ``--scan reduced``    - saturation restricted to a 5-residue alphabet that
  spans charge, hydrophobicity and size (default A/D/K/L/W).

With ``--fragments`` the resulting protein set is codon-optimized, given Golden
Gate adapters derived from a destination plasmid, padded to the vendor minimum,
verified (translation, forbidden sites, overhangs, windowed GC, in-silico
assembly) and written out as plate-upload sheets.

Examples
--------
# Alanine scan of a whole protein, protein sequences only:
python mutational_scan.py --protein-file myprotein.fa --scan alanine

# Site-saturation of three loops, straight through to orderable DNA:
python mutational_scan.py --protein-file myprotein.fa --scan saturation \\
    --positions 31-35,50-65,96-110 --fragments --destination FP01

# Five-residue property scan, 384-well plates:
python mutational_scan.py --protein-file myprotein.fa --scan reduced \\
    --alphabet adklw --fragments --destination FP01 --plate-size 384
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import scan_common as sc  # noqa: E402


def main() -> int:
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    src = p.add_argument_group("input")
    src.add_argument("--protein", default=None, help="Protein sequence, single-letter.")
    src.add_argument("--protein-file", default=None,
                     help="FASTA or plain-text file holding one protein sequence.")
    src.add_argument("--name", default=None, help="Parent name used in variant IDs.")

    sc.add_scan_args(p, default_positions_help="Default: every residue.")
    sc.add_library_args(p)
    sc.add_output_args(p)
    args = p.parse_args()

    name, protein = sc.load_protein(args)
    return sc.report(sc.run(args, name, protein, args.positions), args)


if __name__ == "__main__":
    raise SystemExit(main())
