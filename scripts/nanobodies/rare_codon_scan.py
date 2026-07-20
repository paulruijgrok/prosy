#!/usr/bin/env python3
"""Scan coding sequences for E. coli rare codons and tandem rare-codon runs.

Rare codons are decoded by low-abundance tRNAs in E. coli; isolated ones are
usually tolerated, but **tandem runs** (two or more in a row) — or rare codons
clustered near the 5' end — are associated with ribosome stalling, frameshifting
and reduced yield. This tool flags them in a nanobody mapping CSV's `Coding_DNA`.

The rare set is the classic tRNA-limited group that expression strains (e.g.
Rosetta) supplement: AGA, AGG, CGA, CGG (Arg), ATA (Ile), CTA (Leu), CCC (Pro),
GGA (Gly). GGG (Gly) is borderline and reported separately as "watch".

Example:
    python rare_codon_scan.py \
        --input "../../Working folder/260720_NanobodyMuts_gc75/260626_nanobody_mapping.csv"
"""
from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

# Classic tRNA-limited rare codons in E. coli (amino acid in comment).
RARE = {
    "AGA": "R", "AGG": "R", "CGA": "R", "CGG": "R",
    "ATA": "I", "CTA": "L", "CCC": "P", "GGA": "G",
}
WATCH = {"GGG": "G"}  # borderline-low usage, reported but not counted as rare


def codons_of(dna: str) -> list[str]:
    dna = dna.upper()
    return [dna[i:i + 3] for i in range(0, len(dna) - len(dna) % 3, 3)]


def rare_runs(codons: list[str], min_len: int = 2) -> list[tuple[int, list[str]]]:
    """Return (start_residue_1based, [codons]) for each maximal run of >=min_len
    consecutive rare codons."""
    runs = []
    i = 0
    n = len(codons)
    while i < n:
        if codons[i] in RARE:
            j = i
            while j < n and codons[j] in RARE:
                j += 1
            if j - i >= min_len:
                runs.append((i + 1, codons[i:j]))
            i = j
        else:
            i += 1
    return runs


def scan(rows: list[dict], seq_col: str) -> dict:
    per_variant = []
    total_codons = total_rare = total_watch = 0
    rare_by_codon: dict[str, int] = {}
    variants_with_runs = 0
    total_runs = 0
    longest = 0

    for r in rows:
        seq = (r.get(seq_col) or "").upper()
        if not seq:
            continue
        codons = codons_of(seq)
        total_codons += len(codons)
        rare_positions = [(k + 1, c) for k, c in enumerate(codons) if c in RARE]
        for _, c in rare_positions:
            rare_by_codon[c] = rare_by_codon.get(c, 0) + 1
        total_rare += len(rare_positions)
        total_watch += sum(1 for c in codons if c in WATCH)
        runs = rare_runs(codons)
        if runs:
            variants_with_runs += 1
            total_runs += len(runs)
            longest = max(longest, max(len(c) for _, c in runs))
        per_variant.append({
            "id": r.get("NbID_parent") or r.get("NbID") or "",
            "well": r.get("Well", ""),
            "n_rare": len(rare_positions),
            "runs": runs,
        })

    return {
        "n_variants": len(per_variant),
        "total_codons": total_codons,
        "total_rare": total_rare,
        "total_watch": total_watch,
        "rare_by_codon": rare_by_codon,
        "variants_with_runs": variants_with_runs,
        "total_runs": total_runs,
        "longest_run": longest,
        "per_variant": per_variant,
    }


def print_report(name: str, s: dict, show_runs: bool) -> None:
    pct = 100 * s["total_rare"] / s["total_codons"] if s["total_codons"] else 0
    print(f"=== {name} ===")
    print(f"variants: {s['n_variants']}   codons: {s['total_codons']}")
    print(f"rare codons: {s['total_rare']} ({pct:.2f}% of codons)   "
          f"GGG (watch): {s['total_watch']}")
    if s["rare_by_codon"]:
        by = ", ".join(f"{c}={n}({RARE[c]})" for c, n in
                       sorted(s["rare_by_codon"].items(), key=lambda kv: -kv[1]))
        print(f"by codon: {by}")
    print(f"tandem runs (>=2 in a row): {s['total_runs']} "
          f"across {s['variants_with_runs']} variant(s); longest run = "
          f"{s['longest_run']} codons")
    if show_runs and s["total_runs"]:
        for v in s["per_variant"]:
            for start, cod in v["runs"]:
                print(f"    {v['id']:14} {v['well']:4} residues {start}-"
                      f"{start + len(cod) - 1}: {'-'.join(cod)}")
    print()


def main(argv=None) -> int:
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--input", required=True, nargs="+",
                   help="One or more mapping CSVs (labelled by parent dir name).")
    p.add_argument("--seq-column", default="Coding_DNA")
    p.add_argument("--show-runs", action="store_true",
                   help="List every tandem run with its residue positions.")
    args = p.parse_args(argv)

    for path in args.input:
        rows = list(csv.DictReader(open(path, newline="")))
        label = Path(path).parent.name or Path(path).name
        print_report(label, scan(rows, args.seq_column), args.show_runs)
    return 0


if __name__ == "__main__":
    sys.exit(main())
