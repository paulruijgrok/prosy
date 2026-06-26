#!/usr/bin/env python3
"""Build a 96-well plate of synthesis-ready nanobody DNA fragments.

Task (260626_NanobodyMuts)
--------------------------
For 12 parent nanobodies (7 from ADL1: Nb01-Nb07; 5 from KRU1: Nb50-53, Nb58)
take the top-7 recommended point mutations each, codon-optimize every variant
for E. coli while avoiding BsaI sites, append Golden Gate flanks, pad to a
minimum length, and lay everything out one parent-per-column on a 96-well plate
(parent in row A, its 7 mutants in rows B-H).

Outputs (written next to the inputs):
  * <stamp>_nanobody_plate.xlsx / .xls  - vendor plate upload (Well/Name/Sequence)
  * design_seqs_nanobodies.csv          - intermediate (Name/Sequence/5'/3')
  * <stamp>_nanobody_mapping.csv        - full provenance table

Run with DNAChisel installed to reproduce codon_optimize.py exactly. Without it
the script falls back to a dependency-free optimizer (still BsaI-clean and
translation-correct) so the pipeline remains runnable/testable.
"""

from __future__ import annotations

import argparse
import sys
from dataclasses import dataclass, field
from pathlib import Path

# Make 'prosy' importable when run directly from a checkout.
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from prosy.core import codon as codon_mod  # noqa: E402
from prosy.core import io  # noqa: E402
from prosy.core.cloning import Flanks, add_flanks, pad_to_length  # noqa: E402
from prosy.core.constraints import ConstraintSet  # noqa: E402
from prosy.core.optimize import optimize_cds  # noqa: E402
from prosy.core.plate import PLATE_96  # noqa: E402
from prosy.core.sequence import (  # noqa: E402
    PointMutation,
    apply_mutation,
    translate,
    validate_protein,
)

# --------------------------------------------------------------------------- #
# Task configuration                                                          #
# --------------------------------------------------------------------------- #

# Parent selection, in plate-column order (col 1..12).
ADL1_PARENTS = ["Nb01", "Nb02", "Nb03", "Nb04", "Nb05", "Nb06", "Nb07"]
KRU1_PARENTS = ["Nb50", "Nb51", "Nb52", "Nb53", "Nb58"]
PARENT_ORDER = ADL1_PARENTS + KRU1_PARENTS

MUTATIONS_PER_PARENT = 7
NEW_ID_START = 73          # mutant IDs start at Nb73
NEW_ID_PREFIX = "Nb"
SPECIES = "e_coli"
AVOID_ENZYMES = ["BsaI"]
MIN_LENGTH = 300

# Golden Gate flanks (BsaI), appended to every coding region.
FLANKS = Flanks(
    five_prime="TGTATCGGTCTCgAGGA",
    three_prime="TATTCCgGAGACCTCTAGT",
)

# TSV files per set.
TSV_BY_SET = {
    "ADL1": "recommendations_ADL1.tsv",
    "KRU1": "recommendations_KRU1.tsv",
}


@dataclass
class Variant:
    new_id: str | None          # Nb73.. for mutants; parent's own ID for parents
    parent_id: str              # e.g. Nb01
    descriptive_id: str         # Nb01 (parent) or Nb01_T27F (mutant)
    is_parent: bool
    set_name: str
    parent_tsv_id: str
    mutation: str | None
    score: int | None
    protein: str
    well: str = ""
    dna_coding: str = ""
    dna_final: str = ""

    @property
    def plate_name(self) -> str:
        # Plate "Name" = new NbID + descriptive parent-style name.
        if self.is_parent:
            return self.parent_id
        return f"{self.new_id} {self.descriptive_id}"


# --------------------------------------------------------------------------- #
# Loading & matching                                                          #
# --------------------------------------------------------------------------- #


def load_parents(xlsx_path: Path) -> dict[str, dict]:
    """Return {NbID: {set, alias, protein}} for the selected parents."""
    rows = io.read_xlsx_sheet(xlsx_path, sheet="Nanobodies")
    by_id: dict[str, dict] = {}
    for r in rows:
        nid = r.get("Nanobody ID")
        if nid in PARENT_ORDER:
            by_id[nid] = {
                "set": (r.get("Set name") or "").strip(),
                "alias": r.get("Alias"),
                "protein": validate_protein(str(r.get("AA Sequence")).strip()),
            }
    missing = [p for p in PARENT_ORDER if p not in by_id]
    if missing:
        raise SystemExit(f"Parents not found in sheet: {missing}")
    return by_id


def load_recommendations(folder: Path) -> dict[str, dict[str, list[dict]]]:
    """Return {set_name: {full_protein_seq: [mutation rows]}}."""
    out: dict[str, dict[str, list[dict]]] = {}
    for set_name, fname in TSV_BY_SET.items():
        rows = io.read_delimited(folder / fname, delimiter="\t")
        by_seq: dict[str, list[dict]] = {}
        for r in rows:
            by_seq.setdefault(str(r["sequence"]).strip(), []).append(r)
        out[set_name] = by_seq
    return out


def top_mutations(rows: list[dict], n: int) -> list[dict]:
    """Top-n mutation rows by score (desc), tie-broken by position (asc)."""
    ranked = sorted(rows, key=lambda m: (-int(m["score"]), int(m["position"])))
    return ranked[:n]


# --------------------------------------------------------------------------- #
# Build variants                                                              #
# --------------------------------------------------------------------------- #


def build_variants(parents: dict[str, dict], recs: dict) -> list[Variant]:
    """Assemble parent + mutant variants in plate-fill order.

    Order: parent by parent (column order); within each parent the 7 mutants in
    score rank. New Nb73+ IDs are assigned across mutants in this same order.
    """
    variants: list[Variant] = []
    next_new_id = NEW_ID_START

    for parent_id in PARENT_ORDER:
        meta = parents[parent_id]
        set_name = meta["set"]
        protein = meta["protein"]
        by_seq = recs[set_name]
        if protein not in by_seq:
            raise SystemExit(
                f"{parent_id}: parent sequence not found in {set_name} TSV (no match)."
            )
        rows = by_seq[protein]
        tsv_id = str(rows[0]["seq_id"])
        muts = top_mutations(rows, MUTATIONS_PER_PARENT)
        if len(muts) < MUTATIONS_PER_PARENT:
            raise SystemExit(
                f"{parent_id}: only {len(muts)} mutations available, need {MUTATIONS_PER_PARENT}."
            )

        # Parent variant (row A of its column).
        variants.append(Variant(
            new_id=parent_id, parent_id=parent_id, descriptive_id=parent_id,
            is_parent=True, set_name=set_name, parent_tsv_id=tsv_id,
            mutation=None, score=None, protein=protein,
        ))

        # Mutant variants (rows B-H).
        for m in muts:
            mut = PointMutation(wt=m["wt"], position=int(m["position"]), mut=m["mut"])
            mutated = apply_mutation(protein, mut, check_wt=True)
            new_id = f"{NEW_ID_PREFIX}{next_new_id}"
            next_new_id += 1
            variants.append(Variant(
                new_id=new_id, parent_id=parent_id,
                descriptive_id=f"{parent_id}_{mut}",
                is_parent=False, set_name=set_name, parent_tsv_id=tsv_id,
                mutation=str(mut), score=int(m["score"]), protein=mutated,
            ))
    return variants


def assign_wells(variants: list[Variant]) -> None:
    """Column-major: each parent block (8 variants) fills one plate column."""
    if len(variants) != PLATE_96.size:
        raise SystemExit(f"Expected {PLATE_96.size} variants, got {len(variants)}.")
    block = PLATE_96.rows  # 8
    for col, start in enumerate(range(0, len(variants), block)):
        for row in range(block):
            variants[start + row].well = PLATE_96.well(row, col)


# --------------------------------------------------------------------------- #
# Optimization                                                                #
# --------------------------------------------------------------------------- #


def optimize_all(variants: list[Variant], backend, *, seed: int) -> None:
    constraints = ConstraintSet(avoid_enzymes=AVOID_ENZYMES)
    pad_constraints = ConstraintSet(
        avoid_enzymes=AVOID_ENZYMES, max_homopolymer=4,
        forbid_low_complexity=True, gc_bounds=(0.30, 0.70),
    )
    for i, v in enumerate(variants):
        result = optimize_cds(
            v.protein, species=SPECIES, constraints=constraints, backend=backend,
        )
        v.dna_coding = result.dna
        with_flanks = add_flanks(v.dna_coding, FLANKS)
        v.dna_final = pad_to_length(
            with_flanks, MIN_LENGTH, constraints=pad_constraints,
            species=SPECIES, seed=seed + i,
        )


# --------------------------------------------------------------------------- #
# Output                                                                       #
# --------------------------------------------------------------------------- #


def write_outputs(variants: list[Variant], out_dir: Path, stamp: str) -> dict:
    out_dir.mkdir(parents=True, exist_ok=True)

    # 1) Plate upload sheet.
    plate_records = [
        {"well": v.well, "name": v.plate_name, "sequence": v.dna_final}
        for v in variants
    ]
    xlsx_path = out_dir / f"{stamp}_nanobody_plate.xlsx"
    io.write_plate_xlsx(plate_records, xlsx_path)
    xls_path = io.convert_to_xls(xlsx_path)

    # 2) Design intermediate (analogous to design_seqs_miniprotein.csv).
    design_rows = [
        {
            "Name": v.descriptive_id,
            "Sequence": v.protein,
            "5' nucleotides": FLANKS.five_prime,
            "3' nucleotides": FLANKS.three_prime,
        }
        for v in variants
    ]
    design_path = out_dir / "design_seqs_nanobodies.csv"
    io.write_csv(design_rows, design_path, io.DESIGN_COLUMNS)

    # 3) Full provenance mapping.
    mapping_rows = [
        {
            "Well": v.well,
            "NbID": v.new_id,
            "NbID_parent": v.descriptive_id,
            "Parent": v.parent_id,
            "Set": v.set_name,
            "Parent_TSV_ID": v.parent_tsv_id,
            "Mutation": v.mutation or "",
            "Score": "" if v.score is None else v.score,
            "Type": "parent" if v.is_parent else "mutant",
            "Protein": v.protein,
            "Protein_len": len(v.protein),
            "Coding_DNA": v.dna_coding,
            "Final_DNA": v.dna_final,
            "Final_len": len(v.dna_final),
        }
        for v in variants
    ]
    mapping_cols = list(mapping_rows[0].keys())
    mapping_path = out_dir / f"{stamp}_nanobody_mapping.csv"
    io.write_csv(mapping_rows, mapping_path, mapping_cols)

    return {
        "plate_xlsx": xlsx_path,
        "plate_xls": xls_path,
        "design": design_path,
        "mapping": mapping_path,
    }


# --------------------------------------------------------------------------- #
# Verification                                                                 #
# --------------------------------------------------------------------------- #


def verify(variants: list[Variant]) -> list[str]:
    """Return a list of problem strings; empty means all checks pass."""
    problems: list[str] = []
    from prosy.core import enzymes

    bsai = enzymes.patterns_for(AVOID_ENZYMES)
    wells_seen = set()
    new_ids = []

    for v in variants:
        # translation round-trip on the coding region
        if translate(v.dna_coding) != v.protein:
            problems.append(f"{v.descriptive_id}: coding DNA does not translate to protein.")
        # coding region BsaI-clean
        cd = v.dna_coding.upper()
        for pat in bsai:
            if pat in cd:
                problems.append(f"{v.descriptive_id}: BsaI site {pat} in coding region.")
        # flanks intact in final
        f5 = FLANKS.five_prime.upper()
        f3 = FLANKS.three_prime.upper()
        if f5 not in v.dna_final.upper() or f3 not in v.dna_final.upper():
            problems.append(f"{v.descriptive_id}: flanks not found intact in final sequence.")
        # coding present in final
        if cd not in v.dna_final.upper():
            problems.append(f"{v.descriptive_id}: coding region not preserved in final sequence.")
        # min length
        if len(v.dna_final) < MIN_LENGTH:
            problems.append(f"{v.descriptive_id}: final length {len(v.dna_final)} < {MIN_LENGTH}.")
        wells_seen.add(v.well)
        if not v.is_parent:
            new_ids.append(v.new_id)

    if len(wells_seen) != PLATE_96.size:
        problems.append(f"Expected {PLATE_96.size} unique wells, got {len(wells_seen)}.")
    expected_ids = [f"{NEW_ID_PREFIX}{NEW_ID_START + i}" for i in range(len(new_ids))]
    if new_ids != expected_ids:
        problems.append("New NbID sequence is not contiguous from Nb73.")
    return problems


# --------------------------------------------------------------------------- #
# Main                                                                         #
# --------------------------------------------------------------------------- #


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--data-dir", type=Path,
        default=Path(__file__).resolve().parents[2] / "Working folder" / "260626_NanobodyMuts",
        help="Folder containing the xlsx + TSV inputs (and where outputs are written).",
    )
    parser.add_argument("--stamp", default="260626", help="Filename prefix for outputs.")
    parser.add_argument(
        "--backend", choices=["auto", "dnachisel", "highest_frequency"], default="auto",
        help="Codon backend. 'auto' uses DNAChisel if installed.",
    )
    parser.add_argument("--seed", type=int, default=1, help="RNG seed (fallback + padding).")
    args = parser.parse_args()

    data_dir = args.data_dir
    xlsx = data_dir / "Malaria DX plasmids.xlsx"

    prefer = None if args.backend == "auto" else args.backend
    backend = codon_mod.get_backend(prefer, seed=args.seed)
    print(f"Codon backend: {backend.name}")

    parents = load_parents(xlsx)
    recs = load_recommendations(data_dir)
    variants = build_variants(parents, recs)
    assign_wells(variants)
    optimize_all(variants, backend, seed=args.seed)

    problems = verify(variants)
    out = write_outputs(variants, data_dir, args.stamp)

    print(f"\nParents: {len(PARENT_ORDER)}  Variants: {len(variants)} "
          f"(mutants: {sum(not v.is_parent for v in variants)})")
    print(f"Plate xlsx : {out['plate_xlsx']}")
    print(f"Plate xls  : {out['plate_xls'] or '(LibreOffice not available - xlsx only)'}")
    print(f"Design CSV : {out['design']}")
    print(f"Mapping CSV: {out['mapping']}")

    if problems:
        print("\nVERIFICATION FAILED:")
        for p in problems:
            print("  -", p)
        sys.exit(1)
    print("\nAll verification checks passed.")


if __name__ == "__main__":
    main()
