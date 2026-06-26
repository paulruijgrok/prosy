# scripts/nanobodies

Task scripts for nanobody DNA-fragment design.

## make_nanobody_plate.py

Builds a 96-well plate of synthesis-ready nanobody fragments: 12 parent
nanobodies, each with its top-7 recommended point mutations, one parent per
plate column (parent in row A, mutants in rows B–H).

**Inputs** (in `--data-dir`, default `Working folder/260626_NanobodyMuts/`):

- `Malaria DX plasmids.xlsx` — sheet `Nanobodies` (parent IDs ↔ AA sequences).
- `recommendations_ADL1.tsv`, `recommendations_KRU1.tsv` — proposed mutations.

Parents are matched to the TSVs by full protein sequence. Mutations are ranked
by `score` (desc, ties broken by lowest position) and the top 7 kept.

**Selection:** ADL1 → Nb01–Nb07; KRU1 → Nb50, Nb51, Nb52, Nb53, Nb58.

**IDs:** each mutant gets a new sequential `NbID` starting at `Nb73` (assigned in
plate-fill order), plus a descriptive `NbID_parent` like `Nb01_T27F`. Parents
keep their own ID. Plate `Name` = `"<NbID> <NbID_parent>"` for mutants, the bare
ID for parents.

**Processing per variant:** codon-optimize for *E. coli* avoiding BsaI →
append Golden Gate flanks (`TGTATCGGTCTCgAGGA` / `TATTCCgGAGACCTCTAGT`) → pad to
a 300 bp minimum (BsaI-/repeat-avoiding) if shorter.

**Outputs** (written to `--data-dir`):

- `<stamp>_nanobody_plate.xlsx` (+ `.xls` if LibreOffice is present) — vendor
  upload sheet: `Well Position`, `Name`, `Sequence`.
- `design_seqs_nanobodies.csv` — intermediate, same schema as
  `design_seqs_miniprotein.csv`.
- `<stamp>_nanobody_mapping.csv` — full provenance (well, both IDs, parent, set,
  TSV id, mutation, score, protein, coding + final DNA, lengths).

### Run

```bash
# Production (reproduces codon_optimize.py): requires dnachisel
python scripts/nanobodies/make_nanobody_plate.py

# Force the dependency-free fallback (no dnachisel needed)
python scripts/nanobodies/make_nanobody_plate.py --backend highest_frequency
```

The script self-verifies (96 unique wells, translation round-trip, BsaI-clean
coding regions, intact flanks, contiguous Nb73+ IDs, min length) and exits
non-zero if any check fails.
