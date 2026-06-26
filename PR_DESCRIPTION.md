# ProSy: core package + nanobody plate script

Initial commit of **ProSy**, a toolkit for turning protein sequences into
synthesis-ready DNA, plus the first task script.

## What's included

### `prosy/core/` — general, reusable toolkit
Organism- and cloning-agnostic building blocks. Enzymes, codon tables, GC
targets, flanks and assembly choices are all parameters.

| Module | Purpose |
|---|---|
| `sequence` | Genetic code, translation, reverse-complement, point-mutation parse/apply, validation. |
| `enzymes` | Restriction-enzyme registry (BsaI, BsmBI, BbsI, SapI, EcoRI, …) with both-strand patterns + Type IIS cut metadata. |
| `constraints` | Declarative `ConstraintSet` → forbidden patterns + GC bounds. |
| `codon` | Codon-usage tables and two interchangeable optimization backends. |
| `optimize` | `optimize_cds()` — protein → synthesis-ready CDS. |
| `cloning` | Flank/adapter application + minimum-length padding (input preserved verbatim). |
| `plate` | 96-/384-well coordinate helpers, column- or row-major fill. |
| `io` | Read/write CSV/TSV/XLSX, flexible column matching, vendor plate-sheet writer. |

### Codon backends
- **DNAChisel** (default when installed) reproduces the original
  `codon_optimize.py`: reverse-translate → enforce translation → avoid
  restriction sites → codon-optimize for a species.
- **Highest-frequency fallback** (no dependencies) keeps the pipeline runnable
  and testable without DNAChisel.

### `scripts/nanobodies/make_nanobody_plate.py` — first task script
Builds a 96-well plate of nanobody fragments: 12 parents (ADL1 Nb01–Nb07; KRU1
Nb50–Nb53, Nb58), each with its top-7 recommended mutations, one parent per
column. Codon-optimizes for *E. coli* avoiding BsaI, appends Golden Gate flanks,
pads to a 300 bp minimum, and emits a vendor plate sheet, a design intermediate,
and a full provenance mapping. Self-verifies on every run.

## Tests
`tests/test_core.py` — 9 unit tests covering the genetic code, mutations, the
enzyme registry, the fallback optimizer (translation round-trip + BsaI-clean),
padding, and plate layout.

```bash
python -m pytest tests/        # or: python tests/test_core.py
```

## Notes
- Project data (`Working folder/`) and the original reference script
  (`Inspiration code/`) are intentionally git-ignored.
- Requires Python ≥ 3.10. Runtime dep: `openpyxl`; optional: `dnachisel`
  (`pip install -e ".[optimize]"`).
