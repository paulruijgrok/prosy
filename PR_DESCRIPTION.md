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
| `optimize` | `optimize_cds()` and `build_fragment()` (optimize → flank → pad). |
| `plate` | 96-/384-well coordinate helpers, column- or row-major fill. |
| `layout` | Place groups (parent + N mutants) onto a plate; row/column packing, 96/384, alignment checks. |
| `platemap` | Dependency-free SVG plate-map renderer; optional PNG via matplotlib/LibreOffice. |
| `io` | Read/write CSV/TSV/XLSX, flexible column matching, vendor plate-sheet writer. |

### Codon backends
- **DNAChisel** (default when installed) reproduces the original
  `codon_optimize.py`: reverse-translate → enforce translation → avoid
  restriction sites → codon-optimize for a species.
- **Highest-frequency fallback** (no dependencies) keeps the pipeline runnable
  and testable without DNAChisel.

### `scripts/nanobodies/` — task scripts
- `nanobody_common.py`: shared library (load parents, match mutations by
  sequence, build variants, lay out, optimize, write, verify).
- `make_nanobody_plate.py`: one general CLI. Each parent + its top-N mutations is
  a group that fills a plate column or row; groups pack onto a 96- or 384-well
  plate. Codon-optimizes for *E. coli* avoiding BsaI, appends Golden Gate flanks,
  pads to a minimum length, and emits a vendor plate sheet, a design intermediate
  and a full provenance mapping (all `--stamp`-prefixed). Self-verifies on every
  run. Supports `--orientation row|column`, `--mutations-per-parent`,
  configurable parent lists and `--plate-size 96|384`.

## Tests
- `tests/test_core.py` — genetic code, mutations, enzyme registry, fallback
  optimizer (translation round-trip + BsaI-clean), padding, plate coordinates.
- `tests/test_layout.py` — layout engine: column/row, 96/384, the 16×5 and 24×3
  packings, overflow and straddle handling.

```bash
python -m pytest tests/        # or: python tests/test_core.py
```

## Optional extras
Core needs only `openpyxl`. Feature groups: `.[optimize]` (DNAChisel),
`.[viz]` (matplotlib PNG maps), `.[bio]` (Biopython, python-codon-tables,
Plateo), `.[all]`.

## Notes
- Project data (`Working folder/`) and the original reference script
  (`Inspiration code/`) are intentionally git-ignored.
- Requires Python ≥ 3.10. Runtime dep: `openpyxl`.
