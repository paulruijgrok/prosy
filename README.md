# ProSy

Tools to turn protein sequences into ready-to-order DNA sequences.

ProSy is split in two:

- **`prosy/core/`** — a general, organism- and cloning-agnostic toolkit. Enzymes,
  codon tables, GC targets, flanks and assembly choices are all parameters, so
  the same core serves any host or downstream cloning strategy.
- **`scripts/`** — task-specific scripts that wire the core together for concrete
  jobs encountered in daily practice (starting with `scripts/nanobodies/`).

## Core modules

| Module | Purpose |
|---|---|
| `sequence` | Genetic code, translation, reverse-complement, point-mutation parsing/applying, validation. |
| `enzymes` | Restriction-enzyme registry (BsaI, BsmBI, BbsI, SapI, EcoRI, …) with both-strand patterns and Type IIS cut metadata. |
| `constraints` | Declarative `ConstraintSet` → flat list of forbidden patterns + GC bounds for a backend. |
| `codon` | Codon-usage tables and two interchangeable optimization backends. |
| `optimize` | `optimize_cds()` — high-level protein → synthesis-ready CDS. |
| `cloning` | Flank/adapter application and minimum-length padding (input preserved verbatim). |
| `optimize` | Also `build_fragment()` — the optimize → flank → pad composition. |
| `plate` | 96-/384-well coordinate helpers, column- or row-major fill. |
| `layout` | Place groups (parent + N mutants) onto a plate; row/column packing, 96/384, alignment checks. |
| `platemap` | Dependency-free SVG plate-map renderer (colour-coded, labeled); optional PNG via matplotlib/LibreOffice. |
| `io` | Read/write CSV/TSV/XLSX, flexible column matching, vendor plate-sheet writer. |

## Codon backends

`optimize_cds()` selects a backend automatically:

- **DNAChisel** (default when installed) reproduces the original
  `codon_optimize.py`: reverse-translate → enforce translation → avoid
  restriction sites → codon-optimize for a species via a global solver.
- **Highest-frequency fallback** (no dependencies) picks the most frequent
  synonymous codon per residue and greedily swaps codons to clear forbidden
  patterns. It keeps the whole pipeline runnable and testable without
  DNAChisel, but is not a substitute for DNAChisel on the production run.

Install the production backend with:

```bash
pip install -e ".[optimize]"
```

## Optional extras

The core needs only `openpyxl`. Optional feature groups:

- `.[optimize]` — DNAChisel (production codon backend).
- `.[viz]` — matplotlib (PNG plate maps; SVG needs nothing).
- `.[bio]` — Biopython, python-codon-tables, Plateo (wider synbio ecosystem:
  sequence I/O, codon tables for many organisms, plate graphics / robot picklists).
- `.[all]` — everything above plus pytest.

## Quick start

```python
from prosy.core.constraints import ConstraintSet
from prosy.core.optimize import optimize_cds

res = optimize_cds("QVQLVESGGG...", species="e_coli",
                   constraints=ConstraintSet(avoid_enzymes=["BsaI"]))
print(res.dna, res.backend)
```

## Tests

```bash
python -m pytest tests/        # or: python tests/test_core.py
```
