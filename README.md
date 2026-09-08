# ProSy

Tools to turn protein sequences into ready-to-order DNA sequences: codon
optimization under synthesis constraints, mutational-scan library design, Golden
Gate adapter design against a real destination plasmid, and plate-ready order
files — with every fragment verified in silico before it can be ordered.

ProSy is split in two:

- **`prosy/core/`** — a general, organism- and cloning-agnostic toolkit. Enzymes,
  codon tables, GC targets, flanks and assembly choices are all parameters, so
  the same core serves any host or downstream cloning strategy.
- **`scripts/`** — task-specific scripts that wire the core together for concrete
  jobs encountered in daily practice (`scripts/scan/`, `scripts/nanobodies/`).

## Quick start

```bash
pip install -e ".[optimize,dev]"
python -m pytest tests/ -q

# Alanine-scan a nanobody's CDRs and build fragments for Golden Gate into FP01:
python scripts/nanobodies/nanobody_scan.py --parent Nb01 --scan alanine \
    --fragments --gc-window 50 --gc-max 0.72 --out-dir /tmp/nb01
```

That writes a plate-upload `.xlsx`, a design CSV and a full mapping CSV, and
prints `All checks passed.` only if every fragment translated correctly, carried
the right overhangs and assembled into the expected fusion protein in silico.

The same three scans on any protein, with the whole sequence as the default
target set:

```bash
python scripts/scan/mutational_scan.py --protein-file myprotein.fa --scan saturation \
    --positions 31-35,50-65 --fragments --destination FP01
```

## Installation

Python ≥ 3.10. The core needs only `openpyxl`; the optional feature groups are:

- `.[optimize]` — DNAChisel, the production codon backend. **Required for
  windowed GC targets**; the dependency-free fallback does not solve them.
- `.[viz]` — matplotlib (PNG plate maps and GC plots; SVG needs nothing) and
  fpdf2 (codon-view PDFs).
- `.[bio]` — Biopython, python-codon-tables, Plateo.
- `.[all]` — everything above plus pytest.

`.xlsx` → legacy `.xls` conversion is best-effort and needs LibreOffice on the
PATH; without it the `.xlsx` is kept and nothing fails.

## Pipelines

### Mutational scans → `docs/mutational_scans.md`

Alanine scan, site-saturation and reduced-alphabet (charge/hydrophobicity)
scans of any protein. Positions are explicit — a list or a `"31-35,96-110"`
spec — and default to the whole sequence. Input is a protein sequence or FASTA;
output is a variant CSV, optionally continuing straight into DNA.
Entry point: `scripts/scan/mutational_scan.py`.

### Nanobody scans → `docs/mutational_scans.md`

The same three scans with nanobody defaults: parents are looked up by ID in the
lab spreadsheet, and the target set defaults to the CDRs, annotated by motif
(`prosy.core.antibody`, IMGT boundaries matching the spreadsheet's own CDR
columns). `--region cdr3 / framework / all` and `--positions` override it.
Entry point: `scripts/nanobodies/nanobody_scan.py`.

### Scan → Golden Gate fragments → `docs/golden_gate.md`

Given a destination plasmid and a Type IIS enzyme, ProSy digests the plasmid,
derives the required overhangs, designs the adapters, codon-optimizes each
variant under those constraints, pads to the vendor minimum, and verifies the
result — including a full simulated assembly whose ORF is translated and
compared to the expected fusion protein. `data/plasmids/FP01.fa` is the worked
example: two outward-facing BsaI sites, a `TATG`/`GGTT` overhang pair, and an
`M–nanobody–GS–HaloTag–FLAG–SNAP` product.
Add `--fragments --destination <plasmid>` to either scan script.

### Nanobody mutation plates → `scripts/nanobodies/make_nanobody_plate.py`

The original pipeline: parents plus their top-N recommended point mutations
laid out as plate columns or rows, with plate maps. Predates the scan modules
and still drives the `260626`/`260720` plate sets.

## Core modules

| Module | Purpose |
|---|---|
| `sequence` | Genetic code, translation, reverse-complement, point-mutation parsing/applying, validation. |
| `enzymes` | Restriction-enzyme registry (BsaI, BsmBI, BbsI, SapI, EcoRI, …) with both-strand patterns and Type IIS cut metadata. |
| `constraints` | Declarative `ConstraintSet` → flat list of forbidden patterns + GC bounds (global or per sliding window) for a backend. |
| `codon` | Codon-usage tables and two interchangeable optimization backends. |
| `optimize` | `optimize_cds()` — protein → synthesis-ready CDS; `build_fragment()` — the optimize → flank → pad composition. |
| `cloning` | Flank/adapter application and minimum-length padding (input preserved verbatim). |
| `scan` | Alanine / site-saturation / reduced-alphabet scans over explicit positions; degenerate-codon helpers. |
| `antibody` | Motif-anchored CDR and framework annotation for VHH/VH domains (IMGT, Kabat, extended). |
| `goldengate` | Type IIS digestion, destination analysis, adapter design and in-silico assembly. |
| `library` | Variants → verified, plated, orderable fragments; mapping/design/plate outputs. |
| `plate` | 96-/384-well coordinate helpers, column- or row-major fill. |
| `layout` | Place groups (parent + N mutants) onto a plate; row/column packing, 96/384, alignment checks. |
| `platemap` | Dependency-free SVG plate-map renderer (colour-coded, labeled); optional PNG via matplotlib/LibreOffice. |
| `io` | Read/write CSV/TSV/XLSX, flexible column matching, vendor plate-sheet writer. |

## Codon backends

`optimize_cds()` selects a backend automatically:

- **DNAChisel** (default when installed) reproduces the original
  `codon_optimize.py`: reverse-translate → enforce translation → avoid
  restriction sites → codon-optimize for a species via a global solver. It is
  also the only backend that solves windowed GC and fixed flank contexts.
- **Highest-frequency fallback** (no dependencies) picks the most frequent
  synonymous codon per residue and greedily swaps codons to clear forbidden
  patterns. It keeps the whole pipeline runnable and testable without
  DNAChisel, but is not a substitute for it on a production run. Windowed GC is
  accepted for interface parity and *not* enforced — `library.verify_library()`
  catches that rather than letting an over-GC fragment reach an order.

## Batch processing

Not yet implemented. The scan CLIs are single-parent, single-scan; running a
panel of parents means a shell loop today. A per-dataset-config batch runner
with fail-isolation, per-run logging and resume is the next piece of work — see
*Future work* in `PR_DESCRIPTION.md`.

## Project status / repo map

| Path | What it is |
|---|---|
| `prosy/core/` | Stable. Covered by `tests/`. |
| `scripts/scan/` | New; the generic scan and fragment CLIs. |
| `scripts/nanobodies/` | Task scripts: `nanobody_scan.py` (new), `make_nanobody_plate.py` and its helpers (established). |
| `data/plasmids/` | Destination vectors, resolved by bare name from `--destination`. |
| `docs/` | One page per pipeline. |
| `tests/` | `python -m pytest tests/` — runs green without DNAChisel installed. |
| `journals/` | Per-session work log. |
| `Working folder/` | Dated run outputs; not code. |

## Tests

```bash
python -m pytest tests/        # or: python tests/test_core.py
```

The suite must pass without any optional dependency installed. Tests that need
DNAChisel are skipped rather than failed.
