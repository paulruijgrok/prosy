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
    --fragments --out-dir /tmp/nb01
```

That writes a plate-upload `.xlsx`, a design CSV and a full mapping CSV, and
prints `All checks passed.` only if every fragment translated correctly, carried
the right overhangs, cleared the vendor's manufacturability limits, and
assembled into the expected fusion protein in silico.

The same three scans on any protein, with the whole sequence as the default
target set:

```bash
python scripts/scan/mutational_scan.py --protein-file myprotein.fa --scan saturation \
    --positions 31-35,50-65 --fragments --destination FP01
```

**→ [docs/cookbook.md](docs/cookbook.md) has tested commands for every pipeline**
— scans, nanobody plates, design-set plates, the inspection tools, the Python
API, and a troubleshooting table.

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
| `synthesis` | Vendor manufacturability: repeat coverage, windowed repeat density and GC, homopolymers; named `SynthesisProfile`s and the gate that fails a run. |
| `plate` | 96-/384-well coordinate helpers, column- or row-major fill. |
| `layout` | Place groups (parent + N mutants) onto a plate; row/column packing, 96/384, alignment checks. |
| `platemap` | Dependency-free SVG plate-map renderer (colour-coded, labeled); optional PNG via matplotlib/LibreOffice. |
| `io` | Read/write CSV/TSV/XLSX, flexible column matching, vendor plate-sheet writer. |

## Synthesis manufacturability — read before ordering DNA

Plain codon optimization is **not** orderable. `use_best_codon` picks the single
most frequent codon for each residue, so every Ala is `GCG` and every Leu
`CTG`, and repeated 8-mers end up covering ~70% of a fragment against a typical
vendor limit of 40%. This is the expected output of maximal CAI optimization,
not bad luck — it rejected a whole 84-fragment order in September 2026, and the
96-well nanobody plate built the same way had 64 wells over the limit.

Every script that produces DNA therefore defaults to
`--synthesis-profile vendor-standard`, which constrains repeats, GC and rare
codons at design time and then **gates the run** on the resulting measurements:

| | plain optimization | `vendor-standard` |
|---|---|---|
| 8-mer repeat coverage (both strands) | 58–78% | 0–8% |
| overall GC | 62–66% | 55–56% |
| rare codons | 0% | ~0.3%, no tandem runs |
| fragments a vendor would reject | most | none |

`prosy.core.synthesis` measures what the vendor measures, so a breach fails the
run instead of surfacing days later at the order desk. Use
`--no-check-synthesis` to inspect a borderline set, and `--synthesis-profile
none` for the old unconstrained behaviour — neither is order-ready. A different
vendor is a new `SynthesisProfile`, not a new code path. See
[docs/golden_gate.md](docs/golden_gate.md).

**Reproducibility:** results are deterministic per `--seed`, and a sequence does
not change when other designs are added to or removed from the batch. Orders
placed before 2026-09-21 predate this and cannot be regenerated — their output
files are the record.

## Codon backends

`optimize_cds()` selects a backend automatically:

- **DNAChisel** (default when installed) reproduces the original
  `codon_optimize.py`: reverse-translate → enforce translation → avoid
  restriction sites → codon-optimize for a species via a global solver. It is
  also the only backend that solves windowed GC, fixed flank contexts, k-mer
  uniqueness and the rare-codon floor — i.e. everything a synthesis profile
  asks for.
- **Highest-frequency fallback** (no dependencies) picks the most frequent
  synonymous codon per residue and greedily swaps codons to clear forbidden
  patterns. It keeps the whole pipeline runnable and testable without
  DNAChisel, but is not a substitute for it on a production run. The
  global-solver constraints are accepted for interface parity and *not*
  enforced — the synthesis gate catches that rather than letting an
  unmanufacturable fragment reach an order.

## Batch processing

Not yet implemented. The scan CLIs are single-parent, single-scan; running a
panel of parents means a shell loop today. A per-dataset-config batch runner
with fail-isolation, per-run logging and resume is the next piece of work — see
*Future work* in `PR_DESCRIPTION.md`.

## Project status / repo map

| Path | What it is |
|---|---|
| `prosy/core/` | Stable. Covered by `tests/`. |
| `scripts/scan/` | The generic scan and fragment CLIs. |
| `scripts/designs/` | Design-set plates: a CSV of designed sequences → verified, plated fragments. |
| `scripts/nanobodies/` | Task scripts: `nanobody_scan.py` (new), `make_nanobody_plate.py` and its helpers (established). |
| `data/plasmids/` | Destination vectors, resolved by bare name from `--destination`. |
| `docs/` | `cookbook.md` (commands for everything), plus one page per pipeline. |
| `tests/` | `python -m pytest tests/` — runs green without DNAChisel installed. |
| `journals/` | Per-session work log. |
| `Working folder/` | Dated run outputs; not code. |

## Tests

```bash
python -m pytest tests/        # or: python tests/test_core.py
```

The suite must pass without any optional dependency installed. Tests that need
DNAChisel are skipped rather than failed.
