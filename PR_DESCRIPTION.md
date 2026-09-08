# Mutational-scan library design and Golden Gate assembly, worked through FP01

## Summary

ProSy could optimize *a* protein into *a* fragment. It could not design a
*set* — an alanine scan, a saturation library — and it had no model of the
vector those fragments were going into: the Golden Gate adapters were literal
strings copied into a script's defaults (`DEFAULT_FLANK_5 = "TGTATCGGTCTCgAGGA"`),
with nothing checking they matched any real plasmid.

This PR adds both halves and joins them: four new core modules that go from one
protein sequence to a verified, plated synthesis order, and a worked example
that designs alanine, saturation, reduced-alphabet and NNK libraries of Nb01
for Golden Gate cloning into FP01.

The adapters are no longer trusted — they are *derived* from the destination
plasmid and then proven by assembling the plasmid in silico and translating its
ORF. Notably, the previously hard-coded flanks (`AGGA`/`TTCC` overhangs) do
**not** match FP01, whose BsaI sites demand `TATG`/`GGTT`; the new pipeline
finds that from the sequence rather than from a comment.

## Changes

### New core modules

- **`prosy/core/scan.py`** — mutational scans over an explicit position set
  (`resolve_positions`: an int list, a `"31-35,96-110"` spec, or `None` = the
  whole sequence, with `exclude` and bounds checking). One engine, three named
  scans: `alanine_scan` (native-Ala positions go to Gly by default, or are
  skipped), `saturation_scan` (19 defined variants per position, custom
  alphabets allowed) and `reduced_alphabet_scan` (5-residue property alphabets;
  presets `adklw` / `adrfs` / `gdkvy`, default `adklw` = small / negative /
  positive / aliphatic / aromatic). Plus IUPAC degenerate-codon helpers
  (`expand_degenerate`, `degenerate_residues`) for library-scale saturation.
- **`prosy/core/antibody.py`** — dependency-free, motif-anchored CDR/framework
  annotation for VHH/VH domains. Three boundary conventions (`imgt` default,
  `kabat`, `extended`); raises `AnnotationError` rather than guessing.
- **`prosy/core/goldengate.py`** — Type IIS digestion, destination analysis,
  adapter design and assembly simulation, all in top-strand coordinates so that
  ligation is literally concatenation of top strands.
  `destination_ends()` identifies the backbone as the digest fragment with no
  remaining recognition site and reports the overhangs an insert must carry;
  `design_insert_flanks()` builds the adapters (picking a spacer base that adds
  no extra site); `assemble()` returns the circular product or raises
  `CloningError` on a mismatched, ambiguous or non-circularising design.
- **`prosy/core/library.py`** — variants → verified, plated fragments.
  Optimize (with a per-protein cache), flank, pad, verify, lay out across as
  many plates as needed, write mapping / design / vendor-plate files.
  `build_degenerate_library()` builds one fragment per position with a
  degenerate codon and checks *every* concrete codon it can take against the
  forbidden patterns.

### New scripts

- **`scripts/scan/scan_common.py`** — shared argument groups, input loading, run
  orchestration and reporting for both CLIs.
- **`scripts/scan/mutational_scan.py`** — generic: any protein, default target
  set is the whole sequence, `--fragments` to continue into DNA.
- **`scripts/nanobodies/nanobody_scan.py`** — the same three scans with nanobody
  defaults: parents looked up by ID in the lab spreadsheet, target set defaults
  to the CDRs (`--region cdr|cdr1|cdr2|cdr3|framework|all`, overridden by an
  explicit `--positions`), destination defaults to FP01.

### Data and docs

- **`data/plasmids/FP01.fa`** — the destination vector, resolvable by bare name
  from `--destination`.
- **`docs/mutational_scans.md`**, **`docs/golden_gate.md`** — one page per
  pipeline; `README.md` restructured to the standard layout (summary, quick
  start, installation, pipelines overview with links, batch processing, repo
  map).

## Verification

- **Full suite: 94 passed**, green with and without DNAChisel installed
  (`tests/test_scan.py`, `test_antibody.py`, `test_goldengate.py`,
  `test_library.py`, plus the pre-existing tests).
- **CDR annotation is pinned against the lab's own calls.** The IMGT boundaries
  reproduce the `CDR 1/2/3` columns of the Nanobodies sheet exactly for Nb01,
  Nb02, Nb50 and Nb58; the `kabat` scheme independently reproduces the canonical
  H1 = 31–35 / H2 = 50–65 for Nb01.
- **Cut geometry is tested from first principles** — that a cut leaves its
  overhang on the downstream fragment and not the upstream one, on both strands.
- **FP01 is characterised from the sequence, not from prior belief**: exactly two
  outward-facing BsaI sites, a 1437 bp chloramphenicol stuffer, a 3788 bp
  backbone, and overhangs `TATG` / `GGTT` — the `ATG` completing the vector's
  `CATATG` start context and the 3' end opening the GS linker.
- **The design is proven by assembly, not by inspection.** Assembling FP01 with
  an Nb01 fragment and translating the product's ORF yields
  `M–<Nb01>–GSGSGSGSGS–HaloTag–GSGSGSGSGS–DYKDDDDK–GSGSGSGSGS–SNAP`, and the
  product contains zero BsaI sites. Negative controls are tested too: a
  one-base overhang change and an internal BsaI site both raise `CloningError`.
- **Worked example, all four runs passing every check**
  (`Working folder/260901_Nb01_scans/`, DNAChisel backend, 50 bp / 72% GC cap):

  | run | variants | fragments | plates |
  |---|---|---|---|
  | alanine scan, all CDRs | 31 | 31 × 397 bp (all 31 assembled in silico) | 1 |
  | reduced alphabet `adklw`, all CDRs | 144 | 144 × 397 bp | 2 |
  | site-saturation, CDR3 | 286 | 286 × 397 bp | 3 |
  | NNK library, CDR3 | 15 positions | 15 × 397 bp | 1 |

- **The degenerate-codon safety check caught a real problem**: `NNK` at Nb01
  position 104 can spell `GAGACC` (BsaI, reverse strand) together with its
  neighbours. It is now repaired automatically by a single synonymous swap in a
  neighbouring codon, and reported rather than shipped if no swap works.
- **Hygiene sweep**: `ruff check --select F401,F841` clean, `vulture
  --min-confidence 80` clean, no TODO/FIXME markers, no commented-out code.

## Known limitation surfaced (not introduced)

The dependency-free fallback backend accepts `gc_window` for interface parity
but does not enforce it. Building the FP01 library with the fallback and a 50 bp
/ 72% cap produces windows at 80%. This PR does not fix the fallback — that is a
global-solver problem — but it does make the failure loud:
`verify_library()` reports every breaching fragment, the CLI exits non-zero, and
`tests/test_library.py::test_windowed_gc_cap_needs_dnachisel` pins both halves
(fallback fails verification, DNAChisel passes). The README and
`docs/golden_gate.md` now say so explicitly.

## Future work

- **Batch runner.** The scan CLIs are single-parent, single-scan; a panel means a
  shell loop today. Needs the standard per-dataset-config runner with
  fail-isolation, per-run logging and resume.
- **Fold `nanobody_common.py` onto `library.py`.** Its `optimize_all` / `verify`
  / `write_outputs` are now a special case of the general pipeline
  (`make_nanobody_plate.py` predates it). Deliberately left untouched here to
  keep the established plate scripts byte-stable; worth doing as its own PR with
  a regeneration diff of the 260626/260720 plate sets.
- **Plate maps for scan libraries.** `library.write_library()` does not yet emit
  the SVG/PNG plate maps that `platemap.py` provides for the mutation plates —
  colouring by scanned CDR or by property class would be the natural mapping.
- **Multi-insert assemblies.** `assemble()` already chains an arbitrary number of
  fragments by overhang, but nothing upstream designs a multi-part assembly
  (e.g. splitting a long CDS across two fragments with an internal junction).
- **Windowed GC in the fallback backend**, so the whole pipeline is honest
  without DNAChisel rather than merely loud about it.
- **More destinations.** `data/plasmids/` holds only FP01; other Golden Gate
  destinations should be added as they come into use.
