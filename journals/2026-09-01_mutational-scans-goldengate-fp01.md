# 2026-09-01 — Mutational scans + Golden Gate design, worked through FP01

## What changed and why

Two capabilities were missing and they turned out to be one problem.

1. ProSy could optimize *a* protein into *a* fragment, but had no notion of a
   designed *set* of variants — no alanine scan, no saturation library.
2. It had no model of the destination vector. The Golden Gate adapters were
   literal strings in a script's defaults
   (`DEFAULT_FLANK_5 = "TGTATCGGTCTCgAGGA"`, `DEFAULT_FLANK_3 =
   "GGTTCCgGAGACCTCTAGT"`), with nothing tying them to a plasmid.

Adding (1) without (2) would have produced hundreds of fragments carrying
unverified adapters. So the session built both and joined them: scans produce
proteins, the vector produces the adapter requirements, and the pipeline proves
the combination by assembling the plasmid in silico.

**The hard-coded flanks do not fit FP01.** Working the geometry out from the
sequence: FP01's two outward-facing BsaI sites demand overhangs `TATG` (5') and
`GGTT` (3'). The old defaults give `AGGA` and (bottom-strand) `GGAA`. Whatever
vector those were for, it was not this one — which is exactly the class of error
this PR is meant to make impossible.

## Design decisions worth remembering

**Top-strand coordinates for all Type IIS reasoning.** The whole `goldengate`
module rests on one restatement of the cut geometry: *a cut at top-strand index
`p` severs the top strand before `p` and the bottom strand before `p + k`*.
Therefore the fragment to the right of `p` owns the overhang `seq[p:p+k]` in its
own top strand and the fragment to the left does not contain those bases at all
— so ligation is plain concatenation of top strands, valid exactly when one
fragment's right sticky end equals the next one's left sticky end. Digestion,
adapter design and assembly all fall out of that single rule; there is no
separate bookkeeping of "which strand overhangs where".

**The backbone is identified structurally, not positionally.** It is the digest
fragment with no remaining recognition site — which is what makes a drop-out
destination directional and self-selecting in the first place. If a plasmid does
not have exactly one such fragment, `destination_ends()` refuses rather than
picking the largest.

**FP01's start codon lives in the overhang.** The 5' overhang `TATG` completes
the vector's `CATATG` (pET RBS context), so the insert encodes the payload
**without** a leading Met and `len(coding)` must be a multiple of 3 to keep the
downstream GS linker in frame. The 3' overhang `GGTT` opens that linker
(`…GGT TCC GGT TCT GGT…`). Product:
`M–<nanobody>–GSGSGSGSGS–HaloTag–GSGSGSGSGS–DYKDDDDK–GSGSGSGSGS–SNAP`.

**CDR annotation defaults to IMGT because the lab already does.** The motif
anchors (`C1`, the FR2 `WxRQ` Trp, the FR3 `R[FVLIA][TASV][ILVM][ST]` motif,
`C2`, the FR4 `WGxG` Trp) support three conventions; IMGT reproduces the `CDR
1/2/3` columns of the Nanobodies sheet exactly for Nb01/Nb02/Nb50/Nb58, so the
scan targets match the annotations already in use. That agreement is pinned by
a test.

**Position selection is explicit everywhere in the core.** `resolve_positions()`
takes a list, a `"31-35,96-110"` spec or `None` (= whole sequence). Domain
defaults live in the caller, not the engine — which is why "CDRs by default"
lives in `nanobody_scan.py` and not in `scan.py`.

## Two things the verification caught

- **NNK at Nb01 position 104 can spell `GAGACC`** (BsaI on the reverse strand)
  together with its neighbouring codons. Every concrete expansion of a
  degenerate codon is now checked, and a single synonymous swap in a
  neighbouring codon repairs it automatically; if no swap works, the position is
  reported rather than shipped.
- **The fallback codon backend silently ignores `gc_window`.** Building the FP01
  library with it and a 50 bp / 72% cap gives windows at 80%. Not fixed here
  (it's a global-solver problem, which is DNAChisel's job), but made loud:
  `verify_library()` fails those fragments, the CLI exits non-zero, and a
  skip-if-no-DNAChisel test pins both halves of the behaviour.

## Files touched (uncommitted)

New:

- `prosy/core/scan.py` — position resolution + one scan engine; `alanine_scan`,
  `saturation_scan`, `reduced_alphabet_scan`, degenerate-codon helpers.
- `prosy/core/antibody.py` — motif-anchored CDR/framework annotation (imgt /
  kabat / extended).
- `prosy/core/goldengate.py` — `find_cuts`, `digest`, `destination_ends`,
  `design_insert_flanks`, `check_fragment`, `assemble`.
- `prosy/core/library.py` — `LibraryConfig`, `build_library`,
  `build_degenerate_library`, `verify_library`, `orf_protein`, `write_library`.
- `scripts/scan/scan_common.py`, `scripts/scan/mutational_scan.py`,
  `scripts/nanobodies/nanobody_scan.py`.
- `data/plasmids/FP01.fa`; `docs/mutational_scans.md`, `docs/golden_gate.md`.
- `tests/test_scan.py`, `test_antibody.py`, `test_goldengate.py`,
  `test_library.py`.

Modified:

- `prosy/core/__init__.py` — export the four new modules.
- `README.md` — restructured to the standard layout (quick start, installation,
  pipelines overview with links per pipeline, batch processing, repo map).
- `PR_DESCRIPTION.md` — replaced with this PR's description.

## Exact commands

```bash
cd /path/to/ProSy
pip install -e ".[optimize,dev]"

# Tests + hygiene sweep
python -m pytest tests/ -q
ruff check prosy scripts tests --select F401,F841
vulture prosy scripts --min-confidence 80

# The four worked-example runs (Nb01 -> FP01, BsaI, 50 bp / 72% GC cap)
OUT="Working folder/260901_Nb01_scans"
python scripts/nanobodies/nanobody_scan.py --parent Nb01 --scan alanine \
    --fragments --gc-window 50 --gc-max 0.72 --assembly-checks -1 \
    --out-dir "$OUT" --stamp 260901_Nb01_ala_cdr
python scripts/nanobodies/nanobody_scan.py --parent Nb01 --scan reduced --alphabet adklw \
    --fragments --gc-window 50 --gc-max 0.72 \
    --out-dir "$OUT" --stamp 260901_Nb01_reduced_cdr
python scripts/nanobodies/nanobody_scan.py --parent Nb01 --scan saturation --region cdr3 \
    --fragments --gc-window 50 --gc-max 0.72 \
    --out-dir "$OUT" --stamp 260901_Nb01_ssm_cdr3
python scripts/nanobodies/nanobody_scan.py --parent Nb01 --region cdr3 --degenerate NNK \
    --fragments --out-dir "$OUT" --stamp 260901_Nb01_nnk_cdr3

# Generic (non-antibody) equivalent, whole sequence by default
python scripts/scan/mutational_scan.py --protein-file myprotein.fa --scan saturation \
    --positions 31-35,50-65 --fragments --destination FP01 --out-dir /tmp/out
```

Inspecting FP01 directly:

```python
from prosy.core.goldengate import destination_ends, design_insert_flanks, assemble
from prosy.core.library import orf_protein

ends = destination_ends(fp01, "BsaI")        # ('TATG', 'GGTT')
flanks = design_insert_flanks(ends, "BsaI")  # TGTATCGGTCTCATATG / GGTTAGAGACCTCTAGT
orf_protein(assemble(fp01, [fragment], "BsaI"))
```

## Verification run

`python -m pytest tests/` → **94 passed** (green without DNAChisel installed).
`ruff check --select F401,F841` clean, `vulture --min-confidence 80` clean.

Worked example in `Working folder/260901_Nb01_scans/` — Nb01, FP01, BsaI,
DNAChisel, 50 bp / 72% GC cap, all four runs reporting *All checks passed*:

| run | count | fragments |
|---|---|---|
| alanine scan, all CDRs | 31 variants | 31 × 397 bp, 1 plate — all 31 assembled in silico |
| reduced alphabet `adklw`, all CDRs | 144 variants | 144 × 397 bp, 2 plates |
| site-saturation, CDR3 | 286 variants | 286 × 397 bp, 3 plates |
| NNK, CDR3 | 15 positions | 15 × 397 bp, 1 plate |

## Open threads

- Nothing is committed yet; the working tree also carried pre-existing
  uncommitted work (`scripts/nanobodies/make_db_rows.py`,
  `tests/test_make_db_rows.py`, `scripts/nanobodies/README.md`) from an earlier
  session — untouched here, but it needs sorting before a commit.
- `nanobody_common.py` is now a special case of `library.py` and should be
  folded onto it in its own PR, with a regeneration diff of the 260626/260720
  plate sets as the check.
- No batch runner yet; a panel of parents is a shell loop today.
- `library.write_library()` emits no plate maps, though `platemap.py` is right
  there — colouring by scanned CDR or by property class is the obvious mapping.
