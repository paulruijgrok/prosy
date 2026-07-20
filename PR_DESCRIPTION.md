# Add codon-view and rare-codon inspection tools; tighten GC default to 72%

## Summary

Follow-up to the windowed GC cap. Once the cap was in place we needed to (a) read
the resulting sequences at the codon level and (b) confirm the cap wasn't degrading
codon quality (e.g. swapping preferred glycine codons for the rare `GGA`). This PR
adds two inspection tools and, on the strength of what they showed, tightens the
production GC cap from 75% to **72%** — verified to introduce no rare codons.

## Changes

- **`scripts/nanobodies/codon_view.py`** — new. Codon-aligned view of the in-frame
  `Coding_DNA`: each codon with the amino acid beneath it, a residue-number ruler,
  wrapped into blocks. Writes three formats — monospaced **text**, colour-coded
  **HTML**, and colour-coded **PDF** (via `fpdf2`, easiest to share) — amino acids
  grouped by chemical class, and verifies every sequence translates back to its
  `Protein` column. Flags `--style pipe|space`, `--codons-per-line`, `--no-pdf`.
- **`scripts/nanobodies/rare_codon_scan.py`** — new. Flags E. coli rare codons
  (AGA/AGG/CGA/CGG Arg, ATA Ile, CTA Leu, CCC Pro, GGA Gly; GGG as "watch") and
  reports tandem runs (≥2 consecutive rare codons) per plate.
- **`make_nanobody_plate.py` / `nanobody_common.py`** — production default
  `--gc-max` changed from 75 to **72** (`RunConfig.gc_max = 0.72`).
- **`pyproject.toml`** — `fpdf2>=2.7` added to the `viz` and `all` extras.
- **Tests** — `tests/test_codon_view.py`, `tests/test_rare_codon_scan.py`.
- **Docs** — `scripts/nanobodies/README.md` documents both inspection tools, the
  three codon-view outputs, and the 72% default.

## Verification

- Full test suite passes (`pytest`, **35 passed**), including the new tool tests.
- Rare-codon scan across the June 26 plate, the 74% run, and the new 72% run:
  **0 rare codons and 0 tandem runs in all three** — the GC cap costs nothing on
  codon quality (it sheds GC via common synonyms, never rare codons).
- Regenerated the 72% order: worst 50 bp window **72.0%**; all 96 variants
  translate cleanly; codon-view text/HTML/PDF render correctly (first PDF page
  inspected visually).
- Investigated the GC-dense glycine run at residues 8–10: under a 50 bp window
  cap it stays the optimal `GGC-GGC-GGC`; tightening 74%→72% for Nb01 changed a
  single neighbouring Ala codon (`GCG→GCA`), not the glycine run.

## Future work

- The rare-codon set is the classic Rosetta-supplemented group; `GGG` is treated
  as "watch" only. The set / thresholds could be made configurable if needed.
- The codon-view PDF uses Latin-1 core fonts (fine for ASCII sequence data); a
  Unicode TTF would be needed to render any non-Latin-1 metadata.
