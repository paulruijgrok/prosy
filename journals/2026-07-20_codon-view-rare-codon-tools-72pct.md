# 2026-07-20 — Codon-view + rare-codon inspection tools, GC default → 72%

Second session of the day (follows
`2026-07-20_windowed-gc-cap-nanobody-plate.md`). Focus: tools to *read* and
*quality-check* the GC-capped sequences, and a decision to tighten the cap.

## What changed and why

After the windowed GC cap landed, two needs emerged: (1) read the ordered
sequences at the codon level (which codon codes which residue), and (2) confirm
the cap wasn't trading preferred codons for rare ones (notably glycine
`GGC → GGA`). Built two inspection tools, used them to investigate, and on the
evidence tightened the production cap 75% → 72%.

Key findings from the investigation:

- **Rare codons: none.** Across the June 26 plate, the 74% run, and the 72% run,
  the rare-codon scan found **0 rare codons and 0 tandem runs** (11,464 codons
  each). DNAChisel's `CodonOptimize` objective sheds GC via common synonyms
  (e.g. Gly `GGC → GGT`, the preferred low-GC codon), never rare ones.
- **Residue 8–10 glycine run is not troublesome under a window cap.** Three
  consecutive glycines are 100% GC locally (every Gly codon is `GG*`), but a
  50 bp *window* average doesn't force breaking them — they stay optimal
  `GGC-GGC-GGC` even at 72%. Tightening Nb01 74% → 72% changed exactly one
  neighbouring Ala codon (`GCG → GCA`). Local (short-window) GC is the only thing
  that would touch residues 8–10, and even then the pick would be `GGT` not `GGA`.
- On that basis, **72% is the new production default** (still introduces no rare
  codons; worst 50 bp window sits at exactly 72.0% = 36/50 GC bases).

## Commits

- `4665c71` — Add codon-view and rare-codon inspection tools; tighten GC default
  to 72%.
- `e42e44f` — Remove unused imports in core (cloning `codon_mod`, layout
  `PLATE_384`, platemap `field`) — pre-existing ruff F401, separate hygiene commit.

Both on `main`; **pushed? pending** (push from a terminal with SSH — sandbox
can't reach GitHub).

## Files touched

- `scripts/nanobodies/codon_view.py` (new) — codon-aligned view of `Coding_DNA`:
  codon + amino acid beneath, residue ruler, wrapped blocks. Writes **text +
  colour-coded HTML + colour-coded PDF** (PDF via `fpdf2`); verifies translation
  vs `Protein`. Flags `--style pipe|space`, `--codons-per-line`, `--no-pdf`,
  `--out-{txt,html,pdf}`.
- `scripts/nanobodies/rare_codon_scan.py` (new) — flags E. coli rare codons
  (AGA/AGG/CGA/CGG/ATA/CTA/CCC/GGA; GGG = watch) and tandem runs (≥2 in a row).
- `scripts/nanobodies/make_nanobody_plate.py`, `nanobody_common.py` — default
  `--gc-max` 75 → 72 (`RunConfig.gc_max = 0.72`).
- `pyproject.toml` — `fpdf2>=2.7` added to `viz` and `all` extras.
- `tests/test_codon_view.py`, `tests/test_rare_codon_scan.py` (new).
- `scripts/nanobodies/README.md` — documents both tools, three codon-view
  outputs, 72% default.
- `prosy/core/{cloning,layout,platemap}.py` — removed unused imports (e42e44f).

## Exact commands (run in the `DataAnalysis` conda env)

```bash
conda activate DataAnalysis
cd "/Users/paulruijgrok/Documents/Claude/Projects/ProSy"
python -m pytest tests/ -q                       # 35 passed

cd scripts/nanobodies
# Production order at 72% (inputs copied into the new dir first)
python make_nanobody_plate.py --backend dnachisel --gc-window 50 --gc-max 72 \
  --stamp 260626 --data-dir "../../Working folder/260720_NanobodyMuts_gc72"

# Codon views (text + HTML + PDF) and rare-codon scan
python codon_view.py --input "../../Working folder/260720_NanobodyMuts_gc72/260626_nanobody_mapping.csv"
python rare_codon_scan.py --show-runs \
  --input "../../Working folder/260626_NanobodyMuts/260626_nanobody_mapping.csv" \
          "../../Working folder/260720_NanobodyMuts_gc72/260626_nanobody_mapping.csv"
```

Output dirs: `Working folder/260720_NanobodyMuts_gc72/` (production, 72%),
`..._gc75/` (74% interim), `260626_NanobodyMuts/` (June 26 original) — all have
`*_codons.{txt,html,pdf}`.

## Open items / next steps

- **Push `main`** to GitHub (`git push origin main`) — commits `4665c71`,
  `e42e44f` are local only.
- Rare-codon set is the classic Rosetta-supplemented group; `GGG` is watch-only.
  Make the set / thresholds configurable if a project needs it.
- Codon-view PDF uses Latin-1 core fonts (fine for ASCII sequence data); embed a
  Unicode TTF if non-Latin-1 metadata ever needs rendering.
- `gc75` (74%) dir is now superseded by `gc72` as production; keep or prune.
