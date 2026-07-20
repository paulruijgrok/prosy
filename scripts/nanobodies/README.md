# scripts/nanobodies

Task scripts for nanobody DNA-fragment design.

- **`nanobody_common.py`** — shared library: load parents from the lab sheet,
  match recommended mutations by sequence, build parent+mutant variants, lay them
  onto a plate, run the per-variant synthesis pipeline, write outputs, verify.
- **`make_nanobody_plate.py`** — one general CLI over that library. Each parent +
  its top-N mutations forms a group that fills a plate column or row; groups pack
  onto a 96- or 384-well plate.
- **`gc_sliding_window.py`** — inspection tool: sliding-window GC over a mapping
  CSV → long-format CSV + overlay plot (flags windows outside bounds).
- **`codon_view.py`** — inspection tool: codon-aligned view of each variant's
  coding DNA with the amino acid under each codon, as plain text + colour-coded HTML.

## make_nanobody_plate.py

**Inputs** (in `--data-dir`, default `Working folder/260626_NanobodyMuts/`):

- `Malaria DX plasmids.xlsx` — sheet `Nanobodies` (parent IDs ↔ AA sequences).
- `recommendations_ADL1.tsv`, `recommendations_KRU1.tsv` — proposed mutations.

Parents are matched to the TSVs by full protein sequence. Mutations are ranked by
`score` (desc, ties broken by lowest position) and the top N kept. Each variant
is codon-optimized (E. coli, BsaI-avoiding by default), flanked, and padded to a
minimum length.

By default the codon optimizer also holds the local GC content at or below 72%
over every 50 bp sliding window (`--gc-window`/`--gc-max`). The coding region is
optimized *inside its fixed flank context*, so the cap holds across the
flank/coding junction at the 5' end (a known GC hot spot), not just within the
coding region. The run self-verifies this on the final fragment. Set
`--gc-window 0` to disable and reproduce the pre-cap behaviour.

**IDs:** each mutant gets a new sequential `NbID` (default series from `Nb73`,
assigned in plate-fill order) plus a descriptive `NbID_parent` like `Nb01_T27F`.
Parents keep their own ID. Plate `Name` = `"<NbID> <NbID_parent>"` for mutants.

**Outputs** (prefixed by `--stamp`, so runs don't overwrite each other):

- `<stamp>_nanobody_plate.xlsx` (+ `.xls` if LibreOffice is present)
- `<stamp>_design_seqs_nanobodies.csv` — intermediate (same schema as the
  miniprotein design file)
- `<stamp>_nanobody_mapping.csv` — full provenance
- `<stamp>_nanobody_platemap.svg` (and/or `.png`) — colour-coded plate map

### Options

| Flag | Default | Purpose |
|---|---|---|
| `--adl1` / `--kru1` | Nb01–Nb07 / Nb50,51,52,53,58 | Parent IDs per set, in plate order. |
| `--orientation` | `column` | `column` = one group per column; `row` = one group per row. |
| `--plate-size` | `96` | `96` or `384`. |
| `--plate-map` | `svg` | Plate-map image: `none`/`svg`/`png`/`both` (SVG needs no deps). |
| `--mutations-per-parent` | `7` | Top-N mutations by score (group size = N+1). |
| `--allow-straddle` | off | Permit groups to cross row/column boundaries. |
| `--flank-5` / `--flank-3` | `TGTATCGGTCTCgAGGA` / `GGTTCCgGAGACCTCTAGT` | Flanks. |
| `--species` | `e_coli` | Codon-optimization species. |
| `--avoid-enzymes` | `BsaI` | Enzyme sites to avoid in coding regions. |
| `--min-length` | `300` | Minimum fragment length (bp). |
| `--gc-window` | `50` | Sliding-window width (bp) for the local GC cap; `0` disables. |
| `--gc-max` | `72` | Max GC%% allowed in any `--gc-window` (over a 50 bp window, exactly ≤36 GC bases). |
| `--gc-min` | `0` | Min GC%% required in any `--gc-window` (`0` = no floor). |
| `--new-id-start` / `--new-id-prefix` | `73` / `Nb` | New mutant ID series. |
| `--backend` | `auto` | `auto` / `dnachisel` / `highest_frequency`. |
| `--stamp` | `260626` | Output filename prefix. |
| `--seed` | `1` | RNG seed. |

### Examples

```bash
# Original column plate: 12 parents x 7 mutations, one parent per column
python make_nanobody_plate.py

# Row plate: 5 ADL1 + 3 KRU1, 11 mutations each, one parent per row
python make_nanobody_plate.py --orientation row \
  --adl1 Nb01 Nb02 Nb03 Nb04 Nb05 --kru1 Nb50 Nb51 Nb58 \
  --mutations-per-parent 11 --stamp 260626_row

# 16 parents x 5 mutations, 2 groups per row
python make_nanobody_plate.py --orientation row --mutations-per-parent 5 \
  --adl1 Nb01 Nb02 Nb03 Nb04 Nb05 Nb06 Nb07 Nb08 \
  --kru1 Nb50 Nb51 Nb52 Nb53 Nb54 Nb55 Nb56 Nb57 --stamp 16x5

# 24 parents x 3 mutations, 2 groups per column
python make_nanobody_plate.py --orientation column --mutations-per-parent 3 \
  --adl1 Nb01 Nb02 Nb03 Nb04 Nb05 Nb06 Nb07 Nb08 Nb09 Nb10 Nb11 Nb12 \
  --kru1 Nb50 Nb51 Nb52 Nb53 Nb54 Nb55 Nb56 Nb57 Nb58 Nb59 Nb60 Nb61 \
  --stamp 24x3

# 384-well: 48 parents x 7 mutations, one parent per column
python make_nanobody_plate.py --plate-size 384 --stamp 384plate \
  --adl1 <up to 48 ids ...> --kru1 <...>
```

The script self-verifies every run (unique wells, translation round-trip,
clean coding regions, intact flanks, contiguous ID series, min length) and exits
non-zero if any check fails.

## Inspection tools

Both read a `<stamp>_nanobody_mapping.csv` and write outputs next to it.

### gc_sliding_window.py

Sliding-window GC content across each variant's DNA. Writes a long-format CSV
(`id, well, window_start, window_end, gc`) and an overlay plot; prints a summary
flagging any variant with a window outside `[--low, --high]`% GC.

```bash
python gc_sliding_window.py --input <mapping.csv> --window 50 \
  --out gc_windows.csv --plot gc_profile.png        # --column Coding_DNA to skip flanks
```

### codon_view.py

Codon-aligned view of the in-frame coding sequence (`Coding_DNA`) for manual
reading: each codon with the amino acid it encodes printed directly beneath it,
plus a residue-number ruler, wrapped into blocks. Writes **three** files — a
monospaced plain-text file, a colour-coded HTML file, and a colour-coded PDF
(easiest to share) — with amino acids grouped by chemical class, and verifies
every sequence translates back to its `Protein` column.

```bash
python codon_view.py --input <mapping.csv>     # -> <mapping>_codons.txt + .html + .pdf
```

The PDF needs `fpdf2` (`pip install fpdf2`, or `pip install -e ".[viz]"`); if it's
not installed the PDF is skipped with a note and the text/HTML are still written.

| Flag | Default | Purpose |
|---|---|---|
| `--seq-column` | `Coding_DNA` | In-frame coding-DNA column to render. |
| `--codons-per-line` | `20` | Codons per wrapped block. |
| `--style` | `pipe` | Text layout: `pipe` (`\|CAG\|GTG\|`) or `space` (`CAG GTG`). |
| `--out-txt` / `--out-html` / `--out-pdf` | next to input | Override output paths. |
| `--no-pdf` | off | Skip the PDF (text + HTML only). |
| `--limit` | `0` | Render only the first N variants (0 = all). |
