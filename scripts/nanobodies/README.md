# scripts/nanobodies

Task scripts for nanobody DNA-fragment design.

- **`nanobody_common.py`** — shared library: load parents from the lab sheet,
  match recommended mutations by sequence, build parent+mutant variants, lay them
  onto a plate, run the per-variant synthesis pipeline, write outputs, verify.
- **`make_nanobody_plate.py`** — one general CLI over that library. Each parent +
  its top-N mutations forms a group that fills a plate column or row; groups pack
  onto a 96- or 384-well plate.

## make_nanobody_plate.py

**Inputs** (in `--data-dir`, default `Working folder/260626_NanobodyMuts/`):

- `Malaria DX plasmids.xlsx` — sheet `Nanobodies` (parent IDs ↔ AA sequences).
- `recommendations_ADL1.tsv`, `recommendations_KRU1.tsv` — proposed mutations.

Parents are matched to the TSVs by full protein sequence. Mutations are ranked by
`score` (desc, ties broken by lowest position) and the top N kept. Each variant
is codon-optimized (E. coli, BsaI-avoiding by default), flanked, and padded to a
minimum length.

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
