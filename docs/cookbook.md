# Cookbook

Copy-pasteable commands for every pipeline in the repo. Every command here has
been run; paths are real paths in this repo.

Run everything from the repo root, in the `DataAnalysis` conda env (it has
DNAChisel, which is required for the synthesis constraints):

```bash
cd ~/Documents/Claude/Projects/ProSy
conda activate DataAnalysis
python -m pytest tests/ -q          # sanity check: expect all green
```

All DNA-producing scripts default to `--synthesis-profile vendor-standard` and
gate the run on it — see [the README](../README.md#synthesis-manufacturability--read-before-ordering-dna).
Nothing below needs to ask for that.

---

## 1. Mutational scans — proteins only

Produces `<stamp>_variants.csv`. No DNA, no destination vector needed.

```bash
# Alanine scan of every residue
python scripts/scan/mutational_scan.py \
    --protein-file "Working folder/260917_designs_enzyme/designs01.fasta" \
    --name designs01 --scan alanine \
    --out-dir /tmp/scan --stamp ala

# Site-saturation (19 variants per position) over two loops
python scripts/scan/mutational_scan.py \
    --protein-file "Working folder/260917_designs_enzyme/designs01.fasta" \
    --name designs01 --scan saturation --positions 100-110,150-155 \
    --out-dir /tmp/scan --stamp ssm

# Reduced alphabet: 5 residues spanning charge/hydrophobicity/size
python scripts/scan/mutational_scan.py \
    --protein-file "Working folder/260917_designs_enzyme/designs01.fasta" \
    --name designs01 --scan reduced --alphabet adklw --positions 100-120 \
    --out-dir /tmp/scan --stamp red
```

A protein can be given inline instead of as a file with `--protein MKT…`.
Positions default to the whole sequence; `--exclude` subtracts from the set.

## 2. Scans → orderable DNA

Add `--fragments` and a destination. Adapters are derived from the plasmid.

```bash
python scripts/scan/mutational_scan.py \
    --protein-file "Working folder/260917_designs_enzyme/designs01.fasta" \
    --name designs01 --scan alanine --positions 100-110 \
    --fragments \
    --destination "Working folder/260917_designs_enzyme/gg002.txt" \
    --flank-5 TGTATCGGTCTCGAGGA --flank-3 GGTTCCGGAGACCTCTAGT \
    --expect-prefix MSG --expect-suffix GSHHHHHH \
    --out-dir /tmp/scan --stamp ala_dna
```

> **The `--flank-*` overrides are not optional for gg002.** Its auto-derived 3'
> adapter omits a `GG` that is *coding* — it spells `GGT` with the vector's `T`
> and keeps the downstream His6 in frame. `--expect-prefix/--expect-suffix`
> make the ORF check catch exactly that. FP01 needs no overrides.

Writes the variants CSV plus `<stamp>_library_mapping.csv`,
`<stamp>_design_seqs.csv` and one `<stamp>_plate[_pN].xlsx` per plate.

## 3. Nanobody scans — CDRs by default

Same three scans; parents come from the lab spreadsheet by ID, and the target
set defaults to the CDRs (IMGT, matching the sheet's own CDR columns).
Destination defaults to FP01.

```bash
# Alanine scan of all three CDRs, through to fragments
python scripts/nanobodies/nanobody_scan.py \
    --parent Nb01 --scan alanine --fragments \
    --out-dir /tmp/nbscan --stamp nb01_ala

# Saturation of CDR3 only
python scripts/nanobodies/nanobody_scan.py \
    --parent Nb01 --scan saturation --region cdr3 --fragments \
    --out-dir /tmp/nbscan --stamp nb01_ssm3

# NNK library: one degenerate fragment per CDR3 position
python scripts/nanobodies/nanobody_scan.py \
    --parent Nb01 --region cdr3 --degenerate NNK --fragments \
    --out-dir /tmp/nbscan --stamp nb01_nnk
```

`--region` takes `cdr` (default), `cdr1`/`cdr2`/`cdr3`, `framework` or `all`;
an explicit `--positions` always wins. `--protein …` skips the spreadsheet.

## 4. Nanobody mutation plates

The original pipeline: each parent plus its top-N recommended mutations fills a
plate column or row.

```bash
DATA="Working folder/260720_NanobodyMuts_gc72"

# 12 parents x 7 mutations, one parent per column (the 260626 layout)
python scripts/nanobodies/make_nanobody_plate.py \
    --data-dir "$DATA" --stamp demo_col --plate-map svg

# Row layout: 5 SetA + 3 SetB parents, 11 mutations each
python scripts/nanobodies/make_nanobody_plate.py \
    --data-dir "$DATA" --stamp demo_row --orientation row \
    --REDA Nb01 Nb02 Nb03 Nb04 Nb05 --REDA Nb50 Nb51 Nb58 \
    --mutations-per-parent 11 --plate-map svg
```

Needs `plasmid_database.xlsx` and `recommendations_<SET>.tsv` in `--data-dir`;
outputs are written there too.

## 5. Design-set plates

A CSV of designed sequences → a verified, plated order. Two layout modes.

```bash
DIR="Working folder/260917_designs_enzyme"

# One design group per plate row (oversized groups spill, balanced 7+7)
python scripts/designs/make_design_plate.py \
    --csv "$DIR/sequences.csv" --destination "$DIR/gg002.txt" \
    --flank-5 TGTATCGGTCTCGAGGA --flank-3 GGTTCCGGAGACCTCTAGT \
    --expect-prefix MSG --expect-suffix GSHHHHHH \
    --out-dir /tmp/designs --stamp by_group

# Input order preserved, column-major (8-channel pipette friendly)
python scripts/designs/make_design_plate.py \
    --csv "$DIR/sequences_bypLDDT.csv" --destination "$DIR/gg002.txt" \
    --layout sequential --orientation column \
    --carry-columns selection_bucket complex_plddt \
    --flank-5 TGTATCGGTCTCGAGGA --flank-3 GGTTCCGGAGACCTCTAGT \
    --expect-prefix MSG --expect-suffix GSHHHHHH \
    --check-gc-min 0.40 --check-gc-max 0.58 \
    --out-dir /tmp/designs --stamp by_plddt
```

The second command is the one that produced the accepted `260921_designs_v3`
order. `--carry-columns` copies extra CSV columns into the mapping so the plate
traces back to the scores. Column names are configurable with `--id-column`,
`--group-column`, `--sequence-column`.

## 6. Inspecting a finished run

All four take a run's mapping CSV.

```bash
MAP="Working folder/260720_NanobodyMuts_gc72/260626_nanobody_mapping.csv"

# Sliding-window GC: long CSV + overlay plot, flagging out-of-bounds windows
python scripts/nanobodies/gc_sliding_window.py \
    --input "$MAP" --window 50 --high 72 --plot /tmp/gc.png --out /tmp/gc.csv

# E. coli rare codons and tandem runs (accepts several runs to compare)
python scripts/nanobodies/rare_codon_scan.py --input "$MAP" --show-runs

# Codon-aligned view: text + colour-coded HTML (+ PDF unless --no-pdf)
python scripts/nanobodies/codon_view.py \
    --input "$MAP" --limit 4 --out-txt /tmp/codons.txt \
    --out-html /tmp/codons.html --no-pdf
```

Turning a run into rows for the lab spreadsheet:

```bash
python scripts/nanobodies/make_db_rows.py \
    --data-dir "Working folder/260720_NanobodyMuts_gc72" --stamp 260626 \
    --out /tmp/db_rows.xlsx
```

## 7. From Python

Analysing a destination vector, then simulating a reaction into it. Runs
as-is from the repo root:

```python
import csv
from prosy.core.goldengate import destination_ends, design_insert_flanks, assemble
from prosy.core.library import orf_protein
from prosy.core import synthesis, antibody

def load(path):
    return "".join(l.strip() for l in open(path) if not l.startswith(">")).upper()

# What does the vector demand of an insert?
fp01 = load("data/plasmids/FP01.fa")
ends = destination_ends(fp01, "BsaI")
print(ends.five_overhang, ends.three_overhang)        # TATG GGTT
print(design_insert_flanks(ends, "BsaI"))             # the adapters to use

# Simulate the reaction with a real fragment and read the product back.
row = next(iter(csv.DictReader(open(
    "Working folder/260917_designs_enzyme/260921_designs_v3_mapping.csv"))))
gg002 = load("Working folder/260917_designs_enzyme/gg002.txt")
product = assemble(gg002, [row["Final_DNA"]], "BsaI")
print(len(product), orf_protein(product)[:12])        # 6249 MSGAAEVTEVEV

# Would a vendor accept that fragment?
report = synthesis.check(row["Final_DNA"])            # or check(seq, my_spec)
print(report.ok, report.summary())
# True 1035 bp  GC 55.9%  repeat8 2.3%  dens 8.9%  GC20max 80%

# Where are a nanobody's CDRs?
ann = antibody.annotate("QVQLVESGGGLVQAGGSLRLREDACTED_SEQUENCEAPGKQREWVATFT"
                        "SSGDANYADSVKGRFTISRDNAKSTVYLQMNSLKPEDTAVYYCNADVYGWGY"
                        "TSYSDYWGQGTQVTVSS")
print(ann.describe())      # CDR1 26-33 REDACTED  CDR2 51-57 REDACTE  CDR3 96-110 …
print(ann.cdr_positions("CDR3"))
```

## 8. Recipes

**Add a destination plasmid.** Drop a FASTA or bare sequence into
`data/plasmids/`. `--destination <name>` then resolves it by bare name; a path
or raw DNA also work.

**Check what a vector requires** before committing to a design:

```bash
python -c "
import sys; sys.path.insert(0,'.')
from prosy.core.goldengate import destination_ends
s=''.join(l.strip() for l in open('data/plasmids/FP01.fa') if not l.startswith('>'))
e=destination_ends(s.upper(),'BsaI')
print('overhangs', e.five_overhang, e.three_overhang)
print('backbone', len(e.backbone), 'bp; drop-out', len(e.dropout), 'bp')
print('upstream ', e.upstream[-40:])
print('downstream', e.downstream[:40])
"
```

**Target a different vendor.** Add a `SynthesisProfile` in
`prosy/core/synthesis.py` with that vendor's thresholds and select it with
`--synthesis-profile <name>`. Don't hand-tune flags per run — that is how the
September 2026 rejection happened.

**Compare against unconstrained output** (diagnosis only, never for ordering):

```bash
python scripts/nanobodies/nanobody_scan.py --parent Nb01 --scan alanine \
    --fragments --synthesis-profile none --out-dir /tmp/legacy --stamp legacy
```

**Reproduce an order.** Re-running the same command with the same `--seed`
gives the same sequences, and a sequence does not change if other designs are
added to or removed from the set. Orders placed before 2026-09-21 predate the
seeding fix — for those, the delivered mapping CSV and FASTA are the record.

**Speed.** ~0.3 s per 333 aa sequence with the full constraint set, so an
84-design plate takes about half a minute. `--assembly-checks 0` skips the
in-silico assembly if you only want sequences; `-1` assembles all of them.

## 9. Troubleshooting

| Symptom | Cause and fix |
|---|---|
| `no ladder step satisfied the constraints` | The protein cannot meet the profile even with the soft fallback. Inspect with `--synthesis-profile none` and loosen a specific constraint. |
| `repeated 8-mers cover …% (limit 40%)` | The gate doing its job. Almost always means the run used `--synthesis-profile none`, or the fallback codon backend ran because DNAChisel is not installed. |
| `assembled ORF does not contain the design` | Wrong adapters for the vector. Check `destination_ends()` output against the flanks in use. |
| `assembled ORF is …, expected …GSHHHHHH` | Reading-frame slip at the 3' junction — for gg002, the missing `GG` (see §2). |
| `Expected exactly one …-free backbone fragment` | The plasmid's Type IIS sites don't face outward from a stuffer, or there are more than two. |
| Windowed GC breaches with no explanation | The dependency-free fallback backend is running. Install `.[optimize]`. |
| `Parents not found in sheet` | `--data-dir` lacks `plasmid_database.xlsx`, or the IDs aren't in the `Nanobodies` sheet. |
