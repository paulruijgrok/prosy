# Mutational scans

Turning one protein into a designed set of single mutants, and (optionally) into
synthesis-ready DNA. The protein-level work lives in `prosy.core.scan`; the
antibody-specific position defaults live in `prosy.core.antibody`.

## Position selection

Every scan takes the same `positions` argument and resolves it through
`scan.resolve_positions()`:

| `positions` | meaning |
|---|---|
| `None` (or omitted) | **every residue** of the protein |
| `"31-35,50-65,96-110"` | a spec string of ranges and singletons |
| `[31, 32, 33]` | an explicit list of 1-based positions |

`exclude` takes the same forms and is subtracted afterwards. Positions are
bounds-checked against the protein, so a typo fails immediately instead of
quietly scanning the wrong residues.

## The three scans

### `alanine_scan(protein, positions=None, ...)`

One variant per position, wild type → Ala.

Positions whose wild type *is already* the substitute residue cannot be scanned
with it. By default they are mutated to Gly instead (`native_substitute="G"`);
pass `native_substitute=None` to skip them. `substitute="S"` (or any residue)
turns this into a serine scan, etc.

### `saturation_scan(protein, positions=None, ...)`

Site-saturation: every non-wild-type residue at each position. With the default
20-residue alphabet that is 19 explicitly defined protein sequences per
position — each one an orderable, fully specified clone.

For a library instead of defined clones, use degenerate codons:
`scan.degenerate_saturation_positions()` reports the positions and what a codon
encodes (`NNK` → 20 residues + a stop; `NDT` → 12 residues, no stop), and
`library.build_degenerate_library()` builds one DNA fragment per position.

### `reduced_alphabet_scan(protein, positions=None, alphabet="adklw")`

Saturation restricted to a small alphabet that spans charge, hydrophobicity and
size. Presets:

| preset | residues | rationale |
|---|---|---|
| `adklw` (default) | A D K L W | small-neutral, negative, positive, aliphatic, aromatic-bulky |
| `adrfs` | A D R F S | small, negative, positive (Arg), aromatic, polar-neutral |
| `gdkvy` | G D K V Y | tiny/flexible, negative, positive, β-branched, aromatic-polar |

Any explicit residue string (`alphabet="ADKLW"`) also works. Each variant is
labelled with the property class of the residue substituted in
(`ScanVariant.category`), so results can be grouped by property rather than by
residue identity.

Cost per position: 5 variants (4 where the wild type is already in the
alphabet), against 19 for full saturation.

## CDR annotation

`prosy.core.antibody.annotate()` locates the conserved VH/VHH framework
landmarks by motif and reads the CDRs off as the spans between them. No
dependencies, no numbering server.

Three boundary conventions are available; CDR3 is identical in all of them:

| scheme | CDR1 | CDR2 |
|---|---|---|
| `imgt` (default) | `C1+4 .. W2-3` | `W2+15 .. FR3-9` |
| `kabat` | `C1+9 .. W2-1` | `W2+14 .. FR3-1` |
| `extended` | `C1+4 .. W2-1` | `W2+14 .. FR3-1` |

`imgt` is the default because it reproduces the CDR calls already recorded in
the `CDR 1/2/3` columns of the lab's nanobody spreadsheet —
`tests/test_antibody.py` pins that agreement for four parents.

Annotation raises `AnnotationError` rather than guessing, so a truncated or
non-antibody sequence surfaces immediately and the caller can fall back to
explicit positions.

## Command line

Generic, any protein — default positions are the whole sequence:

```bash
python scripts/scan/mutational_scan.py --protein-file myprotein.fa --scan alanine
python scripts/scan/mutational_scan.py --protein-file myprotein.fa --scan saturation \
    --positions 31-35,50-65
python scripts/scan/mutational_scan.py --protein-file myprotein.fa --scan reduced \
    --alphabet adklw
```

Nanobodies — same scans, default positions are the CDRs:

```bash
python scripts/nanobodies/nanobody_scan.py --parent Nb01 --scan alanine
python scripts/nanobodies/nanobody_scan.py --parent Nb01 --scan saturation --region cdr3
python scripts/nanobodies/nanobody_scan.py --parent Nb01 --scan reduced --region all
```

`--region` accepts `cdr` (default), `cdr1`/`cdr2`/`cdr3`, `framework` or `all`.
An explicit `--positions` always overrides `--region`.

Add `--fragments` to either script to continue into DNA — see
[docs/golden_gate.md](golden_gate.md).

## Outputs

`<stamp>_variants.csv` — one row per variant: `Name`, `Parent`, `Scan`,
`Mutation`, `Position`, `WT`, `Mut`, `Category`, `Length`, `Protein`.

With `--fragments`, the DNA outputs described in
[docs/golden_gate.md](golden_gate.md) are written alongside it.
