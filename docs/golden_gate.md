# Golden Gate: from a protein set to ordered fragments

`prosy.core.goldengate` reads a destination plasmid and works out what an insert
must look like; `prosy.core.library` builds, verifies and plates the fragments.

## The one rule everything rests on

All Type IIS reasoning here happens in **top-strand coordinates**:

> A cut at top-strand index `p` severs the top strand before `p` and the bottom
> strand before `p + k` (`k` = overhang length).

So the fragment to the *right* of `p` carries the overhang `seq[p : p+k]` in its
own top strand, and the fragment to the *left* does not contain those bases at
all. Ligation is therefore plain concatenation of top strands, valid exactly
when one fragment's right-hand sticky end equals the next one's left-hand
sticky end. `digest()`, `assemble()` and `design_insert_flanks()` are all
direct consequences of that.

## Reading a destination

```python
from prosy.core.goldengate import destination_ends, design_insert_flanks

ends = destination_ends(fp01_sequence, "BsaI")     # circular by default
ends.five_overhang, ends.three_overhang            # ('TATG', 'GGTT')
ends.upstream, ends.downstream                     # vector context either side
```

The backbone is identified as the digest fragment that no longer carries a
recognition site — the sites leave with the drop-out stuffer, which is what
makes the reaction directional and self-selecting. If a plasmid does not have
exactly one such fragment, `destination_ends()` raises `CloningError` rather
than guessing.

## Designing adapters

```python
flanks = design_insert_flanks(ends, "BsaI")
# 5': TGTATC GGTCTC A TATG        3': GGTT A GAGACC TCTAGT
```

Layout of the ordered fragment:

```
outer_5 + SITE + N + five_overhang + <coding> + three_overhang + N + rc(SITE) + outer_3
```

The overhangs are part of the adapters, so `<coding>` is exactly the sequence
that sits between the two vector junctions. The single spacer base is chosen so
that each adapter contains exactly one recognition site.

## FP01 specifically

FP01 carries two outward-facing BsaI sites around a chloramphenicol stuffer
(1437 bp), leaving a 3788 bp backbone and this pair of overhangs:

| end | overhang | what it is |
|---|---|---|
| 5' | `TATG` | the `ATG` completes the vector's `CATATG` start context — **the start codon lives in the overhang** |
| 3' | `GGTT` | opens the GS linker: the junction reads `…GGT TCC GGT TCT GGT…` |

Two consequences for the insert:

1. The coding region encodes the payload **without a leading Met** — the vector
   supplies it. For a nanobody, `<coding>` is simply the nanobody CDS.
2. `len(<coding>)` must be a multiple of 3 to keep the downstream linker,
   HaloTag, FLAG and SNAP domains in frame.

The assembled product expresses:

```
M–<nanobody>–GSGSGSGSGS–HaloTag–GSGSGSGSGS–DYKDDDDK–GSGSGSGSGS–SNAP-tag
```

## Proving the design

`assemble()` simulates the reaction and returns the circular product; the
library verifier translates its ORF and checks the expected protein is in
there. This is the check that catches a wrong overhang, a frame slip or a site
that survived codon optimization — before anything is ordered.

```python
from prosy.core.goldengate import assemble
from prosy.core.library import orf_protein

product = assemble(fp01_sequence, [fragment], "BsaI")
orf_protein(product)          # 'MQVQLVESGGG…SNAP…'
```

`assemble()` raises `CloningError` when the overhangs do not chain into a single
closed circle, when an overhang is ambiguous, or when the product still contains
a recognition site.

## What the library pipeline verifies

Every fragment, every run:

- the coding region translates back to the intended protein;
- no forbidden site (the cloning enzyme plus anything in `--avoid`) inside the
  coding region;
- the coding region and both adapters survive padding intact;
- the excised insert's overhangs match what the destination demands;
- sliding-window GC stays inside `--gc-min`/`--gc-max` when `--gc-window` is set;
- no duplicate wells;
- and, for `--assembly-checks` fragments (`-1` = all), a full in-silico assembly
  whose ORF contains the expected protein.

Nothing is written as "passed" unless every one of those holds; the CLI exits
non-zero and lists the failures otherwise.

> **Windowed GC needs DNAChisel.** The dependency-free fallback backend enforces
> GC over the whole coding region only. Verification catches its failures rather
> than letting an over-GC fragment through, but for production runs with
> `--gc-window`, install `.[optimize]`.

## Degenerate-codon libraries

`library.build_degenerate_library()` optimizes the parent CDS once and replaces
one codon per fragment with a degenerate codon (`NNK`, `NDT`, …). Every concrete
codon the degenerate symbol can take is checked against the forbidden patterns;
when one of them would create a site, a single synonymous swap in a neighbouring
codon is used to remove the risk, and if no swap works the position is reported
rather than shipped. (This is not hypothetical: `NNK` at Nb01 position 104 can
spell `GAGACC` with its neighbours, and gets repaired automatically.)

## Outputs

| file | contents |
|---|---|
| `<stamp>_library_mapping.csv` | one row per fragment: plate, well, name, mutation, category, lengths, GC, protein, coding DNA, final DNA |
| `<stamp>_design_seqs.csv` | `Name` / `Sequence` / `5' nucleotides` / `3' nucleotides` |
| `<stamp>_plate[_pN].xlsx` | vendor upload sheet, `Well Position` / `Name` / `Sequence`; one file per 96- or 384-well plate |

## Adding a destination plasmid

Drop a FASTA (or plain sequence) file into `data/plasmids/`. `--destination
<name>` resolves bare names against that directory, and also accepts a path or
raw DNA on the command line.
