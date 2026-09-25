# 2026-09-21 — Design plates for gg002, and beating a vendor's complexity score

## What changed and why

Two threads, one week. First: turn lists of designed enzyme sequences into
ordered plates for Golden Gate into **gg002**. Second — the substantial one —
the first order was **rejected by the DNA manufacturer** on complexity, and
fixing that turned out to be a flaw in how we codon-optimize, not a quirk of
those particular sequences.

### gg002 and the enzyme design designs

`gg002` is a pET-type destination: T7 → lacO → RBS → `CATATG`, two
outward-facing BsaI sites around a 430 bp stuffer, 5244 bp backbone. Required
insert overhangs are `AGGA` / `TTCC` — which is what ProSy's historical
`DEFAULT_FLANK_5/3` were built for, so those flanks are correct *here* (they do
not fit FP01, whose sites demand `TATG`/`GGTT`; the two vectors are not
interchangeable and the flanks are not either).

**The `GG` in `GGTTCCgGAGACCTCTAGT` is coding, not padding.** It pairs with the
vector's `T` to spell `GGT` (Gly) and is what keeps the downstream His6 in
frame. Dropping it — which is exactly what the auto-derived adapter from
`design_insert_flanks()` does — yields an ORF reading `…PPPPPPLMRSGC` (349 aa)
instead of `…GSHHHHHH` (344 aa): a silently untagged protein that every other
check would pass. That near-miss is why `--expect-prefix/--expect-suffix` now
compares the *whole* assembled ORF rather than just asserting the design is
somewhere inside it. Confirmed the check fires by deliberately running with the
auto adapters.

Product: `M–SG–<design>–GS–HHHHHH–stop`, 6249 bp.

### The vendor rejection, and its actual cause

The manufacturer flagged, for most of the 84 sequences:

1. repeated 8-mers covering ~70% of the sequence (limit 40%)
2. >88% of a 90 bp window being repeat
3. a 20 bp window at 95% GC (limit 90%)

Reproducing their numbers exactly was the first job, and the key detail is that
**repeat coverage counts both strands** — a k-mer and its reverse complement
are the same repeat. Single-stranded counting gives 59.2% for the sequence they
scored at 70.2%, which would have quietly passed the 40% limit and sent us
chasing the wrong thing.

The cause is not the designs. It is `CodonOptimize(method="use_best_codon")`:
93.4% of codons were the single most frequent E. coli codon, so every Ala is
`GCG` and every Leu `CTG`, and 8-mers recur constantly. **~70% repeat coverage
is the expected output of maximal CAI optimization**, not bad luck. Any future
order built the old way would have been rejected the same way.

### The fix, and the trade-off it forced

`UniquifyAllKmers(k=8, include_reverse_complement=True)` as a hard constraint
is the main lever. But breaking repeats means using non-optimal synonyms, which
pulls in rare codons — so `AvoidRareCodons(0.10, species="e_coli")` holds that
line. Together: repeats 70% → ~3%, rare codons 0 → 0.27% of codons with **zero
tandem runs** (well under natural E. coli levels; no Rosetta needed).

**k=8 uniqueness is not always satisfiable.** For 8 of 84, the protein carries a
repeated peptide motif that forces a repeated 8-mer whatever the codons.
Raising to k=9 left them at 45–56% — still rejectable. The answer was to keep a
hard k=9 *and* add k=8 uniqueness as a weighted **objective**: 2.3–5.2%. Hence
the relaxation ladder `k=8 → (k=9 + soft k=8) → (k=10 + soft k=8) → soft only`,
with the last step never infeasible. Which step each design landed on is
recorded in the mapping CSV.

### Round three: overall GC

v2 passed complexity but most sequences exceeded the vendor's 58% overall-GC
guidance (57.8–62.5%). There was never any need for that: the theoretical GC
floor for these proteins is ~41% even keeping the rare-codon guard. 60% was
again just `use_best_codon` preferring GC-rich codons. A whole-fragment
`EnforceGCContent(0.45–0.56)` band fixed it.

Two judgement calls worth keeping:

- **Aim at 0.56, not 0.58.** DNAChisel satisfies a cap by sitting against it, so
  targeting 0.58 delivers ~57.9% — passing with no margin if the vendor rounds
  differently or scores the insert alone rather than the whole fragment. Both
  readings now pass (fragment 55.7–55.9%, CDS 55.8–56.0%).
- **Raise the 20 bp GC *floor* to 25%.** On the first v3 attempt 54 of 84 pressed
  against a 20% floor: pushing GC down globally creates AT-rich holes, the same
  problem in the other direction. At a 25% floor everything is still feasible,
  and fewer designs needed the relaxed ladder step (6, down from 13).

**v3 was submitted and passed with no complexity flags at all.**

## Design decisions worth remembering

- **Measure what the vendor measures, then gate on it.** `prosy.core.synthesis`
  exists so the pipeline computes the vendor's own metrics and *fails the run*
  rather than reporting success and letting a rejectable fragment reach an
  order. The regression test carries the actual rejected sequence and asserts
  the vendor's three published numbers, so a drift in the metric breaks the
  build rather than the order.
- **Hard constraint first, soft objective as fallback.** A hard constraint gives
  a guarantee; a soft one never fails. The ladder gets the guarantee where it is
  achievable and the best-effort result where it is not — and the measured gate,
  not the ladder step, decides whether the run passes.
- **Layout: `layout_lines()` for one variable-size group per row** (used for the
  51-sequence set, where the 14-member group balances 7+7), and sequential fill
  for a pre-ranked list (the 84-sequence pLDDT set). Column-major sequential is
  the 8-channel-pipette layout.
- `sequences_bypLDDT.csv` is **not** globally sorted by pLDDT — it is
  design_group → selection_bucket (RS/R/S) → `complex_plddt` descending within
  each of 12 blocks. File order was preserved as asked; worth knowing before
  reading the plate as a ranking.

## Exact commands

```bash
cd ~/Documents/Claude/Projects/ProSy
conda activate DataAnalysis
python -m pytest tests/ -q          # 119 passed

# The final, accepted order (v3)
python scripts/designs/make_design_plate.py \
    --csv "data/runs/enzyme_designs/sequences_bypLDDT.csv" \
    --destination "data/runs/enzyme_designs/gg002.txt" \
    --layout sequential --orientation column \
    --carry-columns selection_bucket calpha_plddt ni_plddt complex_plddt \
    --flank-5 TGTATCGGTCTCGAGGA --flank-3 GGTTCCGGAGACCTCTAGT \
    --expect-prefix MSG --expect-suffix GSHHHHHH \
    --gc-window 50 --gc-max 0.68 --gc-min 0.30 \
    --gc-window-2 20 --gc-window-2-max 0.85 --gc-window-2-min 0.25 \
    --gc-global-max 0.56 --gc-global-min 0.45 \
    --max-homopolymer-a 8 --max-homopolymer-g 5 \
    --unique-kmer 8 --min-codon-frequency 0.10 \
    --check-synthesis --max-repeat-fraction 0.40 --max-gc-20 0.90 \
    --check-gc-min 0.40 --check-gc-max 0.58 \
    --out-dir "data/runs/enzyme_designs" --stamp 260921_designs_v3
```

Checking any sequence against the vendor rule set:

```python
from prosy.core import synthesis
r = synthesis.check(fragment)      # or synthesis.check(seq, my_spec)
r.ok, r.problems, r.summary()
```

## Results

| metric | v1 | v2 | v3 (accepted) | limit |
|---|---|---|---|---|
| overall GC | 62.1–65.9% | 57.8–62.5% | **55.7–55.9%** | <58% |
| 8-mer repeat coverage | 57.9–78.3% | 0.0–6.0% | **0.0–7.6%** | <40% |
| densest 90 bp window | up to 100% | 26.7% | **26.7%** | <88% |
| max 20 bp GC | up to 100% | 85% | **85%** | <90% |
| min 20 bp GC | — | 10% | **25%** | — |
| rare codons | 0% | 0.21% | **0.27%**, 0 tandem | — |
| would be rejected | 84/84 | 83/84 | **0/84** | — |

All three versions: translations correct, zero BsaI in any coding region, all
products 6249 bp, every ORF exactly `MSG + design + GSHHHHHH`.

## Files touched (uncommitted)

New:

- `prosy/core/synthesis.py` — vendor manufacturability metrics (`SynthesisSpec`,
  `repeat_coverage`, `repeat_fraction`, `densest_repeat_window`,
  `extreme_gc_window`, `check`).
- `scripts/designs/make_design_plate.py` — CSV of designs → verified, plated,
  orderable fragments; `--layout group-rows|sequential`, the synthesis
  constraint flags and the relaxation ladder.
- `tests/test_synthesis.py` (regression pin on the rejected sequence),
  `tests/test_layout_lines.py`.

Modified:

- `prosy/core/constraints.py` — `gc_windows`, `max_homopolymer_by_base`,
  `unique_kmer_size`, `soft_unique_kmer_size`/`boost`, `min_codon_frequency`,
  `gc_bands()`.
- `prosy/core/codon.py` — DnaChiselBackend honours all of the above plus
  `codon_method`; `UniquifyAllKmers` / `AvoidRareCodons` wired in, and the soft
  form added as an objective. Signature extended on the Protocol and both
  backends (the fallback accepts and does not enforce — documented).
- `prosy/core/optimize.py` — threads the new fields and `codon_method`.
- `prosy/core/layout.py` — `layout_lines()` for variable-size one-group-per-line
  packing. (Fixed an inverted rows/cols capacity check found by its tests.)
- `prosy/core/platemap.py` — PNG now rasterises the canonical SVG instead of the
  matplotlib reimplementation (which distorted badly with long legend labels);
  legend wraps instead of running off the edge.
- `prosy/core/__init__.py` — export `synthesis`.

## Open items / next steps

- ~~**Nothing is committed.**~~ Committed as `33f9be1` (+ journal `5ff2a57`).
- ~~**Retro-fit the synthesis gate onto the other pipelines.**~~ Done — see
  `2026-09-21b_synthesis-profiles-across-all-scripts.md`. The nanobody plate
  turned out to have 64/96 wells a vendor would have rejected.
- ~~**Make the v3 settings the default.**~~ Done: `--synthesis-profile
  vendor-standard` is the default in all four entry points, and reproduces the
  v3 order with no tuning flags.
- **Fold `nanobody_common.py` onto `library.py`** (carried over; still true).
- Batch runner still missing; a panel of parents is a shell loop.
- `library.write_library()` still emits no plate maps.
- Vendor rules not yet modelled: secondary structure / hairpins, terminal GC,
  and whatever else sits behind the aggregate "complexity score".
