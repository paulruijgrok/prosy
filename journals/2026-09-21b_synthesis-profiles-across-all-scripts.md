# 2026-09-21 (second session) — One synthesis profile across every script

Follow-up to `2026-09-21_design-plates-synthesis-manufacturability.md`, which
closed with: *"the bug is in the shared default, not in the design script."*
This session acts on that.

## What changed and why

The vendor rejection was never really about the designs designs. It was about
`use_best_codon` being the default everywhere, so **every** pipeline in the repo
was producing the same repeat-heavy sequences — `make_design_plate.py` had been
fixed in place, while `make_nanobody_plate.py` and both scan CLIs were still
building orders the old way.

Measured on the 96-well nanobody plate (same parents, same mutations), the old
default would have had **64 of 96 wells rejected**:

| | legacy default | vendor-standard |
|---|---|---|
| 8-mer repeat coverage | 16.7–47.4% | **0.0–8.3%** |
| densest 90 bp window | 70.0% | **18.9%** |
| overall GC | 56.2–61.4% | **54.5–56.0%** |
| rare codons | 0.00% | 0.02% |
| wells a vendor would reject | **64 / 96** | **0 / 96** |

So the nanobody plates carried the same latent defect, just less extreme —
repeat coverage scales with length, and a 121 aa nanobody is a third of a 333 aa
design.

## Design decisions worth remembering

- **A named profile, not a long command line.** `SynthesisProfile` in
  `prosy.core.synthesis` bundles the design-time constraints *and* the
  acceptance spec they target. `vendor-standard` is the default in all four
  entry points; `none` reproduces the legacy behaviour. The settings that got
  an order accepted are now what you get by doing nothing, which is the actual
  fix for "this recurs by omission".
- **The gate is opt-out.** `--no-check-synthesis`, never `--check-synthesis`. A
  run that forgets to ask for checking is still checked. `--check-synthesis`
  survives on `make_design_plate.py` as a deprecated no-op so the recorded v3
  command line keeps working.
- **The profile pairs constraints with the spec they aim at.** The GC design cap
  (0.56) sits deliberately below the acceptance limit (0.58): a solver satisfies
  a cap by sitting against it, so aiming at the limit leaves no margin. Keeping
  both numbers in one object is what stops them drifting apart.
- **Explicit flags replace, never stack.** An explicit `--gc-window 50` removes
  the profile's 50 bp band rather than adding a second one. Two overlapping
  `EnforceGCContent` constraints are mathematically redundant but change
  DNAChisel's search path — this showed up as 8 relaxed designs instead of 6 and
  a different (still valid) sequence set. Fixed in all three call sites.
- **Shared ladder.** `build_fragment_with_ladder()` in `optimize.py` is now the
  single implementation of the "try strictest, degrade to soft objective" retry;
  `make_design_plate.py` had its own copy and no longer does.

## The reproducibility bug found along the way

Checking that the refactor had not changed output revealed that **the pipeline
was never reproducible**: two identical invocations produced 84/84 different
sequences. DNAChisel draws on the *global* `random` (and numpy) state, and
`DnaChiselBackend` never seeded it, so each sequence depended on how much RNG
earlier sequences had consumed.

`DnaChiselBackend` now seeds both generators per call from ``seed`` plus a CRC32
of that call's own inputs, and restores the caller's state afterwards. Two
properties follow, both pinned by tests:

1. runs are reproducible;
2. a sequence does not change when other designs are added to or removed from
   the batch (because the seed derives from the input, not the call order).

**Consequence: the accepted `260921_designs_v3_*` sequences cannot be
regenerated.** They were produced by the unseeded code. The delivered mapping
CSV, FASTA and plate file are the record — do not re-run that command expecting
the same output. The same applies to every order placed before this commit.

A related discovery: a *different* seed often produces the *same* sequence.
DNAChisel converges to one answer for a well-determined problem, so the earlier
chaos came purely from state accumulating across calls, not from seed
sensitivity. The test asserts validity for any seed rather than difference.

## What the gate caught in the test suite

Two `test_library` tests began failing once the gate was on by default — the
gate working, not a regression. They use the dependency-free fallback backend,
which accepts the profile but cannot enforce it. Rather than silence them, the
real behaviour is now pinned: the fallback leaves >15% repeat coverage despite
`unique_kmer_size=8` being requested, and whether that breaches a vendor limit
depends on length. A 121 aa nanobody squeaks under 40%; the 333 aa designs did
not. **The gate, not the request, is what decides.**

## Exact commands

```bash
cd ~/Documents/Claude/Projects/ProSy
conda activate DataAnalysis
python -m pytest tests/ -q          # 143 passed
ruff check prosy scripts tests --select F401,F841
vulture prosy scripts --min-confidence 80

# All four entry points, defaults (profile on, gate on):
python scripts/scan/mutational_scan.py --protein-file p.fa --scan alanine \
    --fragments --destination gg002
python scripts/nanobodies/nanobody_scan.py --parent Nb01 --scan alanine --fragments
python scripts/nanobodies/make_nanobody_plate.py --data-dir <dir> --stamp <stamp>
python scripts/designs/make_design_plate.py --csv designs.csv --destination gg002

# Legacy behaviour, for comparison only:
... --synthesis-profile none
```

The v3 order is now reproducible from the profile alone — no tuning flags:

```bash
python scripts/designs/make_design_plate.py \
    --csv "Working folder/260917_designs_enzyme/sequences_bypLDDT.csv" \
    --destination "Working folder/260917_designs_enzyme/gg002.txt" \
    --layout sequential --orientation column \
    --flank-5 TGTATCGGTCTCGAGGA --flank-3 GGTTCCGGAGACCTCTAGT \
    --expect-prefix MSG --expect-suffix GSHHHHHH \
    --check-gc-min 0.40 --check-gc-max 0.58
```

(verified to produce the same sequences as the fully explicit v3 command, which
in turn matches the pre-refactor reference).

## Files touched (uncommitted)

New:

- `tests/test_synthesis_profiles.py` — profiles, ladder structure, and the four
  reproducibility properties.
- `tests/test_cli_consistency.py` — every fragment-producing script exposes
  `--synthesis-profile` and `--no-check-synthesis`. A guard against the defaults
  drifting apart again, which is what caused the rejection.

Modified:

- `prosy/core/synthesis.py` — `SynthesisProfile`, `PROFILES`
  (`vendor-standard`, `none`), `get_profile()`, `constraint_ladder()`.
- `prosy/core/codon.py` — `DnaChiselBackend(seed=...)` and the
  `_deterministic()` context manager; `get_backend` passes the seed through.
- `prosy/core/optimize.py` — `build_fragment_with_ladder()`.
- `prosy/core/library.py` — `LibraryConfig.synthesis_profile`/`check_synthesis`,
  `constraint_ladder()`, the gate in `verify_library()`, repeat/GC/ladder
  columns in the mapping.
- `scripts/nanobodies/nanobody_common.py`, `make_nanobody_plate.py`,
  `scripts/scan/scan_common.py`, `scripts/designs/make_design_plate.py` — the
  shared flags, reporting and ladder.
- `tests/test_library.py` — fallback-backend expectations made explicit.

## Open items / next steps

- Nothing committed yet.
- **README does not yet describe profiles.** It still presents `use_best_codon`
  as the codon backend story with no mention of the repeat trap. Highest-value
  docs work; a `docs/synthesis.md` page is probably warranted too.
- Fold `nanobody_common.py` onto `library.py` — now more attractive, since both
  carry near-identical ladder/gate wiring.
- Batch runner still missing; `library.write_library()` still emits no plate maps.
- Vendor rules not modelled: secondary structure/hairpins, terminal GC.
- Only one vendor profile exists. A second vendor means a second
  `SynthesisProfile`, which is the shape the code now expects.
