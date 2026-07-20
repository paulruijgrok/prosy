# Add sliding-window GC cap to nanobody DNA fragment design

## Summary

The 260626 nanobody plate fragments had a 5' GC hot spot: a 50 bp sliding-window
scan of the final (flanked) sequences peaked at **84% GC** near the flank/coding
junction — a synthesis risk. The root cause was that the codon optimizer applied
**no GC constraint at all** to the coding region (only the rarely-used padding was
bounded 0.30–0.70).

This PR adds a windowed GC cap to the codon-optimization path and, crucially,
enforces it **across the flank/coding junction** rather than only within the CDS,
by optimizing the coding region embedded in its fixed flank context.

## Changes

- **`prosy/core/codon.py`** — `DnaChiselBackend.optimize` now embeds the CDS in
  its fixed left/right flank context (flanks frozen via `AvoidChanges`;
  `AvoidPattern`/`EnforceTranslation`/`CodonOptimize` scoped to the coding
  location) and applies `EnforceGCContent(window=…)`; it returns just the CDS
  slice. `gc_window`/`left_context`/`right_context` added to the `CodonBackend`
  Protocol and the highest-frequency fallback (fallback accepts but does not
  enforce windowing). Removed a pre-existing unused import.
- **`prosy/core/constraints.py`** — new `gc_window` field on `ConstraintSet`.
- **`prosy/core/optimize.py`** — `optimize_cds` and `build_fragment` thread
  `gc_window` + flank context; `build_fragment` passes the flanks as fixed
  context into `optimize_cds`.
- **`scripts/nanobodies/nanobody_common.py`** — `RunConfig` gains
  `gc_window`/`gc_max`/`gc_min`; `optimize_all` builds the windowed constraint;
  `verify()` adds a `max_window_gc()` self-check that fails the run if any final
  fragment breaches the cap.
- **`scripts/nanobodies/make_nanobody_plate.py`** — new `--gc-window` (50),
  `--gc-max` (75), `--gc-min` (0) flags; run summary prints the active cap.
- **`scripts/nanobodies/gc_sliding_window.py`** — new standalone tool:
  sliding-window GC over a mapping CSV → long-format CSV + overlay plot (mean
  track, out-of-bounds flagging).
- Docs (`README.md`, `scripts/nanobodies/README.md`), a new dnachisel-guarded
  windowed-GC test, and a session journal.

## Verification

- Full test suite passes (`pytest`, **20 passed**), including the new
  `test_windowed_gc_cap_with_flank_context`.
- Regenerated the order at `--gc-max 75` and independently rechecked with
  `gc_sliding_window.py`: worst 50 bp window **84% → 74%** across all 96
  variants; translations preserved, flanks intact, coding regions BsaI-clean,
  all pipeline self-checks pass. Visually confirmed the 5' peak is flattened in
  the profile plot.

## Future work

- The highest-frequency **fallback backend does not enforce windowed GC** (it
  accepts the params but ignores them); `verify()`'s windowed check would fail a
  fallback-backed run. Fine while production uses DNAChisel.
- **Padding** is only globally bounded (0.30–0.70), not covered by the windowed
  optimization; `verify()` still checks the whole final fragment, so a
  pad/flank-junction breach would surface loudly.
