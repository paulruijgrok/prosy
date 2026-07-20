# 2026-07-20 — Sliding-window GC cap for nanobody DNA fragments

## What changed and why

The 260626 nanobody plate order had a GC hot spot at the 5' end: a 50 bp
sliding-window scan of the final (flanked) fragments peaked at **84% GC** near
the flank/coding junction. High local GC at the 5' end risks synthesis failures,
so we wanted to cap the local GC content.

Root cause: the coding-region optimizer applied **no GC constraint at all** to
the CDS — only the padding had GC bounds (0.30–0.70, and padding is usually
absent for full-length nanobodies anyway). So the E. coli codon optimizer was
free to produce GC-rich runs at the start of the CDS.

Fix: added a **windowed GC constraint** to the codon-optimization path, enforced
by DNAChisel's `EnforceGCContent(window=...)`. Crucially, the coding region is
now optimized **inside its fixed flank context** (flanks held immutable via
`AvoidChanges`), so the cap holds across the flank/coding junction — the actual
hot spot — not just within the CDS. Forbidden patterns (BsaI) still apply to the
coding region only, so the flanks keep their intended Type IIS sites.

Target evolved during the session: first "below 76%" (→ windows land at 76.0%),
then tightened to "below 75%". Because a 50 bp window's GC moves in 2% steps,
`--gc-max 75` resolves to **≤74%** in practice. Default is now 75%.

## Files touched (committed as `c080363` on `feature/prosy-core`)

- `prosy/core/codon.py` — `DnaChiselBackend.optimize` rewritten to embed the CDS
  in fixed left/right flank context, add `gc_window`, apply `AvoidChanges` to
  flanks and `AvoidPattern`/`CodonOptimize`/`EnforceTranslation` scoped to the
  coding location; returns the CDS slice. Interface (`gc_window`, `left_context`,
  `right_context`) added to the `CodonBackend` Protocol and the fallback
  `HighestFrequencyBackend` (fallback accepts but does not enforce windowing).
  Also removed a pre-existing unused import (`SYNONYMOUS_CODONS`).
- `prosy/core/constraints.py` — added `gc_window` field to `ConstraintSet`.
- `prosy/core/optimize.py` — `optimize_cds` and `build_fragment` thread
  `gc_window` + flank context; `build_fragment` passes the flanks as fixed
  context into `optimize_cds`.
- `scripts/nanobodies/nanobody_common.py` — `RunConfig` gains
  `gc_window`/`gc_max`/`gc_min`; `optimize_all` builds the windowed
  `ConstraintSet`; `verify()` gains a `max_window_gc()` self-check that fails the
  run if any final fragment breaches the cap.
- `scripts/nanobodies/make_nanobody_plate.py` — new CLI flags `--gc-window` (50),
  `--gc-max` (75), `--gc-min` (0); summary line prints the active cap.
- `scripts/nanobodies/gc_sliding_window.py` — **new** standalone tool: sliding
  window GC over a mapping CSV, long-format CSV output + overlay plot with mean
  track and flagging of windows outside bounds.
- `tests/test_core.py` — new `test_windowed_gc_cap_with_flank_context`
  (dnachisel-guarded) asserting every window (junction included) respects the cap.
- `README.md`, `scripts/nanobodies/README.md` — documented the cap.

## Exact commands (run in the `DataAnalysis` conda env, which has dnachisel/biopython/openpyxl)

```bash
conda activate DataAnalysis
cd "/Users/paulruijgrok/Documents/Claude/Projects/ProSy"

# Full test suite (pre-commit gate)
python -m pytest tests/ -q

# Regenerate the order at <75% GC into a NEW dir (inputs copied there first)
cd scripts/nanobodies
python make_nanobody_plate.py --backend dnachisel --gc-window 50 --gc-max 75 \
  --stamp 260626 --data-dir "../../Working folder/260720_NanobodyMuts_gc75"

# Independent verification + profile plot
python gc_sliding_window.py \
  --input "../../Working folder/260720_NanobodyMuts_gc75/260626_nanobody_mapping.csv" \
  --window 50 \
  --out   "../../Working folder/260720_NanobodyMuts_gc75/260626_gc_windows_50bp.csv" \
  --plot  "../../Working folder/260720_NanobodyMuts_gc75/260626_gc_profile_50bp.png"
```

## Results / verification

- Full suite: **20 passed** (incl. new windowed-GC test).
- Independent recompute of the regenerated order: worst 50 bp window
  **84% → 74%** across all 96 variants; translation preserved, flanks intact,
  coding regions BsaI-clean, all pipeline self-checks pass.
- Output dir: `Working folder/260720_NanobodyMuts_gc75/` (original
  `260626_NanobodyMuts/` left untouched).

## Open items / next steps

- **Committed and pushed** as `c080363` on `feature/prosy-core` (rebased onto
  origin's `12b234b` README tweak; original pre-rebase hash was `06b0334`).
  Still on the feature branch — not yet merged to `main`; open a PR when ready.
- **Fallback backend does not enforce windowed GC** — `HighestFrequencyBackend`
  accepts `gc_window`/context but ignores them. Fine while production uses
  DNAChisel, but the pipeline's `verify()` windowed check would fail a
  fallback-backed run. Consider either a best-effort windowed nudge in the
  fallback or a clearer guard.
- **Padding not covered by the windowed optimization** — the cap is enforced on
  flank+CDS; padding (rare for full-length nanobodies, since ~396 bp > 300 min)
  is only globally bounded 0.30–0.70. `verify()` does check the whole final
  fragment, so a pad/flank-junction breach would surface loudly — but worth a
  note if min-length is ever raised.
- **`gc76` → `gc75` rename**: the output dir was created as `..._gc76` then
  renamed to `..._gc75` after tightening the cap. Confirm nothing downstream
  hardcodes the old name.
- Consider a lower `--gc-min` floor experiment if any 3' AT-rich dips matter
  (currently no floor; min window sits ~38–44%).
