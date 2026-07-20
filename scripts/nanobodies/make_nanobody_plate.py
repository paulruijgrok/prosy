#!/usr/bin/env python3
"""Build a plate of synthesis-ready nanobody DNA fragments (general).

One parameterised entry point for nanobody plates. Each parent nanobody plus its
top-N recommended point mutations forms a group that fills either a plate column
or a plate row; groups pack onto a 96- or 384-well plate.

Defaults reproduce the original 260626 column plate (12 parents x 7 mutations,
column orientation). Change ``--orientation``, ``--mutations-per-parent``, the
parent lists and ``--plate-size`` for other layouts.

Examples
--------
# Original column plate (12 parents, 7 mutations each), one parent per column:
python make_nanobody_plate.py

# Row plate: 5 ADL1 + 3 KRU1 parents, 11 mutations each, one parent per row:
python make_nanobody_plate.py --orientation row \\
    --adl1 Nb01 Nb02 Nb03 Nb04 Nb05 --kru1 Nb50 Nb51 Nb58 \\
    --mutations-per-parent 11 --stamp 260626_row

# 16 parents x 5 mutations, 2 groups per row:
python make_nanobody_plate.py --orientation row --mutations-per-parent 5 \\
    --adl1 ... (8 ids) --kru1 ... (8 ids) --stamp 260626_16x5

Run with DNAChisel installed to reproduce codon_optimize.py exactly; otherwise a
dependency-free fallback keeps the pipeline runnable.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

_HERE = Path(__file__).resolve()
sys.path.insert(0, str(_HERE.parents[2]))   # repo root -> 'prosy'
sys.path.insert(0, str(_HERE.parent))       # this dir -> 'nanobody_common'

import nanobody_common as nb  # noqa: E402
from prosy.core import codon as codon_mod  # noqa: E402
from prosy.core import layout as layout_mod  # noqa: E402
from prosy.core.cloning import Flanks  # noqa: E402

# Default parent selection (reproduces the original column plate).
DEFAULT_ADL1 = ["Nb01", "Nb02", "Nb03", "Nb04", "Nb05", "Nb06", "Nb07"]
DEFAULT_KRU1 = ["Nb50", "Nb51", "Nb52", "Nb53", "Nb58"]

DEFAULT_FLANK_5 = "TGTATCGGTCTCgAGGA"
DEFAULT_FLANK_3 = "GGTTCCgGAGACCTCTAGT"
XLSX_NAME = "Malaria DX plasmids.xlsx"


def main() -> None:
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    p.add_argument("--data-dir", type=Path,
                   default=_HERE.parents[2] / "Working folder" / "260626_NanobodyMuts",
                   help="Folder with the xlsx + TSV inputs (outputs written here).")
    p.add_argument("--stamp", default="260626", help="Output filename prefix.")
    # Parent selection
    p.add_argument("--adl1", nargs="*", default=DEFAULT_ADL1, metavar="NbID",
                   help="ADL1 parent IDs (in plate order).")
    p.add_argument("--kru1", nargs="*", default=DEFAULT_KRU1, metavar="NbID",
                   help="KRU1 parent IDs (in plate order). Appended after ADL1.")
    # Layout
    p.add_argument("--orientation", choices=layout_mod.ORIENTATIONS, default="column",
                   help="'column' = one group per column; 'row' = one group per row.")
    p.add_argument("--plate-size", type=int, choices=[96, 384], default=96)
    p.add_argument("--plate-map", choices=["none", "svg", "png", "both"], default="svg",
                   help="Plate-map image export (SVG needs no deps; PNG uses matplotlib/LibreOffice).")
    p.add_argument("--mutations-per-parent", type=int, default=7,
                   help="Top-N mutations per parent by score (group = parent + N).")
    p.add_argument("--allow-straddle", action="store_true",
                   help="Allow groups to straddle row/column boundaries.")
    # Synthesis knobs
    p.add_argument("--flank-5", default=DEFAULT_FLANK_5)
    p.add_argument("--flank-3", default=DEFAULT_FLANK_3)
    p.add_argument("--species", default="e_coli")
    p.add_argument("--avoid-enzymes", nargs="*", default=["BsaI"], metavar="ENZYME")
    p.add_argument("--min-length", type=int, default=300)
    p.add_argument("--gc-window", type=int, default=50, metavar="BP",
                   help="Sliding-window width (bp) for the local GC cap; 0 disables (default: 50).")
    p.add_argument("--gc-max", type=float, default=75.0, metavar="PCT",
                   help="Max GC%% allowed in any --gc-window (default: 75; "
                        "a 50 bp window resolves this to <=74%% in practice).")
    p.add_argument("--gc-min", type=float, default=0.0, metavar="PCT",
                   help="Min GC%% required in any --gc-window (default: 0 = no floor).")
    p.add_argument("--new-id-start", type=int, default=73)
    p.add_argument("--new-id-prefix", default="Nb")
    p.add_argument("--backend", choices=["auto", "dnachisel", "highest_frequency"],
                   default="auto")
    p.add_argument("--seed", type=int, default=1)
    args = p.parse_args()

    parent_order = list(args.adl1) + list(args.kru1)
    cfg = nb.RunConfig(
        flanks=Flanks(five_prime=args.flank_5, three_prime=args.flank_3),
        species=args.species, avoid_enzymes=args.avoid_enzymes,
        min_length=args.min_length,
        gc_window=args.gc_window if args.gc_window > 0 else None,
        gc_max=args.gc_max / 100.0, gc_min=args.gc_min / 100.0,
        mutations_per_parent=args.mutations_per_parent,
        new_id_start=args.new_id_start, new_id_prefix=args.new_id_prefix,
        orientation=args.orientation, plate_size=args.plate_size,
        plate_map=args.plate_map,
    )

    # Fail early with a clear layout description if it doesn't fit.
    try:
        layout_mod.layout_groups(
            n_groups=len(parent_order), group_size=cfg.group_size,
            plate=cfg.plate_size, orientation=cfg.orientation,
            require_aligned=not args.allow_straddle,
        )
    except layout_mod.LayoutError as e:
        raise SystemExit(f"Layout error: {e}")

    prefer = None if args.backend == "auto" else args.backend
    backend = codon_mod.get_backend(prefer, seed=args.seed)

    print(f"Codon backend: {backend.name}")
    print(layout_mod.describe_layout(
        len(parent_order), cfg.group_size, plate=cfg.plate_size, orientation=cfg.orientation))
    print(f"Parents ({len(parent_order)}): {', '.join(parent_order)}")
    print(f"Species: {cfg.species}  Avoid: {cfg.avoid_enzymes or '(none)'}  "
          f"Min length: {cfg.min_length}")
    if cfg.gc_window:
        floor = f", min {cfg.gc_min*100:.0f}%" if cfg.gc_min > 0 else ""
        print(f"GC cap: max {cfg.gc_max*100:.0f}%{floor} over {cfg.gc_window} bp windows")
    else:
        print("GC cap: (none)")
    print(f"Flanks: 5'={cfg.flanks.five_prime}  3'={cfg.flanks.three_prime}")
    print(f"Mutations/parent: {cfg.mutations_per_parent}  "
          f"New IDs from: {cfg.new_id_prefix}{cfg.new_id_start}")

    variants, problems, out = nb.run_pipeline(
        data_dir=args.data_dir, xlsx_name=XLSX_NAME, parent_order=parent_order,
        tsv_by_set=nb.DEFAULT_TSV_BY_SET, stamp=args.stamp, backend=backend,
        cfg=cfg, seed=args.seed,
    )

    print(f"\nVariants: {len(variants)} "
          f"(mutants: {sum(not v.is_parent for v in variants)})")
    print(f"Plate xlsx : {out['plate_xlsx']}")
    print(f"Plate xls  : {out['plate_xls'] or '(LibreOffice not available - xlsx only)'}")
    print(f"Design CSV : {out['design']}")
    print(f"Mapping CSV: {out['mapping']}")
    if out.get("map_svg"):
        print(f"Plate map  : {out['map_svg']}")
    if out.get("map_png"):
        print(f"Plate map  : {out['map_png']}")

    if problems:
        print("\nVERIFICATION FAILED:")
        for pr in problems:
            print("  -", pr)
        sys.exit(1)
    print("\nAll verification checks passed.")


if __name__ == "__main__":
    main()
