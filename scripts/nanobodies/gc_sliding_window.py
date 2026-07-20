#!/usr/bin/env python3
"""Sliding-window GC content for designed nanobody DNA fragments.

Reads a nanobody mapping CSV (as written by make_nanobody_plate.py) and, for
each variant, slides a fixed-width window along the chosen DNA column and
reports the GC% in each window. Useful for spotting local GC extremes,
especially near the 5' end.

Example:
    conda activate DataAnalysis
    python gc_sliding_window.py \
        --input "../../Working folder/260626_NanobodyMuts/260626_nanobody_mapping.csv" \
        --window 50 --step 1 --out gc_profiles.csv
"""
import argparse
import csv
import sys


def gc_percent(seq: str) -> float:
    """GC% of a sequence (case-insensitive). Empty -> 0.0."""
    if not seq:
        return 0.0
    return 100.0 * sum(c in "GCgc" for c in seq) / len(seq)


def sliding_gc(seq: str, window: int, step: int):
    """Yield (window_start, gc_percent) for each window along seq.

    window_start is 0-based, inclusive. Sequences shorter than `window`
    yield a single window covering the whole sequence.
    """
    n = len(seq)
    if n <= window:
        yield 0, gc_percent(seq)
        return
    for start in range(0, n - window + 1, step):
        yield start, gc_percent(seq[start:start + window])


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--input", required=True,
                   help="Nanobody mapping CSV (needs an ID column and the sequence column).")
    p.add_argument("--column", default="Final_DNA",
                   help="DNA column to scan (default: Final_DNA; use Coding_DNA to skip flanks).")
    p.add_argument("--id-column", default="NbID", help="Variant ID column (default: NbID).")
    p.add_argument("--window", type=int, default=50, help="Window width in bp (default: 50).")
    p.add_argument("--step", type=int, default=1, help="Step between windows in bp (default: 1).")
    p.add_argument("--out", default="gc_profiles.csv",
                   help="Long-format output CSV: id, well, window_start, window_end, gc (default: gc_profiles.csv).")
    p.add_argument("--low", type=float, default=25.0,
                   help="Flag any window below this GC%% in the summary (default: 25).")
    p.add_argument("--high", type=float, default=75.0,
                   help="Flag any window above this GC%% in the summary (default: 75).")
    p.add_argument("--plot", nargs="?", const="gc_profiles.png", default=None,
                   metavar="PATH",
                   help="Also save a GC-track plot (overlay of all variants + mean). "
                        "Optional path (default: gc_profiles.png).")
    p.add_argument("--flank-5-len", type=int, default=17,
                   help="Length of the 5' flank, marked on the plot (default: 17).")
    args = p.parse_args(argv)

    if args.window < 1 or args.step < 1:
        p.error("--window and --step must be >= 1")

    with open(args.input, newline="") as fh:
        rows = list(csv.DictReader(fh))
    if not rows:
        p.error(f"No data rows in {args.input}")
    for col in (args.id_column, args.column):
        if col not in rows[0]:
            p.error(f"Column {col!r} not found. Available: {', '.join(rows[0])}")

    has_well = "Well" in rows[0]
    flagged = []
    written = 0
    profiles = []  # (id, well, starts, vals) for optional plotting
    with open(args.out, "w", newline="") as out:
        w = csv.writer(out)
        w.writerow(["id", "well", "window_start", "window_end", "gc"])
        for r in rows:
            seq = r[args.column] or ""
            starts, vals = [], []
            for start, gc in sliding_gc(seq, args.window, args.step):
                end = min(start + args.window, len(seq))
                w.writerow([r[args.id_column], r.get("Well", "") if has_well else "",
                            start, end, f"{gc:.2f}"])
                starts.append(start)
                vals.append(gc)
                written += 1
            profiles.append((r[args.id_column], r.get("Well", ""), starts, vals))
            if vals:
                lo, hi = min(vals), max(vals)
                if lo < args.low or hi > args.high:
                    flagged.append((r[args.id_column], r.get("Well", ""), lo, hi))

    print(f"Wrote {written} windows for {len(rows)} variants -> {args.out}")
    print(f"Window={args.window} bp, step={args.step} bp, column={args.column}")
    if flagged:
        print(f"\n{len(flagged)} variant(s) with a window outside "
              f"[{args.low:.0f}, {args.high:.0f}]% GC:")
        for vid, well, lo, hi in flagged:
            print(f"  {vid:10} {well:4}  min={lo:5.1f}  max={hi:5.1f}")
    else:
        print(f"\nAll variants stay within [{args.low:.0f}, {args.high:.0f}]% GC.")

    if args.plot is not None:
        plot_profiles(profiles, args)
    return 0


def plot_profiles(profiles, args):
    """Overlay every variant's GC track (thin/translucent) plus the mean track.

    x-axis is window-start position (5' end at left); the 5' flank region and
    the low/high GC bounds are marked for reference.
    """
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print("\n[plot skipped] matplotlib not installed. "
              "Install with: pip install matplotlib")
        return

    profiles = [p for p in profiles if p[3]]
    if not profiles:
        print("\n[plot skipped] no data to plot.")
        return

    fig, ax = plt.subplots(figsize=(11, 5))
    for _id, _well, starts, vals in profiles:
        ax.plot(starts, vals, color="#4C72B0", alpha=0.15, linewidth=0.7)

    # Mean track over the common length (windows shared by all variants).
    min_len = min(len(v) for _, _, _, v in profiles)
    if min_len:
        xs = profiles[0][2][:min_len]
        means = [sum(v[i] for _, _, _, v in profiles) / len(profiles)
                 for i in range(min_len)]
        ax.plot(xs, means, color="#C44E52", linewidth=2.2, label="mean")

    ax.axhline(args.low, color="grey", linestyle=":", linewidth=1)
    ax.axhline(args.high, color="grey", linestyle=":", linewidth=1)
    if args.flank_5_len > 0:
        ax.axvspan(0, args.flank_5_len, color="#DD8452", alpha=0.12,
                   label=f"5' flank ({args.flank_5_len} bp)")

    ax.set_xlabel(f"Window start position (bp) — {args.window} bp window, 5' at left")
    ax.set_ylabel("GC content (%)")
    ax.set_title(f"Sliding-window GC across {len(profiles)} variants "
                 f"({args.column})")
    ax.set_ylim(0, 100)
    ax.margins(x=0)
    ax.legend(loc="upper right", fontsize=9)
    fig.tight_layout()
    fig.savefig(args.plot, dpi=150)
    print(f"\nSaved plot -> {args.plot}")


if __name__ == "__main__":
    sys.exit(main())
