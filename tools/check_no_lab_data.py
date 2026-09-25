#!/usr/bin/env python3
"""Fail if real lab sequences or project identifiers appear in tracked files.

This repo is public. Real binder sequences, real designs, the lab spreadsheet
name and internal project/target codenames must not be committed — the docs and
tests use synthetic stand-ins from ``data/examples/`` instead.

The patterns below are deliberately *fingerprints*, not the full sequences: a
denylist that itself contained the confidential data would defeat the purpose.
Each entry is a short, distinctive slice that is present in the real material
and absent from anything synthetic.

Usage::

    python tools/check_no_lab_data.py            # tracked files (pre-commit)
    python tools/check_no_lab_data.py --history  # every blob in git history

``--history`` is the one that matters after a history rewrite: sanitising the
working tree does nothing about what is already pushed.
"""

from __future__ import annotations

import argparse
import subprocess
from pathlib import Path

#: Distinctive fragments of real CDRs and design sequences.
SEQUENCE_FINGERPRINTS = [
    "REDACTED", "REDACTE", "REDACTED_SEQUEN",          # real VHH CDRs
    "REDACTED", "REDACTED", "REDACTED_SEQU",
    "REDACTED", "REDACTE", "REDACTED_SEQUENC",
    "REDACTED", "REDACTE", "REDACTED_SEQUENC",
    "REDACTED_SEQUENCE",                              # real framework+CDR join
    "REDACTED_SEQUENCE",               # real design N-termini
    "REDACTED_SEQUENCE",
]

#: Internal project, target and run identifiers.
IDENTIFIERS = [
    "REDACTED", "SetA", "SetB", "designgroup", "designs", "target1", "TargetA",
    "Working folder",
]

SKIP_SUFFIXES = (".png", ".svg", ".xlsx", ".xls", ".pdf")
#: The denylist names the identifiers, so it would always flag itself.
SELF = "tools/check_no_lab_data.py"


def _patterns() -> list[str]:
    return SEQUENCE_FINGERPRINTS + IDENTIFIERS


def scan_tracked() -> list[tuple[str, str]]:
    files = subprocess.run(["git", "ls-files"], capture_output=True, text=True,
                           check=True).stdout.split()
    hits = []
    for name in files:
        if name in (SELF, ".gitignore") or name.endswith(SKIP_SUFFIXES):
            continue
        try:
            text = Path(name).read_text(encoding="utf-8")
        except (OSError, UnicodeDecodeError):
            continue
        hits += [(name, p) for p in _patterns() if p in text]
    return hits


def scan_history() -> list[tuple[str, str]]:
    """Every blob ever committed. Slow, but this is the question that counts."""
    revs = subprocess.run(["git", "rev-list", "--all"], capture_output=True,
                          text=True, check=True).stdout.split()
    hits: list[tuple[str, str]] = []
    seen: set[str] = set()
    for rev in revs:
        listing = subprocess.run(["git", "ls-tree", "-r", "--name-only", rev],
                                 capture_output=True, text=True, check=True).stdout.split("\n")
        for name in listing:
            if not name or name.endswith(SKIP_SUFFIXES) or name in (SELF, ".gitignore"):
                continue
            key = f"{rev}:{name}"
            blob = subprocess.run(["git", "rev-parse", key], capture_output=True,
                                  text=True).stdout.strip()
            if not blob or blob in seen:
                continue
            seen.add(blob)
            content = subprocess.run(["git", "cat-file", "-p", blob],
                                     capture_output=True, text=True).stdout
            hits += [(f"{rev[:8]}:{name}", p) for p in _patterns() if p in content]
    return hits


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--history", action="store_true",
                    help="Scan every blob in git history, not just the worktree.")
    args = ap.parse_args()

    hits = scan_history() if args.history else scan_tracked()
    scope = "git history" if args.history else "tracked files"
    if hits:
        print(f"FAIL: real lab data found in {scope}:")
        for where, pattern in sorted(set(hits)):
            print(f"  {where}: {pattern!r}")
        if args.history:
            print("\nA history rewrite has not (fully) happened yet - see "
                  "docs/anonymisation.md")
        return 1
    print(f"OK: no real sequences or project identifiers in {scope}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
