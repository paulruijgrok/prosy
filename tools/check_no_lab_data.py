#!/usr/bin/env python3
"""Fail if real lab sequences or project identifiers appear in the repo.

This repo is public. Real binder sequences, real designs, the lab spreadsheet
name and internal project/target codenames must not be committed — the docs and
tests use synthetic stand-ins from ``data/examples/`` instead.

**The denylist is stored as hashes, not as the strings themselves.** A guard
that spelled out the confidential identifiers would publish exactly what it is
meant to keep out; that was the first version of this file, and it leaked.
Instead each secret is recorded as ``(length, truncated SHA-256)``, and
scanning slides a window of each known length over the text. Nothing here
reveals what is being looked for, and the file can safely scan itself.

Adding a new secret::

    python -c "import hashlib; print(hashlib.sha256(b'SECRET').hexdigest()[:16])"

then drop the digest into the dict keyed by ``len('SECRET')``.

Usage::

    python tools/check_no_lab_data.py            # working tree (runs in pytest)
    python tools/check_no_lab_data.py --history  # every blob in git history

``--history`` is the one that matters after a rewrite: sanitising the working
tree does nothing about what is already pushed.
"""

from __future__ import annotations

import argparse
import hashlib
import subprocess
from pathlib import Path

#: length -> truncated SHA-256 of distinctive slices of real CDRs / designs.
SEQUENCE_HASHES: dict[int, set[str]] = {
    7: {"5bee85a680029d93", "7799a67252da1df7", "7d1f78535edbd56b"},
    8: {"1d662f6c0a7d4bfd", "2dfdd42c8917ee8f", "37646ff561046c5b",
        "69a64b3fb01bc85b", "7457e617d1ad91cf"},
    13: {"6bd55e830beed9da"},
    15: {"f0e19f8e7ee2f758"},
    16: {"6e6b8d3e5a7fd4f6", "fb5ee1c281cea460"},
    19: {"b4d4c74e6eb25561"},
    34: {"2dbb348e0f3a8e13", "8f5e31b3b595b19b"},
}

#: length -> truncated SHA-256 of internal project / target / run identifiers.
#: (Generic words like the name of an ignored folder are deliberately absent:
#: they carry no information and must remain usable in .gitignore.)
IDENTIFIER_HASHES: dict[int, set[str]] = {
    4: {"0dd3f837957da1f5", "b53fbe03936d0109", "d1b513c09422db88"},
    5: {"accfdd41a090a47d", "b871a3d59fad6362"},
    7: {"a4ea198c32796984"},
    11: {"467595504f1d3ae9"},
}

SKIP_SUFFIXES = (".png", ".svg", ".xlsx", ".xls", ".pdf", ".ico", ".zip")


def _digest(text: str) -> str:
    return hashlib.sha256(text.upper().encode()).hexdigest()[:16]


def _scan_text(text: str) -> set[str]:
    """Labels of any secret found in ``text`` (generic, so the report is safe)."""
    found: set[str] = set()
    # Sequences may be line-wrapped, so also scan with whitespace removed.
    variants = [text.upper(), "".join(text.split()).upper()]
    for label, table in (("sequence", SEQUENCE_HASHES), ("identifier", IDENTIFIER_HASHES)):
        for length, digests in table.items():
            for body in variants:
                for i in range(len(body) - length + 1):
                    if _digest(body[i : i + length]) in digests:
                        found.add(f"{label} (len {length})")
                        break
                if f"{label} (len {length})" in found:
                    break
    return found


def _git(*args: str) -> str:
    return subprocess.run(["git", *args], capture_output=True, text=True,
                          check=True).stdout


def scan_worktree() -> list[tuple[str, str]]:
    """Tracked *and* untracked-but-not-ignored files.

    Tracked-only was the earlier bug: a new file is invisible to the guard
    until it is committed, which is one commit too late.
    """
    names = _git("ls-files", "--cached", "--others", "--exclude-standard").split("\n")
    hits: list[tuple[str, str]] = []
    for name in names:
        if not name or name.endswith(SKIP_SUFFIXES):
            continue
        try:
            text = Path(name).read_text(encoding="utf-8")
        except (OSError, UnicodeDecodeError):
            continue
        hits += [(name, label) for label in _scan_text(text)]
    return hits


def scan_history() -> list[tuple[str, str]]:
    """Every blob ever committed. Slow, but this is the question that counts."""
    hits: list[tuple[str, str]] = []
    seen: set[str] = set()
    for rev in _git("rev-list", "--all").split():
        for line in _git("ls-tree", "-r", rev).split("\n"):
            if not line:
                continue
            meta, _, name = line.partition("\t")
            if name.endswith(SKIP_SUFFIXES):
                continue
            blob = meta.split()[2]
            if blob in seen:
                continue
            seen.add(blob)
            content = subprocess.run(["git", "cat-file", "-p", blob],
                                     capture_output=True, text=True).stdout
            hits += [(f"{rev[:8]}:{name}", label) for label in _scan_text(content)]
    return hits


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--history", action="store_true",
                    help="Scan every blob in git history, not just the worktree.")
    args = ap.parse_args()

    hits = scan_history() if args.history else scan_worktree()
    scope = "git history" if args.history else "the working tree"
    if hits:
        print(f"FAIL: real lab data found in {scope}:")
        for where, label in sorted(set(hits)):
            print(f"  {where}: {label}")
        if args.history:
            print("\nHistory has not been rewritten yet - see docs/anonymisation.md")
        else:
            print("\nUse the synthetic sequences in data/examples/ instead.")
        return 1
    print(f"OK: no real sequences or project identifiers in {scope}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
