#!/usr/bin/env python3
"""Emit a git-filter-repo --replace-text list for purging real lab data.

The denylist in ``check_no_lab_data.py`` is stored as hashes, so it cannot tell
you *what* to redact. This script closes that loop: it slides the same hash
windows over every blob in history and, where one matches, recovers the exact
substring from that blob. The result is a replacement list containing precisely
the secrets that are actually still present - no plaintext denylist anywhere,
no guessing at case variants (they are captured as found), and no
over-redaction of synthetic sequences that merely look similar.

An earlier version harvested every long biological string instead. That
over-matched badly: it redacted the synthetic test sequences too and broke 37
tests. Matching known hashes is exact.

**The output contains the confidential strings.** Write it outside the repo,
use it, delete it::

    python tools/build_redactions.py > /tmp/redactions.txt
    git filter-repo --replace-text /tmp/redactions.txt --force
    python tools/check_no_lab_data.py --history        # must print OK
    rm /tmp/redactions.txt

See docs/anonymisation.md for the full procedure.
"""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from check_no_lab_data import (  # noqa: E402
    IDENTIFIER_HASHES,
    SEQUENCE_HASHES,
    SKIP_SUFFIXES,
    _digest,
)


def historical_blobs():
    """Every distinct blob in history, as ``(name, text)``."""
    revs = subprocess.run(["git", "rev-list", "--all"], capture_output=True,
                          text=True, check=True).stdout.split()
    seen: set[str] = set()
    for rev in revs:
        listing = subprocess.run(["git", "ls-tree", "-r", rev], capture_output=True,
                                 text=True, check=True).stdout.split("\n")
        for line in listing:
            if not line:
                continue
            meta, _, name = line.partition("\t")
            if name.endswith(SKIP_SUFFIXES):
                continue
            blob = meta.split()[2]
            if blob in seen:
                continue
            seen.add(blob)
            yield name, subprocess.run(["git", "cat-file", "-p", blob],
                                       capture_output=True, text=True).stdout


def _stripped_positions(text: str) -> tuple[str, list[int]]:
    """Return (whitespace-stripped uppercase text, raw-position for each char).

    Used to find secrets that span line-breaks — e.g. a FASTA sequence wrapped
    at 60 chars per line, or a Python string literal split across two source
    lines.  The returned positions map each character in the stripped string
    back to its index in the original raw text, so the corresponding raw span
    (including any intervening whitespace / punctuation) can be extracted.
    """
    chars: list[str] = []
    positions: list[int] = []
    for i, ch in enumerate(text):
        if not ch.isspace():
            chars.append(ch.upper())
            positions.append(i)
    return "".join(chars), positions


def find_secrets() -> dict[str, str]:
    """``{exact substring found: replacement}`` for everything still in history.

    Scans both the raw text (finds secrets on a single line) and the
    whitespace-stripped text (finds secrets that span line-breaks, e.g. a FASTA
    wrapped at 60 chars per line, or a Python string literal split across two
    source lines).  The raw span — including any intervening whitespace /
    punctuation — is used as the filter-repo literal key so the replacement is
    exact.
    """
    found: dict[str, str] = {}
    tables = ((SEQUENCE_HASHES, "REDACTED_SEQUENCE"), (IDENTIFIER_HASHES, "REDACTED"))
    for _name, text in historical_blobs():
        stripped, positions = _stripped_positions(text)
        for table, replacement in tables:
            for length, digests in table.items():
                # --- raw scan (single-line occurrences) -----------------------
                for i in range(len(text) - length + 1):
                    window = text[i : i + length]
                    if _digest(window) in digests:
                        found.setdefault(window, replacement[:length] or replacement)

                # --- stripped scan (cross-whitespace occurrences) -------------
                for i in range(len(stripped) - length + 1):
                    if _digest(stripped[i : i + length]) in digests:
                        raw_start = positions[i]
                        raw_end = positions[i + length - 1] + 1
                        raw_span = text[raw_start:raw_end]
                        if raw_span not in found:
                            span_len = len(raw_span)
                            repl = (replacement * (span_len // len(replacement) + 1))[
                                :span_len
                            ]
                            found[raw_span] = repl
    return found


def main() -> int:
    secrets = find_secrets()
    if not secrets:
        print("# Nothing left to redact - history is already clean.", file=sys.stderr)
        return 0
    for secret, replacement in sorted(secrets.items(), key=lambda kv: -len(kv[0])):
        print(f"literal:{secret}==>{replacement}")
    print(f"# {len(secrets)} string(s) found in history. CONTAINS SECRETS - "
          "delete this file after use.", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
