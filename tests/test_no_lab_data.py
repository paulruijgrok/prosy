"""The repo is public: real lab sequences and identifiers must stay out of it.

Guards the anonymisation done in 2026-09. `tools/check_no_lab_data.py --history`
is the companion check for what is already pushed, which this cannot see.
"""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path

_ROOT = Path(__file__).resolve().parents[1]


def test_no_real_lab_data_in_tracked_files():
    result = subprocess.run(
        [sys.executable, str(_ROOT / "tools" / "check_no_lab_data.py")],
        capture_output=True, text=True, cwd=_ROOT, timeout=300)
    assert result.returncode == 0, result.stdout + result.stderr


def test_example_data_is_labelled_synthetic():
    for name in ("nb01.fa", "designed_enzyme.fa"):
        text = (_ROOT / "data" / "examples" / name).read_text(encoding="utf-8")
        assert "synthetic" in text.lower(), f"{name} must say it is synthetic"


def test_run_outputs_are_gitignored():
    ignored = (_ROOT / ".gitignore").read_text(encoding="utf-8")
    assert "data/runs/" in ignored, "real run inputs/outputs must not be committable"
