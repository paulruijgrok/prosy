"""Every fragment-producing entry point must expose the same synthesis controls.

The vendor rejection happened because one pipeline's default differed from
another's. These tests are a guard against that drifting apart again: if a new
script builds DNA, it gets the profile flag and the gate too.
"""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path

import pytest

_ROOT = Path(__file__).resolve().parents[1]

#: Scripts whose output can end up in a synthesis order.
FRAGMENT_SCRIPTS = [
    "scripts/scan/mutational_scan.py",
    "scripts/nanobodies/nanobody_scan.py",
    "scripts/nanobodies/make_nanobody_plate.py",
    "scripts/designs/make_design_plate.py",
]


def _help(script: str) -> str:
    result = subprocess.run(
        [sys.executable, str(_ROOT / script), "--help"],
        capture_output=True, text=True, timeout=120, cwd=_ROOT,
    )
    assert result.returncode == 0, result.stderr
    return result.stdout


@pytest.mark.parametrize("script", FRAGMENT_SCRIPTS)
def test_exposes_the_synthesis_profile(script):
    text = _help(script)
    assert "--synthesis-profile" in text
    assert "vendor-standard" in text


@pytest.mark.parametrize("script", FRAGMENT_SCRIPTS)
def test_gate_is_on_by_default_and_can_be_turned_off(script):
    # An opt-*out* flag is the point: a run that forgets to ask for checking
    # is still checked.
    assert "--no-check-synthesis" in _help(script)


if __name__ == "__main__":
    raise SystemExit(pytest.main([__file__, "-q"]))
