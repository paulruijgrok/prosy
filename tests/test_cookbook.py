"""Keep the commands in docs/cookbook.md honest.

Documented commands rot silently: a flag gets renamed and the docs keep saying
the old name until someone pastes it and gets an argparse error. That already
happened once — ``--expect-prefix`` was documented for the scan CLIs before it
existed there.

Two tiers:

* **Lint (always runs, ~1 s).** Parse every fenced block, check the scripts
  exist, check every ``--flag`` is actually accepted by the script it is passed
  to, and check the Python blocks compile. This catches renamed and invented
  flags without running any optimization.
* **Smoke (opt-in, minutes).** Actually execute every bash block, in a sandbox
  copy of the repo so nothing writes into the real ``Working folder/``. Enable
  with ``PROSY_RUN_COOKBOOK=1 pytest tests/test_cookbook.py``. Run it before
  committing changes to the cookbook or to any CLI's arguments.

Mark a block that must not be executed with ``<!-- cookbook:skip -->`` on the
line before its opening fence (used for the install block).
"""

from __future__ import annotations

import ast
import os
import re
import shutil
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path

import pytest

_ROOT = Path(__file__).resolve().parents[1]
COOKBOOK = _ROOT / "docs" / "cookbook.md"

SKIP_MARKER = "<!-- cookbook:skip -->"
_FENCE = re.compile(r"^```(\w+)?\s*$")
#: Data the cookbook uses that is gitignored, so absent from a fresh clone.
_OPTIONAL_PREFIX = "Working folder/"


@dataclass
class Block:
    language: str
    body: str
    line: int
    skip: bool

    @property
    def commands(self) -> list[str]:
        """Logical shell commands, with line continuations joined."""
        joined = self.body.replace("\\\n", " ")
        return [ln.strip() for ln in joined.splitlines()
                if ln.strip() and not ln.strip().startswith("#")]


def parse_blocks(markdown: str) -> list[Block]:
    blocks: list[Block] = []
    lines = markdown.splitlines()
    i = 0
    while i < len(lines):
        match = _FENCE.match(lines[i])
        if not match:
            i += 1
            continue
        language = match.group(1) or ""
        skip = i > 0 and lines[i - 1].strip() == SKIP_MARKER
        start = i + 1
        j = start
        while j < len(lines) and not _FENCE.match(lines[j]):
            j += 1
        blocks.append(Block(language, "\n".join(lines[start:j]), start + 1, skip))
        i = j + 1
    return blocks


@pytest.fixture(scope="module")
def blocks() -> list[Block]:
    assert COOKBOOK.exists(), f"missing {COOKBOOK}"
    return parse_blocks(COOKBOOK.read_text(encoding="utf-8"))


@pytest.fixture(scope="module")
def bash_blocks(blocks) -> list[Block]:
    return [b for b in blocks if b.language == "bash"]


# --------------------------------------------------------------------------- #
# Lint                                                                         #
# --------------------------------------------------------------------------- #


def test_the_cookbook_has_runnable_blocks(bash_blocks):
    runnable = [b for b in bash_blocks if not b.skip]
    assert len(runnable) >= 8, "cookbook lost its command blocks?"


def _script_invocations(command: str) -> list[tuple[str, list[str]]]:
    """``[(script path, flags)]`` for each ``python scripts/….py …`` call."""
    out = []
    for match in re.finditer(r"python3?\s+(scripts/\S+\.py)((?:\s+\S+)*)", command):
        flags = re.findall(r"(--[a-z0-9][a-z0-9-]*)", match.group(2))
        out.append((match.group(1), flags))
    return out


def test_every_referenced_script_exists(bash_blocks):
    missing = {script for b in bash_blocks for cmd in b.commands
               for script, _ in _script_invocations(cmd)
               if not (_ROOT / script).exists()}
    assert not missing, f"cookbook references non-existent scripts: {sorted(missing)}"


@pytest.fixture(scope="module")
def help_text() -> dict[str, str]:
    cache: dict[str, str] = {}

    def get(script: str) -> str:
        if script not in cache:
            result = subprocess.run(
                [sys.executable, str(_ROOT / script), "--help"],
                capture_output=True, text=True, timeout=120, cwd=_ROOT)
            assert result.returncode == 0, f"{script} --help failed: {result.stderr}"
            cache[script] = result.stdout
        return cache[script]

    return get


def test_every_documented_flag_is_accepted(bash_blocks, help_text):
    """The check that would have caught --expect-prefix on the scan CLIs."""
    bad: list[str] = []
    for block in bash_blocks:
        for command in block.commands:
            for script, flags in _script_invocations(command):
                if not (_ROOT / script).exists():
                    continue
                usage = help_text(script)
                bad += [f"line {block.line}: {script} has no {flag}"
                        for flag in flags if flag not in usage]
    assert not bad, "cookbook uses flags that do not exist:\n  " + "\n  ".join(bad)


def test_python_blocks_compile(blocks):
    for block in blocks:
        if block.language != "python" or block.skip:
            continue
        try:
            ast.parse(block.body)
        except SyntaxError as exc:  # pragma: no cover - only on a broken doc
            pytest.fail(f"python block at line {block.line} does not parse: {exc}")


def test_repo_paths_referenced_by_the_cookbook_exist(bash_blocks, blocks):
    """Paths under version control must exist; lab data may legitimately not."""
    text = COOKBOOK.read_text(encoding="utf-8")
    referenced = set(re.findall(r'"?((?:data|scripts|prosy|tests)/[\w./-]+)"?', text))
    missing = sorted(p for p in referenced if not (_ROOT / p).exists())
    assert not missing, f"cookbook references missing repo paths: {missing}"


def test_lab_data_paths_are_flagged_as_optional():
    """`Working folder/` is gitignored, so those commands cannot run on a fresh
    clone. The cookbook has to say so, or a new user hits confusing failures."""
    text = COOKBOOK.read_text(encoding="utf-8")
    if _OPTIONAL_PREFIX in text:
        assert "gitignored" in text or "not in the repo" in text, (
            "cookbook uses 'Working folder/' paths without warning that they are "
            "lab data absent from a fresh clone")


# --------------------------------------------------------------------------- #
# Smoke (opt-in)                                                               #
# --------------------------------------------------------------------------- #

run_cookbook = pytest.mark.skipif(
    os.environ.get("PROSY_RUN_COOKBOOK") != "1",
    reason="set PROSY_RUN_COOKBOOK=1 to actually execute the cookbook commands",
)


@pytest.fixture(scope="module")
def sandbox(tmp_path_factory) -> Path:
    """A copy of the repo, so executed commands cannot touch real lab data."""
    dest = tmp_path_factory.mktemp("cookbook_repo")
    for name in ("prosy", "scripts", "data", "docs"):
        shutil.copytree(_ROOT / name, dest / name)
    working = _ROOT / "Working folder"
    if working.exists():
        shutil.copytree(working, dest / "Working folder",
                        ignore=shutil.ignore_patterns("*.png", "*.svg", "*.pdf"))
    return dest


@run_cookbook
def test_every_bash_block_runs(bash_blocks, sandbox):
    failures: list[str] = []
    for block in bash_blocks:
        if block.skip:
            continue
        result = subprocess.run(
            ["bash", "-euo", "pipefail", "-c", block.body],
            capture_output=True, text=True, timeout=1800, cwd=sandbox,
            env={**os.environ, "PYTHONPATH": str(sandbox)},
        )
        if result.returncode != 0:
            tail = (result.stderr or result.stdout).strip().splitlines()[-4:]
            failures.append(f"line {block.line} exited {result.returncode}:\n    "
                            + "\n    ".join(tail))
    assert not failures, "cookbook commands failed:\n" + "\n".join(failures)


if __name__ == "__main__":
    raise SystemExit(pytest.main([__file__, "-q"]))
