"""Tests for scripts/nanobodies/make_db_rows.py."""

from __future__ import annotations

import csv
import sys
from pathlib import Path

import pytest

_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(_ROOT))
sys.path.insert(0, str(_ROOT / "scripts" / "nanobodies"))

openpyxl = pytest.importorskip("openpyxl")

import make_db_rows as mdr  # noqa: E402

# Two toy parents, 20 aa each, with CDRs at known offsets.
P1 = "QVQLVESGGGAAADEFGHIK"
P2 = "QVQLQESGGGCCCDEFGHIK"

SHEET_HEADER = ["DB ID", "Nanobody ID", "", "Target", "Target spec", "Source",
                "Set name", "Length (AA)", "AA Sequence", "CDR 1", "CDR 2",
                "CDR 3", "Alias", "Old Alias", "Parent", "Linked plasmids",
                "Set rank", "On target affinity", "Notes"]
SHEET_SUB = ["", "Short", "Long"] + [""] * 16


def _make_sheet(path: Path) -> None:
    wb = openpyxl.Workbook()
    ws = wb.active
    ws.title = "Nanobodies"
    ws.append(SHEET_HEADER)
    ws.append(SHEET_SUB)
    ws.append([1, "Nb01", None, "PfLDH", "recb", "P", "ADL1", 20, P1,
               "AAA", "DEF", "GHIK", "P4102_95", None, None, None, 1, None, None])
    ws.append([2, "Nb02", None, "PfLDH", "recb", "Y", "KRU1", 20, P2,
               "CCC", "DEF", "GHIK", None, None, None, None, 2, None, None])
    wb.save(path)


def _make_mapping(path: Path) -> None:
    rows = [
        # Nb01 parent, plus a mutation inside CDR 1 (A11D) and one outside (Q1E).
        dict(Well="A1", NbID="Nb01", NbID_parent="Nb01", Parent="Nb01", Set="ADL1",
             Mutation="", Score="", Type="parent", Protein=P1),
        dict(Well="A2", NbID="Nb03", NbID_parent="Nb01_A11D", Parent="Nb01",
             Set="ADL1", Mutation="A11D", Score="7", Type="mutant",
             Protein=P1[:10] + "D" + P1[11:]),
        dict(Well="A3", NbID="Nb04", NbID_parent="Nb01_Q1E", Parent="Nb01",
             Set="ADL1", Mutation="Q1E", Score="4", Type="mutant",
             Protein="E" + P1[1:]),
        dict(Well="B1", NbID="Nb05", NbID_parent="Nb02_C11W", Parent="Nb02",
             Set="KRU1", Mutation="C11W", Score="3", Type="mutant",
             Protein=P2[:10] + "W" + P2[11:]),
    ]
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)


@pytest.fixture()
def run_dir(tmp_path: Path) -> Path:
    _make_sheet(tmp_path / "Malaria DX plasmids.xlsx")
    _make_mapping(tmp_path / "t_nanobody_mapping.csv")
    return tmp_path


def _read_out(path: Path) -> tuple[list[str], list[dict]]:
    ws = openpyxl.load_workbook(path).active
    rows = list(ws.iter_rows(values_only=True))
    cols = [c for c in rows[0]]
    return cols, [dict(zip(cols, r)) for r in rows[1:]]


def _run(run_dir: Path, *extra: str) -> tuple[list[str], list[dict]]:
    assert mdr.main(["--data-dir", str(run_dir), "--stamp", "t", *extra]) == 0
    return _read_out(run_dir / "t_nanobody_db_rows.xlsx")


def test_columns_mirror_the_sheet(run_dir: Path) -> None:
    cols, _ = _run(run_dir)
    assert cols == ["DB ID", "Nanobody ID", "Long", "Target", "Target spec",
                    "Source", "Set name", "Length (AA)", "AA Sequence", "CDR 1",
                    "CDR 2", "CDR 3", "Alias", "Old Alias", "Parent",
                    "Linked plasmids", "Set rank", "On target affinity", "Notes"]


def test_formula_columns_left_blank(run_dir: Path) -> None:
    _, rows = _run(run_dir)
    for r in rows:
        for col in mdr.FORMULA_COLUMNS:
            assert r[col] in (None, ""), f"{col} must stay empty for the sheet's formula"


def test_mutants_only_by_default_with_continuing_db_ids(run_dir: Path) -> None:
    _, rows = _run(run_dir)
    assert [r["Nanobody ID"] for r in rows] == ["Nb03", "Nb04", "Nb05"]
    assert [r["DB ID"] for r in rows] == [3, 4, 5]  # sheet's max DB ID (2) + 1


def test_include_all_adds_the_parents_in_plate_order(run_dir: Path) -> None:
    _, rows = _run(run_dir, "--include", "all")
    assert [r["Nanobody ID"] for r in rows] == ["Nb01", "Nb03", "Nb04", "Nb05"]


def test_sequences_and_inherited_metadata_come_from_the_run(run_dir: Path) -> None:
    _, rows = _run(run_dir)
    by_id = {r["Nanobody ID"]: r for r in rows}
    assert by_id["Nb03"]["AA Sequence"] == P1[:10] + "D" + P1[11:]
    assert by_id["Nb03"]["Parent"] == "Nb01"
    assert by_id["Nb03"]["Target"] == "PfLDH"
    assert by_id["Nb03"]["Source"] == "P"
    assert by_id["Nb03"]["Set name"] == "ADL1"
    assert by_id["Nb05"]["Source"] == "Y"
    assert by_id["Nb05"]["Set name"] == "KRU1"


def test_alias_extends_the_parent_alias_and_is_blank_without_one(run_dir: Path) -> None:
    _, rows = _run(run_dir)
    by_id = {r["Nanobody ID"]: r for r in rows}
    assert by_id["Nb03"]["Alias"] == "P4102_95_A11D"
    assert by_id["Nb05"]["Alias"] in (None, "")  # Nb02 has no alias


def test_cdr_offsets_mode_reflects_a_mutation_inside_a_cdr(run_dir: Path) -> None:
    _, rows = _run(run_dir)
    by_id = {r["Nanobody ID"]: r for r in rows}
    # CDR 1 spans positions 11-13; A11D changes its first residue.
    assert by_id["Nb03"]["CDR 1"] == "DAA"
    assert by_id["Nb04"]["CDR 1"] == "AAA"      # Q1E is outside every CDR
    assert by_id["Nb05"]["CDR 1"] == "WCC"


def test_cdr_inherit_and_blank_modes(run_dir: Path) -> None:
    _, rows = _run(run_dir, "--cdr-mode", "inherit")
    assert {r["CDR 1"] for r in rows} == {"AAA", "CCC"}
    _, rows = _run(run_dir, "--cdr-mode", "blank")
    assert all(r["CDR 1"] in (None, "") for r in rows)


def test_db_id_start_override(run_dir: Path) -> None:
    _, rows = _run(run_dir, "--db-id-start", "500")
    assert [r["DB ID"] for r in rows] == [500, 501, 502]


def test_missing_parent_in_sheet_fails_loudly(run_dir: Path) -> None:
    path = run_dir / "t_nanobody_mapping.csv"
    text = path.read_text().replace("Nb02", "Nb99")
    path.write_text(text)
    with pytest.raises(SystemExit, match="Nb99"):
        mdr.main(["--data-dir", str(run_dir), "--stamp", "t"])


def test_ambiguous_cdr_is_copied_verbatim_and_warned() -> None:
    parent = "AAABAAACDE"          # "AAA" occurs twice
    cdrs, warns = mdr.cdrs_at_parent_offsets(
        parent, {"CDR 1": "AAA", "CDR 2": "CDE", "CDR 3": ""}, "AAABAAACDF")
    assert cdrs["CDR 1"] == "AAA"
    assert cdrs["CDR 2"] == "CDF"
    assert cdrs["CDR 3"] == ""
    assert any("more than once" in w for w in warns)


def test_cdr_absent_from_parent_is_copied_verbatim_and_warned() -> None:
    cdrs, warns = mdr.cdrs_at_parent_offsets(
        "AAABBB", {"CDR 1": "ZZZ", "CDR 2": "", "CDR 3": ""}, "AAABBC")
    assert cdrs["CDR 1"] == "ZZZ"
    assert any("not found" in w for w in warns)
