"""Tests for layout.layout_lines: one variable-size group per row/column."""

from __future__ import annotations

import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from prosy.core.layout import LayoutError, layout_lines  # noqa: E402

# The 260917 enzyme design set: six groups, one of them larger than a row.
DESIGN_SET_SIZES = [12, 14, 5, 9, 5, 6]


def test_each_group_starts_on_a_fresh_row():
    groups = layout_lines([5, 3, 7])
    assert [g.wells[0] for g in groups] == ["A1", "B1", "C1"]
    assert groups[0].wells == ["A1", "A2", "A3", "A4", "A5"]
    assert groups[1].wells == ["B1", "B2", "B3"]
    assert len(groups[2].wells) == 7


def test_no_two_groups_share_a_row():
    groups = layout_lines(DESIGN_SET_SIZES)
    rows_per_group = [{w[0] for w in g.wells} for g in groups]
    for i, a in enumerate(rows_per_group):
        for b in rows_per_group[i + 1:]:
            assert not a & b


def test_oversized_group_is_balanced_across_rows_by_default():
    groups = layout_lines(DESIGN_SET_SIZES)
    fourteen = groups[1].wells
    assert len(fourteen) == 14
    assert [w for w in fourteen if w[0] == "B"] == [f"B{i}" for i in range(1, 8)]
    assert [w for w in fourteen if w[0] == "C"] == [f"C{i}" for i in range(1, 8)]


def test_fill_rows_packs_to_capacity_instead():
    groups = layout_lines(DESIGN_SET_SIZES, balance_overflow=False)
    fourteen = groups[1].wells
    assert len([w for w in fourteen if w[0] == "B"]) == 12
    assert [w for w in fourteen if w[0] == "C"] == ["C1", "C2"]


def test_design_set_occupies_seven_rows_in_order():
    groups = layout_lines(DESIGN_SET_SIZES)
    used = sorted({w[0] for g in groups for w in g.wells})
    assert used == list("ABCDEFG")           # row H stays free
    assert [len(g.wells) for g in groups] == DESIGN_SET_SIZES
    assert sum(len(g.wells) for g in groups) == 51


def test_wells_are_unique_across_the_whole_layout():
    wells = [w for g in layout_lines(DESIGN_SET_SIZES) for w in g.wells]
    assert len(wells) == len(set(wells))


def test_column_orientation_gives_each_group_its_own_column():
    groups = layout_lines([8, 3, 10], orientation="column")
    assert groups[0].wells == [f"{r}1" for r in "ABCDEFGH"]
    assert groups[1].wells == ["A2", "B2", "C2"]
    # 10 in an 8-row column spills to a second column, balanced 5/5.
    assert groups[2].wells == [f"{r}3" for r in "ABCDE"] + [f"{r}4" for r in "ABCDE"]


def test_384_well_plate():
    groups = layout_lines([24, 30], plate=384)
    assert len(groups[0].wells) == 24 and groups[0].wells[-1] == "A24"
    assert {w[0] for w in groups[1].wells} == {"B", "C"}


def test_too_many_rows_is_an_error():
    with pytest.raises(LayoutError, match="only 8"):
        layout_lines([1] * 9)


def test_capacity_counts_rows_spanned_not_wells():
    layout_lines([12] * 5 + [25])           # 5 rows + ceil(25/12)=3 -> exactly 8
    with pytest.raises(LayoutError):
        layout_lines([12] * 6 + [25])       # 6 + 3 = 9 rows


def test_non_positive_group_size_rejected():
    with pytest.raises(LayoutError):
        layout_lines([5, 0])


def test_bad_orientation_rejected():
    with pytest.raises(LayoutError):
        layout_lines([5], orientation="diagonal")


if __name__ == "__main__":
    raise SystemExit(pytest.main([__file__, "-q"]))
