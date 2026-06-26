"""Tests for the plate layout engine, incl. the natural extension configs."""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from prosy.core import layout
from prosy.core.plate import PLATE_96, PLATE_384


def _all_wells(groups):
    return [w for g in groups for w in g.wells]


def test_column_plate_12x8():
    groups = layout.layout_groups(12, 8, plate=96, orientation="column")
    assert len(groups) == 12
    assert groups[0].wells[:3] == ["A1", "B1", "C1"]
    assert groups[0].wells[-1] == "H1"
    assert groups[1].wells[0] == "A2"
    assert len(set(_all_wells(groups))) == 96


def test_row_plate_8x12():
    groups = layout.layout_groups(8, 12, plate=96, orientation="row")
    assert len(groups) == 8
    assert groups[0].wells[:3] == ["A1", "A2", "A3"]
    assert groups[0].wells[-1] == "A12"
    assert groups[1].wells[0] == "B1"
    assert len(set(_all_wells(groups))) == 96


def test_16_groups_of_6_two_per_row():
    # 16 parents x (1+5) = 96, row-major, 2 groups per row
    groups = layout.layout_groups(16, 6, plate=96, orientation="row")
    assert len(groups) == 16
    assert groups[0].wells == ["A1", "A2", "A3", "A4", "A5", "A6"]
    assert groups[1].wells == ["A7", "A8", "A9", "A10", "A11", "A12"]
    assert groups[2].wells[0] == "B1"
    assert len(set(_all_wells(groups))) == 96


def test_24_groups_of_4_two_per_column():
    # 24 parents x (1+3) = 96, column-major, 2 groups per column
    groups = layout.layout_groups(24, 4, plate=96, orientation="column")
    assert len(groups) == 24
    assert groups[0].wells == ["A1", "B1", "C1", "D1"]
    assert groups[1].wells == ["E1", "F1", "G1", "H1"]
    assert groups[2].wells[0] == "A2"
    assert len(set(_all_wells(groups))) == 96


def test_384_column():
    # 48 parents x (1+7) = 384 on a 384-well plate, column orientation
    groups = layout.layout_groups(48, 8, plate=384, orientation="column")
    assert len(groups) == 48
    assert len(set(_all_wells(groups))) == 384
    assert groups[0].wells[0] == "A1" and groups[0].wells[-1] == "H1"


def test_384_row_24_per_row():
    groups = layout.layout_groups(16, 24, plate=384, orientation="row")
    assert len(groups) == 16
    assert groups[0].wells[0] == "A1" and groups[0].wells[-1] == "A24"
    assert groups[1].wells[0] == "B1"


def test_overflow_raises():
    try:
        layout.layout_groups(13, 8, plate=96, orientation="column")
    except layout.LayoutError as e:
        assert "exceeds" in str(e)
    else:
        raise AssertionError("overflow not detected")


def test_straddle_requires_opt_in():
    # group size 5 does not tile 8 (column) -> error unless allowed
    try:
        layout.layout_groups(10, 5, plate=96, orientation="column", require_aligned=True)
    except layout.LayoutError as e:
        assert "straddle" in str(e)
    else:
        raise AssertionError("straddle not flagged")
    groups = layout.layout_groups(10, 5, plate=96, orientation="column",
                                  require_aligned=False)
    assert len(groups) == 10


def test_resolve_plate():
    assert layout.resolve_plate(96) is PLATE_96
    assert layout.resolve_plate("384") is PLATE_384


def _run_all():
    fns = [v for k, v in sorted(globals().items()) if k.startswith("test_") and callable(v)]
    for fn in fns:
        fn()
        print(f"  ok  {fn.__name__}")
    print(f"\n{len(fns)}/{len(fns)} layout tests passed")


if __name__ == "__main__":
    _run_all()
