"""Plate layout engine.

Places *groups* of wells (e.g. one parent + its M mutants) onto a plate in a
chosen fill order. The same engine covers all the practical combinations:

* group filling one column   (8 wells)  -> 12 groups on a 96-well plate
* group filling one row      (12 wells) -> 8 groups on a 96-well plate
* several groups per row/column (e.g. group size 6 -> 2 per row -> 16 groups)
* 384-well plates (16x24) with the same rules

A "group" here is just a contiguous run of wells of a given size, taken in the
plate's fill order (row-major or column-major). Higher layers decide what each
well holds.
"""

from __future__ import annotations

from dataclasses import dataclass

from prosy.core.plate import PLATE_96, PLATE_384, PlateFormat, resolve_plate_format

ORIENTATIONS = ("column", "row")


@dataclass(frozen=True)
class GroupLayout:
    """The wells assigned to one group, in fill order."""

    index: int          # 0-based group number, in fill order
    wells: list[str]


class LayoutError(ValueError):
    pass


def resolve_plate(plate: PlateFormat | str | int) -> PlateFormat:
    """Accept a PlateFormat, a size (96/384) or a name ('96'/'384')."""
    try:
        return resolve_plate_format(plate)
    except ValueError as e:
        raise LayoutError(str(e)) from e


def _fill_dimension(plate: PlateFormat, orientation: str) -> int:
    """Length of the run that fills first: rows-per-column (column-major) or
    cols-per-row (row-major)."""
    return plate.rows if orientation == "column" else plate.cols


def layout_groups(
    n_groups: int,
    group_size: int,
    *,
    plate: PlateFormat | str | int = PLATE_96,
    orientation: str = "column",
    require_aligned: bool = True,
) -> list[GroupLayout]:
    """Assign ``n_groups`` contiguous groups of ``group_size`` wells.

    Parameters
    ----------
    orientation:
        ``"column"`` fills down each column first (A1,B1,...,H1,A2,...);
        ``"row"`` fills across each row first (A1,A2,...,A12,B1,...).
    require_aligned:
        If True (default), require that ``group_size`` divides the fill
        dimension (plate rows for column-major, plate cols for row-major) so no
        group straddles a row/column boundary. Set False to allow straddling.
    """
    if orientation not in ORIENTATIONS:
        raise LayoutError(f"orientation must be one of {ORIENTATIONS}")
    plate = resolve_plate(plate)
    if n_groups <= 0 or group_size <= 0:
        raise LayoutError("n_groups and group_size must be positive.")

    needed = n_groups * group_size
    if needed > plate.size:
        raise LayoutError(
            f"{n_groups} groups x {group_size} wells = {needed} exceeds plate "
            f"capacity {plate.size}."
        )

    fill_dim = _fill_dimension(plate, orientation)
    if require_aligned and fill_dim % group_size != 0 and group_size % fill_dim != 0:
        raise LayoutError(
            f"group_size {group_size} does not tile the {orientation} dimension "
            f"({fill_dim}); groups would straddle boundaries. Pass "
            f"require_aligned=False to allow this."
        )

    wells = plate.wells(order=orientation)
    groups: list[GroupLayout] = []
    for g in range(n_groups):
        start = g * group_size
        groups.append(GroupLayout(index=g, wells=wells[start : start + group_size]))
    return groups


def describe_layout(
    n_groups: int,
    group_size: int,
    *,
    plate: PlateFormat | str | int = PLATE_96,
    orientation: str = "column",
) -> str:
    """One-line human summary, e.g. '8 groups x 12 wells, row-major on 96-well
    (2 group(s) per row)'."""
    plate = resolve_plate(plate)
    fill_dim = _fill_dimension(plate, orientation)
    per_line = fill_dim / group_size if group_size <= fill_dim else None
    unit = "row" if orientation == "row" else "column"
    if group_size == fill_dim:
        shape = f"1 group per {unit}"
    elif per_line and per_line == int(per_line):
        shape = f"{int(per_line)} group(s) per {unit}"
    else:
        shape = "groups straddle boundaries"
    return (
        f"{n_groups} groups x {group_size} wells, {orientation}-major on "
        f"{plate.size}-well ({shape})"
    )
