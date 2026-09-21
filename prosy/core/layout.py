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

from prosy.core.plate import PLATE_96, PlateFormat, resolve_plate_format

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


def layout_lines(
    sizes: list[int],
    *,
    plate: PlateFormat | str | int = PLATE_96,
    orientation: str = "row",
    balance_overflow: bool = True,
) -> list[GroupLayout]:
    """Give every group its own line (row or column), whatever its size.

    Unlike :func:`layout_groups`, group sizes may differ. Each group starts on a
    fresh line and never shares one with another group, so a plate map reads as
    "one group per row". A group larger than the line length spills onto as many
    further lines as it needs.

    Parameters
    ----------
    sizes:
        Number of wells each group needs, in plate order.
    orientation:
        ``"row"`` (default) gives each group its own row and fills left to
        right; ``"column"`` gives each group its own column, filling downwards.
    balance_overflow:
        For a group that spans several lines, spread its members as evenly as
        possible across them (14 over two 12-well rows -> 7 + 7). Set False to
        fill each line to capacity before starting the next (-> 12 + 2).

    Raises
    ------
    LayoutError
        If the groups need more lines than the plate has.
    """
    if orientation not in ORIENTATIONS:
        raise LayoutError(f"orientation must be one of {ORIENTATIONS}")
    plate = resolve_plate(plate)
    if any(s <= 0 for s in sizes):
        raise LayoutError("Every group size must be positive.")

    # A "line" is a row (row-major) or a column (column-major): `line_len` wells
    # long, and there are `n_lines` of them on the plate.
    line_len = _fill_dimension(plate, orientation)
    n_lines = plate.rows if orientation == "row" else plate.cols

    lines_needed = sum(-(-s // line_len) for s in sizes)  # ceil division
    if lines_needed > n_lines:
        unit = "row" if orientation == "row" else "column"
        raise LayoutError(
            f"{len(sizes)} groups need {lines_needed} {unit}(s) "
            f"(sizes {sizes}, {line_len} wells per {unit}) but the "
            f"{plate.size}-well plate has only {n_lines}."
        )

    all_wells = plate.wells(order=orientation)
    groups: list[GroupLayout] = []
    line = 0
    for index, size in enumerate(sizes):
        spans = -(-size // line_len)
        if balance_overflow and spans > 1:
            base, extra = divmod(size, spans)
            per_line = [base + (1 if i < extra else 0) for i in range(spans)]
        else:
            per_line = [min(line_len, size - i * line_len) for i in range(spans)]
        wells: list[str] = []
        for offset, count in enumerate(per_line):
            start = (line + offset) * line_len
            wells.extend(all_wells[start : start + count])
        groups.append(GroupLayout(index=index, wells=wells))
        line += spans
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
