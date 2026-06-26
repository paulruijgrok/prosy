"""Microplate coordinate helpers (96-well by default).

Wells are addressed as ``<Row><Col>`` (e.g. ``A1``). Rows are letters A.. and
columns are 1-based integers. Supports column-major or row-major fill order.
"""

from __future__ import annotations

from dataclasses import dataclass

_ROW_LETTERS = "ABCDEFGHIJKLMNOP"


@dataclass(frozen=True)
class PlateFormat:
    rows: int = 8
    cols: int = 12

    @property
    def size(self) -> int:
        return self.rows * self.cols

    def row_letter(self, row_idx: int) -> str:
        if not 0 <= row_idx < self.rows:
            raise ValueError(f"Row index {row_idx} out of range for {self.rows} rows.")
        return _ROW_LETTERS[row_idx]

    def well(self, row_idx: int, col_idx: int) -> str:
        """``row_idx``/``col_idx`` are 0-based; returns e.g. ``A1``."""
        if not 0 <= col_idx < self.cols:
            raise ValueError(f"Col index {col_idx} out of range for {self.cols} cols.")
        return f"{self.row_letter(row_idx)}{col_idx + 1}"

    def wells(self, order: str = "column") -> list[str]:
        """All wells in fill order. ``column`` = down each column first (A1,B1,
        ...,H1,A2,...); ``row`` = across each row first (A1,A2,...,A12,B1,...).
        """
        if order == "column":
            return [
                self.well(r, c)
                for c in range(self.cols)
                for r in range(self.rows)
            ]
        if order == "row":
            return [
                self.well(r, c)
                for r in range(self.rows)
                for c in range(self.cols)
            ]
        raise ValueError("order must be 'column' or 'row'")


PLATE_96 = PlateFormat(8, 12)
PLATE_384 = PlateFormat(16, 24)


def resolve_plate_format(plate: "PlateFormat | str | int") -> PlateFormat:
    """Accept a PlateFormat, a size (96/384) or a name ('96'/'384')."""
    if isinstance(plate, PlateFormat):
        return plate
    key = str(plate).strip().lower()
    if key in {"96", "96-well", "plate_96"}:
        return PLATE_96
    if key in {"384", "384-well", "plate_384"}:
        return PLATE_384
    raise ValueError(f"Unknown plate {plate!r}; use 96, 384 or a PlateFormat.")
