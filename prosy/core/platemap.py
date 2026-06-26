"""Plate-map rendering.

A dependency-free SVG renderer draws a plate grid with one cell per well, each
labeled and colour-coded by category (e.g. parent vs mutant, or by set). SVG is
scalable, renders anywhere, and needs nothing installed.

``write_png`` is optional: it uses matplotlib if available, else LibreOffice to
convert the SVG, else returns None. For richer lab graphics and robot picklists,
the optional Plateo package ([bio] extra) is the recommended route.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

from prosy.core.plate import PlateFormat, resolve_plate_format

# A readable categorical palette (fills) with matching darker strokes.
_PALETTE = [
    ("#dbeafe", "#3b82f6"),  # blue
    ("#dcfce7", "#22c55e"),  # green
    ("#fef9c3", "#eab308"),  # yellow
    ("#fee2e2", "#ef4444"),  # red
    ("#f3e8ff", "#a855f7"),  # purple
    ("#ffedd5", "#f97316"),  # orange
    ("#cffafe", "#06b6d4"),  # cyan
    ("#fce7f3", "#ec4899"),  # pink
]
_EMPTY_FILL = "#f8fafc"
_EMPTY_STROKE = "#cbd5e1"


@dataclass
class Cell:
    """Contents of one well for rendering."""

    top: str = ""          # primary label (e.g. NbID)
    bottom: str = ""       # secondary label (e.g. mutation)
    category: str = ""     # drives colour grouping
    emphasize: bool = False  # thicker border (e.g. parents)


def _esc(s: str) -> str:
    return (str(s).replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;"))


def _colors_for(categories: list[str]) -> dict[str, tuple[str, str]]:
    mapping: dict[str, tuple[str, str]] = {}
    i = 0
    for c in categories:
        if c not in mapping:
            mapping[c] = _PALETTE[i % len(_PALETTE)]
            i += 1
    return mapping


def render_svg(
    cells: dict[str, Cell],
    *,
    plate: PlateFormat | str | int = 96,
    title: str | None = None,
    cell_w: int = 70,
    cell_h: int = 50,
) -> str:
    """Render a plate map as an SVG string. ``cells`` maps well -> Cell."""
    plate = resolve_plate_format(plate)
    pad = 16
    label_gutter = 26          # space for row letters / column numbers
    title_h = 30 if title else 0
    legend_h = 30
    grid_w = plate.cols * cell_w
    grid_h = plate.rows * cell_h
    width = pad * 2 + label_gutter + grid_w
    height = pad * 2 + title_h + label_gutter + grid_h + legend_h

    ordered_cats: list[str] = []
    for w in plate.wells(order="row"):
        c = cells.get(w)
        if c and c.category and c.category not in ordered_cats:
            ordered_cats.append(c.category)
    palette = _colors_for(ordered_cats)

    x0 = pad + label_gutter
    y0 = pad + title_h + label_gutter
    parts: list[str] = []
    parts.append(
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" '
        f'viewBox="0 0 {width} {height}" font-family="Helvetica,Arial,sans-serif">'
    )
    parts.append(f'<rect width="{width}" height="{height}" fill="white"/>')
    if title:
        parts.append(
            f'<text x="{pad}" y="{pad + 18}" font-size="18" font-weight="bold" '
            f'fill="#0f172a">{_esc(title)}</text>'
        )

    # Column numbers
    for c in range(plate.cols):
        cx = x0 + c * cell_w + cell_w / 2
        parts.append(
            f'<text x="{cx:.1f}" y="{y0 - 8}" font-size="11" fill="#475569" '
            f'text-anchor="middle">{c + 1}</text>'
        )
    # Row letters
    for r in range(plate.rows):
        cy = y0 + r * cell_h + cell_h / 2
        parts.append(
            f'<text x="{x0 - 8}" y="{cy + 4:.1f}" font-size="11" fill="#475569" '
            f'text-anchor="end">{plate.row_letter(r)}</text>'
        )

    # Wells
    for r in range(plate.rows):
        for c in range(plate.cols):
            well = plate.well(r, c)
            cell = cells.get(well)
            x = x0 + c * cell_w
            y = y0 + r * cell_h
            if cell and cell.category:
                fill, stroke = palette[cell.category]
            else:
                fill, stroke = _EMPTY_FILL, _EMPTY_STROKE
            sw = 2.4 if (cell and cell.emphasize) else 1.0
            parts.append(
                f'<rect x="{x + 2}" y="{y + 2}" width="{cell_w - 4}" height="{cell_h - 4}" '
                f'rx="6" fill="{fill}" stroke="{stroke}" stroke-width="{sw}"/>'
            )
            if cell and cell.top:
                weight = "bold" if cell.emphasize else "normal"
                parts.append(
                    f'<text x="{x + cell_w / 2:.1f}" y="{y + cell_h / 2 - 2:.1f}" '
                    f'font-size="10" font-weight="{weight}" fill="#0f172a" '
                    f'text-anchor="middle">{_esc(cell.top)}</text>'
                )
            if cell and cell.bottom:
                parts.append(
                    f'<text x="{x + cell_w / 2:.1f}" y="{y + cell_h / 2 + 11:.1f}" '
                    f'font-size="9" fill="#475569" text-anchor="middle">'
                    f'{_esc(cell.bottom)}</text>'
                )

    # Legend
    ly = y0 + grid_h + 18
    lx = x0
    for cat in ordered_cats:
        fill, stroke = palette[cat]
        parts.append(
            f'<rect x="{lx}" y="{ly - 10}" width="12" height="12" rx="3" '
            f'fill="{fill}" stroke="{stroke}"/>'
        )
        parts.append(
            f'<text x="{lx + 18}" y="{ly}" font-size="11" fill="#334155">{_esc(cat)}</text>'
        )
        lx += 22 + 8 * len(cat) + 18
    parts.append("</svg>")
    return "\n".join(parts)


def write_svg(cells: dict[str, Cell], path: str | Path, **kwargs) -> Path:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(render_svg(cells, **kwargs), encoding="utf-8")
    return path


def write_png(cells: dict[str, Cell], path: str | Path, **kwargs) -> Path | None:
    """Write a PNG if possible (matplotlib, else LibreOffice from SVG). Else None."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    svg = render_svg(cells, **kwargs)

    # 1) matplotlib path (renders the SVG-equivalent grid natively).
    try:
        return _png_via_matplotlib(cells, path, **kwargs)
    except Exception:
        pass

    # 2) LibreOffice conversion from a temporary SVG.
    import shutil
    import subprocess

    soffice = shutil.which("libreoffice") or shutil.which("soffice")
    if soffice:
        svg_path = path.with_suffix(".svg")
        svg_path.write_text(svg, encoding="utf-8")
        try:
            subprocess.run(
                [soffice, "--headless", "--convert-to", "png", "--outdir",
                 str(path.parent), str(svg_path)],
                check=True, capture_output=True, timeout=120,
            )
            if path.exists():
                return path
        except (subprocess.SubprocessError, OSError):
            return None
    return None


def _png_via_matplotlib(cells, path, *, plate=96, title=None, **_ignore) -> Path:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.patches as mpatches
    import matplotlib.pyplot as plt

    pf = resolve_plate_format(plate)
    ordered_cats: list[str] = []
    for w in pf.wells(order="row"):
        c = cells.get(w)
        if c and c.category and c.category not in ordered_cats:
            ordered_cats.append(c.category)
    palette = _colors_for(ordered_cats)

    fig_w = max(6.0, pf.cols * 0.7)
    fig_h = max(4.0, pf.rows * 0.6)
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    ax.set_xlim(0, pf.cols)
    ax.set_ylim(0, pf.rows)
    ax.invert_yaxis()
    ax.set_xticks([c + 0.5 for c in range(pf.cols)])
    ax.set_xticklabels([str(c + 1) for c in range(pf.cols)], fontsize=8)
    ax.set_yticks([r + 0.5 for r in range(pf.rows)])
    ax.set_yticklabels([pf.row_letter(r) for r in range(pf.rows)], fontsize=8)
    ax.xaxis.tick_top()
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.tick_params(length=0)

    for r in range(pf.rows):
        for c in range(pf.cols):
            cell = cells.get(pf.well(r, c))
            if cell and cell.category:
                fill, stroke = palette[cell.category]
            else:
                fill, stroke = _EMPTY_FILL, _EMPTY_STROKE
            lw = 2.0 if (cell and cell.emphasize) else 0.6
            ax.add_patch(mpatches.FancyBboxPatch(
                (c + 0.06, r + 0.06), 0.88, 0.88,
                boxstyle="round,pad=0.0,rounding_size=0.08",
                facecolor=fill, edgecolor=stroke, linewidth=lw))
            if cell and cell.top:
                ax.text(c + 0.5, r + 0.42, cell.top, ha="center", va="center",
                        fontsize=6.5, fontweight="bold" if cell.emphasize else "normal")
            if cell and cell.bottom:
                ax.text(c + 0.5, r + 0.66, cell.bottom, ha="center", va="center",
                        fontsize=5.5, color="#475569")
    if title:
        ax.set_title(title, fontsize=12, fontweight="bold", loc="left", pad=18)
    if ordered_cats:
        handles = [mpatches.Patch(facecolor=palette[c][0], edgecolor=palette[c][1],
                                  label=c) for c in ordered_cats]
        ax.legend(handles=handles, loc="upper center",
                  bbox_to_anchor=(0.5, -0.04), ncol=len(ordered_cats), fontsize=8,
                  frameon=False)
    fig.tight_layout()
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return path
