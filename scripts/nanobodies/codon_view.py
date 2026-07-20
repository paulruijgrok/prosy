#!/usr/bin/env python3
"""Codon-aligned view of coding DNA sequences for easy manual inspection.

For each variant in a nanobody mapping CSV, render the in-frame coding sequence
(``Coding_DNA``) as blocks of codons with the amino acid each codon encodes
printed directly beneath it, plus a residue-number ruler. Two output formats are
written: a monospaced plain-text file and a colour-coded HTML file (amino acids
grouped by chemical class).

Example:
    conda activate DataAnalysis
    python codon_view.py \
        --input "../../Working folder/260720_NanobodyMuts_gc75/260626_nanobody_mapping.csv"

Layout styles (``--style``):

    pipe (default)          space
    |CAG|GTG|CAG|           CAG GTG CAG
    | Q | V | Q |            Q   V   Q
"""
from __future__ import annotations

import argparse
import csv
import html as html_mod
import sys
from pathlib import Path

_HERE = Path(__file__).resolve()
sys.path.insert(0, str(_HERE.parents[2]))  # repo root -> 'prosy'

from prosy.core.sequence import CODON_TABLE, translate  # noqa: E402

# Amino-acid chemical classes -> (label, hex colour) for the HTML view.
AA_CLASS = {
    **{a: "hydrophobic" for a in "AVLIMFW"},
    **{a: "special" for a in "GP"},
    **{a: "polar" for a in "STNQCY"},
    **{a: "positive" for a in "KRH"},
    **{a: "negative" for a in "DE"},
    "*": "stop",
}
CLASS_COLOR = {
    "hydrophobic": "#8a6d3b",   # brown/tan
    "special": "#6b6b6b",       # grey (Gly/Pro)
    "polar": "#2e7d32",         # green
    "positive": "#1565c0",      # blue
    "negative": "#c62828",      # red
    "stop": "#000000",
    "other": "#333333",
}


def codons_of(dna: str) -> list[str]:
    """Split an in-frame DNA string into 3-nt codons (trailing partial dropped)."""
    dna = dna.upper()
    return [dna[i:i + 3] for i in range(0, len(dna) - len(dna) % 3, 3)]


def aa_of(codon: str) -> str:
    return CODON_TABLE.get(codon, "?")


def _ruler(n: int, start_res: int, offset: int, stride: int, width: int) -> str:
    """Residue-number ruler aligned over codon fields (numbers every 5th + first)."""
    line = [" "] * (offset + n * stride + 8)
    for k in range(n):
        res = start_res + k + 1
        if k == 0 or res % 5 == 0:
            s = str(res)
            center = offset + k * stride + width // 2
            pos = center - (len(s) - 1) // 2
            for j, ch in enumerate(s):
                if 0 <= pos + j < len(line):
                    line[pos + j] = ch
    return "".join(line).rstrip()


def render_text_block(codons: list[str], aas: list[str], start_res: int,
                      style: str) -> list[str]:
    """One wrapped block (ruler + codon line + amino-acid line) as text lines."""
    if style == "space":
        codon_line = " ".join(codons)
        aa_line = " ".join(f"{a:^3}" for a in aas)
        ruler = _ruler(len(codons), start_res, offset=0, stride=4, width=3)
    else:  # pipe
        codon_line = "|" + "|".join(codons) + "|"
        aa_line = "|" + "|".join(f"{a:^3}" for a in aas) + "|"
        ruler = _ruler(len(codons), start_res, offset=1, stride=4, width=3)
    return [ruler, codon_line, aa_line]


def variant_text(rec: dict, seq: str, codons_per_line: int, style: str,
                 mismatch: str | None) -> str:
    codons = codons_of(seq)
    aas = [aa_of(c) for c in codons]
    header = _variant_header(rec, len(codons), mismatch)
    lines = [header]
    for i in range(0, len(codons), codons_per_line):
        block = render_text_block(codons[i:i + codons_per_line],
                                  aas[i:i + codons_per_line], i, style)
        lines.extend(block)
        lines.append("")  # blank line between blocks
    return "\n".join(lines).rstrip() + "\n"


def _variant_header(rec: dict, n_codons: int, mismatch: str | None) -> str:
    bits = []
    for key in ("NbID_parent", "NbID", "Well"):
        if rec.get(key):
            bits.append(str(rec[key]))
    flag = f"   [!] {mismatch}" if mismatch else ""
    return f">>> {'  '.join(bits)}  ({n_codons} codons){flag}"


def _translation_mismatch(seq: str, protein: str | None) -> str | None:
    """Return a human message if translation != given protein, else None."""
    if not protein:
        return None
    try:
        got = translate(seq)
    except Exception as e:  # noqa: BLE001
        return f"cannot translate: {e}"
    got = got.rstrip("*")
    want = protein.strip().upper().rstrip("*")
    if got != want:
        return "translation does not match Protein column"
    return None


# --------------------------------------------------------------------------- #
# HTML                                                                         #
# --------------------------------------------------------------------------- #

def _span(text: str, aa: str) -> str:
    cls = AA_CLASS.get(aa, "other")
    color = CLASS_COLOR.get(cls, CLASS_COLOR["other"])
    return f'<span style="color:{color}">{html_mod.escape(text)}</span>'


def variant_html(rec: dict, seq: str, codons_per_line: int,
                 mismatch: str | None) -> str:
    codons = codons_of(seq)
    aas = [aa_of(c) for c in codons]
    header = html_mod.escape(_variant_header(rec, len(codons), mismatch).lstrip("> "))
    parts = [f'<div class="variant"><div class="hdr'
             f'{" bad" if mismatch else ""}">{header}</div><pre>']
    for i in range(0, len(codons), codons_per_line):
        cs = codons[i:i + codons_per_line]
        as_ = aas[i:i + codons_per_line]
        ruler = html_mod.escape(_ruler(len(cs), i, offset=1, stride=4, width=3))
        codon_line = "|" + "|".join(_span(c, a) for c, a in zip(cs, as_)) + "|"
        aa_line = "|" + "|".join(_span(f"{a:^3}", a) for a in as_) + "|"
        parts.append(ruler + "\n" + codon_line + "\n" + aa_line + "\n")
    parts.append("</pre></div>")
    return "".join(parts)


def _legend_html() -> str:
    items = []
    seen = {}
    for aa, cls in AA_CLASS.items():
        seen.setdefault(cls, []).append(aa)
    labels = {"hydrophobic": "Hydrophobic (AVLIMFW)", "special": "Gly/Pro",
              "polar": "Polar (STNQCY)", "positive": "Positive (KRH)",
              "negative": "Negative (DE)", "stop": "Stop"}
    for cls, label in labels.items():
        items.append(f'<span style="color:{CLASS_COLOR[cls]}">&#9632;</span> '
                     f'{html_mod.escape(label)}')
    return " &nbsp; ".join(items)


def _hex_to_rgb(h: str) -> tuple[int, int, int]:
    h = h.lstrip("#")
    return int(h[0:2], 16), int(h[2:4], 16), int(h[4:6], 16)


def render_pdf(records: list[tuple[dict, str, str | None]], out_pdf: Path,
               codons_per_line: int, title: str) -> bool:
    """Write a colour-coded PDF mirroring the HTML view. Returns False (with a
    message) if fpdf2 is not installed."""
    try:
        from fpdf import FPDF
    except ImportError:
        print("PDF skipped: fpdf2 not installed. Install with: pip install fpdf2")
        return False

    grey = (110, 110, 110)
    black = (34, 34, 34)
    red = (198, 40, 40)
    pdf = FPDF(orientation="P", unit="pt", format="A4")
    pdf.set_auto_page_break(True, margin=28)
    pdf.set_margins(24, 28, 24)
    pdf.add_page()
    fs, lh = 7.5, 9.5  # font size / line height (pt)

    def latin1(s: str) -> str:
        # Core PDF fonts are Latin-1 only; swap common non-Latin-1 punctuation.
        return s.replace("—", "-").replace("→", "->").encode("latin-1", "replace").decode("latin-1")

    pdf.set_font("Courier", "B", 11)
    pdf.set_text_color(*black)
    pdf.multi_cell(0, 14, latin1(title))
    # legend: each class label printed in its own colour (no glyphs needed)
    pdf.set_font("Courier", "", 8)
    for cls, label in (("hydrophobic", "Hydrophobic AVLIMFW"), ("polar", "Polar STNQCY"),
                       ("positive", "Positive KRH"), ("negative", "Negative DE"),
                       ("special", "Gly/Pro"), ("stop", "Stop")):
        pdf.set_text_color(*_hex_to_rgb(CLASS_COLOR[cls]))
        pdf.write(12, f"{label}   ")
    pdf.ln(20)

    def seg(text: str, rgb: tuple[int, int, int]) -> None:
        pdf.set_text_color(*rgb)
        pdf.write(lh, text)

    for rec, seq, mism in records:
        codons = codons_of(seq)
        aas = [aa_of(c) for c in codons]
        pdf.set_font("Courier", "B", 8)
        pdf.set_text_color(*(red if mism else black))
        flag = f"   [!] {mism}" if mism else ""
        pdf.multi_cell(0, 11, latin1(_variant_header(rec, len(codons), mism).lstrip("> ") + flag))
        pdf.set_font("Courier", "", fs)
        for i in range(0, len(codons), codons_per_line):
            cs, as_ = codons[i:i + codons_per_line], aas[i:i + codons_per_line]
            seg(_ruler(len(cs), i, offset=1, stride=4, width=3) + "\n", grey)
            for c, a in zip(cs, as_):
                seg("|", grey)
                seg(c, _hex_to_rgb(CLASS_COLOR[AA_CLASS.get(a, "other")]))
            seg("|\n", grey)
            for a in as_:
                seg("|", grey)
                seg(f"{a:^3}", _hex_to_rgb(CLASS_COLOR[AA_CLASS.get(a, "other")]))
            seg("|\n", grey)
            pdf.ln(lh * 0.4)
        pdf.ln(6)

    pdf.output(str(out_pdf))
    return True


def build_html(title: str, blocks: list[str]) -> str:
    return f"""<!DOCTYPE html>
<html><head><meta charset="utf-8"><title>{html_mod.escape(title)}</title>
<style>
  body {{ font-family: sans-serif; margin: 24px; color: #222; }}
  h1 {{ font-size: 18px; }}
  .legend {{ margin: 8px 0 20px; font-size: 13px; }}
  .variant {{ margin: 0 0 18px; }}
  .hdr {{ font-weight: 600; font-size: 13px; margin-bottom: 2px; }}
  .hdr.bad {{ color: #c62828; }}
  pre {{ margin: 0; font-family: "DejaVu Sans Mono", Menlo, Consolas, monospace;
         font-size: 13px; line-height: 1.35; white-space: pre; }}
</style></head>
<body>
<h1>{html_mod.escape(title)}</h1>
<div class="legend">{_legend_html()}</div>
{''.join(blocks)}
</body></html>
"""


# --------------------------------------------------------------------------- #
# Driver                                                                       #
# --------------------------------------------------------------------------- #

def main(argv=None) -> int:
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--input", required=True, help="Nanobody mapping CSV.")
    p.add_argument("--seq-column", default="Coding_DNA",
                   help="In-frame coding-DNA column (default: Coding_DNA).")
    p.add_argument("--protein-column", default="Protein",
                   help="Protein column used to verify translation (default: Protein).")
    p.add_argument("--codons-per-line", type=int, default=20,
                   help="Codons per wrapped block (default: 20).")
    p.add_argument("--style", choices=["pipe", "space"], default="pipe",
                   help="Text layout: 'pipe' (default) or 'space'.")
    p.add_argument("--out-txt", type=Path, default=None,
                   help="Text output path (default: <input>_codons.txt).")
    p.add_argument("--out-html", type=Path, default=None,
                   help="HTML output path (default: <input>_codons.html).")
    p.add_argument("--out-pdf", type=Path, default=None,
                   help="PDF output path (default: <input>_codons.pdf). Needs fpdf2.")
    p.add_argument("--no-pdf", action="store_true",
                   help="Skip the PDF (text + HTML only).")
    p.add_argument("--limit", type=int, default=0,
                   help="Only render the first N variants (0 = all).")
    args = p.parse_args(argv)

    rows = list(csv.DictReader(open(args.input, newline="")))
    if not rows:
        p.error(f"No rows in {args.input}")
    if args.seq_column not in rows[0]:
        p.error(f"Column {args.seq_column!r} not found. Available: {', '.join(rows[0])}")
    if args.limit:
        rows = rows[:args.limit]

    stem = Path(args.input).with_suffix("")
    out_txt = args.out_txt or Path(f"{stem}_codons.txt")
    out_html = args.out_html or Path(f"{stem}_codons.html")
    out_pdf = args.out_pdf or Path(f"{stem}_codons.pdf")
    title = f"Codon view — {Path(args.input).name}"

    text_blocks, html_blocks, records, mismatches = [], [], [], 0
    for rec in rows:
        seq = (rec.get(args.seq_column) or "").upper()
        if not seq:
            continue
        mism = _translation_mismatch(seq, rec.get(args.protein_column))
        if mism:
            mismatches += 1
        text_blocks.append(variant_text(rec, seq, args.codons_per_line, args.style, mism))
        html_blocks.append(variant_html(rec, seq, args.codons_per_line, mism))
        records.append((rec, seq, mism))

    out_txt.write_text(f"{title}\n{'=' * len(title)}\n\n" + "\n".join(text_blocks))
    out_html.write_text(build_html(title, html_blocks))

    print(f"Rendered {len(text_blocks)} variants")
    print(f"Text : {out_txt}")
    print(f"HTML : {out_html}")
    if not args.no_pdf:
        if render_pdf(records, out_pdf, args.codons_per_line, title):
            print(f"PDF  : {out_pdf}")
    if mismatches:
        print(f"WARNING: {mismatches} variant(s) failed translation check (flagged [!]).")
    else:
        print("All variants translate cleanly to their Protein column.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
