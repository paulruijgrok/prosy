"""Unit tests for scripts/nanobodies/codon_view.py.
Run: python -m pytest tests/ (or python tests/test_codon_view.py).
"""

import sys
from pathlib import Path

_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(_ROOT))
sys.path.insert(0, str(_ROOT / "scripts" / "nanobodies"))

import codon_view as cv  # noqa: E402
from prosy.core.sequence import CODON_TABLE, translate  # noqa: E402

CDS = "CAGGTGCAGCTGGTGGAAAGC"  # Q V Q L V E S


def test_codons_split_and_translate_roundtrip():
    codons = cv.codons_of(CDS)
    assert codons == ["CAG", "GTG", "CAG", "CTG", "GTG", "GAA", "AGC"]
    assert "".join(cv.aa_of(c) for c in codons) == translate(CDS)
    assert cv.aa_of("CAG") == CODON_TABLE["CAG"] == "Q"


def test_codons_ignore_trailing_partial():
    assert cv.codons_of("CAGGT") == ["CAG"]  # trailing 2 nt dropped


def test_text_block_lines_align():
    codons = cv.codons_of(CDS)
    aas = [cv.aa_of(c) for c in codons]
    ruler, codon_line, aa_line = cv.render_text_block(codons, aas, 0, "pipe")
    assert len(codon_line) == len(aa_line)          # codon/AA lines aligned
    assert codon_line.startswith("|") and codon_line.count("|") == len(codons) + 1
    assert " Q " in aa_line and "CAG" in codon_line
    assert ruler.lstrip().startswith("1")


def test_space_style_alignment():
    codons = cv.codons_of(CDS)
    aas = [cv.aa_of(c) for c in codons]
    ruler, codon_line, aa_line = cv.render_text_block(codons, aas, 0, "space")
    assert len(codon_line) == len(aa_line)
    assert "|" not in codon_line


def test_translation_mismatch_detection():
    assert cv._translation_mismatch(CDS, "QVQLVES") is None
    assert cv._translation_mismatch(CDS, "QVQLVES*") is None  # trailing stop ok
    assert cv._translation_mismatch(CDS, "QVQLVES", ) is None
    assert cv._translation_mismatch(CDS, "MMMMMMM") is not None


def test_variant_text_has_header_and_blocks():
    rec = {"NbID_parent": "Nb01_T27F", "NbID": "Nb73", "Well": "B1", "Protein": translate(CDS)}
    out = cv.variant_text(rec, CDS, codons_per_line=4, style="pipe",
                          mismatch=None)
    assert ">>> Nb01_T27F  Nb73  B1" in out
    assert "(7 codons)" in out
    assert "|CAG|GTG|CAG|CTG|" in out   # first block, 4 codons per line


def test_variant_html_is_colored():
    rec = {"NbID": "Nb01", "Protein": translate(CDS)}
    frag = cv.variant_html(rec, CDS, codons_per_line=20, mismatch=None)
    assert "<span" in frag and "CAG" in frag
    assert 'class="variant"' in frag


def test_render_pdf_optional(tmp_path):
    # PDF needs fpdf2; when present it writes a non-empty file, else returns False.
    records = [({"NbID": "Nb01", "Protein": translate(CDS)}, CDS, None)]
    out = tmp_path / "codons.pdf"
    ok = cv.render_pdf(records, out, codons_per_line=20, title="Codon view — test")
    if ok:
        assert out.exists() and out.stat().st_size > 0
    else:
        assert not out.exists()


def _run_all():
    fns = [v for k, v in sorted(globals().items()) if k.startswith("test_") and callable(v)]
    for fn in fns:
        fn()
        print(f"  ok  {fn.__name__}")
    print(f"\n{len(fns)} tests passed")


if __name__ == "__main__":
    _run_all()
