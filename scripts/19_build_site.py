"""Build the GitHub Pages story page docs/index.html and its assets.

Fills scripts/templates/story.html with the values in manuscript/latex/numbers.tex (script 18, from
results/tables), so the page and the manuscript quote the same numbers. Copies the figures and the
manuscript PDF/DOCX into docs/, and writes docs/hub_thumb.png (320 x 180) for the AstroBotany hub card.
"""

import re
import shutil
from string import Template

from paths import FIGURES, ROOT

DOCS = ROOT / "docs"
TEMPLATE = ROOT / "scripts" / "templates" / "story.html"
NUMBERS = ROOT / "manuscript" / "latex" / "numbers.tex"
STEM = "spaceflight-metabolism-natural-variation"
SUP = str.maketrans("0123456789-", "⁰¹²³⁴⁵⁶⁷⁸⁹⁻")


def html_value(v):
    """LaTeX macro body -> text: '= 0.0281' -> '= 0.0281', '\\le 10^{-4}' -> '≤ 10⁻⁴', '7{,}790' -> '7,790'."""
    v = v.replace("{,}", ",").replace("\\le ", "≤ ")
    v = re.sub(r"\\times\s*10\^\{(-?\d+)\}", lambda m: "×10" + m.group(1).translate(SUP), v)
    v = re.sub(r"10\^\{(-?\d+)\}", lambda m: "10" + m.group(1).translate(SUP), v)
    return v.replace("-", "−") if re.fullmatch(r"-?[\d.]+", v) else v


def numbers():
    out = {}
    for line in NUMBERS.read_text().splitlines():
        m = re.match(r"\\newcommand\{\\(\w+)\}\{(.*)\}\s+%", line)
        if m:
            out[m.group(1)] = html_value(m.group(2))
    return out


def thumb():
    from PIL import Image
    im = Image.open(FIGURES / "fig04_geneset_tests.png").convert("RGB")
    w, h = im.size
    tw = min(w, int(h * 16 / 9))
    th = int(tw * 9 / 16)
    left, top = (w - tw) // 2, (h - th) // 2
    im.crop((left, top, left + tw, top + th)).resize((320, 180), Image.LANCZOS).save(DOCS / "hub_thumb.png")


def main():
    (DOCS / "figures").mkdir(parents=True, exist_ok=True)
    page = Template(TEMPLATE.read_text()).substitute(numbers())
    (DOCS / "index.html").write_text(page)
    for f in sorted(FIGURES.glob("fig0*.png")):
        shutil.copyfile(f, DOCS / "figures" / f.name)
    for ext in ("pdf", "docx"):
        src = ROOT / "manuscript" / "latex" / f"main.{ext}"
        if src.exists():
            shutil.copyfile(src, DOCS / f"{STEM}.{ext}")
    (DOCS / ".nojekyll").write_text("")
    thumb()
    print(f"  docs/index.html, {len(list((DOCS / 'figures').glob('*.png')))} figures, hub_thumb.png")


if __name__ == "__main__":
    main()
