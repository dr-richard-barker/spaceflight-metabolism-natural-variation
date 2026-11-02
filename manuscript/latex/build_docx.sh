#!/usr/bin/env bash
# Build main.docx and supplementary.docx from the LaTeX sources with pandoc, for co-authors without
# LaTeX. Run after latexmk has built main.pdf (the SI's cross-references are read from main.aux).
# Word cannot display PDF figures, so the PNG copies from results/figures/ are swapped in;
# numbers.tex and supp_tables.tex are inlined so the macros and tables expand; citations are
# rendered by pandoc's citeproc (author-date), not by sn-nature.bst.
set -euo pipefail
cd "$(dirname "$0")"
# Pin the timestamps pandoc writes into the .docx, so an unchanged source rebuilds byte-identically.
export SOURCE_DATE_EPOCH=1735689600

tmp=$(mktemp -d)
trap 'rm -rf "$tmp"' EXIT
mkdir -p "$tmp/figures"
cp ../../results/figures/*.png "$tmp/figures/"
cp references.bib "$tmp/"
for f in main supplementary; do
python3 - "$tmp" "$f" <<'PY'
import re, sys, pathlib
tmp, name = pathlib.Path(sys.argv[1]), sys.argv[2]
src = pathlib.Path(f"{name}.tex").read_text()
src = src.replace("\\input{numbers}", pathlib.Path("numbers.tex").read_text())
if name == "supplementary":
    src = src.replace("\\input{supp_tables}", pathlib.Path("supp_tables.tex").read_text())
    # xr references into the manuscript, which pandoc cannot follow: resolve them from main.aux.
    labels = dict(re.findall(r"\\newlabel\{([^}]+)\}\{\{([^}]*)\}", pathlib.Path("main.aux").read_text()))
    src = re.sub(r"\\ref\{([^}]+)\}", lambda m: labels.get(m.group(1), m.group(0)), src)
src = re.sub(r"(\\includegraphics(\[[^]]*\])?\{[^}]*)\.pdf\}", r"\1.png}", src)
# pandoc drops tables whose column spec uses a \newcolumntype; Word wraps text anyway.
src = src.replace("\\newcolumntype{L}{>{\\raggedright\\arraybackslash}X}", "")
src = re.sub(r"\\begin\{tabularx\}\{\\textwidth\}\{@\{\}l(L+)@\{\}\}",
             lambda m: "\\begin{tabularx}{\\textwidth}{@{}l" + "X" * len(m.group(1)) + "@{}}", src)
# sn-jnl's \abstract{...} is unknown to pandoc: turn it into the standard environment.
if "\\abstract{" in src:
    i = src.index("\\abstract{") + len("\\abstract{")
    depth, j = 1, i
    while depth:
        depth += {"{": 1, "}": -1}.get(src[j], 0)
        j += 1
    src = src[:i - len("\\abstract{")] + "\\begin{abstract}" + src[i:j - 1] + "\\end{abstract}" + src[j:]
(tmp / f"{name}.tex").write_text(src)
PY
(cd "$tmp" && pandoc "$f.tex" --citeproc --bibliography=references.bib \
    --resource-path=.:figures -o "$f.docx")
cp "$tmp/$f.docx" .
done
echo "wrote main.docx and supplementary.docx"
