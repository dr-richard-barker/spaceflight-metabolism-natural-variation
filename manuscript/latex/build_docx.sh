#!/usr/bin/env bash
# Build main.docx from main.tex with pandoc, for co-authors without LaTeX. Run after latexmk.
# Word cannot display PDF figures, so the PNG copies from results/figures/ are swapped in;
# numbers.tex is inlined so the macros expand; citations are rendered by pandoc's citeproc
# (author-date), not by sn-nature.bst.
set -euo pipefail
cd "$(dirname "$0")"
# Pin the timestamps pandoc writes into the .docx, so an unchanged source rebuilds byte-identically.
export SOURCE_DATE_EPOCH=1735689600

tmp=$(mktemp -d)
trap 'rm -rf "$tmp"' EXIT
mkdir -p "$tmp/figures"
cp ../../results/figures/*.png "$tmp/figures/"
cp references.bib "$tmp/"
python3 - "$tmp" <<'PY'
import re, sys, pathlib
tmp = pathlib.Path(sys.argv[1])
src = pathlib.Path("main.tex").read_text()
src = src.replace("\\input{numbers}", pathlib.Path("numbers.tex").read_text())
src = re.sub(r"(\\includegraphics(\[[^]]*\])?\{[^}]*)\.pdf\}", r"\1.png}", src)
# sn-jnl's \abstract{...} is unknown to pandoc: turn it into the standard environment.
i = src.index("\\abstract{") + len("\\abstract{")
depth, j = 1, i
while depth:
    depth += {"{": 1, "}": -1}.get(src[j], 0)
    j += 1
src = src[:i - len("\\abstract{")] + "\\begin{abstract}" + src[i:j - 1] + "\\end{abstract}" + src[j:]
(tmp / "main.tex").write_text(src)
PY
(cd "$tmp" && pandoc main.tex --citeproc --bibliography=references.bib \
    --resource-path=.:figures -o main.docx)
cp "$tmp/main.docx" .
echo "wrote main.docx"
