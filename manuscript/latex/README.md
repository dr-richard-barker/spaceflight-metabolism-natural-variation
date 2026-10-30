# Manuscript (npj Microgravity style)

```bash
python3 scripts/18_manuscript_numbers.py      # numbers.tex, from results/tables
cd manuscript/latex && latexmk -pdf main.tex   # main.pdf
./build_docx.sh                                # main.docx (pandoc, author-date citations)
```

`numbers.tex` is generated: every number in the prose is a macro read from a results table, and
the comment beside each macro names the table. Edit the script, never the file.

`sn-jnl.cls` and `sn-nature.bst` are the Springer Nature template, copied already patched from
the first study's manuscript (`Photorespiration_multiomics_microgravity/manuscript/latex`, whose
README explains the two class patches TeX Live 2026 needs). `.latexmkrc` makes latexmk run BibTeX.

`references.bib`: every DOI entry was generated from its CrossRef record; run
`python3 scripts/check_references.py` to re-check titles, years and first authors against CrossRef.

## TODO before submission

Placeholders in `main.tex`; nothing was invented to fill them.

- [ ] Affiliation, city, country; co-authors if any; ORCID
- [ ] Author contributions, funding, competing interests
- [ ] Use of AI tools statement (the analysis code and first draft were prepared with an AI
      assistant; the journal's policy decides the wording)
- [ ] Zenodo DOI in Data availability
- [ ] Caldana et al. (2011) genotype: not stated in the abstract or the Wu supplement; the text
      does not call it Col-0 until the paper is checked
- [ ] Wu et al. (2018) stress condition: taken from the supplement ("32 C darkness"); confirm
      against the main text
