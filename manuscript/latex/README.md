# Manuscript (npj Microgravity style)

```bash
python3 scripts/18_manuscript_numbers.py      # numbers.tex, from results/tables
python3 scripts/29_supplementary_tables.py     # supp_tables.tex (Supplementary Tables S2-S14)
cd manuscript/latex && latexmk -pdf main.tex   # main.pdf
latexmk -pdf supplementary.tex                 # supplementary.pdf (reads main.aux for cross-references)
./build_docx.sh                                # main.docx, supplementary.docx (pandoc, author-date citations)
```

The response to the reviewer (`manuscript/response/`, built by its `build_response.py` from the same
`numbers.tex`) quotes confidential review comments and is git-ignored; it lives only in the working copy.

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
- [ ] Author contributions, acknowledgements with funding (npj Microgravity allows no separate
      funding statement), competing interests
- [ ] Use of AI tools statement, now a placeholder in Methods as npj Microgravity requires (the
      analysis code and first draft were prepared with an AI
      assistant; the journal's policy decides the wording)
- [ ] Zenodo DOI in Data availability
- [x] Caldana et al. (2011) genotype: Col-0 per the authors' ArrayExpress deposit E-MTAB-375 and
      Wu et al. (accepted manuscript, lines 476-478); the Caldana paper itself was not accessible
- [x] Wu et al. (2018) stress condition: darkness, 32 C, 75 % humidity, 1,280 min (accepted
      manuscript, lines 512-515)
