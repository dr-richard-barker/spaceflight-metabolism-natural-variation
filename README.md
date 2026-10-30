# Natural variation in *Arabidopsis* light/dark metabolism tests a spaceflight carbon-starvation model

**A gas-transport model of the BRIC-LED-001 flight (OSD-522) predicted that microgravity starves
lit plants of carbon, so lit spaceflight transcriptomes should carry a carbon-starvation signature
and dark ones should not. Measured light/dark metabolism from hundreds of natural accessions
confirms the model's metabolic picture in outline and corrects it in detail. Pre-registered
tests on six flight contrasts find the starvation signal strongest in the dark, where the model
predicts none.**

[**Story page**](https://dr-richard-barker.github.io/spaceflight-metabolism-natural-variation/)
· [manuscript PDF](manuscript/latex/main.pdf) · [manuscript DOCX](manuscript/latex/main.docx)
· [figures](results/figures) · [tables](results/tables) · [pre-registration](prereg/PREDICTIONS.md)

![Pre-registered light/dark gene-set tests](results/figures/fig04_geneset_tests.png)

## What is new

The first study ([Photorespiration_multiomics_microgravity](https://github.com/dr-richard-barker/Photorespiration_multiomics_microgravity),
read here at commit 515226d) predicted the OSD-522 metabolome from a model, because NASA's Open
Science Data Repository holds no plant metabolome. This study tests those predictions, and the
model's light/dark claim, with measured natural-variation data:

- **Zhu et al. 2024**: 259 accessions with leaf metabolomes before and after 6 days of darkness
- **Wu et al. 2018**: 301 accessions with paired control and darkness-stress (32 °C) values, and
  their metabolite GWAS loci per condition
- **Caldana et al. 2011**: a light/dark time course (from the Wu supplement)
- **MetaboLights**: Col-0 control groups from six GC-MS studies, for replicate variance

## Findings

1. **The model's metabolite directions mostly hold, its grouping does not.** Of the 6 predicted
   compounds that measured darkness can test, 5 agree (maltose down; asparagine, valine, leucine,
   isoleucine up). Glutamate, which the model raises with them, falls; glutamine does not move;
   leucine responds 2.7–2.8 times as strongly as valine or asparagine. 41 robust darkness
   responses were not predicted, including malate and 2-oxoglutarate (both fall).
   [T03](results/tables/T03_model_vs_reference.tsv), [fig01](results/figures/fig01_model_vs_reference.png)
2. **A measured metabolite layer leaves PaintOmics unchanged.** Swapping the 21 predicted
   compounds for 96 measured ones kept the same 6 significant pathways (rank correlation 0.95),
   and raw hub p-values were identical for every shared compound: the gene layer sets the result.
   [T04–T07](results/tables), [fig02](results/figures/fig02_jobA_vs_job1.png)
3. **The starvation signal is not light-specific.** 18 of 36 pre-registered outcomes matched,
   but no lit contrast gave a significant starvation signal, while the starvation (DIN) and
   branched-chain amino acid catabolism sets moved strongly in three dark contrasts, down in
   BRIC-PDFU seedlings and up in the leaves of dark-grown CARA seedlings. GWAS-locus gene sets showed nothing
   light-specific. [T08](results/tables/T08_geneset_tests.tsv), [fig04](results/figures/fig04_geneset_tests.png)
4. **The predicted metabolite shifts are too small for a feasible panel.** With 6 plants per arm,
   GC-MS detects about 1.48-fold changes at best; the model's 0.13 log₂ needs 92–1,342 plants per
   arm. [T10](results/tables/T10_gcms_power.tsv), [fig05](results/figures/fig05_gcms_power.png)

Labels used throughout: the PaintOmics metabolite layer here is a **ground carbon-starvation
reference, not flight measurements**.

## Pre-registration

`prereg/` (gene sets, genome order, settings, predicted outcomes) was committed in 858ca21,
before any gene-set statistic was computed on a flight contrast; the tests came in 0d07084.
`scripts/16_gwas_light_dark_tests.py` reads `prereg/` and will not run without it. The DIN set and
the compound comparison are flagged there as re-tests and consistency checks, since the first
study's results and the reference had been seen.

## Layout

| Path | Contents |
|---|---|
| `data/atlas/`, `data/photoresp/` | inputs read with `git show` from the two sibling repositories (`PROVENANCE.md`: commit, path, SHA-256) |
| `data/osdr/` | per-contrast GeneLab statistics for OSD-38, OSD-321, OSD-678 light and dark, OSD-427 |
| `data/reference/` | KEGG name map and reviewed overrides, MetaboLights sources and derived SDs |
| `data/paintomics_upload/` | the three PaintOmics bundles and the layer key |
| `prereg/` | the frozen pre-registration |
| `scripts/` | numbered analysis scripts; `run_all.sh` runs them in order |
| `results/tables/`, `results/figures/` | T01–T10 (see `MANIFEST.tsv`), fig01–fig05 |
| `results/paintomics_raw/` | recovered PaintOmics job records (guest jobs expire after 7 days) |
| `manuscript/latex/` | manuscript source, generated `numbers.tex`, PDF and DOCX |
| `docs/` | the GitHub Pages story page |

## Reproduce

```bash
./run_all.sh          # every table, figure and manuscript number, byte-identical
./run_all.sh --pdf    # also the manuscript PDF and DOCX
python3 scripts/check_references.py   # every DOI against CrossRef
```

Python 3.9 with pandas, numpy, scipy, matplotlib and PyDESeq2 0.4.12. Raw third-party files
(the Wu supplement, KEGG lists, MetaboLights study files) are fetched into the git-ignored
`data/cache/`; only derived values are committed. PaintOmics jobs are not re-submitted: their
records are vendored.

## Data sources and licences

- NASA OSDR: OSD-522, OSD-38, OSD-321, OSD-678, OSD-427 (public domain, CC0)
- Zhu et al. 2024 (*Sci Data*, doi:10.1038/s41597-024-03694-2) and Wu et al. 2018 (*Mol Plant*,
  doi:10.1016/j.molp.2017.08.012), as derived tables in the
  [Arabidopsis Atlas](https://github.com/dr-richard-barker/arabidopsis-atlas) (commit 69b21c6)
- Caldana et al. 2011 (*Plant J*, doi:10.1111/j.1365-313x.2011.04640.x), via Wu et al.'s
  supplementary tables 7–8; only derived ratios are committed
- MetaboLights MTBLS22, 41, 42, 363, 2289, 4066 (EMBL-EBI terms of use); only derived SDs are
  committed
- KEGG REST API (compound names, ChEBI conversion, gene lists, pathway membership)

Code: MIT (`LICENSE`).

## Open items

- Zenodo deposit (manual)
- Manuscript authorship, affiliation, funding and AI-use statements (`manuscript/latex/README.md`)
- Caldana genotype and the Wu stress condition, to confirm against the papers' main text
