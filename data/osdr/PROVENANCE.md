# Per-contrast gene statistics (NASA OSDR, CC0)

Extracted by `scripts/00_vendor_inputs.py` from the GeneLab differential-expression tables the first study downloaded (its `07_hardware_ladder.py`, URL pattern
`https://osdr.nasa.gov/geode-py/ws/studies/{osd}/download?source=datamanager&file={f}`). Each file keeps TAIR, the matched flight-vs-ground `log2fc`, `stat` (Wald), `pvalue` and `padj` columns named in `data/photoresp/study_registry.tsv`, plus the normalised counts of the matched samples where their names state the arm. OSD-522 is not here: it is computed from raw counts (data/atlas/).

| File | Contrast | GeneLab file | Genes | Samples kept | Source SHA-256 (first 16) |
|---|---|---|---|---|---|
| `contrast_OSD-38.csv.gz` | OSD-38 | `GLDS-38_rna_seq_differential_expression.csv` | 27723 | flight 3 / ground 3 | `64cee38cf11e622f` |
| `contrast_OSD-321.csv.gz` | OSD-321 | `GLDS-321_rna_seq_differential_expression_rRNArm_GLbulkRNAseq.csv` | 25505 | none | `adaf833780b85aa0` |
| `contrast_OSD-678-light.csv.gz` | OSD-678-light | `GLDS-612_rna_seq_differential_expression_rRNArm_GLbulkRNAseq.csv` | 24474 | flight 3 / ground 3 | `52d2fb0aae764db4` |
| `contrast_OSD-678-dark.csv.gz` | OSD-678-dark | `GLDS-612_rna_seq_differential_expression_rRNArm_GLbulkRNAseq.csv` | 24474 | flight 3 / ground 3 | `52d2fb0aae764db4` |
| `contrast_OSD-427.csv.gz` | OSD-427 | `GLDS-427_rna_seq_differential_expression.csv` | 24863 | none | `1d011c49a8ce8e0d` |

## OSD-38 ground contrasts

From the same GeneLab table as `contrast_OSD-38.csv.gz` (the first study's cached `dge_OSD-38.csv`), by `scripts/20_osd38_inputs.py`. GC = BRIC-PDFU ground control; RNAlat = plates in closed containers at Ohio University, fixed in RNAlater like GC (Basu et al. 2017).

| File | Contrast | Genes | Source SHA-256 (first 16) |
|---|---|---|---|
| `contrast_OSD-38-hardware.csv.gz` | BRIC-PDFU (ground) vs closed container: `(GC)v(RNAlat)` columns | 27723 | `64cee38cf11e622f` |
