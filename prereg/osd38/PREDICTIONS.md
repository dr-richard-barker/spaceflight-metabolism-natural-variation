# Pre-registered predictions: OSD-38 (BRIC-20) addendum

Written 2026-10-30 by `scripts/21_prereg_osd38.py`, after the OSD-38 inputs were fetched and protein coverage counted (`data/reference/osd38_protein_coverage.tsv`, which lists detected loci only), and **before any OSD-38 protein fold change or hardware-contrast statistic was computed**. The commit adding this folder precedes the commit adding the tests.

## Already known when this was written

- OSD-38 transcriptome, flight vs ground: the first study reported every one of these sets for the OSD-38 transcriptome by Mann-Whitney (its T11), and this study's T08 tested DIN and BCAA degradation there; rna_flight rows are re-tests, not blind.
- The OSDR processed proteome tables are truncated at 65,535 rows (see `scripts/20_osd38_inputs.py`).
- OSD-38 seedlings grew on 0.5x MS with 1 % sucrose (OSDR protocol), so carbon starvation in the classical sense is not expected in either arm.

## Layers

- `rna_flight`: OSD-38 RNA-seq, flight vs ground (both sealed BRIC-PDFU, dark) (primary)
- `membrane_flight`: OSD-38 membrane iTRAQ proteome, flight vs ground (primary)
- `soluble_flight`: OSD-38 soluble iTRAQ proteome, flight vs ground (truncated: top 437 groups) (secondary)
- `rna_hardware`: OSD-38 RNA-seq, ground: BRIC-PDFU vs closed container (both dark) (primary)

## Predictions

Flight layers: the first study's model predicts **no signal for any set** in dark flight. Hardware layer: hypoxia-responsive and fermentation sets **up** in BRIC-PDFU vs container, conditional on the BRIC holding less gas per seedling (the container volume is not reported); no prediction for the other sets (predicted none).

| Layer | Set | Set size | Measured proteins | Predicted |
|---|---|---|---|---|
| rna_flight | `din_first_study` | 5 |  | none |
| rna_flight | `bcaa_degradation` | 53 |  | none |
| rna_flight | `photosynthesis_apparatus` | 31 |  | none |
| rna_flight | `fermentation` | 6 |  | none |
| rna_flight | `hypoxia_responsive` | 3 |  | none |
| rna_flight | `ath00630` | 80 |  | none |
| rna_flight | `ath00710` | 70 |  | none |
| rna_flight | `ath00500` | 181 |  | none |
| rna_flight | `ath00010` | 120 |  | none |
| membrane_flight | `din_first_study` | 5 | 1 | none |
| membrane_flight | `bcaa_degradation` | 53 | 17 | none |
| membrane_flight | `photosynthesis_apparatus` | 31 | 0 | none |
| membrane_flight | `fermentation` | 6 | 1 | none |
| membrane_flight | `hypoxia_responsive` | 3 | 0 | none |
| membrane_flight | `ath00630` | 80 | 41 | none |
| membrane_flight | `ath00710` | 70 | 38 | none |
| membrane_flight | `ath00500` | 181 | 28 | none |
| membrane_flight | `ath00010` | 120 | 60 | none |
| soluble_flight | `din_first_study` | 5 | 0 | none |
| soluble_flight | `bcaa_degradation` | 53 | 7 | none |
| soluble_flight | `photosynthesis_apparatus` | 31 | 0 | none |
| soluble_flight | `fermentation` | 6 | 0 | none |
| soluble_flight | `hypoxia_responsive` | 3 | 0 | none |
| soluble_flight | `ath00630` | 80 | 30 | none |
| soluble_flight | `ath00710` | 70 | 33 | none |
| soluble_flight | `ath00500` | 181 | 13 | none |
| soluble_flight | `ath00010` | 120 | 36 | none |
| rna_hardware | `din_first_study` | 5 |  | none |
| rna_hardware | `bcaa_degradation` | 53 |  | none |
| rna_hardware | `photosynthesis_apparatus` | 31 |  | none |
| rna_hardware | `fermentation` | 6 |  | up |
| rna_hardware | `hypoxia_responsive` | 3 |  | up |
| rna_hardware | `ath00630` | 80 |  | none |
| rna_hardware | `ath00710` | 70 |  | none |
| rna_hardware | `ath00500` | 181 |  | none |
| rna_hardware | `ath00010` | 120 |  | none |

## Settings

- psm_filter: Score > 20 and all six iTRAQ reporter intensities > 0
- loci: AGI locus from Sequence Id; loci with at least 2 distinct peptides (the paper's rule)
- soluble_rows: protein groups 1-437 (the table is cut inside group 438)
- membrane_rows: block 1, rows 1-48,799 (groups 1-2,948); block 2 excluded
- normalisation: log2; each channel centred on its median; each PSM centred on its own mean
- summary: protein value per channel = median over its PSMs
- contrast: log2FC = mean(flight S1-S3) - mean(ground GC1-GC3)
- test: two-sided Student's t, 3 vs 3 (the paper's test)
- relevant_for_paintomics: p <= 0.05 and |log2FC| > 0.2 (the paper's rule)
- Set statistic: signed: mean gene/protein statistic in the set minus the mean over all measured features.
- Null: genomic rotation over the measured features in prereg/gene_order.tsv.gz order, 10000 offsets, seed 20261030.
- Minimum measured members: 3.
- Multiplicity: Benjamini-Hochberg over every testable primary row (rna_flight, membrane_flight, rna_hardware); soluble_flight rows are reported with their own BH and are not used for verdicts.
- Verdict: observed = up/down when q < 0.05, else none; MATCH when observed equals the prediction.
- Concordance: descriptive: Spearman correlation of protein and transcript log2FC over shared loci, OSD-38 (membrane, soluble) and OSD-522 (the first study's proteome, UniProt mapped to AGI by KEGG conv/uniprot/ath).
- PaintOmics: job D: OSD-38 genes (log2FC; relevant = adjusted p < 0.05) + membrane proteins (relevant = p <= 0.05 and |log2FC| > 0.2) + the reference metabolite layer, same settings as job A; compared descriptively with job A (OSD-522, lit BRIC-LED).

## Frozen files

| File | SHA-256 |
|---|---|
| `gene_sets.tsv` | `743c02179c8d7d3bcb3ecf0ab53fa9d7903eb2b4d1b7618cbbd0396817857efd` |
| `settings.json` | `e735ddf2d6f8a749d18f11ce4afc8d4bf7cfc2e90f4c91979b4a1a66b9e8f640` |
| `PREDICTIONS.tsv` | `84dff2bb1f69398821e4bf844fbc34dcaad3a02753c18acd92829726e342a939` |
