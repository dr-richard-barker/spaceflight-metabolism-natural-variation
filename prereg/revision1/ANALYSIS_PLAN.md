# Revision 1: analysis plan for the reviewer-requested analyses

Written 2026-11-02, after the reviewer's comments and **before any of the analyses below was run**.
These analyses are **post hoc**: they answer a reviewer and were not part of the pre-registrations
in `prereg/` (commit 858ca21) or `prereg/osd38/` (commit eb005ec). The commit that adds this folder
precedes the commits that add their scripts and outputs. Every setting is in `settings.json`.

What stays fixed: the pre-registered reference (T02), the pre-registered test outputs (T08, T12) and
the T03 verdict column must regenerate byte-identically. New analyses add tables and columns; they
do not overwrite pre-registered results.

## R1. The darkness reference (reviewer points 1-7)
- **Strata.** The pre-registered consensus rule is applied within each stratum: all three sources
  (the pre-registered reference), dark only (Zhu, Caldana) and dark+heat (Wu: 32 C darkness).
  Each compound gets a per-source direction and qualification, so the table shows where the
  32 C source alone drives a call.
- **Sensitivity.** Accession agreement at 60, 70, 75, 80 and 90 %, crossed with the Caldana window
  at 60, 180 and 360 min or later. Fixed: n >= 20 accessions, 90 % for a robust single source,
  Caldana p < 0.05. Reported: category counts; whether Leu, Ile and Val stay up and malate and
  2-oxoglutarate stay down; Glu and Gln; the T03 verdict counts; and how many robust responses
  the model did not predict.
- **Auditability.** KEGG mapping counts (0, 1 or more candidates; overrides; exclusions; mapped
  fraction by Wu class); coverage (model predictions scorable and comparable; robust responses
  the model predicted); and a `verdict_detail` column in T03. That column splits "untested" into
  *not measured* and *measured, no directional response*, without changing the verdict column.

## R2. Gene-set diagnostics (points 11, 12)
- **Rotation calibration.** On each of the six flight contrasts: 200 random gene sets and 200
  contiguous genomic blocks at sizes 5, 50, 300, 500 and 2000. Each is tested with 999 rotation
  offsets, and the blocks also against a naive random-gene null. Reported: empirical type I error
  at 0.05 and QQ plots. Seed 20261102.
- **Exact sample permutation.** OSD-38, OSD-678 light and OSD-678 dark, 3 vs 3, all 20 labellings.
  The p-values are coarse (minimum 0.1) and are reported as such.
- **Pipeline comparability.** PyDESeq2 on OSD-38 raw counts against GeneLab's Wald statistics.
- **DIN.** Per-gene statistics in each contrast, and a leave-one-out set test.

## R3. OSD-678 flight x light interaction (point 9)
GLDS-612 rRNA-removed RSEM counts, 12 Col-0 samples. edgeR TMM and limma-voom, design ~ flight * light.
Rotation tests on the interaction moderated t for eleven sets, with BH over the family. The model's
expectation is that the starvation sets (DIN, BCAA degradation) have a positive interaction (a larger
flight effect in light); it makes no claim for the others.

## R4. Hypoxia positive control (point 15)
OSD-136 (GeneLab, low O2 at 50 kPa vs O2 restored; Ws-0 roots) and Mustroph et al. 2009 GSE14502
(2 h hypoxia vs control, total mRNA, shoot and root; limma on the RMA series matrix). The sets tested
are hypoxia-responsive, fermentation and DIN, with the rotation null and BH over the family.
**Pass criterion:** a set is hypoxia-sensitive if it is up at q < 0.05 in at least one GSE14502
contrast. The BRIC-vs-container result for a set counts as informative only if that set passes.

## R5. OSD-38 proteome (point 14)
- limma moderated t on the pre-registered protein quantification, membrane fraction, with rotation
  tests for BCAA degradation and the four KEGG carbon pathways.
- A partitioning check on proteins quantified in both fractions.
- A query of the OSDR file versions.
- Raw iTRAQ reprocessing is not attempted in this revision.

## R6. PaintOmics controls (points 8, 9, 16)
Jobs E-K, as defined in `settings.json`:
- metabolite-layer ablations: job A with no metabolites, and with effect-scaled values;
- layer ablations of the hubs: genes only and proteins only, for jobs A and D;
- two DEG-count-matched versions of job B, downsampled to job C's relevant-gene count.

Comparisons: pathways, rank correlation and hub sets. The hub-overlap null is hypergeometric,
plus a KEGG-degree test and a degree-matched permutation.

## R7. Power (point 13)
Rename the 1.5-fold column; show Bonferroni over the 29-compound panel; add a provenance table for
the MetaboLights studies.
