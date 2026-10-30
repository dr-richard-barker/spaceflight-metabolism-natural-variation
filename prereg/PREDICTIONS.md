# Pre-registered predictions

Written 2026-10-30 by `scripts/11_prereg.py`, after the light/dark reference (T01, T02) and the gene sets were built and **before any gene-set statistic was computed on any flight contrast in this repository**. The commit that adds this file precedes the commit that adds the tests; `git log` shows the order.

## What was already known when this was written

- The first study's results: its T07 (OSD-522 transcriptome and proteome against KEGG and marker sets) and T11 (the same sets in all six contrasts), including the DIN set. The `din_first_study` rows below are therefore **re-tests with stronger nulls, not blind predictions**.
- T02, the light/dark metabolite reference. The compound comparison (T03) is a consistency check, not a blind test.
- Nothing about `bcaa_degradation` or the four Wu sets in any flight contrast.

## The model's light/dark prediction

The first study's model puts a microgravity boundary-layer deficit on net assimilation whenever the plants photosynthesise (7.3-8.8 % for a rosette at 1000 umol m-2 s-1; 0.4-2.4 % at OSD-522's 60 umol m-2 s-1, its T13). In darkness there is no assimilation to lose. So: a carbon-starvation signal in **lit** flight contrasts, **none** in dark ones.

- Lit: OSD-522, OSD-678-light, OSD-427
- Dark: OSD-38, OSD-321, OSD-678-dark

## Gene sets (frozen in `gene_sets.tsv`)

| Set | Genes | Definition | Statistic | Lit | Dark |
|---|---|---|---|---|---|
| `din_first_study` | 5 | First study's DIN set (KEGG description regex \bDIN\d+;) | signed | up | none |
| `bcaa_degradation` | 53 | KEGG ath00280 valine, leucine and isoleucine degradation | signed | up | none |
| `wu_stress_only` | 491 | Wu 2018 GWAS regions found only under stress (32 C darkness) | unsigned | enriched | none |
| `wu_dark_metabolite` | 540 | Wu 2018 loci of metabolites with a directional light/dark reference | unsigned | enriched | none |
| `wu_control_only` | 301 | Wu 2018 GWAS regions found only under control (light) | unsigned | none | none |
| `wu_shared` | 1954 | Wu 2018 GWAS regions found under both conditions | unsigned | none | none |

Basis for each prediction:

- `din_first_study`: model: assimilation deficit in light only; RE-TEST, the first study already reported this set in all six contrasts (its T07, T11).
- `bcaa_degradation`: reference: branched-chain amino acids are the strongest darkness response in all three sources; their catabolic route should rise with a carbon deficit.
- `wu_stress_only`: loci that shape metabolism only under 32 C darkness should be perturbed if lit flight is starvation-like; weak prior (a locus gene need not be responsive).
- `wu_dark_metabolite`: loci of metabolites with a directional darkness response; weak prior.
- `wu_control_only`: comparator: control-condition loci, about which the model makes no claim.
- `wu_shared`: comparator: loci found in both conditions, about which the model makes no claim.

## Tests

- Gene statistic: GeneLab Wald stat (Stat_ column) for the five vendored contrasts; for OSD-522 the PyDESeq2 Wald stat recomputed from raw counts as the first study did (6 flight vs 6 ground shoots, genes with total count >= 10).
- Signed set statistic: mean gene statistic in the set minus the mean over all measured genes. Unsigned: mean |gene statistic| in the set minus the mean |statistic| over all measured genes.
- Primary null, OSD-522: sample permutation: all 924 ways to label 6 of the 12 samples as flight; gene statistic a Welch t on log2(median-of-ratios normalised count + 1).
- Primary null, the other five: genomic rotation: the set's gene positions shifted together by a random offset along the genome (KEGG order, prereg/gene_order.tsv.gz), preserving set size and clustering; 10000 offsets. Seed 20261030.
- p-values: two-sided (1 + #{|null| >= |obs|}) / (B + 1) for rotation; #{|null| >= |obs|} / 924 for the sample permutation, which includes the observed labelling; one-sided same with null >= obs, for the unsigned statistic.
- Multiplicity: Benjamini-Hochberg across the 36 primary tests (6 sets x 6 contrasts), q < 0.05.
- Verdict: observed = up/down (signed) or enriched (unsigned) when BH q < 0.05, else none; MATCH when observed equals the prediction.
- Secondary, reported only: rotation null for OSD-522 as well; two-sided Mann-Whitney of set log2FC against all other genes (the first study's test); the non-primary statistic. Reported, not used for verdicts.

## Compound predictions against the reference (T03)

- Scored: predictions whose driver is A or starvation. a T06 prediction whose driver is net assimilation (A) or carbon starvation is compared with the T02 consensus direction: agrees if the signs match, disagrees if they differ, untested if T02 has no directional status for that KEGG ID.
- Not comparable: Vo- and RuBP-driven predictions: darkness stops oxygenation and carboxylation together, so a darkness reference cannot test them.
- T02 was built and inspected on 2026-10-30 before this file was written, so this comparison is a consistency check, not a blind test.

## PaintOmics with the data-derived layer

- Layer: all T02 compounds; value = direction x 0.13052 (the model's carbon-deficit magnitude, T06) for robust and single, 0 for none and conflicting; relevant = robust only. Labelled everywhere: *ground carbon-starvation reference, not flight measurements*.
- Job A: OSD-522 gene and protein files byte-identical to job 1 (SHA-256 checked) + layer. Job B: OSD-678-light genes (GeneLab log2FC; relevant = adjusted p < 0.05) + layer. Job C: OSD-678-dark genes, same rule + layer. Databases KEGG + MapMan, AI off.
- Comparison: descriptive: pathways at combined-p BH FDR < 0.05 gained or lost, rank shifts, hubs at FDR < 0.05.
- Mechanistic check: hub raw p-values for compounds relevant in both job 1 and job A are expected to be identical, because hub significance depends on the gene layer and the KEGG network, not on the metabolite values; FDR may differ because more compounds enter the family.

## GC-MS panel power

- Variance: Col-0 wild-type control groups of MetaboLights GC-MS leaf or rosette studies; within-group SD of log2 intensity per metabolite, median across groups.
- Test: two-sided two-sample t, equal n per arm, noncentral t; alpha 0.05, power 0.8.
- Effects (log2): model_Q60_rosette_max = log2(1 - 0.0237); model_Q1000_T06 = 0.13052; half_fold = 0.585; two_fold = 1.0.
- Multiplicity: unadjusted and Bonferroni across the panel; BH lies between them.

## Frozen files

| File | SHA-256 |
|---|---|
| `gene_sets.tsv` | `6fabfa17ec7f96a88a09a3f24961c4c7c31a9af09e06d5dd6dcfa5929fa1f0d0` |
| `gene_order.tsv.gz` | `c64f61b3ddfcf1b99ab3e7f2b753731d8378e2fd3e2492d7e83c0c1328a7cb77` |
| `settings.json` | `350d6ba095bea1de40be9fc7546519b3c612fe15a8a7b8f7d6347870b150d664` |
| `PREDICTIONS.tsv` | `bfb2af94f99ef6b407ac9a0db349fb0fca677f2b7234dfa6e1e0d4c94ce6973e` |
