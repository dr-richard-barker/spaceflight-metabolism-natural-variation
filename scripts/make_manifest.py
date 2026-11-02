"""Write results/tables/MANIFEST.tsv, one row per derived table, and fail on any undescribed table.

Adapted from the first study's scripts/make_manifest.py (MIT, same author): a new output cannot
appear without saying what it holds and which script made it.
"""

import sys

from paths import TABLES

DESCRIBED = {
    "T01_light_dark_sources.tsv": ("scripts/10_light_dark_reference.py",
        "Per-source light/dark response per compound (Zhu darkness, Wu stress, Caldana dark vs light): "
        "effect, fraction of accessions moving the same way, spread, p, KEGG ID, qualifies",
        "data/atlas (Zhu, Wu); Wu 2018 supplementary tables 7-8 (Caldana time course)"),
    "T02_light_dark_reference.tsv": ("scripts/10_light_dark_reference.py",
        "Consensus light/dark reference per KEGG compound: status, direction, standardised effect",
        "T01"),
    "T03_model_vs_reference.tsv": ("scripts/12_model_vs_reference.py",
        "The first study's 21 predicted compounds against the reference (pre-registered rule), plus "
        "directional reference compounds the model did not predict",
        "data/photoresp/T06_predicted_compounds.tsv, T02"),
    "T04_paintomics_pathways.tsv": ("scripts/15_paintomics_tables.py",
        "All 232 pathways in PaintOmics jobs 1, A, B, C: combined and per-omic p, BH q, rank",
        "data/photoresp/job1_*, results/paintomics_raw/job{A,B,C}_*"),
    "T05_paintomics_hubs.tsv": ("scripts/15_paintomics_tables.py",
        "Metabolite hub analysis per compound and job: best adjusted p, raw p, DE neighbours, step",
        "as T04"),
    "T06_layer_matching.tsv": ("scripts/15_paintomics_tables.py",
        "KEGG IDs PaintOmics used for the reference layer in jobs A-C against the intended IDs",
        "as T04; data/paintomics_upload/metabolomics_key.tsv"),
    "T07_paintomics_summary.tsv": ("scripts/15_paintomics_tables.py",
        "Job-level counts, gained and lost pathways and hubs, the pre-registered hub check",
        "T04, T05"),
    "T08_geneset_tests.tsv": ("scripts/16_gwas_light_dark_tests.py",
        "36 pre-registered gene-set tests: prediction, observed, verdict, primary p and q, z, and "
        "secondary nulls (rotation, permutation, Mann-Whitney)",
        "prereg/, data/osdr, data/atlas GLDS-522 counts"),
    "T09_geneset_summary.tsv": ("scripts/16_gwas_light_dark_tests.py",
        "MATCH counts per gene set in lit and dark contrasts", "T08"),
    "T11_osd38_proteome.tsv": ("scripts/22_osd38_tests.py",
        "OSD-38 protein log2FC, flight vs ground, per iTRAQ fraction (membrane block 1; soluble groups 1-437 "
        "of the truncated OSDR table): t, p, PSMs, relevant at the paper's rule",
        "OSDR GLDS-38 WyattLab soluble/membrane CSVs (data/reference/osd38_sources.tsv)"),
    "T12_osd38_geneset_tests.tsv": ("scripts/22_osd38_tests.py",
        "Pre-registered OSD-38 gene-set tests: nine sets x four layers (flight transcripts, membrane and "
        "soluble proteins; ground BRIC vs container transcripts), prediction, observed, verdict, p, q, z",
        "prereg/osd38/, T11, data/osdr/contrast_OSD-38*.csv.gz"),
    "T13_osd38_concordance.tsv": ("scripts/22_osd38_tests.py",
        "Transcript-protein concordance (OSD-38, OSD-522), the paper's DE-count check, abundance-matched "
        "robustness of the membrane sets, and the exploratory flight-vs-hardware comparison",
        "T11, data/photoresp, data/osdr"),
    "T10_gcms_power.tsv": ("scripts/17_gcms_power.py",
        "Per-compound GC-MS replicate SD (Col-0 controls, MetaboLights) and plants per arm for 80% power, "
        "unadjusted and Bonferroni over the panel",
        "data/reference/metabolights_group_sds.tsv"),
    # Revision 1 (post hoc, prereg/revision1)
    "T14_reference_strata.tsv": ("scripts/23_reference_strata.py",
        "Per compound: each source's effect, n, support and qualification, and the consensus with all three "
        "sources, darkness only (Zhu, Caldana) and dark+heat (Wu); model verdicts in both", "T01, T02, T03"),
    "T15_reference_sensitivity.tsv": ("scripts/23_reference_strata.py",
        "Consensus over accession agreement (60-90%) x Caldana window (>=60/180/360 min), all and dark-only: "
        "category counts, key calls, model verdict counts", "T01 inputs"),
    "T16_kegg_mapping_report.tsv": ("scripts/23_reference_strata.py",
        "How measured names reached KEGG IDs: candidates, rule matches, reviewed overrides, mapped fraction by "
        "Wu class", "data/reference/kegg_map.tsv, kegg_overrides.tsv"),
    "T17_paintomics_controls.tsv": ("scripts/28_paintomics_controls.py",
        "PaintOmics control jobs E-K2 against their parents (pathways, rank correlation, hubs, raw hub p), "
        "per-omic p in job A, relevant compounds without a hub result, BH ties", "results/paintomics_raw/"),
    "T18_hub_overlap_null.tsv": ("scripts/28_paintomics_controls.py",
        "Lit (A) vs dark (D) hub overlap against hypergeometric and reaction-count-matched nulls; KEGG "
        "reaction count per compound and hub status per job", "T05 inputs; KEGG link/reaction (cached)"),
    "T19_osd678_interaction.tsv": ("scripts/24_osd678_interaction.py",
        "OSD-678 Col-0 flight x light (edgeR TMM, limma-voom): rotation tests on the interaction and marginal "
        "moderated t for eleven sets", "OSDR GLDS-612 rRNArm RSEM counts (data/reference/revision1_sources.tsv)"),
    "T20_rotation_calibration.tsv": ("scripts/25_geneset_diagnostics.py",
        "Type I error of the rotation and a naive null on random sets and genomic blocks, per contrast and "
        "size", "data/osdr, data/atlas GLDS-522 counts, prereg/gene_order.tsv.gz"),
    "T21_exact_permutation.tsv": ("scripts/25_geneset_diagnostics.py",
        "Pre-registered sets under exact 3 v 3 sample permutation (OSD-38, OSD-678 light/dark), and OSD-38 "
        "re-tested on PyDESeq2 statistics from raw counts", "data/osdr; OSDR GLDS-38 RSEM counts"),
    "T22_din_genes.tsv": ("scripts/25_geneset_diagnostics.py",
        "DIN genes per contrast, the set with each gene left out, and per-sample values in the 3 v 3 "
        "contrasts", "data/osdr, T08"),
    "T23_osd38_proteome_checks.tsv": ("scripts/27_osd38_proteome_checks.py",
        "OSD-38 membrane sets under limma, membrane-soluble partition check, truncation and filter "
        "accounting, OSDR file listing", "T11 inputs, T12, data/reference/osd38_osdr_versions.tsv"),
    "T24_hypoxia_positive_control.tsv": ("scripts/26_hypoxia_positive_control.py",
        "Hypoxia, fermentation and DIN sets in OSD-136 and GSE14502 (positive control), and whether the "
        "BRIC-vs-container result is informative", "OSDR OSD-136 DE table; GEO GSE14502, GPL198; T12"),
    "T25_hypoxia_signature.tsv": ("scripts/26_hypoxia_positive_control.py",
        "Addendum: an empirical GSE14502 hypoxia signature in OSD-136 and the OSD-38 comparator and flight "
        "contrasts; ADH1 and PDC1 gene by gene", "as T24; data/osdr/contrast_OSD-38*.csv.gz"),
    "T26_metabolights_provenance.tsv": ("scripts/17_gcms_power.py",
        "Per MetaboLights study: organism, variant, tissue, group rule, groups and n, platform, percent-of-"
        "maximum check", "MetaboLights ISA files (data/reference/metabolights_sources.tsv)"),
}


def main():
    present = sorted(p.name for p in TABLES.glob("T*.tsv"))
    missing = [p for p in present if p not in DESCRIBED]
    if missing:
        sys.exit(f"undescribed tables: {missing}")
    lines = ["table\tbytes\tscript\tcontents\tinputs"]
    for name in present:
        script, what, inputs = DESCRIBED[name]
        lines.append(f"{name}\t{(TABLES / name).stat().st_size}\t{script}\t{what}\t{inputs}")
    (TABLES / "MANIFEST.tsv").write_text("\n".join(lines) + "\n")
    print(f"  MANIFEST.tsv: {len(present)} tables")


if __name__ == "__main__":
    main()
