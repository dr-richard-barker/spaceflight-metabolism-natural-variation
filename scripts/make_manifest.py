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
        "Per-compound GC-MS replicate SD (Col-0 controls, MetaboLights) and plants per arm for 80% power",
        "data/reference/metabolights_group_sds.tsv"),
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
