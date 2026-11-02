"""Build the PaintOmics upload bundles for jobs A, B and C (data/paintomics_upload/).

The metabolite layer is the pre-registered one (prereg/settings.json, paintomics): every T02 compound,
valued direction x 0.13052 (the first study's carbon-deficit magnitude, T06) when its light/dark
response is robust or single and 0 otherwise, with the robust compounds marked relevant. It is a
GROUND CARBON-STARVATION REFERENCE, NOT FLIGHT MEASUREMENTS: under the first study's hypothesis that
lit flight is a mild carbon deficit, it is what that deficit would do to these compounds' direction.

Compounds are submitted under their primary KEGG name (the first name of the KEGG COMPOUND entry,
from rest.kegg.jp/list/compound), so PaintOmics' KEGG-derived name matcher sees the name KEGG itself
uses; scripts/15_paintomics_tables.py checks every match against the intended ID.

  A  OSD-522 genes and proteins, byte-identical to the first study's job 1 files, + layer
  B  OSD-678-light genes (GeneLab log2FC; relevant = adjusted p < 0.05) + layer
  C  OSD-678-dark genes, same rule, + layer
  D  OSD-38 (BRIC-20, dark sealed BRIC-PDFU) genes, same rule, + membrane proteins from T11 (relevant at the
     paper's rule, as pre-registered in prereg/osd38) + layer

Revision 1 controls (post hoc, prereg/revision1 R6) reuse these files, plus two more inputs:
  F  metabolomics_scaled_values.tab: direction x |reference_z|, rescaled so the median over directional
     compounds is 0.13052; same names and relevant list as the layer
  K  osd678_light_genes_relevant_k1.tab, _k2.tab: job B's relevant genes randomly downsampled to job
     C's relevant count (seeds 1 and 2)
"""

import hashlib
import json
import shutil

import numpy as np
import pandas as pd

from genesets import kegg
from paths import DATA, OSDR, PHOTORESP, PREREG, TABLES

UPLOAD = DATA / "paintomics_upload"
GAIN = 0.13052


def primary_names():
    out = {}
    for ln in kegg("list/compound", "list_compound.tsv").splitlines():
        cid, names = ln.split("\t")
        out[cid.replace("cpd:", "")] = names.split(";")[0].strip()
    return out


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def layer():
    rule = json.loads((PREREG / "settings.json").read_text())["paintomics"]
    assert "0.13052" in rule["layer"] and "relevant = robust only" in rule["layer"]
    ref = pd.read_csv(TABLES / "T02_light_dark_reference.tsv", sep="\t")
    names = primary_names()
    ref["kegg_name"] = ref.kegg_id.map(names)
    assert ref.kegg_name.notna().all(), ref[ref.kegg_name.isna()]
    assert ref.kegg_name.is_unique, ref[ref.kegg_name.duplicated(keep=False)]
    directional = ref.status.isin(["robust", "single"])
    ref["value"] = (ref.direction * GAIN).where(directional, 0.0)
    ref["relevant"] = ref.status == "robust"
    return ref


def write_layer(ref):
    vals = UPLOAD / "metabolomics_values.tab"
    vals.write_text("#compound\tFLT_vs_GC\n" + "".join(
        f"{r.kegg_name}\t{r.value:.6f}\n" for r in ref.itertuples()))
    (UPLOAD / "metabolomics_relevant.tab").write_text(
        "".join(f"{n}\n" for n in ref.loc[ref.relevant, "kegg_name"]))
    ref[["kegg_id", "name", "kegg_name", "status", "direction", "value", "relevant"]].to_csv(
        UPLOAD / "metabolomics_key.tsv", sep="\t", index=False)


def gene_layer(contrast, tag):
    d = pd.read_csv(OSDR / f"contrast_{contrast}.csv.gz", usecols=["TAIR", "log2fc", "padj"])
    (UPLOAD / f"{tag}_values.tab").write_text("#geneID\tFLT_vs_GC\n" + "".join(
        f"{r.TAIR}\t{r.log2fc:.6f}\n" for r in d.itertuples()))
    rel = d.loc[d.padj < 0.05, "TAIR"]
    (UPLOAD / f"{tag}_relevant.tab").write_text("".join(f"{g}\n" for g in rel))
    return len(d), len(rel)


def main():
    UPLOAD.mkdir(parents=True, exist_ok=True)
    ref = layer()
    write_layer(ref)
    for f in ("gene_expression_values.tab", "gene_expression_relevant.tab",
              "proteomics_values.tab", "proteomics_relevant.tab"):
        shutil.copyfile(PHOTORESP / f, UPLOAD / f"osd522_{f}")
        assert sha(PHOTORESP / f) == sha(UPLOAD / f"osd522_{f}")
    b = gene_layer("OSD-678-light", "osd678_light_genes")
    c = gene_layer("OSD-678-dark", "osd678_dark_genes")
    d = gene_layer("OSD-38", "osd38_genes")
    mem = pd.read_csv(TABLES / "T11_osd38_proteome.tsv", sep="\t").query("fraction == 'membrane'")
    (UPLOAD / "osd38_membrane_proteins_values.tab").write_text("#proteinID\tFLT_vs_GC\n" + "".join(
        f"{r.locus}\t{r.log2fc:.6f}\n" for r in mem.itertuples()))
    (UPLOAD / "osd38_membrane_proteins_relevant.tab").write_text(
        "".join(f"{r.locus}\n" for r in mem[mem.relevant].itertuples()))
    r6 = json.loads((PREREG / "revision1" / "settings.json").read_text())["R6_paintomics"]["jobs"]
    assert "0.13052" in r6["F"] and "seed 1" in r6["K1"]
    d_ = ref.status.isin(["robust", "single"])
    scale = GAIN / ref.loc[d_, "reference_z"].abs().median()
    scaled = (ref.direction * ref.reference_z.abs() * scale).where(d_, 0.0)
    (UPLOAD / "metabolomics_scaled_values.tab").write_text("#compound\tFLT_vs_GC\n" + "".join(
        f"{n}\t{v:.6f}\n" for n, v in zip(ref.kegg_name, scaled)))
    rel_b = sorted((UPLOAD / "osd678_light_genes_relevant.tab").read_text().split())
    n_c = len((UPLOAD / "osd678_dark_genes_relevant.tab").read_text().split())
    for seed in (1, 2):
        pick = sorted(np.random.default_rng(seed).choice(rel_b, n_c, replace=False))
        (UPLOAD / f"osd678_light_genes_relevant_k{seed}.tab").write_text("".join(f"{g}\n" for g in pick))
    lines = [
        "PaintOmics upload bundles, built by scripts/13_build_layer.py",
        "",
        "Metabolite layer (all jobs): GROUND CARBON-STARVATION REFERENCE, NOT FLIGHT MEASUREMENTS.",
        f"  source:   results/tables/T02_light_dark_reference.tsv ({len(ref)} KEGG compounds)",
        f"  value:    direction x {GAIN} for robust and single responses, 0 otherwise",
        f"  relevant: robust ({int(ref.relevant.sum())} compounds)",
        "  names:    primary KEGG COMPOUND name (metabolomics_key.tsv maps name -> KEGG ID)",
        "",
        "Job A: osd522_* gene and protein files, byte-identical to the first study's job 1 upload",
        "       (data/photoresp/, SHA-256 checked) + the metabolite layer",
        f"Job B: osd678_light_genes_* ({b[0]} genes, {b[1]} at GeneLab adjusted p < 0.05) + layer",
        f"Job C: osd678_dark_genes_* ({c[0]} genes, {c[1]} at GeneLab adjusted p < 0.05) + layer",
        f"Job D: osd38_genes_* ({d[0]} genes, {d[1]} at GeneLab adjusted p < 0.05) + "
        f"osd38_membrane_proteins_* ({len(mem)} proteins, AGI loci; {int(mem.relevant.sum())} at p <= 0.05 and "
        "|log2FC| > 0.2, from T11) + layer",
        "",
        "Revision 1 controls (post hoc, prereg/revision1 R6):",
        "Job E: job A without the metabolite layer",
        f"Job F: job A with metabolomics_scaled_values.tab (direction x |reference_z| x {scale:.5f}; median "
        f"|value| over directional compounds {GAIN})",
        "Job G: job A genes + layer;  Job H: job A proteins + layer",
        "Job I: job D genes + layer;  Job J: job D membrane proteins + layer",
        f"Job K1, K2: job B with osd678_light_genes_relevant_k1/k2.tab ({n_c} of {len(rel_b)} relevant genes, "
        "random, seeds 1 and 2) + layer",
    ]
    (UPLOAD / "METHOD.txt").write_text("\n".join(lines) + "\n")
    print("\n".join("  " + ln for ln in lines))


if __name__ == "__main__":
    main()
