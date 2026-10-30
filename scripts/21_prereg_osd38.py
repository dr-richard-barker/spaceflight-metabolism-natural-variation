"""Pre-registration addendum for OSD-38 (BRIC-20): frozen before any OSD-38 protein or hardware statistic.

Writes prereg/osd38/ (gene_sets.tsv, settings.json, PREDICTIONS.tsv, PREDICTIONS.md). The tests in
scripts/22_osd38_tests.py read these files and refuse to run without them. The genome order is the
one already frozen in prereg/gene_order.tsv.gz.

Layers (all dark-grown, 3-day Col-0 seedlings):
  rna_flight        OSD-38 RNA-seq, flight vs ground, both in sealed BRIC-PDFU (GeneLab Wald stat)
  membrane_flight   OSD-38 membrane iTRAQ proteome, flight vs ground (block 1 of the truncated table)
  soluble_flight    OSD-38 soluble iTRAQ proteome, flight vs ground; SECONDARY, because the table is
                    truncated to its 437 top-ranked protein groups
  rna_hardware      OSD-38 RNA-seq on the ground, BRIC-PDFU vs plates in closed containers
"""

import hashlib
import json
import re

import pandas as pd

import genesets
from paths import CACHE, PREREG, REF

OUT = PREREG / "osd38"
MARKER = {  # the first study's own rules (Photorespiration_multiomics_microgravity scripts/genesets.py @515226d)
    "photosynthesis_apparatus": r"photosystem I{1,2} subunit|light harvesting complex|chlorophyll A/B binding",
    "fermentation": r"pyruvate decarboxylase|alcohol dehydrogenase 1;|lactate dehydrogenase",
    "hypoxia_responsive": r"[Hh]ypoxia-responsive",
}
PATHWAYS = {"ath00630": "Glyoxylate & dicarboxylate metabolism", "ath00710": "Carbon fixation",
            "ath00500": "Starch & sucrose metabolism", "ath00010": "Glycolysis / gluconeogenesis"}
LAYERS = {
    "rna_flight": "OSD-38 RNA-seq, flight vs ground (both sealed BRIC-PDFU, dark)",
    "membrane_flight": "OSD-38 membrane iTRAQ proteome, flight vs ground",
    "soluble_flight": "OSD-38 soluble iTRAQ proteome, flight vs ground (truncated: top 437 groups)",
    "rna_hardware": "OSD-38 RNA-seq, ground: BRIC-PDFU vs closed container (both dark)",
}
PRIMARY = ["rna_flight", "membrane_flight", "rna_hardware"]
# Hardware: the first study's enclosure mass balance (its vendored LunarLeaf-CFD T8: a dark BRIC goes from
# 21 % to 5 % O2 in 6.5 days and to 1 % CO2 in 9.4 h, for ITS assumed canister and rosette) is the only
# basis for a prediction, and only conditionally: IF the BRIC-PDFU holds less gas per seedling than the
# Ohio container (not reported), O2 falls further in the BRIC, so hypoxia and fermentation sets rise.
HARDWARE = {"hypoxia_responsive": "up", "fermentation": "up"}
KNOWN_RNA_FLIGHT = ("the first study reported every one of these sets for the OSD-38 transcriptome by "
                    "Mann-Whitney (its T11), and this study's T08 tested DIN and BCAA degradation there; "
                    "rna_flight rows are re-tests, not blind")


def sets():
    fro = pd.read_csv(PREREG / "gene_sets.tsv", sep="\t")
    out = {k: set(fro[fro.set == k].gene) for k in ("din_first_study", "bcaa_degradation")}
    ann = []
    for ln in genesets.kegg("list/ath", "list_ath.tsv").splitlines():
        p = ln.split("\t")
        if len(p) >= 4:
            ann.append((p[0].replace("ath:", ""), p[3]))
    ann = pd.DataFrame(ann, columns=["agi", "desc"])
    for k, rx in MARKER.items():
        out[k] = set(ann[ann.desc.str.contains(rx, regex=True)].agi)
    for pid in PATHWAYS:
        out[pid] = genesets.kegg_pathway(pid)
    return {k: {g for g in v if re.fullmatch(r"AT[1-5CM]G\d{5}", g)} for k, v in out.items()}


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    s = sets()
    pd.DataFrame([(k, g) for k, v in s.items() for g in sorted(v)], columns=["set", "gene"]).to_csv(
        OUT / "gene_sets.tsv", sep="\t", index=False)
    cov = pd.read_csv(REF / "osd38_protein_coverage.tsv", sep="\t")
    measured = {"membrane_flight": set(cov[cov.fraction == "membrane"].locus),
                "soluble_flight": set(cov[cov.fraction == "soluble"].locus)}
    settings = {
        "written": "2026-10-30",
        "layers": LAYERS, "primary_layers": PRIMARY,
        "protein_quantification": {
            "psm_filter": "Score > 20 and all six iTRAQ reporter intensities > 0",
            "loci": "AGI locus from Sequence Id; loci with at least 2 distinct peptides (the paper's rule)",
            "soluble_rows": "protein groups 1-437 (the table is cut inside group 438)",
            "membrane_rows": "block 1, rows 1-48,799 (groups 1-2,948); block 2 excluded",
            "normalisation": "log2; each channel centred on its median; each PSM centred on its own mean",
            "summary": "protein value per channel = median over its PSMs",
            "contrast": "log2FC = mean(flight S1-S3) - mean(ground GC1-GC3)",
            "test": "two-sided Student's t, 3 vs 3 (the paper's test)",
            "relevant_for_paintomics": "p <= 0.05 and |log2FC| > 0.2 (the paper's rule)",
        },
        "gene_statistic": {"rna_flight": "GeneLab Wald stat (data/osdr/contrast_OSD-38.csv.gz)",
                           "rna_hardware": "GeneLab Wald stat, (GC)v(RNAlat) (data/osdr/contrast_OSD-38-hardware.csv.gz)",
                           "membrane_flight": "protein Student's t", "soluble_flight": "protein Student's t"},
        "set_statistic": "signed: mean gene/protein statistic in the set minus the mean over all measured features",
        "null": "genomic rotation over the measured features in prereg/gene_order.tsv.gz order, 10000 offsets, seed 20261030",
        "min_measured": 3,
        "fdr": "Benjamini-Hochberg over every testable primary row (rna_flight, membrane_flight, rna_hardware); "
               "soluble_flight rows are reported with their own BH and are not used for verdicts",
        "verdict": "observed = up/down when q < 0.05, else none; MATCH when observed equals the prediction",
        "concordance": "descriptive: Spearman correlation of protein and transcript log2FC over shared loci, "
                       "OSD-38 (membrane, soluble) and OSD-522 (the first study's proteome, UniProt mapped to "
                       "AGI by KEGG conv/uniprot/ath)",
        "paintomics": "job D: OSD-38 genes (log2FC; relevant = adjusted p < 0.05) + membrane proteins "
                      "(relevant = p <= 0.05 and |log2FC| > 0.2) + the reference metabolite layer, same "
                      "settings as job A; compared descriptively with job A (OSD-522, lit BRIC-LED)",
    }
    (OUT / "settings.json").write_text(json.dumps(settings, indent=2) + "\n")

    rows = []
    for layer in LAYERS:
        for name, genes in s.items():
            n = len(genes & measured[layer]) if layer in measured else None
            if layer == "rna_hardware":
                pred, basis = HARDWARE.get(name, "none"), (
                    "conditional on the BRIC-PDFU holding less gas per seedling than the container "
                    "(first study's enclosure mass balance)" if name in HARDWARE else
                    "no mechanism in the first study's model for this set in a dark enclosure")
            else:
                pred, basis = "none", "first study's model: no assimilation to lose in the dark"
                if layer == "rna_flight":
                    basis += "; RE-TEST: " + KNOWN_RNA_FLIGHT
            rows.append({"layer": layer, "set": name, "n_set": len(genes), "n_measured_protein": n,
                         "predicted": pred, "primary": layer in PRIMARY, "basis": basis})
    pred = pd.DataFrame(rows)
    pred.to_csv(OUT / "PREDICTIONS.tsv", sep="\t", index=False)

    md = ["# Pre-registered predictions: OSD-38 (BRIC-20) addendum", "",
          f"Written {settings['written']} by `scripts/21_prereg_osd38.py`, after the OSD-38 inputs were "
          "fetched and protein coverage counted (`data/reference/osd38_protein_coverage.tsv`, which "
          "lists detected loci only), and **before any OSD-38 protein fold change or hardware-contrast "
          "statistic was computed**. The commit adding this folder precedes the commit adding the tests.",
          "", "## Already known when this was written", "",
          f"- OSD-38 transcriptome, flight vs ground: {KNOWN_RNA_FLIGHT}.",
          "- The OSDR processed proteome tables are truncated at 65,535 rows (see `scripts/20_osd38_inputs.py`).",
          "- OSD-38 seedlings grew on 0.5x MS with 1 % sucrose (OSDR protocol), so carbon starvation in the "
          "classical sense is not expected in either arm.",
          "", "## Layers", ""]
    md += [f"- `{k}`: {v}{' (primary)' if k in PRIMARY else ' (secondary)'}" for k, v in LAYERS.items()]
    md += ["", "## Predictions", "",
           "Flight layers: the first study's model predicts **no signal for any set** in dark flight. "
           "Hardware layer: hypoxia-responsive and fermentation sets **up** in BRIC-PDFU vs container, "
           "conditional on the BRIC holding less gas per seedling (the container volume is not reported); "
           "no prediction for the other sets (predicted none).", "",
           "| Layer | Set | Set size | Measured proteins | Predicted |", "|---|---|---|---|---|"]
    for r in pred.itertuples():
        n = "" if pd.isna(r.n_measured_protein) else int(r.n_measured_protein)
        md.append(f"| {r.layer} | `{r.set}` | {r.n_set} | {n} | {r.predicted} |")
    md += ["", "## Settings", ""]
    md += [f"- {k}: {v}" for k, v in settings["protein_quantification"].items()]
    md += [f"- Set statistic: {settings['set_statistic']}.", f"- Null: {settings['null']}.",
           f"- Minimum measured members: {settings['min_measured']}.", f"- Multiplicity: {settings['fdr']}.",
           f"- Verdict: {settings['verdict']}.", f"- Concordance: {settings['concordance']}.",
           f"- PaintOmics: {settings['paintomics']}.", "", "## Frozen files", "", "| File | SHA-256 |", "|---|---|"]
    md += [f"| `{f}` | `{sha(OUT / f)}` |" for f in ("gene_sets.tsv", "settings.json", "PREDICTIONS.tsv")]
    (OUT / "PREDICTIONS.md").write_text("\n".join(md) + "\n")
    print(f"  prereg/osd38: {len(s)} sets x {len(LAYERS)} layers; protein rows with >= 3 measured: "
          f"{int((pred.n_measured_protein >= 3).sum())}")


if __name__ == "__main__":
    main()
