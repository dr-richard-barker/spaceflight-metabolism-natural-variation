"""Pre-registration: freeze the gene sets, the predictions and every test setting before any test runs.

Writes prereg/gene_sets.tsv, prereg/gene_order.tsv.gz, prereg/settings.json, prereg/PREDICTIONS.tsv
and prereg/PREDICTIONS.md. The test scripts read these frozen files and refuse to run if
prereg/settings.json is missing, so a setting cannot be changed after a result is seen without the
change showing in git.

Where the predictions come from. The first study's model (Photorespiration_multiomics_microgravity,
fvcb.py with LunarLeaf-CFD) puts a microgravity boundary-layer deficit on net assimilation whenever
the plants photosynthesise: 7.3-8.8 % for a rosette at 1000 umol m-2 s-1 and 0.4-2.4 % at OSD-522's
60 umol m-2 s-1 (its T13). In darkness there is no assimilation to lose, so the model predicts a
carbon-starvation signal in lit flight contrasts and none in dark ones. The light/dark reference
(T02) says what carbon starvation measurably does to leaf metabolites; the gene sets carry that
into gene space.
"""

import gzip
import hashlib
import json
import re

import pandas as pd

import genesets
from paths import CACHE, PHOTORESP, PREREG

LIT = ["OSD-522", "OSD-678-light", "OSD-427"]
DARK = ["OSD-38", "OSD-321", "OSD-678-dark"]
# set -> (primary statistic, prediction in lit contrasts, prediction in dark contrasts, basis)
PREDICT = {
    "din_first_study": ("signed", "up", "none",
                        "model: assimilation deficit in light only; RE-TEST, the first study already "
                        "reported this set in all six contrasts (its T07, T11)"),
    "bcaa_degradation": ("signed", "up", "none",
                         "reference: branched-chain amino acids are the strongest darkness response in "
                         "all three sources; their catabolic route should rise with a carbon deficit"),
    "wu_stress_only": ("unsigned", "enriched", "none",
                       "loci that shape metabolism only under 32 C darkness should be perturbed if lit "
                       "flight is starvation-like; weak prior (a locus gene need not be responsive)"),
    "wu_dark_metabolite": ("unsigned", "enriched", "none",
                           "loci of metabolites with a directional darkness response; weak prior"),
    "wu_control_only": ("unsigned", "none", "none",
                        "comparator: control-condition loci, about which the model makes no claim"),
    "wu_shared": ("unsigned", "none", "none",
                  "comparator: loci found in both conditions, about which the model makes no claim"),
}
SETTINGS = {
    "written": "2026-10-30",
    "contrasts": {"lit": LIT, "dark": DARK},
    "gene_statistic": "GeneLab Wald stat (Stat_ column) for the five vendored contrasts; for OSD-522 "
                      "the PyDESeq2 Wald stat recomputed from raw counts as the first study did "
                      "(6 flight vs 6 ground shoots, genes with total count >= 10)",
    "set_statistics": {
        "signed": "mean gene statistic in the set minus the mean over all measured genes",
        "unsigned": "mean |gene statistic| in the set minus the mean |statistic| over all measured genes",
    },
    "primary_statistic": {k: v[0] for k, v in PREDICT.items()},
    "primary_null": {
        "OSD-522": "sample permutation: all 924 ways to label 6 of the 12 samples as flight; gene "
                   "statistic a Welch t on log2(median-of-ratios normalised count + 1)",
        "other": "genomic rotation: the set's gene positions shifted together by a random offset "
                 "along the genome (KEGG order, prereg/gene_order.tsv.gz), preserving set size and "
                 "clustering; 10000 offsets",
    },
    "rotation_draws": 10000,
    "seed": 20261030,
    "min_total_count_osd522": 10,
    "p_two_sided": "(1 + #{|null| >= |obs|}) / (B + 1) for rotation; #{|null| >= |obs|} / 924 for "
                   "the sample permutation, which includes the observed labelling",
    "p_one_sided": "same with null >= obs, for the unsigned statistic",
    "fdr": 0.05,
    "fdr_family": "Benjamini-Hochberg across the 36 primary tests (6 sets x 6 contrasts)",
    "verdict": "observed = up/down (signed) or enriched (unsigned) when BH q < 0.05, else none; "
               "MATCH when observed equals the prediction",
    "secondary": "rotation null for OSD-522 as well; two-sided Mann-Whitney of set log2FC against "
                 "all other genes (the first study's test); the non-primary statistic. Reported, not "
                 "used for verdicts.",
    "compound_rule": {
        "scored_drivers": ["A", "starvation"],
        "rule": "a T06 prediction whose driver is net assimilation (A) or carbon starvation is "
                "compared with the T02 consensus direction: agrees if the signs match, disagrees if "
                "they differ, untested if T02 has no directional status for that KEGG ID",
        "not_comparable": "Vo- and RuBP-driven predictions: darkness stops oxygenation and "
                          "carboxylation together, so a darkness reference cannot test them",
        "blind": False,
        "note": "T02 was built and inspected on 2026-10-30 before this file was written, so this "
                "comparison is a consistency check, not a blind test",
    },
    "paintomics": {
        "layer": "all T02 compounds; value = direction x 0.13052 (the model's carbon-deficit "
                 "magnitude, T06) for robust and single, 0 for none and conflicting; relevant = "
                 "robust only",
        "label": "ground carbon-starvation reference, not flight measurements",
        "job_A": "OSD-522 gene and protein files byte-identical to job 1 (SHA-256 checked) + layer",
        "job_B": "OSD-678-light genes (GeneLab log2FC; relevant = adjusted p < 0.05) + layer",
        "job_C": "OSD-678-dark genes, same rule + layer",
        "databases": ["KEGG", "MapMan"],
        "ai": "off",
        "comparison": "descriptive: pathways at combined-p BH FDR < 0.05 gained or lost, rank "
                      "shifts, hubs at FDR < 0.05",
        "mechanistic_check": "hub raw p-values for compounds relevant in both job 1 and job A are "
                             "expected to be identical, because hub significance depends on the "
                             "gene layer and the KEGG network, not on the metabolite values; FDR "
                             "may differ because more compounds enter the family",
    },
    "power": {
        "variance": "Col-0 wild-type control groups of MetaboLights GC-MS leaf or rosette studies; "
                    "within-group SD of log2 intensity per metabolite, median across groups",
        "test": "two-sided two-sample t, equal n per arm, noncentral t",
        "alpha": 0.05,
        "power": 0.8,
        "effects_log2": {"model_Q60_rosette_max": "log2(1 - 0.0237)",
                         "model_Q1000_T06": 0.13052, "half_fold": 0.585, "two_fold": 1.0},
        "multiplicity": "unadjusted and Bonferroni across the panel; BH lies between them",
    },
}


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def gene_order():
    rows = []
    for ln in (CACHE / "kegg" / "list_ath.tsv").read_text().splitlines():
        p = ln.split("\t")
        m = re.fullmatch(r"ath:(AT[1-5]G\d{5})", p[0])
        pos = re.match(r"([1-5]):(?:complement\()?(?:join\()?(\d+)", p[2]) if len(p) > 2 else None
        if m and pos:
            rows.append((m.group(1), int(pos.group(1)), int(pos.group(2))))
    df = pd.DataFrame(rows, columns=["gene", "chrom", "start"]).sort_values(["chrom", "start", "gene"])
    out = PREREG / "gene_order.tsv.gz"
    with gzip.GzipFile(out, "wb", mtime=0) as gz:
        gz.write(df.to_csv(sep="\t", index=False).encode())
    return len(df)


def main():
    sets = genesets.freeze()
    n_order = gene_order()
    q60 = pd.read_csv(PHOTORESP / "T13_light_sensitivity.tsv", sep="\t", comment="#")
    d = q60[(q60.quantity == "deficit") & (q60.Q_umol_m2_s == 60) & (q60.scale == "rosette")]
    assert float(d.deficit_max_percent.iloc[0]) == 2.37, "T13 Q=60 rosette deficit changed"
    (PREREG / "settings.json").write_text(json.dumps(SETTINGS, indent=2) + "\n")

    rows = []
    for name, (stat, lit, dark, basis) in PREDICT.items():
        for c in LIT + DARK:
            rows.append({"set": name, "contrast": c, "light": "lit" if c in LIT else "dark",
                         "statistic": stat, "predicted": lit if c in LIT else dark,
                         "n_genes": len(sets[name]), "basis": basis})
    pred = pd.DataFrame(rows)
    pred.to_csv(PREREG / "PREDICTIONS.tsv", sep="\t", index=False)

    lines = [
        "# Pre-registered predictions",
        "",
        f"Written {SETTINGS['written']} by `scripts/11_prereg.py`, after the light/dark reference "
        "(T01, T02) and the gene sets were built and **before any gene-set statistic was computed on "
        "any flight contrast in this repository**. The commit that adds this file precedes the "
        "commit that adds the tests; `git log` shows the order.",
        "",
        "## What was already known when this was written",
        "",
        "- The first study's results: its T07 (OSD-522 transcriptome and proteome against KEGG and "
        "marker sets) and T11 (the same sets in all six contrasts), including the DIN set. The "
        "`din_first_study` rows below are therefore **re-tests with stronger nulls, not blind "
        "predictions**.",
        "- T02, the light/dark metabolite reference. The compound comparison (T03) is a consistency "
        "check, not a blind test.",
        "- Nothing about `bcaa_degradation` or the four Wu sets in any flight contrast.",
        "",
        "## The model's light/dark prediction",
        "",
        "The first study's model puts a microgravity boundary-layer deficit on net assimilation "
        "whenever the plants photosynthesise (7.3-8.8 % for a rosette at 1000 umol m-2 s-1; "
        "0.4-2.4 % at OSD-522's 60 umol m-2 s-1, its T13). In darkness there is no assimilation to "
        "lose. So: a carbon-starvation signal in **lit** flight contrasts, **none** in dark ones.",
        "",
        f"- Lit: {', '.join(LIT)}",
        f"- Dark: {', '.join(DARK)}",
        "",
        "## Gene sets (frozen in `gene_sets.tsv`)",
        "",
        "| Set | Genes | Definition | Statistic | Lit | Dark |",
        "|---|---|---|---|---|---|",
    ]
    for name, (stat, lit, dark, basis) in PREDICT.items():
        lines.append(f"| `{name}` | {len(sets[name])} | {genesets.DESCRIPTION[name]} | {stat} | "
                     f"{lit} | {dark} |")
    lines += ["", "Basis for each prediction:", ""]
    lines += [f"- `{k}`: {v[3]}." for k, v in PREDICT.items()]
    s = SETTINGS
    lines += [
        "",
        "## Tests",
        "",
        f"- Gene statistic: {s['gene_statistic']}.",
        f"- Signed set statistic: {s['set_statistics']['signed']}. Unsigned: "
        f"{s['set_statistics']['unsigned']}.",
        f"- Primary null, OSD-522: {s['primary_null']['OSD-522']}.",
        f"- Primary null, the other five: {s['primary_null']['other']}. Seed {s['seed']}.",
        f"- p-values: two-sided {s['p_two_sided']}; one-sided {s['p_one_sided']}.",
        f"- Multiplicity: {s['fdr_family']}, q < {s['fdr']}.",
        f"- Verdict: {s['verdict']}.",
        f"- Secondary, reported only: {s['secondary']}",
        "",
        "## Compound predictions against the reference (T03)",
        "",
        f"- Scored: predictions whose driver is {' or '.join(s['compound_rule']['scored_drivers'])}. "
        f"{s['compound_rule']['rule']}.",
        f"- Not comparable: {s['compound_rule']['not_comparable']}.",
        f"- {s['compound_rule']['note']}.",
        "",
        "## PaintOmics with the data-derived layer",
        "",
        f"- Layer: {s['paintomics']['layer']}. Labelled everywhere: *{s['paintomics']['label']}*.",
        f"- Job A: {s['paintomics']['job_A']}. Job B: {s['paintomics']['job_B']}. Job C: "
        f"{s['paintomics']['job_C']}. Databases {' + '.join(s['paintomics']['databases'])}, AI "
        f"{s['paintomics']['ai']}.",
        f"- Comparison: {s['paintomics']['comparison']}.",
        f"- Mechanistic check: {s['paintomics']['mechanistic_check']}.",
        "",
        "## GC-MS panel power",
        "",
        f"- Variance: {s['power']['variance']}.",
        f"- Test: {s['power']['test']}; alpha {s['power']['alpha']}, power {s['power']['power']}.",
        "- Effects (log2): " + "; ".join(f"{k} = {v}" for k, v in s['power']['effects_log2'].items())
        + ".",
        f"- Multiplicity: {s['power']['multiplicity']}.",
        "",
        "## Frozen files",
        "",
        "| File | SHA-256 |",
        "|---|---|",
    ]
    for f in ("gene_sets.tsv", "gene_order.tsv.gz", "settings.json", "PREDICTIONS.tsv"):
        lines.append(f"| `{f}` | `{sha(PREREG / f)}` |")
    (PREREG / "PREDICTIONS.md").write_text("\n".join(lines) + "\n")
    print(f"  prereg: {len(sets)} sets, {len(pred)} predictions, {n_order} ordered genes")


if __name__ == "__main__":
    main()
