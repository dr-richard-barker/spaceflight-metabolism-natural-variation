"""Gene sets for the light/dark tests, built from KEGG and the Wu 2018 GWAS loci, never from memory.

  din_first_study     the first study's carbon-starvation set, rebuilt with its own rule (KEGG gene
                      descriptions matching r"\\bDIN\\d+;"), so its result can be re-tested here
  bcaa_degradation    KEGG ath00280 (valine, leucine and isoleucine degradation): the catabolic
                      route of the branched-chain amino acids, the strongest darkness response in
                      every reference source (T02)
  wu_control_only     genes in Wu 2018 GWAS regions found only under the control condition (light)
  wu_stress_only      genes in regions found only under the stress condition (32 C darkness)
  wu_shared           genes in regions found under both
  wu_dark_metabolite  genes in Wu loci (either condition) whose traits include a peak of a
                      metabolite with a directional light/dark reference (T02 robust or single)

Wu regions and loci come from the atlas (data/atlas/wu_condition_regions.tsv, wu_gwas_loci.csv); a
region's genes are the union of its loci's genes. KEGG is read through the REST API and cached in
data/cache/kegg. `freeze()` writes the membership to prereg/gene_sets.tsv; the tests read that
frozen file, so a later KEGG release cannot change a pre-registered set.
"""

import re
import urllib.request

import pandas as pd

from paths import ATLAS, CACHE, PREREG, REF, TABLES

KEGG = "https://rest.kegg.jp/{}"
UA = {"User-Agent": "spaceflight-metabolism-natural-variation/1.0"}
DIN_RULE = r"\bDIN\d+;"
DIRECTIONAL = ("robust", "single")
FROZEN = PREREG / "gene_sets.tsv"
DESCRIPTION = {
    "din_first_study": "First study's DIN set (KEGG description regex \\bDIN\\d+;)",
    "bcaa_degradation": "KEGG ath00280 valine, leucine and isoleucine degradation",
    "wu_control_only": "Wu 2018 GWAS regions found only under control (light)",
    "wu_stress_only": "Wu 2018 GWAS regions found only under stress (32 C darkness)",
    "wu_shared": "Wu 2018 GWAS regions found under both conditions",
    "wu_dark_metabolite": "Wu 2018 loci of metabolites with a directional light/dark reference",
}


def kegg(endpoint, name):
    path = CACHE / "kegg" / name
    if not path.exists():
        path.parent.mkdir(parents=True, exist_ok=True)
        with urllib.request.urlopen(urllib.request.Request(KEGG.format(endpoint), headers=UA),
                                    timeout=120) as r:
            path.write_text(r.read().decode())
    return path.read_text()


def din_first_study():
    out = set()
    for line in kegg("list/ath", "list_ath.tsv").splitlines():
        parts = line.split("\t")
        if len(parts) >= 4 and re.search(DIN_RULE, parts[3]):
            out.add(parts[0].replace("ath:", ""))
    return out


def kegg_pathway(pid):
    text = kegg(f"link/ath/path:{pid}", f"link_{pid}.tsv")
    return {ln.split("\t")[1].replace("ath:", "") for ln in text.splitlines() if "\t" in ln}


def wu_region_sets():
    loci = pd.read_csv(ATLAS / "wu_gwas_loci.csv")
    genes = {r.locus_id: set(str(r.genes).split(";")) for r in loci.itertuples()}
    reg = pd.read_csv(ATLAS / "wu_condition_regions.tsv", sep="\t")
    label = {"control": "wu_control_only", "stress": "wu_stress_only", "control+stress": "wu_shared"}
    sets = {v: set() for v in label.values()}
    for r in reg.itertuples():
        for locus in r.loci.split(";"):
            sets[label[r.conditions]] |= genes[locus]
    return sets


def wu_dark_metabolite_loci():
    """Loci whose traits include a peak of a Wu metabolite with a directional reference."""
    ref = pd.read_csv(TABLES / "T02_light_dark_reference.tsv", sep="\t")
    src = pd.read_csv(TABLES / "T01_light_dark_sources.tsv", sep="\t")
    keep = set(ref.loc[ref.status.isin(DIRECTIONAL), "kegg_id"])
    names = set(src.loc[(src.source == "wu_stress") & src.kegg_id.isin(keep), "name"])
    ident = pd.read_csv(ATLAS / "wu_identified_metabolites.csv")
    ident = ident[ident.metabolite_name.isin(names)]
    peaks = set(ident[["peak_control_pos", "peak_control_neg", "peak_stress_pos",
                       "peak_stress_neg"]].stack().dropna())
    loci = pd.read_csv(ATLAS / "wu_gwas_loci.csv")
    hit = loci[[bool(set(t.split(";")) & peaks) for t in loci.traits]]
    return set(";".join(hit.genes).split(";")), hit, sorted(names)


def build():
    sets = {"din_first_study": din_first_study(), "bcaa_degradation": kegg_pathway("ath00280")}
    sets.update(wu_region_sets())
    sets["wu_dark_metabolite"] = wu_dark_metabolite_loci()[0]
    return {k: {g for g in v if re.fullmatch(r"AT[1-5CM]G\d{5}", g)} for k, v in sets.items()}


def freeze():
    sets = build()
    rows = [(k, g) for k, v in sets.items() for g in sorted(v)]
    PREREG.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows, columns=["set", "gene"]).to_csv(FROZEN, sep="\t", index=False)
    return sets


def load():
    df = pd.read_csv(FROZEN, sep="\t")
    return {k: set(g.gene) for k, g in df.groupby("set", sort=False)}


if __name__ == "__main__":
    s = build()
    for k, v in s.items():
        print(f"  {k}: {len(v)} genes — {DESCRIPTION[k]}")
    _, hit, names = wu_dark_metabolite_loci()
    print(f"  wu_dark_metabolite: {len(hit)} loci ({hit.condition.value_counts().to_dict()}) "
          f"from {len(names)} metabolites: {names}")
    din, bcaa = s["din_first_study"], s["bcaa_degradation"]
    print(f"  overlap DIN x ath00280: {sorted(din & bcaa)}")
    for a in ("wu_control_only", "wu_stress_only", "wu_shared", "wu_dark_metabolite"):
        print(f"  overlap {a} x DIN {len(s[a] & din)}, x ath00280 {len(s[a] & bcaa)}")
