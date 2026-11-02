"""Checks on the OSD-38 proteome result (T23).

Revision 1, post hoc: the settings are in prereg/revision1/settings.json (R5), committed before this
script ran. The pre-registered quantification (scripts/20_osd38_inputs.py) and tests (T12) are
unchanged.

  versions   query the OSDR file API for the two processed proteome tables and record what it lists
             (data/reference/osd38_osdr_versions.tsv); network, run once
  build      T23:
             - limma moderated t (lmFit + eBayes, scripts/r/limma_two_group.R) on the membrane
               protein-by-channel matrix, and the five membrane sets re-tested on it with the rotation
               null, next to their T12 results
             - proteins quantified in both fractions: Spearman of membrane against soluble log2FC, and
               the carbon-pathway members among them, soluble fraction
             - accounting: protein groups and loci kept at each step of the truncation rule and filters
             - the OSDR file listing from `versions`
"""

import argparse
import json
import subprocess
import urllib.request

import numpy as np
import pandas as pd
from scipy import stats

from paths import CACHE, PREREG, REF, ROOT, TABLES

R5 = json.loads((PREREG / "revision1" / "settings.json").read_text())["R5_osd38_proteome"]
VERSIONS = REF / "osd38_osdr_versions.tsv"
API = "https://osdr.nasa.gov/osdr/data/osd/files/38{}"
SETS = ["bcaa_degradation", "ath00630", "ath00710", "ath00500", "ath00010"]
CARBON = SETS[1:]


def load(name):
    import importlib.util
    spec = importlib.util.spec_from_file_location(name.replace(".py", ""), ROOT / "scripts" / name)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


T16 = load("16_gwas_light_dark_tests.py")
INPUTS = load("20_osd38_inputs.py")
ORDER = pd.read_csv(PREREG / "gene_order.tsv.gz", sep="\t").gene.tolist()


def versions():
    rows = []
    for q in ("", "?all_files=true"):
        req = urllib.request.Request(API.format(q), headers=INPUTS.UA)
        with urllib.request.urlopen(req, timeout=120) as r:
            d = json.load(r)
        for s in d["studies"].values():
            for f in s["study_files"]:
                if "_proteome_WyattLab_" in f["file_name"] and f["file_name"].endswith((".csv", ".xlsx")):
                    rows.append({"query": q or "(default)", "file_name": f["file_name"], "file_size": f["file_size"],
                                 "date_created_utc": pd.to_datetime(f["date_created"], unit="s").strftime("%Y-%m-%d"),
                                 "date_updated": f.get("date_updated") or "", "files_listed": len(s["study_files"])})
    t = pd.DataFrame(rows)
    t.to_csv(VERSIONS, sep="\t", index=False)
    print(t.to_string(index=False))


def limma_membrane():
    prot, n = INPUTS.matrix("membrane")
    work = CACHE / "osd38" / "limma"
    work.mkdir(parents=True, exist_ok=True)
    prot.to_csv(work / "membrane_expr.tsv", sep="\t")
    grp = ["treated" if c in INPUTS.FLIGHT else "control" for c in prot.columns]
    pd.DataFrame({"sample": prot.columns, "group": grp}).to_csv(work / "membrane_design.tsv", sep="\t", index=False)
    subprocess.run(["Rscript", str(ROOT / "scripts" / "r" / "limma_two_group.R"), str(work / "membrane_expr.tsv"),
                    str(work / "membrane_design.tsv"), str(work / "membrane_limma.tsv")], check=True)
    return pd.read_csv(work / "membrane_limma.tsv", sep="\t").set_index("feature")


def accounting():
    rows = []
    raw = pd.read_csv(INPUTS.OSD38 / INPUTS.FILES["soluble"])
    raw["locus"] = raw["Sequence Id"].str.extract(r"^(AT[1-5CM]G\d{5})", expand=False)
    kept = raw[raw["Protein Group"] <= INPUTS.SOLUBLE_LAST_COMPLETE]
    rows += [("soluble", "data rows in the OSDR table", len(raw)),
             ("soluble", "last protein group present (cut inside it)", int(raw["Protein Group"].max())),
             ("soluble", f"protein groups kept (1-{INPUTS.SOLUBLE_LAST_COMPLETE}, complete)", int(kept["Protein Group"].nunique())),
             ("soluble", "AGI loci in the kept groups", int(kept.locus.nunique())),
             ("soluble", "loci that appear in more than one kept group",
              int((kept.groupby("locus")["Protein Group"].nunique() > 1).sum())),
             ("soluble", "loci after the score, channel and two-peptide filters", int(INPUTS.psms("soluble").locus.nunique()))]
    mem = pd.read_csv(INPUTS.OSD38 / INPUTS.FILES["membrane"])
    b1 = mem.iloc[:INPUTS.MEMBRANE_BLOCK1_ROWS].assign(
        locus=lambda d: d["Sequence Id"].str.extract(r"^(AT[1-5CM]G\d{5})", expand=False))
    rows += [("membrane", "data rows in the OSDR table", len(mem)),
             ("membrane", "rows in block 1 (used)", INPUTS.MEMBRANE_BLOCK1_ROWS),
             ("membrane", "protein groups in block 1", int(b1["Protein Group"].nunique())),
             ("membrane", "highest protein group number in block 1", int(b1["Protein Group"].max())),
             ("membrane", "AGI loci in block 1", int(b1.locus.nunique())),
             ("membrane", "loci after the score, channel and two-peptide filters", int(INPUTS.psms("membrane").locus.nunique()))]
    return [{"kind": "accounting", "fraction": f, "item": i, "n": v} for f, i, v in rows]


def build():
    sets = pd.read_csv(PREREG / "osd38" / "gene_sets.tsv", sep="\t").groupby("set").gene.apply(set)
    t12 = pd.read_csv(TABLES / "T12_osd38_geneset_tests.tsv", sep="\t")
    t12 = t12[t12.layer == "membrane_flight"].set_index("set")
    lim = limma_membrane()
    u = [g for g in ORDER if g in lim.index]
    v = lim.loc[u, "t"].to_numpy(float)
    rng = np.random.default_rng(T16.SEED)
    rows = []
    for k in SETS:
        m = np.isin(u, list(sets[k]))
        so, _ = T16.set_stats(v, m)
        rs, _ = T16.rotation(v, m, rng)
        rows.append({"kind": "membrane set, limma moderated t", "set": k, "n": int(m.sum()), "signed": so,
                     "p": T16.p_two(rs, so, False), "z": (so - rs.mean()) / rs.std(),
                     "T12_p": t12.loc[k, "p"], "T12_z": t12.loc[k, "z"], "T12_observed": t12.loc[k, "observed"]})
    t = pd.DataFrame(rows)
    t["q"] = T16.bh(t.p)
    student = INPUTS.quantify().set_index(["fraction", "locus"])
    rows = [{"kind": "membrane, limma vs Student t", "item": "Spearman of protein t", "n": len(lim),
             "rho": stats.spearmanr(lim.t, student.loc["membrane"].t.reindex(lim.index)).statistic}]
    mem, sol = student.loc["membrane"], student.loc["soluble"]
    both = mem.index.intersection(sol.index)
    rows.append({"kind": "partition", "item": "proteins quantified in both fractions: Spearman of log2FC",
                 "n": len(both), "rho": stats.spearmanr(mem.log2fc[both], sol.log2fc[both]).statistic})
    for k in CARBON:
        g = sorted(set(both) & sets[k])
        if not g:
            rows.append({"kind": "partition", "item": f"{k}: members in both fractions", "n": 0})
            continue
        rows.append({"kind": "partition", "item": f"{k}: members in both fractions", "set": k, "n": len(g),
                     "membrane_mean_log2fc": mem.log2fc[g].mean(), "soluble_mean_log2fc": sol.log2fc[g].mean(),
                     "membrane_down": int((mem.log2fc[g] < 0).sum()), "soluble_down": int((sol.log2fc[g] < 0).sum()),
                     "members": ";".join(f"{x}:{mem.log2fc[x]:.2f}/{sol.log2fc[x]:.2f}" for x in g)})
    for k in CARBON:
        g = sorted(set(sol.index) & sets[k])
        rows.append({"kind": "soluble members", "item": f"{k}: members quantified in the soluble fraction", "set": k,
                     "n": len(g), "soluble_mean_log2fc": sol.log2fc[g].mean() if g else np.nan,
                     "soluble_down": int((sol.log2fc[g] < 0).sum())})
    t = pd.concat([t, pd.DataFrame(rows), pd.DataFrame(accounting())], ignore_index=True)
    if VERSIONS.exists():
        ver = pd.read_csv(VERSIONS, sep="\t")
        for r in ver.itertuples():
            t.loc[len(t)] = {"kind": "OSDR file listing", "item": f"{r.query}: {r.file_name}, {r.file_size} bytes, "
                             f"created {r.date_created_utc}, updated '{r.date_updated if isinstance(r.date_updated, str) else ''}'",
                             "n": r.files_listed}
    lead = ["kind", "set", "item", "n", "signed", "p", "q", "z", "T12_p", "T12_z", "T12_observed", "rho"]
    t = t[lead + [c for c in t.columns if c not in lead]]
    t.to_csv(TABLES / "T23_osd38_proteome_checks.tsv", sep="\t", index=False, float_format="%.4g")
    with pd.option_context("display.max_colwidth", 90, "display.width", 250):
        print(t.drop(columns=["members"], errors="ignore").to_string(index=False))


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("step", choices=["versions", "build"])
    {"versions": versions, "build": build}[ap.parse_args().step]()
