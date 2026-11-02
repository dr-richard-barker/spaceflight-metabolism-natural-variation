"""OSD-678 Col-0: is the flight response different in light and in darkness? (T19)

Revision 1, post hoc: the settings are in prereg/revision1/settings.json (R3), committed before this
script ran. The pre-registered tests treat the OSD-678 light and dark arms as two separate contrasts
(T08); this asks the reviewer's question directly, with one model for the 12 Col-0 samples:
edgeR TMM + limma-voom, design flight x light (scripts/r/limma_voom_contrasts.R). The interaction
moderated t, (flight - ground in light) - (flight - ground in dark), is tested with the genomic
rotation null for the eleven sets, BH over the family. The flight-in-light and flight-in-dark
statistics from the same fit are tested too, for context, with their own BH.

Set statistic: signed (mean t in the set minus the mean over all genes) for the KEGG, DIN and marker
sets; unsigned for the two Wu GWAS sets, as in their pre-registration (prereg/settings.json).
"""

import json
import subprocess
import sys

import numpy as np
import pandas as pd
from scipy import stats

import fetch_record
import genesets
from paths import CACHE, PREREG, ROOT, TABLES

R3 = json.loads((PREREG / "revision1" / "settings.json").read_text())["R3_osd678_interaction"]
COUNTS = "GLDS-612_rna_seq_RSEM_Unnormalized_Counts_rRNArm_GLbulkRNAseq.csv"
WORK = CACHE / "osd678" / "limma"


def load(name):
    import importlib.util
    spec = importlib.util.spec_from_file_location(name.replace(".py", ""), ROOT / "scripts" / name)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


T16 = load("16_gwas_light_dark_tests.py")
ORDER = pd.read_csv(PREREG / "gene_order.tsv.gz", sep="\t").gene.tolist()


def sets():
    a = genesets.load()
    b = pd.read_csv(PREREG / "osd38" / "gene_sets.tsv", sep="\t").groupby("set", sort=False).gene.apply(set).to_dict()
    out = {}
    for k in R3["sets"]:
        out[k] = a[k] if k in a else b[k]
    return out


def fit():
    counts = pd.read_csv(fetch_record.osdr("OSD-678", COUNTS), index_col=0)
    col = [c for c in counts.columns if "_Col-0_" in c]
    assert len(col) == 12, col
    col = sorted(col)
    des = pd.DataFrame({"sample": col,
                        "flight": ["flight" if c.startswith("Flight_") else "ground" for c in col],
                        "light": ["light" if "_Light_" in c else "dark" for c in col]})
    assert (des.groupby(["flight", "light"]).size() == 3).all()
    x = counts[col]
    x = x[x.sum(axis=1) >= T16.MIN_COUNT]
    WORK.mkdir(parents=True, exist_ok=True)
    x.to_csv(WORK / "counts.tsv", sep="\t")
    des.to_csv(WORK / "design.tsv", sep="\t", index=False)
    subprocess.run(["Rscript", str(ROOT / "scripts" / "r" / "limma_voom_contrasts.R"), str(WORK / "counts.tsv"),
                    str(WORK / "design.tsv"), str(WORK / "limma.tsv")], check=True)
    return pd.read_csv(WORK / "limma.tsv", sep="\t").set_index("gene")


def check(res):
    out = {}
    for arm in ("light", "dark"):
        gl = T16.contrast(f"OSD-678-{arm}")
        common = res.index.intersection(gl.index)
        out[arm] = (stats.spearmanr(res.loc[common, f"flight_{arm}_t"], gl.loc[common, "stat"]).statistic, len(common))
        print(f"  limma flight-in-{arm} t vs GeneLab Wald: Spearman {out[arm][0]:.3f} over {len(common)} genes")
    return out


def selftest(v):
    rng = np.random.default_rng(1)
    v = v.copy()
    planted, control = np.zeros(len(v), bool), np.zeros(len(v), bool)
    planted[rng.choice(len(v), 30, replace=False)] = True
    control[rng.choice(np.flatnonzero(~planted), 30, replace=False)] = True
    v[planted] += 2.0
    pp = T16.p_two(T16.rotation(v, planted, rng)[0], T16.set_stats(v, planted)[0], False)
    pc = T16.p_two(T16.rotation(v, control, rng)[0], T16.set_stats(v, control)[0], False)
    ok = pp < 0.01 and pc > 0.05
    print(f"  selftest (interaction t): planted p={pp:.4f}, control p={pc:.3f} -> {'PASS' if ok else 'FAIL'}")
    if not ok:
        sys.exit("selftest failed")


def tests(res, gs, rho):
    unsigned = {k for k, s in T16.S["primary_statistic"].items() if s == "unsigned"}
    u = [g for g in ORDER if g in res.index]
    rows = []
    for term in ("interaction", "flight_light", "flight_dark"):
        v = res.loc[u, f"{term}_t"].to_numpy(float)
        if term == "interaction":
            selftest(v)
        rng = np.random.default_rng(T16.SEED)
        for k, genes in gs.items():
            m = np.isin(u, list(genes))
            so, uo = T16.set_stats(v, m)
            rs, ru = T16.rotation(v, m, rng)
            st = "unsigned" if k in unsigned else "signed"
            p = T16.p_one(ru, uo, False) if st == "unsigned" else T16.p_two(rs, so, False)
            z = (uo - ru.mean()) / ru.std() if st == "unsigned" else (so - rs.mean()) / rs.std()
            rows.append({"term": term, "role": "primary" if term == "interaction" else "context", "set": k,
                         "statistic": st, "n_set": len(genes), "n_used": int(m.sum()), "signed": so,
                         "unsigned": uo, "p": p, "z": z})
    t = pd.DataFrame(rows)
    t["q"] = np.nan
    for term in t.term.unique():
        sel = t.term == term
        t.loc[sel, "q"] = T16.bh(t.loc[sel, "p"])
    t["direction"] = np.where(t.q >= 0.05, "none",
                              np.where(t.statistic == "unsigned", np.where(t.unsigned > 0, "enriched", "depleted"),
                                       np.where(t.signed > 0, "up", "down")))
    for arm in ("light", "dark"):
        t[f"check_spearman_{arm}"] = rho[arm][0]
    t["genes_tested"] = len(u)
    t.to_csv(TABLES / "T19_osd678_interaction.tsv", sep="\t", index=False, float_format="%.4g")
    print(t.pivot(index="set", columns="term", values="z").round(2).to_string())
    print(t[t.term == "interaction"][["set", "statistic", "n_used", "p", "q", "z", "direction"]].to_string(index=False))
    return t


if __name__ == "__main__":
    res = fit()
    rho = check(res)
    tests(res, sets(), rho)
