"""Can these gene sets see hypoxia at all? A positive control for the BRIC-PDFU comparison (T24).

Revision 1, post hoc: the settings are in prereg/revision1/settings.json (R4), committed before this
script ran. OSD-38's ground comparison of the BRIC-PDFU with a closed container found no hypoxia or
fermentation signal (T12, layer rna_hardware). That null only means something if the same sets, with
the same test, detect real hypoxia. Two hypoxia datasets:

  OSD-136      GeneLab's microarray DE table: Ws-0 roots at 50 kPa total pressure, without against with
               O2 restored to 21 kPa (limma t, GeneLab's)
  GSE14502     Mustroph et al. 2009, 7-day Col-0 (p35S:HF-RPL18B) seedlings, 2 h hypoxia against
               control, total mRNA; limma moderated t on the authors' RMA series matrix, shoot (3 v 3)
               and root (3 v 4) separately (scripts/r/limma_two_group.R)

Probes map to genes through GPL198 (GEO annotation; probes naming exactly one AGI; per AGI the probe
with the highest mean expression over the samples used) and, for OSD-136, GeneLab's own TAIR column
with the same one-mapping, highest-mean rule. Sets: hypoxia-responsive, fermentation and DIN. Genomic
rotation null as in T08; BH over the nine rows. Pass criterion, fixed beforehand: a set is
hypoxia-sensitive if it is up at q < 0.05 in at least one GSE14502 contrast, and the BRIC-vs-container
result for that set counts as informative only then.

T25 is the R4b addendum, written after R4 ran: an empirical hypoxia signature (genes more than 2-fold
up, or down, at q < 0.05 in both GSE14502 contrasts) tested in OSD-136 and in the OSD-38 hardware and
flight contrasts, with ADH1 and PDC1 gene by gene.
"""

import gzip
import io
import json
import subprocess

import numpy as np
import pandas as pd

import fetch_record
import genesets
from paths import CACHE, OSDR, PREREG, ROOT, TABLES

R4 = json.loads((PREREG / "revision1" / "settings.json").read_text())["R4_hypoxia_positive_control"]
GSE = CACHE / "gse14502"
OSD136_STAT = "T.stat_(50,000 pascal)v(50,000 pascal + O2 to 21,000 pascal)"
GSE_GROUPS = {"GSE14502_shoot": "shoot", "GSE14502_root": "root"}


def load(name):
    import importlib.util
    spec = importlib.util.spec_from_file_location(name.replace(".py", ""), ROOT / "scripts" / name)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


T16 = load("16_gwas_light_dark_tests.py")
ORDER = pd.read_csv(PREREG / "gene_order.tsv.gz", sep="\t").gene.tolist()
assert OSD136_STAT in R4["datasets"]["OSD-136"]


def sets():
    a = genesets.load()
    b = pd.read_csv(PREREG / "osd38" / "gene_sets.tsv", sep="\t").groupby("set", sort=False).gene.apply(set).to_dict()
    return {k: a[k] if k in a else b[k] for k in R4["sets"]}


def osd136():
    d = pd.read_csv(fetch_record.osdr("OSD-136", "GLDS-136_array_differential_expression_GLmicroarray.csv"),
                    usecols=["TAIR", "ProbesetID", "count_ENSEMBL_mappings", "All.mean", OSD136_STAT])
    d = d[(d.count_ENSEMBL_mappings == 1) & d.TAIR.str.fullmatch(r"AT[1-5CM]G\d{5}", na=False)]
    d = d.sort_values("All.mean", ascending=False).drop_duplicates("TAIR")
    return d.set_index("TAIR")[OSD136_STAT].rename("t")


def gse14502():
    m = fetch_record.get("https://ftp.ncbi.nlm.nih.gov/geo/series/GSE14nnn/GSE14502/matrix/GSE14502_series_matrix.txt.gz",
                         GSE / "GSE14502_series_matrix.txt.gz")
    a = fetch_record.get("https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPLnnn/GPL198/annot/GPL198.annot.gz",
                         GSE / "GPL198.annot.gz")
    lines = gzip.open(m, "rt").read().splitlines()
    titles = [t.strip('"') for t in next(x for x in lines if x.startswith("!Sample_title")).split("\t")[1:]]
    b = next(i for i, x in enumerate(lines) if x.startswith("!series_matrix_table_begin"))
    e = next(i for i, x in enumerate(lines) if x.startswith("!series_matrix_table_end"))
    x = pd.read_csv(io.StringIO("\n".join(lines[b + 1:e])), sep="\t", index_col=0)
    x.columns = titles
    al = gzip.open(a, "rt").read().splitlines()
    hb = next(i for i, x in enumerate(al) if x.startswith("!platform_table_begin"))
    ann = pd.read_csv(io.StringIO("\n".join(al[hb + 1:-1])), sep="\t", dtype=str).set_index("ID")
    agi = ann.Platform_ORF.fillna("").str.upper()
    agi = agi[agi.str.fullmatch(r"AT[1-5CM]G\d{5}")]
    out = {}
    work = GSE / "limma"
    work.mkdir(parents=True, exist_ok=True)
    for name, tissue in GSE_GROUPS.items():
        cols = [c for c in titles if c.startswith(f"{tissue}_") and c.endswith(tuple(f"_total_rep{i}" for i in range(1, 5)))]
        grp = ["treated" if "_hypoxia_" in c else "control" for c in cols]
        y = x.loc[x.index.intersection(agi.index), cols]
        y = y.assign(agi=agi[y.index], mean=y.mean(axis=1)).sort_values("mean", ascending=False)
        y = y.drop_duplicates("agi").set_index("agi")[cols]
        y.to_csv(work / f"{name}_expr.tsv", sep="\t")
        pd.DataFrame({"sample": cols, "group": grp}).to_csv(work / f"{name}_design.tsv", sep="\t", index=False)
        subprocess.run(["Rscript", str(ROOT / "scripts" / "r" / "limma_two_group.R"), str(work / f"{name}_expr.tsv"),
                        str(work / f"{name}_design.tsv"), str(work / f"{name}_limma.tsv")], check=True)
        r = pd.read_csv(work / f"{name}_limma.tsv", sep="\t").set_index("feature")
        r["q"] = T16.bh(r.P)
        out[name] = (r, f"{grp.count('treated')} v {grp.count('control')}")
    return out


def test(t, gs, rng):
    u = [g for g in ORDER if g in t.index]
    v = t.loc[u].to_numpy(float)
    rows = []
    for k, genes in gs.items():
        m = np.isin(u, list(genes))
        so, _ = T16.set_stats(v, m)
        rs, _ = T16.rotation(v, m, rng)
        rows.append({"set": k, "n_set": len(genes), "n_used": int(m.sum()), "signed": so,
                     "p": T16.p_two(rs, so, False), "z": (so - rs.mean()) / rs.std(),
                     "members": ";".join(f"{g}:{t[g]:.2f}" for g in sorted(genes) if g in t.index)})
    return rows, len(u)


def run():
    gs = sets()
    gse = gse14502()
    data = {"OSD-136": (osd136(), "3 v 3")} | {k: (r.t, n) for k, (r, n) in gse.items()}
    rows = []
    for name, (t, n) in data.items():
        r, genes = test(t, gs, np.random.default_rng(T16.SEED))
        rows += [x | {"dataset": name, "arms": n, "genes_tested": genes} for x in r]
    out = pd.DataFrame(rows)
    out["q"] = T16.bh(out.p)
    out["direction"] = np.where(out.q >= 0.05, "none", np.where(out.signed > 0, "up", "down"))
    passed = {k for k in gs if ((out.set == k) & out.dataset.isin(GSE_GROUPS) & (out.direction == "up")).any()}
    t12 = pd.read_csv(TABLES / "T12_osd38_geneset_tests.tsv", sep="\t")
    hw = t12[(t12.layer == "rna_hardware") & t12.set.isin(gs)]
    hw = pd.DataFrame({"dataset": "OSD-38 BRIC-PDFU vs closed container (T12)", "set": hw.set, "n_used": hw.n_measured,
                       "signed": hw.signed, "p": hw.p, "z": hw.z, "q": hw.q, "direction": hw.observed})
    out = pd.concat([out, hw], ignore_index=True)
    out["hypoxia_sensitive"] = out.set.isin(passed)
    out["role"] = np.where(out.dataset.str.startswith("OSD-38"),
                           np.where(out.hypoxia_sensitive, "informative", "not informative"), "positive control")
    lead = ["dataset", "set", "role", "hypoxia_sensitive", "direction", "p", "q", "z", "signed", "n_set", "n_used", "arms"]
    out = out[lead + [c for c in out.columns if c not in lead]]
    out[["n_set", "n_used", "genes_tested"]] = out[["n_set", "n_used", "genes_tested"]].astype("Int64")
    out.to_csv(TABLES / "T24_hypoxia_positive_control.tsv", sep="\t", index=False, float_format="%.4g")
    print(out[["dataset", "set", "role", "direction", "p", "q", "z", "n_used"]].to_string(index=False))
    signature(gse, data["OSD-136"][0])
    return out


MARKERS = {"AT1G77120": "ADH1", "AT4G33070": "PDC1"}


def signature(gse, t136):
    """R4b (addendum, written after R4 ran): an empirical hypoxia signature from GSE14502, tested in
    OSD-136 and in the OSD-38 hardware and flight contrasts (T25)."""
    r4b = json.loads((PREREG / "revision1" / "settings.json").read_text())["R4b_hypoxia_signature"]
    assert "both the shoot and the root" in r4b["signature"].lower()
    sh, rt = (gse[k][0] for k in GSE_GROUPS)
    both = sh.index.intersection(rt.index)
    sig = {"hypoxia_up": {g for g in both if sh.logFC[g] > 1 and rt.logFC[g] > 1 and sh.q[g] < 0.05 and rt.q[g] < 0.05},
           "hypoxia_down": {g for g in both if sh.logFC[g] < -1 and rt.logFC[g] < -1 and sh.q[g] < 0.05 and rt.q[g] < 0.05}}
    targets = {"OSD-136": t136}
    for c, lab in (("OSD-38-hardware", "OSD-38 BRIC-PDFU vs closed container"), ("OSD-38", "OSD-38 flight vs ground")):
        d = pd.read_csv(OSDR / f"contrast_{c}.csv.gz").dropna(subset=["TAIR", "stat"]).drop_duplicates("TAIR").set_index("TAIR")
        targets[lab] = d
    rows = []
    for lab, d in targets.items():
        t = d if isinstance(d, pd.Series) else d.stat
        rr, genes = test(t, sig, np.random.default_rng(T16.SEED))
        rows += [{"kind": "signature", "dataset": lab, "genes_tested": genes, **{k: v for k, v in x.items() if k != "members"}}
                 for x in rr]
        if not isinstance(d, pd.Series):
            for g, name in MARKERS.items():
                rows.append({"kind": "gene", "dataset": lab, "set": name, "n_set": 1, "n_used": 1,
                             "signed": d.stat[g], "log2fc": d.log2fc[g], "padj": d.padj[g]})
    for name, (r, n) in gse.items():
        for g, m in MARKERS.items():
            rows.append({"kind": "gene", "dataset": name, "set": m, "n_set": 1, "n_used": 1, "signed": r.t[g],
                         "log2fc": r.logFC[g], "padj": r.q[g]})
    out = pd.DataFrame(rows)
    sel = out.kind == "signature"
    out.loc[sel, "q"] = T16.bh(out.loc[sel, "p"])
    out.loc[sel, "direction"] = np.where(out.loc[sel, "q"] >= 0.05, "none",
                                         np.where(out.loc[sel, "signed"] > 0, "up", "down"))
    out[["n_set", "n_used", "genes_tested"]] = out[["n_set", "n_used", "genes_tested"]].astype("Int64")
    out.to_csv(TABLES / "T25_hypoxia_signature.tsv", sep="\t", index=False, float_format="%.4g")
    print(f"  R4b signature: {len(sig['hypoxia_up'])} hypoxia-up and {len(sig['hypoxia_down'])} hypoxia-down genes "
          f"(GSE14502 shoot and root, |log2FC| > 1, q < 0.05)")
    print(out[["kind", "dataset", "set", "n_used", "signed", "p", "q", "z", "direction", "log2fc", "padj"]].to_string(index=False))
    return out


if __name__ == "__main__":
    run()
