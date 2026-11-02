"""What drives the PaintOmics results: layer ablations, DEG-count matching and a hub-overlap null (T17, T18).

Revision 1, post hoc: the jobs and comparisons are in prereg/revision1/settings.json (R6), committed
before the jobs were submitted. Jobs A-D and T04-T07 are unchanged; the control jobs are compared with
their parents here.

  E  A without the metabolite layer          F  A with effect-scaled metabolite values
  G  A genes + layer (no proteins)           H  A proteins + layer (no genes)
  I  D genes + layer (no proteins)           J  D membrane proteins + layer (no genes)
  K1, K2  B with its relevant genes downsampled to C's count (seeds 1, 2)

T17 (long format): per control job, pathways at BH q < 0.05 against the parent, Spearman of combined p,
hubs at FDR < 0.05 against the parent, and whether raw hub p-values match the parent's; plus how the
combined pathway p tracks each omic's p in job A, which relevant compounds never reached the hub test,
and the BH step-up ties among job B's hubs.
T18: the overlap of the lit (A) and dark (D) hub sets against a hypergeometric null over the compounds
both jobs tested; each compound's KEGG reaction count (REST link/reaction, cached); a Mann-Whitney test
of reaction count in hubs against non-hubs per job; and the overlap under hub labels shuffled within
reaction-count tertiles.
"""

import gzip
import json

import numpy as np
import pandas as pd
from scipy import stats

from genesets import kegg
from paths import DATA, PAINTOMICS_RAW, PREREG, TABLES

R6 = json.loads((PREREG / "revision1" / "settings.json").read_text())["R6_paintomics"]
ALPHA = 0.05
IDS = {"A": "6kBOS7AN5K", "B": "JfVZ1L37PY", "C": "wgT3IPi7U2", "D": "Es6j23q40o", "E": "K3ju047ZWB",
       "F": "L46ds55AiY", "G": "Y72ky4006o", "H": "N2XgxyWlJn", "I": "sD7V1nZ4t4", "J": "C3kVBAm3N5",
       "K1": "Y0V56yFW12", "K2": "7725j3CkOp"}
PARENT = {"E": "A", "F": "A", "G": "A", "H": "A", "I": "D", "J": "D", "K1": "B", "K2": "B"}
assert set(PARENT) == set(R6["jobs"])
DISPLAY = {"C02504": "2-Isopropylmalate"}   # PaintOmics shows KEGG's first name, alpha-Isopropylmalate
SEED, PERMS = 20261102, 10_000


def load(tag):
    with gzip.open(PAINTOMICS_RAW / f"job{tag}_{IDS[tag]}_full.json.gz", "rt") as fh:
        return json.load(fh)


def bh(p):
    p = np.asarray(p, float)
    o = np.argsort(p)
    q = p[o] * len(p) / np.arange(1, len(p) + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    out = np.empty_like(q)
    out[o] = np.minimum(q, 1)
    return out


def pathways(job):
    rows = [{"pathway": (("KEGG" if p.get("source") == "KEGG" else "MapMan"), p["name"]),
             "p": p["combinedSignificancePvalues"]["Fisher"][0], **(p.get("globalOmicPvalues") or {})}
            for p in job["pathwaysInfo"]]
    d = pd.DataFrame(rows).set_index("pathway")
    d["q"] = bh(d.p)
    return d


def hubs(job):
    h = pd.DataFrame((job.get("hubAnalysisResult") or {}).values())
    if h.empty:
        return h, pd.DataFrame(columns=["pvalue_adjust", "pvalue", "step"])
    return h, h.sort_values(["pvalue_adjust", "pvalue"]).groupby("name").first()[["pvalue_adjust", "pvalue", "step"]]


def names(jobs):
    out = {cid: n for j in jobs.values() for cid, n in (j.get("mappingComp") or {}).items()}
    return out | DISPLAY


def degree(cids):
    """KEGG reactions per compound (rest.kegg.jp/link/reaction, cached in data/cache/kegg)."""
    cids = sorted(cids)
    links = []
    for i in range(0, len(cids), 10):
        chunk = cids[i:i + 10]
        links += kegg("link/reaction/" + "+".join(f"cpd:{c}" for c in chunk),
                      f"link_reaction_{chunk[0]}_{chunk[-1]}.tsv").splitlines()
    n = pd.Series([ln.split("\t")[0].replace("cpd:", "") for ln in links if "\t" in ln]).value_counts()
    return n.reindex(cids).fillna(0).astype(int)


def t17(jobs, nm):
    P = {k: pathways(j) for k, j in jobs.items()}
    H = {k: hubs(j) for k, j in jobs.items()}
    sig = {k: set(p.index[p.q < ALPHA]) for k, p in P.items()}
    hubsig = {k: set(b.index[b.pvalue_adjust < ALPHA]) for k, (_, b) in H.items()}
    lab = lambda s: "; ".join(sorted(nm.get(c, c) for c in s)) or "none"
    pw = lambda s: "; ".join(sorted(f"{d}: {n}" for d, n in s)) or "none"
    rows = []

    def add(job, versus, quantity, value):
        rows.append({"job": job, "job_id": IDS.get(job, ""), "versus": versus, "quantity": quantity, "value": value})

    for c, par in PARENT.items():
        add(c, par, "what", R6["jobs"][c])
        add(c, par, "pathways at BH q < 0.05 (job; parent)", f"{len(sig[c])}; {len(sig[par])}")
        add(c, par, "pathways gained", pw(sig[c] - sig[par]))
        add(c, par, "pathways lost", pw(sig[par] - sig[c]))
        common = P[c].index.intersection(P[par].index)
        add(c, par, "Spearman rho of combined p", round(stats.spearmanr(P[c].p[common], P[par].p[common]).statistic, 4))
        add(c, par, "hub compounds tested (job; parent)", f"{len(H[c][1])}; {len(H[par][1])}")
        add(c, par, "hubs at FDR < 0.05 (job; parent)", f"{len(hubsig[c])}; {len(hubsig[par])}")
        add(c, par, "hubs gained", lab(hubsig[c] - hubsig[par]))
        add(c, par, "hubs lost", lab(hubsig[par] - hubsig[c]))
        if not H[c][0].empty:
            sh = (pd.concat([H[c][0].assign(job="x"), H[par][0].assign(job="y")])
                  .pivot_table(index=["name", "step"], columns="job", values="pvalue").dropna())
            same = bool(np.allclose(sh.x, sh.y, rtol=0, atol=1e-12)) if len(sh) else None
            add(c, par, "raw hub p identical to the parent for shared compound-steps", f"{same} ({len(sh)} compound-steps)")
    for k in ("K1", "K2"):
        add(k, "C", "hubs at FDR < 0.05 (job; C)", f"{len(hubsig[k])}; {len(hubsig['C'])}")
        add(k, "C", "hubs shared with C", lab(hubsig[k] & hubsig["C"]))
        add(k, "C", "pathways at BH q < 0.05 (job; C)", f"{len(sig[k])}; {len(sig['C'])}")

    a = P["A"]
    for omic in ("Gene expression", "Proteomics", "Metabolomics"):
        x = a[omic].dropna()
        add("A", "", f"pathways with a {omic} p", len(x))
        add("A", "", f"Spearman rho of combined p with the {omic} p", round(stats.spearmanr(a.p[x.index], x).statistic, 4))
    key = pd.read_csv(DATA / "paintomics_upload" / "metabolomics_key.tsv", sep="\t")
    relevant = set(key.kegg_id[key.relevant])
    for k in ("A", "B", "C", "D"):
        tested = set(H[k][1].index)
        add(k, "", "relevant compounds (robust) in the layer", len(relevant))
        add(k, "", "relevant compounds hub-tested", len(relevant & tested))
        add(k, "", "relevant compounds with no hub result", lab(relevant - tested))
    b = H["B"][1]
    s = b[b.pvalue_adjust < ALPHA].sort_values("pvalue")
    add("B", "", "hubs at FDR < 0.05: raw p range", f"{s.pvalue.min():.4g}-{s.pvalue.max():.4g}")
    add("B", "", "hubs at FDR < 0.05: distinct adjusted p values", "; ".join(f"{v:.4g}" for v in sorted(s.pvalue_adjust.unique())))
    t = pd.DataFrame(rows)
    t.to_csv(TABLES / "T17_paintomics_controls.tsv", sep="\t", index=False)
    with pd.option_context("display.max_colwidth", 140, "display.width", 250):
        print(t.to_string(index=False))
    return H, hubsig


def t18(H, hubsig, nm):
    rng = np.random.default_rng(SEED)
    rows = []
    tested = {k: set(H[k][1].index) for k in ("A", "B", "C", "D")}
    deg = degree(set().union(*tested.values()))
    cand = sorted(tested["A"] & tested["D"])
    a, d = hubsig["A"] & set(cand), hubsig["D"] & set(cand)
    obs, N = len(a & d), len(cand)
    hg = stats.hypergeom(N, len(a), len(d))
    rows.append({"kind": "overlap A-D", "test": "hypergeometric over compounds hub-tested in both", "n": N,
                 "hubs_A": len(a), "hubs_D": len(d), "observed": obs, "expected": len(a) * len(d) / N,
                 "p": hg.sf(obs - 1)})
    tert = pd.qcut(deg[cand].rank(method="first"), 3, labels=False).to_numpy()
    la, ld = np.isin(cand, list(a)), np.isin(cand, list(d))
    null = np.empty(PERMS, int)
    for i in range(PERMS):
        pa, pd_ = la.copy(), ld.copy()
        for t in range(3):
            m = tert == t
            pa[m], pd_[m] = rng.permutation(la[m]), rng.permutation(ld[m])
        null[i] = int((pa & pd_).sum())
    rows.append({"kind": "overlap A-D", "test": "hub labels shuffled within KEGG reaction-count tertiles", "n": N,
                 "hubs_A": len(a), "hubs_D": len(d), "observed": obs, "expected": null.mean(),
                 "p": (1 + np.sum(null >= obs)) / (PERMS + 1)})
    for k in ("A", "B", "C", "D"):
        hub, non = [c for c in tested[k] if c in hubsig[k]], [c for c in tested[k] if c not in hubsig[k]]
        rows.append({"kind": "reaction count", "test": f"job {k}: Mann-Whitney, hubs vs other hub-tested compounds",
                     "n": len(tested[k]), "hubs": len(hub), "median_reactions_hubs": float(deg[hub].median()),
                     "median_reactions_others": float(deg[non].median()),
                     "p": stats.mannwhitneyu(deg[hub], deg[non], alternative="two-sided").pvalue if hub and non else np.nan})
    rows.append({"kind": "reaction count", "test": "Spearman of reaction count with -log10 best hub p, job A",
                 "n": len(tested["A"]), "rho": stats.spearmanr(deg[sorted(tested["A"])],
                                                               -np.log10(H["A"][1].pvalue[sorted(tested["A"])])).statistic})
    for c in sorted(deg.index, key=lambda c: -deg[c]):
        rows.append({"kind": "compound", "kegg_id": c, "compound": nm.get(c, c), "reactions": int(deg[c]),
                     **{f"hub_{k}": ("hub" if c in hubsig[k] else "tested" if c in tested[k] else "")
                        for k in ("A", "B", "C", "D")}})
    t = pd.DataFrame(rows)
    t.to_csv(TABLES / "T18_hub_overlap_null.tsv", sep="\t", index=False, float_format="%.4g")
    with pd.option_context("display.width", 250):
        print(t[t.kind != "compound"].dropna(axis=1, how="all").to_string(index=False))
    return t


if __name__ == "__main__":
    jobs = {k: load(k) for k in IDS}
    nm = names(jobs)
    H, hubsig = t17(jobs, nm)
    t18(H, hubsig, nm)
