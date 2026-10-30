"""OSD-38 (BRIC-20) tests, read from the frozen prereg/osd38/ (T11, T12, T13, fig06).

  T11  protein log2 fold changes, flight vs ground, per fraction (method in prereg/osd38/settings.json)
  T12  the pre-registered gene-set tests: nine sets in four layers (OSD-38 flight transcriptome, membrane
       and soluble proteomes, and the ground BRIC-PDFU vs closed-container transcriptome)
  T13  descriptive comparisons: transcript-protein concordance (OSD-38 and OSD-522), the paper's own DE
       count as a check on the quantification, and (exploratory, not pre-registered) how the flight
       transcriptome relates to the ground hardware effect

Nulls reuse scripts/16_gwas_light_dark_tests.py (genomic rotation, same seed and draws). Before testing,
a shift is planted into 15 measured membrane proteins and 15 untouched proteins kept as a control; the
script stops unless the rotation finds the first and not the second.
"""

import importlib.util
import json
import sys

import numpy as np
import pandas as pd
from scipy import stats

import style
from genesets import kegg
from paths import DATA, FIGURES, OSDR, PHOTORESP, PREREG, ROOT, TABLES

P38 = PREREG / "osd38"
if not (P38 / "settings.json").exists():
    sys.exit("prereg/osd38/settings.json missing: run scripts/21_prereg_osd38.py and commit it first")
S = json.loads((P38 / "settings.json").read_text())


def load(name):
    spec = importlib.util.spec_from_file_location(name.replace(".py", ""), ROOT / "scripts" / name)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


T16 = load("16_gwas_light_dark_tests.py")
INPUTS = load("20_osd38_inputs.py")
ORDER = pd.read_csv(PREREG / "gene_order.tsv.gz", sep="\t").gene.tolist()
MIN = S["min_measured"]
PRETTY = {"din_first_study": "DIN starvation markers", "bcaa_degradation": "BCAA degradation (ath00280)",
          "photosynthesis_apparatus": "Photosystem & light harvesting", "fermentation": "Fermentation",
          "hypoxia_responsive": "Hypoxia-responsive", "ath00630": "Glyoxylate & dicarboxylate (ath00630)",
          "ath00710": "Carbon fixation (ath00710)", "ath00500": "Starch & sucrose (ath00500)",
          "ath00010": "Glycolysis / gluconeogenesis (ath00010)"}


def layers(prot):
    rna = pd.read_csv(OSDR / "contrast_OSD-38.csv.gz", usecols=["TAIR", "log2fc", "stat"]).dropna()
    hw = pd.read_csv(OSDR / "contrast_OSD-38-hardware.csv.gz", usecols=["TAIR", "log2fc", "stat"]).dropna()
    out = {"rna_flight": rna.set_index("TAIR").stat, "rna_hardware": hw.set_index("TAIR").stat}
    for frac in ("membrane", "soluble"):
        p = prot[prot.fraction == frac].set_index("locus")
        out[f"{frac}_flight"] = p.t
    return out


def ordered(values):
    u = [g for g in ORDER if g in values.index]
    return u, values.loc[u].to_numpy()


def selftest(values):
    rng = np.random.default_rng(1)
    u, v = ordered(values)
    v = v.copy()
    planted, control = np.zeros(len(v), bool), np.zeros(len(v), bool)
    planted[rng.choice(len(v), 15, replace=False)] = True
    control[rng.choice(np.flatnonzero(~planted), 15, replace=False)] = True
    v[planted] += 3.0
    pp = T16.p_two(T16.rotation(v, planted, rng)[0], T16.set_stats(v, planted)[0], False)
    pc = T16.p_two(T16.rotation(v, control, rng)[0], T16.set_stats(v, control)[0], False)
    ok = pp < 0.01 and pc > 0.05
    print(f"  selftest (membrane proteins): planted p={pp:.4f}, control p={pc:.3f} -> {'PASS' if ok else 'FAIL'}")
    if not ok:
        sys.exit("selftest failed")


def tests(lay):
    sets = pd.read_csv(P38 / "gene_sets.tsv", sep="\t").groupby("set", sort=False).gene.apply(set)
    pred = pd.read_csv(P38 / "PREDICTIONS.tsv", sep="\t")
    rows = []
    for layer, values in lay.items():
        u, v = ordered(values)
        rng = np.random.default_rng(T16.SEED)
        for name in pred[pred.layer == layer].set:
            m = np.isin(u, list(sets[name]))
            row = {"layer": layer, "set": name, "n_measured": int(m.sum())}
            if m.sum() >= MIN:
                obs = T16.set_stats(v, m)[0]
                null = T16.rotation(v, m, rng)[0]
                row.update({"signed": obs, "p": T16.p_two(null, obs, False), "z": (obs - null.mean()) / null.std()})
            rows.append(row)
    t = pd.DataFrame(rows).merge(pred[["layer", "set", "predicted", "primary"]], on=["layer", "set"])
    t["q"] = np.nan
    for primary in (True, False):
        sel = t.primary.eq(primary) & t.p.notna()
        t.loc[sel, "q"] = T16.bh(t.loc[sel, "p"])
    t["observed"] = np.where(t.p.isna(), "too few", np.where(t.q < 0.05, np.where(t.signed > 0, "up", "down"), "none"))
    t["verdict"] = np.where(t.p.isna(), "n/a", np.where(t.observed == t.predicted, "MATCH", "MISMATCH"))
    cols = ["layer", "set", "primary", "predicted", "observed", "verdict", "p", "q", "z", "signed", "n_measured"]
    t = t[cols]
    t.to_csv(TABLES / "T12_osd38_geneset_tests.tsv", sep="\t", index=False, float_format="%.4g")
    print(t.to_string(index=False))
    return t


def uniprot_to_agi():
    conv = {}
    for ln in kegg("conv/uniprot/ath", "conv_uniprot_ath.tsv").splitlines():
        a, u = ln.split("\t")
        conv.setdefault(u.replace("up:", ""), set()).add(a.replace("ath:", ""))
    return {u: next(iter(a)) for u, a in conv.items() if len(a) == 1}


def robustness(prot, draws=5000):
    """Not pre-registered: membrane sets against random sets matched on PSM-count decile, because iTRAQ
    biases can track abundance and the carbon-pathway enzymes are abundant proteins."""
    mem = prot[prot.fraction == "membrane"].set_index("locus")
    sets = pd.read_csv(P38 / "gene_sets.tsv", sep="\t").groupby("set").gene.apply(set)
    dec = pd.qcut(np.log(mem.psms), 10, labels=False, duplicates="drop").to_numpy()
    t, rng, rows = mem.t.to_numpy(), np.random.default_rng(T16.SEED), []
    rows.append({"comparison": "OSD-38 membrane: Spearman of protein t with log PSM count",
                 "kind": "robustness, not pre-registered", "n": len(mem),
                 "rho": stats.spearmanr(t, np.log(mem.psms)).correlation, "p": np.nan})
    for name in ("ath00630", "ath00710", "ath00500", "ath00010", "bcaa_degradation"):
        m = mem.index.isin(list(sets[name]))
        bins = pd.Series(dec[m]).value_counts()
        null = np.array([t[np.concatenate([rng.choice(np.flatnonzero(dec == b), n, replace=False)
                                           for b, n in bins.items()])].mean() for _ in range(draws)])
        obs = t[m].mean()
        p = (1 + np.sum(np.abs(null - null.mean()) >= abs(obs - null.mean()))) / (draws + 1)
        rows.append({"comparison": f"OSD-38 membrane {name}: mean t against PSM-decile-matched random sets "
                                   f"(observed {obs:.2f}, null mean {null.mean():.2f})",
                     "kind": "robustness, not pre-registered", "n": int(m.sum()), "rho": np.nan, "p": p})
    return rows


def concordance(prot, lay):
    rows = []
    rna38 = pd.read_csv(OSDR / "contrast_OSD-38.csv.gz", usecols=["TAIR", "log2fc"]).dropna().set_index("TAIR").log2fc
    for frac in ("membrane", "soluble"):
        p = prot[prot.fraction == frac].set_index("locus").log2fc
        common = p.index.intersection(rna38.index)
        r = stats.spearmanr(p[common], rna38[common])
        rows.append({"comparison": f"OSD-38 {frac} protein vs transcript log2FC", "kind": "pre-registered, descriptive",
                     "n": len(common), "rho": r.correlation, "p": r.pvalue})
    conv = uniprot_to_agi()
    p522 = pd.read_csv(PHOTORESP / "proteomics_values.tab", sep="\t", index_col=0).iloc[:, 0]
    p522.index = p522.index.map(conv)
    p522 = p522[p522.index.notna()].groupby(level=0).mean()
    rna522 = pd.read_csv(PHOTORESP / "gene_expression_values.tab", sep="\t", index_col=0).iloc[:, 0]
    common = p522.index.intersection(rna522.index)
    r = stats.spearmanr(p522[common], rna522[common])
    rows.append({"comparison": "OSD-522 protein vs transcript log2FC (first study's upload)",
                 "kind": "pre-registered, descriptive", "n": len(common), "rho": r.correlation, "p": r.pvalue})
    mem = prot[prot.fraction == "membrane"]
    n_rule = int(((mem.p <= 0.05) & (mem.log2fc.abs() > 0.2)).sum())
    rows.append({"comparison": "OSD-38 membrane proteins at the paper's rule (p <= 0.05, |log2FC| > 0.2); "
                               "Kruse et al. 2020 report 103", "kind": "quantification check", "n": n_rule,
                 "rho": np.nan, "p": np.nan})
    rows += robustness(prot)
    fl, hw = lay["rna_flight"], lay["rna_hardware"]
    common = fl.index.intersection(hw.index)
    r = stats.spearmanr(fl[common], hw[common])
    rows.append({"comparison": "OSD-38 flight vs ground (Wald) against BRIC vs container (Wald), all genes",
                 "kind": "exploratory, not pre-registered", "n": len(common), "rho": r.correlation, "p": r.pvalue})
    c = pd.DataFrame(rows)
    c.to_csv(TABLES / "T13_osd38_concordance.tsv", sep="\t", index=False, float_format="%.4g")
    print(c.to_string(index=False))
    return c, p522, rna522, rna38


def figure(t, prot, rna38, lay):
    style.apply()
    import matplotlib.pyplot as plt
    from matplotlib.colors import TwoSlopeNorm

    order_l = ["rna_flight", "membrane_flight", "soluble_flight", "rna_hardware"]
    label_l = ["Transcripts\nflight vs ground", "Membrane proteins\nflight vs ground",
               "Soluble proteins*\nflight vs ground", "Transcripts\nBRIC vs container"]
    sets = list(dict.fromkeys(t.set))
    fig = plt.figure(figsize=(7.2, 6.4))
    gs = fig.add_gridspec(2, 2, height_ratios=[1.25, 1], hspace=0.55, wspace=0.35)
    a = fig.add_subplot(gs[0, :])
    z = t.pivot(index="set", columns="layer", values="z").reindex(index=sets, columns=order_l)
    lim = max(4.0, float(np.nanmax(np.abs(z.to_numpy()))))
    a.imshow(z.to_numpy(), cmap="RdBu_r", norm=TwoSlopeNorm(0, -lim, lim), aspect="auto")
    for i, s in enumerate(sets):
        for j, l in enumerate(order_l):
            r = t[(t.set == s) & (t.layer == l)].iloc[0]
            if r.observed == "too few":
                a.text(j, i, f"n = {r.n_measured}", ha="center", va="center", fontsize=6, color=style.GREY)
                continue
            ink = "white" if abs(r.z) > lim * 0.55 else style.INK
            a.text(j, i - 0.12, r.observed, ha="center", va="center", fontsize=6.6,
                   fontweight="bold" if r.q < 0.05 else "normal", color=ink)
            a.text(j, i + 0.24, f"{'=' if r.verdict == 'MATCH' else '≠'} pred. {r.predicted}", ha="center",
                   va="center", fontsize=5.6, color=ink)
            if r.verdict == "MISMATCH" and r.q < 0.05:
                a.add_patch(plt.Rectangle((j - 0.47, i - 0.46), 0.94, 0.92, fill=False, lw=1.3, ec=style.INK))
    a.set_xticks(range(len(order_l)))
    a.set_xticklabels(label_l, fontsize=6.8)
    a.set_yticks(range(len(sets)))
    a.set_yticklabels([PRETTY.get(x, x) for x in sets], fontsize=6.8)
    a.grid(False)
    style.panel(a, "a", "Pre-registered OSD-38 tests (z against the rotation null; bold: q < 0.05)")

    b = fig.add_subplot(gs[1, 0])
    mem = prot[prot.fraction == "membrane"].set_index("locus").log2fc
    c38 = mem.index.intersection(rna38.index)
    sets_ = pd.read_csv(P38 / "gene_sets.tsv", sep="\t").groupby("set").gene.apply(set)
    carbon = set().union(*(sets_[k] for k in ("ath00630", "ath00710", "ath00500", "ath00010")))
    inc = c38[c38.isin(list(carbon))]
    b.scatter(rna38[c38], mem[c38], s=4, color=style.GREY_LIGHT, lw=0, label=f"all membrane proteins ({len(c38)})")
    b.scatter(rna38[inc], mem[inc], s=12, color=style.VERMILION, lw=0, label=f"four carbon pathways ({len(inc)})")
    b.axhline(0, color=style.GREY, lw=0.6)
    b.axvline(0, color=style.GREY, lw=0.6)
    b.set_xlim(-4, 4)
    b.set_ylim(-1.6, 1.6)
    b.set_xlabel("Transcript log2FC, flight vs ground")
    b.set_ylabel("Membrane protein log2FC")
    b.legend(loc="upper left", fontsize=6, markerscale=1.5, frameon=True, framealpha=0.9, edgecolor="none")
    style.panel(b, "b", "OSD-38 transcripts against proteins")

    c = fig.add_subplot(gs[1, 1])
    fl, hw = lay["rna_flight"], lay["rna_hardware"]
    common = fl.index.intersection(hw.index)
    c.hexbin(hw[common].clip(-15, 15), fl[common].clip(-15, 15), gridsize=45, bins="log", cmap="Greys", mincnt=1)
    c.axhline(0, color=style.GREY, lw=0.6)
    c.axvline(0, color=style.GREY, lw=0.6)
    c.set_xlabel("BRIC vs container, ground (Wald)")
    c.set_ylabel("Flight vs ground (Wald)")
    c.grid(False)
    style.panel(c, "c", "Flight against the hardware effect*")
    fig.text(0.01, 0.005, "* Soluble proteome truncated to its 437 top-ranked groups in OSDR. Panel c is exploratory "
                          "(not pre-registered). OSD-38: 3-day dark-grown Col-0 seedlings, BRIC-PDFU.",
             fontsize=6.2, color=style.GREY)
    style.save(fig, "fig06_osd38", str(FIGURES))


def main():
    prot = INPUTS.quantify()
    prot["relevant"] = (prot.p <= 0.05) & (prot.log2fc.abs() > 0.2)
    prot.to_csv(TABLES / "T11_osd38_proteome.tsv", sep="\t", index=False, float_format="%.5g")
    print(f"  T11: {prot.groupby('fraction').size().to_dict()} proteins")
    lay = layers(prot)
    selftest(lay["membrane_flight"])
    t = tests(lay)
    c, p522, rna522, rna38 = concordance(prot, lay)
    figure(t, prot, rna38, lay)


if __name__ == "__main__":
    main()
