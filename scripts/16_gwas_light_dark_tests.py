"""Test the pre-registered light/dark predictions on the six flight contrasts (T08, T09, fig04).

Everything this script decides is read from prereg/ (settings.json, PREDICTIONS.tsv, gene_sets.tsv,
gene_order.tsv.gz), committed before this script existed. It refuses to run without them.

  OSD-522      Wald statistics recomputed from the raw RSEM counts with PyDESeq2, as the first study
               did (its log2FC is checked to match the vendored job 1 upload). Primary null: all 924
               ways to label 6 of the 12 samples as flight, gene statistic a Welch t on
               log2(median-of-ratios normalised count + 1).
  other five   GeneLab Wald statistics (data/osdr/). Primary null: genomic rotation, the set's
               positions shifted together along the genome, 10000 offsets.

Before testing, `selftest()` plants a known shift and an unshifted control set into real data and
checks that both nulls find the first and not the second; if either fails, nothing is written.
"""

import itertools
import json
import sys

import numpy as np
import pandas as pd
from scipy import stats

import genesets
import style
from paths import ATLAS, FIGURES, OSDR, PHOTORESP, PREREG, TABLES

SETTINGS = PREREG / "settings.json"
if not SETTINGS.exists():
    sys.exit("prereg/settings.json missing: run scripts/11_prereg.py and commit it first")
S = json.loads(SETTINGS.read_text())
LIT, DARK = S["contrasts"]["lit"], S["contrasts"]["dark"]
B, SEED, MIN_COUNT, FDR = S["rotation_draws"], S["seed"], S["min_total_count_osd522"], S["fdr"]


# ---------- gene statistics ----------

def osd522():
    """PyDESeq2 Wald stats and the normalised log counts for the sample permutation."""
    from pydeseq2.dds import DeseqDataSet
    from pydeseq2.ds import DeseqStats

    counts = pd.read_csv(ATLAS / "GLDS-522_rna_seq_RSEM_Unnormalized_Counts_GLbulkRNAseq.csv", index_col=0)
    rs = pd.read_csv(ATLAS / "GLDS-522_rna_seq_bulkRNASeq_v2_runsheet.csv").set_index("Sample Name")
    cond = rs.loc[counts.columns, "Factor Value[Spaceflight]"].map(
        {"Space Flight": "SpaceFlight", "Ground Control": "GroundControl"})
    assert cond.notna().all() and (cond == "SpaceFlight").sum() == 6
    counts = counts[counts.sum(axis=1) >= MIN_COUNT]
    dds = DeseqDataSet(counts=counts.T.astype(int), metadata=pd.DataFrame({"condition": cond}),
                       design_factors="condition", refit_cooks=True, quiet=True)
    dds.deseq2()
    st = DeseqStats(dds, contrast=["condition", "SpaceFlight", "GroundControl"], quiet=True)
    st.summary()
    res = st.results_df.rename(columns={"log2FoldChange": "log2fc"})[["log2fc", "stat"]].dropna()

    first = pd.read_csv(PHOTORESP / "gene_expression_values.tab", sep="\t", index_col=0).iloc[:, 0]
    common = res.index.intersection(first.index)
    gap = float((res.loc[common, "log2fc"] - first[common]).abs().max())
    assert len(common) == len(first) and gap < 1e-5, f"OSD-522 log2FC differs from job 1 by {gap}"

    sf = dds.obsm["size_factors"]
    logn = np.log2(counts.div(sf, axis=1) + 1)
    flight = (cond == "SpaceFlight").to_numpy()
    print(f"  OSD-522: {len(res)} genes, log2FC matches the first study's upload (max gap {gap:.1e})")
    return res, logn, flight


def contrast(name):
    d = pd.read_csv(OSDR / f"contrast_{name}.csv.gz", usecols=["TAIR", "log2fc", "stat"])
    return d.dropna().drop_duplicates("TAIR").set_index("TAIR")


# ---------- nulls ----------

def set_stats(values, mask):
    return values[mask].mean() - values.mean(), np.abs(values[mask]).mean() - np.abs(values).mean()


def rotation(values, mask, rng):
    """Null set statistics with the set's positions rotated together along the gene order."""
    n = len(values)
    idx = np.flatnonzero(mask)
    offsets = rng.integers(1, n, size=B)
    absv, mv, mabs = np.abs(values), values.mean(), np.abs(values).mean()
    signed, unsigned = np.empty(B), np.empty(B)
    for lo in range(0, B, 1000):
        ix = (idx[:, None] + offsets[None, lo:lo + 1000]) % n
        signed[lo:lo + 1000] = values[ix].mean(axis=0) - mv
        unsigned[lo:lo + 1000] = absv[ix].mean(axis=0) - mabs
    return signed, unsigned


def welch_t(x, flight):
    a, b = x[:, flight], x[:, ~flight]
    va, vb = a.var(axis=1, ddof=1) / a.shape[1], b.var(axis=1, ddof=1) / b.shape[1]
    with np.errstate(divide="ignore", invalid="ignore"):
        t = (a.mean(axis=1) - b.mean(axis=1)) / np.sqrt(va + vb)
    return np.nan_to_num(t)


def labellings(n=12, k=6):
    for combo in itertools.combinations(range(n), k):
        lab = np.zeros(n, bool)
        lab[list(combo)] = True
        yield lab


def permutation(logn, flight, masks):
    """Set statistics under all 924 labellings; also returns the observed ones."""
    x = logn.to_numpy()
    out = {k: ([], []) for k in masks}
    for lab in labellings():
        t = welch_t(x, lab)
        for k, m in masks.items():
            s, u = set_stats(t, m)
            out[k][0].append(s)
            out[k][1].append(u)
    obs = welch_t(x, flight)
    return {k: (np.array(s), np.array(u), *set_stats(obs, masks[k])) for k, (s, u) in out.items()}


def p_two(null, obs, include_obs):
    hits = np.sum(np.abs(null) >= abs(obs) - 1e-12)
    return hits / len(null) if include_obs else (1 + hits) / (len(null) + 1)


def p_one(null, obs, include_obs):
    hits = np.sum(null >= obs - 1e-12)
    return hits / len(null) if include_obs else (1 + hits) / (len(null) + 1)


def bh(p):
    p = np.asarray(p, float)
    o = np.argsort(p)
    q = p[o] * len(p) / np.arange(1, len(p) + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    out = np.empty_like(q)
    out[o] = np.minimum(q, 1)
    return out


# ---------- planted-signal check ----------

def selftest(order, logn, flight):
    rng = np.random.default_rng(1)
    d = contrast("OSD-678-light")
    u = [g for g in order if g in d.index]
    v = d.loc[u, "stat"].to_numpy().copy()
    planted, control = np.zeros(len(v), bool), np.zeros(len(v), bool)
    planted[rng.choice(len(v), 60, replace=False)] = True
    control[rng.choice(np.flatnonzero(~planted), 60, replace=False)] = True
    v[planted] += 2.0
    ps, pu = (p_two(rotation(v, planted, rng)[0], set_stats(v, planted)[0], False),
              p_two(rotation(v, control, rng)[0], set_stats(v, control)[0], False))
    x = logn.copy()
    genes = rng.choice(len(x), 60, replace=False)
    ctrl = rng.choice(np.setdiff1d(np.arange(len(x)), genes), 60, replace=False)
    x.iloc[genes, np.flatnonzero(flight)] += 1.0
    m_p, m_c = np.zeros(len(x), bool), np.zeros(len(x), bool)
    m_p[genes], m_c[ctrl] = True, True
    perm = permutation(x, flight, {"planted": m_p, "control": m_c})
    pp = p_two(perm["planted"][0], perm["planted"][2], True)
    pc = p_two(perm["control"][0], perm["control"][2], True)
    ok = ps < 0.01 and pu > 0.05 and pp < 0.01 and pc > 0.05
    print(f"  selftest: rotation planted p={ps:.4f} control p={pu:.3f}; "
          f"permutation planted p={pp:.4f} control p={pc:.3f} -> {'PASS' if ok else 'FAIL'}")
    if not ok:
        sys.exit("selftest failed: the tests cannot find a planted signal, or find one in a control set")


# ---------- run ----------

def run():
    sets = genesets.load()
    order = pd.read_csv(PREREG / "gene_order.tsv.gz", sep="\t").gene.tolist()
    pred = pd.read_csv(PREREG / "PREDICTIONS.tsv", sep="\t")
    res522, logn, flight = osd522()
    selftest(order, logn, flight)

    rows = []
    perm = permutation(logn, flight, {k: logn.index.isin(list(v)) for k, v in sets.items()})
    for c in LIT + DARK:
        d = res522 if c == "OSD-522" else contrast(c)
        u = [g for g in order if g in d.index]
        v = d.loc[u, "stat"].to_numpy()
        rng = np.random.default_rng(SEED)
        for name, genes in sets.items():
            m = np.isin(u, list(genes))
            s_obs, u_obs = set_stats(v, m)
            rs, ru = rotation(v, m, rng)
            inside = d.log2fc[d.index.isin(list(genes))]
            outside = d.log2fc[~d.index.isin(list(genes))]
            mw = stats.mannwhitneyu(inside, outside, alternative="two-sided").pvalue
            row = {"set": name, "contrast": c, "light": "lit" if c in LIT else "dark",
                   "n_set": len(genes), "n_used": int(m.sum()),
                   "signed": s_obs, "unsigned": u_obs,
                   "rot_p_signed": p_two(rs, s_obs, False), "rot_p_unsigned": p_one(ru, u_obs, False),
                   "rot_z_signed": (s_obs - rs.mean()) / rs.std(),
                   "rot_z_unsigned": (u_obs - ru.mean()) / ru.std(),
                   "mw_p": mw, "mw_shift": inside.median() - outside.median()}
            if c == "OSD-522":
                ps, pu, so, uo = perm[name]
                row.update({"perm_signed": so, "perm_unsigned": uo,
                            "perm_p_signed": p_two(ps, so, True), "perm_p_unsigned": p_one(pu, uo, True),
                            "perm_z_signed": (so - ps.mean()) / ps.std(),
                            "perm_z_unsigned": (uo - pu.mean()) / pu.std()})
            rows.append(row)
    t = pd.DataFrame(rows).merge(pred[["set", "contrast", "statistic", "predicted"]], on=["set", "contrast"])
    t["null"] = np.where(t.contrast == "OSD-522", "permutation", "rotation")
    pick = lambda r, what: r[f"{'perm' if r.null == 'permutation' else 'rot'}_{what}_{r.statistic}"]
    t["p"] = [pick(r, "p") for _, r in t.iterrows()]
    t["z"] = [pick(r, "z") for _, r in t.iterrows()]
    t["q"] = bh(t.p)
    sign = np.where(t.statistic == "signed", t.signed, t.unsigned)
    t["observed"] = np.where(t.q >= FDR, "none",
                             np.where(t.statistic == "unsigned", "enriched",
                                      np.where(sign > 0, "up", "down")))
    t["observed"] = np.where((t.statistic == "unsigned") & (t.q < FDR) & (t.unsigned < 0), "depleted",
                             t.observed)
    t["verdict"] = np.where(t.observed == t.predicted, "MATCH", "MISMATCH")
    lead = ["set", "contrast", "light", "statistic", "null", "predicted", "observed", "verdict", "p", "q", "z",
            "n_set", "n_used"]
    t = t[lead + [c for c in t.columns if c not in lead]]
    t.to_csv(TABLES / "T08_geneset_tests.tsv", sep="\t", index=False, float_format="%.4g")

    summ = (t.groupby(["set", "light"]).verdict.apply(lambda v: f"{(v == 'MATCH').sum()}/{len(v)}")
            .unstack().reindex(list(sets))[["lit", "dark"]])
    summ["predicted_lit"] = [pred[(pred.set == s) & (pred.light == "lit")].predicted.iloc[0] for s in summ.index]
    summ["predicted_dark"] = [pred[(pred.set == s) & (pred.light == "dark")].predicted.iloc[0] for s in summ.index]
    summ.rename(columns={"lit": "match_lit", "dark": "match_dark"}).to_csv(
        TABLES / "T09_geneset_summary.tsv", sep="\t")
    print(t[["set", "contrast", "predicted", "observed", "verdict", "p", "q", "z"]].to_string(index=False))
    print(f"  T08: {len(t)} primary tests, {int((t.verdict == 'MATCH').sum())} MATCH")
    return t


def figure(t):
    style.apply()
    import matplotlib.pyplot as plt
    from matplotlib.colors import TwoSlopeNorm

    sets = list(dict.fromkeys(t.set))
    cols = LIT + DARK
    z = t.pivot(index="set", columns="contrast", values="z").reindex(index=sets, columns=cols)
    fig, ax = plt.subplots(figsize=(7.0, 3.5))
    lim = max(4.0, float(np.nanmax(np.abs(z.to_numpy()))))
    im = ax.imshow(z.to_numpy(), cmap="RdBu_r", norm=TwoSlopeNorm(0, -lim, lim), aspect="auto")
    for i, s in enumerate(sets):
        for j, c in enumerate(cols):
            r = t[(t.set == s) & (t.contrast == c)].iloc[0]
            ink = "white" if abs(r.z) > lim * 0.55 else style.INK
            ax.text(j, i - 0.13, f"{r.observed}", ha="center", va="center", fontsize=6.8,
                    fontweight="bold" if r.q < FDR else "normal", color=ink)
            ax.text(j, i + 0.23, f"{'=' if r.verdict == 'MATCH' else '≠'} pred. {r.predicted}", ha="center", va="center",
                    fontsize=5.6, color=ink, style="normal" if r.verdict == "MATCH" else "italic")
            if r.verdict == "MISMATCH" and r.q < FDR:
                ax.add_patch(plt.Rectangle((j - 0.47, i - 0.46), 0.94, 0.92, fill=False, lw=1.4,
                                           ec=style.INK, zorder=3))
    ax.set_xticks(range(len(cols)))
    ax.set_xticklabels(cols, fontsize=7)
    ax.set_yticks(range(len(sets)))
    ax.set_yticklabels([f"{s} ({t[t.set == s].statistic.iloc[0]})" for s in sets], fontsize=7)
    ax.axvline(len(LIT) - 0.5, color=style.INK, lw=1.2)
    ax.text((len(LIT) - 1) / 2, -0.62, "lit flight (model: starvation signal)", ha="center", va="bottom",
            fontsize=7.3, color=style.INK, fontweight="bold")
    ax.text(len(LIT) + (len(DARK) - 1) / 2, -0.62, "dark flight (model: no signal)", ha="center", va="bottom",
            fontsize=7.3, color=style.INK, fontweight="bold")
    ax.grid(False)
    cb = fig.colorbar(im, ax=ax, fraction=0.03, pad=0.02)
    cb.set_label("Set statistic vs its null (z)", fontsize=7)
    ax.set_title("Pre-registered gene-set tests (bold: BH q < 0.05; box: significant and contrary "
                 "to the prediction)", loc="left", fontsize=8.2, pad=20)
    style.save(fig, "fig04_geneset_tests", str(FIGURES))


if __name__ == "__main__":
    figure(run())
