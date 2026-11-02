"""Diagnostics of the pre-registered gene-set tests (T20-T22, figS2, figS3).

Revision 1, post hoc: the settings are in prereg/revision1/settings.json (R2), committed before this
script ran. Nothing here changes T08.

- T20: calibration of the genomic-rotation null. On each flight contrast, random gene sets and
  contiguous genomic blocks of five sizes are tested with the rotation null, and against a naive
  random-gene null. A calibrated test rejects about 5 % of them at alpha 0.05.
- T21: the six pre-registered sets under exact sample permutation (3 v 3, all 20 labellings) for
  OSD-38 and both OSD-678 arms; and, for OSD-38, the same sets re-tested with the rotation null on
  PyDESeq2 statistics recomputed from GeneLab's raw counts, to compare the two pipelines.
- T22: the five DIN genes one by one in every contrast, the DIN set test with each gene left out, and
  each DIN gene's log2 normalised count per sample in the three 3 v 3 contrasts.

The rotation statistic for every offset is computed at once by FFT (circular cross-correlation of
the set indicator with the gene statistics); it is the statistic of scripts/16 `rotation`, checked
against it below.
"""

import importlib.util
import json
import warnings

import numpy as np
import pandas as pd
from scipy import stats

import fetch_record
import genesets
import style
from paths import FIGURES, OSDR, PREREG, ROOT, TABLES

warnings.filterwarnings("ignore", category=RuntimeWarning, module="pydeseq2")
R2 = json.loads((PREREG / "revision1" / "settings.json").read_text())["R2_geneset_diagnostics"]
CAL = R2["rotation_calibration"]
EXACT = R2["exact_sample_permutation"]


def load(name):
    spec = importlib.util.spec_from_file_location(name.replace(".py", ""), ROOT / "scripts" / name)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


T16 = load("16_gwas_light_dark_tests.py")
PRIMARY = T16.S["primary_statistic"]
ORDER = pd.read_csv(PREREG / "gene_order.tsv.gz", sep="\t").gene.tolist()
CONTRASTS = CAL["contrasts"]
assert sorted(CONTRASTS) == sorted(T16.LIT + T16.DARK)


def ordered(d):
    u = [g for g in ORDER if g in d.index]
    return u, d.loc[u, "stat"].to_numpy(float)


def rotated(v, mask):
    """Set mean minus the mean of all genes, at every circular offset; index 0 is the observed set."""
    n = len(v)
    s = np.fft.irfft(np.conj(np.fft.rfft(mask.astype(float))) * np.fft.rfft(v), n)
    return s / mask.sum() - v.mean()


def check_fft(v):
    rng = np.random.default_rng(0)
    m = np.zeros(len(v), bool)
    m[rng.choice(len(v), 40, replace=False)] = True
    ref, _ = T16.rotation(v, m, np.random.default_rng(5))
    off = np.random.default_rng(5).integers(1, len(v), size=T16.B)
    gap = np.abs(rotated(v, m)[off] - ref).max()
    assert gap < 1e-9, f"FFT rotation differs from scripts/16 rotation by {gap}"


# ---------- T20: rotation calibration ----------

def calibration(vals):
    rng = np.random.default_rng(CAL["seed"])
    B, per = CAL["offsets"], CAL["sets_per_size_and_type"]
    rows = []
    for c in CONTRASTS:
        v = vals[c]
        n = len(v)
        for k in CAL["set_sizes"]:
            naive = np.array([v[rng.choice(n, k, replace=False)].mean() for _ in range(B)]) - v.mean()
            for kind in CAL["set_types"]:
                for _ in range(per):
                    m = np.zeros(n, bool)
                    if kind == "random genes":
                        m[rng.choice(n, k, replace=False)] = True
                    else:
                        m[(rng.integers(n) + np.arange(k)) % n] = True
                    r = rotated(v, m)
                    off = rng.integers(1, n, size=B)
                    rows.append({"contrast": c, "set_type": kind, "set_size": k,
                                 "p_rotation": T16.p_two(r[off], r[0], False),
                                 "p_naive": T16.p_two(naive, r[0], False)})
    p = pd.DataFrame(rows)
    out = []
    for (c, kind, k), g in p.groupby(["contrast", "set_type", "set_size"], sort=False):
        for null in ("rotation", "naive"):
            x = int((g[f"p_{null}"] < 0.05).sum())
            ci = stats.binomtest(x, len(g)).proportion_ci(0.95, method="wilson")
            out.append({"contrast": c, "set_type": kind, "set_size": k, "null": null, "n_sets": len(g),
                        "rejected_at_0.05": x, "type1_error": x / len(g), "ci_low": ci.low, "ci_high": ci.high})
    t = pd.DataFrame(out)
    t.to_csv(TABLES / "T20_rotation_calibration.tsv", sep="\t", index=False, float_format="%.4g")
    s = t.groupby(["set_type", "null"]).type1_error.agg(["min", "median", "max"])
    print("  T20 type I error at 0.05 (over contrasts and sizes):\n" + s.round(3).to_string())
    return t, p


# ---------- T21: exact permutation and pipeline check ----------

def welch(x, lab):
    a, b = x[:, lab], x[:, ~lab]
    num = a.mean(axis=1) - b.mean(axis=1)
    den = np.sqrt(a.var(axis=1, ddof=1) / a.shape[1] + b.var(axis=1, ddof=1) / b.shape[1])
    return np.divide(num, den, out=np.zeros_like(num), where=den > 0)


def exact(c, sets):
    d = pd.read_csv(OSDR / f"contrast_{c}.csv.gz").dropna(subset=["TAIR"]).drop_duplicates("TAIR").set_index("TAIR")
    fl = [k for k in d.columns if k.startswith("flight:")]
    gr = [k for k in d.columns if k.startswith("ground:")]
    assert len(fl) == len(gr) == 3, c
    x = np.log2(d[fl + gr].to_numpy(float) + 1)
    keep = np.isfinite(x).all(axis=1)
    x, genes = x[keep], d.index[keep]
    obs_lab = np.array([True] * 3 + [False] * 3)
    labs = list(T16.labellings(6, 3))
    assert len(labs) == 20 and any((lab == obs_lab).all() for lab in labs)
    masks = {k: genes.isin(list(v)) for k, v in sets.items()}
    null = {k: ([], []) for k in sets}
    for lab in labs:
        t = welch(x, lab)
        for k, m in masks.items():
            s, u = T16.set_stats(t, m)
            null[k][0].append(s)
            null[k][1].append(u)
    t_obs = welch(x, obs_lab)
    out = {}
    for k, m in masks.items():
        so, uo = T16.set_stats(t_obs, m)
        ns, nu = np.array(null[k][0]), np.array(null[k][1])
        out[k] = {"n_genes": len(genes), "n_used_perm": int(m.sum()), "perm_signed": so, "perm_unsigned": uo,
                  "perm_p_signed": T16.p_two(ns, so, True), "perm_p_unsigned": T16.p_one(nu, uo, True)}
    return out


def osd38_pydeseq2():
    from pydeseq2.dds import DeseqDataSet
    from pydeseq2.ds import DeseqStats

    counts = pd.read_csv(fetch_record.osdr("OSD-38", "GLDS-38_rna_seq_RSEM_Unnormalized_Counts_GLbulkRNAseq.csv"),
                         index_col=0)
    fetch_record.osdr("OSD-38", "GLDS-38_rna_seq_bulkRNASeq_v2_runsheet.csv")
    cols = pd.read_csv(OSDR / "contrast_OSD-38.csv.gz", nrows=0).columns
    fl = [k.split(":", 1)[1] for k in cols if k.startswith("flight:")]
    gr = [k.split(":", 1)[1] for k in cols if k.startswith("ground:")]
    counts = counts[fl + gr]
    counts = counts[counts.sum(axis=1) >= T16.MIN_COUNT]
    cond = pd.Series(["SpaceFlight"] * len(fl) + ["GroundControl"] * len(gr), index=fl + gr)
    dds = DeseqDataSet(counts=counts.T.astype(int), metadata=pd.DataFrame({"condition": cond}),
                       design_factors="condition", refit_cooks=True, quiet=True)
    dds.deseq2()
    st = DeseqStats(dds, contrast=["condition", "SpaceFlight", "GroundControl"], quiet=True)
    st.summary()
    return st.results_df.rename(columns={"log2FoldChange": "log2fc"})[["log2fc", "stat"]].dropna()


def t21(sets):
    t08 = pd.read_csv(TABLES / "T08_geneset_tests.tsv", sep="\t").set_index(["set", "contrast"])
    rows = []
    for c in EXACT["contrasts"]:
        for k, r in exact(c, sets).items():
            prim = PRIMARY[k]
            rows.append({"set": k, "contrast": c, "statistic": prim, **r,
                         "perm_p": r[f"perm_p_{prim}"], "T08_rotation_p": t08.loc[(k, c), "p"],
                         "T08_q": t08.loc[(k, c), "q"], "T08_observed": t08.loc[(k, c), "observed"]})
    t = pd.DataFrame(rows)

    res = osd38_pydeseq2()
    gl = T16.contrast("OSD-38")
    common = res.index.intersection(gl.index)
    rho_stat = stats.spearmanr(res.loc[common, "stat"], gl.loc[common, "stat"]).statistic
    rho_fc = stats.spearmanr(res.loc[common, "log2fc"], gl.loc[common, "log2fc"]).statistic
    u, v = ordered(res)
    for k, genes in sets.items():
        m = np.isin(u, list(genes))
        so, uo = T16.set_stats(v, m)
        rs, ru = T16.rotation(v, m, np.random.default_rng(T16.SEED))
        p = T16.p_two(rs, so, False) if PRIMARY[k] == "signed" else T16.p_one(ru, uo, False)
        z = (so - rs.mean()) / rs.std() if PRIMARY[k] == "signed" else (uo - ru.mean()) / ru.std()
        i = (t.set == k) & (t.contrast == "OSD-38")
        t.loc[i, ["pydeseq2_rotation_p", "pydeseq2_rotation_z", "T08_z"]] = [p, z, t08.loc[(k, "OSD-38"), "z"]]
    t.loc[t.contrast == "OSD-38", ["pipeline_spearman_stat", "pipeline_spearman_log2fc", "pipeline_genes"]] = \
        [rho_stat, rho_fc, len(common)]
    t["pipeline_genes"] = t.pipeline_genes.astype("Int64")   # a count: keep it out of %.4g
    t.to_csv(TABLES / "T21_exact_permutation.tsv", sep="\t", index=False, float_format="%.4g")
    print(t[["set", "contrast", "statistic", "perm_p", "T08_rotation_p", "T08_observed"]].to_string(index=False))
    print(f"  OSD-38 PyDESeq2 (6 samples) vs GeneLab: Spearman stat {rho_stat:.3f}, log2FC {rho_fc:.3f} over {len(common)} genes")
    print(t[t.contrast == "OSD-38"][["set", "T08_rotation_p", "pydeseq2_rotation_p", "T08_z", "pydeseq2_rotation_z"]].to_string(index=False))
    return t


# ---------- T22: DIN genes ----------

def t22(sets, res522, logn, flight):
    genes = sorted(sets["din_first_study"])
    t08 = pd.read_csv(TABLES / "T08_geneset_tests.tsv", sep="\t").set_index(["set", "contrast"])
    variants = {"none": genes, **{g: [x for x in genes if x != g] for g in genes}}
    perm = T16.permutation(logn, flight, {k: logn.index.isin(v) for k, v in variants.items()})
    rows = []
    for c in CONTRASTS:
        d = res522 if c == "OSD-522" else T16.contrast(c)
        for g in genes:
            rows.append({"kind": "gene", "contrast": c, "gene": g,
                         "stat": d.stat.get(g, np.nan), "log2fc": d.log2fc.get(g, np.nan)})
        u, v = ordered(d)
        for left, members in variants.items():
            m = np.isin(u, members)
            so, _ = T16.set_stats(v, m)
            rs, _ = T16.rotation(v, m, np.random.default_rng(T16.SEED))
            row = {"kind": "set" if left == "none" else "leave one out", "contrast": c,
                   "gene": "" if left == "none" else left, "n_used": int(m.sum()), "signed": so,
                   "rotation_p": T16.p_two(rs, so, False), "rotation_z": (so - rs.mean()) / rs.std()}
            if c == "OSD-522":
                ps, _, po, _ = perm[left]
                row |= {"perm_p": T16.p_two(ps, po, True), "perm_z": (po - ps.mean()) / ps.std()}
            rows.append(row)
    for c in EXACT["contrasts"]:
        d = pd.read_csv(OSDR / f"contrast_{c}.csv.gz").drop_duplicates("TAIR").set_index("TAIR")
        for col in [k for k in d.columns if k.startswith(("flight:", "ground:"))]:
            for g in genes:
                rows.append({"kind": "sample", "contrast": c, "gene": g, "sample": col,
                             "log2_norm_count": np.log2(d.loc[g, col] + 1)})
    t = pd.DataFrame(rows)
    full = t[t.kind == "set"].set_index("contrast")
    for c in CONTRASTS:
        want = t08.loc[("din_first_study", c), "p"]
        got = full.loc[c, "perm_p" if c == "OSD-522" else "rotation_p"]
        assert f"{got:.4g}" == f"{want:.4g}", f"DIN {c}: {got} vs T08 {want}"   # T08 is written at %.4g
    t.to_csv(TABLES / "T22_din_genes.tsv", sep="\t", index=False, float_format="%.4g")
    g = t[t.kind == "gene"].pivot(index="gene", columns="contrast", values="stat")[CONTRASTS]
    print("  T22 DIN gene Wald statistics:\n" + g.round(2).to_string())
    lo = t[t.kind.isin(["set", "leave one out"])].assign(p=lambda x: x.perm_p.where(x.contrast == "OSD-522", x.rotation_p))
    print(lo.pivot(index="gene", columns="contrast", values="p")[CONTRASTS].to_string(float_format=lambda x: f"{x:.4f}"))
    return t


# ---------- figures ----------

def fig_s2(t20, p):
    style.apply()
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D

    fig, (a, b, c) = plt.subplots(1, 3, figsize=(7.2, 2.7), gridspec_kw={"wspace": 0.42})
    rnd, blk = CAL["set_types"]
    cond = [(rnd, "rotation", style.BLUE, "random sets, rotation"),
            (blk, "rotation", style.GREEN, "genomic blocks, rotation"),
            (blk, "naive", style.VERMILION, "genomic blocks, naive null")]
    sizes = CAL["set_sizes"]
    pos = {k: i for i, k in enumerate(sizes)}
    for kind, null, colr, lab in cond:
        g = t20[(t20.set_type == kind) & (t20.null == null)]
        for cn, gc in g.groupby("contrast"):
            a.plot(gc.set_size.map(pos), gc.type1_error, color=colr, lw=0.6, alpha=0.4)
        med = g.groupby("set_size").type1_error.median()
        a.plot(med.index.map(pos), med.values, color=colr, lw=1.6, marker="o", ms=3, label=lab)
    a.axhline(0.05, color=style.INK, lw=0.8, ls=":")
    a.set_xticks(range(len(sizes)))
    a.set_xticklabels([str(s) for s in sizes])
    a.set_xlabel("Set size (genes)")
    a.set_ylabel("Rejected at α = 0.05")
    style.panel(a, "a", "Type I error")

    for ax, sel, title, letter in ((b, cond[:2], "Rotation null", "b"), (c, cond[2:], "Blocks, naive null by size", "c")):
        if letter == "b":
            for kind, null, colr, lab in sel:
                x = np.sort(p[p.set_type == kind][f"p_{null}"].to_numpy())
                e = (np.arange(1, len(x) + 1) - 0.5) / len(x)
                ax.plot(-np.log10(e), -np.log10(x), color=colr, lw=1.2)
        else:
            shades = plt.cm.Oranges(np.linspace(0.35, 0.95, len(sizes)))
            for k, colr in zip(sizes, shades):
                x = np.sort(p[(p.set_type == blk) & (p.set_size == k)].p_naive.to_numpy())
                e = (np.arange(1, len(x) + 1) - 0.5) / len(x)
                ax.plot(-np.log10(e), -np.log10(x), color=colr, lw=1.2, label=f"{k}")
            ax.legend(title="block size", fontsize=5.8, title_fontsize=6, loc="upper left")
        lim = 3.3
        ax.plot([0, lim], [0, lim], color=style.INK, lw=0.7, ls=":")
        ax.set_xlim(0, lim)
        ax.set_ylim(0, lim)
        ax.set_xlabel("Expected −log10 p")
        ax.set_ylabel("Observed −log10 p")
        style.panel(ax, letter, title)
    fig.legend(handles=[Line2D([], [], color=colr, lw=1.6, label=lab) for _, _, colr, lab in cond],
               loc="lower center", ncol=3, bbox_to_anchor=(0.5, -0.1), fontsize=6.4)
    fig.text(0.005, -0.15, f"{CAL['sets_per_size_and_type']} sets per size and type in each of the six flight contrasts; "
             f"{CAL['offsets']} rotation offsets or naive draws per test. Thin lines in a: single contrasts.",
             fontsize=6.0, color=style.GREY)
    style.save(fig, "figS2_rotation_calibration", FIGURES)


def fig_s3(t22):
    style.apply()
    import matplotlib.pyplot as plt
    from matplotlib.colors import TwoSlopeNorm

    g = t22[t22.kind == "gene"]
    genes = sorted(g.gene.unique())
    z = g.pivot(index="gene", columns="contrast", values="stat").reindex(index=genes, columns=CONTRASTS)
    fig, (a, b) = plt.subplots(1, 2, figsize=(7.2, 2.9), gridspec_kw={"width_ratios": [1.1, 1], "wspace": 0.5})
    lim = float(np.nanmax(np.abs(z.to_numpy())))
    im = a.imshow(z.to_numpy(), cmap="RdBu_r", norm=TwoSlopeNorm(0, -lim, lim), aspect="auto")
    for i in range(len(genes)):
        for j in range(len(CONTRASTS)):
            v = z.iloc[i, j]
            a.text(j, i, "n/a" if pd.isna(v) else f"{v:.1f}", ha="center", va="center", fontsize=6,
                   color="white" if abs(v) > lim * 0.6 else style.INK)
    a.set_xticks(range(len(CONTRASTS)))
    a.set_xticklabels(CONTRASTS, rotation=35, ha="right", fontsize=6.4)
    a.set_yticks(range(len(genes)))
    a.set_yticklabels(genes, fontsize=6.6)
    a.axvline(2.5, color=style.INK, lw=1)
    a.grid(False)
    fig.colorbar(im, ax=a, fraction=0.04, pad=0.02).set_label("Wald statistic", fontsize=6.5)
    style.panel(a, "a", "DIN genes, flight vs ground")

    s = t22[t22.kind.isin(["set", "leave one out"])].copy()
    s["z"] = s.perm_z.where(s.contrast == "OSD-522", s.rotation_z)
    for j, c in enumerate(CONTRASTS):
        lo = s[(s.contrast == c) & (s.kind == "leave one out")]
        full = s[(s.contrast == c) & (s.kind == "set")].z.iloc[0]
        b.scatter(np.full(len(lo), j) + np.linspace(-0.18, 0.18, len(lo)), lo.z, s=10, color=style.GREY, zorder=2)
        b.plot([j - 0.3, j + 0.3], [full, full], color=style.INK, lw=1.6, zorder=3)
    b.axhline(0, color=style.GREY, lw=0.7)
    b.axvline(2.5, color=style.INK, lw=1)
    b.set_xticks(range(len(CONTRASTS)))
    b.set_xticklabels(CONTRASTS, rotation=35, ha="right", fontsize=6.4)
    b.set_ylabel("DIN set statistic vs null (z)")
    style.panel(b, "b", "Leave one gene out")
    fig.text(0.005, -0.1, "Left of the vertical line: lit flights; right: dark flights. b: bar = all five genes; dots = "
             "one gene left out. OSD-522 uses the 924-labelling sample permutation, the others the rotation null.",
             fontsize=6.0, color=style.GREY)
    style.save(fig, "figS3_din_genes", FIGURES)


if __name__ == "__main__":
    sets = genesets.load()
    res522, logn, flight = T16.osd522()
    vals = {c: ordered(res522 if c == "OSD-522" else T16.contrast(c))[1] for c in CONTRASTS}
    check_fft(vals["OSD-38"])
    t20, p = calibration(vals)
    t21(sets)
    t22 = t22(sets, res522, logn, flight)
    fig_s2(t20, p)
    fig_s3(t22)
