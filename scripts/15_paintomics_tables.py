"""Tables and figures from the vendored PaintOmics records: jobs A, B, C and the first study's job 1.

  T04  every pathway: combined Fisher p, BH q within each job, rank, and per-omic p, for all four jobs
  T05  metabolite hubs: best adjusted p per compound and job, with the pre-registered mechanistic check
       (raw hub p identical in job 1 and job A for every compound and step the two share)
  T06  layer matching: the KEGG IDs PaintOmics used for the submitted names, against the intended ones
  T07  summary: pathway counts at FDR < 0.05, gained and lost, rank agreement, light vs dark
  fig02 job A vs job 1; fig03 job B (OSD-678 light) vs job C (OSD-678 dark)

Record formats follow the first study's scripts/06_paintomics_tables.py: a pathway's combined p is
combinedSignificancePvalues.Fisher, per-omic p is globalOmicPvalues, and a hub compound's result is its
best (lowest) pvalue_adjust over PaintOmics' steps. Job A's metabolite layer is a ground
carbon-starvation reference, not flight measurements; its gene and protein layers are job 1's.
"""

import gzip
import json

import numpy as np
import pandas as pd
from scipy import stats

import style
from paths import DATA, FIGURES, PAINTOMICS_RAW, PHOTORESP, TABLES

JOBS = {"job1": PHOTORESP / "job1_oK2h66740Q_full.json.gz",
        "A": PAINTOMICS_RAW / "jobA_6kBOS7AN5K_full.json.gz",
        "B": PAINTOMICS_RAW / "jobB_JfVZ1L37PY_full.json.gz",
        "C": PAINTOMICS_RAW / "jobC_wgT3IPi7U2_full.json.gz"}
OMICS = {"gene": "Gene expression", "protein": "Proteomics", "metabolite": "Metabolomics"}
ALPHA = 0.05


def load(path):
    with gzip.open(path, "rt") as fh:
        return json.load(fh)


def bh(p):
    p = np.asarray(p, float)
    o = np.argsort(p)
    q = p[o] * len(p) / np.arange(1, len(p) + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    out = np.empty_like(q)
    out[o] = np.minimum(q, 1)
    return out


def pathways(job, tag):
    rows = []
    for p in job["pathwaysInfo"]:
        g = p.get("globalOmicPvalues") or {}
        rows.append({"db": "KEGG" if p.get("source") == "KEGG" else "MapMan", "pathway": p["name"],
                     f"{tag}_p": p["combinedSignificancePvalues"]["Fisher"][0],
                     **{f"{tag}_p_{k}": g.get(v) for k, v in OMICS.items()}})
    d = pd.DataFrame(rows)
    d[f"{tag}_q"] = bh(d[f"{tag}_p"])
    d[f"{tag}_rank"] = d[f"{tag}_p"].rank(method="min").astype(int)
    return d.set_index(["db", "pathway"])


def hubs(job, tag):
    h = pd.DataFrame(job["hubAnalysisResult"].values())
    best = h.sort_values(["pvalue_adjust", "pvalue"]).groupby("name").first()
    return h.assign(job=tag), best[["pvalue_adjust", "pvalue", "DEN", "step"]].add_prefix(f"{tag}_")


def main():
    jobs = {k: load(v) for k, v in JOBS.items()}
    names = {cid: n for j in jobs.values() for cid, n in (j.get("mappingComp") or {}).items()}

    # T04 pathways
    t4 = pd.concat([pathways(j, k) for k, j in jobs.items()], axis=1).reset_index()
    t4 = t4.dropna(axis=1, how="all").sort_values("A_p")
    t4.to_csv(TABLES / "T04_paintomics_pathways.tsv", sep="\t", index=False, float_format="%.4g")

    # T05 hubs and the mechanistic check
    raw, best = zip(*(hubs(j, k) for k, j in jobs.items()))
    t5 = pd.concat(best, axis=1)
    t5.insert(0, "compound", [names.get(c, "") for c in t5.index])
    allraw = pd.concat(raw)
    shared = (allraw[allraw.job.isin(["job1", "A"])].pivot_table(index=["name", "step"], columns="job",
                                                                 values="pvalue").dropna())
    identical = bool(np.allclose(shared.job1, shared.A, rtol=0, atol=1e-12))
    t5["in_job1_and_A"] = t5.index.isin(shared.index.get_level_values(0))
    t5.index.name = "kegg_id"
    t5 = t5.sort_values("A_pvalue_adjust")
    t5.to_csv(TABLES / "T05_paintomics_hubs.tsv", sep="\t", float_format="%.4g")

    # T06 layer matching
    key = pd.read_csv(DATA / "paintomics_upload" / "metabolomics_key.tsv", sep="\t")
    rows = []
    for tag in ("A", "B", "C"):
        used = set(jobs[tag]["exprssionMetabolites"])
        intended = set(key.kegg_id)
        rows += [{"job": tag, "kegg_id": k, "status": "intended, used" if k in used else "intended, NOT used"}
                 for k in sorted(intended)]
        rows += [{"job": tag, "kegg_id": k, "status": "used, not intended"} for k in sorted(used - intended)]
    t6 = pd.DataFrame(rows)
    t6["compound"] = t6.kegg_id.map(names).fillna("")
    t6.to_csv(TABLES / "T06_layer_matching.tsv", sep="\t", index=False)

    # T07 summary
    sig = {k: set(t4.loc[t4[f"{k}_q"] < ALPHA, "pathway"]) for k in JOBS}
    hubsig = {k: set(t5.index[t5[f"{k}_pvalue_adjust"] < ALPHA]) for k in JOBS}
    rho_1a = stats.spearmanr(t4.job1_p, t4.A_p).correlation
    rho_bc = stats.spearmanr(t4.B_p, t4.C_p).correlation
    s = [("pathways tested", *(len(t4) for _ in JOBS)),
         ("pathways at BH q < 0.05", *(len(sig[k]) for k in JOBS)),
         ("hub compounds tested", *(int(t5[f"{k}_pvalue_adjust"].notna().sum()) for k in JOBS)),
         ("hub compounds at FDR < 0.05", *(len(hubsig[k]) for k in JOBS))]
    t7 = pd.DataFrame(s, columns=["quantity", *JOBS])
    extra = pd.DataFrame([
        ("A vs job1: Spearman rho of combined p", round(rho_1a, 4)),
        ("A vs job1: pathways gained at q < 0.05", "; ".join(sorted(sig["A"] - sig["job1"])) or "none"),
        ("A vs job1: pathways lost at q < 0.05", "; ".join(sorted(sig["job1"] - sig["A"])) or "none"),
        ("A vs job1: hubs gained", "; ".join(sorted(names.get(c, c) for c in hubsig["A"] - hubsig["job1"]))),
        ("A vs job1: hubs lost", "; ".join(sorted(names.get(c, c) for c in hubsig["job1"] - hubsig["A"]))),
        ("A vs job1: raw hub p identical for shared compounds and steps (pre-registered check)",
         f"{identical} ({len(shared)} compound-steps)"),
        ("B vs C: Spearman rho of combined p", round(rho_bc, 4)),
        ("B only (light) at q < 0.05", "; ".join(sorted(sig["B"] - sig["C"])) or "none"),
        ("C only (dark) at q < 0.05", "; ".join(sorted(sig["C"] - sig["B"])) or "none"),
        ("B and C at q < 0.05", "; ".join(sorted(sig["B"] & sig["C"])) or "none"),
        ("hubs B only (light)", "; ".join(sorted(names.get(c, c) for c in hubsig["B"] - hubsig["C"])) or "none"),
        ("hubs C only (dark)", "; ".join(sorted(names.get(c, c) for c in hubsig["C"] - hubsig["B"])) or "none"),
        ("hubs B and C", "; ".join(sorted(names.get(c, c) for c in hubsig["B"] & hubsig["C"])) or "none"),
        ("layer: intended IDs not used (A)", int((t6[(t6.job == 'A')].status == "intended, NOT used").sum())),
        ("layer: used IDs not intended (A)", int((t6[(t6.job == 'A')].status == "used, not intended").sum())),
    ], columns=["quantity", "value"])
    with open(TABLES / "T07_paintomics_summary.tsv", "w") as fh:
        fh.write("# Jobs: job1 = first study oK2h66740Q (OSD-522 genes + proteins + 21 predicted compounds); "
                 "A = 6kBOS7AN5K (same genes + proteins + reference layer); B = JfVZ1L37PY (OSD-678 light "
                 "genes + layer); C = wgT3IPi7U2 (OSD-678 dark genes + layer). The reference layer is a ground "
                 "carbon-starvation reference, not flight measurements.\n")
        t7.to_csv(fh, sep="\t", index=False)
        fh.write("\n")
        extra.to_csv(fh, sep="\t", index=False)
    print(t7.to_string(index=False))
    print(extra.to_string(index=False, max_colwidth=160))
    figures(t4, t5, names)


def figures(t4, t5, names):
    style.apply()
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D

    def scatter(ax, x, y, qx, qy, labels, xl, yl):
        lx, ly = -np.log10(x), -np.log10(y)
        both, onlyx, onlyy = (qx < ALPHA) & (qy < ALPHA), (qx < ALPHA) & (qy >= ALPHA), (qx >= ALPHA) & (qy < ALPHA)
        rest = ~(both | onlyx | onlyy)
        ax.scatter(lx[rest], ly[rest], s=9, color=style.GREY_LIGHT, edgecolor=style.GREY, lw=0.3)
        ax.scatter(lx[both], ly[both], s=16, color=style.INK)
        ax.scatter(lx[onlyx], ly[onlyx], s=16, color=style.BLUE)
        ax.scatter(lx[onlyy], ly[onlyy], s=16, color=style.VERMILION)
        m = max(lx.max(), ly.max()) * 1.05
        ax.plot([0, m], [0, m], color=style.GREY, lw=0.8, ls="--", zorder=0)
        ax.set_xlim(0, m)
        ax.set_ylim(0, m)
        for sel, key in ((onlyx, lx - ly), (onlyy, ly - lx)):
            idx = np.flatnonzero(sel)
            for i in idx[np.argsort(-key[idx])][:4]:
                ax.annotate(labels[i][:30], (lx[i], ly[i]), fontsize=5.6, xytext=(4, 2),
                            textcoords="offset points", color=style.INK)
        ax.set_xlabel(f"{xl}: -log10 combined p")
        ax.set_ylabel(f"{yl}: -log10 combined p")
        return [Line2D([], [], marker="o", ls="", color=style.INK, label="q < 0.05 in both"),
                Line2D([], [], marker="o", ls="", color=style.BLUE, label=f"{xl} only"),
                Line2D([], [], marker="o", ls="", color=style.VERMILION, label=f"{yl} only"),
                Line2D([], [], marker="o", ls="", color=style.GREY_LIGHT, mec=style.GREY, label="neither")]

    def hubpanel(ax, a, b, la, lb):
        h = t5[[f"{a}_pvalue_adjust", f"{b}_pvalue_adjust"]].copy()
        h = h[(h < ALPHA).any(axis=1)]
        h["name"] = [names.get(c, c) for c in h.index]
        h = h.sort_values(f"{b}_pvalue_adjust")
        y = np.arange(len(h))[::-1]
        for col, colr, mk in ((a, style.BLUE, "o"), (b, style.VERMILION, "D")):
            v = h[f"{col}_pvalue_adjust"]
            ax.scatter(-np.log10(v[v.notna()]), y[v.notna().to_numpy()], s=22, color=colr, marker=mk, zorder=3)
            ax.scatter(np.zeros(int(v.isna().sum())), y[v.isna().to_numpy()], s=22, facecolor="white",
                       edgecolor=colr, marker=mk, zorder=3)
        for yi, (pa, pb) in zip(y, zip(h[f"{a}_pvalue_adjust"].fillna(1), h[f"{b}_pvalue_adjust"].fillna(1))):
            ax.plot([-np.log10(pa), -np.log10(pb)], [yi, yi], color=style.GREY_LIGHT, lw=1.2, zorder=1)
        ax.axvline(-np.log10(ALPHA), color=style.GREY, ls="--", lw=0.8)
        ax.set_yticks(y)
        ax.set_yticklabels([n[:30] for n in h.name], fontsize=6.2)
        ax.set_xlabel("Hub -log10 FDR (best step)")
        ax.grid(axis="y", visible=False)
        return [Line2D([], [], marker="o", ls="", color=style.BLUE, label=la),
                Line2D([], [], marker="D", ls="", color=style.VERMILION, label=lb),
                Line2D([], [], marker="o", ls="", mfc="white", mec=style.GREY, label="not entered in that job's hub test")]

    for stem, (a, b, la, lb, title) in {
        "fig02_jobA_vs_job1": ("job1", "A", "job 1 (predicted layer)", "job A (reference layer)",
                               "OSD-522: model-predicted vs reference metabolite layer, same genes and proteins"),
        "fig03_light_vs_dark": ("B", "C", "OSD-678 light", "OSD-678 dark",
                                "OSD-678 leaves: light vs dark flight, same reference layer"),
    }.items():
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(7.2, 4.4), gridspec_kw={"width_ratios": [1.15, 1]})
        h1 = scatter(ax1, t4[f"{a}_p"].to_numpy(), t4[f"{b}_p"].to_numpy(), t4[f"{a}_q"].to_numpy(),
                     t4[f"{b}_q"].to_numpy(), t4.pathway.tolist(), la, lb)
        style.panel(ax1, "a", "Pathways (232; BH q within job)")
        ax1.legend(handles=h1, loc="upper left", bbox_to_anchor=(0, -0.13), ncol=2, fontsize=6.2)
        h2 = hubpanel(ax2, a, b, la, lb)
        style.panel(ax2, "b", "Metabolite hubs at FDR < 0.05 in either")
        ax2.legend(handles=h2, loc="upper left", bbox_to_anchor=(-0.1, -0.13), ncol=1, fontsize=6.2)
        fig.suptitle(title, x=0.01, ha="left", fontsize=9, fontweight="bold", y=0.995)
        fig.subplots_adjust(wspace=0.6, top=0.88)
        fig.text(0.01, -0.12, "Reference layer: ground carbon-starvation responses (T02), not flight measurements. "
                 "Gene and protein layers: NASA OSDR.", fontsize=6.4, color=style.GREY, ha="left")
        style.save(fig, stem, str(FIGURES))


if __name__ == "__main__":
    main()
