"""Write manuscript/latex/supp_tables.tex: Supplementary Tables S2-S14, typeset from results/tables.

Like numbers.tex, the file is generated and never edited by hand, so the SI cannot drift from the
tables it shows. Supplementary Table S1 (where each source fact is stated) is written in
supplementary.tex itself, because it cites papers rather than tables.
"""

import numpy as np
import pandas as pd

from paths import DATA, ROOT, TABLES

OUT = ROOT / "manuscript" / "latex" / "supp_tables.tex"
SETS = {"din_first_study": "DIN", "bcaa_degradation": "BCAA degradation", "wu_stress_only": "Wu stress-only",
        "wu_dark_metabolite": "Wu darkness metabolites", "wu_control_only": "Wu control-only", "wu_shared": "Wu shared",
        "photosynthesis_apparatus": "Photosystem", "fermentation": "Fermentation",
        "hypoxia_responsive": "Hypoxia-responsive", "ath00630": "Glyoxylate/dicarboxylate", "ath00710": "Carbon fixation",
        "ath00500": "Starch and sucrose", "ath00010": "Glycolysis", "hypoxia_up": "Hypoxia-up signature",
        "hypoxia_down": "Hypoxia-down signature"}


def t(name, **kw):
    return pd.read_csv(TABLES / name, sep="\t", comment="#", **kw)


def esc(x):
    s = str(x)
    for a, b in (("\\", "\\textbackslash{}"), ("&", "\\&"), ("%", "\\%"), ("_", "\\_"), ("#", "\\#"),
                 ("$", "\\$"), ("{", "\\{"), ("}", "\\}"), ("~", "\\textasciitilde{}"), ("^", "\\^{}"),
                 ("<", "\\textless{}"), (">", "\\textgreater{}")):
        s = s.replace(a, b)
    return s


def num(v, d=2):
    if v is None or (isinstance(v, float) and np.isnan(v)):
        return "--"
    if isinstance(v, (bool, np.bool_)):
        return "yes" if v else "no"
    if isinstance(v, (int, np.integer)):
        return f"{int(v):,}".replace(",", "{,}")
    v = float(v)
    if v != 0 and abs(v) < 1e-3:
        m, e = f"{v:.1e}".split("e")
        return f"${m}\\times10^{{{int(e)}}}$"
    s = f"{v:.{d}f}"
    return f"${s}$" if s.startswith("-") else s


def table(label, caption, cols, rows, spec, size="\\scriptsize"):
    """A longtable, so the long ones break across pages."""
    head = " & ".join(cols) + " \\\\"
    out = [f"\\begingroup{size}\\setlength{{\\tabcolsep}}{{3pt}}",
           f"\\begin{{longtable}}{{{spec}}}",
           f"\\caption{{{caption}}}\\label{{{label}}}\\\\",
           "\\toprule", head, "\\midrule", "\\endfirsthead",
           f"\\multicolumn{{{len(cols)}}}{{l}}{{\\textit{{(continued)}}}}\\\\", "\\toprule", head, "\\midrule", "\\endhead",
           "\\bottomrule", "\\endfoot"]
    out += [" & ".join(r) + " \\\\" for r in rows]
    out += ["\\end{longtable}", "\\endgroup", ""]
    return "\n".join(out)


def call(status, d):
    if status in ("absent", "none", "conflicting"):
        return {"absent": "--", "none": "none", "conflicting": "confl."}[status]
    arrow = "$\\uparrow$" if d > 0 else "$\\downarrow$"
    return f"{status} {arrow}"


def s2():
    d = t("T14_reference_strata.tsv")
    rows = []
    for r in d.itertuples():
        src = []
        for so in ("zhu", "caldana", "wu_stress"):
            e = getattr(r, f"{so}_effect")
            if pd.isna(e):
                src.append("--")
                continue
            sup = getattr(r, f"{so}_support")
            q = "*" if getattr(r, f"{so}_qualifies") else ""
            src.append(num(e) + q + " / " + (num(sup, 3) if so == "caldana" else f"{100 * sup:.0f}\\,\\%") +
                       f" ({int(getattr(r, f'{so}_n'))})")
        rows.append([esc(r.name[:30]), r.kegg_id, *src, call(r.status_all, r.direction_all),
                     call(r.status_dark_only, r.direction_dark_only), "yes" if r.dark_heat_only else ""])
    return table("tab:s2", "\\textbf{Reference consensus by stratum (T14).} Per source: effect (log$_2$; * the source "
                 "qualifies) and support (Zhu, Wu: share of accessions moving with the median, accessions in brackets; "
                 "Caldana: Welch $p$, values in brackets). Consensus with all three sources (= T02) and with darkness "
                 "only (Zhu, Caldana); ``dark+heat only'': directional only because of Wu's 32\\,\\textdegree C stress.",
                 ["Compound", "KEGG", "Zhu", "Caldana", "Wu (32\\,\\textdegree C)", "All three",
                  "Dark only", "Heat only"], rows, "@{}llrrrlll@{}", "\\tiny")


def s3():
    d = t("T15_reference_sensitivity.tsv")
    rows = []
    arrow = lambda v: v.replace(" up", " $\\uparrow$").replace(" down", " $\\downarrow$").replace("absent", "--")
    for _, r in d.iterrows():
        rows.append([("all three" if r.stratum == "all" else "dark only") + (" (pre-reg.)" if r.preregistered else ""),
                     f"{100 * r.agreement:.0f}", str(r.caldana_min_minutes), str(r.robust), str(r.single),
                     *[arrow(r[k]) for k in ("Leu", "Ile", "Val", "Malate", "2-OG", "Glu", "Gln")],
                     f"{r.model_agrees}/{r.model_disagrees}", str(r.robust_not_predicted)])
    return table("tab:s3", "\\textbf{Threshold sensitivity of the reference (T15).} Accession agreement (\\%) and the "
                 "start of the Caldana window (min). Key calls: robust or single, with direction. Model: agrees/disagrees "
                 "among the scorable predictions; last column, robust responses the model did not predict.",
                 ["Stratum", "Agr.", "Cald.", "Rob.", "Sing.", "Leu", "Ile", "Val", "Malate", "2-OG", "Glu", "Gln",
                  "Model", "Not pred."], rows, "@{}lrrrrlllllllrr@{}", "\\tiny")


def s4():
    d = t("T16_kegg_mapping_report.tsv").fillna("")
    rows = [[esc(r.section), esc(r.item), str(r.n), str(r.of).replace(".0", "")] for r in d.itertuples()]
    return table("tab:s4", "\\textbf{How measured names reached KEGG (T16).} Candidates from KEGG \\texttt{find} on each "
                 "name; the rule takes an exact synonym match; overrides are reviewed by hand "
                 "(\\texttt{data/reference/kegg\\_overrides.tsv}). Wu classes are Wu et al.'s.",
                 ["Section", "Item", "$n$", "of"], rows, "@{}llrr@{}")


def s5():
    d = t("T08_geneset_tests.tsv")
    rows = [[SETS[r.set], r.contrast, r.statistic, r.null, r.predicted, r.observed, r.verdict, num(r.signed if r.statistic == "signed" else r.unsigned, 3),
             num(r.z), num(r.p, 3), num(r.q, 3), f"{r.n_used}/{r.n_set}"] for r in d.itertuples()]
    return table("tab:s5", "\\textbf{All 36 pre-registered gene-set tests (T08).} Set statistic: mean gene statistic "
                 "(signed) or mean absolute statistic (unsigned) in the set minus the genome; $z$ and $p$ against the "
                 "null; BH $q$ over the 36. $n$: measured/listed members. Members: \\texttt{prereg/gene\\_sets.tsv} "
                 "(KEGG retrieved 30 October 2026).",
                 ["Set", "Contrast", "Stat.", "Null", "Pred.", "Obs.", "Verdict", "Set stat.", "$z$", "$p$", "$q$",
                  "$n$"], rows, "@{}lllllllrrrrr@{}", "\\tiny")


def s6():
    d = t("T20_rotation_calibration.tsv")
    g = d.groupby(["set_type", "null", "set_size"]).type1_error.agg(["median", "min", "max"]).reset_index()
    rows = [[r.set_type, r.null, str(r.set_size), f"{r['median']:.3f}", f"{r['min']:.3f}--{r['max']:.3f}"]
            for _, r in g.iterrows()]
    return table("tab:s6", "\\textbf{Calibration of the rotation null (T20).} Share of 200 sets rejected at "
                 "$\\alpha=0.05$, median and range over the six flight contrasts. A calibrated null rejects about 0.05.",
                 ["Set type", "Null", "Size", "Median", "Range"], rows, "@{}llrrr@{}")


def s7():
    d = t("T21_exact_permutation.tsv")
    rows = [[SETS[r.set], r.contrast, r.statistic, num(r.perm_p, 1), num(r.T08_rotation_p, 4), r.T08_observed,
             num(r.pydeseq2_rotation_p, 4) if r.contrast == "OSD-38" else ""] for r in d.itertuples()]
    o = d[d.contrast == "OSD-38"].iloc[0]
    return table("tab:s7", "\\textbf{Exact sample permutation for the 3-versus-3 contrasts (T21).} All 20 labellings, "
                 "observed included, so $p\\ge0.1$. Rotation $p$ and outcome as pre-registered (T08). Last column: "
                 f"OSD-38 sets on PyDESeq2 statistics from GeneLab's raw counts (Spearman with GeneLab's "
                 f"{o.pipeline_spearman_stat:.3f} over {int(o.pipeline_genes):,} genes).".replace(",", "{,}"),
                 ["Set", "Contrast", "Stat.", "Perm. $p$", "Rotation $p$", "T08 outcome", "PyDESeq2 $p$"], rows,
                 "@{}lllrrlr@{}")


def s8():
    d = t("T22_din_genes.tsv")
    contrasts = ["OSD-522", "OSD-678-light", "OSD-427", "OSD-38", "OSD-321", "OSD-678-dark"]
    g = d[d.kind == "gene"].pivot(index="gene", columns="contrast", values="stat")[contrasts]
    rows = [[f"{gene} (Wald)", *[num(v, 1) for v in g.loc[gene]]] for gene in g.index]
    lo = d[d.kind.isin(["set", "leave one out"])].copy()
    lo["pp"] = lo.perm_p.where(lo.contrast == "OSD-522", lo.rotation_p)
    lo["gene"] = lo.gene.fillna("")
    for gene in [""] + list(g.index):
        r = lo[lo.gene == gene].set_index("contrast").pp.reindex(contrasts)
        rows.append([("all five: $p$" if gene == "" else f"without {gene}: $p$"), *[num(v, 4) for v in r]])
    s = d[d.kind == "sample"]
    for c in ("OSD-38", "OSD-678-light", "OSD-678-dark"):
        x = s[s.contrast == c].pivot(index="gene", columns="sample", values="log2_norm_count")
        cols = sorted(x.columns, key=lambda k: (not k.startswith("flight"), k))
        for gene in x.index:
            rows.append([f"{c} {gene}", *[f"{x.loc[gene, k]:.2f}" for k in cols]])
    return table("tab:s8", "\\textbf{The five DIN genes (T22).} Top: Wald statistic per contrast; then the set $p$ "
                 "(OSD-522 sample permutation, others rotation) with all five genes and with each left out. Bottom: "
                 "log$_2$(normalised count + 1) per sample in the 3-versus-3 contrasts, three flight then three ground "
                 "(OSDR sample order). In OSD-678 dark, flight replicate 3 carries the rise.",
                 ["", *contrasts], rows, "@{}lrrrrrr@{}", "\\tiny")


def s9():
    d = t("T19_osd678_interaction.tsv")
    rows = [[SETS[r.set], r.term.replace("_", " "), r.statistic, num(r.z, 1), num(r.p, 4), num(r.q, 4), r.direction,
             str(r.n_used)] for r in d.itertuples()]
    return table("tab:s9", "\\textbf{OSD-678 flight $\\times$ light (T19).} One edgeR-TMM/limma-voom model of the twelve "
                 "Col-0 samples; rotation tests on the moderated $t$ of the interaction (primary; BH over eleven sets) "
                 f"and of flight in light and in darkness (context). Marginal limma $t$ against GeneLab's Wald: Spearman "
                 f"{d.check_spearman_light.iloc[0]:.2f} (light), {d.check_spearman_dark.iloc[0]:.2f} (dark).",
                 ["Set", "Term", "Stat.", "$z$", "$p$", "$q$", "Outcome", "$n$"], rows, "@{}lllrrrlr@{}")


def s10():
    d = t("T17_paintomics_controls.tsv", dtype=str, keep_default_na=False)
    rows = [[f"{r.job} {r.job_id}".strip(), r.versus, esc(r.quantity), esc(r.value)] for r in d.itertuples()
            if r.quantity != "what"]
    return table("tab:s10", "\\textbf{PaintOmics control jobs (T17).} E: job A without metabolites; F: A with "
                 "effect-scaled metabolite values; G, H: A with genes or proteins only (plus the layer); I, J: the same "
                 "for D; K1, K2: B with its relevant genes cut at random to C's number. Also the per-omic $p$ in job A "
                 "and the relevant compounds without a hub result.",
                 ["Job", "vs", "Quantity", "Value"], rows, "@{}llp{4.2cm}p{7.6cm}@{}", "\\tiny")


def s11():
    d = t("T05_paintomics_hubs.tsv")
    names = {"C02504": "2-Isopropylmalate"}
    rows = []
    for r in d.itertuples():
        cells = []
        for j in ("A", "B", "C", "D"):
            fdr, p, step = getattr(r, f"{j}_pvalue_adjust"), getattr(r, f"{j}_pvalue"), getattr(r, f"{j}_step")
            cells += ["--"] if pd.isna(fdr) else [f"{num(p, 4)} / {num(fdr, 4)} ({int(step)})"]
        rows.append([esc(names.get(r.kegg_id, str(r.compound))[:30]), r.kegg_id, *cells])
    o = t("T18_hub_overlap_null.tsv")
    tests = o[o.kind != "compound"]
    note = "; ".join(f"{esc(x.test)}: $p={x.p:.2f}$" for x in tests.itertuples() if pd.notna(x.p))
    return table("tab:s11", "\\textbf{Metabolite hubs (T05) and the lit--dark overlap (T18).} Raw $p$ / BH FDR (best step) "
                 "per job: A 6kBOS7AN5K, B JfVZ1L37PY, C wgT3IPi7U2, D Es6j23q40o; -- not tested. Ties in FDR follow "
                 f"from the BH step-up. Overlap and reaction-count tests: {note}.",
                 ["Compound", "KEGG", "A", "B (light)", "C (dark)", "D"], rows, "@{}llllll@{}", "\\tiny")


def s12():
    d = t("T23_osd38_proteome_checks.tsv")
    rows = []
    for r in d.itertuples():
        what = esc(r.item if isinstance(r.item, str) else SETS.get(r.set, r.set))
        val = []
        for k in ("z", "q", "T12_z", "T12_observed", "rho"):
            v = getattr(r, k)
            if isinstance(v, str) or pd.notna(v):
                val.append(f"{k.replace('_', ' ')} {v if isinstance(v, str) else num(v, 3)}")
        for k in ("membrane_mean_log2fc", "soluble_mean_log2fc"):
            v = getattr(r, k)
            if pd.notna(v):
                val.append(f"{k.split('_')[0]} {num(v)}")
        rows.append([esc(r.kind), what[:110], num(int(r.n)) if pd.notna(r.n) else "", "; ".join(val)])
    return table("tab:s12", "\\textbf{OSD-38 proteome checks (T23).} limma moderated $t$ for the membrane sets next to "
                 "the pre-registered Student $t$ (T12), the membrane--soluble partition check, the truncation and filter "
                 "accounting, and OSDR's file listing (2 November 2026).",
                 ["Check", "Item", "$n$", "Result"], rows, "@{}lp{6.3cm}rp{4.6cm}@{}", "\\tiny")


def s13():
    d = t("T24_hypoxia_positive_control.tsv")
    rows = [[esc(r.dataset), SETS[r.set], r.role, num(r.z, 1), num(r.q, 4), r.direction, str(r.n_used)]
            for r in d.itertuples()]
    e = t("T25_hypoxia_signature.tsv")
    for r in e.itertuples():
        if r.kind == "signature":
            rows.append([esc(r.dataset), SETS[r.set], "addendum", num(r.z, 1), num(r.q, 4), r.direction, str(r.n_used)])
        else:
            rows.append([esc(r.dataset), r.set, "gene", f"log$_2$FC {num(r.log2fc)}", num(r.padj, 4), "", ""])
    return table("tab:s13", "\\textbf{Hypoxia positive control (T24) and signature addendum (T25).} OSD-136: Ws roots, "
                 "50\\,kPa without against with O$_2$ restored to 21\\,kPa (GeneLab limma $t$). GSE14502: Col-0 (35S line) "
                 "shoots and roots, 2\\,h hypoxia, total mRNA (limma). A set is hypoxia-sensitive if up at $q<0.05$ in a "
                 "GSE14502 tissue; the BRIC comparator result is informative only for such sets. The addendum signature "
                 "is the genes up (down) more than two-fold at $q<0.05$ in both GSE14502 tissues. Gene rows: $q$ is "
                 "GeneLab's or limma's adjusted $p$.",
                 ["Dataset", "Set", "Role", "$z$", "$q$", "Outcome", "$n$"], rows, "@{}p{4.6cm}llrrlr@{}")


def s14():
    d = t("T26_metabolights_provenance.tsv")
    rows = [[r.study, esc(r.variant), esc(r.organism_part), esc(r.group_rule), f"{r.groups} / {r.samples} ({r.n_per_group})",
             esc(f"{r.chromatography}; {r.instrument}, {r.mass_analyzer}"), esc(r.compounds_whose_study_max_is_100),
             esc(r.other_organisms_in_study)] for r in d.itertuples()]
    return table("tab:s14", "\\textbf{MetaboLights variance sources (T26)}, from each study's ISA sample and assay files. "
                 "Groups / samples (per group). ``Max 100'': compounds whose maximum over the study is exactly 100, i.e. "
                 "values reported as percent of maximum.",
                 ["Study", "Variant", "Tissue", "Groups used", "Groups / n", "Platform", "Max 100", "Other organisms"],
                 rows, "@{}lp{1.9cm}p{1.1cm}p{3.1cm}p{1.2cm}p{2.6cm}rp{1.2cm}@{}", "\\tiny")


def main():
    parts = ["% Generated by scripts/29_supplementary_tables.py from results/tables. Do not edit by hand.", "",
             s2(), s3(), s4(), s5(), s6(), s7(), s8(), s9(), s10(), s11(), s12(), s13(), s14()]
    OUT.write_text("\n".join(parts))
    print(f"  supp_tables.tex: tables S2-S14 ({OUT.stat().st_size:,} bytes)")


if __name__ == "__main__":
    main()
