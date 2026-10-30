"""Write manuscript/latex/numbers.tex: every number the manuscript quotes, read from results/tables.

The prose uses these macros and never types a result by hand, so a re-run that changes a table
changes the text with it. Each macro's source table is named in the comment beside it.
"""

import numpy as np
import pandas as pd

from paths import DATA, PHOTORESP, ROOT, TABLES

OUT = ROOT / "manuscript" / "latex" / "numbers.tex"


def t(name, **kw):
    return pd.read_csv(TABLES / name, sep="\t", comment="#", **kw)


def p(v, rotation=True):
    """p- and q-values with their relation sign ("= 0.0281", "\\le 10^{-4}"), for use as $p \\pX$.
    A rotation p at its floor, 1/(10000 + 1), is reported as the bound it is."""
    if rotation and v <= 1 / 10001 + 1e-12:
        return "\\le 10^{-4}"
    return "= " + _p(v)


def _p(v):
    if v < 1e-3:
        m, e = f"{v:.2e}".split("e")
        return f"{m}\\times10^{{{int(e)}}}"
    return f"{v:.3g}"


def main():
    m = {}
    src = t("T01_light_dark_sources.tsv")
    ref = t("T02_light_dark_reference.tsv")
    m["nSourceRows"] = (len(src), "T01")
    m["nSourceKegg"] = (int((src.kegg_id.notna()).sum()), "T01")
    m["nZhuAcc"] = (int(src[src.source == "zhu"].n.max()), "T01")
    m["nWuAcc"] = (int(src[src.source == "wu_stress"].n.max()), "T01")
    m["nCaldanaMin"] = (int(src[src.source == "caldana"].n.min()), "T01")
    m["nCaldanaMax"] = (int(src[src.source == "caldana"].n.max()), "T01")
    m["nRefCompounds"] = (len(ref), "T02")
    for k in ("robust", "single", "none", "conflicting"):
        m[f"nRef{k.capitalize()}"] = (int((ref.status == k).sum()), "T02")
    m["nRefUp"] = (int((ref.direction > 0).sum()), "T02")
    m["nRefDown"] = (int((ref.direction < 0).sum()), "T02")

    def eff(source, kid, col="effect"):
        r = src[(src.source == source) & (src.kegg_id == kid)]
        return float(r[col].iloc[0])
    m["leuZhu"] = (f"{eff('zhu', 'C00123'):.2f}", "T01")
    m["leuWu"] = (f"{eff('wu_stress', 'C00123'):.2f}", "T01")
    m["leuCaldana"] = (f"{eff('caldana', 'C00123'):.2f}", "T01")
    m["ileZhu"] = (f"{eff('zhu', 'C00407'):.2f}", "T01")
    m["ileWu"] = (f"{eff('wu_stress', 'C00407'):.2f}", "T01")
    m["gluWu"] = (f"{eff('wu_stress', 'C00025'):.2f}", "T01")
    m["gluWuFrac"] = (f"{100 * eff('wu_stress', 'C00025', 'frac_same'):.0f}", "T01")
    m["glnWu"] = (f"{eff('wu_stress', 'C00064'):.2f}", "T01")
    m["glnWuFrac"] = (f"{100 * eff('wu_stress', 'C00064', 'frac_same'):.0f}", "T01")
    m["serZhu"] = (f"{eff('zhu', 'C00065'):.2f}", "T01")
    m["serZhuFrac"] = (f"{100 * eff('zhu', 'C00065', 'frac_same'):.0f}", "T01")
    m["malZhu"] = (f"{eff('zhu', 'C00149'):.2f}", "T01")
    m["oxgZhu"] = (f"{eff('zhu', 'C00026'):.2f}", "T01")
    m["asnCaldana"] = (f"{eff('caldana', 'C00152'):.2f}", "T01")

    z = ref.set_index("kegg_id").reference_z
    m["ratioLeuVal"] = (f"{z['C00123'] / z['C00183']:.1f}", "T02")
    m["ratioLeuAsn"] = (f"{z['C00123'] / z['C00152']:.1f}", "T02")
    m["nModelCompounds"] = (len(pd.read_csv(PHOTORESP / "T06_predicted_compounds.tsv", sep="\t", comment="#")),
                            "first study T06")
    cmp_ = t("T03_model_vs_reference.tsv")
    mm = cmp_[cmp_.in_model]
    m["nComparable"] = (int(mm.verdict.isin(["agrees", "disagrees"]).sum()), "T03")
    m["nScored"] = (int(mm.verdict.isin(["agrees", "disagrees", "untested"]).sum()), "T03")
    for k, name in (("agrees", "Agree"), ("disagrees", "Disagree"), ("untested", "Untested"),
                    ("not comparable", "NotComparable")):
        m[f"nModel{name}"] = (int((mm.verdict == k).sum()), "T03")
    m["nNotPredicted"] = (int((~cmp_.in_model).sum()), "T03")
    m["nNotPredictedRobust"] = (int(((~cmp_.in_model) & (cmp_.status == "robust")).sum()), "T03")

    gs = t("T08_geneset_tests.tsv")
    m["nTests"] = (len(gs), "T08")
    m["nMatch"] = (int((gs.verdict == "MATCH").sum()), "T08")
    m["nLitMatch"] = (int(((gs.light == "lit") & (gs.verdict == "MATCH")).sum()), "T08")
    m["nDarkMatch"] = (int(((gs.light == "dark") & (gs.verdict == "MATCH")).sum()), "T08")
    m["nLitTests"] = (int((gs.light == "lit").sum()), "T08")
    m["nSigLitPredicted"] = (int(((gs.light == "lit") & (gs.q < 0.05) & (gs.verdict == "MATCH")
                                  & (gs.predicted != "none")).sum()), "T08")
    g = lambda s, c: gs[(gs.set == s) & (gs.contrast == c)].iloc[0]
    for key, s, c in (("DinFiveTwoTwo", "din_first_study", "OSD-522"), ("DinLight", "din_first_study", "OSD-678-light"),
                      ("DinFourTwoSeven", "din_first_study", "OSD-427"), ("BcaaLight", "bcaa_degradation", "OSD-678-light"),
                      ("DinThirtyEight", "din_first_study", "OSD-38"), ("DinThreeTwoOne", "din_first_study", "OSD-321"),
                      ("DinDark", "din_first_study", "OSD-678-dark"), ("BcaaThreeTwoOne", "bcaa_degradation", "OSD-321"),
                      ("BcaaDark", "bcaa_degradation", "OSD-678-dark"), ("SharedLight", "wu_shared", "OSD-678-light"),
                      ("DarkMetFiveTwoTwo", "wu_dark_metabolite", "OSD-522")):
        r = g(s, c)
        m[f"p{key}"] = (p(r.p), "T08")
        m[f"q{key}"] = (p(r.q), "T08")
        m[f"z{key}"] = (f"{r.z:.1f}".replace("-", "-"), "T08")
    r = g("din_first_study", "OSD-522")
    m["pDinFiveTwoTwoMW"] = (p(r.mw_p), "T08")
    m["pDinFiveTwoTwoRot"] = (p(r.rot_p_signed), "T08")
    first = pd.read_csv(PHOTORESP / "T07_falsification_osd522.tsv", sep="\t")
    m["pDinFirstStudy"] = (p(float(first[(first.set == "carbon_starvation_DIN") &
                                           (first.layer == "transcriptome")].p.iloc[0])), "first study T07")
    sets = pd.read_csv(ROOT / "prereg" / "gene_sets.tsv", sep="\t").groupby("set").size()
    for k, name in (("din_first_study", "Din"), ("bcaa_degradation", "Bcaa"), ("wu_control_only", "WuControl"),
                    ("wu_stress_only", "WuStress"), ("wu_shared", "WuShared"), ("wu_dark_metabolite", "WuDarkMet")):
        m[f"nSet{name}"] = (f"{int(sets[k]):,}".replace(",", "{,}"), "prereg/gene_sets.tsv")

    summ = t("T07_paintomics_summary.tsv", nrows=4).set_index("quantity")
    for job, key in (("job1", "One"), ("A", "A"), ("B", "B"), ("C", "C")):
        m[f"nPath{key}"] = (int(summ.loc["pathways at BH q < 0.05", job]), "T07")
        m[f"nHub{key}"] = (int(summ.loc["hub compound at FDR < 0.05".replace("compound", "compounds"), job]), "T07")
        m[f"nHubTested{key}"] = (int(summ.loc["hub compounds tested", job]), "T07")
    path = t("T04_paintomics_pathways.tsv")
    from scipy import stats
    m["rhoOneA"] = (f"{stats.spearmanr(path.job1_p, path.A_p).correlation:.2f}", "T04")
    m["rhoBC"] = (f"{stats.spearmanr(path.B_p, path.C_p).correlation:.2f}", "T04")
    m["nPathTotal"] = (len(path), "T04")
    row = lambda name: path[path.pathway == name].iloc[0]
    gly = row("Glyoxylate and dicarboxylate metabolism")
    m["rankGlyOne"] = (int(gly.job1_rank), "T04")
    m["rankGlyA"] = (int(gly.A_rank), "T04")
    m["pPhotoC"] = (p(row("Photosynthesis").C_p, rotation=False), "T04")
    m["qPhotoB"] = (p(row("Photosynthesis").B_q, rotation=False), "T04")
    m["pBcaaDegC"] = (p(row("Valine, leucine and isoleucine degradation").C_p, rotation=False), "T04")
    m["pBcaaDegB"] = (p(row("Valine, leucine and isoleucine degradation").B_p, rotation=False), "T04")
    m["pErA"] = (p(row("Protein processing in endoplasmic reticulum").A_p, rotation=False), "T04")
    hubs = t("T05_paintomics_hubs.tsv")
    m["hubLightFDR"] = (f"{hubs.B_pvalue_adjust[hubs.B_pvalue_adjust < 0.05].max():.3f}", "T05")
    row = [ln for ln in (TABLES / "T07_paintomics_summary.tsv").read_text().splitlines() if "identical" in ln][0]
    assert "\tTrue (" in row, row
    m["nSharedHubSteps"] = (int(row.split("(")[-1].split()[0]), "T07")
    up = pd.read_csv(DATA / "paintomics_upload" / "metabolomics_key.tsv", sep="\t")
    m["nLayer"] = (len(up), "data/paintomics_upload/metabolomics_key.tsv")
    m["nLayerRelevant"] = (int(up.relevant.sum()), "data/paintomics_upload/metabolomics_key.tsv")
    for c, key in (("OSD-678-light", "Light"), ("OSD-678-dark", "Dark")):
        d = pd.read_csv(DATA / "osdr" / f"contrast_{c}.csv.gz", usecols=["padj"])
        m[f"nDE{key}"] = (f"{int((d.padj < 0.05).sum()):,}".replace(",", "{,}"), f"data/osdr/contrast_{c}")

    de = [int((pd.read_csv(DATA / "osdr" / f"contrast_{c}.csv.gz", usecols=["padj"]).padj < 0.05).sum())
          for c in ("OSD-678-light", "OSD-678-dark")]
    m["ratioDE"] = (f"{de[0] / de[1]:.1f}", "data/osdr")
    pw = t("T10_gcms_power.tsv", index_col=0)
    grp = pd.read_csv(DATA / "reference" / "metabolights_group_sds.tsv", sep="\t")
    m["nPanel"] = (len(pw), "T10")
    m["nPanelGroups"] = (int(grp[["study", "group"]].drop_duplicates().shape[0]), "data/reference/metabolights_group_sds.tsv")
    m["nPanelStudies"] = (int(grp.study.nunique()), "data/reference/metabolights_group_sds.tsv")
    m["sdMin"] = (f"{pw.sd_log2.min():.2f}", "T10")
    m["sdMax"] = (f"{pw.sd_log2.max():.2f}", "T10")
    m["nQThousandMin"] = (int(pw.n_Q1000.min()), "T10")
    m["nQThousandMax"] = (f"{int(pw.n_Q1000.max()):,}".replace(",", "{,}"), "T10")
    m["nQSixtyMin"] = (f"{int(pw.n_Q60.min()):,}".replace(",", "{,}"), "T10")
    m["detectMin"] = (f"{pw.detectable_log2_at_n6.min():.2f}", "T10")
    m["detectMinFold"] = (f"{2 ** pw.detectable_log2_at_n6.min():.2f}", "T10")
    m["nTwoFoldMax"] = (int(pw.n_two_fold.max()), "T10")
    m["qSixtyLogTwo"] = (f"{abs(np.log2(1 - 0.0237)):.3f}", "first study T13")

    # OSD-38 (BRIC-20) addendum
    prot = t("T11_osd38_proteome.tsv")
    m["nProtMem"] = (f"{int((prot.fraction == 'membrane').sum()):,}".replace(",", "{,}"), "T11")
    m["nProtSol"] = (int((prot.fraction == "soluble").sum()), "T11")
    g38 = t("T12_osd38_geneset_tests.tsv")
    r38 = lambda layer, s_: g38[(g38.layer == layer) & (g38.set == s_)].iloc[0]
    for layer, lk in (("rna_flight", "Rna"), ("membrane_flight", "Mem"), ("rna_hardware", "Hw")):
        for s_, sk in (("ath00630", "Gly"), ("ath00710", "Fix"), ("ath00500", "Sta"), ("ath00010", "Glyc"),
                       ("bcaa_degradation", "Bcaa"), ("photosynthesis_apparatus", "Photo"), ("din_first_study", "Din"),
                       ("hypoxia_responsive", "Hyp"), ("fermentation", "Fer")):
            r = r38(layer, s_)
            if pd.notna(r.p):
                m[f"z{lk}{sk}"] = (f"{r.z:.1f}", "T12")
                m[f"q{lk}{sk}"] = (p(r.q), "T12")
            m[f"n{lk}{sk}"] = (int(r.n_measured), "T12")
    mem_carbon = g38[(g38.layer == "membrane_flight") & g38.set.isin(["ath00630", "ath00710", "ath00500", "ath00010"])]
    m["qMemCarbonMin"] = (_p(mem_carbon.q.min()), "T12")
    m["qMemCarbonMax"] = (_p(mem_carbon.q.max()), "T12")
    rna_carbon = g38[(g38.layer == "rna_flight") & g38.set.isin(["ath00630", "ath00710", "ath00500", "ath00010"])]
    m["qRnaCarbonMax"] = (_p(rna_carbon.q.max()), "T12")
    m["nMemCarbonDown"] = (int(((mem_carbon.q < 0.05) & (mem_carbon.signed < 0)).sum()), "T12")
    m["nRnaCarbonUp"] = (int(((rna_carbon.q < 0.05) & (rna_carbon.signed > 0)).sum()), "T12")
    prim = g38[g38.primary & g38.p.notna()]
    m["nOsdTests"] = (len(prim), "T12")
    m["nOsdMatch"] = (int((prim.verdict == "MATCH").sum()), "T12")
    c38 = t("T13_osd38_concordance.tsv")
    row = lambda pat: c38[c38.comparison.str.contains(pat, regex=False)].iloc[0]
    m["rhoMem"] = (f"{row('OSD-38 membrane protein vs transcript').rho:.2f}", "T13")
    m["rhoSol"] = (f"{row('OSD-38 soluble protein vs transcript').rho:.2f}", "T13")
    m["rhoFiveTwoTwoProt"] = (f"{row('OSD-522 protein vs transcript').rho:.2f}", "T13")
    m["nFiveTwoTwoProt"] = (f"{int(row('OSD-522 protein vs transcript').n):,}".replace(",", "{,}"), "T13")
    m["nPaperRule"] = (int(row("paper's rule").n), "T13")
    m["rhoFlightHw"] = (f"{row('against BRIC vs container').rho:.2f}", "T13")
    m["rhoPsm"] = (f"{row('Spearman of protein t with log PSM count').rho:.2f}", "T13")
    rob = c38[c38.comparison.str.contains("PSM-decile-matched") & ~c38.comparison.str.contains("bcaa")]
    m["pRobMin"] = (_p(rob.p.min()), "T13")
    m["pRobMax"] = (_p(rob.p.max()), "T13")
    cov = pd.read_csv(DATA / "reference" / "osd38_protein_coverage.tsv", sep="\t")
    src = pd.read_csv(DATA / "reference" / "osd38_sources.tsv", sep="\t")
    m["nOsdRows"] = (f"{int(src.data_rows.iloc[0]):,}".replace(",", "{,}"), "data/reference/osd38_sources.tsv")
    summ5 = t("T07_paintomics_summary.tsv", nrows=4).set_index("quantity")
    m["nPathD"] = (int(summ5.loc["pathways at BH q < 0.05", "D"]), "T07")
    m["nHubD"] = (int(summ5.loc["hub compounds at FDR < 0.05", "D"]), "T07")
    t7rows = dict(ln.split("\t", 1) for ln in (TABLES / "T07_paintomics_summary.tsv").read_text().splitlines()
                  if ln.startswith("A vs D"))
    m["rhoAD"] = (f"{float(t7rows['A vs D (lit BRIC-LED vs dark BRIC-PDFU): Spearman rho of combined p']):.2f}", "T07")
    hubs5 = t("T05_paintomics_hubs.tsv")
    both = hubs5[(hubs5.A_pvalue_adjust < 0.05) & (hubs5.D_pvalue_adjust < 0.05)]
    m["nHubAD"] = (len(both), "T05")
    m["nPathAD"] = (int(((path.A_q < 0.05) & (path.D_q < 0.05)).sum()), "T04")
    dge = pd.read_csv(DATA / "osdr" / "contrast_OSD-38.csv.gz", usecols=["padj"])
    m["nDEThirtyEight"] = (f"{int((dge.padj < 0.05).sum()):,}".replace(",", "{,}"), "data/osdr/contrast_OSD-38")

    lines = ["% Generated by scripts/18_manuscript_numbers.py. Do not edit by hand."]
    for k, (v, s) in m.items():
        lines.append(f"\\newcommand{{\\{k}}}{{{v}}}  % {s}")
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text("\n".join(lines) + "\n")
    print(f"  numbers.tex: {len(m)} macros")


if __name__ == "__main__":
    main()
