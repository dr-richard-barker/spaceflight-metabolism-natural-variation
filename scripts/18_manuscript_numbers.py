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

    revision(m)
    lines = ["% Generated by scripts/18_manuscript_numbers.py. Do not edit by hand."]
    for k, (v, s) in m.items():
        lines.append(f"\\newcommand{{\\{k}}}{{{v}}}  % {s}")
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text("\n".join(lines) + "\n")
    print(f"  numbers.tex: {len(m)} macros")


def num(v):
    return f"{int(v):,}".replace(",", "{,}")


def kv(name):
    """A long-format table's (job, versus, quantity) -> value lookup."""
    d = t(name, dtype=str, keep_default_na=False)
    return {(r.job, r.versus, r.quantity): r.value for r in d.itertuples()}


def revision(m):
    """Revision 1 (post hoc, prereg/revision1): T14-T26."""
    # R1 reference strata, sensitivity, KEGG mapping, coverage
    s14 = t("T14_reference_strata.tsv")
    d = s14[s14.status_all.isin(["robust", "single"])]
    m["nDirectional"] = (len(d), "T14")
    m["nDarkHeatOnly"] = (int(d.dark_heat_only.sum()), "T14")
    m["nFlipDarkOnly"] = (int(s14.flips_dark_only.sum()), "T14")
    for k in ("robust", "single", "absent"):
        m[f"nDarkOnly{k.capitalize()}"] = (int((s14.status_dark_only == k).sum()), "T14")
    glu = s14.set_index("kegg_id").loc["C00025"]
    m["gluCaldana"] = (f"{glu.caldana_effect:.2f}", "T14")
    m["gluCaldanaP"] = (f"{glu.caldana_support:.2f}", "T14")
    g15 = t("T15_reference_sensitivity.tsv")
    m["nGrid"] = (int((g15.stratum == "all").sum()), "T15")
    m["nGridRows"] = (len(g15), "T15")
    m["nGridBcaaUp"] = (int(g15.bcaa_up.sum()), "T15")
    m["nGridMalOxgDown"] = (int(g15.malate_2og_down.sum()), "T15")
    for st, key in (("all", "All"), ("dark_only", "Dark")):
        g = g15[g15.stratum == st]
        m[f"robust{key}Min"], m[f"robust{key}Max"] = (int(g.robust.min()), "T15"), (int(g.robust.max()), "T15")
        m[f"notPred{key}Min"] = (int(g.robust_not_predicted.min()), "T15")
        m[f"notPred{key}Max"] = (int(g.robust_not_predicted.max()), "T15")
        m[f"nGluDisagree{key}"] = (int((g.model_disagrees > 0).sum()), "T15")
        m[f"nGluUp{key}"] = (int(g.Glu.str.endswith(" up").sum()), "T15")
    k16 = t("T16_kegg_mapping_report.tsv").fillna("")
    v = lambda sec, item: int(k16[(k16.section == sec) & (k16.item == item)].n.iloc[0])
    o = lambda sec, item: int(float(k16[(k16.section == sec) & (k16.item == item)]["of"].iloc[0]))
    m["nKeggNames"] = (v("names", "measured names sent to KEGG (Zhu and Wu)"), "T16")
    for lab, key in (("0", "Zero"), ("1", "One"), (">1", "Multi")):
        m[f"nKegg{key}"] = (v("KEGG find candidates", f"{lab} candidate(s)"), "T16")
    m["nKeggRule"] = (v("rule match (exact synonym)", "mapped by the rule"), "T16")
    m["nOverrides"] = (v("reviewed overrides", "rows in data/reference/kegg_overrides.tsv"), "T16")
    m["nOverRescued"] = (v("reviewed overrides", "rescued (rule found no exact synonym)"), "T16")
    m["nOverCorrected"] = (v("reviewed overrides", "corrected (rule picked another entry)"), "T16")
    m["nKeggMapped"] = (v("after overrides", "measured names with a KEGG ID"), "T16")
    m["nOverCaldanaFix"] = (v("reviewed overrides", "Caldana peak only: Wu's given ID corrected"), "T16")
    m["nOverCaldanaExcl"] = (v("reviewed overrides", "Caldana peak only: excluded"), "T16")
    for src_, key in (("zhu", "Zhu"), ("wu_stress", "Wu")):
        m[f"nKeggMapped{key}"] = (v("after overrides", f"from {src_}"), "T16")
        m[f"nKeggNames{key}"] = (o("after overrides", f"from {src_}"), "T16")
    for cls, key in (("Amino acid", "AminoAcid"), ("Organic acid", "OrganicAcid"), ("Glucosinolate", "Glucosinolate"),
                     ("Flavonoid", "Flavonoid"), ("Phenylpropanoid", "Phenylpropanoid")):
        m[f"wuClass{key}"] = (f"{v('Wu class (after overrides)', cls)} of {o('Wu class (after overrides)', cls)}", "T16")
    src = t("T01_light_dark_sources.tsv")
    wu = src[src.source == "wu_stress"]
    m["wuMissMedian"] = (f"{100 * wu.missing_frac.median():.1f}", "T01")
    m["wuMissMax"] = (f"{100 * wu.missing_frac.max():.0f}", "T01")
    m["wuNMin"] = (int(wu.n.min()), "T01")
    m["nWuPossible"] = (int(wu.n_possible.max()), "T01")
    c3 = t("T03_model_vs_reference.tsv")
    mm = c3[c3.in_model]
    m["nNotMeasured"] = (int((mm.verdict_detail == "not measured").sum()), "T03")
    m["nNoDirection"] = (int((mm.verdict_detail == "measured, no directional response").sum()), "T03")
    ref = t("T02_light_dark_reference.tsv")
    robust = set(ref.kegg_id[ref.status == "robust"])
    m["nRobustCovered"] = (len(robust & set(mm.kegg_id)), "T02, T03")

    # R2 gene-set diagnostics
    c20 = t("T20_rotation_calibration.tsv")
    blk, rnd = "contiguous genomic blocks", "random genes"
    sel = lambda kind, null: c20[(c20.set_type == kind) & (c20.null == null)].type1_error
    m["typeOneRotRandom"] = (f"{sel(rnd, 'rotation').median():.3f}", "T20")
    m["typeOneRotBlock"] = (f"{sel(blk, 'rotation').median():.3f}", "T20")
    m["typeOneRotMin"] = (f"{min(sel(rnd, 'rotation').min(), sel(blk, 'rotation').min()):.3f}", "T20")
    m["typeOneRotMax"] = (f"{max(sel(rnd, 'rotation').max(), sel(blk, 'rotation').max()):.3f}", "T20")
    m["typeOneNaiveBlock"] = (f"{sel(blk, 'naive').median():.3f}", "T20")
    m["typeOneNaiveBlockMax"] = (f"{sel(blk, 'naive').max():.3f}", "T20")
    m["nCalSets"] = (num(c20.n_sets.sum() // 2), "T20")
    e21 = t("T21_exact_permutation.tsv")
    r = lambda s_, c: e21[(e21.set == s_) & (e21.contrast == c)].iloc[0]
    m["permDinThirtyEight"] = (f"{r('din_first_study', 'OSD-38').perm_p:.1f}", "T21")
    m["permDinLight"] = (f"{r('din_first_study', 'OSD-678-light').perm_p:.1f}", "T21")
    m["permDinDark"] = (f"{r('din_first_study', 'OSD-678-dark').perm_p:.1f}", "T21")
    m["permBcaaDark"] = (f"{r('bcaa_degradation', 'OSD-678-dark').perm_p:.1f}", "T21")
    m["permSharedLight"] = (f"{r('wu_shared', 'OSD-678-light').perm_p:.1f}", "T21")
    o38 = r("din_first_study", "OSD-38")
    m["pipeRhoStat"] = (f"{o38.pipeline_spearman_stat:.3f}", "T21")
    m["pipeRhoFc"] = (f"{o38.pipeline_spearman_log2fc:.3f}", "T21")
    m["pipeGenes"] = (num(o38.pipeline_genes), "T21")
    m["pDinThirtyEightDeseq"] = (p(o38.pydeseq2_rotation_p), "T21")
    d22 = t("T22_din_genes.tsv")
    lo = d22[d22.kind == "leave one out"].copy()
    lo["pp"] = lo.perm_p.where(lo.contrast == "OSD-522", lo.rotation_p)
    for c, key in (("OSD-38", "ThirtyEight"), ("OSD-321", "ThreeTwoOne"), ("OSD-678-dark", "Dark")):
        m[f"looMax{key}"] = (_p(lo[lo.contrast == c].pp.max()), "T22")

    # R3 OSD-678 interaction
    i19 = t("T19_osd678_interaction.tsv")
    r = lambda s_, term="interaction": i19[(i19.set == s_) & (i19.term == term)].iloc[0]
    for s_, key in (("din_first_study", "Din"), ("bcaa_degradation", "Bcaa"), ("photosynthesis_apparatus", "Photo"),
                    ("ath00710", "Fix"), ("ath00630", "Gly"), ("ath00010", "Glyc"), ("ath00500", "Sta")):
        m[f"zInt{key}"] = (f"{r(s_).z:.1f}", "T19")
        m[f"qInt{key}"] = (p(r(s_).q), "T19")
    m["zLimmaDinDark"] = (f"{r('din_first_study', 'flight_dark').z:.1f}", "T19")
    m["zLimmaPhotoDark"] = (f"{r('photosynthesis_apparatus', 'flight_dark').z:.1f}", "T19")
    m["zLimmaPhotoLight"] = (f"{r('photosynthesis_apparatus', 'flight_light').z:.1f}", "T19")
    m["limmaRhoLight"] = (f"{i19.check_spearman_light.iloc[0]:.2f}", "T19")
    m["limmaRhoDark"] = (f"{i19.check_spearman_dark.iloc[0]:.2f}", "T19")
    m["nIntSig"] = (int(((i19.term == "interaction") & (i19.q < 0.05)).sum()), "T19")

    # R4 hypoxia positive control, R4b signature
    h24 = t("T24_hypoxia_positive_control.tsv")
    r = lambda ds, s_: h24[(h24.dataset == ds) & (h24.set == s_)].iloc[0]
    m["zFerOsdOneThreeSix"] = (f"{r('OSD-136', 'fermentation').z:.1f}", "T24")
    m["zFerShoot"] = (f"{r('GSE14502_shoot', 'fermentation').z:.1f}", "T24")
    m["zFerRoot"] = (f"{r('GSE14502_root', 'fermentation').z:.1f}", "T24")
    m["nHypMeasured"] = (int(r("GSE14502_shoot", "hypoxia_responsive").n_used), "T24")
    m["nHypSet"] = (int(r("GSE14502_shoot", "hypoxia_responsive").n_set), "T24")
    m["zDinRoot"] = (f"{r('GSE14502_root', 'din_first_study').z:.1f}", "T24")
    h25 = t("T25_hypoxia_signature.tsv")
    sg = h25[h25.kind == "signature"]
    r = lambda ds, s_: sg[(sg.dataset == ds) & (sg.set == s_)].iloc[0]
    m["nHypUp"] = (int(r("OSD-136", "hypoxia_up").n_set), "T25")
    m["nHypDown"] = (int(r("OSD-136", "hypoxia_down").n_set), "T25")
    m["zHypUpOsdOneThreeSix"] = (f"{r('OSD-136', 'hypoxia_up').z:.1f}", "T25")
    hw = "OSD-38 BRIC-PDFU vs closed container"
    m["zHypUpHw"] = (f"{r(hw, 'hypoxia_up').z:.1f}", "T25")
    m["qHypUpHw"] = (p(r(hw, "hypoxia_up").q), "T25")
    m["zHypDownHw"] = (f"{r(hw, 'hypoxia_down').z:.1f}", "T25")
    m["zHypUpFlight"] = (f"{r('OSD-38 flight vs ground', 'hypoxia_up').z:.1f}", "T25")
    m["qHypUpFlight"] = (p(r("OSD-38 flight vs ground", "hypoxia_up").q), "T25")
    gn = h25[(h25.kind == "gene") & (h25.dataset == hw)].set_index("set")
    for g_, key in (("ADH1", "Adh"), ("PDC1", "Pdc")):
        m[f"fold{key}Hw"] = (f"{2 ** gn.loc[g_, 'log2fc']:.1f}", "T25")
        m[f"padj{key}Hw"] = (p(gn.loc[g_, "padj"], rotation=False), "T25")

    # R5 OSD-38 proteome checks
    p23 = t("T23_osd38_proteome_checks.tsv")
    lim = p23[p23.kind == "membrane set, limma moderated t"].set_index("set")
    carbon = lim.loc[["ath00630", "ath00710", "ath00500", "ath00010"]]
    m["qLimmaCarbonMin"] = (_p(carbon.q.min()), "T23")
    m["qLimmaCarbonMax"] = (_p(carbon.q.max()), "T23")
    m["zLimmaBcaa"] = (f"{lim.loc['bcaa_degradation', 'z']:.1f}", "T23")
    m["rhoLimmaStudent"] = (f"{p23[p23.kind == 'membrane, limma vs Student t'].rho.iloc[0]:.2f}", "T23")
    part = p23[p23.kind == "partition"]
    allp = part[part.item.str.startswith("proteins quantified in both")].iloc[0]
    m["rhoPartition"] = (f"{allp.rho:.2f}", "T23")
    m["nPartition"] = (int(allp.n), "T23")
    sol = p23[p23.kind == "soluble members"]
    m["solCarbonFcMin"] = (f"{sol.soluble_mean_log2fc.min():.2f}", "T23")
    m["solCarbonFcMax"] = (f"{sol.soluble_mean_log2fc.max():.2f}", "T23")
    acc = p23[p23.kind == "accounting"]
    a_ = lambda frac, item: int(acc[(acc.fraction == frac) & (acc.item.str.startswith(item))].n.iloc[0])
    m["nSolGroups"] = (a_("soluble", "protein groups kept"), "T23")
    m["nSolLoci"] = (a_("soluble", "AGI loci in the kept groups"), "T23")
    m["nSolTwoGroups"] = (a_("soluble", "loci that appear in more than one"), "T23")
    m["nMemGroups"] = (num(a_("membrane", "protein groups in block 1")), "T23")
    m["nMemGroupMax"] = (num(a_("membrane", "highest protein group number")), "T23")

    # R6 PaintOmics controls and hub overlap
    c17 = kv("T17_paintomics_controls.tsv")
    pair = lambda job, par, q: c17[(job, par, q)].split("; ")
    m["nPathE"] = (int(pair("E", "A", "pathways at BH q < 0.05 (job; parent)")[0]), "T17")
    m["rhoAE"] = (f"{float(c17[('E', 'A', 'Spearman rho of combined p')]):.2f}", "T17")
    m["rhoAF"] = (f"{float(c17[('F', 'A', 'Spearman rho of combined p')]):.2f}", "T17")
    m["nHubF"] = (int(pair("F", "A", "hubs at FDR < 0.05 (job; parent)")[0]), "T17")
    m["nSharedHubStepsF"] = (int(c17[("F", "A", "raw hub p identical to the parent for shared compound-steps")]
                                 .split("(")[1].split()[0]), "T17")
    assert c17[("F", "A", "raw hub p identical to the parent for shared compound-steps")].startswith("True")
    for job, par in (("G", "A"), ("H", "A"), ("I", "D"), ("J", "D"), ("K1", "B"), ("K2", "B")):
        key = job.replace("1", "One").replace("2", "Two")
        m[f"nHub{key}"] = (int(pair(job, par, "hubs at FDR < 0.05 (job; parent)")[0]), "T17")
        m[f"nPath{key}"] = (int(pair(job, par, "pathways at BH q < 0.05 (job; parent)")[0]), "T17")
    for om, key in (("Gene expression", "Gene"), ("Proteomics", "Prot"), ("Metabolomics", "Met")):
        m[f"rhoComb{key}"] = (f"{float(c17[('A', '', f'Spearman rho of combined p with the {om} p')]):.2f}", "T17")
    m["nPathMetP"] = (int(c17[("A", "", "pathways with a Metabolomics p")]), "T17")
    m["nRelevantTested"] = (int(c17[("A", "", "relevant compounds hub-tested")]), "T17")
    m["nRelevantNoHub"] = (int(c17[("A", "", "relevant compounds (robust) in the layer")]) - m["nRelevantTested"][0], "T17")
    m["nWuLoci"] = (len(pd.read_csv(DATA / "atlas" / "wu_gwas_loci.csv")), "data/atlas/wu_gwas_loci.csv")
    m["nWuRegions"] = (len(pd.read_csv(DATA / "atlas" / "wu_condition_regions.tsv", sep="\t")),
                       "data/atlas/wu_condition_regions.tsv")
    zhu = set(pd.read_csv(DATA / "atlas" / "darkness_metabolome.csv", usecols=["accession"]).accession)
    wua = set(pd.read_csv(DATA / "atlas" / "wu_leaf_identified_long.csv", usecols=["accession"]).accession)
    m["nAccShared"] = (len(zhu & wua), "data/atlas darkness_metabolome, wu_leaf_identified_long")
    m["nAccUnion"] = (len(zhu | wua), "data/atlas darkness_metabolome, wu_leaf_identified_long")
    m["nAccSum"] = (m["nZhuAcc"][0] + m["nWuAcc"][0], "T01")
    rng = c17[("B", "", "hubs at FDR < 0.05: raw p range")].split("-")
    m["hubRawMinB"], m["hubRawMaxB"] = (_p(float(rng[0])), "T17"), (_p(float(rng[1])), "T17")
    t18 = t("T18_hub_overlap_null.tsv")
    hg = t18[t18.test == "hypergeometric over compounds hub-tested in both"].iloc[0]
    pm = t18[t18.test.fillna("").str.startswith("hub labels shuffled")].iloc[0]
    m["nHubCand"] = (int(hg.n), "T18")
    m["hubOverlapExp"] = (f"{hg.expected:.1f}", "T18")
    m["hubOverlapP"] = (f"{hg.p:.2f}", "T18")
    m["hubOverlapExpPerm"] = (f"{pm.expected:.1f}", "T18")
    m["hubOverlapPPerm"] = (f"{pm.p:.2f}", "T18")
    mw = t18[t18.kind == "reaction count"].dropna(subset=["hubs"])
    m["mwDegreePMin"] = (f"{mw.p.min():.2f}", "T18")

    # R7 power
    pw = t("T10_gcms_power.tsv", index_col=0)
    m["detectMinBonf"] = (f"{pw.detectable_log2_at_n6_bonf.min():.2f}", "T10")
    m["detectMinFoldBonf"] = (f"{2 ** pw.detectable_log2_at_n6_bonf.min():.2f}", "T10")
    m["nQThousandBonfMin"] = (num(pw.n_Q1000_bonf.min()), "T10")
    m["nQThousandBonfMax"] = (num(pw.n_Q1000_bonf.max()), "T10")
    m["nTwoFoldBonfMax"] = (int(pw.n_two_fold_bonf.max()), "T10")
    m["nFoldOnePFiveMax"] = (int(pw.n_fold_1p5.max()), "T10")
    prov = t("T26_metabolights_provenance.tsv").set_index("study")
    m["nMtblsThreeSixThree"] = (int(prov.loc["MTBLS363", "samples"]), "T26")
    m["nMtblsThreeSixThreeGroups"] = (int(prov.loc["MTBLS363", "groups"]), "T26")


if __name__ == "__main__":
    main()
