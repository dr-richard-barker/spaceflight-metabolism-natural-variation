"""How many plants a targeted GC-MS panel on flight material would need (T10, fig05).

Variance comes from Col-0 wild-type control groups in MetaboLights GC-MS leaf studies (settings in
prereg/settings.json, power). Raw study files are fetched to the git-ignored data/cache/metabolights
(EMBL-EBI terms of use); only derived per-compound SDs are committed. Group rule, per study:

  MTBLS22, 41, 42     Col-0 rosette leaves, salt series: the 0 mM group (Source Name rg1)
  MTBLS2289           Genotype "Col-0 WT", Drought stress "water", shoots (its Col-0 tissue; OSD-522
                      is shoots too), one group per timepoint
  MTBLS4066           Sample group "primary_sample_wildtype" (leaves)
  MTBLS363            the Arabidopsis (Col-0) assay, rosette leaves, one group per timepoint
  excluded            MTBLS8105: its values carry commas whose decimal convention cannot be told
                      from the file. MTBLS7, 11, 12, 13 (single-leaf series): all 49-100 leaf samples
                      share one Source Name and the sheets do not record leaf position, so a group
                      would mix leaves of different age rather than replicate plants

The Golm salt-series studies (MTBLS22, 41, 42) report each compound as a percentage of its own maximum; that is a
per-compound linear rescaling, so log2 SDs are unchanged. Zeros are treated as missing. A group needs
at least 3 values. A compound's SD is the median over all groups (the pre-registered rule); per-study
medians are reported alongside.

Compounds are matched to KEGG by KEGG's own ChEBI conversion (conv/compound/chebi) when it lands on a
target ID, else by name against the target's KEGG synonyms with case, punctuation and L-/D- prefixes
ignored (GC-MS does not resolve enantiomers). Targets: the first study's 21 predicted compounds and
the T02 compounds with a directional darkness response. If several rows of one study map to the same
target, the most abundant (highest median) is used.
"""

import argparse
import hashlib
import json
import re
import urllib.parse
import urllib.request

import numpy as np
import pandas as pd
from scipy import optimize, stats

import style
from genesets import kegg
from paths import CACHE, FIGURES, PHOTORESP, PREREG, REF, TABLES

ML = CACHE / "metabolights"
FTP = "https://ftp.ebi.ac.uk/pub/databases/metabolights/studies/public/{sid}/{name}"
UA = {"User-Agent": "spaceflight-metabolism-natural-variation/1.0"}
STUDIES = ["MTBLS22", "MTBLS41", "MTBLS42", "MTBLS2289", "MTBLS4066", "MTBLS363"]
MIN_N = 3
REPLICATES_OSD522 = 6   # the first study's RNA-seq: 6 flight vs 6 ground shoots


def files(sid):
    s = next(ML.glob(f"{sid}/s_*.txt"))
    a = next(f for f in ML.glob(f"{sid}/a_*.txt") if "agave" not in f.name)
    m = next(f for f in ML.glob(f"{sid}/m_*.tsv") if "agave" not in f.name)
    return s, a, m


def fetch():
    """Download the study files listed, with their SHA-256, in data/reference/metabolights_sources.tsv."""
    src = pd.read_csv(REF / "metabolights_sources.tsv", sep="\t")
    for r in src.itertuples():
        out = ML / r.study / r.file
        if not out.exists():
            out.parent.mkdir(parents=True, exist_ok=True)
            req = urllib.request.Request(FTP.format(sid=r.study, name=urllib.parse.quote(r.file)), headers=UA)
            with urllib.request.urlopen(req, timeout=300) as resp:
                out.write_bytes(resp.read())
        digest = hashlib.sha256(out.read_bytes()).hexdigest()
        print(f"  {r.study}/{r.file} {'ok' if digest == r.sha256 else 'CHANGED since 2026-10-30'}")


def groups(sid):
    s, a, _ = files(sid)
    s = pd.read_csv(s, sep="\t", dtype=str)
    a = pd.read_csv(a, sep="\t", dtype=str)
    col = dict(zip(a["Sample Name"], a["MS Assay Name"]))
    f = lambda c: s[c].fillna("").str.strip() if c in s else pd.Series("", index=s.index)
    part = f("Characteristics[Organism part]").str.lower()
    if sid in ("MTBLS22", "MTBLS41", "MTBLS42"):
        g = s[s["Source Name"] == "rg1"].assign(group="0 mM")
    elif sid == "MTBLS2289":
        g = s[(f("Factor Value[Genotype]") == "Col-0 WT") & (f("Factor Value[Drought stress]") == "water")
              & (part == "shoot")].assign(group=lambda d: "day " + d["Factor Value[Timepoint]"])
    elif sid == "MTBLS4066":
        g = s[f("Factor Value[Sample group]") == "primary_sample_wildtype"].assign(group="wild type")
    elif sid == "MTBLS363":
        g = s[f("Characteristics[Organism]").str.contains("Arabidopsis")].assign(
            group=lambda d: d["Factor Value[Elapsed time]"])
    return {k: [col[n] for n in v["Sample Name"] if n in col] for k, v in g.groupby("group")}


def norm(name):
    n = str(name).lower().strip()
    n = re.sub(r"^(l|d|dl|\(s\)|\(r\))-", "", n)
    return re.sub(r"[^a-z0-9]", "", n)


def targets():
    t6 = pd.read_csv(PHOTORESP / "T06_predicted_compounds.tsv", sep="\t", comment="#")
    t2 = pd.read_csv(TABLES / "T02_light_dark_reference.tsv", sep="\t")
    keep = set(t6.kegg_id) | set(t2.loc[t2.status.isin(["robust", "single"]), "kegg_id"])
    syn, primary = {}, {}
    for ln in kegg("list/compound", "list_compound.tsv").splitlines():
        cid, names = ln.split("\t")
        if cid in keep:
            primary[cid] = names.split(";")[0].strip()
            for n in names.split(";"):
                syn.setdefault(norm(n), set()).add(cid)
    chebi = {}
    for ln in kegg("conv/compound/chebi", "conv_compound_chebi.tsv").splitlines():
        c, k = ln.split("\t")
        chebi.setdefault(c.replace("chebi:", "CHEBI:"), set()).add(k.replace("cpd:", ""))
    return keep, syn, primary, chebi, t6, t2


def match(maf, keep, syn, chebi):
    out = []
    for r in maf.itertuples():
        via_chebi = chebi.get(str(r.database_identifier).strip(), set()) & keep
        via_name = syn.get(norm(r.metabolite_identification), set())
        hit, how = (via_chebi, "ChEBI") if via_chebi else (via_name, "name") if via_name else (set(), "")
        if len(hit) == 1:
            out.append((r.Index, next(iter(hit)), how))
    return out


def sds():
    keep, syn, primary, chebi, *_ = targets()
    rows, matches = [], []
    for sid in STUDIES:
        _, _, m = files(sid)
        maf = pd.read_csv(m, sep="\t", dtype=str)
        grp = groups(sid)
        cols = sorted({c for v in grp.values() for c in v if c in maf.columns})
        vals = maf[cols].apply(lambda s: pd.to_numeric(s, errors="coerce")).where(lambda d: d > 0)
        picked = {}
        for idx, kid, how in match(maf, keep, syn, chebi):
            med = float(np.nanmedian(vals.loc[idx])) if vals.loc[idx].notna().any() else -1
            if kid not in picked or med > picked[kid][1]:
                picked[kid] = (idx, med, how, maf.loc[idx, "metabolite_identification"])
        for kid, (idx, _, how, name) in picked.items():
            matches.append({"study": sid, "kegg_id": kid, "compound": primary[kid], "maf_name": name, "via": how})
            for gname, gcols in grp.items():
                x = np.log2(vals.loc[idx, [c for c in gcols if c in vals.columns]].dropna().astype(float))
                if len(x) >= MIN_N:
                    rows.append({"study": sid, "group": gname, "kegg_id": kid, "n": len(x),
                                 "sd_log2": float(x.std(ddof=1))})
        print(f"  {sid}: {len(grp)} groups ({sum(map(len, grp.values()))} samples), {len(picked)} targets matched")
    return pd.DataFrame(rows), pd.DataFrame(matches), primary


def n_per_arm(effect, sd, alpha, power=0.8):
    """Smallest equal n per arm for a two-sided two-sample t-test (noncentral t)."""
    d = abs(effect) / sd

    def pw(n):
        df, nc = 2 * n - 2, d * np.sqrt(n / 2)
        t = stats.t.ppf(1 - alpha / 2, df)
        return stats.nct.sf(t, df, nc) + stats.nct.cdf(-t, df, nc)
    n = 2
    while pw(n) < power:
        n = n * 2 if n < 2 ** 30 else n + 1
        if n > 10 ** 7:
            return np.inf
    lo, hi = max(2, n // 2), n
    while lo < hi:
        mid = (lo + hi) // 2
        lo, hi = (mid + 1, hi) if pw(mid) < power else (lo, mid)
    return lo


def build():
    S = json.loads((PREREG / "settings.json").read_text())["power"]
    # fold_1p5 is log2(1.5) = 0.585; prereg/settings.json (frozen) calls it half_fold, a misleading name
    effects = {"Q60": float(np.log2(1 - 0.0237)), "Q1000": 0.13052, "fold_1p5": 0.585, "two_fold": 1.0}
    assert S["effects_log2"]["model_Q1000_T06"] == effects["Q1000"]
    assert S["effects_log2"]["half_fold"] == effects["fold_1p5"]
    g, m, primary = sds()
    m.to_csv(REF / "metabolights_kegg_matches.tsv", sep="\t", index=False)
    _, _, _, _, t6, t2 = targets()
    per_study = g.groupby(["kegg_id", "study"]).sd_log2.median().unstack()
    t = g.groupby("kegg_id").agg(sd_log2=("sd_log2", "median"), groups=("sd_log2", "size"),
                                 studies=("study", "nunique"), median_n=("n", "median"))
    t.insert(0, "compound", [primary[k] for k in t.index])
    t["in_first_study"] = t.index.isin(t6.kegg_id)
    t["darkness_status"] = t.index.map(dict(zip(t2.kegg_id, t2.status))).fillna("absent")
    t["sd_study_min"], t["sd_study_max"] = per_study.min(axis=1), per_study.max(axis=1)
    m_panel = len(t)
    for k, e in effects.items():
        t[f"n_{k}"] = [n_per_arm(e, s, S["alpha"], S["power"]) for s in t.sd_log2]
        t[f"n_{k}_bonf"] = [n_per_arm(e, s, S["alpha"] / m_panel, S["power"]) for s in t.sd_log2]
    d_min = t.sd_log2.apply(lambda s: optimize.brentq(
        lambda d: n_per_arm(d, s, S["alpha"], S["power"]) - REPLICATES_OSD522 - 0.5, 1e-4, 50))
    t["detectable_log2_at_n6"] = d_min
    t["detectable_log2_at_n6_bonf"] = t.sd_log2.apply(lambda s: optimize.brentq(
        lambda d: n_per_arm(d, s, S["alpha"] / m_panel, S["power"]) - REPLICATES_OSD522 - 0.5, 1e-4, 50))
    t = t.sort_values("sd_log2")
    with open(TABLES / "T10_gcms_power.tsv", "w") as fh:
        fh.write(f"# Samples per arm for 80% power, two-sided two-sample t, alpha {S['alpha']} "
                 f"(and Bonferroni over the {m_panel}-compound panel). Effects (log2): Q60 = log2(1 - 0.0237) "
                 "(first study's T13, rosette, OSD-522 light), Q1000 = 0.13052 (its T06), fold_1p5 = 0.585 "
                 "(1.5-fold; key half_fold in prereg/settings.json), two_fold 1.0.\n# SD: median within-group SD of log2 intensity over Col-0 wild-type control "
                 f"groups in {', '.join(STUDIES)} (MetaboLights; derived values only).\n"
                 "# detectable_log2_at_n6: smallest |log2 fold change| detectable with 6 per arm "
                 "(unadjusted; _bonf: Bonferroni over the panel).\n")
        t.to_csv(fh, sep="\t", float_format="%.4g")
    g.to_csv(REF / "metabolights_group_sds.tsv", sep="\t", index=False, float_format="%.4g")
    print(t[["compound", "sd_log2", "groups", "studies", "n_Q60", "n_Q1000", "n_Q1000_bonf",
             "n_fold_1p5", "n_fold_1p5_bonf", "detectable_log2_at_n6", "detectable_log2_at_n6_bonf"]].to_string())
    provenance()
    figure(t)


def provenance():
    """Per-study provenance of the variance groups, read from the ISA sample and assay files (T26)."""
    rule = {"MTBLS22": "Source Name rg1 (0 mM NaCl)", "MTBLS41": "Source Name rg1 (0 mM NaCl)",
            "MTBLS42": "Source Name rg1 (0 mM NaCl)",
            "MTBLS2289": "Genotype Col-0 WT, Drought stress water, shoot; one group per timepoint",
            "MTBLS4066": "Sample group primary_sample_wildtype",
            "MTBLS363": "Arabidopsis assay only (the study also profiles Agave americana); one group per elapsed time"}
    rows = []
    for sid in STUDIES:
        s, a, m = files(sid)
        s, a = pd.read_csv(s, sep="\t", dtype=str), pd.read_csv(a, sep="\t", dtype=str)
        grp = groups(sid)
        used = {c for v in grp.values() for c in v}
        names = set(a.loc[a["MS Assay Name"].isin(used), "Sample Name"])
        su = s[s["Sample Name"].isin(names)]
        au = a[a["MS Assay Name"].isin(used)]
        one = lambda d, c: "; ".join(sorted(d[c].dropna().unique())) if c in d else "not in ISA files"
        maf = pd.read_csv(m, sep="\t", dtype=str)
        cols = [c for c in maf.columns if c in set(a["MS Assay Name"])]       # every sample in the study's assay
        v = maf[cols].apply(lambda x: pd.to_numeric(x, errors="coerce"))
        mx = v.max(axis=1).dropna()
        other = sorted(set(s["Characteristics[Organism]"].dropna()) - {"Arabidopsis thaliana", "blank", "reference compound"})
        n = [len(x) for x in grp.values()]
        rows.append({"study": sid, "organism": one(su, "Characteristics[Organism]"),
                     "variant": one(su, "Characteristics[Variant]"), "organism_part": one(su, "Characteristics[Organism part]"),
                     "group_rule": rule[sid], "groups": len(grp), "samples": sum(n),
                     "n_per_group": f"{min(n)}-{max(n)}" if min(n) != max(n) else str(n[0]),
                     "chromatography": one(au, "Parameter Value[Chromatography Instrument]"),
                     "instrument": one(au, "Parameter Value[Instrument]"),
                     "mass_analyzer": one(au, "Parameter Value[Mass analyzer]"),
                     "compounds_whose_study_max_is_100": f"{int(np.isclose(mx, 100).sum())} of {len(mx)}",
                     "other_organisms_in_study": "; ".join(other) or "none"})
    t = pd.DataFrame(rows)
    t.to_csv(TABLES / "T26_metabolights_provenance.tsv", sep="\t", index=False)
    print(t.to_string(index=False))


def figure(t):
    style.apply()
    import matplotlib.pyplot as plt

    t = t.sort_values("sd_log2", ascending=False)
    fig, ax = plt.subplots(figsize=(7.0, 5.0))
    y = np.arange(len(t))
    marks = {"n_Q60": (style.PURPLE, "v", "model, OSD-522 light (Q 60): 0.035 log2"),
             "n_Q1000": (style.VERMILION, "o", "model, Q 1000 (first study's T06): 0.13 log2"),
             "n_fold_1p5": (style.BLUE, "s", "1.5-fold (0.585 log2)"),
             "n_two_fold": (style.GREEN, "D", "2-fold (1 log2)")}
    for col, (c, mk, lab) in marks.items():
        ax.scatter(t[col], y, color=c, marker=mk, s=20, label=lab, zorder=3)
        ax.scatter(t[f"{col}_bonf"], y, facecolor="none", edgecolor=c, marker=mk, s=20, lw=0.8, zorder=3)
        ax.hlines(y, t[col], t[f"{col}_bonf"], color=c, lw=0.5, alpha=0.5, zorder=2)
    ax.axvline(REPLICATES_OSD522, color=style.INK, lw=1, ls="--", label="6 per arm (OSD-522 RNA-seq)")
    ax.set_xscale("log")
    ax.set_yticks(y)
    ax.set_yticklabels([f"{c}{' *' if f else ''}  (SD {s:.2f})" for c, f, s in
                        zip(t.compound, t.in_first_study, t.sd_log2)], fontsize=6.5)
    ax.set_xlabel(f"Plants per arm for 80% power (two-sided t; filled: alpha 0.05; open: Bonferroni over {len(t)} compounds)")
    ax.grid(axis="y", visible=False)
    ax.legend(loc="upper center", bbox_to_anchor=(0.45, -0.1), ncol=2, fontsize=6.3)
    style.panel(ax, "", "Targeted GC-MS panel: replicates needed per compound")
    fig.text(0.01, -0.1, "* predicted by the first study. SD: median within-group SD of log2 intensity in Col-0 "
             "control leaves and shoots, MetaboLights GC-MS (T10).", fontsize=6.4, color=style.GREY)
    style.save(fig, "fig05_gcms_power", str(FIGURES))


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("step", choices=["fetch", "build"])
    {"fetch": fetch, "build": build}[ap.parse_args().step]()
