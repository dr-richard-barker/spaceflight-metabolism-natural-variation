"""The darkness reference by stratum, its threshold sensitivity, and its KEGG mapping (T14-T16, figS1).

Revision 1, post hoc: the settings are in prereg/revision1/settings.json (R1), committed before this
script ran. Nothing here changes the pre-registered reference (T02) or the T03 verdicts.

- T14: per compound, each source's effect, n, support and qualification, and the consensus under
  all three sources (= T02), dark only (Zhu, Caldana) and dark+heat (Wu: darkness at 32 C). It shows
  which calls rest on the 32 C source alone, and the model verdict under the dark-only stratum.
- T15: the consensus recomputed over accession agreement x Caldana window, for the all and the
  dark-only strata: category counts, the key amino- and organic-acid calls, and the T03 verdicts.
- T16: how the 112 measured names reached KEGG IDs (candidates, rule matches, reviewed overrides,
  exclusions) and the mapped fraction by Wu metabolite class.
"""

import importlib.util
import json

import numpy as np
import pandas as pd

import style
from paths import ATLAS, FIGURES, PREREG, REF, ROOT, TABLES

S = json.loads((PREREG / "revision1" / "settings.json").read_text())["R1_reference"]


def load(name):
    spec = importlib.util.spec_from_file_location(name.replace(".py", "").replace("-", "_"), ROOT / "scripts" / name)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


R10 = load("10_light_dark_reference.py")
R12 = load("12_model_vs_reference.py")
fixed = S["fixed"]
assert (R10.MIN_N, R10.ROBUST_FRAC, R10.P_CUT) == (fixed["min_accessions"], fixed["robust_single_source_fraction"],
                                                 fixed["caldana_p"]), "revision settings differ from script 10"
PRE = (R10.QUALIFY_FRAC, R10.LATE_MIN)
SOURCES = ("zhu", "caldana", "wu_stress")
SOURCE_LABEL = {"zhu": "Zhu (6 d dark)", "caldana": "Caldana (dark vs light)", "wu_stress": "Wu (dark, 32 °C)"}
STRATUM_LABEL = {"all": "All three (T02)", "dark_only": "Dark only", "dark_heat": "Dark+heat only"}
KEY = {"C00123": "Leu", "C00407": "Ile", "C00183": "Val", "C00149": "Malate", "C00026": "2-OG",
       "C00025": "Glu", "C00064": "Gln"}
assert set(KEY) == {k for v in S["key_calls"].values() for k in v}


def call(ref, kid):
    if kid not in set(ref.kegg_id):
        return "absent"
    r = ref.set_index("kegg_id").loc[kid]
    return r.status + ({1: " up", -1: " down"}.get(int(r.direction), ""))


def strata(src):
    return {name: R10.consensus(src, tuple(use)) for name, use in S["strata"].items()}


# ---------- T14 ----------

def t14(src, st):
    t02 = pd.read_csv(TABLES / "T02_light_dark_reference.tsv", sep="\t")
    a = st["all"]
    assert list(a.kegg_id) == list(t02.kegg_id) and (a.status.values == t02.status.values).all() \
        and (a.direction.values == t02.direction.values).all(), "the all-source stratum must reproduce T02"
    s = src[src.kegg_id.fillna("") != ""].set_index(["kegg_id", "source"])
    model = {name: R12.score(ref).query("in_model").set_index("kegg_id") for name, ref in st.items()}
    rows = []
    for r in a.itertuples():
        row = {"kegg_id": r.kegg_id, "name": r.name.strip()}
        for so in SOURCES:
            if (r.kegg_id, so) in s.index:
                x = s.loc[(r.kegg_id, so)]
                row |= {f"{so}_n": int(x.n), f"{so}_effect": round(float(x.effect), 3),
                        f"{so}_support": float(x.p) if so == "caldana" else round(float(x.frac_same), 4),
                        f"{so}_qualifies": bool(x.qualifies)}
            else:
                row |= {f"{so}_n": np.nan, f"{so}_effect": np.nan, f"{so}_support": np.nan, f"{so}_qualifies": np.nan}
        for name, ref in st.items():
            hit = ref.set_index("kegg_id").loc[r.kegg_id] if r.kegg_id in set(ref.kegg_id) else None
            row[f"status_{name}"] = hit.status if hit is not None else "absent"
            row[f"direction_{name}"] = int(hit.direction) if hit is not None else 0
        directional = {n: row[f"status_{n}"] in ("robust", "single") for n in st}
        row["dark_heat_only"] = directional["all"] and not directional["dark_only"]
        row["flips_dark_only"] = (directional["all"] and directional["dark_only"]
                                  and row["direction_all"] != row["direction_dark_only"])
        row["in_model"] = r.kegg_id in model["all"].index
        row["verdict"] = model["all"].verdict.get(r.kegg_id, "")
        row["verdict_dark_only"] = model["dark_only"].verdict.get(r.kegg_id, "")
        rows.append(row)
    t = pd.DataFrame(rows)
    t.to_csv(TABLES / "T14_reference_strata.tsv", sep="\t", index=False, float_format="%.4g")
    d = t[t.status_all.isin(["robust", "single"])]
    print(f"  T14: {len(t)} compounds; {len(d)} directional in T02, {int(d.dark_heat_only.sum())} of them rest on "
          f"the 32 C source alone, {int(t.flips_dark_only.sum())} flip dark-only; dark-only status "
          f"{t.status_dark_only.value_counts().to_dict()}")
    m = t[t.in_model]
    print(m[["name", "status_all", "status_dark_only", "verdict", "verdict_dark_only"]].to_string(index=False))
    return t


# ---------- T15 ----------

def t15():
    base = {cm: R10.sources(late_min=cm) for cm in S["sensitivity_caldana_min_minutes"]}
    rows = []
    for qf in S["sensitivity_agreement_fraction"]:
        for cm, b in base.items():
            src = b.copy()
            src["qualifies"] = [R10.qualifies(r, qf) for r in src.itertuples()]
            st = strata(src)
            for name in ("all", "dark_only"):
                ref = st[name]
                sc = R12.score(ref)
                m = sc[sc.in_model]
                calls = {lab: call(ref, kid) for kid, lab in KEY.items()}
                c = ref.status.value_counts()
                rows.append({"stratum": name, "agreement": qf, "caldana_min_minutes": cm,
                             "preregistered": (qf, cm) == PRE,
                             **{k: int(c.get(k, 0)) for k in ("robust", "single", "conflicting", "none")},
                             "up": int((ref.direction > 0).sum()), "down": int((ref.direction < 0).sum()),
                             **calls,
                             "bcaa_up": all(calls[k].endswith(" up") for k in ("Leu", "Ile", "Val")),
                             "malate_2og_down": all(calls[k].endswith(" down") for k in ("Malate", "2-OG")),
                             "model_agrees": int((m.verdict == "agrees").sum()),
                             "model_disagrees": int((m.verdict == "disagrees").sum()),
                             "model_untested": int((m.verdict == "untested").sum()),
                             "robust_not_predicted": int(((~sc.in_model) & (sc.status == "robust")).sum())})
    t = pd.DataFrame(rows)
    t.to_csv(TABLES / "T15_reference_sensitivity.tsv", sep="\t", index=False)
    for name, g in t.groupby("stratum", sort=False):
        print(f"  T15 {name}: robust {g.robust.min()}-{g.robust.max()}, BCAA up in {int(g.bcaa_up.sum())}/{len(g)}, "
              f"malate+2-OG down in {int(g.malate_2og_down.sum())}/{len(g)}, Glu {sorted(set(g.Glu))}, "
              f"agrees {g.model_agrees.min()}-{g.model_agrees.max()}, disagrees {g.model_disagrees.min()}-"
              f"{g.model_disagrees.max()}")
    return t


# ---------- T16 ----------

def wu_class():
    ident = pd.read_csv(ATLAS / "wu_identified_metabolites.csv")
    fix = {"Glocoside": "Glucoside", "Amino acid deri": "Amino acid derivative"}
    c = ident.drop_duplicates("metabolite_name").set_index("metabolite_name").metabolite_class
    return c.replace(fix)


def t16(src):
    km = pd.read_csv(REF / "kegg_map.tsv", sep="\t", dtype=str).fillna("")
    ov = pd.read_csv(REF / "kegg_overrides.tsv", sep="\t", dtype=str, comment="#").fillna("")
    final = R10.load_kegg()
    rows = []

    def add(section, item, n, of="", note=""):
        rows.append({"section": section, "item": item, "n": n, "of": of, "note": note})

    add("names", "measured names sent to KEGG (Zhu and Wu)", len(km))
    for so, g in km.groupby("source"):
        add("names", f"from {so}", len(g), len(km))
    ncand = km.candidates.map(lambda c: len([x for x in c.split(";") if x]))
    for lab, mask in (("0", ncand == 0), ("1", ncand == 1), (">1", ncand > 1)):
        add("KEGG find candidates", f"{lab} candidate(s)", int(mask.sum()), len(km),
            "rest.kegg.jp/find/compound on the name and on the name without a trailing bracket")
    mt = km.match.str.replace(r" \(\d+ hits\)", " (several hits, lowest ID kept)", regex=True)
    for lab, n in mt.value_counts().items():
        add("rule match (exact synonym)", lab, int(n), len(km))
    add("rule match (exact synonym)", "mapped by the rule", int((km.kegg_id != "").sum()), len(km))
    in_map = set(km.name)
    rule = dict(zip(km.name, km.kegg_id))
    kind = []
    for r in ov.itertuples():
        if r.name not in in_map:
            kind.append("Caldana peak only: " + ("excluded" if not r.kegg_id else "Wu's given ID corrected"))
        elif not r.kegg_id:
            kind.append("exclusion")
        elif not rule[r.name]:
            kind.append("rescued (rule found no exact synonym)")
        elif rule[r.name] != r.kegg_id:
            kind.append("corrected (rule picked another entry)")
        else:
            kind.append("confirmed")
    ov["kind"] = kind
    add("reviewed overrides", "rows in data/reference/kegg_overrides.tsv", len(ov))
    for lab, n in ov.kind.value_counts().items():
        add("reviewed overrides", lab, int(n), len(ov))
    mapped = {n for n, k in final.items() if k}
    add("after overrides", "measured names with a KEGG ID", len(mapped), len(km))
    for so, g in km.groupby("source"):
        add("after overrides", f"from {so}", int(g.name.isin(mapped).sum()), len(g))
    cal = src[src.source == "caldana"]
    add("after overrides", "Caldana peaks (KEGG IDs given by Wu et al.)", int((cal.kegg_id.fillna("") != "").sum()), len(cal))
    add("after overrides", "Caldana peaks whose given ID an override changed or removed",
        int((cal.kegg_id.fillna("") != cal.kegg_given).sum()), len(cal), "overrides apply by name to every source")
    cls = wu_class()
    wu = km[km.source == "wu_stress"].assign(cls=lambda d: d.name.map(cls).fillna("not in Wu Table 2"))
    for c, g in wu.groupby("cls"):
        add("Wu class (after overrides)", c, int(g.name.isin(mapped).sum()), len(g))
    t = pd.DataFrame(rows)
    t.to_csv(TABLES / "T16_kegg_mapping_report.tsv", sep="\t", index=False)
    print(t.to_string(index=False))
    return t


# ---------- figure ----------

def figure(t14, t15):
    style.apply()
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D

    d = t14[t14.status_all.isin(["robust", "single", "conflicting"]) | t14.in_model].copy()
    d["z"] = pd.read_csv(TABLES / "T02_light_dark_reference.tsv", sep="\t").set_index("kegg_id").reference_z.reindex(d.kegg_id).values
    d = d.sort_values(["direction_all", "z"], ascending=[False, False]).reset_index(drop=True)
    cols = [*SOURCES, "all", "dark_only"]
    head = ["Zhu\n6 d dark", "Caldana\ndark vs\nlight", "Wu\ndark,\n32 °C", "Cons.:\nall\nthree", "Cons.:\ndark\nonly"]

    fig = plt.figure(figsize=(7.2, 9.4))
    gs = fig.add_gridspec(2, 2, width_ratios=[1, 1.25], height_ratios=[1, 0.95], wspace=0.6, hspace=0.62)
    a = fig.add_subplot(gs[:, 0])
    b = fig.add_subplot(gs[0, 1])
    c = fig.add_subplot(gs[1, 1])

    y = np.arange(len(d))[::-1]
    for j, col in enumerate(cols):
        for yi, r in zip(y, d.itertuples()):
            if col in SOURCES:
                eff, q = getattr(r, f"{col}_effect"), getattr(r, f"{col}_qualifies")
                if pd.isna(eff):
                    continue
                up = eff > 0
                filled = bool(q)
            else:
                st, dr = getattr(r, f"status_{col}"), getattr(r, f"direction_{col}")
                if st == "absent":
                    continue
                if dr == 0:
                    a.plot(j, yi, marker="o", ms=3.2, mfc="none", mec=style.GREY, mew=0.7)
                    continue
                up, filled = dr > 0, st == "robust"
            colr = style.BLUE if up else style.VERMILION
            a.plot(j, yi, marker="^" if up else "v", ms=4.2, color=colr, mfc=colr if filled else "none", mew=0.8)
    a.axvline(2.5, color=style.GREY, lw=0.6)
    a.set_xticks(range(len(cols)))
    a.set_xticklabels(head, fontsize=5.6, linespacing=1.0)
    a.set_yticks(y)
    labels = [f"{r.name}{' *' if r.in_model else ''}" for r in d.itertuples()]
    a.set_yticklabels(labels, fontsize=5.0)
    for tl, r in zip(a.get_yticklabels(), d.itertuples()):
        if r.dark_heat_only:
            tl.set_color(style.VERMILION)
        if r.in_model:
            tl.set_fontweight("bold")
    a.set_ylim(-0.8, len(d) - 0.2)
    a.set_xlim(-0.6, len(cols) - 0.4)
    a.grid(axis="x", visible=False)
    a.tick_params(axis="y", length=0, pad=1.5)
    style.panel(a, "a", "Per-source calls")
    a.legend(handles=[Line2D([], [], marker="^", ls="", color=style.BLUE, ms=4.2, label="up"),
                      Line2D([], [], marker="v", ls="", color=style.VERMILION, ms=4.2, label="down"),
                      Line2D([], [], marker="^", ls="", color=style.GREY, mfc="none", ms=4.2,
                             label="hollow: source does not qualify,\nor consensus single-source"),
                      Line2D([], [], marker="o", ls="", color=style.GREY, mfc="none", ms=3.2, label="no direction")],
             loc="upper center", bbox_to_anchor=(0.4, -0.05), fontsize=5.8, ncol=2, handletextpad=0.3,
             columnspacing=0.8)

    # b: robust + single counts over the grid
    win = list(zip(S["sensitivity_caldana_min_minutes"], (style.SKY, style.GREEN, style.BLACK)))
    for name, ls in (("all", "-"), ("dark_only", "--")):
        for cm, colr in win:
            g = t15[(t15.stratum == name) & (t15.caldana_min_minutes == cm)]
            b.plot(g.agreement * 100, g.robust + g.single, ls=ls, color=colr, marker="o", ms=2.6, lw=1.1)
    pre = t15[t15.preregistered & (t15.stratum == "all")].iloc[0]
    b.plot(pre.agreement * 100, pre.robust + pre.single, marker="o", ms=8, mfc="none", mec=style.INK, mew=0.9)
    b.set_xlabel("Accession agreement threshold (%)")
    b.set_ylabel("Directional compounds\n(robust + single)")
    b.set_xticks([v * 100 for v in S["sensitivity_agreement_fraction"]])
    b.legend(handles=[Line2D([], [], color=colr, lw=1.1, label=f"Caldana ≥ {cm} min") for cm, colr in win]
             + [Line2D([], [], color=style.GREY, ls="-", label="all three sources"),
                Line2D([], [], color=style.GREY, ls="--", label="dark only"),
                Line2D([], [], marker="o", ls="", ms=7, mfc="none", mec=style.INK, label="pre-registered")],
             fontsize=5.8, loc="upper center", bbox_to_anchor=(0.5, -0.2), ncol=2)
    style.panel(b, "b", "Threshold sensitivity")

    # c: key calls over the grid
    order = list(KEY.values())
    wins = sorted(S["sensitivity_caldana_min_minutes"], reverse=True)
    grid = pd.concat([t15[(t15.stratum == name) & (t15.caldana_min_minutes == cm)].sort_values("agreement")
                      for name in ("all", "dark_only") for cm in wins])
    shade = {style.BLUE: "#BFD9EC", style.VERMILION: "#F5D0B5"}
    for i, k in enumerate(order):
        yy = len(order) - 1 - i
        for j, v in enumerate(grid[k]):
            colr = style.BLUE if v.endswith(" up") else style.VERMILION if v.endswith(" down") else None
            if colr is None:
                fc = style.GREY_LIGHT if v != "absent" else "white"
            else:
                fc = colr if v.startswith("robust") else shade[colr]
            c.add_patch(plt.Rectangle((j - 0.42, yy - 0.4), 0.84, 0.8, fc=fc, ec="none"))
            if v == "absent":
                c.text(j, yy, "–", ha="center", va="center", fontsize=5.5, color=style.GREY)
    n = len(grid) // 2
    c.axvline(n - 0.5, color=style.INK, lw=0.9)
    c.set_xlim(-0.6, len(grid) - 0.4)
    c.set_ylim(-0.6, len(order) + 0.5)
    c.set_yticks(range(len(order)))
    c.set_yticklabels(order[::-1])
    c.set_xticks(range(len(grid)))
    c.set_xticklabels([f"{int(r.agreement * 100)}" for r in grid.itertuples()], fontsize=4.8, rotation=90)
    c.tick_params(length=0, pad=1)
    per = len(grid) // (2 * len(wins))
    for blk in range(2 * len(wins)):
        x0 = blk * per - 0.5
        if blk and blk != len(wins):
            c.axvline(x0, color=style.GREY, lw=0.4, ymax=len(order) / (len(order) + 1.1))
        c.text(x0 + per / 2, len(order) - 0.35, f"≥{wins[blk % len(wins)]}", ha="center", va="bottom", fontsize=5.4)
    for k0, lab in ((0, "all three sources"), (n, "dark only")):
        c.text(k0 + n / 2 - 0.5, len(order) + 0.3, lab, ha="center", va="bottom", fontsize=6.6, fontweight="bold")
    for j, r in enumerate(grid.itertuples()):
        if r.preregistered:
            c.add_patch(plt.Rectangle((j - 0.5, -0.5), 1.0, len(order), fc="none", ec=style.INK, lw=0.9))
    c.set_xlabel("Agreement (%) per Caldana window (min)")
    c.grid(False)
    c.spines[["left", "bottom"]].set_visible(False)
    c.set_title("c  Key calls over the grid", loc="left", fontweight="bold", fontsize=9.5, pad=18)
    c.legend(handles=[plt.Rectangle((0, 0), 1, 1, fc=style.BLUE, label="robust up"),
                      plt.Rectangle((0, 0), 1, 1, fc=shade[style.BLUE], label="single up"),
                      plt.Rectangle((0, 0), 1, 1, fc=style.VERMILION, label="robust down"),
                      plt.Rectangle((0, 0), 1, 1, fc=shade[style.VERMILION], label="single down"),
                      plt.Rectangle((0, 0), 1, 1, fc=style.GREY_LIGHT, label="no direction"),
                      plt.Rectangle((0, 0), 1, 1, fc="white", ec=style.GREY, lw=0.4, label="absent (–)")],
             loc="upper center", bbox_to_anchor=(0.5, -0.17), ncol=3, fontsize=5.8)
    fig.text(0.005, -0.012, "* model compound.  Orange names: directional only because of the 32 °C dark+heat source (Wu).  "
             "Outlined column in c: the pre-registered setting (75 %, ≥ 360 min).", fontsize=6.0, color=style.GREY)
    style.save(fig, "figS1_reference_strata", FIGURES)


if __name__ == "__main__":
    src = R10.sources()
    st = strata(src)
    a = t14(src, st)
    b = t15()
    t16(src)
    figure(a, b)
