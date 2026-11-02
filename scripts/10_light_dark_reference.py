"""Measured light/dark metabolite reference: what darkness does to Arabidopsis leaf metabolism.

Three independent ground sources, each a light-vs-dark contrast measured on real plants:

  zhu       Zhu et al. 2024: per-accession BLUPs at 0 d and 6 d extended darkness (259 accessions,
            leaves, LC-MS). Response = 6 d - 0 d per accession.
  wu_stress Wu et al. 2018: control vs a stress treatment that the supplement calls the 32 C darkness
            stress condition (about 300 accessions, leaves, LC-MS). Response = stress - control per
            accession, one ion mode per metabolite and unpaired metabolites left out (wu_modes.py).
  caldana   Wu et al. 2018 Supplementary Table 7: LC-MS of the Caldana et al. 2011 time course (one
            genotype; plants grown at 21 C, 150 uE, then kept or moved to 21 C darkness). Response =
            mean of 21 C dark minus 21 C light (labelled 20-D and 20-L in the file) over the
            timepoints >= 360 min, 3 replicates each; Welch t-test. Only peaks Wu named with a KEGG
            ID (Supplementary Table 8) are used. The genotype is not stated in the Wu supplement or
            the Caldana abstract (the paper is not open access): TODO, so it is not called Col-0.

Rule for a source to count (fixed before the results were looked at):
  zhu, wu_stress: |median response| > 0 and at least 75 % of accessions move the same way, n >= 20.
  caldana:        Welch p < 0.05.
A compound's consensus direction is the shared sign of its qualifying sources; it is "conflicting"
if qualifying sources disagree, "robust" if at least two agree or one has >= 90 % of accessions.
The reference value is the mean of the qualifying sources' standardised effects (response divided
by that source's natural spread: the SD across accessions at 0 d / control, or the light-arm SD).

KEGG compound IDs come from the KEGG REST API (`kegg` step), never from memory: a name maps only
when it equals one of the entry's KEGG synonyms after case and punctuation are ignored, or the same
with an L-/D- prefix. Wu's own IDs for the time-course peaks are kept as given. Every candidate and
the chosen ID go to data/reference/kegg_map.tsv; data/reference/kegg_overrides.tsv records reviewed
corrections. Steps:
  python3 scripts/10_light_dark_reference.py fetch      # Wu supplement into data/cache (git-ignored)
  python3 scripts/10_light_dark_reference.py kegg       # KEGG lookups -> data/reference/kegg_map.tsv
  python3 scripts/10_light_dark_reference.py build      # T01, T02
"""

import argparse
import hashlib
import re
import time
import urllib.error
import urllib.parse
import urllib.request

import numpy as np
import pandas as pd
from scipy.stats import ttest_ind

from paths import ATLAS, CACHE, REF, TABLES, ATLAS_REPO
from wu_modes import stress_minus_control, one_mode

WU_CACHE = CACHE / "wu2018"
WU_URL = "https://ars.els-cdn.com/content/image/1-s2.0-S1674205217302423-{f}"
# SHA-256 prefixes from arabidopsis-atlas data/raw/SUPPLEMENTARY_SOURCES.md
WU_FILES = {"mmc7.xlsx": "ab8526766adfdc6b", "mmc8.xlsx": "7c1080b45b9f2b6c"}
KEGG_CACHE = CACHE / "kegg_find"
LATE_MIN = 360
QUALIFY_FRAC, ROBUST_FRAC, MIN_N, P_CUT = 0.75, 0.90, 20, 0.05
UA = {"User-Agent": "spaceflight-metabolism-natural-variation/0.1"}


def key(name):
    return re.sub(r"[^a-z0-9]", "", str(name).lower())


# ---------- fetch ----------

def fetch():
    WU_CACHE.mkdir(parents=True, exist_ok=True)
    for f, expect in WU_FILES.items():
        dest = WU_CACHE / f
        if not dest.exists():
            req = urllib.request.Request(WU_URL.format(f=f), headers={"User-Agent": "Mozilla/5.0"})
            try:
                with urllib.request.urlopen(req, timeout=120) as r:
                    dest.write_bytes(r.read())
            except Exception as e:  # the Elsevier CDN refuses some clients; fall back to the atlas copy
                local = ATLAS_REPO / "data" / "raw" / "wu_2017_supplementary" / f
                print(f"  CDN fetch failed ({e}); copying {local}")
                dest.write_bytes(local.read_bytes())
        got = hashlib.sha256(dest.read_bytes()).hexdigest()[:16]
        assert got == expect, f"{f}: SHA-256 {got} != {expect}"
        print(f"  {f} ok ({got})")


# ---------- sources ----------

def zhu_response():
    z = pd.read_csv(ATLAS / "darkness_metabolome.csv")
    w = z.pivot_table(index=["metabolite_id", "accession"], columns="timepoint", values="intensity")
    possible = w[["0d darkness", "6d darkness"]].notna().any(axis=1).groupby(level=0).sum()
    d = (w["6d darkness"] - w["0d darkness"]).dropna()
    base = w["0d darkness"].dropna().groupby(level=0).std()
    names = z.drop_duplicates("metabolite_id").set_index("metabolite_id").metabolite_name
    rows = []
    for mid, s in d.groupby(level=0):
        med = s.median()
        rows.append({"source": "zhu", "source_id": mid, "name": names[mid], "n": len(s), "effect": med,
                     "frac_same": float((np.sign(s) == np.sign(med)).mean()) if med != 0 else 0.0,
                     "spread": base.get(mid, np.nan), "p": np.nan,
                     "n_possible": int(possible[mid]), "missing_frac": round(1 - len(s) / possible[mid], 4)})
    return pd.DataFrame(rows)


def wu_stress_response():
    wu = pd.read_csv(ATLAS / "wu_leaf_identified_long.csv")
    ident = pd.read_csv(ATLAS / "wu_identified_metabolites.csv")
    sm = stress_minus_control(wu, ident)
    ctrl_sd = one_mode(wu).query("condition == 'control'").groupby("metabolite_name").intensity.std()
    rows = []
    for name, g in sm.groupby(level=0):
        d = g.delta.dropna()
        if d.empty:
            continue
        med = d.median()
        possible = int((g.control.notna() | g.stress.notna()).sum())
        rows.append({"source": "wu_stress", "source_id": name, "name": name, "n": len(d), "effect": med,
                     "frac_same": float((np.sign(d) == np.sign(med)).mean()) if med != 0 else 0.0,
                     "spread": ctrl_sd.get(name, np.nan), "p": np.nan,
                     "n_possible": possible, "missing_frac": round(1 - len(d) / possible, 4)})
    return pd.DataFrame(rows)


def caldana_response(late_min=LATE_MIN):
    m8 = pd.read_excel(WU_CACHE / "mmc8.xlsx", header=1).rename(columns={"Peak ID": "peak"})
    m8["peak"] = m8.peak.astype(str).str.replace(r"[\sX]", "", regex=True)
    m8 = m8[m8["KEGG ID"].astype(str).str.match(r"^C\d{5}$")].drop_duplicates("peak")
    rows = []
    for sheet in pd.ExcelFile(WU_CACHE / "mmc7.xlsx").sheet_names:
        t = pd.read_excel(WU_CACHE / "mmc7.xlsx", sheet_name=sheet, header=1)
        t["peak"] = t.Peak_ID.astype(str).str.replace(r"[\sX]", "", regex=True)
        t = t[t.peak.isin(m8.peak)]

        def cols(cond):
            out = []
            for c in t.columns:
                m = re.match(rf"^{cond}-(\d+)-\d+$", str(c))
                if m and int(m.group(1)) >= late_min:
                    out.append(c)
            return out
        dark, light = cols("20-D"), cols("20-L")
        assert dark and light, f"{sheet}: no 20-D/20-L columns at >= {late_min} min"
        info_by_peak = m8.set_index("peak")
        dark_vals = t[dark].apply(pd.to_numeric, errors="coerce")
        light_vals = t[light].apply(pd.to_numeric, errors="coerce")
        for i, peak in enumerate(t.peak):
            a, b = dark_vals.iloc[i].dropna(), light_vals.iloc[i].dropna()
            if len(a) < 3 or len(b) < 3:
                continue
            info = info_by_peak.loc[peak]
            rows.append({"source": "caldana", "source_id": peak, "name": info.Metabolite,
                         "kegg_given": info["KEGG ID"], "n": len(a) + len(b), "effect": a.mean() - b.mean(),
                         "frac_same": np.nan, "spread": b.std(), "p": ttest_ind(a, b, equal_var=False).pvalue,
                         "n_possible": len(dark) + len(light),
                         "missing_frac": round(1 - (len(a) + len(b)) / (len(dark) + len(light)), 4),
                         "mean_intensity": pd.concat([a, b]).mean()})
    df = pd.DataFrame(rows)
    # one peak per metabolite: the most intense (fixed rule)
    return df.sort_values("mean_intensity", ascending=False).drop_duplicates("kegg_given")


# ---------- KEGG ----------

def kegg_find(name):
    KEGG_CACHE.mkdir(parents=True, exist_ok=True)
    f = KEGG_CACHE / (hashlib.sha1(name.encode()).hexdigest() + ".txt")
    if not f.exists():
        url = "https://rest.kegg.jp/find/compound/" + urllib.parse.quote(name, safe="")
        req = urllib.request.Request(url, headers=UA)
        try:
            with urllib.request.urlopen(req, timeout=60) as r:
                f.write_text(r.read().decode())
        except urllib.error.HTTPError as e:  # KEGG answers 400 to some long or bracketed names
            if e.code != 400:
                raise
            f.write_text("")
        time.sleep(0.4)
    out = []
    for line in f.read_text().splitlines():
        if "\t" in line:
            cid, syn = line.split("\t", 1)
            out.append((cid.replace("cpd:", ""), [s.strip() for s in syn.split(";")]))
    return out


def kegg_map():
    names = pd.concat([zhu_response()[["source", "name"]], wu_stress_response()[["source", "name"]]])
    names = names.drop_duplicates("name")
    rows = []
    for r in names.itertuples():
        base = re.sub(r"\s*\(.*?\)\s*$", "", r.name).strip()   # "Pantothenic acid (Vitamine B5)"
        chosen, match, cands = None, "none", []
        for q in dict.fromkeys([r.name, base]):
            hits = kegg_find(q)
            cands += [c for c, _ in hits]
            for variant, label in ((q, "exact"), ("L-" + q, "L-prefix"), ("D-" + q, "D-prefix")):
                exact = [c for c, syn in hits if key(variant) in {key(s) for s in syn}]
                if not exact and label != "exact":
                    exact = [c for c, syn in kegg_find(variant) if key(variant) in {key(s) for s in syn}]
                if exact:
                    chosen, match = sorted(exact)[0], label + (f" ({len(exact)} hits)" if len(exact) > 1 else "")
                    break
            if chosen:
                break
        rows.append({"name": r.name, "source": r.source, "kegg_id": chosen or "", "match": match,
                     "candidates": ";".join(dict.fromkeys(cands))[:200]})
    df = pd.DataFrame(rows)
    REF.mkdir(parents=True, exist_ok=True)
    df.to_csv(REF / "kegg_map.tsv", sep="\t", index=False)
    print(f"  kegg_map.tsv: {len(df)} names, {int((df.kegg_id != '').sum())} mapped")


def load_kegg():
    m = pd.read_csv(REF / "kegg_map.tsv", sep="\t", dtype=str).fillna("")
    ov = REF / "kegg_overrides.tsv"
    if ov.exists():
        o = pd.read_csv(ov, sep="\t", dtype=str, comment="#").fillna("")
        for r in o.itertuples():
            m.loc[m.name == r.name, ["kegg_id", "match"]] = [r.kegg_id, f"override: {r.reason}"]
    return dict(zip(m.name, m.kegg_id))


def kegg_overrides():
    """Reviewed corrections by name; they also apply to the IDs Wu gave the time-course peaks."""
    ov = REF / "kegg_overrides.tsv"
    o = pd.read_csv(ov, sep="\t", dtype=str, comment="#").fillna("")
    return {n.strip(): k for n, k in zip(o.name, o.kegg_id)}


# ---------- build ----------

def qualifies(r, qualify_frac=QUALIFY_FRAC):
    if r.source == "caldana":
        return bool(r.p < P_CUT)
    return bool(r.n >= MIN_N and r.effect != 0 and r.frac_same >= qualify_frac)


def sources(qualify_frac=QUALIFY_FRAC, late_min=LATE_MIN):
    """Every source response with its KEGG ID, qualification and standardised effect.
    The defaults are the pre-registered settings; scripts/23_reference_strata.py varies them."""
    kmap = load_kegg()
    src = pd.concat([zhu_response(), wu_stress_response(), caldana_response(late_min)], ignore_index=True)
    ovr = kegg_overrides()
    src["kegg_id"] = [ovr.get(r.name.strip(), r.kegg_given if r.source == "caldana" else kmap.get(r.name, ""))
                      for r in src.itertuples()]
    src["qualifies"] = [qualifies(r, qualify_frac) for r in src.itertuples()]
    src["z"] = src.effect / src.spread
    dup = src[src.kegg_id != ""].duplicated(["source", "kegg_id"], keep=False)
    assert not dup.any(), f"one source maps twice to a KEGG ID: {src[src.kegg_id != ''][dup]}"
    return src


def build():
    src = sources()
    src.round(4).to_csv(TABLES / "T01_light_dark_sources.tsv", sep="\t", index=False)
    ref = consensus(src)
    ref.to_csv(TABLES / "T02_light_dark_reference.tsv", sep="\t", index=False, float_format="%.4g")
    counts = ref.status.value_counts().to_dict()
    print(f"  T01: {len(src)} source rows ({src.kegg_id.ne('').sum()} with KEGG IDs, {int(src.qualifies.sum())} qualifying)")
    print(f"  T02: {len(ref)} KEGG compounds; {counts}; up {int((ref.direction > 0).sum())}, down {int((ref.direction < 0).sum())}")


def consensus(src, use=("zhu", "wu_stress", "caldana")):
    """The pre-registered consensus rule over the sources in `use` (the default is all three)."""
    src = src[src.source.isin(use)]
    rows = []
    for kid, g in src[src.kegg_id.fillna("") != ""].groupby("kegg_id"):
        q = g[g.qualifies]
        signs = set(np.sign(q.effect))
        status = ("none" if q.empty else "conflicting" if len(signs) > 1 else
                  "robust" if (len(q) >= 2 or (q.frac_same >= ROBUST_FRAC).any()) else "single")
        direction = 0 if status in ("none", "conflicting") else int(signs.pop())
        rows.append({"kegg_id": kid, "name": g.name.iloc[0], "sources": ",".join(sorted(g.source)),
                     "qualifying": ",".join(sorted(q.source)), "status": status, "direction": direction,
                     "reference_z": round(float(q.z.mean()), 3) if direction else 0.0,
                     **{f"{s}_effect": round(float(g[g.source == s].effect.iloc[0]), 3) if (g.source == s).any() else np.nan
                        for s in ("zhu", "wu_stress", "caldana")},
                     **{f"{s}_frac_same": round(float(g[g.source == s].frac_same.iloc[0]), 3) if (g.source == s).any() and s != "caldana" else np.nan
                        for s in ("zhu", "wu_stress")},
                     "caldana_p": float(g[g.source == "caldana"].p.iloc[0]) if (g.source == "caldana").any() else np.nan})
    return pd.DataFrame(rows).sort_values(["status", "reference_z"], ascending=[True, False])


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("step", choices=["fetch", "kegg", "build"])
    a = ap.parse_args()
    {"fetch": fetch, "kegg": kegg_map, "build": build}[a.step]()
