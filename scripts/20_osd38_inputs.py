"""OSD-38 (BRIC-20) inputs: the flight proteome and two ground contrasts the study also holds.

OSD-38 grew Col-0 seedlings for about 3 days in darkness (etiolated) in sealed BRIC-PDFU canisters, on
the ISS and on the ground (Kruse et al. 2020), with RNA-seq and iTRAQ proteomics of a membrane and a
soluble fraction. The same study holds a ground hardware experiment, BRIC-PDFU against plates in closed
containers at Ohio University (Basu et al. 2017), and a preservation experiment, RNAlater against liquid
nitrogen (Kruse et al. 2017). Steps:

  fetch      the two processed iTRAQ tables from OSDR into data/cache/osd38 (CC0), checked against the
             SHA-256 recorded in data/reference/osd38_sources.tsv
  contrasts  GeneLab statistics for BRIC vs container (GC v RNAlat) from the same GeneLab table as
             data/osdr/contrast_OSD-38.csv.gz
  coverage   which loci each fraction quantifies under the paper's rule (at least two distinct peptides
             with score > 20), WITHOUT computing any flight/ground comparison, so the protein-level
             predictions can be fixed before a fold change is seen
  quantify   protein-level flight/ground log2 fold changes (run only after prereg/osd38/ is committed)

THE PROCESSED TABLES ARE TRUNCATED. Both CSVs (and the .xlsx copies) stop at 65,535 data rows, the old
Excel row limit:
  soluble   rows run in protein-group order, groups 1-438; the cut falls inside group 438. Only groups
            1-437 are used: the highest-ranked soluble proteins, not the soluble proteome.
  membrane  block 1 (rows 1-48,799) covers protein groups 1-2,948 and ends before the cut; a second
            block restarts at group 1 with different PSMs for the same groups and is cut at group 754.
            Only block 1 is used.
"""

import argparse
import gzip
import hashlib
import io
import urllib.parse
import urllib.request

import numpy as np
import pandas as pd

from paths import CACHE, OSDR, PHOTORESP_REPO, REF

OSD38 = CACHE / "osd38"
URL = "https://osdr.nasa.gov/geode-py/ws/studies/OSD-38/download?source=datamanager&file={}"
UA = {"User-Agent": "spaceflight-metabolism-natural-variation/1.0"}
FILES = {"soluble": "GLDS-38_soluble_proteome_WyattLab_soluble.csv",
         "membrane": "GLDS-38_membrane_proteome_WyattLab_membrane.csv"}
FLIGHT = ["S1:114.1 Intensity", "S2:116.1 Intensity", "S3:118.1 Intensity"]
GROUND = ["GC1:115.1 Intensity", "GC2:117.1 Intensity", "GC3:119.1 Intensity"]
MIN_SCORE, MIN_PEPTIDES = 20, 2
SOLUBLE_LAST_COMPLETE = 437
MEMBRANE_BLOCK1_ROWS = 48_799
GENELAB_DGE = PHOTORESP_REPO / "data" / "cache" / "dge_OSD-38.csv"
CONTRASTS = {"OSD-38-hardware": ("GC", "RNAlat", "BRIC-PDFU (ground) vs closed container")}


def sha(b):
    return hashlib.sha256(b).hexdigest()


def fetch():
    src = REF / "osd38_sources.tsv"
    rows = []
    for frac, name in FILES.items():
        out = OSD38 / name
        if not out.exists():
            OSD38.mkdir(parents=True, exist_ok=True)
            with urllib.request.urlopen(urllib.request.Request(URL.format(urllib.parse.quote(name)),
                                                               headers=UA), timeout=600) as r:
                out.write_bytes(r.read())
        b = out.read_bytes()
        rows.append({"fraction": frac, "file": name, "bytes": len(b), "sha256": sha(b),
                     "data_rows": b.count(b"\n") - 1})
    got = pd.DataFrame(rows)
    if src.exists():
        want = pd.read_csv(src, sep="\t")
        for r in got.itertuples():
            ok = r.sha256 == want.set_index("file").loc[r.file, "sha256"]
            print(f"  {r.file}: {'ok' if ok else 'CHANGED since recorded'} ({r.data_rows} data rows)")
    else:
        got.to_csv(src, sep="\t", index=False)
        print(got.to_string(index=False))


def psms(frac):
    d = pd.read_csv(OSD38 / FILES[frac])
    if frac == "soluble":
        d = d[d["Protein Group"] <= SOLUBLE_LAST_COMPLETE]
    else:
        d = d.iloc[:MEMBRANE_BLOCK1_ROWS]
        assert not (d["Protein Group"].diff() < 0).any(), "membrane block 1 is not one ordered block"
    d = d.assign(locus=d["Sequence Id"].str.extract(r"^(AT[1-5CM]G\d{5})", expand=False))
    d = d[(d.Score > MIN_SCORE) & (d[FLIGHT + GROUND] > 0).all(axis=1) & d.locus.notna()]
    peptides = d.groupby("locus")["Peptide Sequence"].nunique()
    return d[d.locus.isin(peptides.index[peptides >= MIN_PEPTIDES])]


def coverage():
    rows = []
    for frac in FILES:
        d = psms(frac)
        for locus, g in d.groupby("locus"):
            rows.append({"fraction": frac, "locus": locus, "psms": len(g),
                         "peptides": g["Peptide Sequence"].nunique()})
    cov = pd.DataFrame(rows)
    cov.to_csv(REF / "osd38_protein_coverage.tsv", sep="\t", index=False)
    print(cov.groupby("fraction").agg(loci=("locus", "size"), psms=("psms", "sum")).to_string())


def contrasts():
    raw = GENELAB_DGE.read_bytes()
    df = pd.read_csv(io.BytesIO(raw), low_memory=False)
    lines = []
    for name, (a, b, what) in CONTRASTS.items():
        tag = f"({a})v({b})"
        cols = {f"Log2fc_{tag}": "log2fc", f"Stat_{tag}": "stat", f"P.value_{tag}": "pvalue",
                f"Adj.p.value_{tag}": "padj"}
        out = df[["TAIR", *cols]].rename(columns=cols)
        for arm, pat in ((a, f"_sl_{a}_Rep"), (b, f"_sl_{b}_Rep")):
            for c in [c for c in df.columns if pat in c]:
                out[f"{arm}:{c}"] = df[c]
        out = out.dropna(subset=["log2fc"]).drop_duplicates("TAIR")
        buf = io.BytesIO()
        with gzip.GzipFile(fileobj=buf, mode="wb", mtime=0) as gz:
            gz.write(out.to_csv(index=False).encode())
        (OSDR / f"contrast_{name}.csv.gz").write_bytes(buf.getvalue())
        lines.append(f"| `contrast_{name}.csv.gz` | {what}: `{tag}` columns | {len(out)} | `{sha(raw)[:16]}` |")
        print(f"  contrast_{name}.csv.gz: {len(out)} genes")
    prov = OSDR / "PROVENANCE.md"
    text = prov.read_text()
    marker = "\n## OSD-38 ground contrasts\n"
    text = text.split(marker)[0].rstrip("\n") + "\n" + marker + (
        "\nFrom the same GeneLab table as `contrast_OSD-38.csv.gz` (the first study's cached "
        "`dge_OSD-38.csv`), by `scripts/20_osd38_inputs.py`. GC = BRIC-PDFU ground control; RNAlat = "
        "plates in closed containers at Ohio University, fixed in RNAlater like GC (Basu et al. 2017).\n\n| File | Contrast | Genes | Source SHA-256 (first 16) |\n"
        "|---|---|---|---|\n" + "\n".join(lines) + "\n")
    prov.write_text(text)


def quantify():
    """Protein log2 fold changes, flight vs ground, per fraction (method fixed in prereg/osd38)."""
    out = []
    for frac in FILES:
        d = psms(frac)
        x = np.log2(d[FLIGHT + GROUND].astype(float))
        x = x - x.median(axis=0)                    # channel (loading) normalisation
        x = x.sub(x.mean(axis=1), axis=0)           # each PSM relative to its own mean
        x["locus"] = d.locus.to_numpy()
        prot = x.groupby("locus").median()
        f, g = prot[FLIGHT].to_numpy(), prot[GROUND].to_numpy()
        from scipy import stats
        t, p = stats.ttest_ind(f, g, axis=1, equal_var=True)
        res = pd.DataFrame({"fraction": frac, "locus": prot.index, "log2fc": f.mean(1) - g.mean(1),
                            "t": t, "p": p, "psms": d.groupby("locus").size().reindex(prot.index).to_numpy()})
        out.append(res)
    return pd.concat(out, ignore_index=True)


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("step", choices=["fetch", "contrasts", "coverage"])
    {"fetch": fetch, "contrasts": contrasts, "coverage": coverage}[ap.parse_args().step]()
