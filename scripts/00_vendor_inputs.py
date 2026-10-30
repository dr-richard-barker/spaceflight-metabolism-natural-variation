"""Vendor every input this study uses, each traced to one exact source.

- data/atlas/      files read with `git show` from arabidopsis-atlas at ATLAS_COMMIT
- data/photoresp/  files read with `git show` from Photorespiration_multiomics_microgravity at
                   PHOTORESP_COMMIT (the first study: its PaintOmics bundle, job record, predictions)
- data/osdr/       per-contrast gene statistics extracted from the GeneLab differential-expression
                   tables the first study downloaded from NASA OSDR (CC0); only the matched
                   flight-vs-ground columns are kept, plus the matched samples' normalised counts
                   where the sample names say which arm they belong to

Each folder gets a PROVENANCE.md listing source, path, commit or URL, bytes and SHA-256. The
sibling repositories are only needed to re-vendor; everything downstream reads data/.
"""

import gzip
import hashlib
import io
import subprocess
import urllib.parse

import pandas as pd

from paths import (ATLAS, ATLAS_COMMIT, ATLAS_REPO, OSDR, PHOTORESP, PHOTORESP_COMMIT,
                   PHOTORESP_REPO)

ATLAS_FILES = [
    "data/processed/metabolome/darkness_metabolome.csv",
    "data/processed/metabolome/wu_leaf_identified_long.csv",
    "data/processed/metabolome/wu_identified_metabolites.csv",
    "data/processed/metabolome/wu_accessions.csv",
    "data/processed/metabolome/wu_gwas_loci.csv",
    "data/processed/metabolome/wu2016_primary_long.csv",
    "data/processed/metabolome/wu2016_primary_metabolites.csv",
    "data/processed/metabolome/accessions_unified.csv",
    "results/metabolome_meta_analysis/wu_condition_regions.tsv",
    "results/metabolome_meta_analysis/zhu_wu_metabolite_matches.tsv",
    "data/raw/GLDS-522_rna_seq_RSEM_Unnormalized_Counts_GLbulkRNAseq.csv",
    "data/raw/GLDS-522_rna_seq_bulkRNASeq_v2_runsheet.csv",
]
PHOTORESP_FILES = [
    "data/paintomics_upload/gene_expression_values.tab",
    "data/paintomics_upload/gene_expression_relevant.tab",
    "data/paintomics_upload/gene_expression_METHOD.txt",
    "data/paintomics_upload/proteomics_values.tab",
    "data/paintomics_upload/proteomics_relevant.tab",
    "data/paintomics_upload/proteomics_METHOD.txt",
    "data/paintomics_upload/metabolomics_values.tab",
    "data/paintomics_upload/metabolomics_relevant.tab",
    "data/study_registry.tsv",
    "results/paintomics_raw/job1_oK2h66740Q_full.json.gz",
    "results/tables/T06_predicted_compounds.tsv",
    "results/tables/T07_falsification_osd522.tsv",
    "results/tables/T09_paintomics_carbon.tsv",
    "results/tables/T10_metabolite_hubs.tsv",
    "results/tables/T11_hardware_ladder.tsv",
]
# GeneLab differential-expression tables (downloaded by the first study's 07_hardware_ladder.py
# into its git-ignored data/cache/). Sample columns are kept only where the name states the arm.
OSDR_URL = "https://osdr.nasa.gov/geode-py/ws/studies/{osd}/download?source=datamanager&file={f}"
SAMPLE_ARMS = {
    "OSD-38": {"flight": "_FLT_Rep", "ground": "_GC_Rep"},
    "OSD-678-light": {"flight": "Flight_Col-0_Light_", "ground": "Ground_control_Col-0_Light_"},
    "OSD-678-dark": {"flight": "Flight_Col-0_Dark_", "ground": "Ground_control_Col-0_Dark_"},
}


def sha(b):
    return hashlib.sha256(b).hexdigest()


def git_show(repo, commit, path):
    return subprocess.run(["git", "-C", str(repo), "show", f"{commit}:{path}"], check=True,
                          capture_output=True).stdout


def vendor_git(repo, commit, files, dest, title):
    dest.mkdir(parents=True, exist_ok=True)
    rows = []
    for path in files:
        blob = git_show(repo, commit, path)
        out = dest / path.split("/")[-1]
        out.write_bytes(blob)
        rows.append(f"| `{out.name}` | `{path}` | {len(blob)} | `{sha(blob)[:16]}` |")
    head = (f"# {title}\n\nRead with `git show {commit}:<path>` from `{repo.name}` by "
            "`scripts/00_vendor_inputs.py`. Do not edit; re-vendor instead.\n\n"
            "| File | Source path | Bytes | SHA-256 (first 16) |\n|---|---|---|---|\n")
    (dest / "PROVENANCE.md").write_text(head + "\n".join(rows) + "\n")
    print(f"  {dest.name}: {len(rows)} files from {repo.name}@{commit}")


def vendor_osdr():
    OSDR.mkdir(parents=True, exist_ok=True)
    reg = pd.read_csv(PHOTORESP / "study_registry.tsv", sep="\t", comment="#")
    reg = reg[(reg["include"] == "yes") & (reg.accession != "OSD-522")]
    rows = []
    for r in reg.itertuples():
        osd = "-".join(r.accession.split("-")[:2])
        src = PHOTORESP_REPO / "data" / "cache" / f"dge_{osd}.csv"
        raw = src.read_bytes()
        df = pd.read_csv(io.BytesIO(raw), low_memory=False)
        tag = r.contrast.split("Log2fc_", 1)[1]
        cols = {f"Log2fc_{tag}": "log2fc", f"Stat_{tag}": "stat", f"P.value_{tag}": "pvalue",
                f"Adj.p.value_{tag}": "padj"}
        missing = [c for c in cols if c not in df.columns]
        assert not missing, f"{r.accession}: missing {missing}"
        out = df[["TAIR", *cols]].rename(columns=cols)
        arms = SAMPLE_ARMS.get(r.accession)
        n_samples = ""
        if arms:
            for arm, pat in arms.items():
                samp = [c for c in df.columns if pat in c]
                assert len(samp) >= 3, f"{r.accession}: {arm} samples {samp}"
                for c in samp:
                    out[f"{arm}:{c}"] = df[c]
            n_samples = " / ".join(f"{a} {sum(1 for c in out.columns if c.startswith(a + ':'))}" for a in arms)
        out = out.dropna(subset=["log2fc"]).drop_duplicates("TAIR")
        name = f"contrast_{r.accession}.csv.gz"
        buf = io.BytesIO()
        with gzip.GzipFile(fileobj=buf, mode="wb", mtime=0) as gz:
            gz.write(out.to_csv(index=False).encode())
        (OSDR / name).write_bytes(buf.getvalue())
        url = OSDR_URL.format(osd=osd, f=urllib.parse.quote(r.dge_file))
        rows.append(f"| `{name}` | {r.accession} | `{r.dge_file}` | {len(out)} | {n_samples or 'none'} | "
                    f"`{sha(raw)[:16]}` |")
        print(f"  {name}: {len(out)} genes, samples {n_samples or 'none'}")
    head = ("# Per-contrast gene statistics (NASA OSDR, CC0)\n\n"
            "Extracted by `scripts/00_vendor_inputs.py` from the GeneLab differential-expression tables "
            "the first study downloaded (its `07_hardware_ladder.py`, URL pattern\n"
            f"`{OSDR_URL}`). Each file keeps TAIR, the matched flight-vs-ground `log2fc`, `stat` (Wald), "
            "`pvalue` and `padj` columns named in `data/photoresp/study_registry.tsv`, plus the "
            "normalised counts of the matched samples where their names state the arm. OSD-522 is "
            "not here: it is computed from raw counts (data/atlas/).\n\n"
            "| File | Contrast | GeneLab file | Genes | Samples kept | Source SHA-256 (first 16) |\n"
            "|---|---|---|---|---|---|\n")
    (OSDR / "PROVENANCE.md").write_text(head + "\n".join(rows) + "\n")


def main():
    vendor_git(ATLAS_REPO, ATLAS_COMMIT, ATLAS_FILES, ATLAS,
               "Vendored from arabidopsis-atlas (natural-variation metabolome, OSD-522 counts)")
    vendor_git(PHOTORESP_REPO, PHOTORESP_COMMIT, PHOTORESP_FILES, PHOTORESP,
               "Vendored from Photorespiration_multiomics_microgravity (the first study)")
    vendor_osdr()


if __name__ == "__main__":
    main()
