"""Submit jobs A, B and C to PaintOmics and recover their full records.

Adapted from the first study's scripts/05_submit_paintomics.py and its pa_recover_job helper
(Photorespiration_multiomics_microgravity at 515226d, MIT, same author). Its notes still hold:
the browser form cannot be driven under automation, so step 1 is the same multipart POST the form
makes to /pa_step1, with omicN_origin set to "client" (empty makes the job register zero omics
while still returning success); /check_job_status is one-shot. Step 2 (accept the metabolite
matches and run the pathway analysis) is the "Next step" button (#runButton) on
https://paintomics.org/?jobID=<id>. The finished job is then read whole with POST /pa_recover_job
and vendored gzipped, because guest jobs are deleted after 7 days.

  python3 scripts/14_submit_paintomics.py submit --job A
  (open https://paintomics.org/?jobID=<id> and press "Next step"; wait for the pathways)
  python3 scripts/14_submit_paintomics.py recover --job A --id <id>

AI interpretation is off (aiConsent false), as in the first study: it would forward the data to a
separate LLM endpoint.
"""

import argparse
import gzip
import json
import sys
import time

import requests

from paths import DATA, PAINTOMICS_RAW

BASE = "https://paintomics.org"
UPLOAD = DATA / "paintomics_upload"
SPECIE = "ath"
DATABASES = ["KEGG", "MapMan"]
UA = {"User-Agent": "spaceflight-metabolism-natural-variation/1.0"}

GENES = {"omic_name": "Gene expression", "file_type": "Gene Expression file",
         "relevant_file_type": "Relevant Genes list", "match_type": "Gene", "enrichment": "genes"}
PROTEINS = {"omic_name": "Proteomics", "file_type": "Proteomic quatification",
            "relevant_file_type": "Relevant proteins list", "match_type": "Gene",
            "enrichment": "features"}
METABOLITES = {"omic_name": "Metabolomics", "file_type": "Metabolomic quatification",
               "relevant_file_type": "Relevant Compound list", "match_type": "Compound",
               "enrichment": "features", "file": "metabolomics_values.tab",
               "relevant_file": "metabolomics_relevant.tab"}
LAYER = "ground carbon-starvation reference metabolite layer (not flight measurements)"
JOBS = {
    "A": ([{**GENES, "file": "osd522_gene_expression_values.tab",
            "relevant_file": "osd522_gene_expression_relevant.tab"},
           {**PROTEINS, "file": "osd522_proteomics_values.tab",
            "relevant_file": "osd522_proteomics_relevant.tab"}, METABOLITES],
          f"OSD-522 BRIC-LED-001 Arabidopsis shoots, flight vs ground: OSDR transcriptome + proteome "
          f"with a {LAYER}"),
    "B": ([{**GENES, "file": "osd678_light_genes_values.tab",
            "relevant_file": "osd678_light_genes_relevant.tab"}, METABOLITES],
          f"OSD-678 CARA Arabidopsis leaves, light, flight vs ground: OSDR transcriptome with a {LAYER}"),
    "C": ([{**GENES, "file": "osd678_dark_genes_values.tab",
            "relevant_file": "osd678_dark_genes_relevant.tab"}, METABOLITES],
          f"OSD-678 CARA Arabidopsis leaves, dark, flight vs ground: OSDR transcriptome with a {LAYER}"),
    "D": ([{**GENES, "file": "osd38_genes_values.tab", "relevant_file": "osd38_genes_relevant.tab"},
           {**PROTEINS, "file": "osd38_membrane_proteins_values.tab",
            "relevant_file": "osd38_membrane_proteins_relevant.tab"}, METABOLITES],
          f"OSD-38 BRIC-20 dark-grown Arabidopsis seedlings, flight vs ground: OSDR transcriptome + membrane "
          f"proteome with a {LAYER}"),
}


def build(omics, description):
    data = [("specie", SPECIE), ("jobDescription", description), ("experimentDesign", ""),
            ("aiConsent", "false")]
    data += [("databases[]", db) for db in DATABASES]
    files, handles = [], []
    for i, om in enumerate(omics):
        for key, value in (("omic_name", om["omic_name"]), ("origin", "client"), ("filelocation", ""),
                           ("file_type", om["file_type"]), ("relevant_origin", "client"),
                           ("relevant_filelocation", ""), ("relevant_file_type", om["relevant_file_type"]),
                           ("design_origin", ""), ("design_filelocation", ""),
                           ("match_type", om["match_type"]), ("enrichment", om["enrichment"])):
            data.append((f"omic{i}_{key}", value))
        for field in ("file", "relevant_file"):
            path = UPLOAD / om[field]
            fh = open(path, "rb")
            handles.append(fh)
            files.append((f"omic{i}_{field}", (path.name, fh, "text/plain")))
    return data, files, handles


def poll(session, job_id, timeout=900):
    started = time.time()
    while time.time() - started < timeout:
        try:
            info = session.post(f"{BASE}/check_job_status/{job_id}", timeout=60).json()
        except Exception as exc:  # noqa: BLE001 - transient polling errors are expected
            print(f"    poll error: {exc}")
            time.sleep(10)
            continue
        if info.get("success") and "geneBasedInputOmics" in info:
            for o in (info.get("geneBasedInputOmics") or []) + (info.get("compoundBasedInputOmics") or []):
                print(f"      {o.get('omicName')}: {o.get('omicSummary')}")
            return info
        status = info.get("status") or info.get("state")
        if status in ("finished", "Finished", "done", "error", "Error", "failed"):
            return info
        time.sleep(8)
    return {"status": "timeout"}


def submit(job):
    omics, desc = JOBS[job]
    s = requests.Session()
    s.headers.update(UA)
    s.get(BASE, timeout=60)
    data, files, handles = build(omics, desc)
    try:
        r = s.post(f"{BASE}/pa_step1", data=data, files=files, timeout=900)
    finally:
        for fh in handles:
            fh.close()
    payload = r.json()
    PAINTOMICS_RAW.mkdir(parents=True, exist_ok=True)
    job_id = payload.get("jobID") or payload.get("jobId")
    print(f"  job {job}: HTTP {r.status_code}, jobID {job_id}")
    if not job_id:
        sys.exit(f"  rejected: {json.dumps(payload)[:800]}")
    final = poll(s, job_id)
    with gzip.GzipFile(PAINTOMICS_RAW / f"job{job}_{job_id}_step1.json.gz", "wb", mtime=0) as gz:
        gz.write(json.dumps(final, indent=1).encode())
    print(f"  next: open {BASE}/?jobID={job_id}, press 'Next step', then recover")


def recover(job, job_id):
    s = requests.Session()
    s.headers.update(UA)
    r = s.post(f"{BASE}/pa_recover_job", data={"jobID": job_id}, timeout=300)
    r.raise_for_status()
    payload = r.json()
    n = len(payload.get("pathwaysInfo") or [])
    if not n:
        sys.exit(f"  job {job_id} has no pathways yet (step 2 not run, or still running)")
    out = PAINTOMICS_RAW / f"job{job}_{job_id}_full.json.gz"
    with gzip.GzipFile(out, "wb", mtime=0) as gz:
        gz.write(json.dumps(payload, sort_keys=True).encode())
    print(f"  job {job} {job_id}: {n} pathways -> {out.name} ({out.stat().st_size:,} B)")


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("step", choices=["submit", "recover"])
    ap.add_argument("--job", choices=sorted(JOBS), required=True)
    ap.add_argument("--id")
    a = ap.parse_args()
    submit(a.job) if a.step == "submit" else recover(a.job, a.id)
