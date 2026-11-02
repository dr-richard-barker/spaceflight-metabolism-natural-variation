"""Downloads for the revision 1 analyses (prereg/revision1): cached in the git-ignored data/cache/ and
checked against the SHA-256 recorded in data/reference/revision1_sources.tsv on first fetch, so a
changed upstream file stops the analysis instead of silently changing a result."""

import hashlib
import urllib.parse
import urllib.request

import pandas as pd

from paths import CACHE, REF

RECORD = REF / "revision1_sources.tsv"
UA = {"User-Agent": "spaceflight-metabolism-natural-variation/1.0"}
OSDR_URL = "https://osdr.nasa.gov/geode-py/ws/studies/{osd}/download?source=datamanager&file={f}"


def get(url, dest):
    if not dest.exists():
        dest.parent.mkdir(parents=True, exist_ok=True)
        with urllib.request.urlopen(urllib.request.Request(url, headers=UA), timeout=600) as r:
            dest.write_bytes(r.read())
    b = dest.read_bytes()
    sha = hashlib.sha256(b).hexdigest()
    key = str(dest.relative_to(CACHE))
    rec = pd.read_csv(RECORD, sep="\t") if RECORD.exists() else pd.DataFrame(columns=["file", "url", "bytes", "sha256"])
    if key in set(rec.file):
        want = rec.set_index("file").loc[key, "sha256"]
        if sha != want:
            raise SystemExit(f"{key}: SHA-256 {sha[:16]} differs from the recorded {want[:16]}; the upstream file changed")
    else:
        rec = pd.concat([rec, pd.DataFrame([{"file": key, "url": url, "bytes": len(b), "sha256": sha}])])
        rec.sort_values("file").to_csv(RECORD, sep="\t", index=False)
        print(f"  recorded {key} ({len(b)} bytes, {sha[:16]})")
    return dest


def osdr(osd, name):
    return get(OSDR_URL.format(osd=osd, f=urllib.parse.quote(name)), CACHE / osd.lower().replace("-", "") / name)
