"""Check every DOI in manuscript/latex/references.bib against its CrossRef record.

Compares title, year and first-author family name; prints one line per DOI and exits non-zero on
any mismatch, so it can gate a push. No e-mail is sent in the request headers.
"""

import json
import re
import sys
import time
import unicodedata
import urllib.parse
import urllib.request

from paths import ROOT

BIB = ROOT / "manuscript" / "latex" / "references.bib"
UA = {"User-Agent": "spaceflight-metabolism-natural-variation/1.0 (reference check)"}


def plain(s):
    s = re.sub(r"\\[a-zA-Z]+\s*", "", s)
    s = re.sub(r"\\.", "", s)
    s = unicodedata.normalize("NFKD", re.sub(r"<[^>]+>", "", s))
    return re.sub(r"[^a-z0-9]", "", s.encode("ascii", "ignore").decode().lower())


def entries():
    for m in re.finditer(r"@\w+\{(\w+),(.*?)\n\}", BIB.read_text(), re.S):
        f = dict(re.findall(r"(\w+)\s*=\s*\{(.*?)\},?\s*\n", m.group(2) + "\n", re.S))
        if "doi" in f:
            yield m.group(1), f


def main():
    bad = 0
    for key, f in entries():
        url = "https://api.crossref.org/works/" + urllib.parse.quote(f["doi"])
        with urllib.request.urlopen(urllib.request.Request(url, headers=UA), timeout=60) as r:
            c = json.load(r)["message"]
        year = (c.get("published-print") or c.get("published-online") or c["issued"])["date-parts"][0][0]
        first = c["author"][0].get("family", "")
        checks = {"title": plain(f["title"]) == plain(c["title"][0]),
                  "year": str(year) == f["year"].strip(),
                  "first author": plain(f["author"].split(",")[0]) == plain(first),
                  # every author, in order; an entry whose note says its list comes from PubMed
                  # (CrossRef carries only the first author) is checked on the first author alone
                  "authors": [plain(a.split(",")[0]) for a in f["author"].split(" and ")]
                             == [plain(a.get("family", "")) for a in c["author"]]
                             or ("PubMed" in f.get("note", "") and len(c["author"]) == 1)}
        ok = all(checks.values())
        bad += not ok
        print(f"  {'ok ' if ok else 'BAD'} {key:<14} {f['doi']}"
              + ("" if ok else "  " + ", ".join(k for k, v in checks.items() if not v)))
        time.sleep(0.3)
    sys.exit(1 if bad else 0)


if __name__ == "__main__":
    main()
