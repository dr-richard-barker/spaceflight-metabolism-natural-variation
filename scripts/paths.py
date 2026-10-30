"""Paths shared by every script."""

from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "data"
CACHE = DATA / "cache"            # git-ignored raw downloads, re-fetched by the scripts
ATLAS = DATA / "atlas"            # vendored from arabidopsis-atlas
PHOTORESP = DATA / "photoresp"    # vendored from Photorespiration_multiomics_microgravity
OSDR = DATA / "osdr"              # per-contrast gene statistics extracted from GeneLab tables
REF = DATA / "reference"          # derived light/dark reference inputs
TABLES = ROOT / "results" / "tables"
FIGURES = ROOT / "results" / "figures"
PAINTOMICS_RAW = ROOT / "results" / "paintomics_raw"
PREREG = ROOT / "prereg"

# Source repositories and the commits every vendored file is read from (via `git show`), so a
# vendored file can always be traced to one exact blob.
ATLAS_REPO = ROOT.parent / "arabidopsis-atlas"
ATLAS_COMMIT = "69b21c6"
PHOTORESP_REPO = ROOT.parent / "Photorespiration_multiomics_microgravity"
PHOTORESP_COMMIT = "515226d"
