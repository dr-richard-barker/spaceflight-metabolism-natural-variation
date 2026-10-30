# Adapted from Photorespiration_multiomics_microgravity scripts/figures/style.py at 515226d (MIT, same author).
#!/usr/bin/env python3
"""One visual system for every figure in the manuscript.

A figure montage only reads as a single piece of work if the palette, type and sizing are
decided once. Colours are Okabe-Ito, which stays distinguishable under all common forms of
colour-vision deficiency and prints legibly in greyscale — the two constraints a journal
figure actually has to meet.

Semantic colours matter more than pretty ones here, so they are named for what they mean:
an enclosure class, an illumination state, a predicted-versus-measured layer. The same
quantity keeps the same colour in every panel of every figure.
"""

from __future__ import annotations

import os

# Without a fixed SOURCE_DATE_EPOCH, matplotlib stamps the current time into every PDF, so
# a rerun that changes nothing still produces a different file. Pinning it makes figure
# output byte-reproducible, which is the whole point of shipping the scripts.
os.environ.setdefault("SOURCE_DATE_EPOCH", "1735689600")  # 2025-01-01T00:00:00Z

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

# ---- Okabe-Ito ---------------------------------------------------------------------
BLACK = "#000000"
ORANGE = "#E69F00"
SKY = "#56B4E9"
GREEN = "#009E73"
YELLOW = "#F0E442"
BLUE = "#0072B2"
VERMILION = "#D55E00"
PURPLE = "#CC79A7"

GREY = "#8A8F98"
GREY_LIGHT = "#D6D9DE"
INK = "#1A1D21"

# ---- semantics ---------------------------------------------------------------------
# Enclosure classes, ordered by how much carbon they let the plant keep (CFD T11).
ENCLOSURE = {"sealed": VERMILION, "tape": ORANGE, "vented": BLUE}
ENCLOSURE_ORDER = ["sealed", "tape", "vented"]

# Illumination — the axis that turned out to separate the data.
LIGHT = {"light": ORANGE, "dark": BLUE}

# Direction of response.
UP = VERMILION
DOWN = BLUE
FLAT = GREY

# Verdicts against the model's blind prediction.
MATCH = GREEN
MISMATCH = VERMILION
NO_SIGNAL = GREY

# Data provenance — measured versus model output. Predicted quantities are always the
# lighter, dashed, secondary treatment so a reader can never mistake one for the other.
MEASURED = INK
PREDICTED = PURPLE

# Omic layers.
LAYER = {"transcriptome": BLUE, "proteome": GREEN, "metabolome": PURPLE}

# Gene-set families, so a set keeps its colour across figures 5 and 7.
SETS = {
    "carbon_starvation_DIN": VERMILION,
    "photorespiration_core": GREY,
    "ath00630": GREY,
    "fermentation": SKY,
    "hypoxia_responsive": SKY,
    "photosynthesis_apparatus": BLUE,
    "rubisco": BLUE,
    "ath00710": BLUE,
    "ath00500": GREEN,
    "ath00010": GREEN,
}

PRETTY = {
    "ath00630": "Glyoxylate & dicarboxylate\n(photorespiration)",
    "ath00710": "Carbon fixation",
    "ath00500": "Starch & sucrose",
    "ath00010": "Glycolysis",
    "photorespiration_core": "Photorespiration\nC2 enzymes",
    "carbon_starvation_DIN": "Carbon starvation\n(DIN family)",
    "fermentation": "Fermentation",
    "hypoxia_responsive": "Hypoxia-responsive",
    "photosynthesis_apparatus": "Photosystem &\nlight harvesting",
    "rubisco": "Rubisco",
}


def apply() -> None:
    """Set the rcParams every figure script shares."""
    plt.rcParams.update({
        "figure.dpi": 110,
        "savefig.dpi": 300,
        "savefig.bbox": "tight",
        "savefig.pad_inches": 0.03,
        "font.family": "sans-serif",
        "font.sans-serif": ["Helvetica Neue", "Helvetica", "Arial", "DejaVu Sans"],
        "font.size": 8,
        "axes.titlesize": 9.5,
        "axes.titleweight": "bold",
        "axes.labelsize": 8.5,
        "axes.edgecolor": "#3A3F45",
        "axes.linewidth": 0.8,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.grid": True,
        "grid.color": GREY_LIGHT,
        "grid.linestyle": ":",
        "grid.linewidth": 0.6,
        "grid.alpha": 0.9,
        "xtick.labelsize": 7.5,
        "ytick.labelsize": 7.5,
        "xtick.color": INK,
        "ytick.color": INK,
        "text.color": INK,
        "axes.labelcolor": INK,
        "legend.fontsize": 7.2,
        "legend.frameon": False,
        "lines.linewidth": 1.6,
        "lines.solid_capstyle": "round",
        "patch.linewidth": 0.6,
    })


def panel(ax, letter: str, title: str) -> None:
    """Bold panel letter plus title, left-aligned above the axes — npj house style."""
    ax.set_title(f"{letter}  {title}", loc="left", fontweight="bold", fontsize=9.5, pad=7)


def save(fig, stem: str, outdir: str) -> list[str]:
    """Write PNG (raster, for the web) and PDF (vector, for the journal)."""
    os.makedirs(outdir, exist_ok=True)
    written = []
    for ext in ("png", "pdf"):
        path = os.path.join(outdir, f"{stem}.{ext}")
        fig.savefig(path)
        written.append(path)
    plt.close(fig)
    for p in written:
        print(f"  wrote {p}")
    return written


def predicted_note(fig, text: str = "Metabolite layer is model output, not measured.") -> None:
    """The standing caveat, in the same place on every figure that shows a prediction."""
    fig.text(0.005, 0.004, text, fontsize=6.4, color=GREY, ha="left", va="bottom")
