"""Set the first study's 21 compound predictions against the measured light/dark reference (T03, fig01).

The rule is the pre-registered one (prereg/settings.json, compound_rule): predictions driven by net
assimilation (A) or carbon starvation are compared with the T02 consensus direction; Vo- and
RuBP-driven predictions are listed as not comparable, because darkness stops oxygenation and
carboxylation together. T02 was seen before the rule was written, so this is a consistency check.
The table also lists every reference compound with a directional response that the model did not
predict.
"""

import json

import numpy as np
import pandas as pd

import style
from paths import FIGURES, PHOTORESP, PREREG, TABLES

BLOCK = {"Vo": "C2 cycle (Vo)", "A": "Photosynthate (A)", "RuBP": "RuBP",
         "starvation": "Starvation amino acids"}


def build():
    rule = json.loads((PREREG / "settings.json").read_text())["compound_rule"]
    model = pd.read_csv(PHOTORESP / "T06_predicted_compounds.tsv", sep="\t", comment="#")
    ref = pd.read_csv(TABLES / "T02_light_dark_reference.tsv", sep="\t").set_index("kegg_id")
    cols = ["status", "direction", "reference_z", "zhu_effect", "wu_stress_effect", "caldana_effect"]
    rows = []
    for r in model.itertuples():
        hit = ref.loc[r.kegg_id] if r.kegg_id in ref.index else None
        directional = hit is not None and hit.status in ("robust", "single")
        if r.driver not in rule["scored_drivers"]:
            verdict = "not comparable"
        elif not directional:
            verdict = "untested"
        else:
            verdict = "agrees" if np.sign(r.log2FC_FLT_vs_GC) == hit.direction else "disagrees"
        rows.append({"kegg_id": r.kegg_id, "compound": r.compound, "block": BLOCK[r.driver],
                     "model_log2FC": r.log2FC_FLT_vs_GC, "in_model": True,
                     **({c: hit[c] for c in cols} if hit is not None else {"status": "absent"}),
                     "verdict": verdict})
    extra = ref[ref.status.isin(["robust", "single"]) & ~ref.index.isin(model.kegg_id)]
    for kid, h in extra.iterrows():
        rows.append({"kegg_id": kid, "compound": h["name"].strip(), "block": "not in model",
                     "model_log2FC": np.nan, "in_model": False, **{c: h[c] for c in cols},
                     "verdict": "not predicted"})
    t = pd.DataFrame(rows)
    t.to_csv(TABLES / "T03_model_vs_reference.tsv", sep="\t", index=False, float_format="%.4g")
    m = t[t.in_model]
    print(f"  T03: {len(m)} model compounds {m.verdict.value_counts().to_dict()}; "
          f"{len(t) - len(m)} directional reference compounds not predicted "
          f"({int((t[~t.in_model].status == 'robust').sum())} robust)")
    return t


def figure(t):
    style.apply()
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D

    m = t[t.in_model].copy()
    order = list(BLOCK.values())
    m["o"] = m.block.map(order.index)
    m = m.sort_values(["o", "compound"]).reset_index(drop=True)
    x = t[~t.in_model & (t.status == "robust")].sort_values("reference_z").reset_index(drop=True)

    fig = plt.figure(figsize=(7.2, 5.6))
    gs = fig.add_gridspec(1, 5, width_ratios=[0.12, 1, 0.42, 0.1, 1], wspace=0.04)
    am, a, av, b = (fig.add_subplot(gs[0, i]) for i in (0, 1, 2, 4))
    ticks = [-50, -10, -1, 0, 1, 10, 50]
    for ax in (a, b):
        ax.set_xscale("symlog", linthresh=1)
        ax.set_xlim(-80, 80)
        ax.set_xticks(ticks)
        ax.set_xticklabels([str(v) for v in ticks])
        ax.axvline(0, color=style.GREY, lw=0.8)
        ax.grid(axis="y", visible=False)
        ax.set_xlabel("Darkness response (z, symlog)")

    y = np.arange(len(m))[::-1]
    ylim = (-0.7, len(m) - 0.1)
    for ax in (am, a, av):
        ax.set_ylim(*ylim)
    verdict_ink = {"agrees": style.MATCH, "disagrees": style.MISMATCH}
    for i, blk in enumerate(order):
        rows = y[m.block == blk]
        for ax in (am, a, av):
            ax.axhspan(rows.min() - 0.5, rows.max() + 0.5, color=style.GREY_LIGHT,
                       alpha=0.45 if i % 2 else 0, zorder=0, lw=0)
        a.text(0.02, rows.max() + 0.42, blk, transform=a.get_yaxis_transform(), fontsize=6.3,
               color=style.GREY, va="center", style="italic")
    for yi, r in zip(y, m.itertuples()):
        sign = np.sign(r.model_log2FC) if abs(r.model_log2FC) > 0.01 else 0
        am.scatter(0.5, yi, marker={1: "^", -1: "v", 0: "o"}[sign], s=30, color=style.PREDICTED, zorder=3)
        if r.status in ("robust", "single"):
            col = style.UP if r.direction > 0 else style.DOWN
            a.plot([0, r.reference_z], [yi, yi], color=col, lw=1.6, solid_capstyle="butt")
            a.scatter(r.reference_z, yi, s=26, color=col, zorder=3,
                      marker="o" if r.status == "robust" else "D")
        av.text(0.06, yi, r.verdict, va="center", fontsize=6.8,
                color=verdict_ink.get(r.verdict, style.GREY),
                fontweight="bold" if r.verdict in verdict_ink else "normal")
    am.set_xlim(0, 1)
    am.set_xticks([0.5])
    am.set_xticklabels(["model"], fontsize=6.8)
    am.set_yticks(y)
    am.set_yticklabels(m.compound, fontsize=6.8)
    am.grid(False)
    a.set_yticks([])
    av.axis("off")
    for ax in (am,):
        ax.spines["bottom"].set_visible(False)
        ax.tick_params(axis="x", length=0)
    style.panel(am, "a", "First study's predictions")

    yb = np.arange(len(x))
    b.barh(yb, x.reference_z, color=[style.UP if d > 0 else style.DOWN for d in x.direction], height=0.62)
    b.set_yticks(yb)
    b.set_yticklabels(x.compound, fontsize=6.4)
    b.yaxis.tick_right()
    b.spines["left"].set_visible(False)
    b.spines["right"].set_visible(True)
    b.set_ylim(-0.7, len(x) - 0.3)
    style.panel(b, "b", "Robust darkness responses not predicted")

    handles = [Line2D([], [], marker="^", ls="", color=style.PREDICTED, label="Model: up"),
               Line2D([], [], marker="v", ls="", color=style.PREDICTED, label="Model: down"),
               Line2D([], [], marker="o", ls="", color=style.PREDICTED, label="Model: about 0"),
               Line2D([], [], marker="o", ls="", color=style.UP, label="Darkness: up"),
               Line2D([], [], marker="o", ls="", color=style.DOWN, label="Darkness: down"),
               Line2D([], [], marker="D", ls="", color=style.INK, label="Single source (else robust)")]
    fig.legend(handles=handles, ncol=3, loc="upper center", bbox_to_anchor=(0.45, 0.02), fontsize=6.8)
    fig.text(0.0, -0.075, "Model column: the first study's predictions (model output). Lines and bars: "
             "the measured ground light/dark reference (T02), not flight data.", fontsize=6.4,
             color=style.GREY, ha="left")
    style.save(fig, "fig01_model_vs_reference", str(FIGURES))


if __name__ == "__main__":
    figure(build())
