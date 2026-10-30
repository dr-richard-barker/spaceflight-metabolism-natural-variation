# Copied from arabidopsis-atlas scripts/wu_modes.py at 69b21c6 (MIT, same author); keep the two in step.
"""One ion mode per Wu 2018 metabolite, shared by every script that reads wu_leaf_identified_long.csv.

Wu et al. (2018) report some metabolites in both positive and negative ion mode, and the two modes'
intensities sit on different levels (7-methylthioheptyl glucosinolate: control positive-mode median
16.7, stress negative-mode median 26.0). Choosing a mode per accession or per condition mixes those
levels inside one column or across a stress - control difference. So the mode is chosen once per
metabolite: a mode measured in both conditions (negative first), else each condition's own mode.

`stress_unpaired` names the metabolites whose stress - control difference is not like for like and
must be left out of stress responses:
  - no ion mode shared by the two conditions;
  - listed on more than one row of Wu Table 2 (D-Galactonic acid as C6H10O7 and C6H12O7; Guanosine
    2':3'-cyclic monophosphate at two retention times), so control and stress can come from
    different peaks.
Within-condition statistics (ranks, percentiles) can still use these metabolites.
"""

import pandas as pd


def mode_choice(wu):
    """(metabolite_name, condition) -> ion mode."""
    avail = wu.groupby(["metabolite_name", "condition"])["mode"].agg(set).to_dict()
    choice = {}
    for name in wu.metabolite_name.unique():
        c, s = avail.get((name, "control"), set()), avail.get((name, "stress"), set())
        shared = [m for m in ("negative", "positive") if m in c & s]
        for cond, have in (("control", c), ("stress", s)):
            if have:
                choice[(name, cond)] = shared[0] if shared else ("negative" if "negative" in have else "positive")
    return choice


def one_mode(wu):
    """The long table restricted to each metabolite's chosen mode: one value per metabolite x accession x condition."""
    choice = mode_choice(wu)
    keep = [choice.get((n, c)) == m for n, c, m in zip(wu.metabolite_name, wu.condition, wu["mode"])]
    out = wu[keep]
    assert not out.duplicated(["metabolite_name", "accession", "condition"]).any(), "one_mode left duplicates"
    return out


def mode_split(wu):
    """Metabolites measured in both conditions but with no ion mode in common."""
    choice = mode_choice(wu)
    return {n for (n, c), m in choice.items()
            if c == "control" and (n, "stress") in choice and choice[(n, "stress")] != m}


def split_rows(ident):
    """Metabolite names that Wu Table 2 lists on more than one row."""
    return set(ident.metabolite_name[ident.metabolite_name.duplicated(keep=False)])


def stress_unpaired(wu, ident):
    return mode_split(wu) | split_rows(ident)


def stress_minus_control(wu, ident):
    """Per accession stress - control on one mode per metabolite, NaN for `stress_unpaired` metabolites.

    Returns a DataFrame indexed by (metabolite_name, accession) with columns control, stress, delta.
    """
    w = one_mode(wu).pivot_table(index=["metabolite_name", "accession"], columns="condition", values="intensity")
    w = w.reindex(columns=["control", "stress"])
    w["delta"] = w.stress - w.control
    bad = w.index.get_level_values(0).isin(sorted(stress_unpaired(wu, ident)))
    w.loc[bad, "delta"] = float("nan")
    return w


if __name__ == "__main__":
    from pathlib import Path
    proc = Path(__file__).resolve().parent.parent / "data" / "processed" / "metabolome"
    wu = pd.read_csv(proc / "wu_leaf_identified_long.csv")
    ident = pd.read_csv(proc / "wu_identified_metabolites.csv")
    print(f"metabolites {wu.metabolite_name.nunique()}; mode split {sorted(mode_split(wu))}; "
          f"split rows {sorted(split_rows(ident))}; rows kept {len(one_mode(wu)):,} of {len(wu):,}")
