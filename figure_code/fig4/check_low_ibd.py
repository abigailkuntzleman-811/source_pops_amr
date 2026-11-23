#!/usr/bin/env python3
"""Diagnose individuals with unusually low within-population IBD in Figure 4.

For individuals of --pop whose average within-population IBD (individual_ibd.csv,
from 01_individual_ibd.py) is below --below cM, compared with everyone else in
the population:

    masked_inpop   Figure 4's value (pileup-masked IBD, unrelated partners)
    raw_inpop      the same from the raw, unmasked hap-ibd output
    any_cM         total cM shared with anyone in the 1000 Genomes data (raw)
    any_partners   how many different people they share any IBD with (raw)
    top_pop        the population they share the most cM with (raw)
    parent_of      children they have in the 1000 Genomes pedigree

If raw_inpop is normal but masked_inpop is not, the masking removed their IBD.
If any_cM is low too, they share little IBD with anyone (technical, e.g.
sequencing batch or phasing). If any_cM is normal but they share it with
another population, they are genuinely less related to the rest of --pop.

Usage:
    python figure_code/fig4/check_low_ibd.py --pop CLM
    sbatch -J lowibd run.script figure_code/fig4/check_low_ibd.py --pop CLM
"""

import argparse
import sys
from pathlib import Path

import pandas as pd

_root = next(p for p in Path(__file__).resolve().parents if (p / "config.py").exists())
sys.path.insert(0, str(_root))
sys.path.insert(0, str(_root / "figure_code" / "fig3"))
from config import CHROMOSOMES, DATA_RAW, IBD_DIR  # noqa: E402

_fig3 = __import__("01_within_population_ibd")


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--pop", default="CLM")
    p.add_argument("--below", type=float, default=1.0, help="cM threshold for 'low' (default 1)")
    args = p.parse_args()

    pop_of = _fig3.load_population_map()
    related = _fig3.load_related_samples()
    ind = pd.read_csv(IBD_DIR / "individual_ibd.csv")
    ind = ind[ind["pop"] == args.pop].set_index("individual")["avg_inpop_ibd"]
    members = set(ind.index)
    low = set(ind[ind < args.below].index)
    print(f"{args.pop}: {len(low)} of {len(members)} unrelated individuals below {args.below} cM")
    if not low:
        return

    inpop = {"raw": pd.Series(0.0, index=sorted(members))}
    any_cm = pd.Series(0.0, index=sorted(members))
    partners = {i: set() for i in members}
    by_pop = {i: {} for i in members}
    for chrom in CHROMOSOMES:
        n = chrom.removeprefix("chr")
        d = pd.read_csv(IBD_DIR / f"kgp_chr{n}_segments.ibd.gz", sep="\t", header=None,
                        names=_fig3.IBD_COLUMNS, usecols=["id1", "id2", "cM"])
        d = d[d["id1"].isin(members) | d["id2"].isin(members)]
        for a, b in (("id1", "id2"), ("id2", "id1")):
            mine = d[d[a].isin(members)]
            any_cm = any_cm.add(mine.groupby(a)["cM"].sum(), fill_value=0)
            same = mine[mine[b].isin(members)]
            inpop["raw"] = inpop["raw"].add(same.groupby(a)["cM"].sum(), fill_value=0)
            for i, grp in mine[mine[a].isin(low)].groupby(a):
                partners[i].update(grp[b])
                for q, cm in grp.groupby(grp[b].map(pop_of))["cM"].sum().items():
                    by_pop[i][q] = by_pop[i].get(q, 0) + cm
            for i, grp in mine[~mine[a].isin(low)].groupby(a):
                partners[i].update(grp[b])
        print(f"  read {chrom}", flush=True)

    ped = pd.read_csv(DATA_RAW / "metadata" / "20130606_g1k_3202_samples_ped_population.txt", sep=r"\s+")
    cols = {c.lower(): c for c in ped.columns}
    kids = pd.concat([ped[cols["fatherid"]], ped[cols["motherid"]]]).value_counts()

    table = pd.DataFrame({
        "masked_inpop": ind,
        "raw_inpop": inpop["raw"] / (len(members) - 1),
        "any_cM": any_cm,
        "any_partners": pd.Series({i: len(s) for i, s in partners.items()}),
        "parent_of": pd.Series({i: int(kids.get(i, 0)) for i in members}),
    })
    table["top_pop"] = pd.Series({i: (max(by_pop[i], key=by_pop[i].get) if by_pop[i] else "-") for i in low})

    print(f"\nlow individuals (below {args.below} cM):")
    print(table.loc[sorted(low)].round(2).to_string())
    rest = table.drop(index=list(low))
    print(f"\nmedian of the other {len(rest)}:")
    print(rest.drop(columns="top_pop").median().round(2).to_string())
    print("\nhow the low individuals' IBD splits by partner population (cM, raw):")
    print(pd.DataFrame(by_pop).T.loc[sorted(low)].fillna(0).round(0).astype(int).to_string())


if __name__ == "__main__":
    main()
