#!/usr/bin/env python3
"""Independent check of the S7 table: recompute a sample of haplotype pairs by hand.

For each sampled pair, straight from the LAI BED file and the IBD segments,
with plain interval intersection (none of 01's code):
    shared_bp  bp where both haplotypes carry ancestry a
    ibd_bp     bp of IBD between the two haplotypes inside that shared ancestry
and compares them with ancestry_aware_ibd_<POP>.csv.gz written by 01.

Usage:
    python figure_code/figS7_ancestry_aware_ibd/02_check_pairs.py           # CLM, 20 pairs
    python figure_code/figS7_ancestry_aware_ibd/02_check_pairs.py --pop PEL --n 30
"""

import argparse
import sys
from pathlib import Path

import pandas as pd

_root = next(p for p in Path(__file__).resolve().parents if (p / "config.py").exists())
sys.path.insert(0, str(_root))
sys.path.insert(0, str(_root / "figure_code" / "fig3"))
from config import IBD_DIR, IBD_MASKED, LAI_BED_DIR, TABLE_DIR  # noqa: E402

IBD_COLUMNS = __import__("01_within_population_ibd").IBD_COLUMNS
ANCESTRIES = ["EUR", "AMR", "AFR"]


def merge(iv):
    out = []
    for s, e in sorted(iv):
        if out and s <= out[-1][1]:
            out[-1][1] = max(out[-1][1], e)
        else:
            out.append([s, e])
    return out


def intersect(a, b):
    """Intersection of two lists of [start, end) intervals (each merged)."""
    out, i, j = [], 0, 0
    while i < len(a) and j < len(b):
        s, e = max(a[i][0], b[j][0]), min(a[i][1], b[j][1])
        if s < e:
            out.append([s, e])
        if a[i][1] < b[j][1]:
            i += 1
        else:
            j += 1
    return out


def length(iv):
    return sum(e - s for s, e in iv)


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--pop", default="CLM")
    p.add_argument("--n", type=int, default=20, help="pairs to check (half with IBD, half random)")
    p.add_argument("--seed", type=int, default=1)
    args = p.parse_args()

    table = pd.read_csv(TABLE_DIR / f"ancestry_aware_ibd_{args.pop}.csv.gz")
    with_ibd = table[table[[f"{a}_ibd_bp" for a in ANCESTRIES]].sum(axis=1) > 0]
    pairs = pd.concat([with_ibd.sample(min(args.n // 2, len(with_ibd)), random_state=args.seed),
                       table.sample(args.n - min(args.n // 2, len(with_ibd)), random_state=args.seed + 1)])

    lai = pd.read_csv(LAI_BED_DIR / f"{args.pop}_lai.bed", sep="\t")
    lai["hap"] = lai["name"] + "_" + lai["haplotype"].astype(str).str.replace("hap", "", regex=False)
    haps = set(pairs["hap1"]) | set(pairs["hap2"])
    lai = lai[lai["hap"].isin(haps)]

    ids = {h.rsplit("_", 1)[0] for h in haps}
    tag = ".masked" if IBD_MASKED else ""
    ibd = []
    for n in range(1, 23):
        d = pd.read_csv(IBD_DIR / f"kgp_chr{n}_segments{tag}.ibd.gz", sep="\t", header=None, names=IBD_COLUMNS)
        ibd.append(d[d["id1"].isin(ids) & d["id2"].isin(ids)])
    ibd = pd.concat(ibd)
    ibd["h1"] = ibd["id1"] + "_" + ibd["hap1"].astype(str)
    ibd["h2"] = ibd["id2"] + "_" + ibd["hap2"].astype(str)

    rows, worst = [], 0
    for r in pairs.itertuples(index=False):
        seg = ibd[((ibd["h1"] == r.hap1) & (ibd["h2"] == r.hap2)) | ((ibd["h1"] == r.hap2) & (ibd["h2"] == r.hap1))]
        out = {"hap1": r.hap1, "hap2": r.hap2}
        for a in ANCESTRIES:
            shared = ibd_bp = 0
            for chrom in sorted(set(lai["chrom"])):
                ta = merge(lai.loc[(lai["hap"] == r.hap1) & (lai["chrom"] == chrom) & (lai["ancestry"] == a),
                                   ["start", "end"]].values.tolist())
                tb = merge(lai.loc[(lai["hap"] == r.hap2) & (lai["chrom"] == chrom) & (lai["ancestry"] == a),
                                   ["start", "end"]].values.tolist())
                both = intersect(ta, tb)
                shared += length(both)
                s = merge(seg.loc[seg["chrom"] == chrom, ["start_bp", "end_bp"]].values.tolist())
                ibd_bp += length(intersect(both, s))
            want_s, want_i = getattr(r, f"{a}_shared_bp"), getattr(r, f"{a}_ibd_bp")
            out[f"{a}_shared"] = f"{shared:,}" + ("" if shared == want_s else f" (table {want_s:,})")
            out[f"{a}_ibd"] = f"{ibd_bp:,}" + ("" if ibd_bp == want_i else f" (table {want_i:,})")
            worst = max(worst, abs(shared - want_s), abs(ibd_bp - want_i))
        rows.append(out)

    print(pd.DataFrame(rows).to_string(index=False))
    if worst == 0:
        print(f"\nOK: all {len(rows)} pairs match the S7 table exactly (shared ancestry and IBD bp, EUR/AMR/AFR).")
    else:
        print(f"\nMISMATCH: largest difference {worst:,} bp; mismatching values are shown with the table's value.")


if __name__ == "__main__":
    main()
