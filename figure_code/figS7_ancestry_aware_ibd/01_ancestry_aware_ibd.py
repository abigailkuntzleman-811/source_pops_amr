#!/usr/bin/env python3
"""Ancestry-aware IBD ratio (supplementary figure S7).

For two haplotypes i and j from different individuals in the same population,
and a continental ancestry a (EUR, AMR or AFR):

    ratio_a(i, j) = |IBD_ij  ∩  A_a,i  ∩  A_a,j|  /  |A_a,i  ∩  A_a,j|

the fraction of the genome where both haplotypes carry ancestry a that is
covered by IBD between them (all in bp). Undefined when the two haplotypes
share no ancestry-a sequence.

Inputs
    IBD       hap-ibd segments, pileup-masked by default (IBD_MASKED), read via
              figure_code/fig3/01_within_population_ibd.read_ibd_segments
    ancestry  $LAI_BED_DIR/<POP>_lai.bed (per-haplotype local ancestry, columns
              chrom start end ancestry name haplotype)
    samples   unrelated individuals only (RELATED_EXCLUDE, e.g. KING --unrelated),
              restricted to those with local ancestry calls

Haplotype 1 and 2 in the IBD file are matched to hap1 and hap2 in the LAI file,
which assumes both were derived from the same phased haplotypes.

Usage:
    python figure_code/figS7_ancestry_aware_ibd/01_ancestry_aware_ibd.py
    python figure_code/figS7_ancestry_aware_ibd/01_ancestry_aware_ibd.py --replot   # redraw only

Output ($TABLE_DIR, and the figure in $FIGURE_DIR):
    ancestry_aware_ibd_<POP>.csv.gz     one row per haplotype pair: shared bp and
                                        IBD bp on shared ancestry, and the ratio,
                                        for each ancestry
    ancestry_aware_ibd_summary.csv      per population and ancestry: pairs, pairs
                                        sharing the ancestry, pairs with no IBD on
                                        it, mean ratio
    ancestry_aware_ibd_ratios.png       distributions of non-zero ratios
"""

import argparse
import sys
from itertools import combinations
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

_root = next(p for p in Path(__file__).resolve().parents if (p / "config.py").exists())
sys.path.insert(0, str(_root))
sys.path.insert(0, str(_root / "figure_code" / "fig3"))
from config import AMR_POPS, CHROMOSOMES, FIGURE_DIR, LAI_BED_DIR, TABLE_DIR  # noqa: E402

_fig3 = __import__("01_within_population_ibd")

ANCESTRIES = ["EUR", "AMR", "AFR"]


def load_lai(pop, keep_ids):
    path = LAI_BED_DIR / f"{pop}_lai.bed"
    if not path.exists():
        sys.exit(f"Missing {path}\nMake the per-population LAI files first (see Table 1 in the README).")
    lai = pd.read_csv(path, sep="\t")
    lai = lai[lai["name"].isin(keep_ids)].copy()
    lai["hap"] = lai["name"] + "_" + lai["haplotype"].astype(str).str.replace("hap", "", regex=False)
    return lai


def ancestry_grid(lai_chrom, haps):
    """Common breakpoint grid for one chromosome and each haplotype's ancestry on it.

    Returns (edges, codes): edges has G+1 positions, codes is an H x G array of
    indices into ANCESTRIES, -1 where a haplotype has no call.
    """
    edges = np.unique(np.concatenate([lai_chrom["start"].to_numpy(), lai_chrom["end"].to_numpy()]))
    codes = np.full((len(haps), len(edges) - 1), -1, dtype=np.int8)
    row = {h: k for k, h in enumerate(haps)}
    code = {a: k for k, a in enumerate(ANCESTRIES)}
    for r in lai_chrom.itertuples(index=False):
        if r.ancestry not in code or r.hap not in row:
            continue
        a, b = np.searchsorted(edges, [r.start, r.end])
        codes[row[r.hap], a:b] = code[r.ancestry]
    return edges, codes


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--pops", nargs="+", default=AMR_POPS)
    p.add_argument("--replot", action="store_true",
                   help="redraw from the saved ratios (ancestry_aware_ibd_ratios.npz) without recomputing")
    args = p.parse_args()

    if args.replot:
        saved = TABLE_DIR / "ancestry_aware_ibd_ratios.npz"
        if not saved.exists():
            sys.exit(f"Missing {saved}; run once without --replot first.")
        d = np.load(saved)
        FIGURE_DIR.mkdir(parents=True, exist_ok=True)
        plot(args.pops, {(k.split("_", 1)[0], k.split("_", 1)[1]): d[k] for k in d.files})
        return

    pop_of = _fig3.load_population_map()
    related = _fig3.load_related_samples()
    ibd = _fig3.read_ibd_segments()
    ibd = ibd[~ibd["id1"].isin(related) & ~ibd["id2"].isin(related) & (ibd["id1"] != ibd["id2"])]

    TABLE_DIR.mkdir(parents=True, exist_ok=True)
    FIGURE_DIR.mkdir(parents=True, exist_ok=True)
    summary, ratios_for_plot = [], {}

    for pop in args.pops:
        unrelated = set(pop_of[pop_of == pop].index) - related
        lai = load_lai(pop, unrelated)
        inds = sorted(set(lai["name"]))
        haps = [f"{i}_{h}" for i in inds for h in (1, 2)]
        H = len(haps)
        idx = {h: k for k, h in enumerate(haps)}
        print(f"{pop}: {len(inds)} unrelated individuals with LAI calls, {H} haplotypes", flush=True)

        shared = np.zeros((len(ANCESTRIES), H, H))      # bp where both carry ancestry a
        on_ibd = np.zeros((len(ANCESTRIES), H, H))      # of which covered by IBD between them

        sub = ibd[ibd["id1"].isin(inds) & ibd["id2"].isin(inds)]
        for chrom in CHROMOSOMES:
            lc = lai[lai["chrom"] == chrom]
            if lc.empty:
                continue
            edges, codes = ancestry_grid(lc, haps)
            width = np.diff(edges).astype(float)
            for a in range(len(ANCESTRIES)):
                B = (codes == a).astype(float)
                shared[a] += (B * width) @ B.T

            for r in sub[sub["chrom"] == chrom].itertuples(index=False):
                i, j = idx[f"{r.id1}_{r.hap1}"], idx[f"{r.id2}_{r.hap2}"]
                lo, hi = np.searchsorted(edges, [r.start_bp, r.end_bp])
                lo, hi = max(lo - 1, 0), min(hi, len(width))
                seg_w = np.clip(np.minimum(edges[lo + 1:hi + 1], r.end_bp)
                                - np.maximum(edges[lo:hi], r.start_bp), 0, None)
                same = codes[i, lo:hi] == codes[j, lo:hi]
                for a in range(len(ANCESTRIES)):
                    bp = float((seg_w * (same & (codes[i, lo:hi] == a))).sum())
                    on_ibd[a, i, j] += bp
                    on_ibd[a, j, i] += bp

        # all haplotype pairs from different individuals
        pairs = [(i, j) for i, j in combinations(range(H), 2) if haps[i].rsplit("_", 1)[0] != haps[j].rsplit("_", 1)[0]]
        I = np.array([i for i, _ in pairs]); J = np.array([j for _, j in pairs])
        table = pd.DataFrame({"hap1": [haps[i] for i in I], "hap2": [haps[j] for j in J]})
        for a, anc in enumerate(ANCESTRIES):
            s, o = shared[a, I, J], on_ibd[a, I, J]
            table[f"{anc}_shared_bp"] = s.astype(np.int64)
            table[f"{anc}_ibd_bp"] = o.astype(np.int64)
            table[f"{anc}_ratio"] = np.where(s > 0, o / np.where(s > 0, s, 1), np.nan)

            has = s > 0
            summary.append({
                "pop": pop, "ancestry": anc, "n_individuals": len(inds),
                "n_haplotype_pairs": len(pairs),
                "pairs_sharing_ancestry": int(has.sum()),
                "pairs_no_ibd_on_shared": int((has & (o == 0)).sum()),
                "pct_no_ibd_of_sharing_pairs": round(100 * (has & (o == 0)).sum() / max(has.sum(), 1), 1),
                "pct_no_ibd_of_all_pairs": round(100 * (o == 0).sum() / len(pairs), 1),
                "mean_ratio_sharing_pairs": float(np.nanmean(table[f"{anc}_ratio"])) if has.any() else np.nan,
            })
            ratios_for_plot[(pop, anc)] = table.loc[table[f"{anc}_ratio"] > 0, f"{anc}_ratio"].to_numpy()
        table.to_csv(TABLE_DIR / f"ancestry_aware_ibd_{pop}.csv.gz", index=False)

    summary = pd.DataFrame(summary)
    summary.to_csv(TABLE_DIR / "ancestry_aware_ibd_summary.csv", index=False)
    print("\n" + summary.to_string(index=False))

    np.savez_compressed(TABLE_DIR / "ancestry_aware_ibd_ratios.npz",
                        **{f"{pop}_{anc}": v for (pop, anc), v in ratios_for_plot.items()})
    plot(args.pops, ratios_for_plot)
    print(f"wrote ancestry_aware_ibd_<POP>.csv.gz and ancestry_aware_ibd_summary.csv to {TABLE_DIR}")


def plot(pops, ratios_for_plot):
    """Histogram panel from the plotting cell of Chapter_1.ipynb (cell 91),
    unchanged apart from reading the ratios computed above and the font sizes.
    Only non-zero ratios are drawn, as in the original."""
    import seaborn as sns

    POPS = pops
    palette = {"EUR": "#1f77b4", "AMR": "#daa520", "AFR": "#7b3f99"}
    data = {}
    for pop in POPS:
        frames = []
        for anc in ANCESTRIES:
            df = pd.DataFrame({"ratio": ratios_for_plot.get((pop, anc), np.array([]))})
            df["ancestry"] = anc
            frames.append(df[["ancestry", "ratio"]])
        data[pop] = pd.concat(frames, ignore_index=True)

    # ---- 2x2 panel of histograms (non-zero only) ----
    fig, axes = plt.subplots(2, 2, figsize=(12, 8), sharex=True, sharey=True)
    for ax, pop in zip(axes.flat, POPS):
        nz = data[pop][data[pop]["ratio"] > 0]
        sns.histplot(
            data=nz,
            x="ratio",
            hue="ancestry",
            hue_order=ANCESTRIES,
            palette=palette,
            bins=60,
            log_scale=(True, False),
            element="step",
            stat="density",
            common_norm=False,
            ax=ax,
            legend=(ax is axes.flat[0]),
        )
        ax.set_title(pop, fontsize=18)
        ax.set_xlabel("")
        ax.set_ylabel("")
        ax.tick_params(axis="both", labelsize=14)
    leg = axes.flat[0].get_legend()
    if leg is not None:
        plt.setp(leg.get_texts(), fontsize=15)
        plt.setp(leg.get_title(), fontsize=15)

    fig.supxlabel("Proportion of shared ancestry between a haplotype pair that is IBD (bp)", fontsize=17)
    fig.supylabel("Density", fontsize=17)
    fig.tight_layout()
    out = FIGURE_DIR / "ancestry_aware_ibd_ratios.png"
    fig.savefig(out, dpi=200, bbox_inches="tight")
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
