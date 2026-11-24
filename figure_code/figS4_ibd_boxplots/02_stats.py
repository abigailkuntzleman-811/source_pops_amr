#!/usr/bin/env python3
"""Figure S4 statistics: each admixed population vs IBS, bootstrapping individuals.

Haplotype pairs are not independent (every individual appears in hundreds of
pairs), so tests that treat pairs as observations (Kruskal-Wallis, Dunn,
Mann-Whitney on pairs) give p values that are far too small. Instead:

  1. statistic per population, from the per-pair table of
     figure_code/fig3/01_within_population_ibd.py (masked IBD, KING set):
       mean_all      mean IBD per haplotype pair, zero pairs included
       mean_nonzero  mean IBD per haplotype pair among pairs with any IBD
  2. bootstrap: resample the population's individuals with replacement and
     recompute the statistic over the pairs among the resampled individuals
     (an individual drawn k times counts k times; pairs of an individual with
     itself are left out); --boot times, independently per population;
  3. difference vs IBS: 95% percentile interval of the bootstrap differences,
     two-sided bootstrap p value, Holm correction over the four populations
     (separately for each statistic).

With --all-pairs, every pair of the five populations is compared instead
(10 comparisons, Holm over all 10).

Usage:
    python figure_code/figS4_ibd_boxplots/02_stats.py
    python figure_code/figS4_ibd_boxplots/02_stats.py --all-pairs

Output: $TABLE_DIR/IBD_boxplots_excludes0s_stats.csv
        $TABLE_DIR/IBD_boxplots_excludes0s_stats_all_pairs.csv   (--all-pairs)
"""

import argparse
import sys
from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd

_root = next(p for p in Path(__file__).resolve().parents if (p / "config.py").exists())
sys.path.insert(0, str(_root))
from config import IBD_DIR, TABLE_DIR  # noqa: E402

REFERENCE = "IBS"
ADMIXED = ["CLM", "MXL", "PEL", "PUR"]
STATS = ["mean_all", "mean_nonzero"]


def pair_matrices(pairs):
    """Per individual pair: summed cM over its haplotype pairs (S), number of
    haplotype pairs (N), and number with any IBD (K); n x n, symmetric."""
    h = pairs["pair"].str.split("_", n=1, expand=True)
    a = h[0].str.rsplit("-", n=1).str[0]
    b = h[1].str.rsplit("-", n=1).str[0]
    inds = sorted(set(a) | set(b))
    idx = {x: i for i, x in enumerate(inds)}
    i, j = a.map(idx).to_numpy(), b.map(idx).to_numpy()
    n = len(inds)
    S, N, K = (np.zeros((n, n)) for _ in range(3))
    cm = pairs["cM"].to_numpy()
    for M, v in ((S, cm), (N, np.ones_like(cm)), (K, (cm > 0).astype(float))):
        np.add.at(M, (i, j), v)
        np.add.at(M, (j, i), v)
    return S, N, K


def statistic(S, N, K, c):
    """Both statistics for individual weights c (counts in a bootstrap draw)."""
    W = np.outer(c, c)
    np.fill_diagonal(W, 0)
    s, n, k = (W * S).sum(), (W * N).sum(), (W * K).sum()
    return {"mean_all": s / n, "mean_nonzero": s / k if k else np.nan}


def holm(p):
    p = np.asarray(p, float)
    order = np.argsort(p)
    adj = np.empty_like(p)
    running = 0.0
    for rank, k in enumerate(order):
        running = max(running, min(1.0, (len(p) - rank) * p[k]))
        adj[k] = running
    return adj


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--boot", type=int, default=10000)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--all-pairs", action="store_true",
                    help="compare every pair of populations, not just each admixed population vs IBS")
    args = ap.parse_args()
    rng = np.random.default_rng(args.seed)

    path = IBD_DIR / "within_population_pairs.csv.gz"
    if not path.exists():
        sys.exit(f"Missing {path}\nRun: python figure_code/fig3/01_within_population_ibd.py")
    pairs = pd.read_csv(path)

    observed, boots, n_ind = {}, {}, {}
    for pop in [REFERENCE] + ADMIXED:
        S, N, K = pair_matrices(pairs[pairs["pop"] == pop])
        n = len(S)
        n_ind[pop] = n
        observed[pop] = statistic(S, N, K, np.ones(n))
        draws = rng.multinomial(n, np.full(n, 1 / n), size=args.boot)
        boots[pop] = pd.DataFrame([statistic(S, N, K, c) for c in draws])
        print(f"{pop}: {n} individuals, mean_all {observed[pop]['mean_all']:.4f}, "
              f"mean_nonzero {observed[pop]['mean_nonzero']:.4f}", flush=True)

    order = [REFERENCE] + ADMIXED
    comparisons = list(combinations(order, 2)) if args.all_pairs else [(REFERENCE, p) for p in ADMIXED]
    rows = []
    for stat in STATS:
        block = []
        for ref, pop in comparisons:
            d = (boots[pop][stat] - boots[ref][stat]).dropna().to_numpy()
            p = min(1.0, 2 * min((d <= 0).mean(), (d >= 0).mean()))
            p = max(p, 1 / (len(d) + 1))
            lo, hi = np.percentile(d, [2.5, 97.5])
            pop_lo, pop_hi = np.percentile(boots[pop][stat].dropna(), [2.5, 97.5])
            block.append({
                "statistic": stat, "population": pop, "compared_with": ref,
                "n_individuals": n_ind[pop],
                "value": observed[pop][stat], "value_ci_low": pop_lo, "value_ci_high": pop_hi,
                "other_value": observed[ref][stat],
                "difference": observed[pop][stat] - observed[ref][stat],
                "diff_ci_low": lo, "diff_ci_high": hi, "p": p,
            })
        for r, adj in zip(block, holm([r["p"] for r in block])):
            r["p_holm"] = adj
        rows += block

    out = pd.DataFrame(rows)
    TABLE_DIR.mkdir(parents=True, exist_ok=True)
    path = TABLE_DIR / ("IBD_boxplots_excludes0s_stats_all_pairs.csv" if args.all_pairs
                        else "IBD_boxplots_excludes0s_stats.csv")
    out.to_csv(path, index=False)
    print("\n" + out.round(4).to_string(index=False))
    print(f"\n{args.boot} bootstrap draws of individuals; p < {1 / (args.boot + 1):.0e} is the floor")
    print(f"wrote {path}")


if __name__ == "__main__":
    main()
