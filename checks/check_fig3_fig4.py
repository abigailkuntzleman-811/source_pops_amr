#!/usr/bin/env python3
"""Check Figure 3 / S4 and Figure 4 / S6 / S8 against Chapter_1.ipynb.

Runs the notebook's own code (cells 37-41, 52 and the regression of cell 62),
copied unchanged except for the input paths, which point at the inputs we
switched to on purpose: the pileup-masked IBD files, the KING unrelated list
(RELATED_EXCLUDE) and the joint-caller <POP>_lai.bed files. Then compares with
what the repository scripts wrote:

  Fig 3 / S4   IBD_DIR/within_population_pairs.csv.gz, within_population_summary.csv
               (figure_code/fig3/01_within_population_ibd.py)
  Fig 4        IBD_DIR/individual_ibd.csv            (figure_code/fig4/01_individual_ibd.py)
               DATA_PROCESSED/ancestry_proportions/<POP>_individual.csv
                                                    (table1/00_ancestry_proportions.py)
               TABLE_DIR/<ANC>_prop_admix_by_ibd_stats.csv
                                                    (figure_code/fig4/02_plot_admix_vs_ibd.py)

The notebook loops over every IBD segment with iterrows, so run it as a job:
    sbatch -J check_fig3_fig4 run.script checks/check_fig3_fig4.py
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import linregress, spearmanr

_root = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(_root))
from config import (AMR_POPS, DATA_PROCESSED, IBD_DIR, LAI_BED_DIR,  # noqa: E402
                    POP_PANEL, RELATED_EXCLUDE, TABLE_DIR)

if RELATED_EXCLUDE is None:
    sys.exit("Set RELATED_EXCLUDE (the KING unrelated_toberemoved.txt file) in config.local.sh")

PED = str(POP_PANEL)
KING = str(RELATED_EXCLUDE)
IBD = str(IBD_DIR / "kgp_chr{i}_segments.masked.ibd.gz")

sub_pops = ['CHB', 'JPT', 'CHS', 'CDX', 'KHV',
            'CEU', 'TSI', 'GBR', 'FIN', 'IBS',
            'GIH', 'PJL', 'BEB', 'STU', 'ITU',
            'MXL', 'PEL', 'PUR', 'CLM']

# ---------------------------------------------------------------- cell 37
'''Function takes a the metadata file and returns 2 dictionaries: one mapping individuals to
    their population, and the other mapping populations to the list of individuals
'''
def get_pop_dict(df, pop_col, ind_col):
    population_dict = {}
    individual_dict = {}
    for _, row in df.iterrows():
        population = row[pop_col]
        individual = row[ind_col]
        if population not in population_dict:
            population_dict[population] = []
        population_dict[population].append(individual)
        individual_dict[individual] = population
    return population_dict, individual_dict
fam = pd.read_csv(PED, sep=' ')
pop_dict, ind_dict = get_pop_dict(fam, 'Population', 'SampleID')

# ---------------------------------------------------------------- cell 38
#Grabbing the maximal set of unrelated individuals
children_to_filter = pd.read_csv(
    KING,
    sep=r'\s+',
    header=None,
    names=['FID', 'sampleID']
)
children_set = set(children_to_filter['sampleID'])

# ---------------------------------------------------------------- cell 39
# all autosomal IBD segments, with related individuals removed
ibd_df = []
for i in range(1, 23):
    ibd_chr = pd.read_csv(
        IBD.format(i=i),
        sep='\t', header=None,
        names=['ID1', 'hap1', 'ID2', 'hap2', 'chr', 'start_bp', 'end_bp', 'CM']
    )
    ibd_chr = ibd_chr[
        ~ibd_chr['ID1'].isin(children_set) &
        ~ibd_chr['ID2'].isin(children_set)
    ].copy()
    ibd_df.append(ibd_chr)

all_ibd = pd.concat(ibd_df, ignore_index=True)
print(f'{len(all_ibd):,} segments after removing related individuals')

# ---------------------------------------------------------------- cell 40
# Grabbing the individuals who are not in the related set
pop_dict_filtered = {
    k: [v for v in vs if v not in children_set]
    for k, vs in pop_dict.items()
}
''' Initializing a dictionary which will map pairs of individuals (unrelated)
    to their shared IBD for within populations
'''
combo_dict = {}
#only grabbing the individuals that are unrelated
for key1, value1 in pop_dict_filtered.items():
    within_key = f'{key1}_{key1}'
    combo_dict[within_key] = {}
    #this grabs all the haplotypes and pairs for all of them
    for i in range(len(value1)):
        for j in range(i + 1, len(value1)):
            for hap1 in ['1', '2']:
                for hap2 in ['1', '2']:
                    id1 = value1[i] + '-' + hap1
                    id2 = value1[j] + '-' + hap2
                    combo_dict[within_key]['_'.join(sorted([id1, id2]))] = 0

#Filling in combo_dict with IBD shared between haplotypes
for i in range(1, 23):
    ibd = pd.read_csv(
        IBD.format(i=i),
        sep='\t', header=None, compression='gzip',
        names=['id1', 'hap1', 'id2', 'hap2', 'chr', 'start_bp', 'end_bp', 'CM']
    )
    print(f'chr{i}', flush=True)
    for _, row in ibd.iterrows():
        #skipping individuals in the related set
        if row['id1'] in children_set or row['id2'] in children_set:
            continue
        pop1, pop2 = ind_dict[row['id1']], ind_dict[row['id2']]
        #combo_dict is only for within-population pairs of IBD, so exclude pairs of pops where they are
        #not equal to each other
        if pop1 != pop2:
            continue
        #creating the key using the same format as in combo_dict
        ind_key = '_'.join(sorted([row['id1'] + '-' + str(row['hap1']),
                                   row['id2'] + '-' + str(row['hap2'])]))
        #creating the population key, and it's all within-pop so just POP1_POP1
        pop_key = f'{pop1}_{pop1}'
        if ind_key in combo_dict.get(pop_key, {}):
            combo_dict[pop_key][ind_key] += float(row['CM'])

# ---------------------------------------------------------------- cell 41
ibd_diag = []
for pop in sub_pops:
    ibd_diag.append(sum(combo_dict[f'{pop}_{pop}'].values())/len(combo_dict[f'{pop}_{pop}'].values()))


# ---------------------------------------------------------------- cell 52
def get_props(pop_df, all_ibd, pop):
    pop_df = pop_df.copy()
    pop_df['seg_len'] = pop_df['end'] - pop_df['start']

    # per-individual ancestry proportions
    ind_ancestry = pop_df.groupby(['name', 'ancestry'])['seg_len'].sum()
    wide = ind_ancestry.unstack(fill_value=0)
    prop = wide.div(wide.sum(axis=1), axis=0)
    wide_prop = wide.add_suffix('_bp').join(prop.add_suffix('_prop'))

    # unrelated individuals only — related ones are already absent from
    # all_ibd, so counting them in the denominator would understate the mean
    pop_set = set(pop_dict_filtered[pop])
    n_pop = len(pop_set)

    pop_ibd = all_ibd[
        all_ibd['ID1'].isin(pop_set) &
        all_ibd['ID2'].isin(pop_set) &
        (all_ibd['ID1'] != all_ibd['ID2'])
    ]

    s1 = pop_ibd.groupby('ID1')['CM'].sum()
    s2 = pop_ibd.groupby('ID2')['CM'].sum()

    # total cM each individual shares with the rest of the population,
    # averaged over the (n - 1) possible partners
    ibd_series = s1.add(s2, fill_value=0) / (n_pop - 1)
    ibd_series.name = 'avg_inpop_ibd'

    wide_prop = wide_prop.loc[wide_prop.index.isin(pop_set)]

    wide_prop = wide_prop.join(ibd_series)
    wide_prop['avg_inpop_ibd'] = wide_prop['avg_inpop_ibd'].fillna(0)

    return wide_prop


# ================================================================ comparisons
problems = []


def report(name, ok, detail=""):
    print(f"  {'OK  ' if ok else 'DIFF'} {name}{(': ' + detail) if detail else ''}")
    if not ok:
        problems.append(name)


print("\n=== Figure 3 / S4: per-pair IBD (combo_dict) and population means (ibd_diag)")
pairs = pd.read_csv(IBD_DIR / "within_population_pairs.csv.gz")
summary = pd.read_csv(IBD_DIR / "within_population_summary.csv").set_index("pop")
for pop, nb_mean in zip(sub_pops, ibd_diag):
    nb = pd.Series(combo_dict[f'{pop}_{pop}'], dtype=float)
    repo = pairs[pairs["pop"] == pop].set_index("pair")["cM"]
    same_keys = set(nb.index) == set(repo.index)
    maxdiff = float((nb - repo.reindex(nb.index)).abs().max()) if same_keys else np.nan
    report(f"{pop} pairs", same_keys and maxdiff < 1e-6,
           f"{len(nb):,} notebook vs {len(repo):,} repo pairs, max |diff| {maxdiff:.2e} cM")
    rm = summary.loc[pop, "mean_cM"]
    report(f"{pop} mean", abs(nb_mean - rm) < 1e-9,
           f"notebook {nb_mean:.6f} vs repo {rm:.6f} cM; n unrelated {len(pop_dict_filtered[pop])} "
           f"vs {summary.loc[pop, 'n_individuals']}")

print("\n=== Figure 4 / S6 / S8: get_props vs repo IBD_bar_i and ancestry proportions")
ind_ibd = pd.read_csv(IBD_DIR / "individual_ibd.csv")
nb_props = {}
for pop in AMR_POPS:
    pop_df = pd.read_csv(LAI_BED_DIR / f"{pop}_lai.bed", sep='\t')
    nb = get_props(pop_df, all_ibd, pop)
    nb_props[pop] = nb
    repo_ibd = ind_ibd[ind_ibd["pop"] == pop].set_index("individual")["avg_inpop_ibd"]
    d = (nb["avg_inpop_ibd"] - repo_ibd.reindex(nb.index)).abs()
    report(f"{pop} IBD_bar_i", d.notna().all() and d.max() < 1e-9,
           f"{len(nb)} individuals, max |diff| {d.max():.2e} cM")
    repo_p = pd.read_csv(DATA_PROCESSED / "ancestry_proportions" / f"{pop}_individual.csv", index_col=0)
    for col in [c for c in nb.columns if c.endswith("_prop")]:
        d = (nb[col] - repo_p[col].reindex(nb.index)).abs()
        report(f"{pop} {col}", d.notna().all() and d.max() < 1e-12, f"max |diff| {d.max():.2e}")

print("\n=== Figure 4 / S6 / S8: regression statistics (cell 62) vs repo stats tables")
for anc in ["EUR", "AMR", "AFR"]:
    stats_path = TABLE_DIR / f"{anc}_prop_admix_by_ibd_stats.csv"
    if not stats_path.exists():
        print(f"  (no {stats_path.name}; run fig4/02_plot_admix_vs_ibd.py --ancestry {anc})")
        continue
    repo = pd.read_csv(stats_path).set_index("pop")
    for pop in AMR_POPS:
        df = nb_props[pop]
        if f"{anc}_prop" not in df.columns:
            continue
        x = df[f"{anc}_prop"].values
        y = df["avg_inpop_ibd"].values
        mask = ~(np.isnan(x) | np.isnan(y))
        x, y = x[mask], y[mask]
        result = linregress(x, y)
        rho, p_spearman = spearmanr(x, y)
        r = repo.loc[pop]
        ok = (len(x) == r["n"] and np.isclose(result.slope, r["slope"])
              and np.isclose(rho, r["spearman_rho"]) and np.isclose(p_spearman, r["spearman_p"]))
        report(f"{anc} {pop}", ok,
               f"n {len(x)} vs {int(r['n'])}, slope {result.slope:.4f} vs {r['slope']:.4f}, "
               f"rho {rho:.4f} vs {r['spearman_rho']:.4f}, p {p_spearman:.3e} vs {r['spearman_p']:.3e}")

print("\nALL IDENTICAL" if not problems else f"\n{len(problems)} DIFFERENCES: {', '.join(problems)}")
