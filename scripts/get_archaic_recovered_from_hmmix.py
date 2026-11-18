import matplotlib.pyplot as plt
import pandas as pd
import pyranges as pr
import argparse
import numpy as np


parser = argparse.ArgumentParser()

parser.add_argument('filename')
parser.add_argument('nean_type')

args = parser.parse_args()

filename = args.filename
nean_type = args.nean_type

rng = np.random.default_rng(42)


hmmix_nean = pd.read_csv(filename)

clm_nean = hmmix_nean[hmmix_nean['pop'] == 'CLM']
pur_nean = hmmix_nean[hmmix_nean['pop'] == 'PUR']
pel_nean = hmmix_nean[hmmix_nean['pop'] == 'PEL']
mxl_nean = hmmix_nean[hmmix_nean['pop'] == 'MXL']

clm_inds = clm_nean['name'].unique()
pur_inds = pur_nean['name'].unique()
pel_inds = pel_nean['name'].unique()
mxl_inds = mxl_nean['name'].unique()

min_number = min(len(clm_inds), len(pur_inds), len(pel_inds), len(mxl_inds))

def get_archaic_recovered_from_hmmix(pop_df, inds):
    chroms = [f'chr{i}' for i in range(1, 23)]
    df = pop_df[
        pop_df['name'].isin(inds)
        & pop_df['chrom'].isin(chroms)
        & pop_df['haplotype'].isin(['hap1', 'hap2'])
    ].rename(columns={'chrom': 'Chromosome', 'start': 'Start', 'end': 'End'})

    gr = pr.PyRanges(df[['Chromosome', 'Start', 'End']])
    merged = gr.merge()

    return merged


pops = ['CLM', 'PUR', 'PEL', 'MXL']
pop_ids = [clm_inds, pur_inds, pel_inds, mxl_inds]
pop_dfs = [clm_nean, pur_nean, pel_nean, mxl_nean]
for id_set in range(len(pop_ids)):
    x = []
    y = []
    for i in range(1, len(pop_ids[id_set]) + 1):
        x.append(i)
        y.append(get_archaic_recovered_from_hmmix(pop_dfs[id_set], pop_ids[id_set][:i]).lengths().sum())
    with open(f'/users/akuntzle/IBD_archaic_shared/data/recovered_archaic_by_admix_prop/{pops[id_set]}_total_{nean_type}_recovered_from_hmmix_pyranges.txt', 'w') as f:
        f.write("x: " + ", ".join(map(str, x)) + "\n")
        f.write("y: " + ", ".join(map(str, y)) + "\n")
    sampled = rng.choice(np.sort(pop_ids[id_set]), size=min_number, replace=False)
    merged = get_archaic_recovered_from_hmmix(pop_dfs[id_set], sampled)
    pd.DataFrame({'chrom': merged.Chromosome, 'start': merged.Start, 'end': merged.End}).to_csv(
        f'/users/akuntzle/IBD_archaic_shared/data/recovered_archaic_by_admix_prop/{pops[id_set]}_n{min_number}_recovered_{nean_type}.bed',
        sep='\t', header=False, index=False)
    with open(f'/users/akuntzle/IBD_archaic_shared/data/recovered_archaic_by_admix_prop/{pops[id_set]}_n{min_number}_recovered_{nean_type}_samples.txt', 'w') as f:
        f.write("\n".join(sorted(sampled)) + "\n")