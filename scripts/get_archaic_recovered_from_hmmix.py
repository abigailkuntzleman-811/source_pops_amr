import matplotlib.pyplot as plt
import pandas as pd
import pyranges as pr
import argparse

parser = argparse.ArgumentParser()

parser.add_argument('filename')
parser.add_argument('nean_type')

args = parser.parse_args()

filename = args.filename
nean_type = args.nean_type


hmmix_nean = pd.read_csv(filename)

clm_nean = hmmix_nean[hmmix_nean['pop'] == 'CLM']
pur_nean = hmmix_nean[hmmix_nean['pop'] == 'PUR']
pel_nean = hmmix_nean[hmmix_nean['pop'] == 'PEL']
mxl_nean = hmmix_nean[hmmix_nean['pop'] == 'MXL']

clm_inds = clm_nean['name'].unique()
pur_inds = pur_nean['name'].unique()
pel_inds = pel_nean['name'].unique()
mxl_inds = mxl_nean['name'].unique()

# def get_archaic_recovered_from_hmmix(pop_df, inds):
#     total_seq = 0
#     for i in range(1, 23):
#         total_tree = IntervalTree()
#         for ind in inds:
#             for hap in [1, 2]:
#                 these_inds = pop_df[(pop_df['name'] == ind) & (pop_df['chrom'] == f'chr{i}') & (pop_df['haplotype'] == f'hap{hap}')]
#                 for _, row in these_inds.iterrows():
#                     total_tree[row['start']:row['end']] = True
#         total_tree.merge_overlaps()
#         total_seq += sum(interval.end - interval.begin for interval in total_tree)

#     return total_seq

# import pyranges as pr

def get_archaic_recovered_from_hmmix(pop_df, inds):
    chroms = [f'chr{i}' for i in range(1, 23)]
    df = pop_df[
        pop_df['name'].isin(inds)
        & pop_df['chrom'].isin(chroms)
        & pop_df['haplotype'].isin(['hap1', 'hap2'])
    ].rename(columns={'chrom': 'Chromosome', 'start': 'Start', 'end': 'End'})

    gr = pr.PyRanges(df[['Chromosome', 'Start', 'End']])
    merged = gr.merge()

    return merged.lengths().sum()


pops = ['CLM', 'PUR', 'PEL', 'MXL']
pop_ids = [clm_inds, pur_inds, pel_inds, mxl_inds]
pop_dfs = [clm_nean, pur_nean, pel_nean, mxl_nean]
for id_set in range(len(pop_ids)):
    x = []
    y = []
    for i in range(1, len(pop_ids[id_set]) + 1):
        x.append(i)
        y.append(get_archaic_recovered_from_hmmix(pop_dfs[id_set], pop_ids[id_set][:i]))
    with open(f'/users/akuntzle/IBD_archaic_shared/data/recovered_archaic_by_admix_prop/{pops[id_set]}_total_{nean_type}_recovered_from_hmmix_pyranges.txt', 'w') as f:
        f.write("x: " + ", ".join(map(str, x)) + "\n")
        f.write("y: " + ", ".join(map(str, y)) + "\n")