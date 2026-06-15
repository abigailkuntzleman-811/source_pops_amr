import matplotlib.pyplot as plt
import pandas as pd
import argparse
import pyranges as pr

parser = argparse.ArgumentParser()
parser.add_argument("amr_pop")  # positional argument
parser.add_argument("continental_ancestry")

args = parser.parse_args()

# Access the filename
amr_pop = str(args.amr_pop)
continental_ancestry = str(args.continental_ancestry)

# /oscar/data/ehuertas/data/sharing/jaz/1kg_bed_files/CLM_lai.bed
def get_pop_dict_lai():
    pop_dict = {}
    # for pop in ['CLM', 'PEL','PUR','MXL']:
    pop_ids = set()
    with open(f'/oscar/data/ehuertas/data/sharing/jaz/1kg_lai_global/{amr_pop}.bedlist', 'r') as file:
        for line in file.readlines():
            line = line.strip().split('_')
            pop_ids.add(line[3])
    pop_dict[amr_pop] = list(pop_ids)
    return(pop_dict)
pop_dict = get_pop_dict_lai()

# children_to_filter = pd.read_csv('/users/akuntzle/IBD_archaic_shared/data/children_to_remove.txt', header=None, names=['sampleID'])
# ibd_df = []
# for i in range(1, 23):
#     ibd_chr = pd.read_csv(f'/users/akuntzle/data/akuntzle/1kgp_highcov_ibd_calls/kgp_chr{i}_segments.ibd.gz',sep='\t',header=None)
#     ibd_chr.columns = ['ID1', 'hap1', 'ID2', 'hap2', 'chr', 'start_bp', 'end_bp', 'CM']
#     children_set = set(children_to_filter['sampleID'])
#     ibd_chr = ibd_chr[
#         (~ibd_chr['ID1'].isin(children_set)) &
#         (~ibd_chr['ID2'].isin(children_set))
#     ].copy()
#     ibd_df.append(ibd_chr)

# all_ibd = pd.concat(ibd_df, ignore_index=True)

# all_ibd[all_ibd['ID1'].isin(pop_dict['MXL']) & all_ibd['ID2'].isin(pop_dict['MXL'])]

# def get_props(pop_df, all_ibd,pop):
#     prop_dict = {}

#     pop_df['seg_len'] = pop_df['end'] - pop_df['start']
    
#     ind_ancestry = (
#         pop_df
#         .groupby(['name', 'ancestry'])['seg_len']
#         .sum()
#     )
#     wide = ind_ancestry.unstack(fill_value=0)
#     prop = wide.div(wide.sum(axis=1), axis=0)
#     # prop = wide.div(6_000_000_000)
#     wide_prop = wide.add_suffix('_bp').join(
#         prop.add_suffix('_prop')
#     )

    
#     clm_set = set(pop_dict[pop])
#     clm_ibd = all_ibd[
#         all_ibd['ID1'].isin(clm_set) &
#         all_ibd['ID2'].isin(clm_set) &
#         (all_ibd['ID1'] != all_ibd['ID2'])
#     ]
    
#     s1 = clm_ibd.groupby('ID1')['CM'].sum()
#     s2 = clm_ibd.groupby('ID2')['CM'].sum()
    
#     ibd_series = s1.add(s2, fill_value=0)
#     ibd_series.name = 'total_inpop_ibd'
    
#     wide_prop = wide_prop.join(ibd_series)
#     return wide_prop

df = pd.read_csv(f'/oscar/data/ehuertas/data/sharing/jaz/1kg_bed_files/{amr_pop}_lai.bed',sep='\t')

ids = df['name'].unique()


def get_pop_archaic_recovered(ind_ids, lai_pop):
    recovered_seq = 0
    #for each chromosomes/haplotype, I want to merge the sequence over all the individuals
    for i in range(1, 23):
        overlaps = []
        for hap in [1,2]:
            for ind in ind_ids:
                chr_bed = pd.read_csv(
                    f'/oscar/data/ehuertas/data/sharing/jaz/1kg_arc_bed_files/chr{i}_{ind}_{hap}.bed',
                    sep='\t'
                )
                # only introgressed regions
                intro_bed = chr_bed[chr_bed['ancestry'] == 'introgressed']

                
                eur_bed = pd.read_csv(f'/oscar/data/ehuertas/data/sharing/jaz/1kg_bed_files/chr{i}_{ind}_{hap}.bed',
                                      sep = '\t'
                                     )
                #continental regions
                eur_bed = eur_bed[eur_bed['ancestry'] == lai_pop]

                intro_bed = intro_bed.rename(columns={
                    "chrom": "Chromosome",
                    "start": "Start",
                    "end": "End"
                })

                eur_bed = eur_bed.rename(columns={
                    "chrom": "Chromosome",
                    "start": "Start",
                    "end": "End"
                })
                intro_pr = pr.PyRanges(intro_bed)
                eur_pr   = pr.PyRanges(eur_bed)

                overlaps.append(intro_pr.intersect(eur_pr))
        overlap = pr.concat(overlaps)  
        #now I merge all the sequence that is archaic and european in the population
        overlap = overlap.merge()
        recovered_seq += (overlap.End - overlap.Start).sum()
        
    return recovered_seq


# amr_prop = pd.concat([clm_prop, pur_prop, pel_prop, mxl_prop])

# def return_bins(pop):
#     bins = []
#     for i in [0.2, 0.3, 0.4, 0.5, 0.6, 0.7]:
#         bins.append(list(amr_prop[(amr_prop[f'{pop}_prop'] >= i) & (amr_prop[f'{pop}_prop'] < i+ .1)].index))
#     return bins

# labels = ['20-30%', '30-40%', '40-50%', '50-60%', '60-70%', '70-80%']
# amr_bins = return_bins('AMR')
# eur_bins = return_bins('EUR')

# amr_ids = [clm_ids, pur_ids, pel_ids, mxl_ids]
# amr_labels = ['CLM','PUR','PEL','MXL']
# for id_set, label in zip(amr_ids, amr_labels):
x = []
y = []
for i in range(1, len(ids)+ 1):
    subset = ids[:i]
    x.append(i)
    y.append(get_pop_archaic_recovered(subset, continental_ancestry))
with open(f'/users/akuntzle/IBD_archaic_shared/data/recovered_archaic_by_admix_prop/{amr_pop}_only_{continental_ancestry}_archaic_recovered.txt', 'w') as f:
    f.write("x: " + ", ".join(map(str, x)) + "\n")
    f.write("y: " + ", ".join(map(str, y)) + "\n")