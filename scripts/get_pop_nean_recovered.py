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

x = []
y = []
for i in range(1, len(ids)+ 1):
    subset = ids[:i]
    x.append(i)
    y.append(get_pop_archaic_recovered(subset, continental_ancestry))
with open(f'/users/akuntzle/IBD_archaic_shared/data/recovered_archaic_by_admix_prop/{amr_pop}_only_{continental_ancestry}_archaic_recovered.txt', 'w') as f:
    f.write("x: " + ", ".join(map(str, x)) + "\n")
    f.write("y: " + ", ".join(map(str, y)) + "\n")