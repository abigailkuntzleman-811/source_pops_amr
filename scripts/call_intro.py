import demes
import demesdraw
from IPython.display import display, SVG
from matplotlib import pyplot as plt
import msprime
import numpy as np
import pandas as pd
from upsetplot import UpSet
from intervaltree import IntervalTree
import stdpopsim
import argparse
import gzip
import tskit


# ID 0: Population(id=0, metadata={'description': 'Equilibrium/root population', 'name': 'Ancestral'})
# ID 1: Population(id=1, metadata={'description': 'Last common ancestor of archaic hominins.', 'name': 'Archaic'})
# ID 2: Population(id=2, metadata={'description': 'Neanderthal', 'name': 'NEA'})
# ID 3: Population(id=3, metadata={'description': 'Denisovan', 'name': 'DEN'})
# ID 4: Population(id=4, metadata={'description': 'Anatomically modern humans', 'name': 'AMH'})
# ID 5: Population(id=5, metadata={'description': 'Bottleneck out-of-Africa population', 'name': 'OOA'})
# ID 6: Population(id=6, metadata={'description': 'Yoruba in Ibadan, Nigeria. Representative of African ancestries.', 'name': 'YRI'})
# ID 7: Population(id=7, metadata={'description': 'Iberian populations in Spain (IBS). Representative of European ancestries.', 'name': 'IBS'})
# ID 8: Population(id=8, metadata={'description': 'Han Chinese in Beijing, China. Representative of East Asian ancestries.', 'name': 'CHB'})
# ID 9: Population(id=9, metadata={'description': 'Indigenous American in Mexico. Representative of Indigenous American ancestries.', 'name': 'MXB'})
# ID 10: Population(id=10, metadata={'description': 'Mexicans in Los Angeles. Representative of the Mexican population.', 'name': 'MXL'})

# pop_ids = {'IBS': 7, 'CHB': 8, 'MXB': 9, 'MXL': 10}

parser = argparse.ArgumentParser()
parser.add_argument("amr_pop")
parser.add_argument("amr_pop_num")
parser.add_argument("repnum")  # positional argument
parser.add_argument("demography")
parser.add_argument("admix_pop_num")
parser.add_argument("output_file")

args = parser.parse_args()

# Access the filename
amr_pop = str(args.amr_pop)
amr_pop_num = int(args.amr_pop_num)
repnum = int(args.repnum)
demo = str(args.demography)
admix_pop_num = int(args.admix_pop_num)
output = str(args.output_file)

tree_seq = tskit.load(f"/users/akuntzle/data/akuntzle/simulations_donor_eur/{demo}/tree_seqs/{amr_pop}_rep{repnum}.ts.gz")

def get_intro(dest, target_populations, ts):
    # Extract the migrations table.
    mig_table = ts.tables.migrations
    # Extract the denisovan migration table.
    intro_mig_table = mig_table[mig_table.dest == dest]
    # Extract target haplotypes in the form of list.
    tar_haplotype = []
    for pop in target_populations:
        tar_haplotype.extend(ts.get_samples(pop))
    tar_haplotype = sorted(tar_haplotype)
    # Extract the introgressing nodes for the given introgression time and tree coordinates.
    intro_nodes = intro_mig_table.node
    intro_left  = intro_mig_table.left
    intro_right = intro_mig_table.right
    intro_time = intro_mig_table.time
    # Intialize a dictionary.
    intro_dicc = {
        'hap_id': [], 't_l': [], 't_r': [], 'i_l': [], 'i_r': [],
        'o_len': [], 'o_perct': [], 'time': [],
    }
    # For every tree with the specified haplotypes...
    for i, tree in enumerate(ts.trees(tracked_samples=tar_haplotype)):
        # For every introgressing node...
        for j, node in enumerate(intro_nodes):
          # If the tree interval is covered, partially or fully, by the migrating node then it is introgressed...
          if ((tree.interval[0] < intro_right[j] and intro_left[j] < tree.interval[1])):
                # Extract the leaves of the introgressing node.
                leaves = set(tree.get_leaves(node))
                # If the introgressed tree has at least one of our target haplotypes...
                if tree.num_tracked_samples(node) > 0:
                    # For every target haplotype...
                    for target in tar_haplotype:
                        # If the haplotype is one of our leaves...
                        if target in leaves:
                            # Extract the coordinates for the tree.
                            t_l, t_r = tree.interval
                            # Calculate the amount and percentage of overlap.
                            o_len = min(t_r, intro_right[j]) - max(t_l, intro_left[j])
                            o_pert = o_len / (t_r - t_l)
                            # Update the dictionary.
                            intro_dicc['hap_id'].append(target)
                            # intro_dicc['t_idx'].append(i)
                            intro_dicc['t_l'].append(t_l)
                            intro_dicc['t_r'].append(t_r)
                            # intro_dicc['t_span'].append(t_r - t_l)
                            # intro_dicc['i_idx'].append(j)
                            intro_dicc['i_l'].append(intro_left[j])
                            intro_dicc['i_r'].append(intro_right[j])
                            # intro_dicc['i_span'].append(intro_right[j] - intro_left[j])
                            intro_dicc['o_len'].append(o_len)
                            intro_dicc['o_perct'].append(o_pert)
                            intro_dicc['time'].append(intro_time[j])
    # Convert the dictionary to a dataframe.
    intro_df = pd.DataFrame(intro_dicc)
    return intro_df

intro_df = get_intro(admix_pop_num, [amr_pop_num], tree_seq)
# pop_ids = {'IBS': 7, 'CHB': 8, 'MXB': 9, 'MXL': 10}
# for the donor population simulation, pop_ids are
# 0 {'description': 'Equilibrium/root population', 'name': 'Ancestral'}
# 1 {'description': 'Last common ancestor of archaic hominins.', 'name': 'Archaic'}
# 2 {'description': 'Neanderthal', 'name': 'NEA'}
# 3 {'description': 'Denisovan', 'name': 'DEN'}
# 4 {'description': 'Anatomically modern humans', 'name': 'AMH'}
# 5 {'description': 'Bottleneck out-of-Africa population', 'name': 'OOA'}
# 6 {'description': 'Yoruba in Ibadan, Nigeria. Representative of African ancestries.', 'name': 'YRI'}
# 7 {'description': 'Iberian populations in Spain (IBS). Representative of European ancestries.', 'name': 'IBS'}
# 8 {'description': 'Han Chinese in Beijing, China. Representative of East Asian ancestries.', 'name': 'CHB'}
# 9 {'description': 'Indigenous American in Mexico. Representative of Indigenous American ancestries.', 'name': 'MXB'}
# 10 {'description': 'Iberian donor population used for admixture into MXL', 'name': 'IBS_DONOR'}
# 11 {'description': 'Peruvians in Lima. Representative of the Peruvian population.', 'name': 'PEL'}

# Define a function to extract contiguous introgressed tracts.
def extract_introgressed_tracts(intro_df):
    # Initialize a dictionary.
    tract_dicc = {
        "hap_id": [],
        "start": [],
        "end": [],
        "length": [],
        "time": [],
    }
    # For every haplotype with introgressed material.
    for hap_id in intro_df.hap_id.unique():
        # Subset the dataframe.
        hap_df = intro_df[intro_df.hap_id == hap_id]
        # Initialize the tract coordinates.
        start, end = None, None
        # For all trees that overlap introgressed material.
        for _, row in hap_df.iterrows():
            # Extract the coordinates.
            left = row.t_l if row.o_perct == 1 else max(row.t_l, row.i_l)
            right = row.t_r if row.o_perct == 1 else min(row.t_r, row.i_r)
            time = row.time
            # If this is the start of a new segment.
            if start is None and end is None:
                # Initialize the start of the segment.
                start = left
                end = right
            # Else-if, this is a continuation of an existing segment.
            elif end == left:
                # Adjust the end coordinate.
                end = right
            # Else the current tract is complete.
            else:
                # Update the dictionary.
                tract_dicc["hap_id"].append(hap_id)
                tract_dicc["time"].append(time)
                tract_dicc["start"].append(start)
                tract_dicc["end"].append(end)
                tract_dicc["length"].append(end - start)
                # Initialize the new tract coordinates.
                start, end = left, right
        # Add the last tract if it exists.
        if start is not None and end is not None:
            tract_dicc["hap_id"].append(hap_id)
            tract_dicc["time"].append(time)
            tract_dicc["start"].append(start)
            tract_dicc["end"].append(end)
            tract_dicc["length"].append(end - start)
    # Convert the dictionary to a dataframe.
    tract_df = pd.DataFrame(tract_dicc)
    return tract_df

# intro_df.to_csv(f'{output}/neanderthal_csvs_stitched/rep{repnum}_unstitched.csv.gz', index=False, compression='gzip')
# admix_df.to_csv(f'{output}/admix_csvs_stitched/rep{repnum}_unstitched_admix.csv.gz', index=False, compression='gzip')

# intro_df = pd.read_csv(f'{output}/neanderthal_csvs_stitched/rep{repnum}_unstitched.csv.gz', index_col=False, compression='gzip')
# admix_df = pd.read_csv(f'{output}/admix_csvs_stitched/rep{repnum}_unstitched_admix.csv.gz', index_col=False, compression='gzip')

intro_df_stitched = extract_introgressed_tracts(intro_df)

intro_df_stitched.to_csv(f'{output}/admix_csvs/{amr_pop}_{admix_pop_num}tracts_rep{repnum}_stitched.csv.gz', index=False, compression='gzip')