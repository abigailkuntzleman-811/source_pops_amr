import pandas as pd
import pyranges as pr
import argparse

# standard model, no donor pop: pop_ids = {'IBS': 7, 'CHB': 8, 'MXB': 9, 'MXL': 10}

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

parser = argparse.ArgumentParser()
parser.add_argument("repnum")
parser.add_argument("demo")
parser.add_argument("pop")
parser.add_argument("eur_admix_pop_num")

args = parser.parse_args()

# Access the filename
repnum = int(args.repnum)
demo = str(args.demo)
pop = str(args.pop)
eur_admix_pop_num = int(args.eur_admix_pop_num)

arc_df = pd.read_csv(f'/users/akuntzle/data/akuntzle/simulations_donor_eur/{demo}/admix_csvs/{pop}_2tracts_rep{repnum}_stitched.csv.gz')
eur_df = pd.read_csv(f'/users/akuntzle/data/akuntzle/simulations_donor_eur/{demo}/admix_csvs/{pop}_{eur_admix_pop_num}tracts_rep{repnum}_stitched.csv.gz')
mxb_df = pd.read_csv(f'/users/akuntzle/data/akuntzle/simulations_donor_eur/{demo}/admix_csvs/{pop}_9tracts_rep{repnum}_stitched.csv.gz')

arc_df["hap_id"] = arc_df["hap_id"].astype(int)
eur_df["hap_id"] = eur_df["hap_id"].astype(int)
mxb_df["hap_id"] = mxb_df["hap_id"].astype(int)

arc_df["Chromosome"] = "chr1"
eur_df["Chromosome"] = "chr1"
mxb_df["Chromosome"] = "chr1"

arc_df = arc_df.rename(columns={"start": "Start", "end": "End"})
eur_df = eur_df.rename(columns={"start": "Start", "end": "End"})
mxb_df = mxb_df.rename(columns={"start": "Start", "end": "End"})

common_eur_haps = set(arc_df.hap_id).intersection(eur_df.hap_id)
common_mxb_haps = set(arc_df.hap_id).intersection(mxb_df.hap_id)

eur_arc = []
mxb_arc = []

for hap in common_eur_haps:
    arc_sub = arc_df[arc_df.hap_id == hap]
    eur_sub = eur_df[eur_df.hap_id == hap]
    
    arc_gr = pr.PyRanges(arc_sub[["Chromosome","Start","End"]])
    eur_gr = pr.PyRanges(eur_sub[["Chromosome","Start","End"]])
    
    # Intersect
    overlap_eur = arc_gr.intersect(eur_gr)
    
    if len(overlap_eur) > 0:
        df_eur = overlap_eur.df.copy()
        df_eur["hap_id"] = hap
        eur_arc.append(df_eur)

for hap in common_mxb_haps:
    mxb_sub = mxb_df[mxb_df.hap_id == hap]
    arc_sub = arc_df[arc_df.hap_id == hap]

    mxb_gr = pr.PyRanges(mxb_sub[["Chromosome","Start","End"]])
    arc_gr = pr.PyRanges(arc_sub[["Chromosome","Start","End"]])

    overlap_mxb = arc_gr.intersect(mxb_gr)

    if len(overlap_mxb) > 0:
        df_mxb = overlap_mxb.df.copy()
        df_mxb["hap_id"] = hap
        mxb_arc.append(df_mxb)


final_eur_arc = pd.concat(eur_arc, ignore_index=True) if eur_arc else pd.DataFrame()
final_mxb_arc = pd.concat(mxb_arc, ignore_index=True) if mxb_arc else pd.DataFrame()

final_eur_arc.to_csv(
    f'/users/akuntzle/data/akuntzle/simulations_donor_eur/{demo}/admix_csvs/{pop}_eur_arc_ancestry_rep{repnum}.csv',
    index=False
)

final_mxb_arc.to_csv(
    f'/users/akuntzle/data/akuntzle/simulations_donor_eur/{demo}/admix_csvs/{pop}_amr_arc_ancestry_rep{repnum}.csv',
    index=False
)