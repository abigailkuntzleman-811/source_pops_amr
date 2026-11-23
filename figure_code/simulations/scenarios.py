"""Simulation scenarios behind Figures 6, S5 and S9.

For each scenario: the folder under SIM_DIR, the demes model in
data/demographic_models/ (with {pop} for the admixed population), and the
European source population(s) whose tracts are extracted. The number in the
tract file names (e.g. MXL_7tracts_...) is that source's population id in the
model, not a count.
"""

SCENARIOS = {
    # standard model: European ancestry enters directly from IBS (id 7)
    "no_bottleneck_donor": dict(demes="{pop_lower}_introgression.yml", sources=["IBS"]),
    # one bottlenecked donor, 10% of IBS (id 10)
    "eur_donor": dict(demes="{pop_lower}_introgression_donor.yml", sources=["IBS_DONOR"]),
    # MXL only: two donors (ids 8, 9), merged into {pop}_merged_89tracts_rep{rep}.csv.gz
    "mxl_double_donor_recent": dict(demes="mxl_introgression_double_donor_recent.yml",
                                    sources=["IBS_split1", "IBS_split2"],
                                    merged="{pop}_merged_89tracts_rep{rep}.csv.gz"),
    # MXL only: six donors (ids 8-13), merged into {pop}_totalIBSmerged_rep{rep}.csv.gz
    "mxl_donor6": dict(demes="mxl_introgression_donor6.yml",
                       sources=[f"IBS_split{k}" for k in range(1, 7)],
                       merged="{pop}_totalIBSmerged_rep{rep}.csv.gz"),
}


def tract_file(pop, source_id, rep):
    return f"{pop}_{source_id}tracts_rep{rep}_stitched.csv.gz"
