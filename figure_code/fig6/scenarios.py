"""The simulation scenarios behind Figure 6.

Each scenario is one demographic model, simulated for SIM_REPS replicates. The
directory names are historical and not always self-describing, so the mapping
to models and to the figure's labels lives here rather than being spread
across the scripts.

`tracts` is the filename of the stitched ancestry tract file within a
scenario's admix_csvs/ directory. The number in those names is a *population
id* in the demographic model, not a count: 7 is IBS in the standard model, 10
is IBS_DONOR in the donor model. Two scenarios read files merged across
several donor populations, which 02_extract_tracts.py does not produce by
itself.
"""

SCENARIOS = {
    "no_bottleneck_donor": dict(
        label="Standard",
        # No donor population: European ancestry enters directly from IBS.
        # The directory name is misleading and kept only for continuity with
        # the existing output.
        demes="mxl_introgression.yml",
        admix_pop="IBS",
        tracts="{pop}_7tracts_rep{rep}_stitched.csv.gz",
        color="#0072B2",
        linestyle="-",
    ),
    "eur_donor": dict(
        label="Single 10% European donor",
        demes="mxl_introgression_donor.yml",
        admix_pop="IBS_DONOR",
        tracts="{pop}_10tracts_rep{rep}_stitched.csv.gz",
        color="#CC79A7",
        linestyle=":",
    ),
    "mxl_double_donor_recent": dict(
        label="Two 10% European donors",
        demes="mxl_introgression_double_donor_recent.yml",
        admix_pop=None,            # merged across two donor populations
        tracts="{pop}_merged_89tracts_rep{rep}.csv.gz",
        color="#009E73",
        linestyle="-.",
    ),
    "mxl_donor6": dict(
        label="Six 10% European donors",
        demes="mxl_introgression_donor6.yml",
        admix_pop=None,            # merged across six donor populations
        tracts="{pop}_totalIBSmerged_rep{rep}.csv.gz",
        color="#D55E00",
        linestyle=(0, (5, 1)),
    ),
}

# Order the layers are revealed in, and plotted in.
PLOT_ORDER = [
    "no_bottleneck_donor",
    "eur_donor",
    "mxl_double_donor_recent",
    "mxl_donor6",
]
