# source_pops_amr

Analysis code for *[PAPER TITLE]*.

**Author:** Abigail Kuntzleman [, co-authors]
**Paper / preprint:** [DOI or link]

This repository is being organised one figure at a time. `figure_code/fig1/` holds the
complete pipeline for **Figure 1**, the heatmap of shared archaic ancestry
between populations. The remaining analyses currently live in
`Chapter_1.ipynb` and will be added as they are tidied.

---

## Data

No large data files are stored here. `figure_code/fig1/00_download_data.sh` retrieves
everything that is public.

### Archaic introgression calls

Published hmmix calls for the 1000 Genomes cohort, in GRCh38 coordinates.

> Skov, L. (2024). *Archaic introgression for HGDP and 1000genomes in hg38*
> (v1.0.0) [Data set]. Zenodo. https://doi.org/10.5281/zenodo.14136628

Licensed CC BY 4.0. Segments were called with
[hmmix](https://github.com/LauritsSkov/Introgression-detection).

These calls are based on the **1000 Genomes 30× high-coverage callset on
GRCh38** (3,202 samples):

> Byrska-Bishop, M., et al. (2022). High-coverage whole-genome sequencing of
> the expanded 1000 Genomes Project cohort including 602 trios. *Cell*,
> 185(18), 3426–3440.e19. https://doi.org/10.1016/j.cell.2022.08.004

The call set covers 2,309 target individuals. African populations are absent
from the targets because they are used to construct the hmmix outgroup.

**Filtering applied here:** segments with mean posterior probability ≥ 0.80,
split by `ND_type` into Neanderthal and Denisovan sets
(`figure_code/fig1/01_filter_hmmix_calls.py`). The 0.80 threshold follows the
recommendation in the Zenodo record.

---

## Setup

```bash
pip install -r requirements.txt
```

Also required: **bedops** ≥ 2.4, which provides `bedops --intersect`.
On an HPC system it may already be available as a module, e.g.
`module load bedops`.

If you are running in an existing virtual environment rather than installing
fresh, activate it before running any of the scripts:

```bash
source /path/to/your/venv/bin/activate
```

All paths live in `config.sh` (shell) and `config.py` (Python).

For settings specific to your machine — where the local ancestry data lives,
which Python environment to activate — copy the example file instead of
editing the tracked ones:

```bash
cp config.local.sh.example config.local.sh
```

`config.local.sh` is not tracked by git, so your paths stay local. Both
`config.sh` and `config.py` read it, and Slurm jobs pick it up automatically,
including activating the environment named by `VENV_ACTIVATE`.

### Running on a Slurm cluster

Submit from the repository root, so the scripts can find `config.sh`:

```bash
cd /path/to/source_pops_amr
mkdir -p logs            # Slurm will not create the log directory itself
sbatch --array=0-7 figure_code/fig2/01_recovered_by_ancestry.script
```

If you must submit from elsewhere, set `REPO_ROOT` to the repository path.

Each `.script` has a commented placeholder near the top for activating a
Python environment; fill that in if your cluster needs it, since a Slurm job
does not inherit your interactive shell. The same applies to any variables you
exported by hand, such as `LAI_DIR` — set them in `config.sh` rather than
relying on the environment, or pass `--export` to `sbatch`.

---

## Figure 1 — shared archaic ancestry heatmap

```bash
bash   figure_code/fig1/00_download_data.sh                    # ~1.5 GB from Zenodo
python figure_code/fig1/01_filter_hmmix_calls.py               # mean_prob >= 0.80
bash   figure_code/fig1/02_segments_to_bed.sh                  # per-haplotype BED files
bash   figure_code/fig1/03_pairwise_sharing.sh neanderthal     # pairwise intersections
python figure_code/fig1/04_heatmap.py neanderthal              # -> results/figures/
```

Step 3 is the expensive one: every haplotype pair within and between 19
populations. On a Slurm cluster run it as an array job instead:

```bash
sbatch --array=1-190 figure_code/fig1/03_pairwise_sharing.sh neanderthal
```

Repeat steps 3–4 with `denisovan` or `all` for the other panels.

### Outputs

| Paper item | Command |
|---|---|
| Figure 1, Neanderthal panel | `python figure_code/fig1/04_heatmap.py neanderthal` |
| Figure 1, Denisovan panel | `python figure_code/fig1/04_heatmap.py denisovan` |

Each run writes the figure and a CSV of the underlying matrix to
`results/figures/`, so the values behind the plot can be checked without
rerunning the pipeline.

---

## Figure 2 — archaic sequence recovered, by continental ancestry

How much unique archaic sequence is recovered as more individuals are sampled,
split by the continental ancestry of the haplotype carrying it. Requires the
local ancestry inference files (see `LAI_DIR` below).

```bash
# one job per population x ancestry; see the .script for a Slurm array
python figure_code/fig2/01_recovered_by_ancestry.py MXL EUR
python figure_code/fig2/01_recovered_by_ancestry.py MXL AMR
# ... CLM, PEL, PUR

python figure_code/fig2/02_plot_recovered_by_ancestry.py        # -> results/figures/
```

On a Slurm cluster, run all eight combinations as an array:

```bash
sbatch --array=0-7 figure_code/fig2/01_recovered_by_ancestry.script
```

The plotting step takes `--ancestries` if you want a different set than the
default `EUR AMR`.

| Paper item | Command |
|---|---|
| Figure 2 | `python figure_code/fig2/02_plot_recovered_by_ancestry.py` |

---

## Running scripts as Slurm jobs

Any Python script here can run on a compute node instead of the login node:

```bash
mkdir -p logs                     # once
sbatch -J figS5 run.script figure_code/fig6/03_plot_sim_ibd_per_pair.py --reps 1000
tail -f logs/figS5_<jobid>.out    # follow progress
```

Everything after the script name is passed to it. Defaults: 12 h, 32 GB;
override on the command line (`sbatch -t 24:00:00 --mem=64G run.script ...`).
Submit with your Python environment active.

## IBD calling

Detects IBD segments in the 1000 Genomes 30x phased panel with hap-ibd. This
is independent of the archaic analyses and feeds Figure 3.

```bash
bash hapIBD_calls/00_setup_resources.sh             # hap-ibd.jar + genetic maps

mkdir -p logs
sbatch --array=1-22 hapIBD_calls/01_call_ibd.script  # one job per autosome
```

`00_setup_resources.sh` downloads hap-ibd and the PLINK-format GRCh38 HapMap
genetic maps, which hap-ibd needs in cM units.

Each chromosome is filtered to common (MAF > 0.05) biallelic SNPs with
duplicate positions removed, then passed to hap-ibd.

> Zhou, Y., Browning, S.R., Browning, B.L. (2020). A fast and simple method
> for detecting identity by descent segments in large-scale data.
> *American Journal of Human Genetics*, 106(4), 426-437.
> https://doi.org/10.1016/j.ajhg.2020.02.010

`02_within_population_ibd.py` sums segment lengths across autosomes for every
within-population haplotype pair. Pairs with no detected IBD count as zero, so
the denominator is the number of possible pairs rather than the number
observed. Related individuals are excluded: set `RELATED_EXCLUDE` to a file of
sample IDs (KING output works), or leave it unset to exclude the samples with
a parent recorded in the 1000 Genomes pedigree.

---

## Figure 3 — shared archaic ancestry versus shared IBD

```bash
bash   figure_code/fig1/03_pairwise_sharing.sh all        # archaic sharing, if not already run
python figure_code/fig3/01_within_population_ibd.py       # aggregate the hap-ibd output
python figure_code/fig3/02_plot_ibd_vs_archaic.py         # -> results/figures/
```

| Paper item | Command |
|---|---|
| Figure 3 | `python figure_code/fig3/02_plot_ibd_vs_archaic.py` |

Also writes `results/tables/archaic_share_vs_ibd.csv` with the plotted values.

---

## Figure 4 — admixture proportion against per-individual IBD

Each point is an individual: their proportion of a continental ancestry
against the average cM they share with each other member of their population.
Needs both the local ancestry data and the hap-ibd output.

```bash
python table1/00_ancestry_proportions.py            # needs LAI_BED_DIR
python figure_code/fig4/01_individual_ibd.py        # needs hapIBD_calls output
python figure_code/fig4/02_plot_admix_vs_ibd.py     # -> results/figures/
```

`--ancestry AMR` or `--ancestry AFR` plots the other ancestries.

| Paper item | Command |
|---|---|
| Figure 4 | `python figure_code/fig4/02_plot_admix_vs_ibd.py` |

The IBD average here is per individual, over the unrelated members of the
population: total cM shared with everyone else, divided by (n - 1). That is a
different quantity from Figure 3, which averages over haplotype pairs.

### IBD pileups

hap-ibd reports huge numbers of short segments in a few regions next to
centromeres (segmental duplications and genetic map errors), e.g. 99% of
chr10's 3.6 million segments overlap 49 Mb. These are masked before any IBD
analysis:

```bash
python hapIBD_calls/02_find_pileups.py     # optional: report pileups and their effect, changes nothing
python hapIBD_calls/03_mask_pileups.py     # write kgp_chr<N>_segments.masked.ibd.gz
```

`03_mask_pileups.py` flags 100 kb bins whose segment coverage exceeds 3x the
genome-wide median, removes segments entirely inside a flagged region, and
trims segments that cross one, keeping trimmed pieces of at least 2 cM (cM
from the genetic map hap-ibd used). The segment handling follows Guo et al.
(2024, Nat Commun 15:2499); the peak threshold does not, for the reason given
in the script's docstring. Figures 3 and 4 read the masked files by default;
set `IBD_MASKED=0` to use the raw hap-ibd output.

### Related individuals

All IBD analyses (Figures 3 and 4, and S6 and S8) use the maximal set of
unrelated individuals from KING (`king --unrelated`), which contains no first-
or second-degree relatives. KING writes the samples to drop to
`unrelated_toberemoved.txt` (two columns, FID and IID); point `RELATED_EXCLUDE`
at that file:

```bash
export RELATED_EXCLUDE="/path/to/unrelated_toberemoved.txt"
```

Each script prints `excluding N samples listed in <file>` when it is in use.
Without `RELATED_EXCLUDE`, the scripts fall back to dropping every sample with
a parent recorded in the 1000 Genomes pedigree, which removes fewer people and
gives slightly different sample sizes.

Analyses that do not use IBD (Figure 2, Table 1) keep every individual,
related ones included.

Spearman's rho and the regression p values are reported in the panel titles
and saved to `results/tables/`.

---

## Figure S7 — ancestry-aware IBD ratio

For each pair of haplotypes from different unrelated individuals, the fraction
of their shared EUR, AMR or AFR ancestry (bp) that is covered by IBD between
them. Uses the pileup-masked IBD, the per-population LAI files and the KING
unrelated set.

```bash
python figure_code/figS7_ancestry_aware_ibd/01_ancestry_aware_ibd.py
```

Writes per-pair tables and `ancestry_aware_ibd_summary.csv` to `TABLE_DIR`
and `ancestry_aware_ibd_ratios.png` to `FIGURE_DIR`. Assumes haplotypes 1 and
2 mean the same thing in the hap-ibd output and the LAI files.

---

## Figure 6 — simulated admixture models: admixture proportion against IBD

Needs the msprime simulations of each model; set `SIM_DIR` in
`config.local.sh` (one folder per model, each with `tree_seqs/` and
`admix_csvs/`).

```bash
python figure_code/fig6/01_plot_donor_ibd.py                # -> results/figures/, one PNG per step
python figure_code/fig6/compare_original_vs_fixed.py        # optional: old vs fixed IBD side by side
python figure_code/fig6/02_plot_donor_vs_standard_all_pops.py   # Figure S9: donor vs standard, 4 populations
python figure_code/fig6/03_plot_sim_ibd_per_pair.py           # Figure S5: IBD per pair, standard model
```

IBD in the simulations is computed by `figure_code/fig6/sim_ibd.py`. tskit
reports shared ancestry in pieces that break wherever the local tree changes,
so touching pieces are joined into tracts before a 2 cM cutoff (matching
hap-ibd), counting only shared ancestors from the last 500 generations.
Applying the cutoff to the raw pieces (tskit's `min_span`) undercounts IBD.

---

## Reviewer analysis — archaic sequence shared vs unique between ancestry backgrounds

Within each admixed population: of the archaic sequence, what percentage is
found on both European and Indigenous American tracts (shared), and what
percentage on only one. At each position with at least 20 haplotypes of each
ancestry, 20 of each are picked at random (100 times, averaged); the null
shuffles the EUR/AMR labels before picking (1000 times). Details in the
script's docstring.

```bash
mkdir -p logs
sbatch --array=0-3 figure_code/archaic_by_background/01_shared_unique_by_background.script
# or one population directly:
python3 figure_code/archaic_by_background/01_shared_unique_by_background.py PEL
```

Output: `$TABLE_DIR/archaic_by_background_<POP>.csv`.

## Table 1 — admixture proportion and shared archaic ancestry

Reuses the pairwise sharing output from Figure 1, so run `figure_code/fig1/00`–`figure_code/fig1/03`
first.

```bash
python table1/00_ancestry_proportions.py          # needs LAI_BED_DIR (see below)
python table1/01_table1_values.py                 # -> results/tables/
```

`--unit individual` (the default) averages over people and reproduces the
published Table 1 exactly. `--unit haplotype` averages over haplotypes
instead, giving the same means and marginally different standard deviations.
The table in the paper reports per-individual values.

### Local ancestry inference — not public

`table1/00_ancestry_proportions.py` needs the joint-caller local ancestry BED
files, which are not publicly available. Point `LAI_BED_DIR` at them:

```bash
export LAI_BED_DIR=/path/to/autosome_lai_bedfiles   # expects <POP>_lai.bed
```

Table 1 uses every individual in those files, related individuals included.
Figure 4 drops related individuals itself.

Without them, that step stops with an explanatory message and
`01_table1_values.py` still produces the shared-ancestry column, which depends
only on the public Skov (2024) calls.

[Add the data availability statement for these files here: who generated them,
and how a reader can obtain them.]

---

## Repository layout

```
config.py, config.sh         all paths and parameters
  fig1/                      Figure 1 pipeline, numbered in run order
  00_download_data.sh        fetch public inputs
  01_filter_hmmix_calls.py   posterior probability filter
  02_segments_to_bed.py/.sh  per-haplotype BED files
  03_pairwise_sharing.sh     pairwise intersections (uses the two helpers)
  pairwise_hmmix.sh          between-population comparisons
  pairwise_within_hmmix.sh   within-population comparisons
  04_heatmap.py              builds the matrix and draws the figure
figure_code/fig2/                        Figure 2, needs ARC_BED_DIR and LAI_BED_DIR
  01_recovered_by_ancestry.py   recovery curves per population and ancestry
  02_plot_recovered_by_ancestry.py  draws the figure
figure_code/fig3/                        Figure 3, needs fig1 and hapIBD_calls output
  01_within_population_ibd.py  within-population IBD per haplotype pair
  02_plot_ibd_vs_archaic.py    archaic sharing against IBD sharing
fig4/                        Figure 4, needs LAI_BED_DIR and hapIBD_calls output
  01_individual_ibd.py         average IBD per individual
  02_plot_admix_vs_ibd.py      admixture proportion against IBD
hapIBD_calls/                IBD calling with hap-ibd
  00_setup_resources.sh      downloads hap-ibd.jar
  01_call_ibd.sh/.script     filter the panel and call IBD, per chromosome
table1/                      Table 1, reuses fig1 steps 00-03
  00_ancestry_proportions.py continental ancestry proportions (needs LAI_BED_DIR)
  01_table1_values.py        assembles the table
data/raw/                    downloaded inputs (not tracked)
data/processed/              intermediates (not tracked)
data/demographic_models/     demes models for the simulation analyses
results/figures/             generated figures
scripts/                     earlier versions, superseded by figure_code/fig1/
Chapter_1.ipynb              full exploratory notebook
```

---

## Citation

```
[BibTeX once the paper has a DOI]
```

## License

Code: MIT (see `LICENSE`). Input data keep the licenses of their original
sources; the Skov (2024) calls are CC BY 4.0.

## Contact

Abigail Kuntzleman — questions welcome via the Issues tab.
