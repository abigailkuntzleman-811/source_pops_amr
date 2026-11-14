# source_pops_amr

Analysis code for *[PAPER TITLE]*.

**Author:** Abigail Kuntzleman [, co-authors]
**Paper / preprint:** [DOI or link]

This repository is being organised one figure at a time. `fig1/` holds the
complete pipeline for **Figure 1**, the heatmap of shared archaic ancestry
between populations. The remaining analyses currently live in
`Chapter_1.ipynb` and will be added as they are tidied.

---

## Data

No large data files are stored here. `fig1/00_download_data.sh` retrieves
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
(`fig1/01_filter_hmmix_calls.py`). The 0.80 threshold follows the
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

All paths live in `config.sh` (shell) and `config.py` (Python). Edit those two
files, or export the same variables, and nothing else needs changing.

---

## Figure 1 — shared archaic ancestry heatmap

```bash
bash   fig1/00_download_data.sh                    # ~1.5 GB from Zenodo
python fig1/01_filter_hmmix_calls.py               # mean_prob >= 0.80
bash   fig1/02_segments_to_bed.sh                  # per-haplotype BED files
bash   fig1/03_pairwise_sharing.sh neanderthal     # pairwise intersections
python fig1/04_heatmap.py neanderthal              # -> results/figures/
```

Step 3 is the expensive one: every haplotype pair within and between 19
populations. On a Slurm cluster run it as an array job instead:

```bash
sbatch --array=1-190 fig1/03_pairwise_sharing.sh neanderthal
```

Repeat steps 3–4 with `denisovan` or `all` for the other panels.

### Outputs

| Paper item | Command |
|---|---|
| Figure 1, Neanderthal panel | `python fig1/04_heatmap.py neanderthal` |
| Figure 1, Denisovan panel | `python fig1/04_heatmap.py denisovan` |

Each run writes the figure and a CSV of the underlying matrix to
`results/figures/`, so the values behind the plot can be checked without
rerunning the pipeline.

---

## Table 1 — admixture proportion and shared archaic ancestry

Reuses the pairwise sharing output from Figure 1, so run `fig1/00`–`fig1/03`
first.

```bash
python table1/00_ancestry_proportions.py          # needs LAI_DIR (see below)
python table1/01_table1_values.py                 # -> results/tables/
```

`--unit individual` (the default) averages over people; `--unit haplotype`
averages over haplotypes. The means are nearly identical, the standard
deviations are not, so state in the methods which one the table reports.

### Local ancestry inference — not public

`table1/00_ancestry_proportions.py` needs local ancestry inference BED files
that are not publicly available. Point `LAI_DIR` at them:

```bash
export LAI_DIR=/path/to/local_ancestry      # expects 1kg_bed_files/<POP>_lai.bed
```

Without them, that step stops with an explanatory message and
`01_table1_values.py` still produces the shared-ancestry column, which depends
only on the public Skov (2024) calls.

[Add the data availability statement for these files here: who generated them,
and how a reader can obtain them.]

---

## Repository layout

```
config.py, config.sh         all paths and parameters
fig1/                        Figure 1 pipeline, numbered in run order
  00_download_data.sh        fetch public inputs
  01_filter_hmmix_calls.py   posterior probability filter
  02_segments_to_bed.py/.sh  per-haplotype BED files
  03_pairwise_sharing.sh     pairwise intersections (uses the two helpers)
  pairwise_hmmix.sh          between-population comparisons
  pairwise_within_hmmix.sh   within-population comparisons
  04_heatmap.py              builds the matrix and draws the figure
table1/                      Table 1, reuses fig1 steps 00-03
  00_ancestry_proportions.py continental ancestry proportions (needs LAI_DIR)
  01_table1_values.py        assembles the table
data/raw/                    downloaded inputs (not tracked)
data/processed/              intermediates (not tracked)
data/demographic_models/     demes models for the simulation analyses
results/figures/             generated figures
scripts/                     earlier versions, superseded by fig1/
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
