# TODO

Grouped by kind of work. Figures come first because the text depends on them.
Figure numbers follow the current manuscript: S4 real IBD boxplots, S5 simulated
IBD boxplots, S6 Indigenous ancestry vs IBD, S7 ancestry-aware IBD ratios,
S8 African ancestry vs IBD, S9 simulated donor vs standard.

## 1. Regenerate figures

- [x] Fig 2: recovered archaic ancestry from EUR/AMR tracts with the new joint-caller LAI calls (already in Overleaf)
- [x] Mask IBD pileups in the hap-ibd output (`hapIBD_calls/03_mask_pileups.py`; 12 regions, 30.8 Mb)
- [ ] Fig 4, S6, S8: rerun with masked IBD, KING unrelated set and new LAI calls, and upload to Overleaf
- [ ] S5: simulated IBD per pair with the corrected (joined) IBD (`figure_code/fig6/03_plot_sim_ibd_per_pair.py`), and upload to Overleaf
- [x] Fig 3: rerun with masked IBD and KING unrelated set
- [ ] S7: run `figure_code/figS7_ancestry_aware_ibd/01_ancestry_aware_ibd.py` (masked IBD, new LAI, KING set)
- [ ] S7: confirm with Jazeps that hap-ibd and the joint caller used the same phased haplotypes (haplotype 1/2 must mean the same thing in both)
- [ ] S4: in the notebook, read the masked IBD files (`kgp_chr{i}_segments.masked.ibd.gz`) and use the KING unrelated set (its current table counts full samples, related individuals included)
- [ ] Fig 6 and S9: check the reruns finished with the IBD fix and the tract-overlap fix, 1000 replicates
- [ ] Upload all regenerated figures to Overleaf (Fig 3, Fig 4, S4–S9)

## 2. Statistics

- [ ] Redo the tests on IBD distributions: real data (S4) with a bootstrap over individuals, simulations (S5) with one value per replicate; decide and state whether zero-IBD pairs are included
- [ ] S9: report the replicate-level test of the slope difference (`sims_IBD_admix_prop_stats.csv`) rather than the pooled interaction p-value
- [ ] Optional: robustness check with IBD segments of at least 3 or 4 cM

## 3. New analyses and additions

- [ ] Reviewer request: shared vs unique archaic ancestry on European vs Indigenous tracts within each population. Script written and tested: run `sbatch --array=0-3 figure_code/archaic_by_background/01_shared_unique_by_background.script` (20 haplotypes per background, 100 random picks, 1000 label shuffles)
- [ ] Add the UpSet plot to the paper
- [ ] Add the UpSet analysis to the methods and supplement (downsample each population to the smallest sample size, then compute all possible intersections)
- [ ] Supplementary table of masked pileup regions (`hapIBD_calls/04_masked_regions_table.py` writes CSV and LaTeX)

## 4. Results and interpretation (after figures are final)

- [ ] Fig 2: interpretation; in the new figure PEL's EUR and AMR curves end nearly equal and PEL's AMR curve is above its EUR curve for most of the range
- [ ] Fig 4 and abstract: after masking, MXL is weakly positive (rho 0.26) and PEL shows no relationship (rho -0.07); the abstract, Results and Fig 4 caption say no correlation for MXL and a weak negative trend for PEL
- [ ] Fig 6 / S9 framing: the multiple-donor simulations are framed as explaining MXL's lack of a positive correlation; revisit now that MXL is weakly positive
- [ ] S6 caption: MXL no longer positive (rho -0.23, n.s.)
- [ ] S8 caption and Results/Discussion: "negatively correlated in all populations" no longer holds for MXL (rho 0.11, n.s.)
- [ ] Fig 3 Results: "MXL has the least amount of IBD sharing of all four admixed groups" no longer holds; PEL is lowest (mean cM per haplotype pair: PEL 0.28, MXL 0.33, CLM 3.67, PUR 5.14)
- [ ] Fig 3 Discussion: "PEL shows levels most similar to those observed in East Asian populations" doesn't really hold; PEL (0.28) is as close to European (GBR 0.29, CEU 0.24) as to East Asian populations (CDX 0.35, CHS 0.22)
- [ ] Fig 3 Results/Discussion: PUR and CLM are now above FIN (2.80), not "comparable"/"more similar" to it
- [ ] Fig 3: "IBD and archaic sharing are not correlated": state the test (Spearman rho -0.37, p = 0.12 across 19 populations; Pearson p = 0.047, driven by PUR and CLM)
- [ ] Results on S5 (standard model does not reproduce real IBD differences) and on simulated PEL in S9, once reruns are in

## 5. Methods and captions

- [ ] IBD methods: describe pileup masking (100 kb bins above 3x median coverage; segment removal and trimming as in Guo et al. 2024) and report the hap-ibd version: 1.0 (15Jun23.92f), default parameters (min-seed 2.0, min-extend 1.0, min-output 2.0 cM, max-gap 1000 bp, min-markers 100, min-mac 2), HapMap GRCh38 PLINK maps linked from the hap-ibd README; results stable from 3x to 30x masking thresholds, and the chr10 pileup appears in all 26 populations
- [ ] IBD_bar_i definition: real data is per individual (N unrelated individuals, IBD_ij summed over both haplotypes of each); "N haploid samples" only fits the simulations
- [ ] Simulation methods: IBD from tskit with touching pieces joined into tracts of at least 2 cM, shared ancestors within 500 generations; admixture proportions from merged (non-overlapping) tracts
- [ ] Clarify where the whole set vs the KING unrelated set is used, with sample sizes: Fig 2 and Table 1 use everyone with LAI calls (CLM 132, MXL 97, PEL 89, PUR 139); IBD analyses use KING unrelated individuals (Fig 4/S6/S8, with LAI: CLM 94, MXL 62, PEL 52, PUR 104)
- [ ] Table 1: European admixture comes from the joint-caller LAI calls while shared Neanderthal ancestry comes from hmmix; resolve or state clearly
- [ ] Methods: resolve the note to Jazeps on PEL counts (89 admixed PEL = 122 total - 33 used as the Indigenous American reference panel)
- [ ] Replicate counts in captions: S9 says 100 (now 1000); check Fig 6; S5 says 1000 (correct)
- [ ] S4 caption: say what the boxes show (zero-IBD pairs excluded, whiskers 1.5x IQR) and the unit of a pair; S5: pairs are haplotype pairs
- [ ] Fix cross-references: `Figure 3\ref{archaic_share_ibd}` prints "Figure 3 3" (Methods and Results); Discussion `Figure 3\label{archaic_share_ibd}` (\label should be \ref); Discussion sentence on EUR admixture vs IBD cites `\ref{eur_amr_recovered}` (Fig 2) instead of Fig 4

## 6. Checks

- [ ] Literature: is the low within-population IBD in PEL and MXL after masking (about 1 cM per partner) consistent with published estimates (e.g. Mooney et al. 2018)?
- [ ] Table 1: new African proportions (CLM 0.15, PUR 0.20) look high next to published estimates; check the joint caller's reference panels
- [ ] Masking used the maps in `data/raw/genetic_maps` (not found in the hap-ibd logs); confirm they are the maps hap-ibd used

## 7. Housekeeping

- [ ] Delete the nine branches already merged into main (could not be done from Claude's session):
      `git push origin --delete explicit-exclusions fig2-pipeline fig4 fig4-pipeline figures-folder ibd-pipeline local-config reorganize-for-publication table1-pipeline`
- [ ] Decide on the `fig6` branch (its tract-overlap fix is now in main; its IBD code uses the old method)
- [ ] Remove the stray folder named `$HOME` in the repo on Oscar: `rm -r '$HOME'` (single quotes)
- [ ] Move the April lowercase `*_only_eur_/amr_archaic_recovered.txt` files out of `recovered_archaic_by_admix_prop`
