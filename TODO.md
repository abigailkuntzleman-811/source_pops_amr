# TODO


Grouped by kind of work. Figures come first because the text depends on them.
Figure numbers follow the current manuscript: S1 demographic models, S2 HMMix
recovered, S3 UpSet, S4 simulated recovered, S5 real IBD boxplots, S6 simulated
IBD boxplots, S7 Indigenous ancestry vs IBD, S8 African ancestry vs IBD, S9
ancestry-aware IBD ratios, S10 simulated donor vs standard. Tables: S1 shared vs
unique by background, S2 masked regions, S3 unrelated individuals.

## 1. Regenerate figures

- [x] Fig 2: recovered archaic ancestry from EUR/AMR tracts with the new joint-caller LAI calls (already in Overleaf)
- [x] Mask IBD pileups in the hap-ibd output (`hapIBD_calls/03_mask_pileups.py`; 12 regions, 30.8 Mb)
- [x] Fig 4, S7, S8: rerun with masked IBD, KING unrelated set and new LAI calls, and upload to Overleaf
- [x] S6: simulated IBD per pair with the corrected (joined) IBD (`figure_code/fig6/03_plot_sim_ibd_per_pair.py`), and upload to Overleaf
- [x] Fig 3: rerun with masked IBD and KING unrelated set, and upload to Overleaf
- [x] S9: run `figure_code/figS7_ancestry_aware_ibd/01_ancestry_aware_ibd.py` (masked IBD, new LAI, KING set), new table, and upload to Overleaf
- [x] S5: masked IBD and KING unrelated set (`figure_code/figS4_ibd_boxplots/01_plot_real_ibd_per_pair.py`), new table, and upload to Overleaf
- [x] Fig 6 and S10: rerun with the repo scripts (IBD fix, tract-overlap fix, 1000 replicates) and upload to Overleaf
- [x] Upload the remaining regenerated figures to Overleaf (Fig 6, S10)
- [x] Fig 6: bigger, legend inside upper left, y axis from 1, new y label; uploaded
- [x] S4: recomputed over 1000 replicates; uploaded
- [x] Fig 1 drawn from code as published (A Neanderthal, B Denisovan; numbers identical to the notebook)
- [x] Style pass: larger axis labels and text on Fig 1, 2, 3, 6, S4, S5, S6, S9, S10; "per pair" lower case; scientific-notation axes on Fig 1, 3, S4

## 2. Statistics

- [x] S5 and S6 statistics redone and in the captions (S5: bootstrap over individuals; S6: one value per replicate, Wilcoxon signed-rank vs IBS from the same replicate; both Holm-corrected; `figS4_ibd_boxplots/02_stats.py`, `fig6/04_sim_ibd_stats.py`)
- [x] S10: rerun with the new band (mean of the 1,000 per-replicate lines +- 95% CI across replicates; previously the CI of one regression on all samples pooled) and re-upload to Overleaf; caption: band = 95% CI across 1,000 replicates; test compares per-replicate slopes (Welch); call the differences differences in slope (cM per 50 Mb per unit EUR proportion), not "inflation of IBD"; PEL donor slope is flat (0.09)
- [x] S10 caption/Results: report the replicate-paired p for the slope difference (not the pooled interaction p), and the slopes themselves (e.g. PUR donor 7.93 vs standard 2.99 cM; with 1000 replicates any real difference gives a tiny p)

## 3. New analyses and additions

- [x] Reviewer request: shared vs unique archaic ancestry on European vs Indigenous tracts: supplementary Table S1 added with caption (archaic Mb now from the same 100 draws as the percentages; percentages unchanged). All archaic only (joint-caller calls have no Neanderthal/Denisovan label)
- [ ] Discuss Table S1 in the Results text and write the response to the reviewer (Methods done) (shared ~18-21% vs ~49-53% expected under the shuffled null, all p < 0.001; European-only largest in every population; noisy local ancestry pushes toward the null)
- [x] Add the UpSet plot to the paper (combined A Neanderthal, B Denisovan, C all; `upset_combined.png`)
- [ ] Discuss the UpSet figure (S3) in the Results text (Methods done)
- [x] Supplementary table of masked pileup regions (12 regions, 30.8 Mb; `hapIBD_calls/04_masked_regions_table.py` writes CSV and LaTeX)

## 4. Results and interpretation (after figures are final)

- [ ] Fig 2: interpretation; in the new figure PEL's EUR and AMR curves end nearly equal and PEL's AMR curve is above its EUR curve for most of the range
- [ ] Fig 4 and abstract: after masking, MXL is weakly positive (rho 0.26) and PEL shows no relationship (rho -0.07); the abstract and Results (caption done) say no correlation for MXL and a weak negative trend for PEL
- [ ] Fig 4 storyline: PUR and MXL now have the same rho (0.26, p = 0.007 and 0.038); only CLM is strongly positive (0.77). The "CLM and PUR positive vs MXL and PEL not" contrast (abstract, Results, Discussion) needs rethinking: the split is now CLM vs the rest
- [ ] Fig 6 / S10 framing: the multiple-donor simulations are framed as explaining MXL's lack of a positive correlation; revisit now that MXL is weakly positive. New Fig 6 mean slopes over 1000 replicates (+- SE): standard 0.35 +- 0.07, one donor 5.73 +- 0.09, two donors 2.67 +- 0.08, six donors 0.71 +- 0.07. Note: simulated IBD is a total over the 99 other haploid samples, not IBD_bar_i, so these slopes are not on the same scale as the real-data slopes in Fig 4 (do not compare them directly). Single replicates are noisy (slope SD 2.1-2.8, pooled R2 < 0.08), worth a clause in the caption
- [x] S7 caption: MXL no longer positive (rho -0.23, n.s.)
- [ ] S8 Results/Discussion (caption done): "negatively correlated in all populations" no longer holds for MXL (rho 0.11, n.s.)
- [ ] Fig 3 Results: "MXL has the least amount of IBD sharing of all four admixed groups" no longer holds; PEL is lowest (mean cM per haplotype pair: PEL 0.28, MXL 0.33, CLM 3.67, PUR 5.14)
- [ ] Fig 3 Discussion: "PEL shows levels most similar to those observed in East Asian populations" doesn't really hold; PEL (0.28) is as close to European (GBR 0.29, CEU 0.24) as to East Asian populations (CDX 0.35, CHS 0.22)
- [ ] Fig 3 Results/Discussion: PUR and CLM are now above FIN (2.80), not "comparable"/"more similar" to it
- [ ] Fig 3: "IBD and archaic sharing are not correlated": state the test (Spearman rho -0.37, p = 0.12 across 19 populations; Pearson p = 0.047, driven by PUR and CLM)
- [ ] S9 figure: re-upload (now drawn with the notebook's seaborn code: density per log10 unit, unfilled; the old version used density per linear unit, which pushed the peaks one to two decades left). Caption: histograms show only pairs with non-zero ratios; y axis is density per log10 unit, each ancestry normalized separately
- [ ] Re-upload Fig 4, S7 and S8 (confidence band restored to the notebook formula, t x residual SD; no star in titles)
- [ ] S9 caption: nearly done (example now 53.4%, letters dropped, non-zero pairs, axes). Left: decide whether to switch the y axis to proportion of pairs per bin (`stat="probability"`, same shapes, easier caption) or keep density per log10 unit; fix `log$_10$` -> `log$_{10}$`; move the no-IBD-anywhere sentence next to the table sentence
- [ ] S9 Results: the EUR < AMR < AFR gradient in pairs with no IBD is gone (now within 1-3 points in every population); "shared African ancestry contains fewer IBD segments more often" still holds but only slightly. Higher EUR ratio in CLM/PUR than MXL/PEL still holds (mean 0.0017-0.0021 vs 0.00016-0.00017)
- [ ] S9: pairs with zero IBD on shared ancestry are almost all pairs with no IBD anywhere (no IBD anywhere / EUR / AMR / AFR: CLM 53.2 / 53.4 / 53.4 / 55.4%, MXL 92.1 / 92.6 / 92.2 / 93.5%, PEL 92.1 / 93.5 / 92.1 / 95.1%, PUR 32.8 / 33.2 / 33.2 / 35.0%), so the zero fractions reflect overall IBD sharing, not ancestry; the non-zero ratio histogram is the ancestry comparison. (S9 checked: 30 pairs each in CLM/MXL/PEL recomputed by hand match exactly; same pairs as Fig 3.)
- [ ] S9 Discussion: "CLM and PUR distribution is bimodal ... wider distribution for MXL" is not visible in the new figure
- [ ] S9 / limitations: the joint-caller local ancestry switches every few 300-SNP windows (~1,350 switches per haplotype on chr10 vs ~10-30 expected), so ancestry at a given position on a given haplotype is noisy; affects S9 most (and Fig 2 / reviewer analysis), not ancestry proportions
- [ ] Results on S6 (standard model does not reproduce real IBD differences) and on simulated PEL in S10, once reruns are in

## 5. Methods and captions

- [x] IBD methods: describe pileup masking (100 kb bins above 3x median coverage; segment removal and trimming as in Guo et al. 2024) and report the hap-ibd version: 1.0 (15Jun23.92f), default parameters (min-seed 2.0, min-extend 1.0, min-output 2.0 cM, max-gap 1000 bp, min-markers 100, min-mac 2), HapMap GRCh38 PLINK maps linked from the hap-ibd README; results stable from 3x to 30x masking thresholds, and the chr10 pileup appears in all 26 populations
- [x] Methods: KING was run with `--unrelated --degree 3` (KING 2.3.2, LD-pruned SNPs from the hap-ibd input), so the unrelated set has no relatives up to **third** degree; the text says "no first or second degree relatives"
- [x] Simulation methods: say Fig 6 / S10 use total IBD each haploid sample shares with the other 99 over 50 Mb (not averaged); real-data IBD_bar_i is now defined per individual in the IBD Methods
- [x] Simulation methods: IBD from tskit with touching pieces joined into tracts of at least 2 cM, shared ancestors within 500 generations; admixture proportions from merged (non-overlapping) tracts
- [x] Clarify where the whole set vs the KING unrelated set is used, with sample sizes: Fig 2 and Table 1 use everyone with LAI calls (CLM 132, MXL 97, PEL 89, PUR 139); IBD analyses use KING unrelated individuals (Fig 4/S7/S8, with LAI: CLM 94, MXL 62, PEL 52, PUR 104)
- [x] Table 1: European admixture comes from the joint-caller LAI calls while shared Neanderthal ancestry comes from hmmix; resolve or state clearly
- [x] Methods: resolve the note to Jazeps on PEL counts (89 admixed PEL = 122 total - 33 used as the Indigenous American reference panel)
- [x] Replicate counts in captions: S10 done; S4 says 100 (now 1000); check Fig 6; S6 says 1000 (correct)
- [x] S5 and S6 captions: boxes, whiskers, haplotype pairs, zero pairs, tables
- [x] S10 caption: add that the y axis is each haploid genome's total IBD with the other 99, over a 50 Mb region (label now just cM); Fig 6 caption done. Re-upload Fig 6 (new across-replicate band, job running) and S10 after --replot for the new label
- [x] Fix cross-references: `Figure 3\ref{archaic_share_ibd}` prints "Figure 3 3" (Methods and Results); Discussion `Figure 3\label{archaic_share_ibd}` (\label should be \ref); Discussion sentence on EUR admixture vs IBD cites `\ref{eur_amr_recovered}` (Fig 2) instead of Fig 4

## 6. Checks

- [x] Table 1: new African proportions (CLM 0.15, PUR 0.20): checked; reference panels are all IBS and all YRI, values accepted
- [x] Masking used the maps in `data/raw/genetic_maps`: checked, same as hap-ibd's maps (chr1/chr10 median difference < 0.001 cM); regions before a map's first position (chr14, 15, 22) are extrapolated by hap-ibd

## 7. Housekeeping

- [ ] Delete the nine branches already merged into main (could not be done from Claude's session):
      `git push origin --delete explicit-exclusions fig2-pipeline fig4 fig4-pipeline figures-folder ibd-pipeline local-config reorganize-for-publication table1-pipeline`
- [ ] Decide on the `fig6` branch (its tract-overlap fix is now in main; its IBD code uses the old method)
- [ ] Remove the stray folder named `$HOME` in the repo on Oscar: `rm -r '$HOME'` (single quotes)
- [ ] Move the April lowercase `*_only_eur_/amr_archaic_recovered.txt` files out of `recovered_archaic_by_admix_prop`
