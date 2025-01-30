[![DOI](https://zenodo.org/badge/826514679.svg)](https://doi.org/10.5281/zenodo.14774593)

This repository contains code used to analyze data reported in the Alberge, Dutta, Poletti, et al. paper. 

## Data availability

Most data from the PCROWD study reported as average are available to other researchers upon request to Dr Ghobrial (mailto: Irene_Ghobrial [at] DFCI[dot]HARVARD[dot]EDU). Sequencing data will be made available shortly on the dbGAP study phsxxxx [New Stury Registered].

## Code used for Figures

* Figure 1A

The code for figure and Fisher's Exact tests are in: `code/2.1bis-create-new-matrix.R`. Metadata are provided in Supplementary Table S1 (List of samples), S2 (List of MutSig2CV drivers), S3 (List of GISTIC2.0 drivers), S4 (Odds ratios and significance of driver enrichments), S9 (List of SV driver candidates from SSVGAR)

* Figure 1BCD

The code for figures and statistical tests is in `code/2.2_compute-MM-like-score.R`.

* Figure 2

Data is provided in Supplementary Table S5, exact lists of somatic mutations are accessible upon request with dbGAP access rights.

* Figure 3ABC

Timing data is provided in Supplementary Table S6, while code to perform the analysis are in `code/5.1_clocklike-timing-EM-code.R` and to reproduce the figures and statistical tests in `code/5.2_clocklike-timing-figures-statistics.R`.

* Figure 3D

Collaboration Dr. Kübler and Dr. Arndt [data to be uploaded].

* Figure 3E

Collaboration Dr. Poletti, Dr. Zamagni, Dr. Terragna, code and figures for the League Models are in [TiMMing](https://github.com/andrea-poletti-unibo/MS_TiMMing) repository. 
Code for Figure 3E is in `code/5.3_clonality-league-comparison-bradley-terry-models-figure.R`. Values from the model are provided in Supplementary Table S7.

* Figure 4A

Collaboration Dr. Coorens, code is in `code/HDP_complete_workflow_TIMCOORENS.R` and `code/deconvolute_hdp_sigs_TIMCOORENS.R`. Frequencies of substitutions are provided in Supplementary Table S8.

* Figure 4BC

Code is in `code/3.2_signatures-corr-models-dup-ccf.R`.

* Figure 5AB

Collaboration SSVGAR Loinaz, code is in companion repository [SSVGAR)](https://github.com/getzlab/SSVGAR), list of drivers is provided in Supplementary Table S9

* Figure 5C

Code is in `code/2.5-SV-circos-MYC.R`

* Figure 6A

Code is in `code/6.1_DIG-noncoding-driver-analysis.R`, detailed list of drivers is provided in Supplementary Table S10.

* Figure 6BC

Code is in `code/6.2_noncoding-ILF2-promoter-zoom.R`.

* Figure 6DE

Code is in `code/6.3_noncoding-ILF2_commpass_expr.R`, detailed participant IDs from the CoMMpass study and their mutational status is provided in Supplementary Table S11 and S12.
