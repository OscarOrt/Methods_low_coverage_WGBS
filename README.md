# Low coverage Whole Genome Bisulfite Sequencing

This repository is intended to provide the code used for the bioinformatic analysis for the methods paper Estimating global methylation and erasure using low-coverage whole genome bisulfite sequencing (WGBS).

The code is implemented in R (v3.4.4) and the following libraries:
- ggplot2 3.0.0
- dplyr

This reposity contains the scripts, the data used for the analysis and the genome used for upstream analysis.

## Components of the analysis

- Plot margin of error theoretical model:
This code uses a theretical model to predict the margin of error in low-coverage WGBS.

- Plot example data bovine blastocyst:
This code uses the theoretical model to estimate the margin of error in control (fibroblasts) and bovine blastocyst samples.

- Plot dispersion bootstraping model:
This code plots the relation between number of CGs sampled and the porcetage of methylation obtained for each sample.

- Plot margin of error bootstraping model:
This code uses the bootstraping sampling data to estimate the margin of error at different CGs sample numbers for two samples.

- Plot margin of error bootstraping model and first CG:
This code adds the results obtained when just the first CG from each reads is analysed to estimate the ME.