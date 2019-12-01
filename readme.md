# Mapping Tumor-Specific Expression QTLs in Impure Tumor Samples

# Author Contributions Checklist Form

## Data

### Abstract

We have used gene expression and genotype data obtained from The Cancer Genome Atlas
(TCGA). All the data are publicly available

### Availability

We applied and were approved for TCGA data access through NCBI dbGap.

All the data we used can be downloaded from NCI Genomic Data Commons (https://gdc.cancer.gov/). While new versions of data have been generated, the version that we used could be downloaded from their legacy archive: https://portal.gdc.cancer.gov/legacyarchive/search/f.

Genotype data were downloaded as CEL file format. Gene expression data were downloaded
as BAM files. We have provided details on data processing in the supplementary materials of
this paper, including the exact command that we used to process the data step by step.

## Code

### Abstract

We have implemented our code in an R package named pTReCASE, whch is available at
https://github.com/Sun-lab/pTReCASE. In addition, we have also provided our codes for data
processing, simulation and real data analysis as a supplementary zip file of this paper.

### Description

In the folder of “data_preparation”, we provided data preparation pipelines to prepare gene
expression and genotype data for real data analysis. Most files have associated .Rout files that show the intermediate output.

In the folder “ex_real_data_analysis”, we provided data as well as R code used in real data
analysis while replacing real genotype data with simulated ones.

In the folder “simulation”, we provided codes as well as results for simulation analysis.

More details are provided in a README.txt file


## Instructions for use

### Reproducibility

All the results of our analyses can be reproduced. We have provided the R code to reproduce
the simulation results. The real data analyses take multiple steps because we have started from raw data. We have included the details of data processing and we provide processed gene
expression data. However, genotype are not allowed to be cannot be shared in a website
without access control. Anyone who is interested in accessing genotype data would need to
download it from TCGA genomic data commons.

