# PRS-calculation_LDpred2_Lassosum!
This codes will help to calculate PRS score using the online available datasets

**Polygenic Risk Scores:**
* PRS is the only method that provides an estimate of genetic liability to a trait at the individual level.
* It is calculated using GWAS summary statistics (Base data) and a target data.
* Calculation - By computing the sum of risk alleles that an individual has, weighted by the risk allele effect sizes as extimated by a GWAS on the phenotype.
* Chatterjee said "PRS of an individual is defined as a quantitative measure of the disease over susceptibility variants"

**PRS R packages:**
1. LDpred2
2. Lassosum

**LDpred2**
LDpred-2 is one of the dedicated PRS programs which is an R package that uses a Bayesian approach to polygenic risk scoring.

You can install LDpred and its dependencies in R with the following command: 

install.packages("remotes")
library(remotes)
remotes::install_github("https://github.com/privefl/bigsnpr.git")

**Lassosum**
lassosum is one of the dedicated PRS programs which is an R package that uses penalised regression (LASSO) in its approach to PRS calculation.

You can install Lassosum and its dependencies in R with the following command:

install.packages(c("devtools","RcppArmadillo", "data.table", "Matrix"), dependencies=TRUE)
library(devtools)
install_github("tshmak/lassosum")
