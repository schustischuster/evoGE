
## Perform comparative analysis of DevSeq and Brawand data

This code allows to reproduce the inter-organ intra-species, inter-species and cross-kingdom analyses of the DevSeq and Brawand [(Brawand et al., 2011)](https://pubmed.ncbi.nlm.nih.gov/22012392/) data sets. 


## Contents

* [Getting Started](#getting-started)
  * [Required Packages](#required-packages)
  * [Data input](#data-input)
* [Data analysis and vizualization](#data-analysis-and-vizualization)
* [Session info](#session-info)


## Getting started


### Required Packages
Install and load the following R packages before running the reproducible scripts:

```R
if (!require(plyr)) install.packages('plyr')
library(plyr)
if (!require(dplyr)) install.packages('dplyr')
library(dplyr)
if (!require(ggplot2)) install.packages('ggplot2')
library(ggplot2)
if (!require(mgcv)) install.packages('mgcv')
library(mgcv)
if (!require(grid)) install.packages('grid')
library(grid)
if (!require(scales)) install.packages('scales')
library(scales)

```

### Data input
Download the entire subdirectory containing the [data](https://github.com/schustischuster/evoGEx/tree/master/compAnalysis/data) folder and [R script](https://github.com/schustischuster/evoGEx/tree/master/compAnalysis/R) to the working directory on your computer, e.g. by using [GitZip](http://kinolien.github.io/gitzip/), and extract the file. Then, set the path for input and output files and source the R scripts: 

```R
in_dir <- "./compAnalysis/data"
out_dir <- "./compAnalysis"

source("exprGenes/R/makeCompAnalysis.R")

```

## Data analysis and vizualization

The following function will load and analyze the DevSeq and Brawand ortholog expression data and generate the plots. 

```R
makeCompAnylsis(dataset = c("Brawand", "DevSeq"), expr_estimation = c("TPM", "counts"), coefficient = c("pearson", "spearman"))

```
To reproduce the plots of this study, run the function with the following parameters:

```R
makeCompAnylsis(dataset = "Brawand", expr_estimation = "TPM", coefficient = "spearman")
makeCompAnylsis(dataset = "DevSeq", expr_estimation = "TPM", coefficient = "spearman")
makeCompAnylsis(dataset = "DevSeq", expr_estimation = "TPM", coefficient = "pearson")

```



The plotting functions will generate the panels for the following figures:
