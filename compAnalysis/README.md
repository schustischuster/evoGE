
## Perform comparative analysis of DevSeq and Brawand data

This code allows to reproduce the intra-species, inter-species and cross-kingdom analyses of the DevSeq and Brawand [(Brawand et al., 2011)](https://pubmed.ncbi.nlm.nih.gov/22012392/) data sets. 


## Contents

* [Getting Started](#getting-started)
  * [Required Packages](#required-packages)
  * [Data input](#data-input)
* [Data analysis and vizualization](#data-analysis-and-vizualization)
  * [Retrieve mapping statistics](#retrieve-mapping-statistics)
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

source("exprGenes/R/getStats.R")
source("exprGenes/R/getExprGenes.R")

```

## Data analysis

### Retrieve mapping statistics

The following function will merge the DevSeq mapping statistics and create a data table for each the _A.thaliana_ data set, the non-ATH data, and the comparative data set. 

```R
getStats()

```

## Visualization

Set the file path for the data generated in the previous steps and source the R script:

```R
in_dir_cd <- "./output/mapping_statistics"

source("exprGenes_plots.R")

```

The plotting functions will generate the panels for the following figures:
