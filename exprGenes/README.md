
## Expressed genes and data statistics

This code allows to summarize the DevSeq data statistics and to reproduce the results of the intra-species transcriptome analyses. 


## Contents

* [Getting Started](#getting-started)
  * [Required Packages](#required-packages)
  * [Data input](#data-input)
* [Data analysis](#data-analysis)
  * [Retrieve mapping statistics](#retrieve-mapping-statistics)
* [Visualization](#visualization)
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
Download the entire subdirectory containing the [data](https://github.com/schustischuster/evoGEx/tree/master/exprGenes/data) folder and [R scripts](https://github.com/schustischuster/evoGEx/tree/master/exprGenes/R) to the working directory on your computer, e.g. by using [GitZip](http://kinolien.github.io/gitzip/), and extract the file. Then, set the path for input and output files and source the R scripts: 

```R
in_dir <- "./exprGenes/data"
out_dir <- "./exprGenes"

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
