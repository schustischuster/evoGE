
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

# Create list of required packages
lib_List <- c("dplyr", "ggplot2", "data.table", "mgcv", "grid", "gtable", "scales", "factoextra", "dendextend")

loadLibrary <- function(x) { 
    if (!require(x, character.only = T)) {
        install.packages('x')
        library(x)
    }
}

# Load packages
invisible(lapply(lib_List, loadLibrary))

```

### Data input
Download the entire subdirectory containing the `data` and `R` folders to the working directory on your computer, e.g. by using [GitZip](http://kinolien.github.io/gitzip/), and extract the file. Then, set the path for input and output files and source the R scripts:  

```R
in_dir <- file.path("exprGenes", "data")
out_dir <- file.path("exprGenes")

source(file.path("exprGenes", "R", "getStats.R"))
source(file.path("exprGenes", "R", "getExprGenes.R"))

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
in_dir_cd <- file.path("exprGenes", "output", "mapping_statistics")

source(file.path("exprGenes", "R", "exprGenes_plots.R"))

```

The plotting functions will generate the panels for the following figures:
