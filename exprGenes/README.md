
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

### Retrieve number of expressed genes

The following function will apply a threshold function based on ERCC spike-ins at different threshold levels. 

```R
getExprGenes(species = c("ATH", "AL", "CR", "ES", "TH", "MT", "BD"), 
             experiment = c("single-species", "comparative"), threshold)

```
To generate the panels of the figures, execute the following function calls:

```R
thresholds <- list(0, 0.01, 0.05, 0.1) # ERCC threshold values are 0 (a fixed TPM threshold of 0.05)
# or perc of expressed spike-ins for 0.01/0.05/0.1

lapply(thresholds, getExprGenes, species = "ATH", experiment = "single-species")
lapply(thresholds, getExprGenes, species = "AL", experiment = "comparative")
lapply(thresholds, getExprGenes, species = "CR", experiment = "comparative")
lapply(thresholds, getExprGenes, species = "ES", experiment = "comparative")
lapply(thresholds, getExprGenes, species = "TH", experiment = "comparative")
lapply(thresholds, getExprGenes, species = "MT", experiment = "comparative")
lapply(thresholds, getExprGenes, species = "BD", experiment = "comparative")

```

## Visualization

Set the file path for the data generated in the previous steps and source the R script:

```R
in_dir_stats <- file.path("exprGenes", "output", "mapping_statistics")
in_dir_expr_genes <- file.path("exprGenes", "output", "expr_genes")

source(file.path("exprGenes", "R", "exprGenes_plots.R"))

```

The plotting functions will generate the panels for the following figures:
