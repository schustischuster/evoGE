## Analyse rates of gene expression evolution of functional groups

This code allows to analyse the rate of gene expression evolution of functionally related genes. This involves estimating the stability of correlations using Monte-Carlo simulations, matching optimal number of control gene sets to each GO slim term category by assessing balance statistics, evaluating the strength of expression conservation in relation to gene expression levels, and applying non-linear negative exponential growth models to estimate the rates of gene expression evolution. 


## Contents

* [Getting Started](#getting-started)
  * [Required Packages](#required-packages)
  * [Data input](#data-input)
* [Data vizualization](#data-vizualization)
* [Session info](#session-info)

---
## Getting started


### Required Packages
Install and load the following R packages before running the reproducible scripts:

```R
# Create list of required packages
lib_List <- c("dplyr", "MatchIt", "gplots", "ggplot2", "scales")

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
Download and extract the entire directory to the working directory on your computer. Then, set the path for input and output files and source the R scripts: 

```R
in_dir <- file.path("evoGE", "GOslim", "data")
out_dir <- file.path("evoGE", "GOslim")

# Source R files
sourceDir <- function(path, trace = TRUE, ...) {
   for (nm in list.files(path, pattern = "[.][RrSsQq]$")) {
      if(trace) cat(nm,":")
      source(file.path(path, nm), ...)
      if(trace) cat("\n")
   }
}
 
sourceDir(path_to_R_files)
```
---
## Data analysis and vizualization

The following function will : 

```R
plotPhyloCore(div_times = c("Median", "Estimated"))

```
</br>

| Arguments  |  |
| :---  | :---  |
| div_times  | Indicates which pairwise divergence time to use. Must be either `"Median"` (median time derived from all studies) or `"Estimated"` (TTOL estimation). For more details, see [Hedges et al., MBE (2015)](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4379413/). |

</br>

To reproduce the results of this study, execute the following function call:

```R
plotPhyloCore(div_times = "Estimated")

```


---
## Session info

```R
sessionInfo()
```
