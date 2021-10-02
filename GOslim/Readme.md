## Analyse rates of gene expression evolution of functional groups

This code allows to analyse the rate of gene expression evolution of functionally related genes. This involves estimating the stability of correlations using Monte-Carlo simulations, matching optimal number of control gene sets to each GO slim term category by assessing balance statistics, evaluating the strength of expression conservation in relation to gene expression levels, and applying non-linear regression models to estimate the rates of gene expression evolution. 


## Contents

* [Getting Started](#getting-started)
  * [Required Packages](#required-packages)
  * [Data input](#data-input)
* [Data analysis and visualization](#data-analysis-and-visualization)
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
path_to_R_files <- file.path("evoGE", "GOslim", "R")

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
## Data analysis and visualization

Sample correlations converge to the population value with increasing sample size, and it has been shown that the sample size should approach n=250 for stable estimates ([Schönbrodt and Perugini, J Res Pers. 2013](https://www.sciencedirect.com/science/article/abs/pii/S0092656613000858)). Since the points of stability (POS) published in this work were based on normal distributions, we implemented Monte-Carlo simulations to retrieve sample size estimates for stable correlations derived from gene expression data.

```R
estimatePOS(nbootstrap, coswidth, bss, ...)

```
</br>

| Arguments  |  |
| :---  | :---  |
| nbootstrap  | Number of sampling trajectories. |
| coswidth  | Indicates the corridor of stability (COS). Values of 0.1/0.15/0.2 correspond to small/medium/large effect sizes, respectively. For more details, see [Cohen, Psychol Bull (1992)](https://pubmed.ncbi.nlm.nih.gov/19565683/). |
| bss  | Indicates the level of confidence for the confidence interval. |

</br>

To reproduce the results of this study, execute the following function call:

```R
estimatePOS(nbootstrap = 1000, coswidth = 0.15, bss = 0.95)

```

Next, we wanted to test wether the evolutionary stability of gene subsets is affected by gene expression levels. We therefore generated subsets of orthologous genes according to quantiles of average expression across either all samples (inter-organ inter-species), or within the same organ across species (intra-organ inter-species). For each quantile set, we then calculated metric pearson distances and fitted non-linear regression models to estimate quantile-specific rates of gene expression evolution.

```R
getExprCons(nquant, qtype = c("base_mean", "organ_spec"), ...)

```
</br>

| Arguments  |  |
| :---  | :---  |
| nquant  | Number of quantiles. |
| qtype  | Type of quantile. Use "base_mean" for quantiles of average gene expression across organs and species, and "organ_spec" for organ-specific (intra-organ) quantiles of average gene expression across species. |

</br>

To reproduce the results of this study, execute the following function call:

```R
getExprCons(nquant = 500, qtype = "base_mean")
getExprCons(nquant = 500, qtype = "organ_spec")

```

---
## Session info

```R
sessionInfo()
```








