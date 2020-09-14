## Make angiosperm phylogeny and ortholog gene plots

This code allows to generate the angiosperm phylogeny and coding/non-coding orthologous gene plots. 


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
lib_List <- c("dplyr", "ggplot2", "ape", "scales", "gtable")

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
in_dir <- file.path("orthologs", "data")
out_dir <- file.path("orthologs")

source(file.path("orthologs", "R", "plotPhyloCore.R"))

```
---
## Data analysis and vizualization

The following function will generate the angiosperm phylogeny and orthologous gene plots:: 

```R
plotPhyloCore(div_times = c("Median", "Estimated"))

```
</br>

| Arguments  |  |
| :---  | :---  |
| div_times  | Indicates which pairwise divergence time to use. Must be either `"Median"` (median time derived from all studies) or `"Estimated"` (TTOL estimation). For more details, see [Hedges et al., MBE (2015)](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4379413/). |

</br>

To reproduce the results of this study, execute the following function calls:

```R
plotPhyloCore(div_times = "Median")
plotPhyloCore(div_times = "Estimated")

```

This will generate the panels for the following figures:


---
## Session info

```R
sessionInfo()
```
