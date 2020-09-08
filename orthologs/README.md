## Make angiosperm phylogeny and orthologous gene plots

This code allows to generate the angiosperm phylogeny and coding and non-coding orthologous gene plots. 


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
lib_List <- c("dplyr", "ggplot2", "ape")

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

The following function will load and analyze the DevSeq and Brawand ortholog expression data and generate the plots: 

```R
plotPhyloCore(div_times = c("Median", "Estimated"))

```
</br>

| Arguments  |  |
| :---  | :---  |
| dataset  | Indicates which data set to use. Can be either `"Brawand"` (mammalian) or `"DevSeq"` (angiosperm) data. |

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
