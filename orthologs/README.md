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
lib_List <- c("dplyr", "ggplot2", "ape", "scales", "gtable", "ggtree")

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
in_dir <- file.path("evoGE", "orthologs", "data")
out_dir <- file.path("evoGE", "orthologs")

source(file.path("evoGE", "orthologs", "R", "plotPhyloCore.R"))

```
---
## Data analysis and vizualization

The following function will generate the angiosperm phylogeny and orthologous gene plots:

* `plotPhyloCore(div_times = c("Median", "Estimated"))`

</br>

| Arguments  |  |
| :---  | :---  |
| div_times  | Indicates which pairwise divergence time to use. Must be either `"Median"` (median time derived from all studies) or `"Estimated"` (TTOL estimation). For more details, see [Hedges et al., MBE (2015)](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4379413/). |

</br>

To reproduce the results of this study, execute the following function call:

```R
plotPhyloCore(div_times = "Estimated")

```

The following function will generate simplified phylogenetic trees for angiosperms and mammals: 

* `plotPhyloComp(div_times = c("Median", "Estimated"))`

This function takes the same arguments described above. TTOL estimation (`"Estimated"`) was used to generate the plots:

```R
plotPhyloComp(div_times = "Estimated")

```


---
## Session info

```R
sessionInfo()
```

```R
#> R version 3.3.3 (2017-03-06)
#> Platform: x86_64-apple-darwin13.4.0 (64-bit)
#> Running under: OS X Mavericks 10.9.5

#> locale:
#> [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

#> attached base packages:
#> [1] stats     graphics     grDevices utils     datasets     methods     base   

#> other attached packages:
#> [1] gtable_0.3.0     scales_0.5.0     ape_5.0     ggplot2_2.2.1     dplyr_0.7.4 

#> loaded via a namespace (and not attached):
#> [1] Rcpp_1.0.7         lattice_0.20-34  assertthat_0.2.1  grid_3.3.3      R6_2.4.1 
#> [6] plyr_1.8.4         nlme_3.1-131     magrittr_1.5      rlang_0.1.6     lazyeval_0.2.1    
#> [11] bindrcpp_0.2      glue_1.2.0       munsell_0.5.0     parallel_3.3.3  pkgconfig_2.0.3    
#> [16] colorspace_1.3-2  bindr_0.1.1      tibble_1.3.4     
 
```
