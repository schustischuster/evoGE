## Make angiosperm phylogeny and ortholog gene plots

This code allows to build the angiosperm and mammalian phylogenies based on divergence times obtained from [TimeTree](http://www.timetree.org/), and to generate plots showing the total number of coding/non-coding orthologous genes for the pairwise species comparisons. 


## Contents

* [Getting Started](#getting-started)
  * [Required Packages](#required-packages)
  * [Data input](#data-input)
* [Data analysis and visualization](#data-analysis-and-visualization)
* [Session info](#session-info)

---
## Getting started


### Required Packages
Install Bioconductor core packages and ggtree:

```R
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install()
BiocManager::install(ggtree)

```

Install and load the following R packages before running the reproducible scripts:

```R
# Required packages
lib_List <- c("dplyr", "ggplot2", "ape", "scales", "gtable", "ggtree")

# Install missing packages
instpack <- lib_List %in% installed.packages()[,"Package"]
if (any(instpack == FALSE)) {
  install.packages(lib_List[!instpack])
}

# Load packages
invisible(lapply(lib_List, library, character.only = TRUE))

```

### Data input
Download and extract the evoGE repository to the working directory on your computer. Then, set the path for input and output files and source the R scripts: 

```R
in_dir <- file.path("evoGE-master", "orthologs", "data")
out_dir <- file.path("evoGE-master", "orthologs")

source(file.path("evoGE-master", "orthologs", "R", "plotPhyloCore.R"))
source(file.path("evoGE-master", "orthologs", "R", "plotPhyloComp.R"))

```
---
## Data analysis and visualization

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
#> [1] stats     graphics     grDevices     utils     datasets     methods     base   

#> other attached packages:
#> [1] ggtree_1.6.11     gtable_0.3.0     scales_0.5.0     ape_5.0     ggplot2_2.2.1     dplyr_0.7.4 

#> loaded via a namespace (and not attached):
#> [1] Rcpp_1.0.7       lattice_0.20-34     tidyr_0.7.2        assertthat_0.2.1     grid_3.3.3      
#> [6] R6_2.4.1         plyr_1.8.4          jsonlite_1.5       nlme_3.1-131         stats4_3.3.3    
#>[11] magrittr_1.5     rlang_0.1.6         lazyeval_0.2.1     bindrcpp_0.2         glue_1.2.0      
#>[16] purrr_0.2.4      munsell_0.5.0       parallel_3.3.3     pkgconfig_2.0.3      colorspace_1.3-2
#>[21] bindr_0.1.1      tibble_1.3.4     
 
```

```R
#> R version 4.1.3 (2022-03-10)
#> Platform: x86_64-apple-darwin17.0 (64-bit)
#> Running under: macOS Monterey 12.3

#> Matrix products: default
#> BLAS:   /Library/Frameworks/R.framework/Versions/4.1/Resources/lib/libRblas.0.dylib
#> LAPACK: /Library/Frameworks/R.framework/Versions/4.1/Resources/lib/libRlapack.dylib

#> locale:
#> [1] de_DE.UTF-8/de_DE.UTF-8/de_DE.UTF-8/C/de_DE.UTF-8/de_DE.UTF-8

#> attached base packages:
#> [1] stats     graphics     grDevices     utils     datasets     methods     base     

#> other attached packages:
#> [1] ggtree_3.2.1        gtable_0.3.0        scales_1.1.1        ape_5.6-2          
#> [5] ggplot2_3.3.5       dplyr_1.0.8         BiocManager_1.30.16

#> loaded via a namespace (and not attached):
#> [1]  Rcpp_1.0.8.2       pillar_1.7.0       compiler_4.1.3     yulab.utils_0.0.4  tools_4.1.3       
#> [6]  aplot_0.1.2        jsonlite_1.8.0     tidytree_0.3.9     lifecycle_1.0.1    tibble_3.1.6      
#> [11] nlme_3.1-155       lattice_0.20-45    pkgconfig_2.0.3    rlang_1.0.2        cli_3.2.0         
#> [16] ggplotify_0.1.0    patchwork_1.1.1    parallel_4.1.3     treeio_1.18.1      withr_2.5.0       
#> [21] generics_0.1.2     vctrs_0.3.8        gridGraphics_0.5-1 grid_4.1.3         tidyselect_1.1.2  
#> [26] glue_1.6.2         R6_2.5.1           fansi_1.0.2        purrr_0.3.4        tidyr_1.2.0       
#> [31] magrittr_2.0.2     ellipsis_0.3.2     colorspace_2.0-3   utf8_1.2.2         lazyeval_0.2.2    
#> [36] munsell_0.5.0      ggfun_0.0.5        crayon_1.5.0 

```
