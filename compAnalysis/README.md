
## Perform comparative analysis of DevSeq and Brawand data

This code allows to reproduce the inter-organ intra-species, inter-species and cross-kingdom analyses of the DevSeq and Brawand [(Brawand et al., 2011)](https://pubmed.ncbi.nlm.nih.gov/22012392/) data sets. 


## Contents

* [Getting Started](#getting-started)
  * [Required Packages](#required-packages)
  * [Data input](#data-input)
* [Data analysis and vizualization](#data-analysis-and-vizualization)
* [Session info](#session-info)


## Getting started


### Required Packages
Install and load the following R packages before running the reproducible scripts:

```R
if (!require(dplyr)) install.packages('dplyr')
library(dplyr)
if (!require(gplots)) install.packages('gplots')
library(gplots)
if (!require(ggplot2)) install.packages('ggplot2')
library(ggplot2)
if (!require(factoextra)) install.packages('factoextra')
library(factoextra)
if (!require(dendextend)) install.packages('dendextend')
library(dendextend)
if (!require(ggbeeswarm)) install.packages('ggbeeswarm')
library(ggbeeswarm)

```

### Data input
Download the entire subdirectory containing the [data](https://github.com/schustischuster/evoGEx/tree/master/compAnalysis/data) folder and [R script](https://github.com/schustischuster/evoGEx/tree/master/compAnalysis/R) to the working directory on your computer, e.g. by using [GitZip](http://kinolien.github.io/gitzip/), and extract the file. Then, set the path for input and output files and source the R scripts: 

```R
in_dir <- "./compAnalysis/data"
out_dir <- "./compAnalysis"

source("compAnalysis/R/makeCompAnalysis.R")

```

## Data analysis and vizualization

The following function will load and analyze the DevSeq and Brawand ortholog expression data and generate the plots. 

```R
makeCompAnylsis(dataset = c("Brawand", "DevSeq"), expr_estimation = c("TPM", "counts"), 
                coefficient = c("pearson", "spearman"), devseq_spec = c("Brassicaceae", "all"))

```
To reproduce the results of this study, execute the following function calls:

```R
makeCompAnylsis(dataset = "Brawand", expr_estimation = "TPM", coefficient = "spearman")
makeCompAnylsis(dataset = "DevSeq", expr_estimation = "TPM", coefficient = "spearman", devseq_spec = "Brassicaceae")
makeCompAnylsis(dataset = "DevSeq", expr_estimation = "TPM", coefficient = "spearman", devseq_spec = "all")
makeCompAnylsis(dataset = "DevSeq", expr_estimation = "TPM", coefficient = "pearson", devseq_spec = "Brassicaceae")
makeCompAnylsis(dataset = "DevSeq", expr_estimation = "TPM", coefficient = "pearson", devseq_spec = "all")

```

This will generate the panels for the following figures:


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
#> [1] stats     graphics  grDevices utils     datasets  methods   base

#> other attached packages:
#> [1] ggbeeswarm_0.6.0  dendextend_1.12.0 factoextra_1.0.5  ggplot2_2.2.1     gplots_3.0.1.1   
#> [6] dplyr_0.7.4      

#> loaded via a namespace (and not attached):
#>  [1] Rcpp_0.12.14       bindr_0.1.1        magrittr_1.5       munsell_0.5.0     
#>  [5] viridisLite_0.3.0  colorspace_1.3-2   R6_2.4.1           rlang_0.1.6       
#>  [9] vipor_0.4.5        plyr_1.8.4         caTools_1.17.1     grid_3.3.3        
#> [13] beeswarm_0.2.3     gtable_0.3.0       KernSmooth_2.23-15 gtools_3.5.0      
#> [17] lazyeval_0.2.1     assertthat_0.2.1   tibble_1.3.4       bindrcpp_0.2      
#> [21] gridExtra_2.3      viridis_0.5.1      bitops_1.0-6       ggrepel_0.7.0     
#> [25] glue_1.2.0         gdata_2.18.0       scales_0.5.0       pkgconfig_2.0.3

```
