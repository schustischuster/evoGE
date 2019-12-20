## Comparative Analysis of ATGE and DevSeq data sets

This code visualizes the results of the AtGenExpress (ATGE) and DevSeq comparative analysis, which can be found [here](https://github.com/schustischuster/ATGE-DevSeq).


## Contents

* [Getting Started](#getting-started)
  * [Required Packages](#required-packages)
  * [Data input](#data-input)
* [Visualization](#visualization)
* [Session info](#session-info)


## Getting started


### Required Packages
Install and load the following R packages before running the reproducible scripts:

```R
if (!require(dplyr)) install.packages('dplyr')
library(dplyr)
if (!require(gplots)) install.packages('gplots')
library(gplots)
if (!require(factoextra)) install.packages('factoextra')
library(factoextra)
if (!require(dendextend)) install.packages('dendextend')
library(dendextend)

```
  
### Data input
Download the [data](https://github.com/schustischuster/evoGEx/tree/master/ATGE-DevSeq/data) folder and R script to the working directory on your computer. Then, set the file path for input and output files and source the script: 

```R
in_dir <- "./data"
out_dir <- "."

# Store plots in /out_dir/output/plots
if (!dir.exists(file.path(out_dir, "output", "plots"))) 
  dir.create(file.path(out_dir, "output", "plots"), recursive = TRUE)

source("DevSeq_ATGE_plots.R")

```

## Visualization

After loading the data and sourcing the R script, run the following commands to generate the plots:

```R
# Boxplot showing pairwise ATGE-DevSeq gene correlations 
plot_Gene_Corr(spearman_RE, pearson_RE, pearson_log2_RE)

# Boxplot showing ATGE-DevSeq sample correlations
plot_Sample_Corr(atge_devseq_spearman, atge_devseq_pearson, atge_devseq_log_pearson)

# Correlation heatmap of merged ATGE-DevSeq data
makeCorrplot(atge_devseq_re, coefficient = "pearson")
makeCorrplot(atge_devseq_re_log, coefficient = "pearson")

# hclust dendrogram of ATGE and DevSeq data
makeDendrogram(atge_re, coefficient = "pearson")
makeDendrogram(devseq_re, coefficient = "pearson")

```

The plotting functions will generate the panels for the following figure:


![ATGE-DevSeq](README_files/ATGE-DevSeq.png)


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
#> [1] parallel  stats4    stats     graphics  grDevices utils     datasets  methods   base    

#> other attached packages:
#> [1] bindrcpp_0.2         rtracklayer_1.34.2   GenomicRanges_1.26.4 GenomeInfoDb_1.10.3  IRanges_2.8.2       
#> [6] S4Vectors_0.12.2     BiocGenerics_0.20.0  dplyr_0.7.4    

#> loaded via a namespace (and not attached):
#> [1] Rcpp_0.12.14               bindr_0.1.1                XVector_0.14.1             magrittr_1.5              
#> [5] zlibbioc_1.20.0            GenomicAlignments_1.10.1   BiocParallel_1.8.2         lattice_0.20-34           
#> [9] R6_2.4.1                   rlang_0.1.6                tools_3.3.3                grid_3.3.3                
#> [13] SummarizedExperiment_1.4.0 Biobase_2.34.0             assertthat_0.2.1           tibble_1.3.4              
#> [17] Matrix_1.2-8               bitops_1.0-6               RCurl_1.95-4.10            glue_1.2.0                
#> [21] Biostrings_2.42.1          Rsamtools_1.26.2           XML_3.98-1.9               pkgconfig_2.0.3  

```
