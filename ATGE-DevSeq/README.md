## Comparative Analysis of ATGE and DevSeq data sets

This code visualizes the results of the AtGenExpress (ATGE) and DevSeq comparative analysis, which can be found [here](https://github.com/schustischuster/ATGE-DevSeq).


## Contents

* [Getting Started](#getting-started)
  * [Required Packages](#required-packages)
  * [Data input](#data-input)
* [Data visualization](#data-visualization)
* [Session info](#session-info)


## Getting started


### Required Packages
Install and load the following R packages before running the reproducible script:

```R

# Required packages
lib_List <- c("dplyr", "gplots", "factoextra", "dendextend")

# Install missing packages
instpack <- lib_List %in% installed.packages()[,"Package"]
if (any(instpack == FALSE)) {
  install.packages(lib_List[!instpack])
}

# Load packages
invisible(lapply(lib_List, library, character.only = TRUE))

```
  
### Data input
Download and extract the evoGE repository to the working directory on your computer. Then, set the path for input and output files and source the R script:

```R
in_dir <- file.path("evoGE", "ATGE-DevSeq", "data")
out_dir <- file.path("evoGE", "ATGE-DevSeq")

source(file.path("evoGE", "ATGE-DevSeq", "R", "DevSeq_ATGE_plots.R"))

```

## Data visualization

After loading the data and sourcing the R script, run the following commands to generate the plots:

```R
# List of genes for plotting
genelist <- c("WUS", "REV", "AP1", "AG", "LFY", "PLT1", "FLC", "PIN1")

# Pairwise gene correlation plots
plotRE(exp_data = devseq_log2_re_vs_atge_log2_re, genelist = genelist)

# Correlation heatmap of combined ATGE-DevSeq data
makeCorrplot(exp_data=atge_devseq_re_log, coefficient="pearson", clustm="complete")

# Boxplots of pairwise sample and gene correlations
plotCor(data = gene_sample_cor)

# hclust dendrogram of ATGE and DevSeq data
makeDendrogram(atge, coefficient = "pearson", clustm="complete")
makeDendrogram(devseq, coefficient = "pearson", clustm="complete")

```

---
## Session info

This code was developed and tested on MacOS X 10.9.5 in R version 3.3.3 and on MacOS 12.3 in R version 4.1.3. 

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
#> [1] dendextend_1.12.0 factoextra_1.0.5  ggplot2_2.2.1     gplots_3.0.1.1    dplyr_0.7.4 

#> loaded via a namespace (and not attached):
#> [1] Rcpp_0.12.14       bindr_0.1.1        magrittr_1.5       munsell_0.5.0      viridisLite_0.3.0 
#> [6] colorspace_1.3-2   R6_2.4.1           rlang_0.1.6        plyr_1.8.4         caTools_1.17.1    
#> [11] grid_3.3.3         gtable_0.3.0       KernSmooth_2.23-15 gtools_3.5.0       lazyeval_0.2.1    
#> [16] assertthat_0.2.1   tibble_1.3.4       bindrcpp_0.2       gridExtra_2.3      viridis_0.5.1     
#> [21] bitops_1.0-6       ggrepel_0.7.0      glue_1.2.0         gdata_2.18.0       scales_0.5.0      
#> [26] pkgconfig_2.0.3  

```
