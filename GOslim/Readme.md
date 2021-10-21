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


* estimatePOS(nbootstrap, coswidth, bss, ...)


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

Next, we wanted to test wether the evolutionary stability of gene subsets is affected by gene expression levels. We therefore generated subsets of orthologous genes according to quantiles of average expression either across all samples (inter-organ inter-species), or within the same organ across species (intra-organ inter-species). For each quantile set, we then calculated metric pearson distances and fitted non-linear regression models to estimate quantile-specific rates of gene expression evolution.


* getExprCons(nquant, qtype, ...)


</br>

| Arguments  |  |
| :---  | :---  |
| nquant  | Number of genes in each quantile. |
| qtype  | Type of quantile. Use "base_mean" for quantiles of average gene expression across organs and species, and "organ_spec" for organ-specific (intra-organ) quantiles of average gene expression across species. |

</br>

To reproduce the results of this study, execute the following function calls:

```R
getExprCons(nquant = 500, qtype = "base_mean")
getExprCons(nquant = 500, qtype = "organ_spec")

```
Now, multiple control genes will be matched to each gene of a GO slim category that is larger than the size threshold (POS) estimated above. The optimal number of control sets will be determined using balance statistics (standardized mean difference and variance ratio). Subsequently, intra-organ distances will be calculated for all species pairs and gene sets, and non-linear regression model will be fitted to the data. Finally, the regression slopes of treatment and control groups of each functional category will be compared using nonparametric statistics (Wilcoxon rank-sum test, permutation test).


* getGOSLIM(aspect, sample_size)


</br>

| Arguments  |  |
| :---  | :---  |
| aspect  | GOslim term aspect. Can be either "biological_process" or "molecular_function" |
| sample_size  | Indicates the minimum number of genes required in a GO slim category. This is the point of stability (POS) determined in the previous step. |

</br>

To reproduce the results of this study, execute the following function calls:

```R
getGOSLIM(aspect = "biological_process", sample_size = 293)
getGOSLIM(aspect = "molecular_function", sample_size = 293)

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
#> [1] stats    graphics    grDevices    utils    datasets    methods    base    

#> other attached packages:
#> [1] gplots_3.0.1.1    MatchIt_4.2.0    scales_0.5.0    ggplot2_2.2.1    dplyr_0.7.4

#> loaded via a namespace (and not attached):
#> [1] Rcpp_1.0.7         gtools_3.5.0       assertthat_0.2.1      bitops_1.0-6      
#> [5] grid_3.3.3         R6_2.4.1           plyr_1.8.4            backports_1.2.1   
#> [9] gtable_0.3.0       magrittr_1.5       KernSmooth_2.23-15    rlang_0.1.6       
#>[13] lazyeval_0.2.1     gdata_2.18.0       bindrcpp_0.2          tools_3.3.3       
#>[17] glue_1.2.0         munsell_0.5.0      pkgconfig_2.0.3       colorspace_1.3-2  
#>[21] caTools_1.17.1     bindr_0.1.1        tibble_1.3.4  

```







