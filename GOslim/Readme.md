## Analyse rates of gene expression evolution of functional groups

This code allows to analyse the rate of gene expression evolution for functionally related genes. This involves estimating the stability of correlations using a Monte-Carlo simulation, matching optimal number of control gene sets to each GO slim term category under the assessment of balance statistics, evaluating the strength of expression conservation in relation to gene expression levels, and applying non-linear regression models to estimate the rates of gene expression evolution for functionally related genes. 


## Contents

* [Getting Started](#getting-started)
  * [Required Packages](#required-packages)
  * [Data input](#data-input)
* [Data analysis and visualization](#data-analysis-and-visualization)
  * [Stability of correlations](#stability-of-correlations)
  * [Relationship between gene expression level and evolutionary conservation](#relationship-between-gene-expression-level-and-evolutionary-conservation)
  * [Rate of expression evolution](#rate-of-expression-evolution)
  * [Generate plots of GO analysis](#generate-plots-of-go-analysis)
  * [Coefficient of variation](#coefficient-of-variation)
  * [Expression levels across functional groups](#expression-levels-across-functional-groups)
* [Session info](#session-info)

---
## Getting started


### Required Packages
Install and load the following R packages before running the reproducible scripts:

```R
# Required packages
lib_List <- c("dplyr", "MatchIt", "gplots", "ggplot2", "scales")

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
in_dir <- file.path("evoGE-master", "GOslim", "data")
out_dir <- file.path("evoGE-master", "GOslim")
path_to_R_files <- file.path("evoGE-master", "GOslim", "R")

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

### Stability of correlations

Sample correlations converge to the population value with increasing sample size, and it has been shown that the sample size should approach n=250 for stable estimates ([Schönbrodt and Perugini, J Res Pers. 2013](https://www.sciencedirect.com/science/article/abs/pii/S0092656613000858)). Since the points of stability (POS) published in this work were based on normal distributions, we implemented a Monte-Carlo simulation to retrieve sample size estimates for stable correlations derived from gene expression data.


* `estimatePOS(nbootstrap, coswidth, clevel, ...)`


| Arguments  |  |
| :---  | :---  |
| nbootstrap  | Number of sampling trajectories. |
| coswidth  | Indicates the corridor of stability (COS). Values of 0.1/0.15/0.2 correspond to small/medium/large effect sizes, respectively. For more details, see [Cohen, Psychol Bull (1992)](https://pubmed.ncbi.nlm.nih.gov/19565683/). |
| clevel  | Indicates the level of confidence for the confidence interval. |


To reproduce the results of this study, execute the following function call:

```R
estimatePOS(nbootstrap = 1000, coswidth = 0.1, clevel = 0.8)

```
### Relationship between gene expression level and evolutionary conservation

Next, we wanted to test whether the evolutionary stability of gene subsets is affected by gene expression levels. We therefore generated subsets of orthologous genes according to quantiles of average expression either across all samples (inter-organ inter-species), or within the same organ across species (intra-organ inter-species). For each quantile set, we then calculated metric Pearson distances and fitted non-linear regression models to estimate quantile-specific rates of gene expression evolution.


* `getExprCons(nquant, qtype, ...)`


| Arguments  |  |
| :---  | :---  |
| nquant  | Number of genes in each quantile. |
| qtype  | Type of quantile. Use "base_mean" for quantiles of average gene expression across organs and species, and "organ_spec" for organ-specific (intra-organ) quantiles of average gene expression across species. |


To reproduce the results of this study, execute the following function calls:

```R
getExprCons(nquant = 500, qtype = "base_mean")
getExprCons(nquant = 500, qtype = "organ_spec")

```
### Rate of expression evolution

Now, multiple control genes will be matched to each gene of a GO slim category that is larger than the size threshold (POS) estimated above. The optimal number of control sets will be determined using balance statistics (standardized mean difference and variance ratio). Subsequently, intra-organ distances will be calculated for all species pairs and gene sets, and non-linear regression model will be fitted to the data. Finally, the regression slopes of treatment and control groups of each functional category will be compared using nonparametric statistics (Wilcoxon rank-sum test, permutation test).


* `getGOSLIM(aspect, sample_size)`


| Arguments  |  |
| :---  | :---  |
| aspect  | GO slim term aspect. Can be either "biological_process" or "molecular_function" |
| sample_size  | Indicates the minimum number of genes required in a GO slim category. This is the point of stability (POS) determined in the previous step. |


To reproduce the results of this study, execute the following function calls:

```R
getGOSLIM(aspect = "biological_process", sample_size = 412)
getGOSLIM(aspect = "molecular_function", sample_size = 412)

```
### Generate plots of GO analysis

To visualize the results of the GO enrichment analysis for genes of the first and last expression quantile, and to plot the test statistics for the relative rates of expression evolution of genes belonging to predefined functional groups, execute the following function (make sure to run getExprCons and getGOSLIM functions first):

```R
plotGOs()

```
### Coefficient of variation

Every gene was classified as evolutionary stable or variable using the coefficient of variation (CV). This coefficient was calculated independently for each organ across species, and then averaged, resulting in a mean coefficient of variation for each gene. Stable and variable genes were matched based on their mean expression level across samples using the "nearest" method and caliper option. Then, the fraction of stable and variable genes was assessed in each functional category. Each category was tested for a significant increase in the fraction of stable or variable genes using a Chi-squared test.

* `getCV(aspect, estimate, sample_size)`

| Arguments  |  |
| :---  | :---  |
| aspect  | GO slim term aspect. Can be either "biological_process" or "molecular_function" |
| estimate  | Expression estimate. Use "VST" for Variance Stabilization Transformed counts, or "TPM" for Transcripts Per Million |
| sample_size  | Indicates the minimum number of genes required in a GO slim category. This is the point of stability (POS) determined in estimatePOS(). |


The following function call was used to generate the results of this study:

```R
getCV(aspect = "biological_process", estimate = "VST", sample_size = 412)

```

### Expression levels across functional groups

Finally, the mean expression level of each gene across organs and species was calculated, and their distribution was visualized for all GO slim categories that showed a significant different proportion of evolutionarily stable and variable genes. Expression values are log2-transformed Transcripts Per Million (TPM).

* `plotGroupEx(sample_size, ...)`

| Arguments  |  |
| :---  | :---  |
| sample_size  | Indicates the minimum number of genes required in a GO slim category. This is the point of stability (POS) determined in estimatePOS(). |

```R
plotGroupEx(sample_size = 412)

```

---
## Session info

This code was developed and tested on MacOS X 10.9.5 in R version 3.3.3 and on MacOS 12.3.1 in R version 4.1.3.

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
#> [1] gplots_3.0.1.1     MatchIt_4.2.0     scales_0.5.0     ggplot2_2.2.1     dplyr_0.7.4

#> loaded via a namespace (and not attached):
#> [1] Rcpp_1.0.7         gtools_3.5.0       assertthat_0.2.1      bitops_1.0-6      
#> [5] grid_3.3.3         R6_2.4.1           plyr_1.8.4            backports_1.2.1   
#> [9] gtable_0.3.0       magrittr_1.5       KernSmooth_2.23-15    rlang_0.1.6       
#>[13] lazyeval_0.2.1     gdata_2.18.0       bindrcpp_0.2          tools_3.3.3       
#>[17] glue_1.2.0         munsell_0.5.0      pkgconfig_2.0.3       colorspace_1.3-2  
#>[21] caTools_1.17.1     bindr_0.1.1        tibble_1.3.4  

```

```R
#> R version 4.1.3 (2022-03-10)
#> Platform: x86_64-apple-darwin17.0 (64-bit)
#> Running under: macOS Monterey 12.3.1

#> Matrix products: default
#> BLAS:   /Library/Frameworks/R.framework/Versions/4.1/Resources/lib/libRblas.0.dylib
#> LAPACK: /Library/Frameworks/R.framework/Versions/4.1/Resources/lib/libRlapack.dylib

#> locale:
#> [1] de_DE.UTF-8/de_DE.UTF-8/de_DE.UTF-8/C/de_DE.UTF-8/de_DE.UTF-8

#> attached base packages:
#> [1] stats     graphics     grDevices     utils     datasets     methods     base     

#> other attached packages:
#> [1] scales_1.1.1     ggplot2_3.3.5     gplots_3.1.1     MatchIt_4.3.4     dplyr_1.0.8  

#> loaded via a namespace (and not attached):
#>  [1] Rcpp_1.0.8.2        pillar_1.7.0        compiler_4.1.3      RColorBrewer_1.1-2
#>  [5] bitops_1.0-7        tools_4.1.3         digest_0.6.29       lifecycle_1.0.1   
#>  [9] tibble_3.1.6        gtable_0.3.0        nlme_3.1-155        lattice_0.20-45   
#> [13] mgcv_1.8-39         pkgconfig_2.0.3     rlang_1.0.2         Matrix_1.4-0      
#> [17] cli_3.2.0           withr_2.5.0         generics_0.1.2      vctrs_0.3.8       
#> [21] gtools_3.9.2        caTools_1.18.2      grid_4.1.3          tidyselect_1.1.2  
#> [25] glue_1.6.2          R6_2.5.1            fansi_1.0.2         purrr_0.3.4       
#> [29] farver_2.1.0        magrittr_2.0.2      backports_1.4.1     ellipsis_0.3.2    
#> [33] splines_4.1.3       colorspace_2.0-3    labeling_0.4.2      utf8_1.2.2        
#> [37] KernSmooth_2.23-20  munsell_0.5.0       crayon_1.5.0

```







