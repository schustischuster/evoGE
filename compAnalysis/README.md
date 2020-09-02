
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
if (!require(mblm)) install.packages('mblm')
library(mblm)
if (!require(lsmeans)) install.packages('lsmeans')
library(lsmeans)
if (!require(rcompanion)) install.packages('rcompanion')
library(rcompanion)

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
	        coefficient = c("pearson", "spearman"), devseq_spec = c("Brassicaceae", "all"), 
                data_norm = c("intra-organ", "inter-organ"))

```

| Arguments  | Description |
| ------------- | ------------- |
| dataset  | Indicates which data set to use. Can be either `"Brawand"` or `"DevSeq"`. |
| expr_estimation  | The expression estimation measure; Must be one of `"TPM"` or `"counts"` (VST). |
| coefficient  | A character string that defines which correlation coefficient will be used; Can be either `"pearson"` or `"spearman"`. |
| devseq_spec  | Use one of the two string options: `"Brassicaceae"` for Brassicaceae-specific analysis of the DevSeq data, or `"all"` to perform analysis on all DevSeq angiosperm species (7). |
| data_norm  | This argument indicates the normalization method that was used for RNA-Seq data normalization; In case of `"intra-organ"`, data was normalized within comparative organs across species, whereas for `"inter-organ"` data was normalized between organs and species. |


To reproduce the results of this study, execute the following function calls:

```R
makeCompAnylsis(dataset="DevSeq", expr_estimation="TPM", coefficient="pearson", spec="Brassicaeae", data_norm="inter-organ")
makeCompAnylsis(dataset="DevSeq", expr_estimation="TPM", coefficient="pearson", spec="all", data_norm="inter-organ")
makeCompAnylsis(dataset="DevSeq", expr_estimation="counts", coefficient="pearson", spec="Brassicaeae", data_norm="inter-organ")
makeCompAnylsis(dataset="DevSeq", expr_estimation="counts", coefficient="pearson", spec="all", data_norm="inter-organ")
makeCompAnylsis(dataset="Brawand", expr_estimation="counts", coefficient="pearson", data_norm="inter-organ")

```
The following function will compare the gene expression divergence times between Angiosperms (DevSeq data set) and Vertebrates (Brawand data set). 

```R
getATDiv <- function(expr_estimation = c("TPM", "counts"), coefficient = c("pearson", "spearman"))

```
To reproduce the results of this study, execute the following function calls:

```R
getATDiv(expr_estimation = "TPM", coefficient = "pearson")
getATDiv(expr_estimation = "counts", coefficient = "pearson")

```
These function calls will generate the panels for the following figures:


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
#> [1] rcompanion_1.11.1 lsmeans_2.27-61   mblm_0.12         ggbeeswarm_0.6.0  dendextend_1.12.0
#> [6] factoextra_1.0.5  ggplot2_2.2.1     gplots_3.0.1.1    dplyr_0.7.4            

#> loaded via a namespace (and not attached):
#>  [1] viridis_0.5.1        viridisLite_0.3.0    splines_3.3.3        BSDA_1.2.0          
#>  [5] gtools_3.5.0         ucminf_1.1-4         assertthat_0.2.1     expm_0.999-2        
#>  [9] stats4_3.3.3         coin_1.2-2           vipor_0.4.5          ggrepel_0.7.0       
#> [13] lattice_0.20-34      quantreg_5.34        glue_1.2.0           minqa_1.2.4         
#> [17] colorspace_1.3-2     sandwich_2.5-1       Matrix_1.2-8         plyr_1.8.4          
#> [21] pkgconfig_2.0.3      SparseM_1.77         EMT_1.1              xtable_1.8-4        
#> [25] mvtnorm_1.0-6        scales_0.5.0         gdata_2.18.0         manipulate_1.0.1    
#> [29] lme4_1.1-15          MatrixModels_0.4-1   tibble_1.3.4         mgcv_1.8-17         
#> [33] car_2.1-6            TH.data_1.0-10       maxLik_1.3-4         nnet_7.3-12         
#> [37] lazyeval_0.2.1       pbkrtest_0.4-7       survival_2.40-1      magrittr_1.5        
#> [41] ordinal_2015.6-28    estimability_1.2     nlme_3.1-131         MASS_7.3-45         
#> [45] WRS2_0.9-2           RVAideMemoire_0.9-69 foreign_0.8-67       class_7.3-14        
#> [49] beeswarm_0.2.3       vegan_2.4-5          tools_3.3.3          multcomp_1.4-8      
#> [53] munsell_0.5.0        cluster_2.0.5        bindrcpp_0.2         ade4_1.7-10         
#> [57] e1071_1.6-8          multcompView_0.1-7   caTools_1.17.1       rlang_0.1.6         
#> [61] grid_3.3.3           nloptr_1.0.4         miscTools_0.6-22     hermite_1.1.1       
#> [65] bitops_1.0-6         boot_1.3-18          DescTools_0.99.23    gtable_0.3.0        
#> [69] codetools_0.2-15     reshape_0.8.7        R6_2.4.1             gridExtra_2.3       
#> [73] zoo_1.8-1            mc2d_0.1-18          nortest_1.0-4        bindr_0.1.1         
#> [77] KernSmooth_2.23-15   permute_0.9-4        modeltools_0.2-22    parallel_3.3.3      
#> [81] Rcpp_0.12.14         lmtest_0.9-35        coda_0.19-1     

```
