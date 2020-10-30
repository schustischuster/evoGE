
## Perform comparative analysis of DevSeq and Brawand data

This code allows to reproduce the inter-organ intra-species, inter-species and cross-kingdom analyses of the DevSeq and Brawand [(Brawand et al., 2011)](https://pubmed.ncbi.nlm.nih.gov/22012392/) data sets. 


## Contents

* [Getting Started](#getting-started)
  * [Required Packages](#required-packages)
  * [Data input](#data-input)
* [Data analysis and vizualization](#data-analysis-and-vizualization)
* [Session info](#session-info)

---
## Getting started


### Required Packages
Install and load the following R packages before running the reproducible scripts:

```R
# Create list of required packages
lib_List <- c("dplyr", "gplots", "ggplot2", "factoextra", "dendextend", "ggbeeswarm", "mblm", "lsmeans", "rcompanion", "devtools")

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
in_dir <- file.path("compAnalysis", "data")
out_dir <- file.path("compAnalysis")

source(file.path("compAnalysis", "R", "makeCompAnalysis.R"))
source(file.path("compAnalysis", "R", "getATDiv.R"))

```
---
## Data analysis and vizualization

The following function will load and analyze the DevSeq and Brawand ortholog expression data and generate the plots: 

```R
makeCompAnylsis(dataset = c("Brawand", "DevSeq"), expr_estimation = c("TPM", "counts"), 
	        coefficient = c("pearson", "spearman"), devseq_spec = c("Brassicaceae", "all"), 
                data_norm = c("intra-organ", "inter-organ"))

```
</br>

| Arguments  |  |
| :---  | :---  |
| dataset  | Indicates which data set to use. Can be either `"Brawand"` (mammalian) or `"DevSeq"` (angiosperm) data. |
| expr_estimation  | The expression estimation measure; Must be one of `"TPM"` or `"counts"` (VST). |
| coefficient  | A character string that defines which correlation coefficient will be used; Can be either `"pearson"` or `"spearman"`. |
| devseq_spec  | Use one of the two string options: `"Brassicaceae"` for Brassicaceae-specific analysis of the DevSeq data, or `"all"` to perform analysis on all DevSeq angiosperm species (7). |
| data_norm  | This argument indicates how the input data was normalized. For `"intra-organ"`, expression data will be loaded that was normalized between comparative organs across species, whereas for `"inter-organ"` data was normalized between organs AND species. |

</br>

For re-analysis of the mammalian data set, biological replicates that showed a sample correlation below 0.85 (Pearson's r) were excluded. To reproduce the results of this study, execute the following function calls:

```R
makeCompAnylsis(dataset="DevSeq", expr_estimation="TPM", coefficient="pearson", spec="Brassicaeae", data_norm="inter-organ")
makeCompAnylsis(dataset="DevSeq", expr_estimation="TPM", coefficient="pearson", spec="all", data_norm="inter-organ")
makeCompAnylsis(dataset="DevSeq", expr_estimation="counts", coefficient="pearson", spec="Brassicaeae", data_norm="inter-organ")
makeCompAnylsis(dataset="DevSeq", expr_estimation="counts", coefficient="pearson", spec="all", data_norm="inter-organ")
makeCompAnylsis(dataset="Brawand", expr_estimation="counts", coefficient="pearson", data_norm="inter-organ")

```

If pearson expression correlation is chosen, it will compute both the metric pearson distance and the expression distance that is based on the stationary Ohrenstein-Ulenbeck model with variable µ-distance. 
To format the ortholog gene expression tables for correct parsing in treeExp2, execute the following command. The results will be saved in ./compAnalysis/output/data.

```R
getTaxoInput()

```

The following function will compare the gene expression divergence rates between Angiosperms (DevSeq data set) and Mammals (Brawand data set): 

```R
getATDiv <- function(expr_estimation = c("TPM", "counts"), coefficient = c("pearson", "spearman"))

```

The arguments of this function are described above. To reproduce the results of this study, execute the following function calls:

```R
getATDiv(expr_estimation = "TPM", coefficient = "pearson")
getATDiv(expr_estimation = "counts", coefficient = "pearson")

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
#> [1] devtools_1.13.4   rcompanion_1.11.1 lsmeans_2.27-61   mblm_0.12         ggbeeswarm_0.6.0
#> [6] dendextend_1.12.0 factoextra_1.0.5  ggplot2_2.2.1     gplots_3.0.1.1    dplyr_0.7.4

#> loaded via a namespace (and not attached):
#>  [1] nlme_3.1-131         bitops_1.0-6         pbkrtest_0.4-7       ordinal_2015.6-28   
#>  [5] tools_3.3.3          R6_2.4.1             vegan_2.4-5          KernSmooth_2.23-15  
#>  [9] vipor_0.4.5          nortest_1.0-4        lazyeval_0.2.1       mgcv_1.8-17         
#> [13] colorspace_1.3-2     permute_0.9-4        ade4_1.7-10          nnet_7.3-12         
#> [17] withr_2.1.2          gridExtra_2.3        quantreg_5.34        hermite_1.1.1       
#> [21] SparseM_1.77         expm_0.999-2         sandwich_2.5-1       caTools_1.17.1      
#> [25] scales_0.5.0         lmtest_0.9-35        mvtnorm_1.0-6        mc2d_0.1-18         
#> [29] multcompView_0.1-7   digest_0.6.13        foreign_0.8-67       minqa_1.2.4         
#> [33] WRS2_0.9-2           pkgconfig_2.0.3      lme4_1.1-15          manipulate_1.0.1    
#> [37] rlang_0.1.6          bindr_0.1.1          zoo_1.8-1            gtools_3.5.0        
#> [41] car_2.1-6            magrittr_1.5         modeltools_0.2-22    Matrix_1.2-8        
#> [45] Rcpp_0.12.14         DescTools_0.99.23    munsell_0.5.0        viridis_0.5.1       
#> [49] ucminf_1.1-4         multcomp_1.4-8       MASS_7.3-45          plyr_1.8.4          
#> [53] grid_3.3.3           parallel_3.3.3       gdata_2.18.0         ggrepel_0.7.0       
#> [57] lattice_0.20-34      splines_3.3.3        EMT_1.1              boot_1.3-18         
#> [61] estimability_1.2     codetools_0.2-15     stats4_3.3.3         glue_1.2.0          
#> [65] nloptr_1.0.4         miscTools_0.6-22     MatrixModels_0.4-1   gtable_0.3.0        
#> [69] reshape_0.8.7        assertthat_0.2.1     coin_1.2-2           xtable_1.8-4        
#> [73] e1071_1.6-8          coda_0.19-1          class_7.3-14         survival_2.40-1     
#> [77] BSDA_1.2.0           viridisLite_0.3.0    tibble_1.3.4         beeswarm_0.2.3      
#> [81] memoise_1.1.0        bindrcpp_0.2         cluster_2.0.5        maxLik_1.3-4        
#> [85] TH.data_1.0-10       RVAideMemoire_0.9-69


```
