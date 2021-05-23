
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
# Load the required packages
lib_List <- c("dplyr", "gplots", "ggplot2", "factoextra", "dendextend", "ggbeeswarm", "lsmeans", "rcompanion", "scales", 
"matrixStats", "ape", "ggtree")

loadLibrary <- function(x) { 
    if (!require(x, character.only = T)) {
        install.packages('x')
        library(x)
    }
}

# Load packages
invisible(lapply(lib_List, loadLibrary))

```

Install TreeExp2 package for phylogenetic transcriptome analysis through devtools:

```R
install.packages('devtools')
devtools::install_github("jingwyang/TreeExp")

# Load package
library('TreeExp')

```

### Data input
Download and extract the entire directory to the working directory on your computer. Then, set the path for input and output files and source the R scripts: 

```R
in_dir <- file.path("evoGE", "compAnalysis", "data")
out_dir <- file.path("evoGE", "compAnalysis")

source(file.path("evoGE", "compAnalysis", "R", "makeCompAnalysis.R"))
source(file.path("evoGE", "compAnalysis", "R", "getATDiv.R"))

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

To format the ortholog gene expression tables for correct parsing in treeExp2, execute the following command. The results will be stored in ./compAnalysis/output/data.

```R
getTaxoInput()

```

The following function will compare the gene expression divergence rates between Angiosperms (DevSeq data set) and Mammals (Brawand data set). If `pearson` expression correlation is chosen, both metric pearson distance and an expression distance under the stationary Ornstein-Uhlenbeck (OU) model with variable optimal expression level [(Yang et al., 2019)](https://pubmed.ncbi.nlm.nih.gov/31609424/) will be estimated. A detailed description of the variable-µ expression distance including case studies can be found [here](https://jingwyang.github.io/TreeExp-Tutorial/). 

```R
getATDiv(expr_estimation = c("TPM", "counts"), coefficient = c("pearson", "spearman"))

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
#> [1] matrixStats_0.52.2 scales_0.5.0       rcompanion_1.11.1  lsmeans_2.27-61   
#> [5] mblm_0.12          ggbeeswarm_0.6.0   dendextend_1.12.0  factoextra_1.0.5  
#> [9] ggplot2_2.2.1      gplots_3.0.1.1     dplyr_0.7.4       

#> loaded via a namespace (and not attached):
#>  [1] viridis_0.5.1        viridisLite_0.3.0    splines_3.3.3        BSDA_1.2.0 
#>  [5] gtools_3.5.0         ucminf_1.1-4         assertthat_0.2.1     expm_0.999-2        
#>  [9] stats4_3.3.3         coin_1.2-2           vipor_0.4.5          ggrepel_0.7.0       
#> [13] lattice_0.20-34      quantreg_5.34        glue_1.2.0           minqa_1.2.4         
#> [17] colorspace_1.3-2     sandwich_2.5-1       Matrix_1.2-8         plyr_1.8.4          
#> [21] pkgconfig_2.0.3      SparseM_1.77         EMT_1.1              xtable_1.8-4        
#> [25] mvtnorm_1.0-6        gdata_2.18.0         manipulate_1.0.1     lme4_1.1-15         
#> [29] MatrixModels_0.4-1   tibble_1.3.4         mgcv_1.8-17          car_2.1-6           
#> [33] TH.data_1.0-10       maxLik_1.3-4         nnet_7.3-12          lazyeval_0.2.1      
#> [37] pbkrtest_0.4-7       survival_2.40-1      magrittr_1.5         ordinal_2015.6-28   
#> [41] estimability_1.2     nlme_3.1-131         MASS_7.3-45          WRS2_0.9-2          
#> [45] RVAideMemoire_0.9-69 foreign_0.8-67       class_7.3-14         beeswarm_0.2.3      
#> [49] vegan_2.4-5          tools_3.3.3          multcomp_1.4-8       munsell_0.5.0       
#> [53] cluster_2.0.5        bindrcpp_0.2         ade4_1.7-10          e1071_1.6-8         
#> [57] multcompView_0.1-7   caTools_1.17.1       rlang_0.1.6          grid_3.3.3          
#> [61] nloptr_1.0.4         miscTools_0.6-22     hermite_1.1.1        bitops_1.0-6        
#> [65] boot_1.3-18          DescTools_0.99.23    gtable_0.3.0         codetools_0.2-15    
#> [69] reshape_0.8.7        R6_2.4.1             gridExtra_2.3        zoo_1.8-1           
#> [73] mc2d_0.1-18          nortest_1.0-4        bindr_0.1.1          KernSmooth_2.23-15  
#> [77] permute_0.9-4        modeltools_0.2-22    parallel_3.3.3       Rcpp_0.12.14        
#> [81] lmtest_0.9-35        coda_0.19-1         



```
