
## Expressed genes and data statistics

This code allows to summarize the DevSeq data statistics and to reproduce the results of the intra-species transcriptome analyses. 


## Contents

* [Getting Started](#getting-started)
  * [Required Packages](#required-packages)
  * [Data input](#data-input)
* [Data analysis](#data-analysis)
  * [Retrieve mapping statistics](#retrieve-mapping-statistics)
  * [Retrieve number of expressed genes and transcripts](#retrieve-number-of-expressed-genes-and-transcripts)
* [Visualization](#visualization)
* [Session info](#session-info)

---
## Getting started


### Required Packages
Install and load the following R packages before running the reproducible scripts:

```R

# List of required packages
lib_List <- c("dplyr", "ggplot2", "data.table", "grid", "gtable", "scales", "factoextra", "dendextend")

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
in_dir <- file.path("evoGE", "exprGenes", "data")
out_dir <- file.path("evoGE", "exprGenes")
path_to_R_files <- file.path("evoGE", "exprGenes", "R")

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
## Data analysis

### Retrieve mapping statistics

The following function will merge the DevSeq mapping statistics and create a data table for each the _A.thaliana_ data set, the non-ATH data, and the comparative data set. 

```R
getStats()

```

### Retrieve number of expressed genes and transcripts

The following function will apply a threshold based on ERCC spike-ins at different threshold levels. ERCC spike-ins are a common set of external RNA controls that allow to measure  both sensitivity (lower limit of detection) and dynamic range of an RNA-Seq experiment. A gene (protein-coding/lncRNA) is considered to be expressed if it's expression value is above the threshold level in at least two out of three biological replicates. 

```R
getExprGenes(species = c("ATH", "AL", "CR", "ES", "TH", "MT", "BD"), 
             experiment = c("single-species", "comparative"), threshold)

```
</br>

| Arguments  |  |
| :---  | :---  |
| species  | A character string that defines the species to be analyzed. Can be one of `"ATH"` (A. thaliana), `"AL"` (A. lyrata), `"CR"` (C. rubella), `"ES"` (E. salsugineum), `"TH"` (T. hassleriana), `"MT"` (M. truncatula) or `"BD"` (B. distachyon). |
| experiment  | Defines which data set to use; `"single-species"` indicates all organs and stages, whereas `"comparative"` only includes organs of the comparative data set. |
| threshold  | A positive number between 0 and 1 defining the threshold. |

</br>

To reproduce the results of this study, execute the following function calls:

```R
thresholds <- list(0, 0.01, 0.05, 0.1)  # ERCC threshold values are 0 (static TPM threshold of 0.5)
                                        # or the 0.01/0.05/0.1 percentile of detected spike-ins

lapply(thresholds, getExprGenes, species = "ATH", experiment = "single-species")
lapply(thresholds, getExprGenes, species = "AL", experiment = "comparative")
lapply(thresholds, getExprGenes, species = "CR", experiment = "comparative")
lapply(thresholds, getExprGenes, species = "ES", experiment = "comparative")
lapply(thresholds, getExprGenes, species = "TH", experiment = "comparative")
lapply(thresholds, getExprGenes, species = "MT", experiment = "comparative")
lapply(thresholds, getExprGenes, species = "BD", experiment = "comparative")

```

To retrieve the number of expressed transcripts, the following function can be used:

```R
getExprTranscripts(species = c("ATH", "AL", "CR", "ES", "TH", "MT", "BD"), 
                   experiment = c("single-species", "comparative"), threshold)

```
It takes the same arguments described above. To reproduce the results of this study, execute the following function calls:

```R
thresholds <- list(0, 0.01, 0.05, 0.1)

lapply(thresholds, getExprTranscripts, species = "ATH", experiment = "single-species")
lapply(thresholds, getExprTranscripts, species = "AL", experiment = "comparative")
lapply(thresholds, getExprTranscripts, species = "CR", experiment = "comparative")
lapply(thresholds, getExprTranscripts, species = "ES", experiment = "comparative")
lapply(thresholds, getExprTranscripts, species = "TH", experiment = "comparative")
lapply(thresholds, getExprTranscripts, species = "MT", experiment = "comparative")
lapply(thresholds, getExprTranscripts, species = "BD", experiment = "comparative")

```


---
## Visualization

To visualize the results from the previous steps, execute the following function:

```R

plotExpr()

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
#> [1] grid      stats     graphics  grDevices utils     datasets  methods   base     

#> other attached packages:
#> [1] dendextend_1.12.0   factoextra_1.0.5    scales_0.5.0        gtable_0.3.0       
#> [5] mgcv_1.8-17         nlme_3.1-131        data.table_1.10.4-3 ggplot2_2.2.1      
#> [9] dplyr_0.7.4        

#> loaded via a namespace (and not attached):
#>  [1] Rcpp_0.12.14      bindr_0.1.1       magrittr_1.5      munsell_0.5.0    
#>  [5] colorspace_1.3-2  viridisLite_0.3.0 lattice_0.20-34   R6_2.4.1         
#>  [9] rlang_0.1.6       plyr_1.8.4        lazyeval_0.2.1    assertthat_0.2.1 
#> [13] tibble_1.3.4      Matrix_1.2-8      bindrcpp_0.2      gridExtra_2.3    
#> [17] viridis_0.5.1     ggrepel_0.7.0     glue_1.2.0        pkgconfig_2.0.3 

```
