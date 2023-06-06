## Coding-coding SAS and non-coding-coding SAS expression analysis

This code allows to reproduce the results of the protein-coding protein-coding sense-antisense (SAS) pair and non-coding protein-coding SAS pair expression analysis. 


## Contents

* [Getting Started](#getting-started)
  * [Required Packages](#required-packages)
  * [Data input](#data-input)
* [Data analysis](#data-analysis)
  * [Retrieve coding-coding gene overlapp and pairwise expression correlation](#retrieve-coding-coding-gene-overlapp-and-pairwise-expression-correlation)
  * [Calculate pairwise non-coding/protein-coding gene correlation](#calculate-pairwise-non-coding-protein-coding-gene-correlation)
  * [Get cisNAT-protein-coding gene overlap length](#get-cisNAT-protein-coding-gene-overlap-length)
  * [Get intergenic distance of neighboring genes](#get-intergenic-distance-of-neighboring-genes)
* [Visualization](#visualization)
* [Session info](#session-info)

---
## Getting started


### Required Packages

Install Bioconductor core packages, GenomicRanges and rtracklayer:

```R
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install()
BiocManager::install("GenomicRanges")
BiocManager::install("rtracklayer")

```

Install and load the following R packages before running the reproducible scripts:

```R
# Required packages
lib_List <- c("plyr", "dplyr", "GenomicRanges", "rtracklayer", "ggplot2", "scales", "mgcv", "data.table", "R.utils")

# Install missing packages
instpack <- lib_List %in% installed.packages()[,"Package"]
if (any(instpack == FALSE)) {
  install.packages(lib_List[!instpack])
}

# Load packages
invisible(lapply(lib_List, library, character.only = TRUE))

```
  
### Data input
Download and extract the entire directory to the working directory on your computer. Then, set the path for input and output files and source the R scripts: 

```R
in_dir <- file.path("evoGE-master", "cisNAT", "data")
out_dir <- file.path("evoGE-master", "cisNAT")
path_to_R_files <- file.path("evoGE-master", "cisNAT", "R")

# Source R files
sourceDir <- function(path, trace = TRUE, ...) {
   for (nm in list.files(path, "^[^plotcisNAT].+[.][RrSsQq]$")) {
      if(trace) cat(nm,":")
      source(file.path(path, nm), ...)
      if(trace) cat("\n")
   }
}
 
sourceDir(path_to_R_files)

```

Unzip the GTF files:

```R
gtf_file_ls <- list.files(file.path(in_dir, "GTF"), pattern = "[.][gz]")

unzipFiles <- function(f) {
  file <- paste(file.path(in_dir, "GTF"), f, sep = "/")
  gunzip(file, remove = FALSE)
}

# Unzip GTF files
lapply(gtf_file_ls, unzipFiles)

```
---
## Data analysis

Data analysis is performed across all species defined in a species list:

```R
species_ls <- list("AT", "AL", "CR", "ES", "TH", "MT", "BD")

```

### Retrieve coding-coding gene overlapp and pairwise expression correlation

The following function will extract all protein-coding protein-coding sense-antisense (SAS) pairs from the GTF file, apply an expression threshold, retrieve maximum and mean expression values for each gene, compute pairwise SAS correlations across all samples, and write the results to a CSV file. The threshold is set as follows: an expression value of both sense and antisense transcript greater than 0.5 TPM in at least two out of three replicates in at least one sample type. 

* `getPcPc(species = c("AT", "AL", "CR", "ES", "TH", "MT", "BD"), experiment = c("single-species", "comparative"), threshold)`


| Arguments  |  |
| :---  | :---  |
| species  | Defines the species to be analyzed. Can be one of "AT" (*Arabidopsis thaliana*), "AL" (*Arabidopsis lyrata*), "CR" (*Capsella rubella*), "ES" (*Eutrema salsugineum*), "TH" (*Tarenaya hassleriana*), "MT" (*Medicago truncatula*), "BD" (*Brachypodium distachyon*).|
| experiment  | Type of experiment. Choose "single-species" to run the analysis on all samples, and "comparative" to limit the analysis to the comparative samples. |
| threshold  | Indicates the TPM threshold above which a gene is considered expressed. |


To generate all data tables used in this study, execute the following function calls: 

```R
getPcPc(species = "AT", experiment = "single-species", threshold = 0.5)
lapply(species_ls, getPcPc, experiment = "comparative", threshold = 0.5)

```

### Calculate pairwise non-coding protein-coding gene correlation

The following function will compute pairwise cis-natural antisense transcript (cisNAT)/protein-coding (PC) gene correlations across all samples, retrieve maximum and mean expression values for each gene, and write the results to a CSV file. A threshold of 0.5 TPM is applied for both cisNAT and PC genes. 

* `getCorNcPc(species = c("AT", "AL", "CR", "ES", "TH", "MT", "BD"), experiment = c("single-species", "comparative"))`

To generate all data tables used in this study, execute the following function calls: 

```R
getCorNcPc(species = "AT", experiment = "single-species")
lapply(species_ls, getCorNcPc, experiment = "comparative")

```

### Get cisNAT-protein-coding gene overlap length

The following function will extract the overlap length between cis-natural antisense transcripts (cisNATs) and protein-coding gene pairs from the species GTF files. The results will be written to a CSV file. 

```R
lapply(species_ls, getNcPcOverlap)

```

### Get intergenic distance of neighboring genes

Numerous studies using different model organisms including _Arabidopsis_, _Caenorhabditis_, _Drosophila_, _human_, and _Saccharomyces_ have shown that neighboring genes tend to be coexpressed, e.g. [Cohen  et al. (2000)](https://www.ncbi.nlm.nih.gov/pubmed/11017073), [Boutanaev et al. (2000)](https://www.ncbi.nlm.nih.gov/pubmed/12478293), [Lercher et al. (2002)](https://www.ncbi.nlm.nih.gov/pubmed/11992122), [Lercher et al. (2003)](https://www.ncbi.nlm.nih.gov/pubmed/12566401), [Williams and Bowles (2004)](https://www.ncbi.nlm.nih.gov/pubmed/15173112). We wanted to test if a similar trend can be found in the DevSeq data set. The following function will extract all neighbouring protein-coding gene pairs and their intergenic distance from the GTF file, apply an expression threshold, compute pairwise Pearson and Spearman expression correlations across all samples, and write the results to a CSV file.

* `getPcPcNO(species = c("AT", "AL", "CR", "ES", "TH", "MT", "BD"), experiment = c("single-species", "comparative"), threshold)`

To generate the data table for _Arabidopsis thaliana_, execute the following function call: 

```R
getPcPcNO("AT", "single-species", 0.5)

```

## Visualization

To visualize the results generated in the previous steps, source the following R script:

```R
source(file.path("evoGE-master", "cisNAT", "R", "plotcisNAT.R"))

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
#> [1] parallel  stats4    stats     graphics  grDevices utils     datasets  methods   base    

#> other attached packages:
#> [1] mgcv_1.8-17          nlme_3.1-131         ggplot2_2.2.1        rtracklayer_1.34.2   GenomicRanges_1.26.4
#> [6] GenomeInfoDb_1.10.3  IRanges_2.8.2        S4Vectors_0.12.2     BiocGenerics_0.20.0  dplyr_0.7.4 
#>[11] plyr_1.8.4

#> loaded via a namespace (and not attached):
#> [1] Rcpp_0.12.14               bindr_0.1.1                XVector_0.14.1             magrittr_1.5              
#> [5] zlibbioc_1.20.0            GenomicAlignments_1.10.1   munsell_0.5.0              BiocParallel_1.8.2        
#> [9] colorspace_1.3-2           lattice_0.20-34            R6_2.4.1                   rlang_0.1.6               
#>[13] tools_3.3.3                grid_3.3.3                 SummarizedExperiment_1.4.0 gtable_0.3.0              
#>[17] Biobase_2.34.0             lazyeval_0.2.1             assertthat_0.2.1           tibble_1.3.4              
#>[21] Matrix_1.2-8               bindrcpp_0.2               bitops_1.0-6               RCurl_1.95-4.10           
#>[25] glue_1.2.0                 scales_0.5.0               Biostrings_2.42.1          Rsamtools_1.26.2          
#>[29] XML_3.98-1.9               pkgconfig_2.0.3  

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
#> [1] stats4    stats     graphics  grDevices utils     datasets  methods   base     

#> other attached packages:
#> [1] R.utils_2.12.2       R.oo_1.25.0          R.methodsS3_1.8.2    data.table_1.14.2    mgcv_1.8-39         
#> [6] nlme_3.1-155         scales_1.1.1         ggplot2_3.3.5        rtracklayer_1.54.0   GenomicRanges_1.46.1
#>[11] GenomeInfoDb_1.30.1  IRanges_2.28.0       S4Vectors_0.32.4     BiocGenerics_0.40.0  dplyr_1.0.8         
#>[16] plyr_1.8.6           BiocManager_1.30.16 

#> loaded via a namespace (and not attached):
#> [1] SummarizedExperiment_1.24.0 tidyselect_1.1.2            purrr_0.3.4                 splines_4.1.3              
#> [5] lattice_0.20-45             colorspace_2.0-3            vctrs_0.3.8                 generics_0.1.2             
#> [9] yaml_2.3.5                  utf8_1.2.2                  XML_3.99-0.14               rlang_1.0.2                
#>[13] pillar_1.7.0                glue_1.6.2                  withr_2.5.0                 BiocParallel_1.28.3        
#>[17] matrixStats_0.61.0          GenomeInfoDbData_1.2.7      lifecycle_1.0.1             zlibbioc_1.40.0            
#>[21] MatrixGenerics_1.6.0        Biostrings_2.62.0           munsell_0.5.0               gtable_0.3.0               
#>[25] restfulr_0.0.15             labeling_0.4.2              Biobase_2.54.0              parallel_4.1.3             
#>[29] fansi_1.0.2                 Rcpp_1.0.8.2                KernSmooth_2.23-20          DelayedArray_0.20.0        
#>[33] XVector_0.34.0              farver_2.1.0                Rsamtools_2.10.0            digest_0.6.29              
#>[37] rjson_0.2.21                BiocIO_1.4.0                grid_4.1.3                  cli_3.2.0                  
#>[41] tools_4.1.3                 bitops_1.0-7                magrittr_2.0.2              RCurl_1.98-1.12            
#>[45] tibble_3.1.6                crayon_1.5.0                pkgconfig_2.0.3             ellipsis_0.3.2             
#>[49] Matrix_1.4-0                R6_2.5.1                    GenomicAlignments_1.30.0    compiler_4.1.3

```
