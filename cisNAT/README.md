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


## Getting started


### Required Packages
Install and load the following R packages before running the reproducible scripts:

```R

# Required packages
lib_List <- c("plyr", "dplyr", "GenomicRanges", "rtracklayer", "ggplot2", "scales", "mgcv", "data.table")

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
   for (nm in list.files(path, pattern = "[.][RrSsQq]$")) {
      if(trace) cat(nm,":")
      source(file.path(path, nm), ...)
      if(trace) cat("\n")
   }
}
 
sourceDir(path_to_R_files)

```

## Data analysis

### Retrieve coding-coding gene overlapp and pairwise expression correlation

The following function will extract all protein-coding protein-coding sense-antisense (SAS) pairs from the GTF file, apply an expression threshold, retrieve maximum and mean expression values for each gene, compute pairwise SAS correlations across all samples, and write the results to a CSV file. The threshold is set as follows: an expression value of both sense and antisense transcript greater than 0.5 TPM in at least two out of three replicates in at least one sample type. 

* `getPcPc(species = c("AT", "AL", "CR", "ES", "TH", "MT", "BD"), experiment = c("single-species", "comparative"), threshold)`

To generate all data tables used in this study, execute the following function calls: 

```R
species_ls <- list("AT", "AL", "CR", "ES", "TH", "MT", "BD")

getPcPc(species = "AT", experiment = "single-species", threshold = 0.5)
lapply(species_ls, getPcPc, experiment = "comparative", threshold = 0.5)

```

### Calculate pairwise non-coding protein-coding gene correlation

The following function will compute pairwise cis-natural antisense transcript (cisNAT)/protein-coding (PC) gene correlations across all samples, retrieve maximum and mean expression values for each gene, and write the results to a CSV file. A threshold of 0.5 TPM is applied for both cisNAT and PC genes. 

* `getCorNcPc(species = c("AT", "AL", "CR", "ES", "TH", "MT", "BD"), experiment = c("single-species", "comparative"))`

To generate all data tables used in this study, execute the following function calls: 

```R
species_ls <- list("AT", "AL", "CR", "ES", "TH", "MT", "BD")

getCorNcPc(species = "AT", experiment = "single-species")
lapply(species_ls, getCorNcPc, experiment = "comparative")

```

### Get cisNAT-protein-coding gene overlap length

The following function will extract the overlap length between cis-natural antisense transcripts (cisNATs) and protein-coding gene pairs from the species GTF files. The results will be written to a CSV file. 

```R
species_ls <- list("AT", "AL", "CR", "ES", "TH", "MT", "BD")

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

Set the file path for the data generated in the previous steps and source the R script:

```R
in_dir_cd <- file.path("evoGE-master", "cisNAT", "output", "overlap_pc_genes")
in_dir_nc <- file.path("evoGE-master", "cisNAT", "output", "overlap_nc_genes")
in_dir_PC_pairs <- file.path("evoGE-master", "cisNAT", "output", "overlap_nc_genes")
in_dir_NAT_cor <- file.path("evoGE-master", "cisNAT", "output", "NAT_expr_cor")

source(file.path("evoGE-master", "cisNAT", "R", "plotcisNAT.R"))

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
