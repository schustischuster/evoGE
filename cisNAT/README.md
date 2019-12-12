## Coding-coding SAS and non-coding-coding SAS expression analysis

This code allows to reproduce the results of the protein-coding protein-coding sense-antisense (SAS) pair and non-coding protein-coding SAS pair expression analysis. 


## Contents

* [Getting Started](#getting-started)
  * [Required Packages](#required-packages)
  * [Data input](#data-input)
* [Data analysis](#data-analysis)
  * [Get coding-coding gene overlapp](#get-coding-coding-gene-overlapp)
  * [Get non-coding-coding gene overlapp](#get-non-coding-coding-gene-overlapp)
* [Visualization](#visualization)
* [Session info](#session-info)


## Getting started


### Required Packages
Install and load the following R packages before running the reproducible scripts:

```R
if (!require(dplyr)) install.packages('dplyr')
library(dplyr)
if (!require(dplyr)) install.packages('dplyr')
library(dplyr)
if (!require(GenomicRanges)) install.packages('GenomicRanges')
library(GenomicRanges)
if (!require(rtracklayer)) install.packages('rtracklayer')
library(rtracklayer)

```
  
### Data input
Download 'data' folder to the working directory on your computer. Then, set file path for input and output files: 

```R
in_dir <- "./data"
out_dir <- "."

```

## Data analysis

### Get coding-coding gene overlapp

The following function will extract all protein-coding protein-coding sense-antisense (SAS) pairs from the GTF file, apply an expression threshold, compute pairwise SAS correlations, and write the results to an CSV file. The threshold is set as follows: an expression value greater 0.5 TPM in at least two out of three replicates in at least one sample type. 

```R
getPcPc(species = "ATH", experiment = "single-species")

```
To generate all data tables used in this study, execute the following function calls: 

```R
getPcPc("ATH", "single-species")
getPcPc("ATH", "comparative")
getPcPc("AL", "single-species")
getPcPc("AL", "comparative")
getPcPc("CR")
getPcPc("ES")
getPcPc("TH")
getPcPc("MT")
getPcPc("BD")

```

### Get non-coding-coding gene overlapp

```R
getNcPc(species = "ATH", experiment = "single-species", threshold = 0.5)

```
