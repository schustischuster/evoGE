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
Download 'data' folder to the working directory on your computer. Then, set file path and input files: 

```R
in_dir <- "./data"
out_dir <- "."

```

## Data analysis

### Get coding-coding gene overlapp

```R
# Run pc-pc SAS analysis
```
### Get non-coding-coding gene overlapp
