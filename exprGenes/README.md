
## Expressed genes and data statistics

This code allows to reproduce the results of the ... 


## Contents

* [Getting Started](#getting-started)
  * [Required Packages](#required-packages)
  * [Data input](#data-input)
* [Data analysis](#data-analysis)
  * [Retrieve sample statistics](#retrieve-sample-statistics)
* [Visualization](#visualization)
* [Session info](#session-info)


## Getting started


### Required Packages
Install and load the following R packages before running the reproducible scripts:

```R
if (!require(plyr)) install.packages('plyr')
library(plyr)
if (!require(dplyr)) install.packages('dplyr')
library(dplyr)
if (!require(ggplot2)) install.packages('ggplot2')
library(ggplot2)
if (!require(mgcv)) install.packages('mgcv')
library(mgcv)
if (!require(grid)) install.packages('grid')
library(grid)

```
