## Make angiosperm phylogeny and orthologous gene plots

This code allows to generate the angiosperm phylogeny and coding and non-coding orthologous gene plots. 


## Contents

* [Getting Started](#getting-started)
  * [Required Packages](#required-packages)
  * [Data input](#data-input)
* [Data vizualization](#data-vizualization)
* [Session info](#session-info)

---
## Getting started


### Required Packages
Install and load the following R packages before running the reproducible scripts:

```R
# Create list of required packages
lib_List <- c("dplyr", "gplots", "ggplot2", "factoextra", "dendextend", "ggbeeswarm", "mblm", "lsmeans", "rcompanion")

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

To reproduce the results of this study, execute the following function calls:

```R
makeCompAnylsis(dataset="DevSeq", expr_estimation="TPM", coefficient="pearson", spec="Brassicaeae", data_norm="inter-organ")

```

This will generate the panels for the following figures:


---
## Session info

```R
sessionInfo()
```
