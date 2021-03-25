# Set working directory
setwd("/Volumes/User/Shared/Christoph_manuscript/DevSeq_paper/Analysis/Analysis_2019/A_thaliana_gene_exression_map")


# Set file path and input files
in_dir <- "./20200401_CS_exprGenes/data"
out_dir <- "./20200401_CS_exprGenes"


# Create list of required packages
lib_List <- c("dplyr", "data.table", "ggplot2", "mgcv", "grid", "gtable", "scales", "factoextra", 
"dendextend" )


loadLibrary <- function(x) { 
    if (!require(x, character.only = T)) {
        install.packages('x')
        library(x)
    }
}

# Load packages
# invisible() will suppress the contents of the list object 
invisible(lapply(lib_List, loadLibrary))
