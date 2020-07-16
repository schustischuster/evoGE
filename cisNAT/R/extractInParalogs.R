# Extract all in-paralog genes from OrthoFinder2 output
# Data input: Orthofinder2 tab-delimited file


#------------------- Load packages, set directories and read sample tables ---------------------


# Install and load packages
if (!require(dplyr)) install.packages('dplyr')
library(dplyr)


# Set file path and input files
in_dir <- "/Volumes/User/Shared/Christoph_manuscript/DevSeq_paper/Analysis/Analysis_2019/A_thaliana_gene_exression_map/20191121_CS_coding_cisNAT_analysis/data"
out_dir <- "/Volumes/User/Shared/Christoph_manuscript/DevSeq_paper/Analysis/Analysis_2019/A_thaliana_gene_exression_map/20191121_CS_coding_cisNAT_analysis"


extractInParalogs <- function(species = c("ATH", "AL", "CR", "ES", "TH", "MT", "BD")) {

    # Show error message if no species is chosen
    if (missing(species))
   
       stop(
       "Please choose one of the available species: 
	   'ATH', 'AL', 'CR', 'ES', 'TH', 'MT', 'BD'",
	   call. = TRUE
       )

    # Show error message if incorrect species is chosen
    if (!is.element(species, c("ATH", "AL", "CR", "ES", "TH", "MT", "BD")))
   
       stop(
       "Please choose one of the available species: 
	   'ATH', 'AL', 'CR', 'ES', 'TH', 'MT', 'BD'",
	   call. = TRUE
       )


    # Set ortholog table input file
    # It contains all orthogroups for AT/AL/CR/ES/TH/MT/BD from the Orthofinder2 output
    orthoTable = file.path(in_dir, "Orthofinder2", "DevSeq_Orthogroups.tsv")

	orthotable <- read.table(orthoTable, sep="\t", header=TRUE)




    #------- Process orthoginder data and remove inparalogs them from expression tables --------


    if (is.element("ATH", species)) {

    	orthotable_spec <- orthotable$Athaliana
        species_id <- "ATH"

    } else if (is.element("AL", species)) {

    	orthotable_spec <- orthotable$Alyrata
        species_id <- "AL"

    } else if (is.element("CR", species)) {

    	orthotable_spec <- orthotable$Crubella
        species_id <- "CR"

    } else if (is.element("ES", species)) {

    	orthotable_spec <- orthotable$Esalsugineum
        species_id <- "ES"

    } else if (is.element("TH", species)) {

    	orthotable_spec <- orthotable$Thassleriana
        species_id <- "TH"

    } else if (is.element("MT", species)) {

    	orthotable_spec <- orthotable$Mtruncatula
        species_id <- "MT"

    } else if (is.element("BD", species)) {

    	orthotable_spec <- orthotable$Bdistachyon
        species_id <- "BD"
    }


    # Get number of inparalog genes greater 0 in selected species 
    orth_num <- sapply((gregexpr(",", orthotable_spec, fixed=TRUE)), function(i) sum(i > 0))
    orth_spec_num <- data.frame(orthotable_spec, orth_num)
    colnames(orth_spec_num) <- c("orthogroup", "orth_greater_1")

    # Remove all orthogroups that only contain one member in target species
    # Note: it's part of an orthogroup because it contains in at least one of the species some inparalogs
    orth_spec_num_great1 <- orth_spec_num[orth_spec_num$orth_greater_1 > 0,]
    all_inparalog_transcripts <- orth_spec_num_great1$orthogroup

    # Simplify orthogroup transcript list to set of genes that have inparalogs in selected species
    all_inparalog_transcripts <-unlist(strsplit(as.character(all_inparalog_transcripts),", "))
    all_inparalog_genes <- gsub("\\..*","", all_inparalog_transcripts)
    all_inparalog_genes <- unique(all_inparalog_genes)

    all_inparalog_genes <- as.data.frame(all_inparalog_genes)
    colnames(all_inparalog_genes) <- paste0(species_id, "_in_paralogs", sep="")

    


    #----------------------------------- Write output to csv -----------------------------------


    fname <- sprintf('%s.csv', paste("All_inparalog_genes", species_id, sep="_"))

    write.table(all_inparalog_genes, file=file.path(out_dir, "output", "cd_gene_pairs", fname), 
		sep=";", dec=".", row.names=FALSE, col.names=TRUE)

}

extractInParalogs(species = "ATH")



