# Extract all in-paralog genes from OrthoFinder2 output
# Data input: Orthofinder2 tab-delimited file


#------------------- Load packages, set directories and read sample tables ---------------------


getInParalogs <- function(species = c("AT", "AL", "CR", "ES", "TH", "MT", "BD")) {

    # Show error message if no species is chosen
    if (missing(species))
   
       stop(
       "Please choose one of the available species: 
	   'AT', 'AL', 'CR', 'ES', 'TH', 'MT', 'BD'",
	   call. = TRUE
       )

    # Show error message if incorrect species is chosen
    if (!is.element(species, c("AT", "AL", "CR", "ES", "TH", "MT", "BD")))
   
       stop(
       "Please choose one of the available species: 
	   'AT', 'AL', 'CR', 'ES', 'TH', 'MT', 'BD'",
	   call. = TRUE
       )


    # Set ortholog table input file
    # It contains all orthogroups for AT/AL/CR/ES/TH/MT/BD from the Orthofinder2 output
    orthoTable = file.path(in_dir, "Orthofinder2", "DevSeq_Orthogroups.tsv")

	orthotable <- read.table(orthoTable, sep = "\t", header = TRUE)



    #------- Process orthoginder data and remove inparalogs them from expression tables --------


    if (is.element("AT", species)) {

    	orthotable_spec <- orthotable$Athaliana
        species_id <- "AT"

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
    orth_num <- sapply((gregexpr(",", orthotable_spec, fixed = TRUE)), function(i) sum(i > 0))
    orth_spec_num <- data.frame(orthotable_spec, orth_num)
    colnames(orth_spec_num) <- c("orthogroup", "num_orthologs")

    # Remove all orthogroups that only contain one member in target species
    # Note: it's part of an orthogroup because it contains in at least one of the species some inparalogs
    orth_spec_num_great1 <- orth_spec_num[orth_spec_num$num_orthologs > 0,]
    all_inparalog_transcripts <- orth_spec_num_great1$orthogroup
    all_inparalog_transcripts <-strsplit(as.character(all_inparalog_transcripts),", ")

    # A) Create data frame with species orthologs and same row length by adding NA values
    all_inparalog_transcripts_l <- lapply(all_inparalog_transcripts, `length<-`, max(lengths(all_inparalog_transcripts)))
    all_inparalog_transcripts_m <- do.call(rbind, all_inparalog_transcripts_l)
    orthogroup_genes <- as.data.frame(gsub("\\..*","", all_inparalog_transcripts_m))
    orthogroup_genes <- cbind(ID_orth = seq(1:nrow(orthogroup_genes)), orthogroup_genes)

    # B) Simplify orthogroup transcripts to a list inparalogs in selected species
    all_inparalog_transcripts <- unlist(all_inparalog_transcripts)
    all_inparalog_genes <- gsub("\\..*","", all_inparalog_transcripts)
    all_inparalog_genes <- as.data.frame(unique(all_inparalog_genes))
    colnames(all_inparalog_genes) <- paste0(species_id, "_inparalogs", sep = "")
    


    #----------------------------------- Write output to csv -----------------------------------

    
    fname_orth <- sprintf('%s.csv', paste("All_orthogroups", species_id, sep = "_"))
    
    fname_simpl <- sprintf('%s.csv', paste("All_inparalog_genes", species_id, sep = "_"))

    write.table(orthogroup_genes, file = file.path(out_dir, "output", "cd_gene_pairs", fname_orth), 
		sep = ";", dec = ".", row.names = FALSE, col.names = TRUE)

    write.table(all_inparalog_genes, file = file.path(out_dir, "output", "cd_gene_pairs", fname_simpl), 
		sep = ";", dec = ".", row.names = FALSE, col.names = TRUE)

}




