# This script loads and processes gene expression correlation tables of sense-antisense (SAS) 
# gene pairs for coding/cisNATs SAS of the DevSeq data set
# It computes maximum expression level and average expression level for both coding and 
# non-coding SAS, and calculates the ratio between coding and non-coding gene expression
# Input expression values are in log2 TPM

#------------------- Load packages, set directories and read sample tables ---------------------



# Install and load packages
if (!require(dplyr)) install.packages('dplyr')
library(dplyr)


# Set file path and input files
in_dir <- "/Volumes/User/Shared/Christoph_manuscript/DevSeq_paper/Analysis/Analysis_2019/A_thaliana_gene_exression_map/20191121_CS_coding_cisNAT_analysis/output/overlap_nc_genes"
out_dir <- "/Volumes/User/Shared/Christoph_manuscript/DevSeq_paper/Analysis/Analysis_2019/A_thaliana_gene_exression_map/20191121_CS_coding_cisNAT_analysis"



getExprRatio <- function() {

	# Read all csv files in input file path
	readTableQuery <- function(path, pattern = "*query_expr_0.5.csv") {
		files = list.files(path, pattern, full.names = TRUE)
		lapply(files, function(x) read.table(x, sep=";", dec=".", header=TRUE, stringsAsFactors=FALSE))
	}
	readTableSubject <- function(path, pattern = "*subject_expr_0.5.csv") {
		files = list.files(path, pattern, full.names = TRUE)
		lapply(files, function(x) read.table(x, sep=";", dec=".", header=TRUE, stringsAsFactors=FALSE))
	}

	query_gene_tables <- readTableQuery(in_dir)
	subject_gene_tables <- readTableSubject(in_dir)


	# Get file names and save them in character vector
	query_gene_tables_list <- as.character(list.files(in_dir, pattern = "*query_expr_0.5.csv"))
	query_gene_tables_names <- gsub('\\.csv$', '', query_gene_tables_list)
	subject_gene_tables_list <- as.character(list.files(in_dir, pattern = "*subject_expr_0.5.csv"))
	subject_gene_tables_names <- gsub('\\.csv$', '', subject_gene_tables_list)


	# Change data frame names in list
	names(query_gene_tables) <- query_gene_tables_names
	names(subject_gene_tables) <- subject_gene_tables_names

	# list2env(query_gene_tables, envir = .GlobalEnv)
	# list2env(subject_gene_tables, envir = .GlobalEnv)

	# Stop function here to allow specific analysis of a single species
    return_list <- list("query_gene_tables" = query_gene_tables, "subject_gene_tables" = subject_gene_tables)
    return(return_list)
	}
	genes_tables <- getExprRatio()
	list2env(genes_tables, envir = .GlobalEnv)




	#------- Generate NAT and coding tables containing average and maximum expression levels --------


	# Filter tables by biotype (NAT | conding), select and rename columns, and extend table names
	NAT_query_tables <- lapply(query_gene_tables, function(x) {
		x <- x  %>% filter(biotype_query == "lnc_exonic_antisense" | biotype_query == "lnc_intronic_antisense")
		x <- dplyr::select(x, -c(id_minus_strand, id_minus, start_minus, end_minus, width_subject, 
			strand_subject, biotype_subject, gene_source_subject))
		x <- dplyr::rename(x, gene_id=id_plus_strand, id=id_plus, start=start_plus, end=end_plus, 
			width=width_query, strand=strand_query, biotype=biotype_query, gene_source=gene_source_query)
	})
	names(NAT_query_tables) <- lapply(names(NAT_query_tables), function(x) {x <- paste0(x, "_NAT")
	})

	NAT_subject_tables <- lapply(subject_gene_tables, function(x) {
		x <- x  %>% filter(biotype_subject == "lnc_exonic_antisense" | biotype_subject == "lnc_intronic_antisense")
		x <- dplyr::select(x, -c(id_plus_strand, id_plus, start_plus, end_plus, width_query, 
			strand_query, biotype_query, gene_source_query))
		x <- dplyr::rename(x, gene_id=id_minus_strand, id=id_minus, start=start_minus, end=end_minus, 
			width=width_subject, strand=strand_subject, biotype=biotype_subject, gene_source=gene_source_subject)
	})
	names(NAT_subject_tables) <- lapply(names(NAT_subject_tables), function(x) {x <- paste0(x, "_NAT")
	})

	list2env(NAT_query_tables, envir = .GlobalEnv)
	list2env(NAT_subject_tables, envir = .GlobalEnv)


	# Create combined "+" and "-" strand NAT lists for each species
	AL_comparative_samples_cd_nc_SAS_cor_wo_pollen_expr_0.5_NAT <- rbind(AL_comparative_samples_cd_nc_SAS_cor_wo_pollen_query_expr_0.5_NAT,
		AL_comparative_samples_cd_nc_SAS_cor_wo_pollen_subject_expr_0.5_NAT)
	ATH_cd_nc_SAS_cor_wo_pollen_expr_0.5_NAT <- rbind(ATH_cd_nc_SAS_cor_wo_pollen_query_expr_0.5_NAT,
		ATH_cd_nc_SAS_cor_wo_pollen_subject_expr_0.5_NAT)
	ATH_comparative_samples_cd_nc_SAS_cor_wo_pollen_expr_0.5_NAT <- rbind(ATH_comparative_samples_cd_nc_SAS_cor_wo_pollen_query_expr_0.5_NAT,
		ATH_comparative_samples_cd_nc_SAS_cor_wo_pollen_subject_expr_0.5_NAT)
	BD_cd_nc_SAS_cor_wo_pollen_expr_0.5_NAT <- rbind(BD_cd_nc_SAS_cor_wo_pollen_query_expr_0.5_NAT,
		BD_cd_nc_SAS_cor_wo_pollen_subject_expr_0.5_NAT)
	CR_cd_nc_SAS_cor_wo_pollen_expr_0.5_NAT <- rbind(CR_cd_nc_SAS_cor_wo_pollen_query_expr_0.5_NAT,
		CR_cd_nc_SAS_cor_wo_pollen_subject_expr_0.5_NAT)
	ES_cd_nc_SAS_cor_wo_pollen_expr_0.5_NAT <- rbind(ES_cd_nc_SAS_cor_wo_pollen_query_expr_0.5_NAT,
		ES_cd_nc_SAS_cor_wo_pollen_subject_expr_0.5_NAT)
	MT_cd_nc_SAS_cor_wo_pollen_expr_0.5_NAT <- rbind(MT_cd_nc_SAS_cor_wo_pollen_query_expr_0.5_NAT,
		MT_cd_nc_SAS_cor_wo_pollen_subject_expr_0.5_NAT)
	TH_cd_nc_SAS_cor_wo_pollen_expr_0.5_NAT <- rbind(TH_cd_nc_SAS_cor_wo_pollen_query_expr_0.5_NAT,
		TH_cd_nc_SAS_cor_wo_pollen_subject_expr_0.5_NAT)



	#-------------- Find 'ATGE_NAT_ID' and 'ATH_cd_nc_SAS_cor_wo_pollen_0.5' overlap ----------------


	# Work only with ATGE protein-coding gene IDs since lncRNA IDs may have changed

	# Extract protein-coding genes from DevSeq cd-nc SAS ATH data
	DevSeq_NAT_SAS_pair_plus_coding <- ATH_cd_nc_SAS_cor_wo_pollen_0.5 %>% filter(
		biotype_query == "protein_coding")
	DevSeq_NAT_plus_coding_ID <- DevSeq_NAT_SAS_pair_plus_coding %>% select(id_plus_strand)
	names(DevSeq_NAT_plus_coding_ID) <- "coding_id"

	DevSeq_NAT_SAS_pair_minus_coding <- ATH_cd_nc_SAS_cor_wo_pollen_0.5 %>% filter(
		biotype_subject == "protein_coding")
	DevSeq_NAT_minus_coding_ID <- DevSeq_NAT_SAS_pair_minus_coding %>% select(id_minus_strand)
	names(DevSeq_NAT_minus_coding_ID) <- "coding_id"

	# Combine plus and minus strand protein ID lists
	DevSeq_NAT_coding_ID <- rbind(DevSeq_NAT_plus_coding_ID, DevSeq_NAT_minus_coding_ID)

	# Select all protein-coding genes from DevSeq cd-nc SAS ATH data that are in ATGE SAS list 
	DevSeq_ATGE_NAT_coding_ID <- DevSeq_NAT_coding_ID[(
		((DevSeq_NAT_coding_ID$coding_id %in% ATGE_NAT_ID$AGI_1) | 
			(DevSeq_NAT_coding_ID$coding_id %in% ATGE_NAT_ID$AGI_2))
		),]

	# Select rows in "ATH_cd_nc_SAS_cor_wo_pollen_0.5" with a protein id in "DevSeq_ATGE_NAT_coding_ID"
	ATH_cd_nc_SAS_cor_wo_pollen_0.5_in_ATGE <- ATH_cd_nc_SAS_cor_wo_pollen_0.5[(
		((ATH_cd_nc_SAS_cor_wo_pollen_0.5$id_plus_strand %in% DevSeq_ATGE_NAT_coding_ID) | 
			(ATH_cd_nc_SAS_cor_wo_pollen_0.5$id_minus_strand %in% DevSeq_ATGE_NAT_coding_ID))
		),]



	#---------------------------------------- Write output -----------------------------------------


	# Write final data tables to csv files and store them in /out_dir/output/data_tables
	if (!dir.exists(file.path(out_dir, "output", "SAS_DevSeq_ATGE"))) 
		dir.create(file.path(out_dir, "output", "SAS_DevSeq_ATGE"), recursive = TRUE)

	write.table(ATH_cd_nc_SAS_cor_wo_pollen_0.5_in_ATGE, 
		file=file.path(out_dir, "output", "SAS_DevSeq_ATGE", "ATH_cd_nc_SAS_cor_wo_pollen_0.5_in_ATGE.csv"), 
		sep=";", dec=".", row.names=FALSE, col.names=TRUE)

}


