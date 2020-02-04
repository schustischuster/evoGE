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


	# Filter tables by biotype (NAT), select and rename columns, and extend table names
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


	# Filter tables by biotype (protein_conding), select and rename columns, and extend table names
	cd_query_tables <- lapply(query_gene_tables, function(x) {
		x <- x  %>% filter(biotype_query == "protein_coding")
		x <- dplyr::select(x, -c(id_minus_strand, id_minus, start_minus, end_minus, width_subject, 
			strand_subject, biotype_subject, gene_source_subject))
		x <- dplyr::rename(x, gene_id=id_plus_strand, id=id_plus, start=start_plus, end=end_plus, 
			width=width_query, strand=strand_query, biotype=biotype_query, gene_source=gene_source_query)
	})
	names(cd_query_tables) <- lapply(names(cd_query_tables), function(x) {x <- paste0(x, "_cd")
	})

	cd_subject_tables <- lapply(subject_gene_tables, function(x) {
		x <- x  %>% filter(biotype_subject == "protein_coding")
		x <- dplyr::select(x, -c(id_plus_strand, id_plus, start_plus, end_plus, width_query, 
			strand_query, biotype_query, gene_source_query))
		x <- dplyr::rename(x, gene_id=id_minus_strand, id=id_minus, start=start_minus, end=end_minus, 
			width=width_subject, strand=strand_subject, biotype=biotype_subject, gene_source=gene_source_subject)
	})
	names(cd_subject_tables) <- lapply(names(cd_subject_tables), function(x) {x <- paste0(x, "_cd")
	})

	list2env(cd_query_tables, envir = .GlobalEnv)
	list2env(cd_subject_tables, envir = .GlobalEnv)


	# Create combined "+" and "-" strand NAT lists for each species
	AL_comparative_samples_cd_nc_SAS_cor_wo_pollen_expr_0.5_cd <- rbind(AL_comparative_samples_cd_nc_SAS_cor_wo_pollen_query_expr_0.5_cd,
		AL_comparative_samples_cd_nc_SAS_cor_wo_pollen_subject_expr_0.5_cd)
	ATH_cd_nc_SAS_cor_wo_pollen_expr_0.5_cd <- rbind(ATH_cd_nc_SAS_cor_wo_pollen_query_expr_0.5_cd,
		ATH_cd_nc_SAS_cor_wo_pollen_subject_expr_0.5_cd)
	ATH_comparative_samples_cd_nc_SAS_cor_wo_pollen_expr_0.5_cd <- rbind(ATH_comparative_samples_cd_nc_SAS_cor_wo_pollen_query_expr_0.5_cd,
		ATH_comparative_samples_cd_nc_SAS_cor_wo_pollen_subject_expr_0.5_cd)
	BD_cd_nc_SAS_cor_wo_pollen_expr_0.5_cd <- rbind(BD_cd_nc_SAS_cor_wo_pollen_query_expr_0.5_cd,
		BD_cd_nc_SAS_cor_wo_pollen_subject_expr_0.5_cd)
	CR_cd_nc_SAS_cor_wo_pollen_expr_0.5_cd <- rbind(CR_cd_nc_SAS_cor_wo_pollen_query_expr_0.5_cd,
		CR_cd_nc_SAS_cor_wo_pollen_subject_expr_0.5_cd)
	ES_cd_nc_SAS_cor_wo_pollen_expr_0.5_cd <- rbind(ES_cd_nc_SAS_cor_wo_pollen_query_expr_0.5_cd,
		ES_cd_nc_SAS_cor_wo_pollen_subject_expr_0.5_cd)
	MT_cd_nc_SAS_cor_wo_pollen_expr_0.5_cd <- rbind(MT_cd_nc_SAS_cor_wo_pollen_query_expr_0.5_cd,
		MT_cd_nc_SAS_cor_wo_pollen_subject_expr_0.5_cd)
	TH_cd_nc_SAS_cor_wo_pollen_expr_0.5_cd <- rbind(TH_cd_nc_SAS_cor_wo_pollen_query_expr_0.5_cd,
		TH_cd_nc_SAS_cor_wo_pollen_subject_expr_0.5_cd)


	# Make sure there are only matching id's in NAT and coding tables
	AL_comparative_samples_cd_nc_SAS_cor_wo_pollen_expr_0.5_NAT <- AL_comparative_samples_cd_nc_SAS_cor_wo_pollen_expr_0.5_NAT[(
		AL_comparative_samples_cd_nc_SAS_cor_wo_pollen_expr_0.5_NAT$id %in% AL_comparative_samples_cd_nc_SAS_cor_wo_pollen_expr_0.5_cd$id),]
	ATH_cd_nc_SAS_cor_wo_pollen_expr_0.5_NAT <- ATH_cd_nc_SAS_cor_wo_pollen_expr_0.5_NAT[(
		ATH_cd_nc_SAS_cor_wo_pollen_expr_0.5_NAT$id %in% ATH_cd_nc_SAS_cor_wo_pollen_expr_0.5_cd$id),]
	ATH_comparative_samples_cd_nc_SAS_cor_wo_pollen_expr_0.5_NAT <- ATH_comparative_samples_cd_nc_SAS_cor_wo_pollen_expr_0.5_NAT[(
		ATH_comparative_samples_cd_nc_SAS_cor_wo_pollen_expr_0.5_NAT$id %in% ATH_comparative_samples_cd_nc_SAS_cor_wo_pollen_expr_0.5_cd$id),]
	BD_cd_nc_SAS_cor_wo_pollen_expr_0.5_NAT <- BD_cd_nc_SAS_cor_wo_pollen_expr_0.5_NAT[(
		BD_cd_nc_SAS_cor_wo_pollen_expr_0.5_NAT$id %in% BD_cd_nc_SAS_cor_wo_pollen_expr_0.5_cd$id),]
	CR_cd_nc_SAS_cor_wo_pollen_expr_0.5_NAT <- CR_cd_nc_SAS_cor_wo_pollen_expr_0.5_NAT[(
		CR_cd_nc_SAS_cor_wo_pollen_expr_0.5_NAT$id %in% CR_cd_nc_SAS_cor_wo_pollen_expr_0.5_cd$id),]
	ES_cd_nc_SAS_cor_wo_pollen_expr_0.5_NAT <- ES_cd_nc_SAS_cor_wo_pollen_expr_0.5_NAT[(
		ES_cd_nc_SAS_cor_wo_pollen_expr_0.5_NAT$id %in% ES_cd_nc_SAS_cor_wo_pollen_expr_0.5_cd$id),]
	MT_cd_nc_SAS_cor_wo_pollen_expr_0.5_NAT <- MT_cd_nc_SAS_cor_wo_pollen_expr_0.5_NAT[(
		MT_cd_nc_SAS_cor_wo_pollen_expr_0.5_NAT$id %in% MT_cd_nc_SAS_cor_wo_pollen_expr_0.5_cd$id),]
	TH_cd_nc_SAS_cor_wo_pollen_expr_0.5_NAT <- TH_cd_nc_SAS_cor_wo_pollen_expr_0.5_NAT[(
		TH_cd_nc_SAS_cor_wo_pollen_expr_0.5_NAT$id %in% TH_cd_nc_SAS_cor_wo_pollen_expr_0.5_cd$id),]


	AL_comparative_samples_cd_nc_SAS_cor_wo_pollen_expr_0.5_cd <- AL_comparative_samples_cd_nc_SAS_cor_wo_pollen_expr_0.5_cd[(
		AL_comparative_samples_cd_nc_SAS_cor_wo_pollen_expr_0.5_cd$id %in% AL_comparative_samples_cd_nc_SAS_cor_wo_pollen_expr_0.5_NAT$id),]
	ATH_cd_nc_SAS_cor_wo_pollen_expr_0.5_cd <- ATH_cd_nc_SAS_cor_wo_pollen_expr_0.5_cd[(
		ATH_cd_nc_SAS_cor_wo_pollen_expr_0.5_cd$id %in% ATH_cd_nc_SAS_cor_wo_pollen_expr_0.5_NAT$id),]
	ATH_comparative_samples_cd_nc_SAS_cor_wo_pollen_expr_0.5_cd <- ATH_comparative_samples_cd_nc_SAS_cor_wo_pollen_expr_0.5_cd[(
		ATH_comparative_samples_cd_nc_SAS_cor_wo_pollen_expr_0.5_cd$id %in% ATH_comparative_samples_cd_nc_SAS_cor_wo_pollen_expr_0.5_NAT$id),]
	BD_cd_nc_SAS_cor_wo_pollen_expr_0.5_cd <- BD_cd_nc_SAS_cor_wo_pollen_expr_0.5_cd[(
		BD_cd_nc_SAS_cor_wo_pollen_expr_0.5_cd$id %in% BD_cd_nc_SAS_cor_wo_pollen_expr_0.5_NAT$id),]
	CR_cd_nc_SAS_cor_wo_pollen_expr_0.5_cd <- CR_cd_nc_SAS_cor_wo_pollen_expr_0.5_cd[(
		CR_cd_nc_SAS_cor_wo_pollen_expr_0.5_cd$id %in% CR_cd_nc_SAS_cor_wo_pollen_expr_0.5_NAT$id),]
	ES_cd_nc_SAS_cor_wo_pollen_expr_0.5_cd <- ES_cd_nc_SAS_cor_wo_pollen_expr_0.5_cd[(
		ES_cd_nc_SAS_cor_wo_pollen_expr_0.5_cd$id %in% ES_cd_nc_SAS_cor_wo_pollen_expr_0.5_NAT$id),]
	MT_cd_nc_SAS_cor_wo_pollen_expr_0.5_cd <- MT_cd_nc_SAS_cor_wo_pollen_expr_0.5_cd[(
		MT_cd_nc_SAS_cor_wo_pollen_expr_0.5_cd$id %in% MT_cd_nc_SAS_cor_wo_pollen_expr_0.5_NAT$id),]
	TH_cd_nc_SAS_cor_wo_pollen_expr_0.5_cd <- TH_cd_nc_SAS_cor_wo_pollen_expr_0.5_cd[(
		TH_cd_nc_SAS_cor_wo_pollen_expr_0.5_cd$id %in% TH_cd_nc_SAS_cor_wo_pollen_expr_0.5_NAT$id),]


	# Create NAT and coding lists
	NAT_list_all_spec <- list("AL_NAT"=AL_comparative_samples_cd_nc_SAS_cor_wo_pollen_expr_0.5_NAT,
		"ATH_all_NAT"=ATH_cd_nc_SAS_cor_wo_pollen_expr_0.5_NAT, "ATH_comp_NAT"=ATH_comparative_samples_cd_nc_SAS_cor_wo_pollen_expr_0.5_NAT,
		"BD_NAT"=BD_cd_nc_SAS_cor_wo_pollen_expr_0.5_NAT, "CR_NAT"=CR_cd_nc_SAS_cor_wo_pollen_expr_0.5_NAT,
		"ES_NAT"=ES_cd_nc_SAS_cor_wo_pollen_expr_0.5_NAT, "MT_NAT"=MT_cd_nc_SAS_cor_wo_pollen_expr_0.5_NAT,
		"TH_NAT"=TH_cd_nc_SAS_cor_wo_pollen_expr_0.5_NAT)

	cd_list_all_spec <- list("AL_cd"=AL_comparative_samples_cd_nc_SAS_cor_wo_pollen_expr_0.5_cd, 
		"ATH_all_cd"=ATH_cd_nc_SAS_cor_wo_pollen_expr_0.5_cd, "ATH_comp_cd"=ATH_comparative_samples_cd_nc_SAS_cor_wo_pollen_expr_0.5_cd, 
		"BD_cd"=BD_cd_nc_SAS_cor_wo_pollen_expr_0.5_cd, "CR_cd"=CR_cd_nc_SAS_cor_wo_pollen_expr_0.5_cd, 
		"ES_cd"=ES_cd_nc_SAS_cor_wo_pollen_expr_0.5_cd, "MT_cd"=MT_cd_nc_SAS_cor_wo_pollen_expr_0.5_cd, 
		"TH_cd"=TH_cd_nc_SAS_cor_wo_pollen_expr_0.5_cd)


	# Compute average expression for replicates
	getAverage <- function(x) {

		# Compute average expression for all elements in list
		calculateAvgExpr <- function(df) {

			# Split data frame by sample replicates into a list
			# then get rowMeans for each subset, simplify output and bind to gene_id column
	
			averaged_replicates <- data.frame(df[1:12],

				sapply(split.default(df[13:ncol(df)], 
					rep(seq_along(df), 
					each = 3, 
					length.out=ncol(df)-12)
					), rowMeans)
				)
		
			return(averaged_replicates)
		}
		y <- calculateAvgExpr(x)

		# define column name string
		colnames_descr <- names(x[1:12])
		colnames_samples <- unique(gsub(".{2}$", "", names(x[13:ncol(x)])))
		colnames(y) <- c(colnames_descr, colnames_samples)
		return(y)
	}

	NAT_list_all_spec_avg <- lapply(NAT_list_all_spec, getAverage)
	cd_list_all_spec_avg <- lapply(cd_list_all_spec, getAverage)


	# Get maximum expression for each gene of each element in list
	getMax <- function(x) {
		x$max <- do.call(pmax, x[13:ncol(x)])
		return(x)
	}

	NAT_list_all_spec_avg_max <- lapply(NAT_list_all_spec_avg, getMax)
	cd_list_all_spec_avg_max <- lapply(cd_list_all_spec_avg, getMax)


	# Compute mean expression for each gene of each element in list
	getMean <- function(df) {
		last_col <- (ncol(df))-1
		avg_df <- data.frame(df[1:ncol(df)], means=rowMeans(df[13:last_col]))
		return(avg_df)
	}

	NAT_list_all_spec_avg_max_avg <- lapply(NAT_list_all_spec_avg_max, getMean)
	cd_list_all_spec_avg_max_avg <- lapply(cd_list_all_spec_avg_max, getMean)


	list2env(NAT_list_all_spec_avg_max_avg, envir = .GlobalEnv)
	list2env(cd_list_all_spec_avg_max_avg, envir = .GlobalEnv)




	#------- Generate NAT tables containing relative (max)NAT-coding gene expression levels --------




















	#---------------------------------------- Write output -----------------------------------------


	# Write final data tables to csv files and store them in /out_dir/output/data_tables
	if (!dir.exists(file.path(out_dir, "output", "SAS_DevSeq_ATGE"))) 
		dir.create(file.path(out_dir, "output", "SAS_DevSeq_ATGE"), recursive = TRUE)

	write.table(ATH_cd_nc_SAS_cor_wo_pollen_0.5_in_ATGE, 
		file=file.path(out_dir, "output", "SAS_DevSeq_ATGE", "ATH_cd_nc_SAS_cor_wo_pollen_0.5_in_ATGE.csv"), 
		sep=";", dec=".", row.names=FALSE, col.names=TRUE)

}


