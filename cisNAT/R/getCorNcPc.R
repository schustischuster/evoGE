# Find overlapping protein-coding coding SAS gene pairs
# Data input: 1) GTF file | 2) Expression_data WITHOUT mito and chloroplast genes
# Analysis can be performed on both whole single species datasets (ATH: 132 samples; AL: 36 samples)
# OR on comparative data sets (27 samples)
# Input sample tables should have the following format:
# id / biotype / source / DEVSEQ_SAMPLE_REPLICATES(between 27 and 132 depending on species)
# GTF annotation: 
# 1. genes
# 2. transcripts
# 3. exon
# 4. CDS
# source and gene_source can be 'araport11' or 'devseq'
# feature.type can be 'gene' / 'transcript' / 'exon' / 'CDS'
# strand can be '+' or '-'
# gene_biotype can be 'protein_coding', 'lnc_exonic_antisense', 'lnc_intronic_antisense', 'lnc_intergenic', 
# 'snoRNA' 'tRNA', 'miRNA' etc.
#---------- for genes ---------
# chromosome / source / feature.type(gene) / start_position / end_position / ??? / strand / gene_id /
# gene_source / gene_biotype / gene_name
#---------- for transcript ---------
# chromosome / source / feature.type(transcript) / start_position / end_position / ??? / strand / gene_id / 
# gene_source / gene_biotype / transcript_id / transcript_source(araport11/DevSeq) / 
# transcript_biotype / gene_name / info
#---------- for exon ---------
# chromosome / source / feature.type(exon) / start_position / end_position / ??? / strand / gene_id / 
# gene_source / gene_biotype / transcript_id / transcript_source (araport11/DevSeq) / 
# transcript_biotype / exon_number / exon_id / gene_name
#---------- for CDS ---------
# chromosome / source / feature.type(CDS) / start_position / end_position / ??? / strand / gene_id / 
# gene_source / gene_biotype / transcript_id / transcript_source (araport11/DevSeq) / 
# transcript_biotype / exon_number / gene_name



#------------------- Load packages, set directories and read sample tables ---------------------


# Install and load packages
if (!require(dplyr)) install.packages('dplyr')
library(dplyr)
if (!require(GenomicRanges)) install.packages('GenomicRanges')
library(GenomicRanges)
if (!require(rtracklayer)) install.packages('rtracklayer')
library(rtracklayer)


# Set file path and input files
in_dir <- "/Volumes/User/Shared/Christoph_manuscript/DevSeq_paper/Analysis/Analysis_2019/A_thaliana_gene_exression_map/20191121_CS_coding_cisNAT_analysis/data"
out_dir <- "/Volumes/User/Shared/Christoph_manuscript/DevSeq_paper/Analysis/Analysis_2019/A_thaliana_gene_exression_map/20191121_CS_coding_cisNAT_analysis"



getCorNcPc <- function(species = c("AT", "AL", "CR", "ES", "TH", "MT", "BD"), 
	experiment = c("single-species", "comparative")) {
	
	# Show error message if no species is chosen
    if (missing(species))
   
       stop(
       "Please choose one of the available species: 
	   'AT', 'AL', 'CR', 'ES', 'TH', 'MT', 'BD'",
	   call. = TRUE
       )

   	# Show error message for ATH and AL if no experiment is chosen
    if (missing(experiment) && (is.element("AT", species) | is.element("AL", species)))
   
       stop(
       "Please choose one of the available experiments: 
	   'single-species', 'comparative'",
	   call. = TRUE
       )

   	if (!is.element(experiment, c("comparative", "single-species")) 
   		&& (is.element("AT", species) | is.element("AL", species)))

   		stop(
       "Please choose one of the available experiments: 
	   'single-species', 'comparative'",
	   call. = TRUE
       )

    # Show error message for CR/ES/TH/MT/BD if no experiment and no species are chosen
    if (missing(experiment) && (!is.element(species, c("CR", "ES", "TH", "MT", "BD"))))
   
       stop(
       "Please choose one of the available species: 
	   'AT', 'AL', 'CR', 'ES', 'TH', 'MT', 'BD'",
	   call. = TRUE
       )


	# Set GTF input gtf file
    if (is.element("AT", species)) {
    	GTFfile = file.path(in_dir, "GTF", "AT_final_annotation.gtf")
        genesTPM = file.path(in_dir, "Expression_data", "AT_lnc_all_antisense_genes_inter_tpm_mat_deseq_sample_ids_extended.csv")
        species_id <- "AT"

    } else if (is.element("AL", species)) {
		GTFfile = file.path(in_dir, "GTF", "AL_final_annotation.gtf")
		genesTPM = file.path(in_dir, "Expression_data", "AL_lnc_all_antisense_genes_inter_tpm_mat_deseq_sample_ids_extended.csv")
		species_id <- "AL"

    } else if (is.element("CR", species)) {
		GTFfile = file.path(in_dir, "GTF", "CR_final_annotation.gtf")
		genesTPM = file.path(in_dir, "Expression_data", "CR_lnc_all_antisense_genes_inter_tpm_mat_deseq_sample_ids_extended.csv")
		species_id <- "CR"

    } else if (is.element("ES", species)) {
		GTFfile = file.path(in_dir, "GTF", "ES_final_annotation.gtf")
		genesTPM = file.path(in_dir, "Expression_data", "ES_lnc_all_antisense_genes_inter_tpm_mat_deseq_sample_ids_extended.csv")
		species_id <- "ES"

    } else if (is.element("TH", species)) {
		GTFfile = file.path(in_dir, "GTF", "TH_final_annotation.gtf")
		genesTPM = file.path(in_dir, "Expression_data", "TH_lnc_all_antisense_genes_inter_tpm_mat_deseq_sample_ids_extended.csv")
		species_id <- "TH"

    } else if (is.element("MT", species)) {
		GTFfile = file.path(in_dir, "GTF", "MT_final_annotation.gtf")
		genesTPM = file.path(in_dir, "Expression_data", "MT_lnc_all_antisense_genes_inter_tpm_mat_deseq_sample_ids_extended.csv")
		species_id <- "MT"

    } else if (is.element("BD", species)) {
		GTFfile = file.path(in_dir, "GTF", "BD_final_annotation.gtf")
		genesTPM = file.path(in_dir, "Expression_data", "BD_lnc_all_antisense_genes_inter_tpm_mat_deseq_sample_ids_extended.csv")
		species_id <- "BD"
    }


	# Import gtf file
	GTF = import.gff(GTFfile, format = "gtf", feature.type = "gene")

	# Read expression data
	all_genes_tpm <- read.table(genesTPM, sep = ";", dec = ".", header = TRUE, stringsAsFactors = FALSE)
	colnames(all_genes_tpm)[1] <- "gene_id"


	# Extract expression data for comparative samples (w/o pollen)
    if ((is.element("AT", species)) && (is.element("comparative", experiment))) {

		all_genes_tpm <- dplyr::select(all_genes_tpm, c(
			gene_id, prt_id, biotype, source, info, 
			Sample_106, # root 5d.1
			Sample_107, # root 5d.2
			Sample_108, # root 5d.3
			Sample_16, # hypocotyl 10d.1
			Sample_17, # hypocotyl 10d.2
			Sample_18, # hypocotyl 10d.3
			Sample_31, # leaf 1+2 7d.1
			Sample_32, # leaf 1+2 7d.2
			Sample_33, # leaf 1+2 7d.3
			Sample_115, # apex veg 7d.1
			Sample_116, # apex veg 7d.2
			Sample_117, # apex veg 7d.3
			Sample_55, # apex inf 21d.1
			Sample_56, # apex inf 21d.2
			Sample_57, # apex inf 21d.3
			Sample_70, # flower stg12 21d.1
			Sample_71, # flower stg12 21d.2
			Sample_72, # flower stg12 21d.3
			Sample_82, # flower stg12 stamen.1
			Sample_83, # flower stg12 stamen.2
			Sample_84, # flower stg12 stamen.3
			Sample_283, # flower stg12 carpel.1
			Sample_284, # flower stg12 carpel.2
			Sample_285)) # flower stg12 carpel.3
		### tibble w/o pollen samles

		species_id <- "AT_comparative_samples"


    } else if ((is.element("AL", species)) && (is.element("comparative", experiment))) {

		all_genes_tpm <- dplyr::select(all_genes_tpm, -c(
			Sample_295, # flower stg11 stamen.1 
			Sample_296, # flower stg11 stamen.2 
			Sample_297, # flower stg11 stamen.3
			Sample_298, # flower early stg12 stamen.1
			Sample_299, # flower early stg12 stamen.2
			Sample_300, # flower early stg12 stamen.3
			Sample_301, # flower late stg12 stamen.1
			Sample_302, # flower late stg12 stamen.2
			Sample_303)) # flower late stg12 stamen.3
		### tibble w/o pollen samles

		species_id <- "AL_comparative_samples"

	}


    # Stop function here to allow specific analysis of a single species
    # return_list <- list("species_id" = species_id, "GTF" = GTF, "all_genes_tpm" = all_genes_tpm, "experiment" = experiment)
    # return(return_list)
    # }
    # return_objects <- getCorNcPc("AT", "single-species") # read in GTF and expression data for A.thaliana
    # list2env(return_objects, envir = .GlobalEnv)



	#------------------- Compute correlation between coding gene and cisNAT  -------------------

	if (species_id == "AT_comparative_samples" | species_id == "AL_comparative_samples") {

		all_genes_tpm <- all_genes_tpm[rowSums(all_genes_tpm[, 6:ncol(all_genes_tpm)]) > 0, ]

	} # Remove all genes that are not expressed in comparative samples

	all_genes_tpm_ls <- split(all_genes_tpm, f = all_genes_tpm$prt_id) 
	# Split data by same protein-coding gene id
	# 4103 protein-coding gene/cisNAT pairs for AT

	all_genes_tpm_ls <- Filter(function(dt) nrow(dt) > 1, all_genes_tpm_ls)
	# Remove all list entries that contain unpaired entries

	length(all_genes_tpm_ls[purrr::map(all_genes_tpm_ls, nrow) > 2])
	# check how many genes have more than one NAT overlapping = 152 for AT

	length(all_genes_tpm_ls[purrr::map(all_genes_tpm_ls, nrow) > 3])
	# 7 genes overlapping 3 cisNATs for AT


	# Compute spearman and pearson correlations
	getCdNcCor <- function(df) {

		df[, 6:ncol(df)] <- log2(df[, 6:ncol(df)] + 1) # log-transform TPM

		if (nrow(df) == 2) {

			cor <- cor.test(as.numeric(df[1, 6:ncol(df)]), as.numeric(df[2, 6:ncol(df)]))$estimate
			max <- apply(X = df[6:ncol(df)], MARGIN = 1, FUN = max)
			avg <- apply(X = df[6:ncol(df)], MARGIN = 1, FUN = mean)
			df_out <- data.frame(df[-1, 1:5], cor = c(cor), maxPC = max[1], maxNC = max[2], 
				meanPC = avg[1], meanNC = avg[2], maxRatio = max[2]/max[1], meanRatio = avg[2]/avg[1])
		
		} else if (nrow(df) == 3) {

			cor1 <- cor.test(as.numeric(df[1, 6:ncol(df)]), as.numeric(df[2, 6:ncol(df)]))$estimate
			cor2 <- cor.test(as.numeric(df[1, 6:ncol(df)]), as.numeric(df[3, 6:ncol(df)]))$estimate
			max <- apply(X = df[6:ncol(df)], MARGIN = 1, FUN = max)
			avg <- apply(X = df[6:ncol(df)], MARGIN = 1, FUN = mean)
			df_out <- data.frame(df[-1, 1:5], cor = c(cor1, cor2), maxPC = rep(max[1], 2), maxNC = c(max[2], max[3]), 
				meanPC = rep(avg[1], 2), meanNC = c(avg[2], avg[3]), maxRatio = c(max[2]/max[1], max[3]/max[1]), 
				meanRatio = c(avg[2]/avg[1], avg[3]/avg[1]))
		
		} else if (nrow(df) == 4) {

			cor1 <- cor.test(as.numeric(df[1, 6:ncol(df)]), as.numeric(df[2, 6:ncol(df)]))$estimate
			cor2 <- cor.test(as.numeric(df[1, 6:ncol(df)]), as.numeric(df[3, 6:ncol(df)]))$estimate
			cor3 <- cor.test(as.numeric(df[1, 6:ncol(df)]), as.numeric(df[4, 6:ncol(df)]))$estimate
			max <- apply(X = df[6:ncol(df)], MARGIN = 1, FUN = max)
			avg <- apply(X = df[6:ncol(df)], MARGIN = 1, FUN = mean)
			df_out <- data.frame(df[-1, 1:5], cor = c(cor1, cor2, cor3), maxPC = rep(max[1], 3), maxNC = c(max[2], max[3], max[4]), 
				meanPC = rep(avg[1], 3), meanNC = c(avg[2], avg[3], avg[4]), maxRatio = c(max[2]/max[1], max[3]/max[1], max[4]/max[1]), 
				meanRatio = c(avg[2]/avg[1], avg[3]/avg[1], avg[4]/avg[1]))
		
		}

		return(df_out)
	}

	# startup message
	message("Computing correlation...")

	cd_nc_cor_ls <- lapply(all_genes_tpm_ls, getCdNcCor)
	cd_nc_cor <- do.call("rbind", cd_nc_cor_ls)




#------------------------- Create final data table and write csv file  -------------------------


	# Create data table containing both strand plus and minus genes and cor values

	strand_minus_overlap_genes_descript = strand_minus_overlap_genes_tpm_th %>% select(
		gene_id, id, gene_source, gene_biotype, seqnames, start, end, strand, width, width_overlap)
	names(strand_minus_overlap_genes_descript) <- c(
		"id_minus_strand", "id_minus", "gene_source_subject", "biotype_subject", "seqnames", 
		"start_minus", "end_minus", "strand_subject", "width_subject", "width_overlap_subject")
	strand_plus_overlap_genes <- dplyr::select(strand_plus_overlap_genes, -c(seqnames))
	names(strand_plus_overlap_genes)[1:9] <- c("id_plus_strand", "id_plus", "gene_source_query", 
		"biotype_query", "start_plus", "end_plus", "strand_query", "width_query", "width_overlap_query")
	overlap_cd_nc_genes_cor <- cbind(strand_plus_overlap_genes, strand_minus_overlap_genes_descript)
	overlap_cd_nc_genes_cor$NAT_overlap_width <- 
		overlap_cd_nc_genes_cor$width_overlap_query + overlap_cd_nc_genes_cor$width_overlap_subject
	overlap_cd_nc_genes_cor = overlap_cd_nc_genes_cor %>% select(
		id_plus_strand, 
		id_plus, 
		seqnames, 
		start_plus, 
		end_plus, 
		width_query, 
		strand_query, 
		biotype_query, 
		gene_source_query,
		id_minus_strand, 
		id_minus, 
		start_minus, 
		end_minus, 
		width_subject, 
		strand_subject, 
		biotype_subject, 
		gene_source_subject, 
		NAT_overlap_width, 
		Spearman, 
		Pearson)


	strand_minus_overlap_genes_descript_wo_pollen = strand_minus_overlap_genes_tpm_th_wo_pollen %>% select(
		gene_id, id, gene_source, gene_biotype, seqnames, start, end, strand, width, width_overlap)
	names(strand_minus_overlap_genes_descript_wo_pollen) <- c(
		"id_minus_strand", "id_minus", "gene_source_subject", "biotype_subject", "seqnames", 
		"start_minus", "end_minus", "strand_subject", "width_subject", "width_overlap_subject")
	strand_plus_overlap_genes_wo_pollen <- dplyr::select(strand_plus_overlap_genes_wo_pollen, -c(seqnames))
	names(strand_plus_overlap_genes_wo_pollen)[1:9] <- c("id_plus_strand", "id_plus", "gene_source_query", 
		"biotype_query", "start_plus", "end_plus", "strand_query", "width_query", "width_overlap_query")
	overlap_cd_nc_genes_cor_wo_pollen <- cbind(strand_plus_overlap_genes_wo_pollen, 
		strand_minus_overlap_genes_descript_wo_pollen)

	overlap_cd_nc_genes_cor_wo_pollen$NAT_overlap_width <- 
		overlap_cd_nc_genes_cor_wo_pollen$width_overlap_query + overlap_cd_nc_genes_cor_wo_pollen$width_overlap_subject
	
	overlap_cd_nc_genes_cor_wo_pollen = overlap_cd_nc_genes_cor_wo_pollen %>% select(
		id_plus_strand, 
		id_plus, 
		seqnames, 
		start_plus, 
		end_plus, 
		width_query, 
		strand_query, 
		biotype_query, 
		gene_source_query,
		id_minus_strand, 
		id_minus, 
		start_minus, 
		end_minus, 
		width_subject, 
		strand_subject, 
		biotype_subject, 
		gene_source_subject, 
		NAT_overlap_width, 
		Spearman, 
		Pearson)


	# Remove coding-coding SAS pairs from data
	overlap_cd_nc_genes_cor <- overlap_cd_nc_genes_cor %>% filter(biotype_query != biotype_subject)
	overlap_cd_nc_genes_cor_wo_pollen <- overlap_cd_nc_genes_cor_wo_pollen %>% filter(biotype_query != biotype_subject)


	# Create a list of all lncRNAs above the threshold (dataset including pollen samples)
	overlap_nc_gene_ids_plus = overlap_cd_nc_genes_cor %>% select(id_plus_strand, biotype_query)
	names(overlap_nc_gene_ids_plus) <- c("id", "biotype")
	overlap_nc_gene_ids_minus = overlap_cd_nc_genes_cor %>% select(id_minus_strand, biotype_subject)
	names(overlap_nc_gene_ids_minus) <- c("id", "biotype")
	NATs_above_threshold <- rbind(overlap_nc_gene_ids_plus, overlap_nc_gene_ids_minus)
	NATs_above_threshold <- filter(NATs_above_threshold, biotype != "protein_coding")
	NATs_above_threshold <- unique(NATs_above_threshold)


	# Set filename
    fname <- sprintf('%s.csv', paste(species_id, "cd_nc_SAS_cor", threshold, sep="_"))
	fname_wo_pollen <- sprintf('%s.csv', paste(species_id, "cd_nc_SAS_cor_wo_pollen", threshold, sep="_"))
	fname_NATs_above_th <- sprintf('%s.csv', paste(species_id, "NATs_above", threshold, "TPM", sep="_"))


	# Write final data tables to csv files and store them in /out_dir/output/data_tables
	if (!dir.exists(file.path(out_dir, "output", "overlap_nc_genes"))) 
		dir.create(file.path(out_dir, "output", "overlap_nc_genes"), recursive = TRUE)
	message("Storing results in: ", file.path("output", "overlap_nc_genes"))

	write.table(overlap_cd_nc_genes_cor, file=file.path(out_dir, "output", "overlap_nc_genes", fname), 
		sep=";", dec=".", row.names=FALSE, col.names=TRUE)
	write.table(overlap_cd_nc_genes_cor_wo_pollen, file=file.path(out_dir, "output", "overlap_nc_genes", fname_wo_pollen), 
		sep=";", dec=".", row.names=FALSE, col.names=TRUE)
	write.table(NATs_above_threshold, file=file.path(out_dir, "output", "overlap_nc_genes", fname_NATs_above_th), 
		sep=";", dec=".", row.names=FALSE, col.names=TRUE)




#-- Create data tables containing expression values for all species w/o pollen and write csv file ---


	# Create data table containing both strand plus and minus genes, their expression data and cor values
	if (threshold==0.5) {

		strand_minus_overlap_genes_wo_pollen <- getCor(
			strand_minus_overlap_genes_tpm_th_wo_pollen, strand_plus_overlap_genes_tpm_th_wo_pollen)

		strand_plus_overlap_genes_descript_wo_pollen = strand_plus_overlap_genes_tpm_th_wo_pollen %>% select(
			gene_id, id, gene_source, gene_biotype, seqnames, start, end, strand, width, width_overlap)
		names(strand_plus_overlap_genes_descript_wo_pollen) <- c(
			"id_plus_strand", "id_plus", "gene_source_query", "biotype_query", "seqnames", 
			"start_plus", "end_plus", "strand_query", "width_query", "width_overlap_query")
		strand_minus_overlap_genes_wo_pollen <- dplyr::select(strand_minus_overlap_genes_wo_pollen, -c(seqnames))
		names(strand_minus_overlap_genes_wo_pollen)[1:9] <- c("id_minus_strand", "id_minus", "gene_source_subject", 
			"biotype_subject", "start_minus", "end_minus", "strand_subject", "width_subject", "width_overlap_subject")
		overlap_cd_nc_genes_cor_wo_pollen_subject_expr <- cbind(strand_minus_overlap_genes_wo_pollen, 
			strand_plus_overlap_genes_descript_wo_pollen)

		overlap_cd_nc_genes_cor_wo_pollen_subject_expr$NAT_overlap_width <- 
			overlap_cd_nc_genes_cor_wo_pollen_subject_expr$width_overlap_query + overlap_cd_nc_genes_cor_wo_pollen_subject_expr$width_overlap_subject
	
		overlap_cd_nc_genes_cor_wo_pollen_subject_expr = overlap_cd_nc_genes_cor_wo_pollen_subject_expr %>% select(
			id_minus_strand, 
			id_minus, 
			seqnames, 
			start_minus, 
			end_minus, 
			width_subject, 
			strand_subject, 
			biotype_subject, 
			gene_source_subject,
			id_plus_strand, 
			id_plus, 
			start_plus, 
			end_plus, 
			width_query, 
			strand_query, 
			biotype_query, 
			gene_source_query, 
			NAT_overlap_width, 
			Spearman, 
			Pearson,
			everything())

		overlap_cd_nc_genes_cor_wo_pollen_subject_expr	<- dplyr::select(
			overlap_cd_nc_genes_cor_wo_pollen_subject_expr, -c(width_overlap_subject, width_overlap_query))

		overlap_cd_nc_genes_cor_wo_pollen_query_expr <- cbind(strand_plus_overlap_genes_wo_pollen, 
			strand_minus_overlap_genes_descript_wo_pollen)

		overlap_cd_nc_genes_cor_wo_pollen_query_expr$NAT_overlap_width <- 
			overlap_cd_nc_genes_cor_wo_pollen_query_expr$width_overlap_query + overlap_cd_nc_genes_cor_wo_pollen_query_expr$width_overlap_subject

		overlap_cd_nc_genes_cor_wo_pollen_query_expr = overlap_cd_nc_genes_cor_wo_pollen_query_expr %>% select(
			id_plus_strand, 
			id_plus, 
			seqnames, 
			start_plus, 
			end_plus, 
			width_query, 
			strand_query, 
			biotype_query, 
			gene_source_query,
			id_minus_strand, 
			id_minus, 
			start_minus, 
			end_minus, 
			width_subject, 
			strand_subject, 
			biotype_subject, 
			gene_source_subject, 
			NAT_overlap_width, 
			Spearman, 
			Pearson,
			everything())

		overlap_cd_nc_genes_cor_wo_pollen_query_expr <- dplyr::select(
			overlap_cd_nc_genes_cor_wo_pollen_query_expr, -c(width_overlap_subject, width_overlap_query))

		# Remove coding-coding SAS pairs from data
		overlap_cd_nc_genes_cor_wo_pollen_query_expr <- overlap_cd_nc_genes_cor_wo_pollen_query_expr %>% filter(biotype_query != biotype_subject)
		overlap_cd_nc_genes_cor_wo_pollen_subject_expr <- overlap_cd_nc_genes_cor_wo_pollen_subject_expr %>% filter(biotype_subject != biotype_query)

		# Set filename
		fname_wo_pollen_query <- sprintf('%s.csv', paste(species_id, "cd_nc_SAS_cor_wo_pollen_query_expr", threshold, sep="_"))
		fname_wo_pollen_subject <- sprintf('%s.csv', paste(species_id, "cd_nc_SAS_cor_wo_pollen_subject_expr", threshold, sep="_"))

		# Write data tables to csv files and store them in /out_dir/output/data_tables
		write.table(overlap_cd_nc_genes_cor_wo_pollen_query_expr, file=file.path(out_dir, "output", "overlap_nc_genes", fname_wo_pollen_query), 
			sep=";", dec=".", row.names=FALSE, col.names=TRUE)
		write.table(overlap_cd_nc_genes_cor_wo_pollen_subject_expr, file=file.path(out_dir, "output", "overlap_nc_genes", fname_wo_pollen_subject), 
			sep=";", dec=".", row.names=FALSE, col.names=TRUE)
	}

}


thresholds <- list(0.5, 2, 5, 10)

lapply(thresholds, getNcPc, species = "ATH", experiment = "single-species")
lapply(thresholds, getNcPc, species = "ATH", experiment = "comparative")
lapply(thresholds, getNcPc, species = "AL", experiment = "comparative")
lapply(thresholds, getNcPc, species = "CR", experiment = "comparative")
lapply(thresholds, getNcPc, species = "ES", experiment = "comparative")
lapply(thresholds, getNcPc, species = "TH", experiment = "comparative")
lapply(thresholds, getNcPc, species = "MT", experiment = "comparative")
lapply(thresholds, getNcPc, species = "BD", experiment = "comparative")

getNcPc("ATH", "single-species", 0)

