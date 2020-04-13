# Find expressed genes in each sample replicate using ERCC spike-in information
# Thresholds: 0 (no ERCC threshold-everything above 0 TPM), 0.01 (weakest expressed spike-in), 
# 0.05 (lowest 5% of expressed spike-ins), 0.1 (lowest 10% of expressed spike-ins)
# Data input: 1) Expression_data WITH spike-in information
# Analysis can be performed on both whole single species datasets (ATH: 132 samples; AL: 36 samples)
# OR on comparative data sets (27 samples)
# Input sample tables should have the following format:
# id / biotype / source / DEVSEQ_SAMPLE_REPLICATES(between 27 and 132 depending on species)




#------------------- Load packages, set directories and read sample tables ---------------------


# Install and load packages
if (!require(dplyr)) install.packages('dplyr')
library(dplyr)
if (!require(data.table)) install.packages('data.table')
library(data.table)


# Set file path and input files
in_dir <- "/Volumes/User/Shared/Christoph_manuscript/DevSeq_paper/Analysis/Analysis_2019/A_thaliana_gene_exression_map/20200113_expressed_genes_CS/data"
out_dir <- "/Volumes/User/Shared/Christoph_manuscript/DevSeq_paper/Analysis/Analysis_2019/A_thaliana_gene_exression_map/20200113_expressed_genes_CS"





# Define function to get expressed genes

getExprGenes <- function(species = c("ATH", "AL", "CR", "ES", "TH", "MT", "BD"), 
	experiment = c("single-species", "comparative"), threshold) {
	
	# Show error message if no species is chosen
    if (missing(species))
   
       stop(
       "Please choose one of the available species: 
	   'ATH', 'AL', 'CR', 'ES', 'TH', 'MT', 'BD'",
	   call. = TRUE
       )

   	# Show error message for ATH and AL if no experiment is chosen
    if (missing(experiment) && (is.element("ATH", species) | is.element("AL", species)))
   
       stop(
       "Please choose one of the available experiments: 
	   'single-species', 'comparative'",
	   call. = TRUE
       )

   	if (!is.element(experiment, c("comparative", "single-species")) 
   		&& (is.element("ATH", species) | is.element("AL", species)))

   		stop(
       "Please choose one of the available experiments: 
	   'single-species', 'comparative'",
	   call. = TRUE
       )

    # Show error message for CR/ES/TH/MT/BD if no experiment and no species are chosen
    if (missing(experiment) && (!is.element(species, c("CR", "ES", "TH", "MT", "BD"))))
   
       stop(
       "Please choose one of the available species: 
	   'ATH', 'AL', 'CR', 'ES', 'TH', 'MT', 'BD'",
	   call. = TRUE
       )

    # Show error message if no threshold is given
    if (missing(threshold))
   
       stop(
       "Please set a spike-in threshold:
       e.g. '0.05' for the lowest 5% of expressed spike-ins",
	   call. = TRUE
       )

    # Add an error if threshold < 0
  	if (threshold < 0)
    	stop(
        "'threshold' must be >= 0",
	   	call. = TRUE
    	)


    # Show startup message
    message("Reading data...")


	# Set GTF input gtf file
    if (is.element("ATH", species)) {
        genesTPM = file.path(in_dir, "Expression_data", "ATH_no_TE_genes_tpm_sample_names.csv")
        species_id <- "ATH"

    } else if (is.element("AL", species)) {
		genesTPM = file.path(in_dir, "Expression_data", "AL_no_TE_genes_tpm_sample_names.csv")
		species_id <- "AL"

    } else if (is.element("CR", species)) {
		genesTPM = file.path(in_dir, "Expression_data", "CR_no_TE_genes_tpm_sample_names.csv")
		species_id <- "CR"

    } else if (is.element("ES", species)) {
		genesTPM = file.path(in_dir, "Expression_data", "ES_no_TE_genes_tpm_sample_names.csv")
		species_id <- "ES"

    } else if (is.element("TH", species)) {
		genesTPM = file.path(in_dir, "Expression_data", "TH_no_TE_genes_tpm_sample_names.csv")
		species_id <- "TH"

    } else if (is.element("MT", species)) {
		genesTPM = file.path(in_dir, "Expression_data", "MT_no_TE_genes_tpm_sample_names.csv")
		species_id <- "MT"

    } else if (is.element("BD", species)) {
		genesTPM = file.path(in_dir, "Expression_data", "BD_no_TE_genes_tpm_sample_names.csv")
		species_id <- "BD"
    }


	# Read expression data
	all_genes_tpm <- read.table(genesTPM, sep=";", dec=".", header=TRUE, stringsAsFactors=FALSE)
	colnames(all_genes_tpm)[1] <- "gene_id"


	# Format expression data and rename pollen samples
    if ((is.element("ATH", species)) && (is.element("comparative", experiment))) {

		all_genes_tpm <- dplyr::select(all_genes_tpm, c(
			gene_id, biotype, source, 
			root_whole_root_5d_1,
			root_whole_root_5d_2,
			root_whole_root_5d_3,
			hypocotyl_10d_1,
			hypocotyl_10d_2,
			hypocotyl_10d_3,
			leaf_1.2_7d_1,
			leaf_1.2_7d_2,
			leaf_1.2_7d_3,
			apex_vegetative_7d_1,
			apex_vegetative_7d_2,
			apex_vegetative_7d_3,
			apex_inflorescence_21d_1,
			apex_inflorescence_21d_2,
			apex_inflorescence_21d_3,
			flower_stg12_21d._1,
			flower_stg12_21d._2,
			flower_stg12_21d._3,
			flower_stg12_stamens_21d._1,
			flower_stg12_stamens_21d._2,
			flower_stg12_stamens_21d._3,
			flowers_mature_pollen_28d_1,
			flowers_mature_pollen_28d_2,
			flowers_mature_pollen_28d_3,
			flower_early_stg12_carpels_21d._1,
			flower_early_stg12_carpels_21d._2,
			flower_early_stg12_carpels_21d._3)) #tibble w/o pollen samles

		colnames(all_genes_tpm)[colnames(all_genes_tpm) == 'flowers_mature_pollen_28d_1'] <- 'flowers_mature_pollen_1'
		colnames(all_genes_tpm)[colnames(all_genes_tpm) == 'flowers_mature_pollen_28d_2'] <- 'flowers_mature_pollen_2'
		colnames(all_genes_tpm)[colnames(all_genes_tpm) == 'flowers_mature_pollen_28d_3'] <- 'flowers_mature_pollen_3'

		species_id <- "ATH_comparative_samples"


    } else if ((is.element("AL", species)) && (is.element("comparative", experiment))) {

		all_genes_tpm <- dplyr::select(all_genes_tpm, -c(
			flower_stg11_stamens_8w.10w.25d_1, 
			flower_stg11_stamens_8w.10w.25d_2, 
			flower_stg11_stamens_8w.10w.25d_3,
			flower_early_stg12_stamens_8w.10w.23d_1,
			flower_early_stg12_stamens_8w.10w.23d_2,
			flower_early_stg12_stamens_8w.10w.23d_3,
			flower_late_stg12_stamens_8w.10w.21d_1,
			flower_late_stg12_stamens_8w.10w.21d_2,
			flower_late_stg12_stamens_8w.10w.21d_3)) #tibble w/o pollen samles

		colnames(all_genes_tpm)[colnames(all_genes_tpm) == 'flowers_mature_pollen_8w.10w.25d_1'] <- 'flowers_mature_pollen_1'
		colnames(all_genes_tpm)[colnames(all_genes_tpm) == 'flowers_mature_pollen_8w.10w.25d_2'] <- 'flowers_mature_pollen_2'
		colnames(all_genes_tpm)[colnames(all_genes_tpm) == 'flowers_mature_pollen_8w.10w.25d_3'] <- 'flowers_mature_pollen_3'

		species_id <- "AL_comparative_samples"


    } else if (is.element("ATH", species)) {
		colnames(all_genes_tpm)[colnames(all_genes_tpm) == 'flowers_mature_pollen_28d_1'] <- 'flowers_mature_pollen_1'
		colnames(all_genes_tpm)[colnames(all_genes_tpm) == 'flowers_mature_pollen_28d_2'] <- 'flowers_mature_pollen_2'
		colnames(all_genes_tpm)[colnames(all_genes_tpm) == 'flowers_mature_pollen_28d_3'] <- 'flowers_mature_pollen_3'

    } else if (is.element("AL", species)) {
    	colnames(all_genes_tpm)[colnames(all_genes_tpm) == 'flowers_mature_pollen_8w.10w.25d_1'] <- 'flowers_mature_pollen_1'
		colnames(all_genes_tpm)[colnames(all_genes_tpm) == 'flowers_mature_pollen_8w.10w.25d_2'] <- 'flowers_mature_pollen_2'
		colnames(all_genes_tpm)[colnames(all_genes_tpm) == 'flowers_mature_pollen_8w.10w.25d_3'] <- 'flowers_mature_pollen_3'

    } else if (is.element("CR", species)) {
    	colnames(all_genes_tpm)[colnames(all_genes_tpm) == 'flowers_mature_pollen_6w.7w.22d_1'] <- 'flowers_mature_pollen_1'
		colnames(all_genes_tpm)[colnames(all_genes_tpm) == 'flowers_mature_pollen_6w.7w.22d_2'] <- 'flowers_mature_pollen_2'
		colnames(all_genes_tpm)[colnames(all_genes_tpm) == 'flowers_mature_pollen_6w.7w.22d_3'] <- 'flowers_mature_pollen_3'

    } else if (is.element("ES", species)) {
		colnames(all_genes_tpm)[colnames(all_genes_tpm) == 'flowers_mature_pollen_7w.8w.17d_1'] <- 'flowers_mature_pollen_1'
		colnames(all_genes_tpm)[colnames(all_genes_tpm) == 'flowers_mature_pollen_7w.8w.17d_2'] <- 'flowers_mature_pollen_2'
		colnames(all_genes_tpm)[colnames(all_genes_tpm) == 'flowers_mature_pollen_7w.8w.17d_3'] <- 'flowers_mature_pollen_3'

    } else if (is.element("TH", species)) {
		colnames(all_genes_tpm)[colnames(all_genes_tpm) == 'flowers_mature_pollen_11w_1'] <- 'flowers_mature_pollen_1'
		colnames(all_genes_tpm)[colnames(all_genes_tpm) == 'flowers_mature_pollen_11w_2'] <- 'flowers_mature_pollen_2'
		colnames(all_genes_tpm)[colnames(all_genes_tpm) == 'flowers_mature_pollen_11w_3'] <- 'flowers_mature_pollen_3'

    } else if (is.element("MT", species)) {
		colnames(all_genes_tpm)[colnames(all_genes_tpm) == 'flowers_mature_pollen_7w_1'] <- 'flowers_mature_pollen_1'
		colnames(all_genes_tpm)[colnames(all_genes_tpm) == 'flowers_mature_pollen_7w_2'] <- 'flowers_mature_pollen_2'
		colnames(all_genes_tpm)[colnames(all_genes_tpm) == 'flowers_mature_pollen_7w_3'] <- 'flowers_mature_pollen_3'

    } else if (is.element("BD", species)) {
		colnames(all_genes_tpm)[colnames(all_genes_tpm) == 'flowers_mature_pollen_32d_1'] <- 'flowers_mature_pollen_1'
		colnames(all_genes_tpm)[colnames(all_genes_tpm) == 'flowers_mature_pollen_32d_2'] <- 'flowers_mature_pollen_2'
		colnames(all_genes_tpm)[colnames(all_genes_tpm) == 'flowers_mature_pollen_32d_3'] <- 'flowers_mature_pollen_3'
    }



    # Stop function here to allow specific analysis of a single species
    # return_list <- list("species_id" = species_id, "all_genes_tpm" = all_genes_tpm, "threshold" = threshold)
    # return(return_list)
    # }
    # return_objects <- getExprGenes(species="ATH", experiment="single-species", 0.05) # read in GTF and expression data for A.thaliana
    # list2env(return_objects, envir = .GlobalEnv)



	# Extract ERCC data
	ERCC <- all_genes_tpm[all_genes_tpm$gene_id %like% "ERCC", ]


	# Get lowest ERCC spike-in after applying threshold
	ERCC_threshold <- rbind(ERCC, cbind(
		data.frame(gene_id="ERCC_threshold"), data.frame(biotype="<NA>"), data.frame(source="<NA>"), 
		as.data.frame(t(
			data.frame(ERCC_threshold = sapply(ERCC[,4:ncol(ERCC)], function(x) {
				ERCC_cutoff <- floor(sum(x>0) - (threshold * sum(x>0))) 
				# use ceiling() instead floor() to round up if neccessary
				ERCC_cutoff_low <- nrow(ERCC) - ERCC_cutoff
				return(ERCC_cutoff_low)
				}
			))
		))
	))


	# Get sample-specific TPM threshold
	ERCC_cutoff <- rbind(ERCC_threshold, cbind(
		data.frame(gene_id="TPM_cutoff"), data.frame(biotype="<NA>"), data.frame(source="<NA>"), 
		as.data.frame(t(
			data.frame(TPM_cutoff = sapply(ERCC_threshold[,4:ncol(ERCC_threshold)], function(x) {
				sort(x[1:nrow(ERCC_threshold)-1], partial=x[nrow(ERCC_threshold)])[x[nrow(ERCC_threshold)]]
				}
			))
		))
	))


	# Bind sample-specific threshold to expression table
	all_genes_tpm_cutoff <- rbind(all_genes_tpm,ERCC_cutoff[nrow(ERCC_cutoff),])











}




test <- getExprGenes(species="ATH", experiment="single-species", threshold=0.05)














# test data
df <- data.frame (ID  = c('data1', 'data2', 'data3', 'data4', 'TPM_cutoff'), 
	source  = c('DevSeq', 'DevSeq', 'DevSeq', 'DevSeq', 'NA'), 
    sample1 = c(2, 1, 18, 3, 3),
    sample2 = c(4, 3, 17, 16, 5),
    sample3 = c(3, 2, 11, 2, 2),
    sample4 = c(22, 9, 11, 35, 11),
    sample5 = c(10, 5, 8, 22, 7),
    sample6 = c(17, 6, 9, 11, 8))






# Define threshold function
		getThreshold <- function(df) {

			# Split data frame by sample replicates into a list then apply threshold for each subset
	
			th_replicates <- do.call(cbind, lapply(split.default(df[2:ncol(df)], #adjust columns
								rep(seq_along(df), each = 1, length.out = ncol(df)-1)),
								function(x) {
									x[x <= x[nrow(df),], ] <- 0;
									x
								}
							))

			# Bind key/id/prt_id/symbol/biotype/source columns to thresholded data frame
			th_replicates <- cbind(df[1], th_replicates)

			return(th_replicates)
		}









	# This is a threshold function that can be applied to expression tables
	# Settings: TPM > ERCC spike-in threshold in at least 2 of 3 replicates in at least one sample
	applyThreshold <- function(express_df) {

		# Add keys to data frame
		key <- seq(1, nrow(express_df), 1)
		express_df <- cbind(as.data.frame(key),express_df)

		getSampleTH <- function(df) {

			# Split data frame by sample replicates into a list then apply threshold for each subset
	
			th_replicates <- do.call(cbind, lapply(split.default(df[4:ncol(df)], #adjust columns
								rep(seq_along(df), each = 1, length.out = ncol(df)-3)),
								function(x) {
									x[x <= x[nrow(df),], ] <- 0;
									x
								}
							))

			# Bind key/id/prt_id/symbol/biotype/source columns to thresholded data frame
			th_replicates <- cbind(df[1:3], th_replicates)

			return(th_replicates)
		}

		df <- getSampleTH(express_df)

		# Define threshold function
		getThreshold <- function(df) {

			# Split data frame by sample replicates into a list then apply threshold for each subset
	
			th_replicates <- do.call(cbind, lapply(split.default(df[4:ncol(df)], #adjust columns
								rep(seq_along(df), each = 3, length.out = ncol(df)-3)), #adjust columns
								function(x) {
									x[rowSums(x > 0) < 2, ] <- 0; 
									x
								}
							))

			# Bind key/id/prt_id/symbol/biotype/source columns to thresholded data frame
			th_replicates <- cbind(df[1:3], th_replicates)

			# Remove all rows that only contain "0"
			th_replicates <- th_replicates[which(rowSums(th_replicates[,-1:-3, drop = FALSE] > 0) > 0),]

			return(th_replicates)
		}

		# Apply threshold to data and extract keys ("key")
		keys_data <- getThreshold(df)
		keys_data <- keys_data[,1:2]
		names(keys_data) <- c("key","ID")

		# Generate thresholded data frame based on keys
		th_df <- merge(keys_data, express_df, by="key")
		th_df <- th_df[order(th_df$key),]
		th_df <- th_df[-1:-2]

		return(th_df)
	}






















# Define function to get overlapping protein-coding genes

getPcPc <- function(species = c("ATH", "AL", "CR", "ES", "TH", "MT", "BD"), 
	experiment = c("single-species", "comparative")) {
	
	# Show error message if no species is chosen
    if (missing(species))
   
       stop(
       "Please choose one of the available species: 
	   'ATH', 'AL', 'CR', 'ES', 'TH', 'MT', 'BD'",
	   call. = TRUE
       )

   	# Show error message for ATH and AL if no experiment is chosen
    if (missing(experiment) && (is.element("ATH", species) | is.element("AL", species)))
   
       stop(
       "Please choose one of the available experiments: 
	   'single-species', 'comparative'",
	   call. = TRUE
       )

   	if (!is.element(experiment, c("comparative", "single-species")) 
   		&& (is.element("ATH", species) | is.element("AL", species)))

   		stop(
       "Please choose one of the available experiments: 
	   'single-species', 'comparative'",
	   call. = TRUE
       )

    # Show error message for CR/ES/TH/MT/BD if no experiment and no species are chosen
    if (missing(experiment) && (!is.element(species, c("CR", "ES", "TH", "MT", "BD"))))
   
       stop(
       "Please choose one of the available species: 
	   'ATH', 'AL', 'CR', 'ES', 'TH', 'MT', 'BD'",
	   call. = TRUE
       )


	# Set GTF input gtf file
    if (is.element("ATH", species)) {
    	GTFfile = file.path(in_dir, "GTF", "AT_final_annotation.gtf")
        genesTPM = file.path(in_dir, "Expression_data", "ATH_no_TE_genes_tpm_sample_names.csv")
        species_id <- "ATH"

    } else if (is.element("AL", species)) {
		GTFfile = file.path(in_dir, "GTF", "AL_final_annotation.gtf")
		genesTPM = file.path(in_dir, "Expression_data", "AL_no_TE_genes_tpm_sample_names.csv")
		species_id <- "AL"

    } else if (is.element("CR", species)) {
		GTFfile = file.path(in_dir, "GTF", "CR_final_annotation.gtf")
		genesTPM = file.path(in_dir, "Expression_data", "CR_no_TE_genes_tpm_sample_names.csv")
		species_id <- "CR"

    } else if (is.element("ES", species)) {
		GTFfile = file.path(in_dir, "GTF", "ES_final_annotation.gtf")
		genesTPM = file.path(in_dir, "Expression_data", "ES_no_TE_genes_tpm_sample_names.csv")
		species_id <- "ES"

    } else if (is.element("TH", species)) {
		GTFfile = file.path(in_dir, "GTF", "TH_final_annotation.gtf")
		genesTPM = file.path(in_dir, "Expression_data", "TH_no_TE_genes_tpm_sample_names.csv")
		species_id <- "TH"

    } else if (is.element("MT", species)) {
		GTFfile = file.path(in_dir, "GTF", "MT_final_annotation.gtf")
		genesTPM = file.path(in_dir, "Expression_data", "MT_no_TE_genes_tpm_sample_names.csv")
		species_id <- "MT"

    } else if (is.element("BD", species)) {
		GTFfile = file.path(in_dir, "GTF", "BD_final_annotation.gtf")
		genesTPM = file.path(in_dir, "Expression_data", "BD_no_TE_genes_tpm_sample_names.csv")
		species_id <- "BD"
    }


	# Import gtf file
	GTF = import.gff(GTFfile, format="gtf", feature.type="gene")

	# Read expression data
	all_genes_tpm <- read.table(genesTPM, sep=";", dec=".", header=TRUE, stringsAsFactors=FALSE)
	colnames(all_genes_tpm)[1] <- "gene_id"


	# Format expression data and rename pollen samples
    if ((is.element("ATH", species)) && (is.element("comparative", experiment))) {

		all_genes_tpm <- dplyr::select(all_genes_tpm, c(
			gene_id, biotype, source, 
			root_whole_root_5d_1,
			root_whole_root_5d_2,
			root_whole_root_5d_3,
			hypocotyl_10d_1,
			hypocotyl_10d_2,
			hypocotyl_10d_3,
			leaf_1.2_7d_1,
			leaf_1.2_7d_2,
			leaf_1.2_7d_3,
			apex_vegetative_7d_1,
			apex_vegetative_7d_2,
			apex_vegetative_7d_3,
			apex_inflorescence_21d_1,
			apex_inflorescence_21d_2,
			apex_inflorescence_21d_3,
			flower_stg12_21d._1,
			flower_stg12_21d._2,
			flower_stg12_21d._3,
			flower_stg12_stamens_21d._1,
			flower_stg12_stamens_21d._2,
			flower_stg12_stamens_21d._3,
			flowers_mature_pollen_28d_1,
			flowers_mature_pollen_28d_2,
			flowers_mature_pollen_28d_3,
			flower_early_stg12_carpels_21d._1,
			flower_early_stg12_carpels_21d._2,
			flower_early_stg12_carpels_21d._3)) #tibble w/o pollen samles

		colnames(all_genes_tpm)[colnames(all_genes_tpm) == 'flowers_mature_pollen_28d_1'] <- 'flowers_mature_pollen_1'
		colnames(all_genes_tpm)[colnames(all_genes_tpm) == 'flowers_mature_pollen_28d_2'] <- 'flowers_mature_pollen_2'
		colnames(all_genes_tpm)[colnames(all_genes_tpm) == 'flowers_mature_pollen_28d_3'] <- 'flowers_mature_pollen_3'

		species_id <- "ATH_comparative_samples"


    } else if ((is.element("AL", species)) && (is.element("comparative", experiment))) {

		all_genes_tpm <- dplyr::select(all_genes_tpm, -c(
			flower_stg11_stamens_8w.10w.25d_1, 
			flower_stg11_stamens_8w.10w.25d_2, 
			flower_stg11_stamens_8w.10w.25d_3,
			flower_early_stg12_stamens_8w.10w.23d_1,
			flower_early_stg12_stamens_8w.10w.23d_2,
			flower_early_stg12_stamens_8w.10w.23d_3,
			flower_late_stg12_stamens_8w.10w.21d_1,
			flower_late_stg12_stamens_8w.10w.21d_2,
			flower_late_stg12_stamens_8w.10w.21d_3)) #tibble w/o pollen samles

		colnames(all_genes_tpm)[colnames(all_genes_tpm) == 'flowers_mature_pollen_8w.10w.25d_1'] <- 'flowers_mature_pollen_1'
		colnames(all_genes_tpm)[colnames(all_genes_tpm) == 'flowers_mature_pollen_8w.10w.25d_2'] <- 'flowers_mature_pollen_2'
		colnames(all_genes_tpm)[colnames(all_genes_tpm) == 'flowers_mature_pollen_8w.10w.25d_3'] <- 'flowers_mature_pollen_3'

		species_id <- "AL_comparative_samples"


    } else if (is.element("ATH", species)) {
		colnames(all_genes_tpm)[colnames(all_genes_tpm) == 'flowers_mature_pollen_28d_1'] <- 'flowers_mature_pollen_1'
		colnames(all_genes_tpm)[colnames(all_genes_tpm) == 'flowers_mature_pollen_28d_2'] <- 'flowers_mature_pollen_2'
		colnames(all_genes_tpm)[colnames(all_genes_tpm) == 'flowers_mature_pollen_28d_3'] <- 'flowers_mature_pollen_3'

    } else if (is.element("AL", species)) {
    	colnames(all_genes_tpm)[colnames(all_genes_tpm) == 'flowers_mature_pollen_8w.10w.25d_1'] <- 'flowers_mature_pollen_1'
		colnames(all_genes_tpm)[colnames(all_genes_tpm) == 'flowers_mature_pollen_8w.10w.25d_2'] <- 'flowers_mature_pollen_2'
		colnames(all_genes_tpm)[colnames(all_genes_tpm) == 'flowers_mature_pollen_8w.10w.25d_3'] <- 'flowers_mature_pollen_3'

    } else if (is.element("CR", species)) {
    	colnames(all_genes_tpm)[colnames(all_genes_tpm) == 'flowers_mature_pollen_6w.7w.22d_1'] <- 'flowers_mature_pollen_1'
		colnames(all_genes_tpm)[colnames(all_genes_tpm) == 'flowers_mature_pollen_6w.7w.22d_2'] <- 'flowers_mature_pollen_2'
		colnames(all_genes_tpm)[colnames(all_genes_tpm) == 'flowers_mature_pollen_6w.7w.22d_3'] <- 'flowers_mature_pollen_3'

    } else if (is.element("ES", species)) {
		colnames(all_genes_tpm)[colnames(all_genes_tpm) == 'flowers_mature_pollen_7w.8w.17d_1'] <- 'flowers_mature_pollen_1'
		colnames(all_genes_tpm)[colnames(all_genes_tpm) == 'flowers_mature_pollen_7w.8w.17d_2'] <- 'flowers_mature_pollen_2'
		colnames(all_genes_tpm)[colnames(all_genes_tpm) == 'flowers_mature_pollen_7w.8w.17d_3'] <- 'flowers_mature_pollen_3'

    } else if (is.element("TH", species)) {
		colnames(all_genes_tpm)[colnames(all_genes_tpm) == 'flowers_mature_pollen_11w_1'] <- 'flowers_mature_pollen_1'
		colnames(all_genes_tpm)[colnames(all_genes_tpm) == 'flowers_mature_pollen_11w_2'] <- 'flowers_mature_pollen_2'
		colnames(all_genes_tpm)[colnames(all_genes_tpm) == 'flowers_mature_pollen_11w_3'] <- 'flowers_mature_pollen_3'

    } else if (is.element("MT", species)) {
		colnames(all_genes_tpm)[colnames(all_genes_tpm) == 'flowers_mature_pollen_7w_1'] <- 'flowers_mature_pollen_1'
		colnames(all_genes_tpm)[colnames(all_genes_tpm) == 'flowers_mature_pollen_7w_2'] <- 'flowers_mature_pollen_2'
		colnames(all_genes_tpm)[colnames(all_genes_tpm) == 'flowers_mature_pollen_7w_3'] <- 'flowers_mature_pollen_3'

    } else if (is.element("BD", species)) {
		colnames(all_genes_tpm)[colnames(all_genes_tpm) == 'flowers_mature_pollen_32d_1'] <- 'flowers_mature_pollen_1'
		colnames(all_genes_tpm)[colnames(all_genes_tpm) == 'flowers_mature_pollen_32d_2'] <- 'flowers_mature_pollen_2'
		colnames(all_genes_tpm)[colnames(all_genes_tpm) == 'flowers_mature_pollen_32d_3'] <- 'flowers_mature_pollen_3'
    }


    
    #--------- Extract protein-coding genes, seperate them by strand and find overlap ----------


	GTF_df = as.data.frame(GTF)

	# Get all protein-coding genes
	GTF_df_cd <- subset(GTF_df, gene_biotype == "protein_coding")

	# Separate genes from plus strand and minus strand
	strand_plus <- subset(GTF_df_cd, strand == "+")
	strand_minus <- subset(GTF_df_cd, strand == "-")


	strand_plus_granges <- makeGRangesFromDataFrame(strand_plus, keep.extra.columns=FALSE, 
		ignore.strand=FALSE, seqinfo=NULL, seqnames.field=c(
		"seqnames", "seqname","chromosome", "chrom","chr", "chromosome_name","seqid"),
		start.field="start", end.field=c("end", "stop"),starts.in.df.are.0based=FALSE)

	strand_minus_granges <- makeGRangesFromDataFrame(strand_minus, keep.extra.columns=FALSE, 
		ignore.strand=FALSE, seqinfo=NULL, seqnames.field=c(
		"seqnames", "seqname","chromosome", "chrom","chr", "chromosome_name","seqid"),
		start.field="start", end.field=c("end", "stop"),starts.in.df.are.0based=FALSE)


	# Find overlaps between genes from plus and minus strands
	overlap_with_strand = findOverlaps(strand_plus_granges, strand_minus_granges, ignore.strand = TRUE)
	overlap_with_strand_df = as.data.frame(overlap_with_strand)
	names(overlap_with_strand_df) = c("key_plus", "key_minus")
	dim(overlap_with_strand) #number of protein-coding genes that overlapp




	#--- Extract overlapping protein coding genes from plus_strand/minus_strand data frames ----


	# Add key to plus and minus strand data frames
	key_plus <- seq(1, nrow(strand_plus), 1)
	strand_plus <- cbind(as.data.frame(key_plus), strand_plus)

	key_minus <- seq(1, nrow(strand_minus), 1)
	strand_minus <- cbind(as.data.frame(key_minus), strand_minus)


	# Generate thresholded data frame based on keys
	strand_plus_overlapp_genes <- merge(overlap_with_strand_df, strand_plus, by="key_plus")
	strand_plus_overlapp_genes$id  <- seq(1, nrow(strand_plus_overlapp_genes), 1)
	strand_plus_overlapp_genes = strand_plus_overlapp_genes %>% select(
		id,
		gene_id,
		seqnames, 
		start, 
		end, 
		strand)

	strand_minus_overlapp_genes <- merge(overlap_with_strand_df, strand_minus, by="key_minus")
	strand_minus_overlapp_genes$id  <- seq(1, nrow(strand_minus_overlapp_genes), 1)
	strand_minus_overlapp_genes = strand_minus_overlapp_genes %>% select(
		id,
		gene_id,
		seqnames, 
		start, 
		end, 
		strand)


	#---------- Merge plus_strand/minus_strand data tables with DevSeq expression data ---------

	strand_plus_overlapp_genes_tpm <- merge(strand_plus_overlapp_genes, all_genes_tpm, by="gene_id")
	strand_plus_overlapp_genes_tpm <- strand_plus_overlapp_genes_tpm[order(strand_plus_overlapp_genes_tpm$id),]
	strand_plus_overlapp_genes_tpm <- dplyr::select(strand_plus_overlapp_genes_tpm, -c(
		seqnames, source))
	strand_minus_overlapp_genes_tpm <- merge(strand_minus_overlapp_genes, all_genes_tpm, by="gene_id")
	strand_minus_overlapp_genes_tpm <- strand_minus_overlapp_genes_tpm[order(strand_minus_overlapp_genes_tpm$id),]
	strand_minus_overlapp_genes_tpm <- dplyr::select(strand_minus_overlapp_genes_tpm, -c(
		seqnames, source))



	#---------------------------- Apply TPM expression threshold -------------------------------


	# This is a threshold function that can be applied to expression tables
	# Settings: TPM > 0.5 in at least 2 of 3 replicates
	applyThreshold <- function(df, threshold) {
  
  		#* Add an error if radius < 0
  		if (threshold < 0)
    		stop(
           	"'threshold' must be >= 0",
	   		call. = TRUE
    		)

		# Add keys to data frame
		key <- seq(1, nrow(df), 1)
		df <- cbind(as.data.frame(key),df)

		# Define threshold function
		getThreshold <- function(df) {

			# Split data frame by sample replicates into a list then apply threshold for each subset
	
			th_replicates <- do.call(cbind, lapply(split.default(df[8:ncol(df)], #adjust columns
								rep(seq_along(df), each = 3, length.out = ncol(df)-7)), #adjust columns
								function(x) {
									x[rowSums(x > threshold) < 2, ] <- 0; 
									x
								}
							))

			# Bind key/id/prt_id/symbol/biotype/source columns to thresholded data frame
			th_replicates <- cbind(df[1:7], th_replicates)

			# Remove all rows that only contain "0"
			th_replicates <- th_replicates[which(rowSums(th_replicates[,-1:-7, drop = FALSE] > 0) > 0),]

			return(th_replicates)
		}

		# Apply threshold to data and extract keys ("key")
		keys_data <- getThreshold(df)
		keys_data <- keys_data[,1:2]
		names(keys_data) <- c("key","ID")

		# Generate thresholded data frame based on keys
		th_df <- merge(keys_data, df, by="key")
		th_df <- th_df[order(th_df$key),]
		th_df <- th_df[-1:-2]

		return(th_df)
	}


	# Apply threshold function
	strand_plus_overlapp_genes_tpm_0.5 <- applyThreshold(strand_plus_overlapp_genes_tpm,0.5)
	strand_minus_overlapp_genes_tpm_0.5 <- applyThreshold(strand_minus_overlapp_genes_tpm,0.5)

	strand_plus_overlapp_genes_tpm_0.5 <- strand_plus_overlapp_genes_tpm_0.5[(
		strand_plus_overlapp_genes_tpm_0.5$id %in% strand_minus_overlapp_genes_tpm_0.5$id),]
	strand_minus_overlapp_genes_tpm_0.5 <- strand_minus_overlapp_genes_tpm_0.5[(
		strand_minus_overlapp_genes_tpm_0.5$id %in% strand_plus_overlapp_genes_tpm_0.5$id),]




	#-------------------------- Get expression table w/o pollen data  --------------------------


	# Remove pollen triplicates from expression table
	strand_plus_overlapp_genes_tpm_0.5_wo_pollen <- dplyr::select(strand_plus_overlapp_genes_tpm_0.5, -c(
		flowers_mature_pollen_1, 
		flowers_mature_pollen_2, 
		flowers_mature_pollen_3)) #tibble w/o pollen samles

	strand_minus_overlapp_genes_tpm_0.5_wo_pollen <- dplyr::select(strand_minus_overlapp_genes_tpm_0.5, -c(
		flowers_mature_pollen_1, 
		flowers_mature_pollen_2, 
		flowers_mature_pollen_3)) #tibble w/o pollen samles



	#------------------- Compute correlation between coding gene and cisNAT  -------------------


	# Compute spearman and pearson correlations
	getCor <- function(df1, df2) {

		df1_col <- ncol(df1)
		df2_col <- ncol(df2)

		# startup message
		message("Computing correlation...")

		df1$Spearman <- sapply(1:nrow(df1), function(i) 
	    	cor(as.numeric(df1[i, 7:df1_col]), as.numeric(df2[i, 7:df2_col]), method=c("spearman")))

		df1$Pearson <- sapply(1:nrow(df1), function(i) 
	    	cor(as.numeric(df1[i, 7:df1_col]), as.numeric(df2[i, 7:df2_col]), method=c("pearson")))

		return(df1)
	}


	strand_plus_overlapp_genes <- getCor(
		strand_plus_overlapp_genes_tpm_0.5, strand_minus_overlapp_genes_tpm_0.5)
	strand_plus_overlapp_genes_wo_pollen <- getCor(
		strand_plus_overlapp_genes_tpm_0.5_wo_pollen, strand_minus_overlapp_genes_tpm_0.5_wo_pollen)




#------------------------- Create final data table and write csv file  -------------------------


	# Create data table containing both strand plus and minus genes and cor values
	strand_minus_overlapp_genes_descript = strand_minus_overlapp_genes_tpm_0.5 %>% select(gene_id, start, end, strand, biotype)
	names(strand_minus_overlapp_genes_descript) <- c("id_minus_strand", "start_minus", "end_minus", "strand_subject", "biotype_subject")
	strand_plus_overlapp_genes <- strand_plus_overlapp_genes[-2]
	names(strand_plus_overlapp_genes)[1:5] <- c("id_plus_strand", "start_plus", "end_plus", "strand_query", "biotype_query")
	overlapp_cd_genes_cor <- cbind(strand_plus_overlapp_genes, strand_minus_overlapp_genes_descript)
	overlapp_cd_genes_cor = overlapp_cd_genes_cor %>% select(
		id_plus_strand, 
		start_plus, 
		end_plus, 
		strand_query, 
		biotype_query, 
		id_minus_strand, 
		start_minus, 
		end_minus, 
		strand_subject, 
		biotype_subject, 
		Spearman, 
		Pearson)


	strand_minus_overlapp_genes_descript_wo_pollen = strand_minus_overlapp_genes_tpm_0.5_wo_pollen %>% select(gene_id, start, end, strand, biotype)
	names(strand_minus_overlapp_genes_descript_wo_pollen) <- c("id_minus_strand", "start_minus", "end_minus", "strand_subject", "biotype_subject")
	strand_plus_overlapp_genes_wo_pollen <- strand_plus_overlapp_genes_wo_pollen[-2]
	names(strand_plus_overlapp_genes_wo_pollen)[1:5] <- c("id_plus_strand", "start_plus", "end_plus", "strand_query", "biotype_query")
	overlapp_cd_genes_cor_wo_pollen <- cbind(strand_plus_overlapp_genes_wo_pollen, 
		strand_minus_overlapp_genes_descript_wo_pollen)
	overlapp_cd_genes_cor_wo_pollen = overlapp_cd_genes_cor_wo_pollen %>% select(
		id_plus_strand, 
		start_plus, 
		end_plus, 
		strand_query, 
		biotype_query, 
		id_minus_strand, 
		start_minus, 
		end_minus, 
		strand_subject, 
		biotype_subject, 
		Spearman, 
		Pearson)


	# Set filename
    fname <- sprintf('%s.csv', paste(species_id, "coding_SAS_cor", sep="_"))
	fname_wo_pollen <- sprintf('%s.csv', paste(species_id, "coding_SAS_cor_wo_pollen", sep="_"))


	# Write final data tables to csv files and store them in /out_dir/output/data_tables
	if (!dir.exists(file.path(out_dir, "output", "overlapp_cd_genes"))) 
		dir.create(file.path(out_dir, "output", "overlapp_cd_genes"), recursive = TRUE)

	write.table(overlapp_cd_genes_cor, file=file.path(out_dir, "output", "overlapp_cd_genes", fname), 
		sep=";", dec=".", row.names=FALSE, col.names=TRUE)
	write.table(overlapp_cd_genes_cor_wo_pollen, file=file.path(out_dir, "output", "overlapp_cd_genes", fname_wo_pollen), 
		sep=";", dec=".", row.names=FALSE, col.names=TRUE)

}



# Execute getPcPc function
getPcPc("ATH", "single-species")
getPcPc("ATH", "comparative")
getPcPc("AL", "single-species")
getPcPc("AL", "comparative")
getPcPc("CR", "comparative")
getPcPc("ES", "comparative")
getPcPc("TH", "comparative")
getPcPc("MT", "comparative")
getPcPc("BD", "comparative")




