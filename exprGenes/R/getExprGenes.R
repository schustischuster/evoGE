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
in_dir <- "/Volumes/User/Shared/Christoph_manuscript/DevSeq_paper/Analysis/Analysis_2019/A_thaliana_gene_exression_map/20200401_CS_exprGenes/data"
out_dir <- "/Volumes/User/Shared/Christoph_manuscript/DevSeq_paper/Analysis/Analysis_2019/A_thaliana_gene_exression_map/20200401_CS_exprGenes"




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
        genesCounts = file.path(in_dir, "Expression_data", "ATH_no_TE_genes_counts_sample_names.csv")
        species_id <- "ATH"

    } else if (is.element("AL", species)) {
		genesTPM = file.path(in_dir, "Expression_data", "AL_no_TE_genes_tpm_sample_names.csv")
		genesCounts = file.path(in_dir, "Expression_data", "AL_no_TE_genes_counts_sample_names.csv")
		species_id <- "AL"

    } else if (is.element("CR", species)) {
		genesTPM = file.path(in_dir, "Expression_data", "CR_no_TE_genes_tpm_sample_names.csv")
		genesCounts = file.path(in_dir, "Expression_data", "CR_no_TE_genes_counts_sample_names.csv")
		species_id <- "CR"

    } else if (is.element("ES", species)) {
		genesTPM = file.path(in_dir, "Expression_data", "ES_no_TE_genes_tpm_sample_names.csv")
		genesCounts = file.path(in_dir, "Expression_data", "ES_no_TE_genes_counts_sample_names.csv")
		species_id <- "ES"

    } else if (is.element("TH", species)) {
		genesTPM = file.path(in_dir, "Expression_data", "TH_no_TE_genes_tpm_sample_names.csv")
		genesCounts = file.path(in_dir, "Expression_data", "TH_no_TE_genes_counts_sample_names.csv")
		species_id <- "TH"

    } else if (is.element("MT", species)) {
		genesTPM = file.path(in_dir, "Expression_data", "MT_no_TE_genes_tpm_sample_names.csv")
		genesCounts = file.path(in_dir, "Expression_data", "MT_no_TE_genes_counts_sample_names.csv")
		species_id <- "MT"

    } else if (is.element("BD", species)) {
		genesTPM = file.path(in_dir, "Expression_data", "BD_no_TE_genes_tpm_sample_names.csv")
		genesCounts = file.path(in_dir, "Expression_data", "BD_no_TE_genes_counts_sample_names.csv")
		species_id <- "BD"
    }


    # Get experiment
    if (is.element("single-species", experiment)) {
    	experiment <- "single-species"
    } else experiment <- "comparative"


	# Read expression data
	all_genes_tpm <- read.table(genesTPM, sep=";", dec=".", header=TRUE, stringsAsFactors=FALSE)
	colnames(all_genes_tpm)[1] <- "gene_id"
	all_genes_counts <- read.table(genesCounts, sep=";", dec=".", header=TRUE, stringsAsFactors=FALSE)
	colnames(all_genes_counts)[1] <- "gene_id"


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

		all_genes_counts <- dplyr::select(all_genes_counts, -c(
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
    # return_list <- list("species_id" = species_id, "experiment" = experiment, "all_genes_tpm" = all_genes_tpm, "all_genes_counts" = all_genes_counts, "threshold" = threshold)
    # return(return_list)
    # }
    # return_objects <- getExprGenes(species="ATH", experiment="single-species", 0.05) # read in GTF and expression data for A.thaliana
    # list2env(return_objects, envir = .GlobalEnv)



#--------------------------- Get and apply sample-specific threshold  --------------------------


    # Show message
    message("Starting analysis...")


	# Extract ERCC data
	ERCC <- all_genes_tpm[all_genes_tpm$gene_id %like% "ERCC", ]


	if (threshold > 0) {

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


	if (threshold == 0) {
		th_values <- data.frame(t(rep(0, ncol(ERCC)-3)))
		names(th_values) <- colnames(all_genes_tpm)[4:ncol(all_genes_tpm)]
		ERCC_cutoff <- data.frame(gene_id="TPM_cutoff", biotype="NA", source="NA", th_values)
		all_genes_tpm_cutoff <- rbind(all_genes_tpm, ERCC_cutoff)
	}




	#* test data
	#* df <- data.frame (ID  = c('data1', 'data2', 'data3', 'data4', 'data5', 'data6', 'data7', 'TPM_cutoff'), 
	#* 	biotype  = c('coding', 'coding', 'coding', 'coding', 'coding', 'coding', 'coding', 'NA'), 
	#* 	source  = c('DevSeq', 'DevSeq', 'DevSeq', 'DevSeq', 'DevSeq', 'Araport', 'DevSeq', 'NA'), 
	#*     sample1 = c(2, 7, 1, 18, 3, 0.1, 0, 3),
	#*     sample2 = c(4, 0, 3, 17, 16, 0.3, 0, 5),
	#*     sample3 = c(3, 4, 2, 11, 2, 0.2, 0, 2),
	#*     sample4 = c(22, 14, 9, 11, 35, 0, 0, 11),
	#*     sample5 = c(10, 8, 5, 8, 22, 0, 0, 7),
	#*     sample6 = c(17, 11, 6, 9, 11, 0.1, 0, 8))




	# This is a threshold function that can be applied to expression tables
	# Settings: TPM > ERCC spike-in threshold in at least 2 of 3 replicates in at least one sample
	applyThreshold <- function(express_df) {

		# Add keys to data frame
		key <- seq(1, nrow(express_df), 1)
		express_df <- cbind(as.data.frame(key),express_df)

		# Replace all values with "0" that are below sample threshold (either 0 or ERCC)
		getSampleTH <- function(df) {

			# Split data frame by sample replicates into a list then apply threshold for each subset
	
			th_replicates <- do.call(cbind, lapply(split.default(df[5:ncol(df)], #adjust columns
								rep(seq_along(df), each = 1, length.out = ncol(df)-4)),
								function(x) {
									x[x <= x[nrow(df),], ] <- 0;
									x
								}
							))

			# Bind key/id/prt_id/symbol/biotype/source columns to thresholded data frame
			th_replicates <- cbind(df[1:4], th_replicates)

			return(th_replicates)
		}

		if (threshold > 0) {
			df <- getSampleTH(express_df)
		} else 
			df <- express_df

		# Define threshold function
		# This function will remove all rows that do not show expression in at least two of three
		# replicates in at least one sample type
		getThreshold <- function(df) {

			# Split data frame by sample replicates into a list then apply threshold for each subset
	
			th_replicates <- do.call(cbind, lapply(split.default(df[5:ncol(df)], #adjust columns
								rep(seq_along(df), each = 3, length.out = ncol(df)-4)), #adjust columns
								function(x) {
									x[rowSums(x > 0) < 2, ] <- 0; 
									x
								}
							))

			# Bind key/id/prt_id/symbol/biotype/source columns to thresholded data frame
			th_replicates <- cbind(df[1:4], th_replicates)

			# Remove all rows that only contain "0"
			th_replicates <- th_replicates[which(rowSums(th_replicates[,-1:-4, drop = FALSE] > 0) > 0),]

			return(th_replicates)
		}

		# Apply threshold to data and extract keys ("key")
		keys_data_repl <- getThreshold(df)
		keys_data <- keys_data_repl[,1:2]
		names(keys_data) <- c("key","ID")

		# Generate thresholded data frame based on keys
		th_df <- merge(keys_data, express_df, by="key")
		th_df <- th_df[order(th_df$key),]
		th_df <- th_df[-1:-2]
		keys_data_repl <- keys_data_repl[-1]

		# Create list of return objects
		return_th_list <- list("express_data_th"=th_df, "express_data_th_replicates"=keys_data_repl)
		return(return_th_list)
	}

    return_objects_th <- applyThreshold(all_genes_tpm_cutoff)
    list2env(return_objects_th, envir = .GlobalEnv)
    #* there are two return objects: 
    #* (1) express_data_th = expression data that only contains genes that show in at least one 
       #* sample type in at least two out of three replicates an expression above the individual 
       #* ERCC spike-in threshold (at defined level of either 0.01/0.05/0.1) or above 0 (if
       #* threshold=0)
    #* (2) express_data_th_replicates = same file as above, but if expression value of replicate
       #* is below the individual ERCC spike-in threshold (or 0), the expression values of the
       #* replicates will be replaced by "0", and if two out of three replicates of one sample
       #* type get "0", the third replicate value also gets "0" (therefore sum of all three
       #* replicates will be zero, which defines a gene that is not expressed [above threshold
       #* level])




#----------- Merge replicates and retrieve number of expressed genes for each sample -----------


    calculateAvgExpr <- function(df) {

	# Split data frame by sample replicates into a list
	# then get rowMeans for each subset and bind averaged data to gene_id/biotype/source column
	
	averaged_replicates <- do.call(cbind, lapply(split.default(df[4:ncol(df)], 
			rep(seq_along(df), 
			each = 3, 
			length.out=ncol(df)-3)
			), rowMeans)
		)

		averaged_replicates <- cbind(df[1:3], averaged_replicates)
		
		return(averaged_replicates)
	}


	express_data_th_avg <- calculateAvgExpr(express_data_th_replicates)

	repl_names <- unique(substr(names(express_data_th_replicates[-1:-3]), 1, nchar(names(express_data_th_replicates[-1:-3]))-2))
	repl_names <- substring(repl_names, 3)

	if ((is.element("ATH", species_id)) && (is.element("single-species", experiment))) { 

		repl_names[-1:-9] <- substring(repl_names[-1:-9], 2)
		repl_names[-1:-23] <- substr(repl_names[-1:-23],1,nchar(repl_names[-1:-23])-1)
		repl_names[10] <- substr(repl_names[10],1,nchar(repl_names[10])-1)
		repl_names <- gsub("apex_inflorescence_28", "apex_inflorescence_28d", repl_names)
		repl_names <- gsub("flowers_mature_polle", "flowers_mature_pollen", repl_names)
	}

	colnames(express_data_th_avg) <- c("gene_id", "biotype", "source", repl_names)


	# Get number of genes expressed in each sample type
	protein_coding_subset <- subset(express_data_th_avg, biotype=="protein_coding")
	lnc_intergenic_subset <- subset(express_data_th_avg, biotype=="lnc_intergenic")
	lnc_antisense_subset <- express_data_th_avg[express_data_th_avg$biotype %like% "antisense", ]

	
	# Remove chloroplast and mito genes from coding gene list
	`%nlike%` = Negate(`%like%`)

	if ((is.element("ATH", species_id))) { 

		protein_coding_subset <- protein_coding_subset[protein_coding_subset$gene_id %nlike% "ATMG", ]
		protein_coding_subset <- protein_coding_subset[protein_coding_subset$gene_id %nlike% "ATCG", ]
	}


	expr_protein_coding <- colSums(protein_coding_subset > 0)
	expr_lnc_antisense <- colSums(lnc_antisense_subset > 0)
	expr_lnc_intergenic <- colSums(lnc_intergenic_subset > 0)


	# Create final list of expressed genes per organ/ sample type
	expressed_genes_th_avg <- rbind(expr_protein_coding, expr_lnc_antisense, expr_lnc_intergenic)
	expressed_genes_th_avg <- subset(expressed_genes_th_avg, select = -c(biotype, source))
	colnames(expressed_genes_th_avg)[1] <- "total_expressed"
	biotype <- c("protein_coding", "lnc_antisense", "lnc_intergenic")
	expressed_genes_th_avg <- as.data.frame(cbind(biotype, expressed_genes_th_avg))




#----------- Calculate spearman correlation for biological replicates of each sample -----------


	# Get rlog count thresholded expression table based on TPM cutoff
	express_data_th_counts <- subset(all_genes_counts, gene_id %in% express_data_th$gene_id)


	# Get replicate correlation
	replCorr <- function(df, coefficient=c("spearman", "pearson"), expr_est=c("tpm", "counts")) {

		if (is.element("pearson", coefficient) && is.element("tpm", expr_est)) {
			df[,4:ncol(df)] <- log2(df[,4:ncol(df)] + 1)
		}

		replicate_corr <- do.call(cbind, lapply(split.default(df[4:ncol(df)], #adjust columns
								rep(seq_along(df), each = 3, length.out = ncol(df)-3)), #adjust columns
									function(x) {
									repl_corr <- data.frame(cbind(
										cor(x[,1],x[,2],method=coefficient),
										cor(x[,1],x[,3],method=coefficient),
										cor(x[,2],x[,3],method=coefficient)));
									repl_corr
								}
							))

		return(replicate_corr)
	}

	
	# Get replicate correlations based on TPM values
	repl_corr_df_tpm <- replCorr(express_data_th, coefficient="pearson", expr_est="tpm")
	colnames(repl_corr_df_tpm) <- rep(repl_names,each=3)

	# Get replicate correlations based on normalized rlog counts
	repl_corr_df_counts <- replCorr(express_data_th_counts, coefficient="pearson", expr_est="counts")
	colnames(repl_corr_df_counts) <- rep(repl_names,each=3)




#-------- Generate replicate expression lists with thresholded genes for each biotype ---------


	# For ERCC-thresholded rlog count data
	protein_coding_th_repl_counts <- subset(all_genes_counts, gene_id %in% protein_coding_subset$gene_id)
	lnc_antisense_th_repl_counts <- subset(all_genes_counts, gene_id %in% lnc_antisense_subset$gene_id)
	lnc_intergenic_th_repl_counts <- subset(all_genes_counts, gene_id %in% lnc_intergenic_subset$gene_id)
	th_genes_counts <- rbind(protein_coding_th_repl_counts, lnc_antisense_th_repl_counts, lnc_intergenic_th_repl_counts)

	# For ERCC-thresholded tpm data - first log2 transform data
	all_genes_tpm_log <- all_genes_tpm
	all_genes_tpm_log[,4:ncol(all_genes_tpm_log)] <- log2(all_genes_tpm_log[,4:ncol(all_genes_tpm_log)] + 1)
	protein_coding_th_repl_tpm <- subset(all_genes_tpm_log, gene_id %in% protein_coding_subset$gene_id)
	lnc_antisense_th_repl_tpm <- subset(all_genes_tpm_log, gene_id %in% lnc_antisense_subset$gene_id)
	lnc_intergenic_th_repl_tpm <- subset(all_genes_tpm_log, gene_id %in% lnc_intergenic_subset$gene_id)
	th_genes_tpm <- rbind(protein_coding_th_repl_tpm, lnc_antisense_th_repl_tpm, lnc_intergenic_th_repl_tpm)




#--------------------------------------- Write csv file  ---------------------------------------


	# Show message
    message("Writing output...")


	# Set filename
    fname_expr_genes <- sprintf('%s.csv', paste(species_id, "expr_genes", threshold, sep="_"))
    fname_repl_corr_tpm <- sprintf('%s.csv', paste(species_id, "repl_corr_tpm", threshold, sep="_"))
    fname_repl_corr_counts <- sprintf('%s.csv', paste(species_id, "repl_corr_counts", threshold, sep="_"))
    fname_th_genes_repl_tpm <- sprintf('%s.csv', paste(species_id, "th_genes_repl_tpm", threshold, sep="_"))
    fname_th_genes_repl_counts <- sprintf('%s.csv', paste(species_id, "th_genes_repl_counts", threshold, sep="_"))


	# Write final data tables to csv files and store them in /out_dir/output/data_tables
	if (!dir.exists(file.path(out_dir, "output", "expr_genes"))) 
		dir.create(file.path(out_dir, "output", "expr_genes"), recursive = TRUE)

	write.table(expressed_genes_th_avg, file=file.path(out_dir, "output", "expr_genes", fname_expr_genes), 
		sep=";", dec=".", row.names=FALSE, col.names=TRUE)
	write.table(repl_corr_df_tpm, file=file.path(out_dir, "output", "expr_genes", fname_repl_corr_tpm), 
		sep=";", dec=".", row.names=FALSE, col.names=TRUE)
	write.table(repl_corr_df_counts, file=file.path(out_dir, "output", "expr_genes", fname_repl_corr_counts), 
		sep=";", dec=".", row.names=FALSE, col.names=TRUE)

	# Only write thresholded expression tables for 0.05 ERCC threshold to file
	if (threshold == 0.05) {
		write.table(th_genes_tpm, file=file.path(out_dir, "output", "expr_genes", fname_th_genes_repl_tpm), 
			sep=";", dec=".", row.names=FALSE, col.names=TRUE)
		write.table(th_genes_counts, file=file.path(out_dir, "output", "expr_genes", fname_th_genes_repl_counts), 
			sep=";", dec=".", row.names=FALSE, col.names=TRUE)
	}

}



thresholds <- list(0, 0.01, 0.05, 0.1) # threshold values are 0 or perc of expressed spike-ins

lapply(thresholds, getExprGenes, species = "ATH", experiment = "single-species")
lapply(thresholds, getExprGenes, species = "AL", experiment = "comparative")
lapply(thresholds, getExprGenes, species = "CR", experiment = "comparative")
lapply(thresholds, getExprGenes, species = "ES", experiment = "comparative")
lapply(thresholds, getExprGenes, species = "TH", experiment = "comparative")
lapply(thresholds, getExprGenes, species = "MT", experiment = "comparative")
lapply(thresholds, getExprGenes, species = "BD", experiment = "comparative")


