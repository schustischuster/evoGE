# Find expressed transcripts in each sample replicate using ERCC spike-in information
# Thresholds: 0 (static 0.5 TPM threshold) OR
# 0.01, 0.05, 0.1 ERCC (percentile of detected spike-ins)
# Data input: Expression_data WITH spike-in information
# Analysis can be performed on both whole single species datasets (ATH: 132 samples; AL: 36 samples)
# OR on comparative data sets (each species: 27 samples)
# Input sample tables should have the following format:
# id / biotype / source / info / DEVSEQ_SAMPLE_REPLICATES(between 27 and 132 depending on species)




#------------------- Load packages, set directories and read sample tables ---------------------


# Define function to get expressed transcripts

getExprTranscripts <- function(species = c("ATH", "AL", "CR", "ES", "TH", "MT", "BD"), 
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
        transcriptsTPM = file.path(in_dir, "Expression_data", "AT_transcripts_complete_table_tpm_sample_names.csv")
        transcriptsPtMt = file.path(in_dir, "Expression_data", "AT_Pt_Mt_Orthologs.csv")
        species_id <- "ATH"

    } else if (is.element("AL", species)) {
		transcriptsTPM = file.path(in_dir, "Expression_data", "AL_transcripts_complete_table_tpm_sample_names.csv")
		transcriptsPtMt = file.path(in_dir, "Expression_data", "AL_Pt_Mt_Orthologs.csv")
		species_id <- "AL"

    } else if (is.element("CR", species)) {
		transcriptsTPM = file.path(in_dir, "Expression_data", "CR_transcripts_complete_table_tpm_sample_names.csv")
		transcriptsPtMt = file.path(in_dir, "Expression_data", "CR_Pt_Mt_Orthologs.csv")
		species_id <- "CR"

    } else if (is.element("ES", species)) {
		transcriptsTPM = file.path(in_dir, "Expression_data", "ES_transcripts_complete_table_tpm_sample_names.csv")
		transcriptsPtMt = file.path(in_dir, "Expression_data", "ES_Pt_Mt_Orthologs.csv")
		species_id <- "ES"

    } else if (is.element("TH", species)) {
		transcriptsTPM = file.path(in_dir, "Expression_data", "TH_transcripts_complete_table_tpm_sample_names.csv")
		transcriptsPtMt = file.path(in_dir, "Expression_data", "TH_Pt_Mt_Orthologs.csv")
		species_id <- "TH"

    } else if (is.element("MT", species)) {
		transcriptsTPM = file.path(in_dir, "Expression_data", "MT_transcripts_complete_table_tpm_sample_names.csv")
		transcriptsPtMt = file.path(in_dir, "Expression_data", "MT_Pt_Mt_Orthologs.csv")
		species_id <- "MT"

    } else if (is.element("BD", species)) {
		transcriptsTPM = file.path(in_dir, "Expression_data", "BD_transcripts_complete_table_tpm_sample_names.csv")
		transcriptsPtMt = file.path(in_dir, "Expression_data", "BD_Pt_Mt_Orthologs.csv")
		species_id <- "BD"
    }


    # Get experiment
    if (is.element("single-species", experiment)) {
    	experiment <- "single-species"
    } else experiment <- "comparative"


	# Read expression data
	all_transcripts_tpm <- read.table(transcriptsTPM, sep=";", dec=".", header=TRUE, stringsAsFactors=FALSE)
	colnames(all_transcripts_tpm)[1] <- "transcript_id"
	transcripts_ptmt <- read.table(transcriptsPtMt, sep=";", dec=".", header=TRUE, stringsAsFactors=FALSE)


	# Format expression data and rename pollen samples
    if ((is.element("ATH", species)) && (is.element("comparative", experiment))) {

		all_transcripts_tpm <- dplyr::select(all_transcripts_tpm, c(
			"transcript_id", "biotype", "source", "info", 
			"root_whole_root_5d_1",
			"root_whole_root_5d_2",
			"root_whole_root_5d_3",
			"hypocotyl_10d_1",
			"hypocotyl_10d_2",
			"hypocotyl_10d_3",
			"leaf_1+2_7d_1",
			"leaf_1+2_7d_2",
			"leaf_1+2_7d_3",
			"apex_vegetative_7d_1",
			"apex_vegetative_7d_2",
			"apex_vegetative_7d_3",
			"apex_inflorescence_21d_1",
			"apex_inflorescence_21d_2",
			"apex_inflorescence_21d_3",
			"flower_stg12_21d+_1",
			"flower_stg12_21d+_2",
			"flower_stg12_21d+_3",
			"flower_stg12_stamens_21d+_1",
			"flower_stg12_stamens_21d+_2",
			"flower_stg12_stamens_21d+_3",
			"flowers_mature_pollen_28d_1",
			"flowers_mature_pollen_28d_2",
			"flowers_mature_pollen_28d_3",
			"flower_early_stg12_carpels_21d+_1",
			"flower_early_stg12_carpels_21d+_2",
			"flower_early_stg12_carpels_21d+_3")) #tibble w/o pollen samples

		species_id <- "ATH_comparative_samples"


    } else if ((is.element("AL", species)) && (is.element("comparative", experiment))) {

		all_transcripts_tpm <- all_transcripts_tpm[, -which(names(all_transcripts_tpm) %in% c(
			"flower_stg11_stamens_8w.10w.25d_1", 
			"flower_stg11_stamens_8w.10w.25d_2", 
			"flower_stg11_stamens_8w.10w.25d_3",
			"flower_early_stg12_stamens_8w.10w.23d_1",
			"flower_early_stg12_stamens_8w.10w.23d_2",
			"flower_early_stg12_stamens_8w.10w.23d_3",
			"flower_late_stg12_stamens_8w.10w.21d_1",
			"flower_late_stg12_stamens_8w.10w.21d_2",
			"flower_late_stg12_stamens_8w.10w.21d_3"))] #tibble w/o pollen samples

    }



    # Stop function here to allow specific analysis of a single species
    # return_list <- list("species_id" = species_id, "experiment" = experiment, "all_transcripts_tpm" = all_transcripts_tpm, "threshold" = threshold, "transcripts_ptmt" = transcripts_ptmt)
    # return(return_list)
    # }
    # return_objects <- getExprTranscripts(species="ATH", experiment="single-species", 0.05) # read in GTF and expression data for A.thaliana
    # list2env(return_objects, envir = .GlobalEnv)



#--------------------------- Get and apply sample-specific threshold  --------------------------


    # Show message
    message("Starting analysis...")


    # Remove info column
    all_transcripts_tpm <- dplyr::select(all_transcripts_tpm, -c(info))
    

	# Extract ERCC data
	ERCC <- all_transcripts_tpm[all_transcripts_tpm$biotype %like% "ERCC", ]


	if (threshold > 0) {

	   # Get sample-specific TPM threshold
	   ERCC_cutoff <- cbind(
		   data.frame(transcript_id="TPM_cutoff"), data.frame(biotype="<NA>"), data.frame(source="<NA>"), 
		   as.data.frame(t(data.frame(
			TPM_cutoff = apply(ERCC[,4:ncol(ERCC)], 2, function(i)quantile(i[i>0], threshold)))
		   ))
	   )

	   # Bind sample-specific threshold to expression table
	   all_transcripts_tpm_cutoff <- rbind(all_transcripts_tpm, ERCC_cutoff)
	}


	# Set static threshold of 0.5 TPM for "0 ERCC threshold"
	if (threshold == 0) {
		th_values <- data.frame(t(rep(0.5, ncol(ERCC)-3)))
		names(th_values) <- colnames(all_transcripts_tpm)[4:ncol(all_transcripts_tpm)]
		ERCC_cutoff <- data.frame(transcript_id="TPM_cutoff", biotype="NA", source="NA", th_values)
		all_transcripts_tpm_cutoff <- rbind(all_transcripts_tpm, ERCC_cutoff)
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
		express_kdf <- cbind(as.data.frame(key),express_df)

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

		df_exp <- getSampleTH(express_kdf)

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
		keys_data_repl <- getThreshold(df_exp)
		keys_data <- keys_data_repl[,1:2]
		names(keys_data) <- c("key","ID")

		# Generate thresholded data frame based on keys
		th_df <- merge(keys_data, express_kdf, by="key")
		th_df <- th_df[order(th_df$key),]
		th_df <- th_df[-1:-2]
		keys_data_repl <- keys_data_repl[-1]

		# Create list of return objects
		return_th_list <- list("express_data_th"=th_df, "express_data_th_replicates"=keys_data_repl)
		return(return_th_list)
	}

    return_objects_th <- applyThreshold(all_transcripts_tpm_cutoff)
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




#----------- Merge replicates and retrieve number of expressed transcripts for each sample -----------


    calculateAvgExpr <- function(df) {

	# Split data frame by sample replicates into a list
	# then get rowMeans for each subset and bind averaged data to transcript_id/biotype/source column
	
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

	
	# Set replicate colnames
	repl_names <- unique(substr(names(express_data_th_replicates[-1:-3]), 1, nchar(names(express_data_th_replicates[-1:-3]))-2))
	repl_names <- gsub("^[^.]*.", "", repl_names)
	repl_names <- substr(repl_names, 1, nchar(repl_names)-2)
	repl_names <- gsub('^\\.|\\.$', '', repl_names)

	colnames(express_data_th_avg) <- c("transcript_id", "biotype", "source", repl_names)


	# Check which biotypes are in df
	express_data_th_avg$biotype[!(duplicated(express_data_th_avg$biotype))]

	
	# Remove Pt and Mt transcripts
	express_data_th_avg <- express_data_th_avg[!(express_data_th_avg$transcript_id %in% transcripts_ptmt$transcript_id),]


	# Get number of transcripts expressed in each sample type
	protein_coding_subset <- subset(express_data_th_avg, biotype=="protein_coding")

	
	# Double-check that no chloroplast and mito transcripts from coding transcript list
	`%nlike%` = Negate(`%like%`)

	if ((is.element("ATH", species_id))) { 

		protein_coding_subset <- protein_coding_subset[protein_coding_subset$transcript_id %nlike% "ATMG", ]
		protein_coding_subset <- protein_coding_subset[protein_coding_subset$transcript_id %nlike% "ATCG", ]
	}


	expr_protein_coding <- colSums(protein_coding_subset > 0)


	# Create final list of expressed genes per organ/ sample type
	expr_protein_coding_wo_bt <- expr_protein_coding[names(expr_protein_coding) %nlike% "biotype"]
	expressed_transcripts_th_avg <- expr_protein_coding_wo_bt[names(expr_protein_coding_wo_bt) %nlike% "source"]
	names(expressed_transcripts_th_avg)[1] <- "total_expressed"
	biotype_df <- data.frame(biotype = "coding_transcripts")
	expressed_transcripts_th_avg <- as.data.frame(cbind(biotype_df, as.data.frame(t(expressed_transcripts_th_avg))))




#--------------------------------------- Write csv file  ---------------------------------------


	# Show message
    message("Writing output...")


	# Set filename
    fname_expr_transcripts <- sprintf('%s.csv', paste(species_id, "expr_coding_transcripts", threshold, sep="_"))


	# Write final data tables to csv files and store them in /out_dir/output/data_tables
	if (!dir.exists(file.path(out_dir, "output", "expr_genes"))) 
		dir.create(file.path(out_dir, "output", "expr_genes"), recursive = TRUE)

	write.table(expressed_transcripts_th_avg, file=file.path(out_dir, "output", "expr_genes", fname_expr_transcripts), 
		sep=";", dec=".", row.names=FALSE, col.names=TRUE)

}



thresholds <- list(0, 0.01, 0.05, 0.1) # ERCC threshold values are 0 (a fixed TPM threshold of 0.5)
# or perc of expressed spike-ins for 0.01/0.05/0.1

lapply(thresholds, getExprTranscripts, species = "ATH", experiment = "single-species")
lapply(thresholds, getExprTranscripts, species = "AL", experiment = "comparative")
lapply(thresholds, getExprTranscripts, species = "CR", experiment = "comparative")
lapply(thresholds, getExprTranscripts, species = "ES", experiment = "comparative")
lapply(thresholds, getExprTranscripts, species = "TH", experiment = "comparative")
lapply(thresholds, getExprTranscripts, species = "MT", experiment = "comparative")
lapply(thresholds, getExprTranscripts, species = "BD", experiment = "comparative")


