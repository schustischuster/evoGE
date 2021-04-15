# Find expressed genes in each sample replicate using ERCC spike-in information
# Thresholds: 0 (no ERCC threshold-everything above 0.5 TPM), 0.01 (weakest expressed spike-in), 
# 0.05 (lowest 5% of expressed spike-ins), 0.1 (lowest 10% of expressed spike-ins)
# Data input: 1) Expression_data WITH spike-in information
# Analysis can be performed on both whole single species datasets (ATH: 132 samples; AL: 36 samples)
# OR on comparative data sets (27 samples)
# Input sample tables should have the following format:
# id / biotype / source / info / DEVSEQ_SAMPLE_REPLICATES(between 27 and 132 depending on species)




#------------------- Load packages, set directories and read sample tables ---------------------


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
        genesTPM = file.path(in_dir, "Expression_data", "AT_genes_complete_table_tpm_with_circRNA_sample_names.csv")
        genesCounts = file.path(in_dir, "Expression_data", "AT_genes_intra_norm_count_mat_vsd_sample_names.csv")
        genesPtMt = file.path(in_dir, "Expression_data", "AT_Pt_Mt_Orthologs.csv")
        species_id <- "ATH"

    } else if (is.element("AL", species)) {
		genesTPM = file.path(in_dir, "Expression_data", "AL_genes_complete_table_tpm_with_circRNA_sample_names.csv")
		genesCounts = file.path(in_dir, "Expression_data", "AL_genes_intra_norm_count_mat_vsd_sample_names.csv")
		genesPtMt = file.path(in_dir, "Expression_data", "AL_Pt_Mt_Orthologs.csv")
		species_id <- "AL"

    } else if (is.element("CR", species)) {
		genesTPM = file.path(in_dir, "Expression_data", "CR_genes_complete_table_tpm_with_circRNA_sample_names.csv")
		genesCounts = file.path(in_dir, "Expression_data", "CR_genes_intra_norm_count_mat_vsd_sample_names.csv")
		genesPtMt = file.path(in_dir, "Expression_data", "CR_Pt_Mt_Orthologs.csv")
		species_id <- "CR"

    } else if (is.element("ES", species)) {
		genesTPM = file.path(in_dir, "Expression_data", "ES_genes_complete_table_tpm_with_circRNA_sample_names.csv")
		genesCounts = file.path(in_dir, "Expression_data", "ES_genes_intra_norm_count_mat_vsd_sample_names.csv")
		genesPtMt = file.path(in_dir, "Expression_data", "ES_Pt_Mt_Orthologs.csv")
		species_id <- "ES"

    } else if (is.element("TH", species)) {
		genesTPM = file.path(in_dir, "Expression_data", "TH_genes_complete_table_tpm_with_circRNA_sample_names.csv")
		genesCounts = file.path(in_dir, "Expression_data", "TH_genes_intra_norm_count_mat_vsd_sample_names.csv")
		genesPtMt = file.path(in_dir, "Expression_data", "TH_Pt_Mt_Orthologs.csv")
		species_id <- "TH"

    } else if (is.element("MT", species)) {
		genesTPM = file.path(in_dir, "Expression_data", "MT_genes_complete_table_tpm_with_circRNA_sample_names.csv")
		genesCounts = file.path(in_dir, "Expression_data", "MT_genes_intra_norm_count_mat_vsd_sample_names.csv")
		genesPtMt = file.path(in_dir, "Expression_data", "MT_Pt_Mt_Orthologs.csv")
		species_id <- "MT"

    } else if (is.element("BD", species)) {
		genesTPM = file.path(in_dir, "Expression_data", "BD_genes_complete_table_tpm_with_circRNA_sample_names.csv")
		genesCounts = file.path(in_dir, "Expression_data", "BD_genes_intra_norm_count_mat_vsd_sample_names.csv")
		genesPtMt = file.path(in_dir, "Expression_data", "BD_Pt_Mt_Orthologs.csv")
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
	genes_ptmt <- read.table(genesPtMt, sep=";", dec=".", header=TRUE, stringsAsFactors=FALSE)


	# Format expression data and rename pollen samples
    if ((is.element("ATH", species)) && (is.element("comparative", experiment))) {

		all_genes_tpm <- dplyr::select(all_genes_tpm, c(
			gene_id, biotype, source, info, 
			root_whole_root_5d_.1.,
			root_whole_root_5d_.2.,
			root_whole_root_5d_.3.,
			hypocotyl_10d_.1.,
			hypocotyl_10d_.2.,
			hypocotyl_10d_.3.,
			leaf_1.2_7d_.1.,
			leaf_1.2_7d_.2.,
			leaf_1.2_7d_.3.,
			apex_vegetative_7d_.1.,
			apex_vegetative_7d_.2.,
			apex_vegetative_7d_.3.,
			apex_inflorescence_21d_.1.,
			apex_inflorescence_21d_.2.,
			apex_inflorescence_21d_.3.,
			flower_stg12_21d._.1.,
			flower_stg12_21d._.2.,
			flower_stg12_21d._.3.,
			flower_stg12_stamens_21d._.1.,
			flower_stg12_stamens_21d._.2.,
			flower_stg12_stamens_21d._.3.,
			flowers_mature_pollen_28d_.1.,
			flowers_mature_pollen_28d_.2.,
			flowers_mature_pollen_28d_.3.,
			flower_early_stg12_carpels_21d._.1.,
			flower_early_stg12_carpels_21d._.2.,
			flower_early_stg12_carpels_21d._.3.)) #tibble w/o pollen samples

		species_id <- "ATH_comparative_samples"


    } else if ((is.element("AL", species)) && (is.element("comparative", experiment))) {

		all_genes_tpm <- dplyr::select(all_genes_tpm, -c(
			flower_stg11_stamens_8w.10w.25d_.1., 
			flower_stg11_stamens_8w.10w.25d_.2., 
			flower_stg11_stamens_8w.10w.25d_.3.,
			flower_early_stg12_stamens_8w.10w.23d_.1.,
			flower_early_stg12_stamens_8w.10w.23d_.2.,
			flower_early_stg12_stamens_8w.10w.23d_.3.,
			flower_late_stg12_stamens_8w.10w.21d_.1.,
			flower_late_stg12_stamens_8w.10w.21d_.2.,
			flower_late_stg12_stamens_8w.10w.21d_.3.)) #tibble w/o pollen samples

    }



    # Stop function here to allow specific analysis of a single species
    # return_list <- list("species_id" = species_id, "experiment" = experiment, "all_genes_tpm" = all_genes_tpm, "all_genes_counts" = all_genes_counts, "threshold" = threshold, "genes_ptmt" = genes_ptmt)
    # return(return_list)
    # }
    # return_objects <- getExprGenes(species="ATH", experiment="single-species", 0.05) # read in GTF and expression data for A.thaliana
    # list2env(return_objects, envir = .GlobalEnv)



#--------------------------- Get and apply sample-specific threshold  --------------------------


    # Show message
    message("Starting analysis...")


    # Remove info column
    all_genes_tpm <- dplyr::select(all_genes_tpm, -c(info))
    

	# Extract ERCC data
	ERCC <- all_genes_tpm[all_genes_tpm$gene_id %like% "ERCC", ]


	if (threshold > 0) {

	   # Get sample-specific TPM threshold
	   ERCC_cutoff <- cbind(
		   data.frame(gene_id="TPM_cutoff"), data.frame(biotype="<NA>"), data.frame(source="<NA>"), 
		   as.data.frame(t(data.frame(
			TPM_cutoff = apply(ERCC[,4:ncol(ERCC)], 2, function(i)quantile(i[i>0], threshold)))
		   ))
	   )

	   # Bind sample-specific threshold to expression table
	   all_genes_tpm_cutoff <- rbind(all_genes_tpm, ERCC_cutoff)
	}


	# Set static threshold of 0.5 TPM for "0 ERCC threshold"
	if (threshold == 0) {
		th_values <- data.frame(t(rep(0.5, ncol(ERCC)-3)))
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

	
	# Set replicate colnames
	repl_names <- unique(substr(names(express_data_th_replicates[-1:-3]), 1, nchar(names(express_data_th_replicates[-1:-3]))-2))
	repl_names <- gsub("^[^.]*.", "", repl_names)
	repl_names <- substr(repl_names, 1, nchar(repl_names)-2)
	repl_names <- gsub('^\\.|\\.$', '', repl_names)

	colnames(express_data_th_avg) <- c("gene_id", "biotype", "source", repl_names)


	# Check which biotypes are in df
	express_data_th_avg$biotype[!(duplicated(express_data_th_avg$biotype))]

	
	# Remove Pt and Mt genes
	express_data_th_avg <- express_data_th_avg[!(express_data_th_avg$gene_id %in% genes_ptmt$gene_id),]


	# Get number of genes expressed in each sample type
	protein_coding_subset <- subset(express_data_th_avg, biotype=="protein_coding")
	lnc_intergenic_subset <- subset(express_data_th_avg, biotype=="lnc_intergenic")
	lnc_antisense_subset <- express_data_th_avg[express_data_th_avg$biotype %like% "antisense", ]
	LTR_subset <- subset(express_data_th_avg, biotype=="LTR_retrotransposon")

	
	# Double-check that no chloroplast and mito genes from coding gene list
	`%nlike%` = Negate(`%like%`)

	if ((is.element("ATH", species_id))) { 

		protein_coding_subset <- protein_coding_subset[protein_coding_subset$gene_id %nlike% "ATMG", ]
		protein_coding_subset <- protein_coding_subset[protein_coding_subset$gene_id %nlike% "ATCG", ]
	}


	expr_protein_coding <- colSums(protein_coding_subset > 0)
	expr_lnc_antisense <- colSums(lnc_antisense_subset > 0)
	expr_lnc_intergenic <- colSums(lnc_intergenic_subset > 0)
	expr_LTR <- colSums(LTR_subset > 0)


	# Create final list of expressed genes per organ/ sample type
	expressed_genes_th_avg <- rbind(expr_protein_coding, expr_lnc_antisense, expr_lnc_intergenic, expr_LTR)
	expressed_genes_th_avg <- subset(expressed_genes_th_avg, select = -c(biotype, source))
	colnames(expressed_genes_th_avg)[1] <- "total_expressed"
	biotype <- substring(rownames(expressed_genes_th_avg), 6)
	expressed_genes_th_avg <- as.data.frame(cbind(biotype, expressed_genes_th_avg))




#--------------- Calculate correlation for biological replicates of each sample ----------------


	# Remove ERCC spike-ins from VST count thresholded expression table
	`%nlike%` = Negate(`%like%`)
	express_data_th_counts <- all_genes_counts[rownames(all_genes_counts) %nlike% "ERCC", ]


	# Get replicate correlation
	replCorr <- function(df, coefficient=c("spearman", "pearson")) {

		replicate_corr <- do.call(cbind, lapply(split.default(df[1:ncol(df)], #adjust columns
								rep(seq_along(df), each = 3, length.out = ncol(df))), #adjust columns
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


	# Get replicate correlations based on normalized VST counts
	repl_corr_df_counts <- replCorr(express_data_th_counts, coefficient="pearson")
	colnames(repl_corr_df_counts) <- colnames(express_data_th_counts)




#-------- Generate replicate expression lists with thresholded genes for each biotype ---------


	# For ERCC-thresholded VST count data
	protein_coding_th_repl_counts <- subset(all_genes_counts, rownames(all_genes_counts) %in% protein_coding_subset$gene_id)
	lnc_antisense_th_repl_counts <- subset(all_genes_counts, rownames(all_genes_counts) %in% lnc_antisense_subset$gene_id)
	lnc_intergenic_th_repl_counts <- subset(all_genes_counts, rownames(all_genes_counts) %in% lnc_intergenic_subset$gene_id)
	circRNA_th_repl_counts <- subset(all_genes_counts, rownames(all_genes_counts) %in% circRNA_subset$gene_id)
	th_genes_counts <- rbind(protein_coding_th_repl_counts, lnc_antisense_th_repl_counts, lnc_intergenic_th_repl_counts, 
		circRNA_th_repl_counts)




#--------------------------------------- Write csv file  ---------------------------------------


	# Show message
    message("Writing output...")


	# Set filename
    fname_expr_genes <- sprintf('%s.csv', paste(species_id, "expr_genes", threshold, sep="_"))
    fname_repl_corr_counts <- sprintf('%s.csv', paste(species_id, "repl_corr_counts", threshold, sep="_"))
    fname_th_genes_repl_counts <- sprintf('%s.csv', paste(species_id, "th_genes_repl_counts", threshold, sep="_"))


	# Write final data tables to csv files and store them in /out_dir/output/data_tables
	if (!dir.exists(file.path(out_dir, "output", "expr_genes"))) 
		dir.create(file.path(out_dir, "output", "expr_genes"), recursive = TRUE)

	write.table(expressed_genes_th_avg, file=file.path(out_dir, "output", "expr_genes", fname_expr_genes), 
		sep=";", dec=".", row.names=FALSE, col.names=TRUE)
	write.table(repl_corr_df_counts, file=file.path(out_dir, "output", "expr_genes", fname_repl_corr_counts), 
		sep=";", dec=".", row.names=FALSE, col.names=TRUE)

	# Only write thresholded expression tables for 0.5 TPM threshold to file
	if (threshold == 0) {
		write.table(th_genes_counts, file=file.path(out_dir, "output", "expr_genes", fname_th_genes_repl_counts), 
			sep=";", dec=".", row.names=TRUE, col.names=TRUE)
	}

}



thresholds <- list(0, 0.01, 0.05, 0.1) # ERCC threshold values are 0 (a fixed TPM threshold of 0.5)
# or perc of expressed spike-ins for 0.01/0.05/0.1

lapply(thresholds, getExprGenes, species = "ATH", experiment = "single-species")
lapply(thresholds, getExprGenes, species = "AL", experiment = "comparative")
lapply(thresholds, getExprGenes, species = "CR", experiment = "comparative")
lapply(thresholds, getExprGenes, species = "ES", experiment = "comparative")
lapply(thresholds, getExprGenes, species = "TH", experiment = "comparative")
lapply(thresholds, getExprGenes, species = "MT", experiment = "comparative")
lapply(thresholds, getExprGenes, species = "BD", experiment = "comparative")


