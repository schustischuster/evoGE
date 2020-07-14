# Compute expression correlation of 10.000 randomized gene pairs in ATH
# Repeat for number of given bootstrap replicates 
# Data input: Expression_data WITHOUT mito and chloroplast genes


#------------------- Load packages, set directories and read sample tables ---------------------


# Install and load packages
if (!require(dplyr)) install.packages('dplyr')
library(dplyr)


# Set file path and input files
in_dir <- "/Volumes/User/Shared/Christoph_manuscript/DevSeq_paper/Analysis/Analysis_2019/A_thaliana_gene_exression_map/20191121_CS_coding_cisNAT_analysis/data"
out_dir <- "/Volumes/User/Shared/Christoph_manuscript/DevSeq_paper/Analysis/Analysis_2019/A_thaliana_gene_exression_map/20191121_CS_coding_cisNAT_analysis"


getRandGeneCor <- function(species = c("ATH", "AL", "CR", "ES", "TH", "MT", "BD"), 
    cor_method = c("Spearman", "Pearson"), experiment = c("single-species", "comparative"), 
    bootstrap_repl) {

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

    # Show error message if no species is chosen
    if (missing(cor_method))
   
       stop(
       "Please choose one of the available cor_method methods: 
	   'Spearman', 'Pearson'",
	   call. = TRUE
       )

    # Show error message if incorrect correlation coefficient is chosen
    if (!is.element(cor_method, c("Spearman", "Pearson")))
   
       stop(
       "Please choose one of the available cor_method methods: 
	   'Spearman', 'Pearson'",
	   call. = TRUE
       )

    # Show error message if no number of bootstrap replicates is given
    if (missing(bootstrap_repl))
   
       stop(
       "Please provide number for bootstrap_repl",
	   call. = TRUE
       )


    # Set expression data input file
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
    # return_list <- list("species_id" = species_id, "all_genes_tpm" = all_genes_tpm, "cor_method" = cor_method, "bootstrap_repl" = bootstrap_repl)
    # return(return_list)
    # }
    # return_objects <- getRandGeneCor(species = "ATH", cor_method = "Pearson", experiment = "single-species", bootstrap_repl = 10) # read in GTF and expression data for A.thaliana
    # list2env(return_objects, envir = .GlobalEnv)



    #------------------- Select protein-coding genes and remove pollen data  -------------------


	# Remove pollen triplicates from expression table
	all_genes_tpm_wo_pollen <- dplyr::select(all_genes_tpm, -c(
		flowers_mature_pollen_1, 
		flowers_mature_pollen_2, 
		flowers_mature_pollen_3)) #tibble w/o pollen samles


	# Extract all protein-coding genes
	all_genes_tpm_wo_pollen_cd <- subset(all_genes_tpm_wo_pollen, biotype == "protein_coding")

	


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
	
			th_replicates <- do.call(cbind, lapply(split.default(df[4:ncol(df)], #adjust columns
								rep(seq_along(df), each = 3, length.out = ncol(df)-3)), #adjust columns
								function(x) {
									x[rowSums(x > threshold) < 2, ] <- 0; 
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
		th_df <- merge(keys_data, df, by="key")
		th_df <- th_df[order(th_df$key),]
		th_df <- th_df[-1:-2]

		return(th_df)
	}


	# Apply threshold function
	all_genes_tpm_0.5 <- applyThreshold(all_genes_tpm_wo_pollen_cd, 0.5)




    #--------------- Select 10.000 random gene pairs and compute correlation -------------------



	randomizeGene <- function(df, cor_method) {

        # Define two sets of each 10.000 genes
        gene_set_1 <- df[sample(nrow(df), 10000, replace = FALSE), ]

        remaining_genes <- df[!(df$gene_id %in% gene_set_1$gene_id),]
        gene_set_2 <- remaining_genes[sample(nrow(remaining_genes), 10000, replace = FALSE), ]


        # Compute pairwise gene correlation
        getCor <- function(df1, df2, method) {

		    df1_col <- ncol(df1)
		    df2_col <- ncol(df2)

		    # startup message
		    message("Computing correlation...")

		    if (method == "Spearman") {

		        cor_out <- sapply(1:nrow(df1), function(i) 
		    	    cor(as.numeric(df1[i, 4:df1_col]), as.numeric(df2[i, 4:df2_col]), method=c("spearman")))
		
		    } else if (method == "Pearson") {

		        # log2-transform TPM values before computing Pearson
		        df1[, 4:df1_col] <- log2(df1[, 4:df1_col] + 1)
		        df2[, 4:df1_col] <- log2(df2[, 4:df2_col] + 1)

		        cor_out <- sapply(1:nrow(df1), function(i) 
		    	    cor(as.numeric(df1[i, 4:df1_col]), as.numeric(df2[i, 4:df2_col]), method=c("pearson")))
		    }
	    }

        cor <- getCor(gene_set_1, gene_set_2, method = cor_method)

        return(cor)
    }

    bs_cor_PCT_pairs <- replicate(bootstrap_repl, randomizeGene(all_genes_tpm_0.5, cor_method), simplify=FALSE)


    # Write bootstrap replicate correlation values result list to data frame
    bs_cor_PCT_pairs_df <- data.frame(matrix(unlist(bs_cor_PCT_pairs), nrow=10000, byrow=F), 
    	                              stringsAsFactors=FALSE)


    # Set colnames
    bs_number <- seq(1:bootstrap_repl)
    bs_tag <- rep("bs_repl_", bootstrap_repl)
    bs_col_names <- paste0(bs_tag, bs_number)

    colnames(bs_cor_PCT_pairs_df) <- bs_col_names




    #----------------------------------- Write output to csv -----------------------------------


	# Set filename
    fname <- sprintf('%s.csv', paste("PCT_corr", species_id, cor_method, bootstrap_repl, "BS_repl", sep="_"))

	# Write final data tables to csv files and store them in /out_dir/output/data_tables
	if (!dir.exists(file.path(out_dir, "output", "cd_gene_pairs"))) 
		dir.create(file.path(out_dir, "output", "cd_gene_pairs"), recursive = TRUE)

	write.table(bs_cor_PCT_pairs_df, file=file.path(out_dir, "output", "cd_gene_pairs", fname), 
		sep=";", dec=".", row.names=FALSE, col.names=TRUE)

}



# Execute getRandGeneCor function
getRandGeneCor(species="ATH", cor_method="Pearson", experiment="single-species", bootstrap_repl=10)






