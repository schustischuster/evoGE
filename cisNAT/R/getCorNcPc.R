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
        genesCount = file.path(in_dir, "Expression_data", "AT_genes_inter_norm_count_mat_vsd_sample_names.csv")
        species_id <- "AT"

    } else if (is.element("AL", species)) {
		GTFfile = file.path(in_dir, "GTF", "AL_final_annotation.gtf")
		genesTPM = file.path(in_dir, "Expression_data", "AL_lnc_all_antisense_genes_inter_tpm_mat_deseq_sample_ids_extended.csv")
		genesCount = file.path(in_dir, "Expression_data", "AL_genes_inter_norm_count_mat_vsd_sample_names.csv")
		species_id <- "AL"

    } else if (is.element("CR", species)) {
		GTFfile = file.path(in_dir, "GTF", "CR_final_annotation.gtf")
		genesTPM = file.path(in_dir, "Expression_data", "CR_lnc_all_antisense_genes_inter_tpm_mat_deseq_sample_ids_extended.csv")
		genesCount = file.path(in_dir, "Expression_data", "CR_genes_inter_norm_count_mat_vsd_sample_names.csv")
		species_id <- "CR"

    } else if (is.element("ES", species)) {
		GTFfile = file.path(in_dir, "GTF", "ES_final_annotation.gtf")
		genesTPM = file.path(in_dir, "Expression_data", "ES_lnc_all_antisense_genes_inter_tpm_mat_deseq_sample_ids_extended.csv")
		genesCount = file.path(in_dir, "Expression_data", "ES_genes_inter_norm_count_mat_vsd_sample_names.csv")
		species_id <- "ES"

    } else if (is.element("TH", species)) {
		GTFfile = file.path(in_dir, "GTF", "TH_final_annotation.gtf")
		genesTPM = file.path(in_dir, "Expression_data", "TH_lnc_all_antisense_genes_inter_tpm_mat_deseq_sample_ids_extended.csv")
		genesCount = file.path(in_dir, "Expression_data", "TH_genes_inter_norm_count_mat_vsd_sample_names.csv")
		species_id <- "TH"

    } else if (is.element("MT", species)) {
		GTFfile = file.path(in_dir, "GTF", "MT_final_annotation.gtf")
		genesTPM = file.path(in_dir, "Expression_data", "MT_lnc_all_antisense_genes_inter_tpm_mat_deseq_sample_ids_extended.csv")
		genesCount = file.path(in_dir, "Expression_data", "MT_genes_inter_norm_count_mat_vsd_sample_names.csv")
		species_id <- "MT"

    } else if (is.element("BD", species)) {
		GTFfile = file.path(in_dir, "GTF", "BD_final_annotation.gtf")
		genesTPM = file.path(in_dir, "Expression_data", "BD_lnc_all_antisense_genes_inter_tpm_mat_deseq_sample_ids_extended.csv")
		genesCount = file.path(in_dir, "Expression_data", "BD_genes_inter_norm_count_mat_vsd_sample_names.csv")
		species_id <- "BD"
    }


	# Import gtf file
	GTF = import.gff(GTFfile, format = "gtf", feature.type = "gene")

	# Read expression data
	all_genes_tpm <- read.table(genesTPM, sep = ";", dec = ".", header = TRUE, stringsAsFactors = FALSE)
	all_genes_count <- read.table(genesCount, sep = ";", dec = ".", header = TRUE, stringsAsFactors = FALSE)
	colnames(all_genes_tpm)[1] <- "gene_id"
	all_genes_count <- tibble::rownames_to_column(all_genes_count, "gene_id")


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

		all_genes_count <- dplyr::select(all_genes_count, c(
			gene_id,  
			root_whole_root_5d_.1., # root 5d.1
			root_whole_root_5d_.2., # root 5d.2
			root_whole_root_5d_.3., # root 5d.3
			hypocotyl_10d_.1., # hypocotyl 10d.1
			hypocotyl_10d_.2., # hypocotyl 10d.2
			hypocotyl_10d_.3., # hypocotyl 10d.3
			leaf_1.2_7d_.1., # leaf 1+2 7d.1
			leaf_1.2_7d_.2., # leaf 1+2 7d.2
			leaf_1.2_7d_.3., # leaf 1+2 7d.3
			apex_vegetative_7d_.1., # apex veg 7d.1
			apex_vegetative_7d_.2., # apex veg 7d.2
			apex_vegetative_7d_.3., # apex veg 7d.3
			apex_inflorescence_21d_.1., # apex inf 21d.1
			apex_inflorescence_21d_.2., # apex inf 21d.2
			apex_inflorescence_21d_.3., # apex inf 21d.3
			flower_stg12_21d._.1., # flower stg12 21d.1
			flower_stg12_21d._.2., # flower stg12 21d.2
			flower_stg12_21d._.3., # flower stg12 21d.3
			flower_stg12_stamens_21d._.1., # flower stg12 stamen.1
			flower_stg12_stamens_21d._.2., # flower stg12 stamen.2
			flower_stg12_stamens_21d._.3., # flower stg12 stamen.3
			flower_early_stg12_carpels_21d._.1., # flower stg12 carpel.1
			flower_early_stg12_carpels_21d._.2., # flower stg12 carpel.2
			flower_early_stg12_carpels_21d._.3.)) # flower stg12 carpel.3
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

		all_genes_count <- dplyr::select(all_genes_count, -c(
			flower_stg11_stamens_8w.10w.25d_.1., # flower stg11 stamen.1 
			flower_stg11_stamens_8w.10w.25d_.2., # flower stg11 stamen.2 
			flower_stg11_stamens_8w.10w.25d_.3., # flower stg11 stamen.3
			flower_early_stg12_stamens_8w.10w.23d_.1., # flower early stg12 stamen.1
			flower_early_stg12_stamens_8w.10w.23d_.2., # flower early stg12 stamen.2
			flower_early_stg12_stamens_8w.10w.23d_.3., # flower early stg12 stamen.3
			flower_late_stg12_stamens_8w.10w.21d_.1., # flower late stg12 stamen.1
			flower_late_stg12_stamens_8w.10w.21d_.2., # flower late stg12 stamen.2
			flower_late_stg12_stamens_8w.10w.21d_.3.)) # flower late stg12 stamen.3
		### tibble w/o pollen samles

		species_id <- "AL_comparative_samples"

	}


    # Stop function here to allow specific analysis of a single species
    # return_list <- list("species_id" = species_id, "GTF" = GTF, "all_genes_tpm" = all_genes_tpm, "all_genes_count" = all_genes_count, "experiment" = experiment)
    # return(return_list)
    # }
    # return_objects <- getCorNcPc("AT", "single-species") # read in GTF and expression data for A.thaliana
    # list2env(return_objects, envir = .GlobalEnv)



	#------------------- Compute correlation between coding gene and cisNAT  -------------------

	if (species_id == "AT_comparative_samples" | species_id == "AL_comparative_samples") {

		all_genes_tpm <- all_genes_tpm[rowSums(all_genes_tpm[, 6:ncol(all_genes_tpm)]) > 0, ]

	} # Remove all genes that are not expressed in comparative samples

	all_genes_tpm[6:ncol(all_genes_tpm)] <- log2(all_genes_tpm[6:ncol(all_genes_tpm)] + 1)

	# Replace TPM expression values by VST counts
	all_genes_count <- merge(all_genes_count, all_genes_tpm, by = "gene_id")

	if (species_id == "AT") {
		
		all_genes_count <- cbind(all_genes_count[c("gene_id", "prt_id", "biotype", "source", "info")], 
			all_genes_count[2:130])
	
	} else if (species_id == "AL") {
		
		all_genes_count <- cbind(all_genes_count[c("gene_id", "prt_id", "biotype", "source", "info")], 
			all_genes_count[2:34])
	
	} else {

		all_genes_count <- cbind(all_genes_count[c("gene_id", "prt_id", "biotype", "source", "info")], 
		all_genes_count[2:25])
	}


	all_genes_tpm_ls <- split(all_genes_tpm, f = all_genes_tpm$prt_id)
	all_genes_count_ls <- split(all_genes_count, f = all_genes_count$prt_id)
	# Split data by same protein-coding gene id
	# 4100 protein-coding gene/cisNAT pairs for AT comparative (tpm)

	all_genes_tpm_ls <- Filter(function(dt) nrow(dt) > 1, all_genes_tpm_ls)
	all_genes_count_ls <- Filter(function(dt) nrow(dt) > 1, all_genes_count_ls)
	# Remove all list entries that contain unpaired entries
	# 3843 protein-coding gene/cisNAT pairs for AT comparative (tpm)


	all_genes_tpm_ls <- Filter(function(db) nrow(db[db$biotype == "protein_coding",]) > 0, all_genes_tpm_ls)
	all_genes_count_ls <- Filter(function(db) nrow(db[db$biotype == "protein_coding",]) > 0, all_genes_count_ls)
	# Remove all list elements that do not contain a protein-coding gene partner


	length(all_genes_tpm_ls[purrr::map(all_genes_tpm_ls, nrow) > 2])
	# check how many genes have more than one NAT overlapping = 150 for AT

	length(all_genes_tpm_ls[purrr::map(all_genes_tpm_ls, nrow) > 3])
	# 7 genes overlapping 3 cisNATs for AT


	# Show message
	message("Calculate correlations...")


	# Compute spearman and pearson correlations
	getCdNcCor <- function(df) {

		# no log-transformation required for Pearson correlation estimation (VST)

		corSP <- function(p) {
				spe_df <- cor.test(as.numeric(df[1, 6:ncol(df)]), as.numeric(df[p, 6:ncol(df)]), method = "spearman")$estimate
				pea_df <- cor.test(as.numeric(df[1, 6:ncol(df)]), as.numeric(df[p, 6:ncol(df)]), method = "pearson")$estimate
				df_out <- rbind(spe_df, pea_df)
			}

		if (nrow(df) == 2) {

			cor_spe <- cor.test(as.numeric(df[1, 6:ncol(df)]), as.numeric(df[2, 6:ncol(df)]), method = "spearman")$estimate
			cor_pea <- cor.test(as.numeric(df[1, 6:ncol(df)]), as.numeric(df[2, 6:ncol(df)]), method = "pearson")$estimate
			max <- apply(X = df[6:ncol(df)], MARGIN = 1, FUN = max)
			avg <- apply(X = df[6:ncol(df)], MARGIN = 1, FUN = mean)
			df_out <- data.frame(df[!grepl("protein_coding", df$biotype), 1:5], Spearman = cor_spe, Pearson = cor_pea , maxPC = max[1], maxNC = max[2], 
				meanPC = avg[1], meanNC = avg[2], maxRatio = max[2]/max[1], meanRatio = avg[2]/avg[1], 
				maxSum = max[2]+max[1], meanSum = avg[2]+avg[1])
		
		} else if (nrow(df) == 3) {

			r_list <- seq(2, 3)

			p_df <- do.call(rbind, lapply(r_list, corSP))

			max <- apply(X = df[6:ncol(df)], MARGIN = 1, FUN = max)
			avg <- apply(X = df[6:ncol(df)], MARGIN = 1, FUN = mean)

			df_out <- data.frame(df[!grepl("protein_coding", df$biotype), 1:5], Spearman = p_df[rownames(p_df) == "spe_df"], Pearson = p_df[rownames(p_df) == "pea_df"], 
				maxPC = rep(max[1], 2), maxNC = max[2:3], meanPC = rep(avg[1], 2), meanNC = avg[2:3], 
				maxRatio = c(max[2]/max[1], max[3]/max[1]), meanRatio = c(avg[2]/avg[1], avg[3]/avg[1]), 
				maxSum = c(max[2]+max[1], max[3]+max[1]), meanSum = c(avg[2]+avg[1], avg[3]+avg[1]))
		
		} else if (nrow(df) == 4) {

			r_list <- seq(2, 4)

			p_df <- do.call(rbind, lapply(r_list, corSP))

			max <- apply(X = df[6:ncol(df)], MARGIN = 1, FUN = max)
			avg <- apply(X = df[6:ncol(df)], MARGIN = 1, FUN = mean)

			df_out <- data.frame(df[!grepl("protein_coding", df$biotype), 1:5], Spearman = p_df[rownames(p_df) == "spe_df"], Pearson = p_df[rownames(p_df) == "pea_df"], 
				maxPC = rep(max[1], 3), maxNC = max[2:4], meanPC = rep(avg[1], 3), meanNC = avg[2:4], 
				maxRatio = c(max[2]/max[1], max[3]/max[1], max[4]/max[1]), meanRatio = c(avg[2]/avg[1], avg[3]/avg[1], avg[4]/avg[1]), 
				maxSum = c(max[2]+max[1], max[3]+max[1], max[4]+max[1]), meanSum = c(avg[2]+avg[1], avg[3]+avg[1], avg[4]+avg[1]))
		
		} else if (nrow(df) == 5) {

			r_list <- seq(2, 5)

			p_df <- do.call(rbind, lapply(r_list, corSP))

			max <- apply(X = df[6:ncol(df)], MARGIN = 1, FUN = max)
			avg <- apply(X = df[6:ncol(df)], MARGIN = 1, FUN = mean)

			df_out <- data.frame(df[!grepl("protein_coding", df$biotype), 1:5], Spearman = p_df[rownames(p_df) == "spe_df"], Pearson = p_df[rownames(p_df) == "pea_df"], 
				maxPC = rep(max[1], 4), maxNC = max[2:5], meanPC = rep(avg[1], 4), meanNC = avg[2:5], 
				maxRatio = c(max[2]/max[1], max[3]/max[1], max[4]/max[1], max[5]/max[1]), meanRatio = c(avg[2]/avg[1], avg[3]/avg[1], avg[4]/avg[1], avg[5]/avg[1]), 
				maxSum = c(max[2]+max[1], max[3]+max[1], max[4]+max[1], max[5]+max[1]), meanSum = c(avg[2]+avg[1], avg[3]+avg[1], avg[4]+avg[1], avg[5]+avg[1]))
		
		} else if (nrow(df) == 6) {

			r_list <- seq(2, 6)

			p_df <- do.call(rbind, lapply(r_list, corSP))

			max <- apply(X = df[6:ncol(df)], MARGIN = 1, FUN = max)
			avg <- apply(X = df[6:ncol(df)], MARGIN = 1, FUN = mean)

			df_out <- data.frame(df[!grepl("protein_coding", df$biotype), 1:5], Spearman = p_df[rownames(p_df) == "spe_df"], Pearson = p_df[rownames(p_df) == "pea_df"], 
				maxPC = rep(max[1], 5), maxNC = max[2:6], meanPC = rep(avg[1], 5), meanNC = avg[2:6], 
				maxRatio = c(max[2]/max[1], max[3]/max[1], max[4]/max[1], max[5]/max[1], max[6]/max[1]), meanRatio = c(avg[2]/avg[1], avg[3]/avg[1], avg[4]/avg[1], avg[5]/avg[1], avg[6]/avg[1]), 
				maxSum = c(max[2]+max[1], max[3]+max[1], max[4]+max[1], max[5]+max[1], max[6]+max[1]), meanSum = c(avg[2]+avg[1], avg[3]+avg[1], avg[4]+avg[1], avg[5]+avg[1], avg[6]+avg[1]))
		
		} else if (nrow(df) == 7) {

			r_list <- seq(2, 7)

			p_df <- do.call(rbind, lapply(r_list, corSP))

			max <- apply(X = df[6:ncol(df)], MARGIN = 1, FUN = max)
			avg <- apply(X = df[6:ncol(df)], MARGIN = 1, FUN = mean)

			df_out <- data.frame(df[!grepl("protein_coding", df$biotype), 1:5], Spearman = p_df[rownames(p_df) == "spe_df"], Pearson = p_df[rownames(p_df) == "pea_df"], 
				maxPC = rep(max[1], 6), maxNC = max[2:7], meanPC = rep(avg[1], 6), meanNC = avg[2:7], 
				maxRatio = c(max[2]/max[1], max[3]/max[1], max[4]/max[1], max[5]/max[1], max[6]/max[1], max[7]/max[1]), meanRatio = c(avg[2]/avg[1], avg[3]/avg[1], avg[4]/avg[1], avg[5]/avg[1], avg[6]/avg[1], avg[7]/avg[1]), 
				maxSum = c(max[2]+max[1], max[3]+max[1], max[4]+max[1], max[5]+max[1], max[6]+max[1], max[7]+max[1]), meanSum = c(avg[2]+avg[1], avg[3]+avg[1], avg[4]+avg[1], avg[5]+avg[1], avg[6]+avg[1], avg[7]+avg[1]))
		
		} else if (nrow(df) == 8) {

			r_list <- seq(2, 8)

			p_df <- do.call(rbind, lapply(r_list, corSP))

			max <- apply(X = df[6:ncol(df)], MARGIN = 1, FUN = max)
			avg <- apply(X = df[6:ncol(df)], MARGIN = 1, FUN = mean)

			df_out <- data.frame(df[!grepl("protein_coding", df$biotype), 1:5], Spearman = p_df[rownames(p_df) == "spe_df"], Pearson = p_df[rownames(p_df) == "pea_df"], 
				maxPC = rep(max[1], 7), maxNC = max[2:8], meanPC = rep(avg[1], 7), meanNC = avg[2:8], 
				maxRatio = c(max[2]/max[1], max[3]/max[1], max[4]/max[1], max[5]/max[1], max[6]/max[1], max[7]/max[1], max[8]/max[1]), meanRatio = c(avg[2]/avg[1], avg[3]/avg[1], avg[4]/avg[1], avg[5]/avg[1], avg[6]/avg[1], avg[7]/avg[1], avg[8]/avg[1]), 
				maxSum = c(max[2]+max[1], max[3]+max[1], max[4]+max[1], max[5]+max[1], max[6]+max[1], max[7]+max[1], max[8]+max[1]), meanSum = c(avg[2]+avg[1], avg[3]+avg[1], avg[4]+avg[1], avg[5]+avg[1], avg[6]+avg[1], avg[7]+avg[1], avg[8]+avg[1]))
		
		}

		return(df_out)
	}

	# startup message
	message("Computing correlation...")

	cd_nc_cor_count <- suppressWarnings(do.call("rbind", lapply(all_genes_count_ls, getCdNcCor)))
	cd_nc_cor_tpm <- suppressWarnings(do.call("rbind", lapply(all_genes_tpm_ls, getCdNcCor)))

	cd_nc_cor_count <- cd_nc_cor_count[cd_nc_cor_count$biotype != "protein_coding",]
	cd_nc_cor_tpm <- cd_nc_cor_tpm[cd_nc_cor_tpm$biotype != "protein_coding",]



	# Show message
	message("Write data...")

	# Write final data tables to csv files and store them in /out_dir/output/plots
	if (!dir.exists(file.path(out_dir, "output", "NAT_expr_cor"))) 
		dir.create(file.path(out_dir, "output", "NAT_expr_cor"), recursive = TRUE)

	fname_count <- sprintf('%s.csv', paste(species_id, "cd_nc_cor_count", sep = "_"))
	fname_tpm <- sprintf('%s.csv', paste(species_id, "cd_nc_cor_tpm", sep = "_"))

	write.table(cd_nc_cor_count, file = file.path(out_dir, "output", "NAT_expr_cor", fname_count), 
		sep = ";", dec = ".", row.names = FALSE, col.names = TRUE)
	write.table(cd_nc_cor_tpm, file = file.path(out_dir, "output", "NAT_expr_cor", fname_tpm), 
		sep = ";", dec = ".", row.names = FALSE, col.names = TRUE)


}

