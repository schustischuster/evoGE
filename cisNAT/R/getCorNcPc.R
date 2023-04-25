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
	all_genes_count <- all_genes_count[-30:-53]
	all_genes_count <- cbind(all_genes_count[c("gene_id", "prt_id", "biotype", "source", "info")], 
		all_genes_count[2:25])


	all_genes_tpm_ls <- split(all_genes_tpm, f = all_genes_tpm$prt_id)
	all_genes_count_ls <- split(all_genes_count, f = all_genes_count$prt_id)
	# Split data by same protein-coding gene id
	# 4100 protein-coding gene/cisNAT pairs for AT comparative (tpm)

	all_genes_tpm_ls <- Filter(function(dt) nrow(dt) > 1, all_genes_tpm_ls)
	all_genes_count_ls <- Filter(function(dt) nrow(dt) > 1, all_genes_count_ls)
	# Remove all list entries that contain unpaired entries
	# 3843 protein-coding gene/cisNAT pairs for AT comparative (tpm)

	length(all_genes_tpm_ls[purrr::map(all_genes_tpm_ls, nrow) > 2])
	# check how many genes have more than one NAT overlapping = 150 for AT

	length(all_genes_tpm_ls[purrr::map(all_genes_tpm_ls, nrow) > 3])
	# 7 genes overlapping 3 cisNATs for AT


	# Show message
	message("Calculate correlations...")


	# Compute spearman and pearson correlations
	getCdNcCor <- function(df) {

		# no log-transformation required for Pearson correlation estimation (VST)

		if (nrow(df) == 2) {

			cor_spe <- cor.test(as.numeric(df[1, 6:ncol(df)]), as.numeric(df[2, 6:ncol(df)]), method = "spearman")$estimate
			cor_pea <- cor.test(as.numeric(df[1, 6:ncol(df)]), as.numeric(df[2, 6:ncol(df)]), method = "pearson")$estimate
			max <- apply(X = df[6:ncol(df)], MARGIN = 1, FUN = max)
			avg <- apply(X = df[6:ncol(df)], MARGIN = 1, FUN = mean)
			df_out <- data.frame(df[!grepl("protein_coding", df$biotype), 1:5], Spearman = cor_spe, Pearson = cor_pea , maxPC = max[1], maxNC = max[2], 
				meanPC = avg[1], meanNC = avg[2], maxRatio = max[2]/max[1], meanRatio = avg[2]/avg[1], 
				maxSum = max[2]+max[1], meanSum = avg[2]+avg[1])
		
		} else if (nrow(df) == 3) {

			cor1_spe <- cor.test(as.numeric(df[1, 6:ncol(df)]), as.numeric(df[2, 6:ncol(df)]), method = "spearman")$estimate
			cor2_spe <- cor.test(as.numeric(df[1, 6:ncol(df)]), as.numeric(df[3, 6:ncol(df)]), method = "spearman")$estimate
			cor1_pea <- cor.test(as.numeric(df[1, 6:ncol(df)]), as.numeric(df[2, 6:ncol(df)]), method = "pearson")$estimate
			cor2_pea <- cor.test(as.numeric(df[1, 6:ncol(df)]), as.numeric(df[3, 6:ncol(df)]), method = "pearson")$estimate
			max <- apply(X = df[6:ncol(df)], MARGIN = 1, FUN = max)
			avg <- apply(X = df[6:ncol(df)], MARGIN = 1, FUN = mean)
			df_out <- data.frame(df[!grepl("protein_coding", df$biotype), 1:5], Spearman = c(cor1_spe, cor2_spe), Pearson = c(cor1_pea, cor2_pea), 
				maxPC = rep(max[1], 2), maxNC = c(max[2], max[3]), meanPC = rep(avg[1], 2), meanNC = c(avg[2], avg[3]), 
				maxRatio = c(max[2]/max[1], max[3]/max[1]), meanRatio = c(avg[2]/avg[1], avg[3]/avg[1]), 
				maxSum = c(max[2]+max[1], max[3]+max[1]), meanSum = c(avg[2]+avg[1], avg[3]+avg[1]))
		
		} else if (nrow(df) == 4) {

			cor1_spe <- cor.test(as.numeric(df[1, 6:ncol(df)]), as.numeric(df[2, 6:ncol(df)]), method = "spearman")$estimate
			cor2_spe <- cor.test(as.numeric(df[1, 6:ncol(df)]), as.numeric(df[3, 6:ncol(df)]), method = "spearman")$estimate
			cor3_spe <- cor.test(as.numeric(df[1, 6:ncol(df)]), as.numeric(df[4, 6:ncol(df)]), method = "spearman")$estimate
			cor1_pea <- cor.test(as.numeric(df[1, 6:ncol(df)]), as.numeric(df[2, 6:ncol(df)]), method = "pearson")$estimate
			cor2_pea <- cor.test(as.numeric(df[1, 6:ncol(df)]), as.numeric(df[3, 6:ncol(df)]), method = "pearson")$estimate
			cor3_pea <- cor.test(as.numeric(df[1, 6:ncol(df)]), as.numeric(df[4, 6:ncol(df)]), method = "pearson")$estimate
			max <- apply(X = df[6:ncol(df)], MARGIN = 1, FUN = max)
			avg <- apply(X = df[6:ncol(df)], MARGIN = 1, FUN = mean)
			df_out <- data.frame(df[!grepl("protein_coding", df$biotype), 1:5], Spearman = c(cor1_spe, cor2_spe, cor3_spe), Pearson = c(cor1_pea, cor2_pea, cor3_pea), 
				maxPC = rep(max[1], 3), maxNC = c(max[2], max[3], max[4]), meanPC = rep(avg[1], 3), meanNC = c(avg[2], avg[3], avg[4]), 
				maxRatio = c(max[2]/max[1], max[3]/max[1], max[4]/max[1]), meanRatio = c(avg[2]/avg[1], avg[3]/avg[1], avg[4]/avg[1]), 
				maxSum = c(max[2]+max[1], max[3]+max[1], max[4]+max[1]), meanSum = c(avg[2]+avg[1], avg[3]+avg[1], avg[4]+avg[1]))
		
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

	fname_count <- sprintf('%s.csv', paste(species_id, experiment, "cd_nc_cor_count", sep = "_"))
	fname_tpm <- sprintf('%s.csv', paste(species_id, experiment, "cd_nc_cor_tpm", sep = "_"))

	write.table(cd_nc_cor_count, file = file.path(out_dir, "output", "NAT_expr_cor", fname_count), 
		sep = ";", dec = ".", row.names = FALSE, col.names = TRUE)
	write.table(cd_nc_cor_tpm, file = file.path(out_dir, "output", "NAT_expr_cor", fname_tpm), 
		sep = ";", dec = ".", row.names = FALSE, col.names = TRUE)



	# Show message
	message("Generate plots...")


	# Write final data tables to csv files and store them in /out_dir/output/plots
	if (!dir.exists(file.path(out_dir, "output", "plots"))) 
		dir.create(file.path(out_dir, "output", "plots"), recursive = TRUE)


	# Generate plots for comparative samples
	# AT single species correlation plot has to be generated independently, because it also
	# contains correlations of neighbouring and overlapping protein-coding genes

	if (species_id != "AT" | species_id != "AL") {


		# prepare data for ggplot2
		cor_AT <- data.frame(
			class = rep("PC_NC", nrow(cd_nc_cor)), 
			cor_NC_PC = cd_nc_cor$cor)



		# Generate plots
		plotCorAT <- function(data) {

			fname <- sprintf('%s.pdf', paste(deparse(substitute(data)), sep="_"))

			x_labels = c("PC_NC" = expression(atop(NA, atop(textstyle('PC/'), textstyle('NC')))))

			#data$cor <- factor(data$class, levels = unique(data$class))

			p <- ggplot(data, aes(x = class, y = cor_NC_PC, color = class)) + 
			geom_boxplot(aes(fill = class), colour = "black", width = 0.44, outlier.shape = NA, 
				size = 0.8, fatten = 2.8, notch = TRUE) + 
			scale_x_discrete(expand = c(0.005, 0), labels = x_labels) + 
			scale_y_continuous(limits = c(-1.05, 1), expand = c(0, 0), breaks = c(-1, -0.5, 0, 0.5, 1))

			q <- p + 
			scale_fill_manual(values = c("PC_NC" = "#f7ddb0")) + 
			theme_classic() + 
			xlab("Sense-Antisense Pair") + ylab("Correlation") + ggtitle("") + 
			theme(text = element_text(size = 23.5), 
				axis.ticks.length = unit(0.2, "cm"), 
				axis.ticks = element_line(colour = "black", size = 0.95), 
				axis.line = element_line(colour = 'black', size = 0.95), 
				plot.margin = unit(c(1, 1, 1, 1), "cm"), 
				axis.title.y = element_text(size = 18.4, margin = margin(t = 0, r = 6.4, b = 0, l = 3.38), 
					colour = "black", face = "plain"), 
				axis.title.x = element_text(size = 18.75, margin = margin(t = 6.5, r = 0, b = 5.75, l = 0), 
					colour = "black", face = "bold"), 
				axis.text.x = element_text(size = 16.25, margin = margin(t = -7, b = 2), colour = "black", 
					angle = 0, vjust = 1, hjust = 0.5), 
				axis.text.y = element_text(size = 16.5, angle = 0, margin = margin(l = 0, r = 1.5), colour = "black"),   
				legend.position = "none")

			ggsave(file = file.path(out_dir, "output", "plots", fname), plot = q, 
				width = 7, height = 5.75, units = c("in"))
		}

      plotCorAT(data = cor_AT)





	}



	#------------------------------------- Generate plots -------------------------------------


	


	










	# Generate plots
      plotMaxExprDist <- function(data, species) {

         if (species == "ACE") {

            p_mwu <- p_mwu[!grepl("AT", p_mwu$species),]

            # Create df for FDR p-value mapping
            mwu_df <- data.frame(
                class = rep(c("coding_non-core", "coding_core", "lncRNA_non-core", "lncRNA_core"), 
                    times = 3), 
                y = rep(c(19.68, 18.25), times = 6),
                label = ifelse(p_mwu$p_value < 1e-07, "****", 

                    c(paste("italic('P =')~", set_scientific(p_mwu$p_value)))), 

                species = rep(c("A.lyrata", "C.rubella", "E.salsugineum"), each = 4)
            )

            # Create df for gem_segments
            h_seg_df <- data.frame(
                x = rep(c(1.105, 2.105, 4.107, 5.107), times = 3), 
                xend = rep(c(3.105, 3.105, 6.107, 6.107), times = 3), 
                y = rep(c(20.08, 18.65, 20.08, 18.65), times = 3), 
                yend = rep(c(20.08, 18.65, 20.08, 18.65), times = 3), 
                species = rep(c("A.lyrata", "C.rubella", "E.salsugineum"), each = 4)
            )

            v_seg_df <- data.frame(
                x = rep(c(1.107, 3.107, 2.107, 3.107, 4.107, 6.107, 5.107, 6.107), times = 3), 
                xend = rep(c(1.107, 3.107, 2.107, 3.107, 4.107, 6.107, 5.107, 6.107), times = 3), 
                y = rep(c(19.64, 19.64, 18.2, 18.2, 19.64, 19.64, 18.2, 18.2), times = 3), 
                yend = rep(c(20.08, 20.08, 18.65, 18.65, 20.08, 20.08, 18.65, 18.65), times = 3), 
                species = rep(c("A.lyrata", "C.rubella", "E.salsugineum"), each = 4)
            )

            # Adjust position of p-value labels
            mwu_df$label <- paste0(mwu_df$label, c("", "              "))

            y_scale <- c(2.9, 21.125)

            plt_mar <- c(0.1, 1.55, 1.7, 0.55)

            stp_mar <- margin(0.25, 0, 0.25, 0, "cm")

         } else if (species == "AT") { 

          p_mwu <- p_mwu[!grepl("AT", p_mwu$species),]

            # Create df for FDR p-value mapping
            mwu_df <- data.frame(
                class = rep(c("coding_non-core", "coding_core", "lncRNA_non-core", "lncRNA_core"), 
                    times = 1), 
                y = rep(c(18.77, 17.6), times = 2),
                label = ifelse(p_mwu$p_value < 1e-07, "****", 

                    c(paste("italic('P =')~", set_scientific(p_mwu$p_value)))), 

                species = rep(c("A.thaliana"), each = 4)
            )

            # Create df for gem_segments
            h_seg_df <- data.frame(
                x = rep(c(1.105, 2.105, 4.107, 5.107), times = 1), 
                xend = rep(c(3.105, 3.105, 6.107, 6.107), times = 1), 
                y = rep(c(19.08, 17.925, 19.08, 17.925), times = 1), 
                yend = rep(c(19.08, 17.925, 19.08, 17.925), times = 1), 
                species = rep(c("A.thaliana"), each = 4)
            )

            v_seg_df <- data.frame(
                x = rep(c(1.107, 3.107, 2.107, 3.107, 4.107, 6.107, 5.107, 6.107), times = 1), 
                xend = rep(c(1.107, 3.107, 2.107, 3.107, 4.107, 6.107, 5.107, 6.107), times = 1), 
                y = rep(c(18.64, 18.64, 17.45, 17.45, 18.64, 18.64, 17.45, 17.45), times = 1), 
                yend = rep(c(19.08, 19.08, 17.9, 17.9, 19.08, 19.08, 17.9, 17.9), times = 1), 
                species = rep(c("A.thaliana"), each = 4)
            )

            # Adjust position of p-value labels
            mwu_df$label <- paste0(mwu_df$label, c("", "              "))

            y_scale <- c(5.5, 19.89)

            plt_mar <- c(0.1, 32.475, 1.7, 0.55)

            stp_mar <- margin(0.24, 0, 0.26, 0, "cm")

         }

         fname <- sprintf('%s.pdf', paste(deparse(substitute(data)), sep="_"))

         x_lab <- c(Root = "Rt", Hypocotyl = "Hc", Leaf = "Lf", Apex_veg = "Av", 
            Apex_inf = "Ai", Flower = "Fl", Stamen = "St", Carpel = "Ca")

         x_labels = c("coding_all" = expression(atop(NA, atop(textstyle('All'), textstyle('PC')))), 
            "coding_non-core" = expression(atop(NA, atop(textstyle('PC w/o'), textstyle('Ortho')))), 
            "coding_core" = expression(atop(NA, atop(textstyle('Ortho'), textstyle('PC')))), 
            "lncRNA_all" = expression(atop(NA, atop(textstyle('All'), textstyle('lnc')))), 
            "lncRNA_non-core" = expression(atop(NA, atop(textstyle('lnc w/o'), textstyle('Ortho')))), 
            "lncRNA_core" = expression(atop(NA, atop(textstyle('Ortho'), textstyle('lnc')))))

         data$species <- gsub("AT", "A.thaliana", data$species)
         data$species <- gsub("AL", "A.lyrata", data$species)
         data$species <- gsub("CR", "C.rubella", data$species)
         data$species <- gsub("ES", "E.salsugineum", data$species)

         data$class <- factor(data$class, levels = unique(data$class))
         data$conservation <- factor(data$conservation, levels = unique(data$conservation))
         data$species <- factor(data$species, levels = unique(data$species))

         p <- ggplot(data, aes(x = class, y = max_expr, color = class)) + geom_flat_violin(aes(fill = class), colour = "black", position = position_nudge(x = -0.037, y = 0), alpha = 1, size = 0.8) + 
         geom_boxplot(aes(fill = class), colour = "black", width = 0.44, outlier.shape = NA, position = position_hnudge(x = 0.25), size = 0.8, fatten = 2.8, notch = TRUE) +
         scale_x_discrete(expand = c(0.005, 0), labels = x_labels) + 
         scale_y_continuous(limits = y_scale, expand = c(0, 0), breaks = c(5,7.5,10,12.5,15,17.5))
         q <- p + 
         scale_fill_manual(values = c("coding_all" = "#f7ddb0", "coding_non-core" = "#edbb5c", 
            "coding_core" = "#e7a007", "lncRNA_all" = "#cdbee5", "lncRNA_non-core" = "#A689CE", 
            "lncRNA_core" = "#8055b8")) + 
         geom_text(data = mwu_df, mapping = aes(x = class, y = y, label = label), size = 9.275, colour = "black", 
            parse = FALSE, hjust = 0.325, vjust = 0) + 
         geom_segment(data = h_seg_df, mapping = aes(x = x, xend = xend, y = y, yend = yend), size = 0.8, colour = "black") + 
         geom_segment(data = v_seg_df, mapping = aes(x = x, xend = xend, y = y, yend = yend), size = 0.8, colour = "black") + 
         theme_classic() + 
         xlab("") + ylab("Maximum expression    \n (VST-normalized counts)     ") + ggtitle("") + 
         theme(text = element_text(size = 23.5), 
            strip.text = element_text(size = 19.5, face = "italic"), 
                strip.text.x = element_text(margin = stp_mar), 
                strip.background = element_rect(colour = 'white', fill = NA, size = 0.1), 
                axis.ticks.length = unit(0.2, "cm"), 
                axis.ticks = element_line(colour = "black", size = 0.95), 
                axis.line = element_line(colour = 'black', size = 0.95), 
                plot.margin = unit(plt_mar, "cm"), 
                axis.title.y = element_text(size = 18.4, margin = margin(t = 0, r = 6.4, b = 0, l = 3.38), 
                    colour = "black", face = "plain"), 
                axis.title.x = element_text(size = 18.75, margin = margin(t = 6.5, r = 0, b = 5.75, l = 0), 
                    colour = "black", face = "bold"), 
                axis.text.x = element_text(size = 16.25, margin = margin(t = -7, b = 2), colour = "black", 
                    angle = 0, vjust = 1, hjust = 0.5), 
                axis.text.y = element_text(size = 16.5, angle = 0, margin = margin(l = 0, r = 1.5), colour = "black"), 
                panel.spacing = unit(0.7, "cm"), 
                panel.grid.major = element_blank(),
                panel.grid.minor.x = element_blank(), 
                panel.grid.minor.y = element_blank(),  
                legend.position ="none")

         q <- q + facet_wrap(~ factor(species, levels = c("A.thaliana", "A.lyrata", "C.rubella", "E.salsugineum")) , nrow = 1, scales = "free_x")

            ggsave(file = file.path(out_dir, "output", "plots", fname), plot = q, 
                width = 20, height = 5.75, units = c("in"))
      }

      plotMaxExprDist(data = max_expr_dist_non_AT, species = "ACE")
      plotMaxExprDist(data = max_expr_dist_AT, species = "AT")

























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

