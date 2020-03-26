# Find non-overlapping neighboring protein-coding genes
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
if (!require(plyr)) install.packages('plyr')
library(plyr)
if (!require(dplyr)) install.packages('dplyr')
library(dplyr)
if (!require(GenomicRanges)) install.packages('GenomicRanges')
library(GenomicRanges)
if (!require(rtracklayer)) install.packages('rtracklayer')
library(rtracklayer)


# Set file path and input files
in_dir <- "/Volumes/User/Shared/Christoph_manuscript/DevSeq_paper/Analysis/Analysis_2019/A_thaliana_gene_exression_map/20191121_CS_coding_cisNAT_analysis/data"
out_dir <- "/Volumes/User/Shared/Christoph_manuscript/DevSeq_paper/Analysis/Analysis_2019/A_thaliana_gene_exression_map/20191121_CS_coding_cisNAT_analysis"



# Define function to get overlapping protein-coding genes

getPcPcNO <- function(species = c("ATH", "AL", "CR", "ES", "TH", "MT", "BD"), 
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


    # Stop function here to allow specific analysis of a single species
    # return_list <- list("species_id" = species_id, "GTF" = GTF, "all_genes_tpm" = all_genes_tpm)
    # return(return_list)
    # }
    # return_objects <- getPcPc("ATH", "single-species") # read in GTF and expression data for A.thaliana
    # list2env(return_objects, envir = .GlobalEnv)



    
    #--------- Extract protein-coding genes, seperate them by strand and find overlap ----------


	GTF_df = as.data.frame(GTF)

	# Get all protein-coding genes
	GTF_df_cd <- subset(GTF_df, gene_biotype == "protein_coding")

	# Remove all chloroplast and mito genes
	# Reason: RNA-seq prep kit used is Ribo_zero kit, which incompletely removes Pt and
	# Mt transcripts, so expression estimates of those genes is likely to be wrong
	GTF_df_cd <- subset(GTF_df_cd, seqnames != "Pt")
	GTF_df_cd <- subset(GTF_df_cd, seqnames != "Mt")

	# Find nearest protein-coding gene for plus and minus strands
	GTF_df_cd_GR <- makeGRangesFromDataFrame(GTF_df_cd, keep.extra.columns=FALSE, 
		ignore.strand=FALSE, seqinfo=NULL, seqnames.field=c(
		"seqnames", "seqname","chromosome", "chrom","chr", "chromosome_name","seqid"),
		start.field="start", end.field=c("end", "stop"),starts.in.df.are.0based=FALSE)

	# Get nearest protein-coding genes and their distance to each other
	# a distance <= 0 indicates overlapping genes on same or opposite strand -> remove them afterwards!
	transcript_pairs <- distanceToNearest(GTF_df_cd_GR, select=c("arbitrary"), ignore.strand=TRUE)
	transcript_pairs_df <- as.data.frame(transcript_pairs)
	transcript_pairs_df <- subset(transcript_pairs_df, distance != 0)



	
	#------- Extract non-overlapping protein coding gene pair descriptions from GTF data --------


	# Add key to plus_minus strand data frames
	queryHits <- seq(1, nrow(GTF_df_cd), 1)
	subjectHits <- seq(1, nrow(GTF_df_cd), 1)
	strand_plus_minus_query <- cbind(as.data.frame(queryHits), GTF_df_cd)
	strand_plus_minus_subject <- cbind(as.data.frame(subjectHits), GTF_df_cd)


	strand_plus_minus_query_genes <- join(transcript_pairs_df, strand_plus_minus_query, by="queryHits")
	strand_plus_minus_query_genes$id  <- seq(1, nrow(strand_plus_minus_query_genes), 1)
	strand_plus_minus_query_genes = strand_plus_minus_query_genes %>% select(
		id,
		subjectHits,
		queryHits,
		distance,
		gene_id,
		seqnames, 
		start, 
		end, 
		strand)

	strand_plus_minus_subject_genes <- join(transcript_pairs_df, strand_plus_minus_subject, by="subjectHits")
	strand_plus_minus_subject_genes$id  <- seq(1, nrow(strand_plus_minus_subject_genes), 1)
	strand_plus_minus_subject_genes = strand_plus_minus_subject_genes %>% select(
		id,
		queryHits,
		subjectHits,
		distance,
		gene_id,
		seqnames, 
		start, 
		end, 
		strand)




	#---------- Merge plus_strand/minus_strand data tables with DevSeq expression data ---------

	strand_plus_minus_query_genes_tpm <- join(strand_plus_minus_query_genes, all_genes_tpm, by="gene_id")
	strand_plus_minus_subject_genes_tpm <- join(strand_plus_minus_subject_genes, all_genes_tpm, by="gene_id")




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
	
			th_replicates <- do.call(cbind, lapply(split.default(df[13:ncol(df)], #adjust columns
								rep(seq_along(df), each = 3, length.out = ncol(df)-12)), #adjust columns
								function(x) {
									x[rowSums(x > threshold) < 2, ] <- 0; 
									x
								}
							))

			# Bind key/id/prt_id/symbol/biotype/source columns to thresholded data frame
			th_replicates <- cbind(df[1:12], th_replicates)

			# Remove all rows that only contain "0"
			th_replicates <- th_replicates[which(rowSums(th_replicates[,-1:-12, drop = FALSE] > 0) > 0),]

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
	strand_plus_minus_query_genes_tpm_0.5 <- applyThreshold(strand_plus_minus_query_genes_tpm,0.5)
	strand_plus_minus_subject_genes_tpm_0.5 <- applyThreshold(strand_plus_minus_subject_genes_tpm,0.5)

	strand_plus_minus_query_genes_tpm_0.5 <- strand_plus_minus_query_genes_tpm_0.5[(
		strand_plus_minus_query_genes_tpm_0.5$id %in% strand_plus_minus_subject_genes_tpm_0.5$id),]
	strand_plus_minus_subject_genes_tpm_0.5 <- strand_plus_minus_subject_genes_tpm_0.5[(
		strand_plus_minus_subject_genes_tpm_0.5$id %in% strand_plus_minus_query_genes_tpm_0.5$id),]



	#---- Split data to sense-sense (same_strand PCT pairs) and sense-antisense (SAS PCT pairs)  ----

	PCT_pairs_query_plus <- subset(strand_plus_minus_query_genes_tpm_0.5, strand == "+")
	PCT_pairs_subject_plus <- subset(strand_plus_minus_subject_genes_tpm_0.5, strand == "+")
	PCT_pairs_query_minus <- subset(strand_plus_minus_query_genes_tpm_0.5, strand == "-")
	PCT_pairs_subject_minus <- subset(strand_plus_minus_subject_genes_tpm_0.5, strand == "-")

	same_strand_PCT_pairs_query_plus <- PCT_pairs_query_plus[(
		PCT_pairs_query_plus$id %in% PCT_pairs_subject_plus$id),]
	same_strand_PCT_pairs_subject_plus <- PCT_pairs_subject_plus[(
		PCT_pairs_subject_plus$id %in% PCT_pairs_query_plus$id),]
	same_strand_PCT_pairs_query_minus <- PCT_pairs_query_minus[(
		PCT_pairs_query_minus$id %in% PCT_pairs_subject_minus$id),]
	same_strand_PCT_pairs_subject_minus <- PCT_pairs_subject_minus[(
		PCT_pairs_subject_minus$id %in% PCT_pairs_query_minus$id),]

	same_strand_PCT_pairs_query <- rbind(same_strand_PCT_pairs_query_plus, same_strand_PCT_pairs_query_minus)
	same_strand_PCT_pairs_subject <- rbind(same_strand_PCT_pairs_subject_plus, same_strand_PCT_pairs_subject_minus)

	SAS_PCT_pairs_query_plus <- PCT_pairs_query_plus[(
		PCT_pairs_query_plus$id %in% PCT_pairs_subject_minus$id),]
	SAS_PCT_pairs_subject_plus <- PCT_pairs_subject_plus[(
		PCT_pairs_subject_plus$id %in% PCT_pairs_query_minus$id),]
	SAS_PCT_pairs_query_minus <- PCT_pairs_query_minus[(
		PCT_pairs_query_minus$id %in% PCT_pairs_subject_plus$id),]
	SAS_PCT_pairs_subject_minus <- PCT_pairs_subject_minus[(
		PCT_pairs_subject_minus$id %in% PCT_pairs_query_plus$id),]

	SAS_PCT_pairs_query <- rbind(SAS_PCT_pairs_query_plus, SAS_PCT_pairs_query_minus)
	SAS_PCT_pairs_subject <- rbind(SAS_PCT_pairs_subject_minus, SAS_PCT_pairs_subject_plus)



	#-------------------------- Get expression table w/o pollen data  --------------------------


	# Remove pollen triplicates from expression table
	same_strand_PCT_pairs_query_wo_pollen <- dplyr::select(same_strand_PCT_pairs_query, -c(
		flowers_mature_pollen_1, 
		flowers_mature_pollen_2, 
		flowers_mature_pollen_3)) #tibble w/o pollen samles

	same_strand_PCT_pairs_subject_wo_pollen <- dplyr::select(same_strand_PCT_pairs_subject, -c(
		flowers_mature_pollen_1, 
		flowers_mature_pollen_2, 
		flowers_mature_pollen_3)) #tibble w/o pollen samles

	SAS_PCT_pairs_query_wo_pollen <- dplyr::select(SAS_PCT_pairs_query, -c(
		flowers_mature_pollen_1, 
		flowers_mature_pollen_2, 
		flowers_mature_pollen_3)) #tibble w/o pollen samles

	SAS_PCT_pairs_subject_wo_pollen <- dplyr::select(SAS_PCT_pairs_subject, -c(
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
	    	cor(as.numeric(df1[i, 12:df1_col]), as.numeric(df2[i, 12:df2_col]), method=c("spearman")))

		# log2-transform TPM values before computing Pearson
		df1[, 12:df1_col] <- log2(df1[, 12:df1_col] + 1)
		df2[, 12:df1_col] <- log2(df2[, 12:df2_col] + 1)

		df1$Pearson <- sapply(1:nrow(df1), function(i) 
	    	cor(as.numeric(df1[i, 12:df1_col]), as.numeric(df2[i, 12:df2_col]), method=c("pearson")))

		return(df1)
	}


	same_strand_PCT_pairs_wo_pollen <- getCor(
		same_strand_PCT_pairs_query_wo_pollen, same_strand_PCT_pairs_subject_wo_pollen)
	SAS_PCT_pairs_wo_pollen <- getCor(
		SAS_PCT_pairs_query_wo_pollen, SAS_PCT_pairs_subject_wo_pollen)




#------------------------- Create final data table and write csv file  -------------------------


	# Create data table containing both strand plus and minus genes and cor values
	same_strand_PCT_pairs_subject_descript = same_strand_PCT_pairs_subject_wo_pollen %>% select(id, gene_id, start, end, strand, biotype)
	names(same_strand_PCT_pairs_subject_descript) <- c("key_subject" ,"id_subject", "start_subject", "end_subject", "strand_subject", "biotype_subject")
	same_strand_PCT_pairs_wo_pollen_descript <- same_strand_PCT_pairs_wo_pollen
	names(same_strand_PCT_pairs_wo_pollen_descript)[1:10] <- c("key_query", "subjectHits", "queryHits", "distance", "id_query", "seqnames", "start_query", "end_query", "strand_query", "biotype_query")
	same_strand_PCT_pairs <- cbind(same_strand_PCT_pairs_wo_pollen_descript, same_strand_PCT_pairs_subject_descript)
	same_strand_PCT_pairs = same_strand_PCT_pairs %>% select(
		id_query, 
		start_query, 
		end_query, 
		strand_query, 
		biotype_query, 
		id_subject, 
		start_subject, 
		end_subject, 
		strand_subject, 
		biotype_subject, 
		distance, 
		seqnames, 
		queryHits, 
		subjectHits, 
		Spearman, 
		Pearson)


	SAS_PCT_pairs_subject_descript = SAS_PCT_pairs_subject_wo_pollen %>% select(id, gene_id, start, end, strand, biotype)
	names(SAS_PCT_pairs_subject_descript) <- c("key_subject" ,"id_subject", "start_subject", "end_subject", "strand_subject", "biotype_subject")
	SAS_PCT_pairs_wo_pollen_descript <- SAS_PCT_pairs_wo_pollen
	names(SAS_PCT_pairs_wo_pollen_descript)[1:10] <- c("key_query", "subjectHits", "queryHits", "distance", "id_query", "seqnames", "start_query", "end_query", "strand_query", "biotype_query")
	SAS_PCT_pairs <- cbind(SAS_PCT_pairs_wo_pollen_descript, SAS_PCT_pairs_subject_descript)
	SAS_PCT_pairs = SAS_PCT_pairs %>% select(
		id_query, 
		start_query, 
		end_query, 
		strand_query, 
		biotype_query, 
		id_subject, 
		start_subject, 
		end_subject, 
		strand_subject, 
		biotype_subject, 
		distance, 
		seqnames, 
		queryHits, 
		subjectHits, 
		Spearman, 
		Pearson)


	# Set filename
    same_strand_fname <- sprintf('%s.csv', paste(species_id, "same_strand_PCT_pairs", sep="_"))
	SAS_fname <- sprintf('%s.csv', paste(species_id, "SAS_PCT_pairs", sep="_"))


	# Write final data tables to csv files and store them in /out_dir/output/data_tables
	if (!dir.exists(file.path(out_dir, "output", "cd_gene_pairs"))) 
		dir.create(file.path(out_dir, "output", "cd_gene_pairs"), recursive = TRUE)

	write.table(same_strand_PCT_pairs, file=file.path(out_dir, "output", "cd_gene_pairs", same_strand_fname), 
		sep=";", dec=".", row.names=FALSE, col.names=TRUE)
	write.table(SAS_PCT_pairs, file=file.path(out_dir, "output", "cd_gene_pairs", SAS_fname), 
		sep=";", dec=".", row.names=FALSE, col.names=TRUE)

}



# Execute getPcPc function
getPcPcNO("ATH", "single-species")
getPcPcNO("ATH", "comparative")
getPcPcNO("AL", "single-species")
getPcPcNO("AL", "comparative")
getPcPcNO("CR", "comparative")
getPcPcNO("ES", "comparative")
getPcPcNO("TH", "comparative")
getPcPcNO("MT", "comparative")
getPcPcNO("BD", "comparative")




