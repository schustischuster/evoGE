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



# Define function to get overlapping non-coding / protein-coding genes

getNcPc <- function(species = c("ATH", "AL", "CR", "ES", "TH", "MT", "BD"), 
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

   	# Add an error if threshold < 0
  	if (threshold < 0)
    	stop(
        "'threshold' must be >= 0",
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

	# Save threshold in variable
	threshold <- threshold


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
    # return_list <- list("species_id" = species_id, "GTF" = GTF, "all_genes_tpm" = all_genes_tpm, "threshold" = threshold)
    # return(return_list)
    # }
    # return_objects <- getNcPc("ATH", "single-species", 0.5) # read in GTF and expression data for A.thaliana
    # list2env(return_objects, envir = .GlobalEnv)



    
    #--------- Extract protein-coding genes, seperate them by strand and find overlap ----------


	GTF_df = as.data.frame(GTF)


	# Get all protein-coding genes
	GTF_df_cd_nc <- GTF_df[GTF_df$gene_biotype %in% c("protein_coding","lnc_exonic_antisense","lnc_intronic_antisense"), ]


	# Separate genes from plus strand and minus strand
	strand_plus <- subset(GTF_df_cd_nc, strand == "+")
	strand_minus <- subset(GTF_df_cd_nc, strand == "-")


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


	# Compute overlap length between cd and nc genes
	overlaps <- pintersect(strand_plus_granges[queryHits(overlap_with_strand)], 
						   strand_minus_granges[subjectHits(overlap_with_strand)], 
						   ignore.strand = TRUE)


	widthOverlap <- width(overlaps)
	widthOverlap_df <- as.data.frame(widthOverlap)
	names(widthOverlap_df) <- "width_overlap"


	# Bind "width_overlap" and "percent_overlap" columns to "overlap_with_strand_df"
	overlap_with_strand_df <- cbind(overlap_with_strand_df, widthOverlap_df)



	#--- Extract overlapping protein coding genes from plus_strand/minus_strand data frames ----


	# Add key to plus and minus strand data frames
	strand_plus$key_plus <- seq(1, nrow(strand_plus), 1)
	strand_minus$key_minus <- seq(1, nrow(strand_minus), 1)


	# Generate strand plus and strand minus GTF data frames with overlaping genes based on keys
	strand_plus_overlap_genes <- merge(overlap_with_strand_df, strand_plus, by="key_plus")
	strand_plus_overlap_genes$id  <- seq(1, nrow(strand_plus_overlap_genes), 1)
	strand_plus_overlap_genes = strand_plus_overlap_genes %>% select(
		id,
		gene_id,
		gene_source,
		gene_biotype,
		seqnames, 
		start, 
		end, 
		strand,
		width,
		width_overlap)

	strand_minus_overlap_genes <- merge(overlap_with_strand_df, strand_minus, by="key_minus")
	strand_minus_overlap_genes <- strand_minus_overlap_genes[order(strand_minus_overlap_genes$key_plus),]
	strand_minus_overlap_genes$id  <- seq(1, nrow(strand_minus_overlap_genes), 1)
	strand_minus_overlap_genes = strand_minus_overlap_genes %>% select(
		id,
		gene_id,
		gene_source,
		gene_biotype,
		seqnames, 
		start, 
		end, 
		strand,
		width,
		width_overlap)



	#---------- Merge plus_strand/minus_strand data tables with DevSeq expression data ---------

	## Replace strand_width by "0" if gene_biotype is "protein-coding" - only keep antisense overlap


	strand_plus_overlap_genes_tpm <- merge(strand_plus_overlap_genes, all_genes_tpm, by="gene_id")
	strand_plus_overlap_genes_tpm <- strand_plus_overlap_genes_tpm[order(strand_plus_overlap_genes_tpm$id),]
	strand_plus_overlap_genes_tpm <- dplyr::select(strand_plus_overlap_genes_tpm, -c(
		biotype, source))
	strand_plus_overlap_genes_tpm <- strand_plus_overlap_genes_tpm %>% mutate(
		width_overlap = ifelse(gene_biotype == "protein_coding", 0, width_overlap))

	strand_minus_overlap_genes_tpm <- merge(strand_minus_overlap_genes, all_genes_tpm, by="gene_id")
	strand_minus_overlap_genes_tpm <- strand_minus_overlap_genes_tpm[order(strand_minus_overlap_genes_tpm$id),]
	strand_minus_overlap_genes_tpm <- dplyr::select(strand_minus_overlap_genes_tpm, -c(
		biotype, source))
	strand_minus_overlap_genes_tpm <- strand_minus_overlap_genes_tpm %>% mutate(
		width_overlap = ifelse(gene_biotype == "protein_coding", 0, width_overlap))



	#-------------------------- Get expression table w/o pollen data  --------------------------


	# Remove pollen triplicates from expression table
	strand_plus_overlap_genes_tpm_wo_pollen <- dplyr::select(strand_plus_overlap_genes_tpm, -c(
		flowers_mature_pollen_1, 
		flowers_mature_pollen_2, 
		flowers_mature_pollen_3)) #tibble w/o pollen samles

	strand_minus_overlap_genes_tpm_wo_pollen <- dplyr::select(strand_minus_overlap_genes_tpm, -c(
		flowers_mature_pollen_1, 
		flowers_mature_pollen_2, 
		flowers_mature_pollen_3)) #tibble w/o pollen samles




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
	
			th_replicates <- do.call(cbind, lapply(split.default(df[12:ncol(df)], #adjust columns
								rep(seq_along(df), each = 3, length.out = ncol(df)-11)), #adjust columns
								function(x) {
									x[rowSums(x > threshold) < 2, ] <- 0; 
									x
								}
							))

			# Bind key/gene_id/id/gene_source/gene_biotype/seqnames/start/end/strand/width/width_overlap columns to thresholded data frame
			th_replicates <- cbind(df[1:11], th_replicates)

			# Remove all rows that only contain "0"
			th_replicates <- th_replicates[which(rowSums(th_replicates[,-1:-11, drop = FALSE] > 0) > 0),]

			return(th_replicates)
		}

		# Apply threshold to data and extract keys ("key")
		keys_data <- getThreshold(df)
		keys_data <- keys_data[,1:2]
		names(keys_data) <- c("key","ID")

		# Generate thresholded data frame based on keys
		th_df <- merge(keys_data, df, by="key")
		th_df <- th_df[-1:-2]

		return(th_df)
	}


	# Apply threshold function for expression tables w/ pollen
	strand_plus_overlap_genes_tpm_th <- applyThreshold(strand_plus_overlap_genes_tpm, threshold)
	strand_minus_overlap_genes_tpm_th <- applyThreshold(strand_minus_overlap_genes_tpm, threshold)


	strand_plus_overlap_genes_tpm_th <- strand_plus_overlap_genes_tpm_th[(
		strand_plus_overlap_genes_tpm_th$id %in% strand_minus_overlap_genes_tpm_th$id),]
	strand_minus_overlap_genes_tpm_th <- strand_minus_overlap_genes_tpm_th[(
		strand_minus_overlap_genes_tpm_th$id %in% strand_plus_overlap_genes_tpm_th$id),]


	# Apply threshold function for expression tables w/o pollen
	strand_plus_overlap_genes_tpm_th_wo_pollen <- applyThreshold(strand_plus_overlap_genes_tpm_wo_pollen, threshold)
	strand_minus_overlap_genes_tpm_th_wo_pollen <- applyThreshold(strand_minus_overlap_genes_tpm_wo_pollen, threshold)


	strand_plus_overlap_genes_tpm_th_wo_pollen <- strand_plus_overlap_genes_tpm_th_wo_pollen[(
		strand_plus_overlap_genes_tpm_th_wo_pollen$id %in% strand_minus_overlap_genes_tpm_th_wo_pollen$id),]
	strand_minus_overlap_genes_tpm_th_wo_pollen <- strand_minus_overlap_genes_tpm_th_wo_pollen[(
		strand_minus_overlap_genes_tpm_th_wo_pollen$id %in% strand_plus_overlap_genes_tpm_th_wo_pollen$id),]




	#------------------- Compute correlation between coding gene and cisNAT  -------------------


	# Compute spearman and pearson correlations
	getCor <- function(df1, df2) {

		df1_col <- ncol(df1)
		df2_col <- ncol(df2)

		# startup message
		message("Computing correlation...")

		df1$Spearman <- sapply(1:nrow(df1), function(i) 
	    	cor(as.numeric(df1[i, 11:df1_col]), as.numeric(df2[i, 11:df2_col]), method=c("spearman")))

		df1$Pearson <- sapply(1:nrow(df1), function(i) 
	    	cor(as.numeric(df1[i, 11:df1_col]), as.numeric(df2[i, 11:df2_col]), method=c("pearson")))

		return(df1)
	}


	strand_plus_overlap_genes <- getCor(
		strand_plus_overlap_genes_tpm_th, strand_minus_overlap_genes_tpm_th)
	strand_plus_overlap_genes_wo_pollen <- getCor(
		strand_plus_overlap_genes_tpm_th_wo_pollen, strand_minus_overlap_genes_tpm_th_wo_pollen)




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




#----- Create ATH_all w/o pollen data tables containing expression values and write csv file ------


	# Create data table containing both strand plus and minus genes, their expression data and cor values
	if ((is.element("ATH", species)) && (is.element("single-species", experiment))) {

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

		overlap_cd_nc_genes_cor_wo_pollen_query_expr	<- dplyr::select(
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



# Execute getNcPc function
getNcPc("ATH", "single-species", 0.5)
getNcPc("ATH", "comparative", 0.5)
getNcPc("AL", "single-species", 0.5)
getNcPc("AL", "comparative", 0.5)
getNcPc("CR", "comparative", 0.5)
getNcPc("ES", "comparative", 0.5)
getNcPc("TH", "comparative", 0.5)
getNcPc("MT", "comparative", 0.5)
getNcPc("BD", "comparative", 0.5)

getNcPc("ATH", "single-species", 2)
getNcPc("ATH", "comparative", 2)
getNcPc("AL", "single-species", 2)
getNcPc("AL", "comparative", 2)
getNcPc("CR", "comparative", 2)
getNcPc("ES", "comparative", 2)
getNcPc("TH", "comparative", 2)
getNcPc("MT", "comparative", 2)
getNcPc("BD", "comparative", 2)

getNcPc("ATH", "single-species", 5)
getNcPc("ATH", "comparative", 5)
getNcPc("AL", "single-species", 5)
getNcPc("AL", "comparative", 5)
getNcPc("CR", "comparative", 5)
getNcPc("ES", "comparative", 5)
getNcPc("TH", "comparative", 5)
getNcPc("MT", "comparative", 5)
getNcPc("BD", "comparative", 5)

getNcPc("ATH", "single-species", 0)

