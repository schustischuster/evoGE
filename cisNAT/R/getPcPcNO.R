# Find non-overlapping neighbouring protein-coding genes
# Data input: 1) GTF file | 2) Expression_data WITHOUT mito and chloroplast genes
# 3) List of orthologous protein-coding genes for each species (defined as orthogroups)
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


#----------------------------------------- Read data -----------------------------------------


getPcPcNO <- function(species = c("AT", "AL", "CR", "ES", "TH", "MT", "BD"), 
	experiment = c("single-species", "comparative"), threshold) {
	
	# Show error message if no species is chosen
    if (missing(species))
   
       stop(
       "Please choose one of the available species: 
	   'AT', 'AL', 'CR', 'ES', 'TH', 'MT', 'BD'",
	   call. = TRUE
       )

   	# Show error message for AT and AL if no experiment is chosen
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

   	# Add an error if threshold < 0
  	if (threshold < 0)
    	stop(
        "'threshold' must be >= 0",
	   	call. = TRUE
    	)


	# Set GTF input gtf file
    if (is.element("AT", species)) {
    	GTFfile = file.path(in_dir, "GTF", "AT_final_annotation.gtf")
        genesCounts = file.path(in_dir, "Expression_data", "AT_genes_inter_norm_count_mat_vsd_sample_names.csv")
        genesTPM = file.path(in_dir, "Expression_data", "AT_genes_inter_norm_tpm_mat_deseq_sample_names.csv")
        species_id <- "AT"

    } else if (is.element("AL", species)) {
		GTFfile = file.path(in_dir, "GTF", "AL_final_annotation.gtf")
		genesCounts = file.path(in_dir, "Expression_data", "AL_genes_inter_norm_count_mat_vsd_sample_names.csv")
		genesTPM = file.path(in_dir, "Expression_data", "AL_genes_inter_norm_tpm_mat_deseq_sample_names.csv")
		species_id <- "AL"

    } else if (is.element("CR", species)) {
		GTFfile = file.path(in_dir, "GTF", "CR_final_annotation.gtf")
		genesCounts = file.path(in_dir, "Expression_data", "CR_genes_inter_norm_count_mat_vsd_sample_names.csv")
		genesTPM = file.path(in_dir, "Expression_data", "CR_genes_inter_norm_tpm_mat_deseq_sample_names.csv")
		species_id <- "CR"

    } else if (is.element("ES", species)) {
		GTFfile = file.path(in_dir, "GTF", "ES_final_annotation.gtf")
		genesCounts = file.path(in_dir, "Expression_data", "ES_genes_inter_norm_count_mat_vsd_sample_names.csv")
		genesTPM = file.path(in_dir, "Expression_data", "ES_genes_inter_norm_tpm_mat_deseq_sample_names.csv")
		species_id <- "ES"

    } else if (is.element("TH", species)) {
		GTFfile = file.path(in_dir, "GTF", "TH_final_annotation.gtf")
		genesCounts = file.path(in_dir, "Expression_data", "TH_genes_inter_norm_count_mat_vsd_sample_names.csv")
		genesTPM = file.path(in_dir, "Expression_data", "TH_genes_inter_norm_tpm_mat_deseq_sample_names.csv")
		species_id <- "TH"

    } else if (is.element("MT", species)) {
		GTFfile = file.path(in_dir, "GTF", "MT_final_annotation.gtf")
		genesCounts = file.path(in_dir, "Expression_data", "MT_genes_inter_norm_count_mat_vsd_sample_names.csv")
		genesTPM = file.path(in_dir, "Expression_data", "MT_genes_inter_norm_tpm_mat_deseq_sample_names.csv")
		species_id <- "MT"

    } else if (is.element("BD", species)) {
		GTFfile = file.path(in_dir, "GTF", "BD_final_annotation.gtf")
		genesCounts = file.path(in_dir, "Expression_data", "BD_genes_inter_norm_count_mat_vsd_sample_names.csv")
		genesTPM = file.path(in_dir, "Expression_data", "BD_genes_inter_norm_tpm_mat_deseq_sample_names.csv")
		species_id <- "BD"
    }


	# Import gtf file
	GTF = import.gff(GTFfile, format = "gtf", feature.type = "gene")

	# Read expression data
	all_genes_counts <- read.table(genesCounts, sep = ";", dec = ".", header = TRUE, stringsAsFactors = FALSE)
	all_genes_tpm <- read.table(genesTPM, sep = ";", dec = ".", header = TRUE, stringsAsFactors = FALSE)

	# Save threshold in variable
	threshold <- threshold


	# Format expression data and rename pollen samples
    if ((is.element("AT", species)) && (is.element("comparative", experiment))) {

		all_genes_counts <- dplyr::select(all_genes_counts, c(
			root_root_tip_5d_.1.,
			root_root_tip_5d_.2.,
			root_root_tip_5d_.3.,
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
			flower_early_stg12_carpels_21d._.1.,
			flower_early_stg12_carpels_21d._.2.,
			flower_early_stg12_carpels_21d._.3.)) #tibble w/o pollen samles


		all_genes_tpm <- dplyr::select(all_genes_tpm, c(
			root_root_tip_5d_.1.,
			root_root_tip_5d_.2.,
			root_root_tip_5d_.3.,
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
			flower_early_stg12_carpels_21d._.1.,
			flower_early_stg12_carpels_21d._.2.,
			flower_early_stg12_carpels_21d._.3.)) #tibble w/o pollen samles


		species_id <- "AT_comp"


    } else if ((is.element("AL", species)) && (is.element("comparative", experiment))) {

		all_genes_counts <- dplyr::select(all_genes_counts, -c(
			flower_stg11_stamens_8w.10w.25d_.1., 
			flower_stg11_stamens_8w.10w.25d_.2., 
			flower_stg11_stamens_8w.10w.25d_.3.,
			flower_early_stg12_stamens_8w.10w.23d_.1.,
			flower_early_stg12_stamens_8w.10w.23d_.2.,
			flower_early_stg12_stamens_8w.10w.23d_.3.,
			flower_late_stg12_stamens_8w.10w.21d_.1.,
			flower_late_stg12_stamens_8w.10w.21d_.2.,
			flower_late_stg12_stamens_8w.10w.21d_.3.)) #tibble w/o pollen samles

		all_genes_tpm <- dplyr::select(all_genes_tpm, -c(
			flower_stg11_stamens_8w.10w.25d_.1., 
			flower_stg11_stamens_8w.10w.25d_.2., 
			flower_stg11_stamens_8w.10w.25d_.3.,
			flower_early_stg12_stamens_8w.10w.23d_.1.,
			flower_early_stg12_stamens_8w.10w.23d_.2.,
			flower_early_stg12_stamens_8w.10w.23d_.3.,
			flower_late_stg12_stamens_8w.10w.21d_.1.,
			flower_late_stg12_stamens_8w.10w.21d_.2.,
			flower_late_stg12_stamens_8w.10w.21d_.3.)) #tibble w/o pollen samles


		species_id <- "AL_comp"


    }



    # return_list <- list("species_id" = species_id, "GTF" = GTF, "all_genes_counts" = all_genes_counts, "all_genes_tpm" = all_genes_tpm, "threshold" = threshold)
    # return(return_list)
    # }
    # return_objects <- getPcPcNO("AT", "single-species", 0.5) # read in GTF and expression data for A.thaliana
    # list2env(return_objects, envir = .GlobalEnv)



    
    #--------- Extract protein-coding genes, seperate them by strand and find overlap ----------


	GTF_df = as.data.frame(GTF)


	# Get all protein-coding genes
	GTF_df_cd <- subset(GTF_df, gene_biotype == "protein_coding")

	# Remove all chloroplast and mito genes
	# Reason: RNA-seq prep kit used is Ribo_zero, which incompletely removes Pt and
	# Mt transcripts, so expression estimates of those genes are likely to be wrong
	GTF_df_cd <- subset(GTF_df_cd, seqnames != "Pt")
	GTF_df_cd <- subset(GTF_df_cd, seqnames != "Mt")

	# Get positional, chromosome and strand informarmation for protein-coding genes
	GTF_df_cd_GR <- makeGRangesFromDataFrame(GTF_df_cd, keep.extra.columns = FALSE, 
		ignore.strand = FALSE, seqinfo = NULL, seqnames.field = "seqnames",
		start.field = "start", end.field = "end", starts.in.df.are.0based = FALSE)

	# Get nearest protein-coding genes and their distance to each other
	# a distance <= 0 indicates overlapping genes on same or opposite strand -> remove them afterwards!
	transcript_pairs <- distanceToNearest(GTF_df_cd_GR, select = "arbitrary", ignore.strand = TRUE)
	transcript_pairs_df <- as.data.frame(transcript_pairs)
	transcript_pairs_df <- subset(transcript_pairs_df, distance != 0)



	#------- Extract non-overlapping protein coding gene pair descriptions from GTF data --------


	# Add key to plus_minus strand data frames
	queryHits <- seq(1, nrow(GTF_df_cd), 1)
	subjectHits <- seq(1, nrow(GTF_df_cd), 1)
	strand_plus_minus_query <- cbind(as.data.frame(queryHits), GTF_df_cd)
	strand_plus_minus_subject <- cbind(as.data.frame(subjectHits), GTF_df_cd)


	strand_plus_minus_query_genes <- join(transcript_pairs_df, strand_plus_minus_query, by = "queryHits")
	strand_plus_minus_query_genes$id  <- seq(1, nrow(strand_plus_minus_query_genes), 1)
	strand_plus_minus_query_genes = strand_plus_minus_query_genes %>% select(
		id,
		queryHits,
		subjectHits,
		distance,
		gene_id,
		seqnames, 
		start, 
		end, 
		strand)

	strand_plus_minus_subject_genes <- join(transcript_pairs_df, strand_plus_minus_subject, by = "subjectHits")
	strand_plus_minus_subject_genes$id  <- seq(1, nrow(strand_plus_minus_subject_genes), 1)
	strand_plus_minus_subject_genes = strand_plus_minus_subject_genes %>% select(
		id,
		subjectHits,
		queryHits,
		distance,
		gene_id,
		seqnames, 
		start, 
		end, 
		strand)



	#---------- Merge plus_strand/minus_strand data tables with DevSeq expression data ---------


	all_genes_counts <- tibble::rownames_to_column(all_genes_counts, "gene_id")
	all_genes_tpm <- tibble::rownames_to_column(all_genes_tpm, "gene_id")


	strand_plus_minus_query_genes_counts <- merge(strand_plus_minus_query_genes, all_genes_counts, by = "gene_id")
	strand_plus_minus_query_genes_counts <- strand_plus_minus_query_genes_counts[order(strand_plus_minus_query_genes_counts$id),]

	strand_plus_minus_subject_genes_counts <- merge(strand_plus_minus_subject_genes, all_genes_counts, by = "gene_id")
	strand_plus_minus_subject_genes_counts <- strand_plus_minus_subject_genes_counts[order(strand_plus_minus_subject_genes_counts$id),]



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

		# Define threshold function
		getThreshold <- function(df) {

			# Split data frame by sample replicates into a list then apply threshold for each subset
	
			th_replicates <- do.call(cbind, lapply(split.default(df[2:ncol(df)], #adjust columns
								rep(seq_along(df), each = 3, length.out = ncol(df)-1)), #adjust columns
								function(x) {
									x[rowSums(x > threshold) < 2, ] <- 0; 
									x
								}
							))

			# Get max and mean expression for each gene
			th_replicates$max <- apply(X = th_replicates, MARGIN = 1, FUN = max)
			th_replicates$avg <- apply(X = th_replicates, MARGIN = 1, FUN = mean)

			# Bind key and gene_id columns to thresholded data frame
			th_replicates <- cbind(df[1], th_replicates)

			# Reorder df
			th_replicates <- cbind(th_replicates[c("gene_id", "max", "avg")], 
				th_replicates[2:(ncol(th_replicates)-2)])

			# Remove all rows that only contain "0"
			th_replicates <- th_replicates[which(rowSums(th_replicates[,-1:-3, drop = FALSE] > 0) > 0),]

			return(th_replicates)
		}

		# Apply threshold to data and extract keys ("key")
		data <- getThreshold(df)
		th_df <- data[1:3]

		return(th_df)
	}


	# Apply threshold function for expression tables w/ pollen
	genes_tpm_th <- applyThreshold(all_genes_tpm, threshold)


	# Remove all gene pairs that show expression below threshold
	strand_plus_minus_query_genes_counts_th <- merge(
		genes_tpm_th, strand_plus_minus_query_genes_counts, by = "gene_id")
	strand_plus_minus_subject_genes_counts_th <- merge(
		genes_tpm_th, strand_plus_minus_subject_genes_counts, by = "gene_id")


	# Remove all gene pairs where only one partner shows expression above threshold
	strand_plus_minus_query_genes_counts_th <- strand_plus_minus_query_genes_counts_th[(
		strand_plus_minus_query_genes_counts_th$id %in% strand_plus_minus_subject_genes_counts_th$id),]
	strand_plus_minus_subject_genes_counts_th <- strand_plus_minus_subject_genes_counts_th[(
		strand_plus_minus_subject_genes_counts_th$id %in% strand_plus_minus_query_genes_counts_th$id),]

	strand_plus_minus_query_genes_counts_th <- strand_plus_minus_query_genes_counts_th[order(strand_plus_minus_query_genes_counts_th$id),]
	strand_plus_minus_subject_genes_counts_th <- strand_plus_minus_subject_genes_counts_th[order(strand_plus_minus_subject_genes_counts_th$id),]



	#------------------- Compute correlation between coding gene and cisNAT  -------------------


	# Compute spearman and pearson correlations
	getCor <- function(df1, df2) {

		df1_col <- ncol(df1)
		df2_col <- ncol(df2)

		# startup message
		message("Computing correlation...")

		df1$Spearman <- sapply(1:nrow(df1), function(i) 
	    	cor(as.numeric(df1[i, 12:df1_col]), as.numeric(df2[i, 12:df2_col]), method = c("spearman")))

		# no log-transformation of values required before computing Pearson (VST counts)

		df1$Pearson <- sapply(1:nrow(df1), function(i) 
	    	cor(as.numeric(df1[i, 12:df1_col]), as.numeric(df2[i, 12:df2_col]), method = c("pearson")))

		return(df1)
	}


	strand_plus_minus_query_genes_cor <- getCor(
		strand_plus_minus_query_genes_counts_th, strand_plus_minus_subject_genes_counts_th)
	strand_plus_minus_subject_genes_cor <- getCor(
		strand_plus_minus_subject_genes_counts_th, strand_plus_minus_query_genes_counts_th)


	# Create a single table containing protein-coding/protein-coding and protein-coding/cisNAT pairs
	all_neighbouring_gene_pairs <- data.frame(
		id1 = strand_plus_minus_query_genes_cor$id,
		gene_id1 = strand_plus_minus_query_genes_cor$gene_id,
		seqnames1 = strand_plus_minus_query_genes_cor$seqnames,
		start1 = strand_plus_minus_query_genes_cor$start,
		end1 = strand_plus_minus_query_genes_cor$end,
		strand1 = strand_plus_minus_query_genes_cor$strand,
		max_expr1 = strand_plus_minus_query_genes_cor$max,
		mean_expr1 = strand_plus_minus_query_genes_cor$avg,
		id2 = strand_plus_minus_subject_genes_cor$id,
		gene_id2 = strand_plus_minus_subject_genes_cor$gene_id,
		seqnames2 = strand_plus_minus_subject_genes_cor$seqnames,
		start2 = strand_plus_minus_subject_genes_cor$start,
		end2 = strand_plus_minus_subject_genes_cor$end,
		strand2 = strand_plus_minus_subject_genes_cor$strand,
		max_expr2 = strand_plus_minus_subject_genes_cor$max,
		mean_expr2 = strand_plus_minus_subject_genes_cor$avg,
		distance = strand_plus_minus_query_genes_cor$distance,
		Spearman = strand_plus_minus_query_genes_cor$Spearman,
		Pearson = strand_plus_minus_query_genes_cor$Pearson)



	# Remove one copy of duplicated gene pairs (based on distance and Pearson cor)
	dup_pairs <- all_neighbouring_gene_pairs[duplicated(all_neighbouring_gene_pairs[,c(17,19)]),]
	all_neighbouring_dedup_gene_pairs <- all_neighbouring_gene_pairs[!all_neighbouring_gene_pairs$id1 %in% dup_pairs$id1, ]
	# 10932 unique gene pairs in AT



#--------------------------------------- Write csv file ---------------------------------------


	# Set filename
    fname <- sprintf('%s.csv', paste(species_id, "cd_cd_NO_cor", threshold, sep = "_"))


	# Write final data tables to csv files and store them in /out_dir/output/data_tables
	if (!dir.exists(file.path(out_dir, "output", "cd_gene_pairs"))) 
		dir.create(file.path(out_dir, "output", "cd_gene_pairs"), recursive = TRUE)
	message("Storing results in: ", file.path("output", "cd_gene_pairs"))

	write.table(all_neighbouring_dedup_gene_pairs, file = file.path(out_dir, "output", "cd_gene_pairs", fname), 
		sep = ";", dec = ".", row.names = FALSE, col.names = TRUE)


}


