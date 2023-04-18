# Find overlapping protein-coding gene pairs
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



#----------------------------------------- Read data -----------------------------------------


getPcPc <- function(species = c("AT", "AL", "CR", "ES", "TH", "MT", "BD"), 
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
	GTF = import.gff(GTFfile, format="gtf", feature.type="gene")

	# Read expression data
	all_genes_counts <- read.table(genesCounts, sep=";", dec=".", header=TRUE, stringsAsFactors=FALSE)
	all_genes_tpm <- read.table(genesTPM, sep=";", dec=".", header=TRUE, stringsAsFactors=FALSE)

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


		species_id <- "AT_comparative_samples"


    } else if ((is.element("AL", species)) && (is.element("comparative", experiment))) {

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


		species_id <- "AL_comparative_samples"


    }



    # return_list <- list("species_id" = species_id, "GTF" = GTF, "all_genes_counts" = all_genes_counts, "all_genes_tpm" = all_genes_tpm, "threshold" = threshold)
    # return(return_list)
    # }
    # return_objects <- getPcPc("AT", "single-species", 0.5) # read in GTF and expression data for A.thaliana
    # list2env(return_objects, envir = .GlobalEnv)



    
    #--------- Extract protein-coding genes, seperate them by strand and find overlap ----------


	GTF_df = as.data.frame(GTF)


	# Get all protein-coding genes
	GTF_df_cd_nc <- GTF_df[GTF_df$gene_biotype %in% c("protein_coding","lnc_exonic_antisense","lnc_intronic_antisense"), ]

	# Remove all chloroplast and mito genes
	# Reason: RNA-seq prep kit used is Ribo_zero, which incompletely removes Pt and
	# Mt transcripts, so expression estimates of those genes are likely to be wrong
	GTF_df_cd_nc <- subset(GTF_df_cd_nc, seqnames != "Pt")
	GTF_df_cd_nc <- subset(GTF_df_cd_nc, seqnames != "Mt")


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


	# Generate strand plus and strand minus GTF data frames with overlapping genes based on keys
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


	all_genes_counts <- tibble::rownames_to_column(all_genes_counts, "gene_id")
	all_genes_tpm <- tibble::rownames_to_column(all_genes_tpm, "gene_id")


	strand_plus_overlap_genes_counts <- merge(strand_plus_overlap_genes, all_genes_counts, by = "gene_id")
	strand_plus_overlap_genes_counts <- strand_plus_overlap_genes_counts[order(strand_plus_overlap_genes_counts$id),]

	strand_minus_overlap_genes_counts <- merge(strand_minus_overlap_genes, all_genes_counts, by = "gene_id")
	strand_minus_overlap_genes_counts <- strand_minus_overlap_genes_counts[order(strand_minus_overlap_genes_counts$id),]



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
	
			th_replicates <- do.call(cbind, lapply(split.default(df[3:ncol(df)], #adjust columns
								rep(seq_along(df), each = 3, length.out = ncol(df)-2)), #adjust columns
								function(x) {
									x[rowSums(x > threshold) < 2, ] <- 0; 
									x
								}
							))

			# Bind key/gene_id/id/gene_source/gene_biotype/seqnames/start/end/strand/width/width_overlap columns to thresholded data frame
			th_replicates <- cbind(df[1:11], th_replicates)

			# Remove all rows that only contain "0"
			th_replicates <- th_replicates[which(rowSums(th_replicates[,-1:-2, drop = FALSE] > 0) > 0),]

			return(th_replicates)
		}

		# Apply threshold to data and extract keys ("key")
		keys_data <- getThreshold(df)
		keys_data <- keys_data[, 1:2]
		names(keys_data) <- c("key", "ID")

		# Generate thresholded data frame based on keys
		th_df <- merge(keys_data, df, by = "key")
		th_df <- th_df[-1:-2]

		return(th_df)
	}


	# Apply threshold function for expression tables w/ pollen
	genes_tpm_th <- applyThreshold(all_genes_tpm, threshold)


	strand_plus_overlap_genes_counts_th <- strand_plus_overlap_genes_counts[(
		strand_plus_overlap_genes_counts$gene_id %in% genes_tpm_th$gene_id),]
	strand_minus_overlap_genes_counts_th <- strand_minus_overlap_genes_counts[(
		strand_minus_overlap_genes_counts$gene_id %in% genes_tpm_th$gene_id),]


	# Remove all gene pairs that show expression below threshold
	strand_plus_overlap_genes_counts_th <- strand_plus_overlap_genes_counts_th[(
		strand_plus_overlap_genes_counts_th$id %in% strand_minus_overlap_genes_counts_th$id),]
	strand_minus_overlap_genes_counts_th <- strand_minus_overlap_genes_counts_th[(
		strand_minus_overlap_genes_counts_th$id %in% strand_plus_overlap_genes_counts_th$id),]



	#------------------- Compute correlation between coding gene and cisNAT  -------------------


	# Compute spearman and pearson correlations
	getCor <- function(df1, df2) {

		df1_col <- ncol(df1)
		df2_col <- ncol(df2)

		# startup message
		message("Computing correlation...")

		df1$Spearman <- sapply(1:nrow(df1), function(i) 
	    	cor(as.numeric(df1[i, 11:df1_col]), as.numeric(df2[i, 11:df2_col]), method=c("spearman")))

		# log2-transform TPM values before computing Pearson
		df1[, 11:df1_col] <- log2(df1[, 11:df1_col] + 1)
		df2[, 11:df1_col] <- log2(df2[, 11:df2_col] + 1)

		df1$Pearson <- sapply(1:nrow(df1), function(i) 
	    	cor(as.numeric(df1[i, 11:df1_col]), as.numeric(df2[i, 11:df2_col]), method=c("pearson")))

		return(df1)
	}


	strand_plus_overlap_genes <- getCor(
		strand_plus_overlap_genes_counts_th, strand_minus_overlap_genes_counts_th)
	strand_minus_overlap_genes <- getCor(
		strand_minus_overlap_genes_counts_th, strand_plus_overlap_genes_counts_th)


	# Create a single table containing protein-coding/protein-coding and protein-coding/cisNAT pairs
	all_overlap_gene_pairs <- data.frame(
		id = strand_plus_overlap_genes$id,
		gene_id1 = strand_plus_overlap_genes$gene_id,
		gene_source1 = strand_plus_overlap_genes$gene_source,
		gene_biotype1 = strand_plus_overlap_genes$gene_biotype,
		seqnames1 = strand_plus_overlap_genes$seqnames,
		start1 = strand_plus_overlap_genes$start,
		end1 = strand_plus_overlap_genes$end,
		strand1 = strand_plus_overlap_genes$strand,
		width1 = strand_plus_overlap_genes$width,
		gene_id2 = strand_minus_overlap_genes$gene_id,
		gene_source2 = strand_minus_overlap_genes$gene_source,
		gene_biotype2 = strand_minus_overlap_genes$gene_biotype,
		seqnames2 = strand_minus_overlap_genes$seqnames,
		start2 = strand_minus_overlap_genes$start,
		end2 = strand_minus_overlap_genes$end,
		strand2 = strand_minus_overlap_genes$strand,
		width2 = strand_minus_overlap_genes$width,
		overlap = strand_minus_overlap_genes$width_overlap,
		Spearman = strand_plus_overlap_genes$Spearman,
		Pearson = strand_plus_overlap_genes$Pearson)


	# Make table of overlapping protein-coding/protein-coding gene pairs
	pc_pc_overlap_gene_pairs <- filter(all_overlap_gene_pairs, 
		gene_biotype1 == "protein_coding" & gene_biotype2 == "protein_coding")



#--------------------------------------- Write csv file ---------------------------------------


	# Set filename
    fname <- sprintf('%s.csv', paste(species_id, "cd_cd_SAS_cor", threshold, sep="_"))


	# Write final data tables to csv files and store them in /out_dir/output/data_tables
	if (!dir.exists(file.path(out_dir, "output", "overlap_pc_genes"))) 
		dir.create(file.path(out_dir, "output", "overlap_pc_genes"), recursive = TRUE)
	message("Storing results in: ", file.path("output", "overlap_pc_genes"))

	write.table(pc_pc_overlap_gene_pairs, file=file.path(out_dir, "output", "overlap_pc_genes", fname), 
		sep=";", dec=".", row.names = FALSE, col.names = TRUE)


}


