# Get overlap length between non-coding/protein-coding gene pairs
# Data input: 1) GTF file | 2) Expression_data WITHOUT mito and chloroplast genes
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


getNcPcOverlap <- function(species = c("AT", "AL", "CR", "ES", "TH", "MT", "BD")) {
	
	# Show error message if no species is chosen
    if (missing(species))
   
       stop(
       "Please choose one of the available species: 
	   'AT', 'AL', 'CR', 'ES', 'TH', 'MT', 'BD'",
	   call. = TRUE
       )

	# Set GTF input gtf file
    if (is.element("AT", species)) {
    	GTFfile = file.path(in_dir, "GTF", "AT_final_annotation.gtf")
        species_id <- "AT"

    } else if (is.element("AL", species)) {
		GTFfile = file.path(in_dir, "GTF", "AL_final_annotation.gtf")
		species_id <- "AL"

    } else if (is.element("CR", species)) {
		GTFfile = file.path(in_dir, "GTF", "CR_final_annotation.gtf")
		species_id <- "CR"

    } else if (is.element("ES", species)) {
		GTFfile = file.path(in_dir, "GTF", "ES_final_annotation.gtf")
		species_id <- "ES"

    } else if (is.element("TH", species)) {
		GTFfile = file.path(in_dir, "GTF", "TH_final_annotation.gtf")
		species_id <- "TH"

    } else if (is.element("MT", species)) {
		GTFfile = file.path(in_dir, "GTF", "MT_final_annotation.gtf")
		species_id <- "MT"

    } else if (is.element("BD", species)) {
		GTFfile = file.path(in_dir, "GTF", "BD_final_annotation.gtf")
		species_id <- "BD"
    }


	# Import gtf file
	GTF = import.gff(GTFfile, format="gtf", feature.type="gene")


    # return_list <- list("species_id" = species_id, "GTF" = GTF)
    # return(return_list)
    # }
    # return_objects <- getNcPcOverlap("AT") # read in GTF for A.thaliana
    # list2env(return_objects, envir = .GlobalEnv)



    
    #--------- Extract protein-coding genes, seperate them by strand and find overlap ----------


	GTF_df = as.data.frame(GTF)


	# Get all cis-natural antisense transcripts (cisNATs)
	GTF_df_cd_nc <- GTF_df[GTF_df$gene_biotype %in% c("protein_coding", "lnc_intronic_antisense", "lnc_exonic_antisense"), ]

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



	#- Extract overlapping cisNAT/protein coding genes from plus_strand/minus_strand data frames -


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

	## only keep cisNAT/protein-coding gene pairs

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
		overlap = strand_minus_overlap_genes$width_overlap)


	# Make sure that table only contains overlapping protein-coding gene pairs
	nc_pc_overlap_gene_pairs <- filter(all_overlap_gene_pairs, 
		!(gene_biotype1 == "protein_coding" & gene_biotype2 == "protein_coding"))



#--------------------------------------- Write csv file ---------------------------------------


	# Set filename
    fname <- sprintf('%s.csv', paste(species_id, "nd_cd_overlap", sep="_"))


	# Write final data tables to csv files and store them in /out_dir/output/data_tables
	if (!dir.exists(file.path(out_dir, "output", "overlap_nc_genes"))) 
		dir.create(file.path(out_dir, "output", "overlap_nc_genes"), recursive = TRUE)
	message("Storing results in: ", file.path("output", "overlap_nc_genes"))

	write.table(nc_pc_overlap_gene_pairs, file = file.path(out_dir, "output", "overlap_nc_genes", fname), 
		sep = ";", dec = ".", row.names = FALSE, col.names = TRUE)


}


