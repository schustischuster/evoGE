# This script loads and processes gene expression correlation table of sense-antisense (SAS) 
# gene pairs for coding/cisNATs SAS of the ATGE data set

#------------------- Load packages, set directories and read sample tables ---------------------


getDevSeq_ATGE <- function() {

	# Read ATGE_NAT ID table
	ATGE_NAT_ID <- read.table(file=file.path(in_dir, "ATGE_NAT", "ATGE_NAT.csv"), sep=";", dec=".", header=TRUE, stringsAsFactors = FALSE)
	ATH_cd_nc_SAS_cor_wo_pollen_0.5 <- read.table(file=file.path(out_dir, "output", "overlapp_nc_genes", "ATH_cd_nc_SAS_cor_wo_pollen_0.5.csv"), sep=";", dec=".", header=TRUE, stringsAsFactors = FALSE)



#-------------- Find 'ATGE_NAT_ID' and 'ATH_cd_nc_SAS_cor_wo_pollen_0.5' overlap ----------------


	# Work only with ATGE protein-coding gene IDs since lncRNA IDs may have changed

	# Extract protein-coding genes from DevSeq cd-nc SAS ATH data
	DevSeq_NAT_SAS_pair_plus_coding <- ATH_cd_nc_SAS_cor_wo_pollen_0.5 %>% filter(
		biotype_query == "protein_coding")
	DevSeq_NAT_plus_coding_ID <- DevSeq_NAT_SAS_pair_plus_coding %>% select(id_plus_strand)
	names(DevSeq_NAT_plus_coding_ID) <- "coding_id"

	DevSeq_NAT_SAS_pair_minus_coding <- ATH_cd_nc_SAS_cor_wo_pollen_0.5 %>% filter(
		biotype_subject == "protein_coding")
	DevSeq_NAT_minus_coding_ID <- DevSeq_NAT_SAS_pair_minus_coding %>% select(id_minus_strand)
	names(DevSeq_NAT_minus_coding_ID) <- "coding_id"

	# Combine plus and minus strand protein ID lists
	DevSeq_NAT_coding_ID <- rbind(DevSeq_NAT_plus_coding_ID, DevSeq_NAT_minus_coding_ID)

	# Select all protein-coding genes from DevSeq cd-nc SAS ATH data that are in ATGE SAS list 
	DevSeq_ATGE_NAT_coding_ID <- DevSeq_NAT_coding_ID[(
		((DevSeq_NAT_coding_ID$coding_id %in% ATGE_NAT_ID$AGI_1) | 
			(DevSeq_NAT_coding_ID$coding_id %in% ATGE_NAT_ID$AGI_2))
		),]

	# Select rows in "ATH_cd_nc_SAS_cor_wo_pollen_0.5" with a protein id in "DevSeq_ATGE_NAT_coding_ID"
	ATH_cd_nc_SAS_cor_wo_pollen_0.5_in_ATGE <- ATH_cd_nc_SAS_cor_wo_pollen_0.5[(
		((ATH_cd_nc_SAS_cor_wo_pollen_0.5$id_plus_strand %in% DevSeq_ATGE_NAT_coding_ID) | 
			(ATH_cd_nc_SAS_cor_wo_pollen_0.5$id_minus_strand %in% DevSeq_ATGE_NAT_coding_ID))
		),]



#---------------------------------------- Write output -----------------------------------------


	# Write final data tables to csv files and store them in /out_dir/output/data_tables
	if (!dir.exists(file.path(out_dir, "output", "SAS_DevSeq_ATGE"))) 
		dir.create(file.path(out_dir, "output", "SAS_DevSeq_ATGE"), recursive = TRUE)

	write.table(ATH_cd_nc_SAS_cor_wo_pollen_0.5_in_ATGE, 
		file=file.path(out_dir, "output", "SAS_DevSeq_ATGE", "ATH_cd_nc_SAS_cor_wo_pollen_0.5_in_ATGE.csv"), 
		sep=";", dec=".", row.names=FALSE, col.names=TRUE, col.names=NA)

}


