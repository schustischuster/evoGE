# This script loads gene expression correlation tables of sense-antisense (SAS) gene pairs for 
# both coding/cisNATs SAS and coding-coding SAS, defines TPM ranges, performs a Wilcox test
# and generates the plots for the sense-antisense pair expression correlation figure

# Expression correlation tables have the following format:
# id_plus_strand / start_plus / end_plus / strand_query / biotype_subject / id_minus_strand 
# / start_minus / end_minus / strand_subject / biotype_subject / Spearman / Pearson


#------------------- Load packages, set directories and read sample tables ---------------------


# Install and load packages
if (!require(dplyr)) install.packages('dplyr')
library(dplyr)
if (!require(ggplot2)) install.packages('ggplot2')
library(ggplot2)
if (!require(mgcv)) install.packages('mgcv')
library(mgcv)


# Set file path and input files
in_dir_cd <- "/Volumes/User/Shared/Christoph_manuscript/DevSeq_paper/Analysis/Analysis_2019/A_thaliana_gene_exression_map/20191121_CS_coding_cisNAT_analysis/output/overlap_cd_genes"
in_dir_nc <- "/Volumes/User/Shared/Christoph_manuscript/DevSeq_paper/Analysis/Analysis_2019/A_thaliana_gene_exression_map/20191121_CS_coding_cisNAT_analysis/output/overlap_nc_genes"
in_dir_ATGE <- "/Volumes/User/Shared/Christoph_manuscript/DevSeq_paper/Analysis/Analysis_2019/A_thaliana_gene_exression_map/20191121_CS_coding_cisNAT_analysis/output/SAS_DevSeq_ATGE"
in_dir_expr <- "/Volumes/User/Shared/Christoph_manuscript/DevSeq_paper/Analysis/Analysis_2019/A_thaliana_gene_exression_map/20191121_CS_coding_cisNAT_analysis/output/NAT_expr_cor"
out_dir <- "/Volumes/User/Shared/Christoph_manuscript/DevSeq_paper/Analysis/Analysis_2019/A_thaliana_gene_exression_map/20191121_CS_coding_cisNAT_analysis"


# Read all csv files in input file path
readTable <- function(path, pattern = "*.csv") {
    files = list.files(path, pattern, full.names = TRUE)
    lapply(files, function(x) read.table(x, sep=";", dec=".", header=TRUE, stringsAsFactors=FALSE))
}

coding_genes_tables <- readTable(in_dir_cd)
NAT_genes_tables <- readTable(in_dir_nc)
NAT_expr_cor <- readTable(in_dir_expr)


# Get file names and save them in character vector
coding_genes_table_list <- as.character(list.files(in_dir_cd, pattern = "*.csv"))
coding_genes_table_names <- gsub('\\.csv$', '', coding_genes_table_list)

NAT_genes_table_list <- as.character(list.files(in_dir_nc, pattern = "*.csv"))
NAT_genes_table_names <- gsub('\\.csv$', '', NAT_genes_table_list)

NAT_expr_cor_list <- as.character(list.files(in_dir_expr, pattern = "*.csv"))
NAT_expr_cor_names <- gsub('\\.csv$', '', NAT_expr_cor_list)


# Change data frame names in list
names(coding_genes_tables) <- coding_genes_table_names
list2env(coding_genes_tables, envir = .GlobalEnv)

names(NAT_genes_tables) <- NAT_genes_table_names
list2env(NAT_genes_tables, envir = .GlobalEnv)

names(NAT_expr_cor) <- NAT_expr_cor_names
list2env(NAT_expr_cor, envir = .GlobalEnv)


# Read ATGE_NAT ID table
ATGE_NAT_ID <- read.table(file=file.path(in_dir_ATGE, "ATH_cd_nc_SAS_cor_wo_pollen_0.5_in_ATGE.csv"), 
	sep=";", dec=".", header=TRUE, stringsAsFactors = FALSE)




#----------------------- Define TPM ranges for all species (w/o pollen)  -----------------------


# ATH all samples
ATH_cd_nc_SAS_cor_wo_pollen_0.5_2 <- ATH_cd_nc_SAS_cor_wo_pollen_0.5 %>% filter(
	!((id_plus_strand %in% ATH_cd_nc_SAS_cor_wo_pollen_2$id_plus_strand) &
		(id_minus_strand %in% ATH_cd_nc_SAS_cor_wo_pollen_2$id_minus_strand)))
ATH_cd_nc_SAS_cor_wo_pollen_2_5 <- ATH_cd_nc_SAS_cor_wo_pollen_2 %>% filter(
	!((id_plus_strand %in% ATH_cd_nc_SAS_cor_wo_pollen_5$id_plus_strand)
	 & (id_minus_strand %in% ATH_cd_nc_SAS_cor_wo_pollen_5$id_minus_strand)))
ATH_cd_nc_SAS_cor_wo_pollen_5_10 <- ATH_cd_nc_SAS_cor_wo_pollen_5 %>% filter(
	!((id_plus_strand %in% ATH_cd_nc_SAS_cor_wo_pollen_10$id_plus_strand)
	 & (id_minus_strand %in% ATH_cd_nc_SAS_cor_wo_pollen_10$id_minus_strand)))
# ATH comparative samples
ATH_comparative_samples_cd_nc_SAS_cor_wo_pollen_0.5_2 <- ATH_comparative_samples_cd_nc_SAS_cor_wo_pollen_0.5 %>% filter(
	!((id_plus_strand %in% ATH_comparative_samples_cd_nc_SAS_cor_wo_pollen_2$id_plus_strand) &
		(id_minus_strand %in% ATH_comparative_samples_cd_nc_SAS_cor_wo_pollen_2$id_minus_strand)))
ATH_comparative_samples_cd_nc_SAS_cor_wo_pollen_2_5 <- ATH_comparative_samples_cd_nc_SAS_cor_wo_pollen_2 %>% filter(
	!((id_plus_strand %in% ATH_comparative_samples_cd_nc_SAS_cor_wo_pollen_5$id_plus_strand)
	 & (id_minus_strand %in% ATH_comparative_samples_cd_nc_SAS_cor_wo_pollen_5$id_minus_strand)))
ATH_comparative_samples_cd_nc_SAS_cor_wo_pollen_5_10 <- ATH_comparative_samples_cd_nc_SAS_cor_wo_pollen_5 %>% filter(
	!((id_plus_strand %in% ATH_comparative_samples_cd_nc_SAS_cor_wo_pollen_10$id_plus_strand)
	 & (id_minus_strand %in% ATH_comparative_samples_cd_nc_SAS_cor_wo_pollen_10$id_minus_strand)))
# AL comparative samples
AL_comparative_samples_cd_nc_SAS_cor_wo_pollen_0.5_2 <- AL_comparative_samples_cd_nc_SAS_cor_wo_pollen_0.5 %>% filter(
	!((id_plus_strand %in% AL_comparative_samples_cd_nc_SAS_cor_wo_pollen_2$id_plus_strand) &
		(id_minus_strand %in% AL_comparative_samples_cd_nc_SAS_cor_wo_pollen_2$id_minus_strand)))
AL_comparative_samples_cd_nc_SAS_cor_wo_pollen_2_5 <- AL_comparative_samples_cd_nc_SAS_cor_wo_pollen_2 %>% filter(
	!((id_plus_strand %in% AL_comparative_samples_cd_nc_SAS_cor_wo_pollen_5$id_plus_strand)
	 & (id_minus_strand %in% AL_comparative_samples_cd_nc_SAS_cor_wo_pollen_5$id_minus_strand)))
AL_comparative_samples_cd_nc_SAS_cor_wo_pollen_5_10 <- AL_comparative_samples_cd_nc_SAS_cor_wo_pollen_5 %>% filter(
	!((id_plus_strand %in% AL_comparative_samples_cd_nc_SAS_cor_wo_pollen_10$id_plus_strand)
	 & (id_minus_strand %in% AL_comparative_samples_cd_nc_SAS_cor_wo_pollen_10$id_minus_strand)))
# CR
CR_cd_nc_SAS_cor_wo_pollen_0.5_2 <- CR_cd_nc_SAS_cor_wo_pollen_0.5 %>% filter(
	!((id_plus_strand %in% CR_cd_nc_SAS_cor_wo_pollen_2$id_plus_strand) &
		(id_minus_strand %in% CR_cd_nc_SAS_cor_wo_pollen_2$id_minus_strand)))
CR_cd_nc_SAS_cor_wo_pollen_2_5 <- CR_cd_nc_SAS_cor_wo_pollen_2 %>% filter(
	!((id_plus_strand %in% CR_cd_nc_SAS_cor_wo_pollen_5$id_plus_strand)
	 & (id_minus_strand %in% CR_cd_nc_SAS_cor_wo_pollen_5$id_minus_strand)))
CR_cd_nc_SAS_cor_wo_pollen_5_10 <- CR_cd_nc_SAS_cor_wo_pollen_5 %>% filter(
	!((id_plus_strand %in% CR_cd_nc_SAS_cor_wo_pollen_10$id_plus_strand)
	 & (id_minus_strand %in% CR_cd_nc_SAS_cor_wo_pollen_10$id_minus_strand)))
# ES
ES_cd_nc_SAS_cor_wo_pollen_0.5_2 <- ES_cd_nc_SAS_cor_wo_pollen_0.5 %>% filter(
	!((id_plus_strand %in% ES_cd_nc_SAS_cor_wo_pollen_2$id_plus_strand) &
		(id_minus_strand %in% ES_cd_nc_SAS_cor_wo_pollen_2$id_minus_strand)))
ES_cd_nc_SAS_cor_wo_pollen_2_5 <- ES_cd_nc_SAS_cor_wo_pollen_2 %>% filter(
	!((id_plus_strand %in% ES_cd_nc_SAS_cor_wo_pollen_5$id_plus_strand)
	 & (id_minus_strand %in% ES_cd_nc_SAS_cor_wo_pollen_5$id_minus_strand)))
ES_cd_nc_SAS_cor_wo_pollen_5_10 <- ES_cd_nc_SAS_cor_wo_pollen_5 %>% filter(
	!((id_plus_strand %in% ES_cd_nc_SAS_cor_wo_pollen_10$id_plus_strand)
	 & (id_minus_strand %in% ES_cd_nc_SAS_cor_wo_pollen_10$id_minus_strand)))
# TH
TH_cd_nc_SAS_cor_wo_pollen_0.5_2 <- TH_cd_nc_SAS_cor_wo_pollen_0.5 %>% filter(
	!((id_plus_strand %in% TH_cd_nc_SAS_cor_wo_pollen_2$id_plus_strand) &
		(id_minus_strand %in% TH_cd_nc_SAS_cor_wo_pollen_2$id_minus_strand)))
TH_cd_nc_SAS_cor_wo_pollen_2_5 <- TH_cd_nc_SAS_cor_wo_pollen_2 %>% filter(
	!((id_plus_strand %in% TH_cd_nc_SAS_cor_wo_pollen_5$id_plus_strand)
	 & (id_minus_strand %in% TH_cd_nc_SAS_cor_wo_pollen_5$id_minus_strand)))
TH_cd_nc_SAS_cor_wo_pollen_5_10 <- TH_cd_nc_SAS_cor_wo_pollen_5 %>% filter(
	!((id_plus_strand %in% TH_cd_nc_SAS_cor_wo_pollen_10$id_plus_strand)
	 & (id_minus_strand %in% TH_cd_nc_SAS_cor_wo_pollen_10$id_minus_strand)))
# MT
MT_cd_nc_SAS_cor_wo_pollen_0.5_2 <- MT_cd_nc_SAS_cor_wo_pollen_0.5 %>% filter(
	!((id_plus_strand %in% MT_cd_nc_SAS_cor_wo_pollen_2$id_plus_strand) &
		(id_minus_strand %in% MT_cd_nc_SAS_cor_wo_pollen_2$id_minus_strand)))
MT_cd_nc_SAS_cor_wo_pollen_2_5 <- MT_cd_nc_SAS_cor_wo_pollen_2 %>% filter(
	!((id_plus_strand %in% MT_cd_nc_SAS_cor_wo_pollen_5$id_plus_strand)
	 & (id_minus_strand %in% MT_cd_nc_SAS_cor_wo_pollen_5$id_minus_strand)))
MT_cd_nc_SAS_cor_wo_pollen_5_10 <- MT_cd_nc_SAS_cor_wo_pollen_5 %>% filter(
	!((id_plus_strand %in% MT_cd_nc_SAS_cor_wo_pollen_10$id_plus_strand)
	 & (id_minus_strand %in% MT_cd_nc_SAS_cor_wo_pollen_10$id_minus_strand)))
# BD
BD_cd_nc_SAS_cor_wo_pollen_0.5_2 <- BD_cd_nc_SAS_cor_wo_pollen_0.5 %>% filter(
	!((id_plus_strand %in% BD_cd_nc_SAS_cor_wo_pollen_2$id_plus_strand) &
		(id_minus_strand %in% BD_cd_nc_SAS_cor_wo_pollen_2$id_minus_strand)))
BD_cd_nc_SAS_cor_wo_pollen_2_5 <- BD_cd_nc_SAS_cor_wo_pollen_2 %>% filter(
	!((id_plus_strand %in% BD_cd_nc_SAS_cor_wo_pollen_5$id_plus_strand)
	 & (id_minus_strand %in% BD_cd_nc_SAS_cor_wo_pollen_5$id_minus_strand)))
BD_cd_nc_SAS_cor_wo_pollen_5_10 <- BD_cd_nc_SAS_cor_wo_pollen_5 %>% filter(
	!((id_plus_strand %in% BD_cd_nc_SAS_cor_wo_pollen_10$id_plus_strand)
	 & (id_minus_strand %in% BD_cd_nc_SAS_cor_wo_pollen_10$id_minus_strand)))




#------------------------ Prepare data for plotting and test statistics ------------------------


# Create list of coding SAS cor tables excluding pollen samples
coding_genes_tables_wo_pollen_list <- list(
	AL_coding_SAS_cor_wo_pollen = AL_coding_SAS_cor_wo_pollen, 
	AL_comp_samples_coding_SAS_cor_wo_pollen = AL_comparative_samples_coding_SAS_cor_wo_pollen, 
	ATH_coding_SAS_cor_wo_pollen = ATH_coding_SAS_cor_wo_pollen,
	ATH_coding_SAS_cor = ATH_coding_SAS_cor, #including pollen samples for ATH
	ATH_comp_samples_coding_SAS_cor_wo_pollen = ATH_comparative_samples_coding_SAS_cor_wo_pollen,
	ATH_comp_samples_coding_SAS_cor = ATH_comparative_samples_coding_SAS_cor, #including pollen samples for ATH
	BD_coding_SAS_cor_wo_pollen = BD_coding_SAS_cor_wo_pollen,
	CR_coding_SAS_cor_wo_pollen = CR_coding_SAS_cor_wo_pollen,
	ES_coding_SAS_cor_wo_pollen = ES_coding_SAS_cor_wo_pollen,
	MT_coding_SAS_cor_wo_pollen = MT_coding_SAS_cor_wo_pollen,
	TH_coding_SAS_cor_wo_pollen = TH_coding_SAS_cor_wo_pollen)


# Create list of non-coding NAT SAS cor tables excluding pollen samples
NAT_genes_tables_wo_pollen_list <- list(
	ATH_cd_nc_SAS_cor_wo_pollen_0.5 = ATH_cd_nc_SAS_cor_wo_pollen_0.5,
	ATH_cd_nc_SAS_cor_wo_pollen_2 = ATH_cd_nc_SAS_cor_wo_pollen_2,
	ATH_cd_nc_SAS_cor_wo_pollen_5 = ATH_cd_nc_SAS_cor_wo_pollen_5,
	ATH_cd_nc_SAS_cor_wo_pollen_0.5_2 = ATH_cd_nc_SAS_cor_wo_pollen_0.5_2,
	ATH_cd_nc_SAS_cor_wo_pollen_2_5 = ATH_cd_nc_SAS_cor_wo_pollen_2_5,
	ATH_cd_nc_SAS_cor_wo_pollen_5_10 = ATH_cd_nc_SAS_cor_wo_pollen_5_10,
	ATH_cd_nc_SAS_cor_wo_pollen_10 = ATH_cd_nc_SAS_cor_wo_pollen_10,
	ATH_cd_nc_SAS_cor_0.5 = ATH_cd_nc_SAS_cor_0.5, #including pollen samples for ATH
	AL_comp_samples_cd_nc_SAS_cor_wo_pollen_0.5 = AL_comparative_samples_cd_nc_SAS_cor_wo_pollen_0.5, 
	AL_comp_samples_cd_nc_SAS_cor_wo_pollen_2 = AL_comparative_samples_cd_nc_SAS_cor_wo_pollen_2, 
	AL_comp_samples_cd_nc_SAS_cor_wo_pollen_5 = AL_comparative_samples_cd_nc_SAS_cor_wo_pollen_5,
	AL_comp_samples_cd_nc_SAS_cor_wo_pollen_0.5_2 = AL_comparative_samples_cd_nc_SAS_cor_wo_pollen_0.5_2,
	AL_comp_samples_cd_nc_SAS_cor_wo_pollen_2_5 = AL_comparative_samples_cd_nc_SAS_cor_wo_pollen_2_5,
	AL_comp_samples_cd_nc_SAS_cor_wo_pollen_5_10 = AL_comparative_samples_cd_nc_SAS_cor_wo_pollen_5_10,
	AL_comp_samples_cd_nc_SAS_cor_wo_pollen_10 = AL_comparative_samples_cd_nc_SAS_cor_wo_pollen_10,
	ATH_comp_samples_cd_nc_SAS_cor_wo_pollen_0.5 = ATH_comparative_samples_cd_nc_SAS_cor_wo_pollen_0.5,
	ATH_comp_samples_cd_nc_SAS_cor_wo_pollen_2 = ATH_comparative_samples_cd_nc_SAS_cor_wo_pollen_2, 
	ATH_comp_samples_cd_nc_SAS_cor_wo_pollen_5 = ATH_comparative_samples_cd_nc_SAS_cor_wo_pollen_5,
	ATH_comp_samples_cd_nc_SAS_cor_wo_pollen_0.5_2 = ATH_comparative_samples_cd_nc_SAS_cor_wo_pollen_0.5_2,
	ATH_comp_samples_cd_nc_SAS_cor_wo_pollen_2_5 = ATH_comparative_samples_cd_nc_SAS_cor_wo_pollen_2_5,
	ATH_comp_samples_cd_nc_SAS_cor_wo_pollen_5_10 = ATH_comparative_samples_cd_nc_SAS_cor_wo_pollen_5_10,
	ATH_comp_samples_cd_nc_SAS_cor_wo_pollen_10 = ATH_comparative_samples_cd_nc_SAS_cor_wo_pollen_10,
	ATH_comp_samples_cd_nc_SAS_cor_0.5 = ATH_comparative_samples_cd_nc_SAS_cor_0.5, #including pollen samples for ATH
	BD_cd_nc_SAS_cor_wo_pollen_0.5 = BD_cd_nc_SAS_cor_wo_pollen_0.5, 
	BD_cd_nc_SAS_cor_wo_pollen_2 = BD_cd_nc_SAS_cor_wo_pollen_2,
	BD_cd_nc_SAS_cor_wo_pollen_5 = BD_cd_nc_SAS_cor_wo_pollen_5,
	BD_cd_nc_SAS_cor_wo_pollen_0.5_2 = BD_cd_nc_SAS_cor_wo_pollen_0.5_2,
	BD_cd_nc_SAS_cor_wo_pollen_2_5 = BD_cd_nc_SAS_cor_wo_pollen_2_5,
	BD_cd_nc_SAS_cor_wo_pollen_5_10 = BD_cd_nc_SAS_cor_wo_pollen_5_10,
	BD_cd_nc_SAS_cor_wo_pollen_10 = BD_cd_nc_SAS_cor_wo_pollen_10,
	CR_cd_nc_SAS_cor_wo_pollen_0.5 = CR_cd_nc_SAS_cor_wo_pollen_0.5,
	CR_cd_nc_SAS_cor_wo_pollen_2 = CR_cd_nc_SAS_cor_wo_pollen_2,
	CR_cd_nc_SAS_cor_wo_pollen_5 = CR_cd_nc_SAS_cor_wo_pollen_5,
	CR_cd_nc_SAS_cor_wo_pollen_0.5_2 = CR_cd_nc_SAS_cor_wo_pollen_0.5_2,
	CR_cd_nc_SAS_cor_wo_pollen_2_5 = CR_cd_nc_SAS_cor_wo_pollen_2_5,
	CR_cd_nc_SAS_cor_wo_pollen_5_10 = CR_cd_nc_SAS_cor_wo_pollen_5_10,
	CR_cd_nc_SAS_cor_wo_pollen_10 = CR_cd_nc_SAS_cor_wo_pollen_10,
	ES_cd_nc_SAS_cor_wo_pollen_0.5 = ES_cd_nc_SAS_cor_wo_pollen_0.5,
	ES_cd_nc_SAS_cor_wo_pollen_2 = ES_cd_nc_SAS_cor_wo_pollen_2,
	ES_cd_nc_SAS_cor_wo_pollen_5 = ES_cd_nc_SAS_cor_wo_pollen_5,
	ES_cd_nc_SAS_cor_wo_pollen_0.5_2 = ES_cd_nc_SAS_cor_wo_pollen_0.5_2,
	ES_cd_nc_SAS_cor_wo_pollen_2_5 = ES_cd_nc_SAS_cor_wo_pollen_2_5,
	ES_cd_nc_SAS_cor_wo_pollen_5_10 = ES_cd_nc_SAS_cor_wo_pollen_5_10,
	ES_cd_nc_SAS_cor_wo_pollen_10 = ES_cd_nc_SAS_cor_wo_pollen_10,
	MT_cd_nc_SAS_cor_wo_pollen_0.5 = MT_cd_nc_SAS_cor_wo_pollen_0.5,
	MT_cd_nc_SAS_cor_wo_pollen_2 = MT_cd_nc_SAS_cor_wo_pollen_2,
	MT_cd_nc_SAS_cor_wo_pollen_5 = MT_cd_nc_SAS_cor_wo_pollen_5,
	MT_cd_nc_SAS_cor_wo_pollen_0.5_2 = MT_cd_nc_SAS_cor_wo_pollen_0.5_2,
	MT_cd_nc_SAS_cor_wo_pollen_2_5 = MT_cd_nc_SAS_cor_wo_pollen_2_5,
	MT_cd_nc_SAS_cor_wo_pollen_5_10 = MT_cd_nc_SAS_cor_wo_pollen_5_10,
	MT_cd_nc_SAS_cor_wo_pollen_10 = MT_cd_nc_SAS_cor_wo_pollen_10,
	TH_cd_nc_SAS_cor_wo_pollen_0.5 = TH_cd_nc_SAS_cor_wo_pollen_0.5,
	TH_cd_nc_SAS_cor_wo_pollen_2 = TH_cd_nc_SAS_cor_wo_pollen_2,
	TH_cd_nc_SAS_cor_wo_pollen_5 = TH_cd_nc_SAS_cor_wo_pollen_5,
	TH_cd_nc_SAS_cor_wo_pollen_0.5_2 = TH_cd_nc_SAS_cor_wo_pollen_0.5_2,
	TH_cd_nc_SAS_cor_wo_pollen_2_5 = TH_cd_nc_SAS_cor_wo_pollen_2_5,
	TH_cd_nc_SAS_cor_wo_pollen_5_10 = TH_cd_nc_SAS_cor_wo_pollen_5_10,
	TH_cd_nc_SAS_cor_wo_pollen_10 = TH_cd_nc_SAS_cor_wo_pollen_10)


# Extract Spearman and Pearson correlations from expression cor tables excluding pollen samples
# Correlation values will be stored in a numeric atomic vector with the original data frame name plus
# attached type of correlation measure, e.g. "AL_coding_SAS_cor_wo_pollen_pearson"

# Coding SAS pairs
coding_genes_tables_wo_pollen_list_spearman <- lapply(coding_genes_tables_wo_pollen_list, function(x) {
	as.numeric(unlist(x[, 11]))})
names(coding_genes_tables_wo_pollen_list_spearman) <- paste(names(
	coding_genes_tables_wo_pollen_list_spearman),"_spearman", sep="")

coding_genes_tables_wo_pollen_list_pearson <- lapply(coding_genes_tables_wo_pollen_list, function(x) {
	as.numeric(unlist(x[, 12]))})
names(coding_genes_tables_wo_pollen_list_pearson) <- paste(names(
	coding_genes_tables_wo_pollen_list_pearson),"_pearson", sep="")


# non-coding NAT / protein-coding SAS pairs
NAT_genes_tables_wo_pollen_list_spearman <- lapply(NAT_genes_tables_wo_pollen_list, function(x) {
	as.numeric(unlist(x[, 19]))})
names(NAT_genes_tables_wo_pollen_list_spearman) <- paste(names(
	NAT_genes_tables_wo_pollen_list_spearman),"_spearman", sep="")

NAT_genes_tables_wo_pollen_list_pearson <- lapply(NAT_genes_tables_wo_pollen_list, function(x) {
	as.numeric(unlist(x[, 20]))})
names(NAT_genes_tables_wo_pollen_list_pearson) <- paste(names(
	NAT_genes_tables_wo_pollen_list_pearson),"_pearson", sep="")


# Combine coding SAS and non-coding NAT / protein-coding SAS lists
SAS_pairs_list_pearson <- append(
	coding_genes_tables_wo_pollen_list_pearson, NAT_genes_tables_wo_pollen_list_pearson)
SAS_pairs_list_spearman <- append(
	coding_genes_tables_wo_pollen_list_spearman, NAT_genes_tables_wo_pollen_list_spearman)


# Show names of all correlation vectors
names(SAS_pairs_list_pearson)
names(SAS_pairs_list_spearman)

list2env(SAS_pairs_list_pearson, envir = .GlobalEnv)
list2env(SAS_pairs_list_spearman, envir = .GlobalEnv)


# Extract pearson and spearman correlations from DevSeq-ATGE expression cor tables
ATH_cd_nc_SAS_cor_wo_pollen_0.5_in_ATGE_spearman <- as.numeric(unlist(
	ATGE_NAT_ID[, 19]))
ATH_cd_nc_SAS_cor_wo_pollen_0.5_in_ATGE_pearson <- as.numeric(unlist(
	ATGE_NAT_ID[, 20]))


# Create "plots" folder in /out_dir/output/plots
if (!dir.exists(file.path(out_dir, "output", "plots"))) 
	dir.create(file.path(out_dir, "output", "plots"), recursive = TRUE)




#--------------- Prepare data for plotting NAT expression intensities and ratios ---------------


# Generate data tables that are classified by cisNAT-coding gene correlation ranges
cd_nc_expr_below_min03_list <- lapply(NAT_expr_cor, function(x) {
	dplyr::filter(x, x$Pearson < -0.3)})
cd_nc_expr_betw_min02_02_list <- lapply(NAT_expr_cor, function(x) {
	dplyr::filter(x, ((x$Pearson > -0.2) & (x$Pearson < 0.2)))})
cd_nc_expr_above_05_list <- lapply(NAT_expr_cor, function(x) {
	dplyr::filter(x, x$Pearson > 0.5)})


# Set names in NAT expression lists
element_names <- names(NAT_expr_cor)
cd_nc_expr_below_min03_list_names <- paste(element_names,"below_min03", sep='_')
cd_nc_expr_betw_min02_02_list_names <- paste(element_names,"betw_min02_02", sep='_')
cd_nc_expr_above_05_list_names <- paste(element_names,"above_05", sep='_')

names(cd_nc_expr_below_min03_list) <- cd_nc_expr_below_min03_list_names
names(cd_nc_expr_betw_min02_02_list) <- cd_nc_expr_betw_min02_02_list_names
names(cd_nc_expr_above_05_list) <- cd_nc_expr_above_05_list_names

list2env(cd_nc_expr_below_min03_list, envir = .GlobalEnv)
list2env(cd_nc_expr_betw_min02_02_list, envir = .GlobalEnv)
list2env(cd_nc_expr_above_05_list, envir = .GlobalEnv)


# Prepare data for ggplot2 density plot w/ expression data below_min03, between_min02_02, above_05
# For both max coding and NAT gene expression
combineExprDataCdNc <- function(below_min03, between_min02_02, above_05) {

	number_values <- 2*(nrow(below_min03)+nrow(between_min02_02)+nrow(above_05))
	species_name = as.data.frame(rep(c(sub("\\_.*", "", deparse(substitute(below_min03)))),each=number_values))
	names(species_name) <- "species"

	class_0 = as.data.frame(rep(c(">05_cd"),each=nrow(above_05)))
	names(class_0) <- "class"
	class_1 = as.data.frame(rep(c(">-02 <02_cd"),each=nrow(between_min02_02)))
	names(class_1) <- "class"
	class_2 = as.data.frame(rep(c("-0.3_cd"),each=nrow(below_min03)))
	names(class_2) <- "class"
	class_3 = as.data.frame(rep(c(">05_nc"),each=nrow(above_05)))
	names(class_3) <- "class"
	class_4 = as.data.frame(rep(c(">-02 <02_nc"),each=nrow(between_min02_02)))
	names(class_4) <- "class"
	class_5 = as.data.frame(rep(c("-0.3_nc"),each=nrow(below_min03)))
	names(class_5) <- "class"

	expr_values_0 = as.data.frame(above_05$max_coding)
	names(expr_values_0) <- "max_expression"
	expr_values_1 = as.data.frame(between_min02_02$max_coding)
	names(expr_values_1) <- "max_expression"
	expr_values_2 = as.data.frame(below_min03$max_coding)
	names(expr_values_2) <- "max_expression"
	expr_values_3 = as.data.frame(above_05$max_NAT)
	names(expr_values_3) <- "max_expression"
	expr_values_4 = as.data.frame(between_min02_02$max_NAT)
	names(expr_values_4) <- "max_expression"
	expr_values_5 = as.data.frame(below_min03$max_NAT)
	names(expr_values_5) <- "max_expression"

	expression_df = data.frame(species_name, rbind(class_0, class_1, class_2, class_3, class_4, class_5), 
		rbind(expr_values_0, expr_values_1, expr_values_2, expr_values_3, expr_values_4, expr_values_5))
	expression_df <- na.omit(expression_df)

	return(expression_df)
}


ATH_all_cd_nc_max_expr_pearson <- combineExprDataCdNc(ATH_all_cd_nc_expr_below_min03, 
	ATH_all_cd_nc_expr_betw_min02_02, ATH_all_cd_nc_expr_above_05)
ATH_comp_cd_nc_max_expr_pearson <- combineExprDataCdNc(ATH_comp_cd_nc_expr_below_min03, 
	ATH_comp_cd_nc_expr_betw_min02_02, ATH_comp_cd_nc_expr_above_05)
AL_cd_nc_max_expr_pearson <- combineExprDataCdNc(AL_cd_nc_expr_below_min03, 
	AL_cd_nc_expr_betw_min02_02, AL_cd_nc_expr_above_05)
CR_cd_nc_max_expr_pearson <- combineExprDataCdNc(CR_cd_nc_expr_below_min03, 
	CR_cd_nc_expr_betw_min02_02, CR_cd_nc_expr_above_05)
ES_cd_nc_max_expr_pearson <- combineExprDataCdNc(ES_cd_nc_expr_below_min03, 
	ES_cd_nc_expr_betw_min02_02, ES_cd_nc_expr_above_05)
TH_cd_nc_max_expr_pearson <- combineExprDataCdNc(TH_cd_nc_expr_below_min03, 
	TH_cd_nc_expr_betw_min02_02, TH_cd_nc_expr_above_05)
MT_cd_nc_max_expr_pearson <- combineExprDataCdNc(MT_cd_nc_expr_below_min03, 
	MT_cd_nc_expr_betw_min02_02, MT_cd_nc_expr_above_05)
BD_cd_nc_max_expr_pearson <- combineExprDataCdNc(BD_cd_nc_expr_below_min03, 
	BD_cd_nc_expr_betw_min02_02, MT_cd_nc_expr_above_05)


# Perform wilcox test on all combinations
wilcoxMaxNAT <- function(x,y) {
	z <- wilcox.test(x[,17], y[,17])$p.value
	z <- formatC(z, format = "e", digits = 0)
}

ATH_all_wilcox_03_02 <- wilcoxMaxNAT(ATH_all_cd_nc_expr_below_min03, ATH_all_cd_nc_expr_betw_min02_02)
ATH_all_wilcox_03_05 <- wilcoxMaxNAT(ATH_all_cd_nc_expr_below_min03, ATH_all_cd_nc_expr_above_05)
ATH_all_wilcox_02_05 <- wilcoxMaxNAT(ATH_all_cd_nc_expr_betw_min02_02, ATH_all_cd_nc_expr_above_05)
ATH_comp_wilcox_03_02 <- wilcoxMaxNAT(ATH_comp_cd_nc_expr_below_min03, ATH_comp_cd_nc_expr_betw_min02_02)
ATH_comp_wilcox_03_05 <- wilcoxMaxNAT(ATH_comp_cd_nc_expr_below_min03, ATH_comp_cd_nc_expr_above_05)
ATH_comp_wilcox_02_05 <- wilcoxMaxNAT(ATH_comp_cd_nc_expr_betw_min02_02, ATH_comp_cd_nc_expr_above_05)
AL_wilcox_03_02 <- wilcoxMaxNAT(AL_cd_nc_expr_below_min03, AL_cd_nc_expr_betw_min02_02)
AL_wilcox_03_05 <- wilcoxMaxNAT(AL_cd_nc_expr_below_min03, AL_cd_nc_expr_above_05)
AL_wilcox_02_05 <- wilcoxMaxNAT(AL_cd_nc_expr_betw_min02_02, AL_cd_nc_expr_above_05)
CR_wilcox_03_02 <- wilcoxMaxNAT(CR_cd_nc_expr_below_min03, CR_cd_nc_expr_betw_min02_02)
CR_wilcox_03_05 <- wilcoxMaxNAT(CR_cd_nc_expr_below_min03, CR_cd_nc_expr_above_05)
CR_wilcox_02_05 <- wilcoxMaxNAT(CR_cd_nc_expr_betw_min02_02, CR_cd_nc_expr_above_05)
ES_wilcox_03_02 <- wilcoxMaxNAT(ES_cd_nc_expr_below_min03, ES_cd_nc_expr_betw_min02_02)
ES_wilcox_03_05 <- wilcoxMaxNAT(ES_cd_nc_expr_below_min03, ES_cd_nc_expr_above_05)
ES_wilcox_02_05 <- wilcoxMaxNAT(ES_cd_nc_expr_betw_min02_02, ES_cd_nc_expr_above_05)
TH_wilcox_03_02 <- wilcoxMaxNAT(TH_cd_nc_expr_below_min03, TH_cd_nc_expr_betw_min02_02)
TH_wilcox_03_05 <- wilcoxMaxNAT(TH_cd_nc_expr_below_min03, TH_cd_nc_expr_above_05)
TH_wilcox_02_05 <- wilcoxMaxNAT(TH_cd_nc_expr_betw_min02_02, TH_cd_nc_expr_above_05)
MT_wilcox_03_02 <- wilcoxMaxNAT(MT_cd_nc_expr_below_min03, MT_cd_nc_expr_betw_min02_02)
MT_wilcox_03_05 <- wilcoxMaxNAT(MT_cd_nc_expr_below_min03, MT_cd_nc_expr_above_05)
MT_wilcox_02_05 <- wilcoxMaxNAT(MT_cd_nc_expr_betw_min02_02, MT_cd_nc_expr_above_05)
BD_wilcox_03_02 <- wilcoxMaxNAT(BD_cd_nc_expr_below_min03, BD_cd_nc_expr_betw_min02_02)
BD_wilcox_03_05 <- wilcoxMaxNAT(BD_cd_nc_expr_below_min03, MT_cd_nc_expr_above_05)
BD_wilcox_02_05 <- wilcoxMaxNAT(BD_cd_nc_expr_betw_min02_02, MT_cd_nc_expr_above_05)


# Function to scatter plot max expression versus pearson correlation
makeScrPlotMaxExpr <- function(data, lim_y, p03_02, p03_05, p02_05,
	plot_title = c("ATH_all", "ATH", "CR", "ES", "TH", "MT", "BD"), yadj) {

	fname <- sprintf('%s.jpg', paste(deparse(substitute(data)), sep="_"))

	n_ac <- nrow(subset(data, class=="-0.3_nc"))
	n_nc <- nrow(subset(data, class==">-02 <02_nc"))
	n_pc <- nrow(subset(data, class==">05_nc"))
	n_ac = paste("n=", n_ac, sep="")
	n_nc = paste("n=", n_nc, sep="")
	n_pc = paste("n=", n_pc, sep="")

	cor03_02 = paste("vs.    ", " p=", p03_02, sep="")
	cor03_05 = paste("vs.    ", " p=", p03_05, sep="")
	cor02_05 = paste("vs.    ", " p=", p02_05, sep="")

	blu = rgb(0, 70, 139, max = 255, alpha = 35)
	gray = rgb(131, 145, 145, max = 255, alpha = 60)
	grn = rgb(94, 200, 100, max = 255, alpha = 75)
	blu_ln = rgb(0, 70, 139, max = 255, alpha = 0)
	gray_ln = rgb(131, 145, 145, max = 255, alpha = 0)
	grn_ln = rgb(73, 180, 60, max = 255, alpha = 0)

	p <- ggplot(data, aes(x=max_expression, group=class, fill=class, colour=class, linetype=class)) +
	geom_density(adjust=1.35, size=1.25) + 
	scale_x_continuous(limits = c(0,12), expand = c(0, 0)) +
	scale_y_continuous(limits = lim_y, expand = c(0, 0)) + 
	annotate("rect", xmin=c(6.5,8.25,6.5,8.25,6.5,8.25,6.5,7.45,8.25,6.5,7.45,8.25,6.5,7.45,8.25), 
		xmax=c(7,8.75,7,8.75,7,8.75,7.075,7.68,8.75,7.075,7.68,8.75,7.075,7.68,8.75), 
		ymin=c(0.4158*yadj,0.4158*yadj,0.3853*yadj,0.3853*yadj,0.3548*yadj,0.3548*yadj,0.3336*yadj,0.3336*yadj,0.3243*yadj,0.3031*yadj,0.3031*yadj,0.2938*yadj,0.2726*yadj,0.2726*yadj,0.2633*yadj), 
		ymax=c(0.4348*yadj,0.4348*yadj,0.4043*yadj,0.4043*yadj,0.3738*yadj,0.3738*yadj,0.3340*yadj,0.3340*yadj,0.3433*yadj,0.3035*yadj,0.3035*yadj,0.3128*yadj,0.2730*yadj,0.2730*yadj,0.2823*yadj), 
		color=c("#49b43c","#00468b","#49b43c","#839191","#00468b","#839191","#49b43c","#49b43c","#49b43c","#00468b","#00468b","#00468b","#839191","#839191","#839191"), 
		size=1.2, fill=c(grn,blu,grn,gray,blu,gray,grn,grn,grn,blu,blu,blu,gray,gray,gray)) + 
	annotate("segment", x=c(6), xend=c(12), y=c(0.248*yadj), yend=c(0.248*yadj), color="grey20", size=0.7) + 
	annotate("segment", x=c(6.03), xend=c(6.03), y=c(0.247*yadj), yend=c(0.449*yadj), color="grey20", size=0.7) + 
	annotate("text", x = -Inf, y = Inf, hjust = -1.644, vjust = 1.8, size=5.5, label = cor03_02) + 
	annotate("text", x = -Inf, y = Inf, hjust = -1.644, vjust = 3.425, size=5.5, label = cor03_05) + 
	annotate("text", x = -Inf, y = Inf, hjust = -1.644, vjust = 5.05, size=5.5, label = cor02_05) + 
	annotate("text", x = 9.066, y = Inf, hjust = 0, vjust = 6.675, size=5.5, label = n_ac) + 
	annotate("text", x = 9.066, y = Inf, hjust = 0, vjust = 8.3, size=5.5, label = n_nc) + 
	annotate("text", x = 9.066, y = Inf, hjust = 0, vjust = 9.925, size=5.5, label = n_pc)
	q <- p + ggtitle(plot_title) + theme_bw() + scale_fill_manual(values = c(gray_ln, blu_ln, grn_ln, gray, blu, grn)) +
		scale_color_manual(values = c("#839191", "#00468b", "#52b540", "#839191", "#00468b", "#52b540")) + xlab("Expression (log2 TPM)") + ylab("Density") + 
		scale_linetype_manual(values = c("dotdash","dotdash","dotdash", "solid","solid","solid")) + 
		theme(text=element_text(size=16), 
  		axis.ticks.length = unit(.3, "cm"),
  		plot.margin = unit(c(5.5, 10.5, 3.0, 3.5), "points"),
  		axis.text.x = element_text(colour = "black", size=16, angle=0, margin = margin(t = 5, r = 0, b = 0, l = 0)), 
  		axis.text.y = element_text(colour = "black", size=16, angle=0, margin = margin(t = 0, r = 5, b = 0, l = 0)),
  		axis.title.x = element_text(colour = "black", size=17.5, margin = margin(t = 14.5, r = 0, b = 0, l = 0)),
  		axis.title.y = element_text(colour = "black", size=17.5, margin = margin(t = 0, r = 12.5, b = 0, l = 0)),
  		plot.title = element_text(colour = "black", size=17.5, margin = margin(t = 17.25, r = 0, b = 8, l = 0), hjust = 0.5),
  		legend.position = "bottom",
  		panel.border = element_rect(colour = "black", fill=NA, size=0.5))

  	ggsave(file = file.path(out_dir, "output", "plots", fname), plot = q,
		scale = 1, width = 4.848485, height = 5.95, units = c("in"), 
		dpi = 825, limitsize = FALSE)
}

makeScrPlotMaxExpr(data=ATH_all_cd_nc_max_expr_pearson, lim_y=c(0,0.45), p03_02=ATH_all_wilcox_03_02, p03_05=ATH_all_wilcox_03_05, p02_05=ATH_all_wilcox_02_05, plot_title="ATH_all", yadj=1)
makeScrPlotMaxExpr(data=ATH_comp_cd_nc_max_expr_pearson, lim_y=c(0,0.545), p03_02=ATH_comp_wilcox_03_02, p03_05=ATH_comp_wilcox_03_05, p02_05=ATH_comp_wilcox_02_05, plot_title="ATH_comp", yadj=1.211)
makeScrPlotMaxExpr(data=AL_cd_nc_max_expr_pearson, lim_y=c(0,0.552), p03_02=AL_wilcox_03_02, p03_05=AL_wilcox_03_05, p02_05=AL_wilcox_02_05, plot_title="AL_", yadj=1.2267)
makeScrPlotMaxExpr(data=CR_cd_nc_max_expr_pearson, lim_y=c(0,0.457), p03_02=CR_wilcox_03_02, p03_05=CR_wilcox_03_05, p02_05=CR_wilcox_02_05, plot_title="CR_", yadj=1.0155)
makeScrPlotMaxExpr(data=ES_cd_nc_max_expr_pearson, lim_y=c(0,0.480), p03_02=ES_wilcox_03_02, p03_05=ES_wilcox_03_05, p02_05=ES_wilcox_02_05, plot_title="ES_", yadj=1.0667)
makeScrPlotMaxExpr(data=TH_cd_nc_max_expr_pearson, lim_y=c(0,0.479), p03_02=TH_wilcox_03_02, p03_05=TH_wilcox_03_05, p02_05=TH_wilcox_02_05, plot_title="TH_", yadj=1.0645)
makeScrPlotMaxExpr(data=MT_cd_nc_max_expr_pearson, lim_y=c(0,0.771), p03_02=MT_wilcox_03_02, p03_05=MT_wilcox_03_05, p02_05=MT_wilcox_02_05, plot_title="MT_", yadj=1.7135)
makeScrPlotMaxExpr(data=BD_cd_nc_max_expr_pearson, lim_y=c(0,0.590), p03_02=BD_wilcox_03_02, p03_05=BD_wilcox_03_05, p02_05=BD_wilcox_02_05, plot_title="BD_", yadj=1.3112)




# Prepare data for ggplot2 density plot w/ expression data below_min03, between_min02_02, above_05
# For coding and NAT gene expression ratio
combineExprDataRatio <- function(below_min03, between_min02_02, above_05) {

	number_values <- (nrow(below_min03)+nrow(between_min02_02)+nrow(above_05))
	species_name = as.data.frame(rep(c(sub("\\_.*", "", deparse(substitute(below_min03)))),each=number_values))
	names(species_name) <- "species"

	class_0 = as.data.frame(rep(c(">05"),each=nrow(above_05)))
	names(class_0) <- "class"
	class_1 = as.data.frame(rep(c(">-02 <02"),each=nrow(between_min02_02)))
	names(class_1) <- "class"
	class_2 = as.data.frame(rep(c("-0.3"),each=nrow(below_min03)))
	names(class_2) <- "class"

	expr_values_0 = as.data.frame(above_05$max_ratio_nc_cd)
	names(expr_values_0) <- "max_expression"
	expr_values_1 = as.data.frame(between_min02_02$max_ratio_nc_cd)
	names(expr_values_1) <- "max_expression"
	expr_values_2 = as.data.frame(below_min03$max_ratio_nc_cd)
	names(expr_values_2) <- "max_expression"

	expression_df = data.frame(species_name, rbind(class_0, class_1, class_2), 
		rbind(expr_values_0, expr_values_1, expr_values_2))
	expression_df <- na.omit(expression_df)

	return(expression_df)
}


ATH_all_cd_nc_max_expr_ratio <- combineExprDataRatio(ATH_all_cd_nc_expr_below_min03, 
	ATH_all_cd_nc_expr_betw_min02_02, ATH_all_cd_nc_expr_above_05)
ATH_comp_cd_nc_max_expr_ratio <- combineExprDataRatio(ATH_comp_cd_nc_expr_below_min03, 
	ATH_comp_cd_nc_expr_betw_min02_02, ATH_comp_cd_nc_expr_above_05)
AL_cd_nc_max_expr_ratio <- combineExprDataRatio(AL_cd_nc_expr_below_min03, 
	AL_cd_nc_expr_betw_min02_02, AL_cd_nc_expr_above_05)
CR_cd_nc_max_expr_ratio <- combineExprDataRatio(CR_cd_nc_expr_below_min03, 
	CR_cd_nc_expr_betw_min02_02, CR_cd_nc_expr_above_05)
ES_cd_nc_max_expr_ratio <- combineExprDataRatio(ES_cd_nc_expr_below_min03, 
	ES_cd_nc_expr_betw_min02_02, ES_cd_nc_expr_above_05)
TH_cd_nc_max_expr_ratio <- combineExprDataRatio(TH_cd_nc_expr_below_min03, 
	TH_cd_nc_expr_betw_min02_02, TH_cd_nc_expr_above_05)
MT_cd_nc_max_expr_ratio <- combineExprDataRatio(MT_cd_nc_expr_below_min03, 
	MT_cd_nc_expr_betw_min02_02, MT_cd_nc_expr_above_05)
BD_cd_nc_max_expr_ratio <- combineExprDataRatio(BD_cd_nc_expr_below_min03, 
	BD_cd_nc_expr_betw_min02_02, MT_cd_nc_expr_above_05)


# Perform wilcox test on all combinations
wilcoxRatioSAS <- function(x,y) {
	z <- wilcox.test(x[,21], y[,21])$p.value
	z <- formatC(z, format = "e", digits = 0)
}

ATH_all_wilcox_03_02_ratio <- wilcoxRatioSAS(ATH_all_cd_nc_expr_below_min03, ATH_all_cd_nc_expr_betw_min02_02)
ATH_all_wilcox_03_05_ratio <- wilcoxRatioSAS(ATH_all_cd_nc_expr_below_min03, ATH_all_cd_nc_expr_above_05)
ATH_all_wilcox_02_05_ratio <- wilcoxRatioSAS(ATH_all_cd_nc_expr_betw_min02_02, ATH_all_cd_nc_expr_above_05)
ATH_comp_wilcox_03_02_ratio <- wilcoxRatioSAS(ATH_comp_cd_nc_expr_below_min03, ATH_comp_cd_nc_expr_betw_min02_02)
ATH_comp_wilcox_03_05_ratio <- wilcoxRatioSAS(ATH_comp_cd_nc_expr_below_min03, ATH_comp_cd_nc_expr_above_05)
ATH_comp_wilcox_02_05_ratio <- wilcoxRatioSAS(ATH_comp_cd_nc_expr_betw_min02_02, ATH_comp_cd_nc_expr_above_05)
AL_wilcox_03_02_ratio <- wilcoxRatioSAS(AL_cd_nc_expr_below_min03, AL_cd_nc_expr_betw_min02_02)
AL_wilcox_03_05_ratio <- wilcoxRatioSAS(AL_cd_nc_expr_below_min03, AL_cd_nc_expr_above_05)
AL_wilcox_02_05_ratio <- wilcoxRatioSAS(AL_cd_nc_expr_betw_min02_02, AL_cd_nc_expr_above_05)
CR_wilcox_03_02_ratio <- wilcoxRatioSAS(CR_cd_nc_expr_below_min03, CR_cd_nc_expr_betw_min02_02)
CR_wilcox_03_05_ratio <- wilcoxRatioSAS(CR_cd_nc_expr_below_min03, CR_cd_nc_expr_above_05)
CR_wilcox_02_05_ratio <- wilcoxRatioSAS(CR_cd_nc_expr_betw_min02_02, CR_cd_nc_expr_above_05)
ES_wilcox_03_02_ratio <- wilcoxRatioSAS(ES_cd_nc_expr_below_min03, ES_cd_nc_expr_betw_min02_02)
ES_wilcox_03_05_ratio <- wilcoxRatioSAS(ES_cd_nc_expr_below_min03, ES_cd_nc_expr_above_05)
ES_wilcox_02_05_ratio <- wilcoxRatioSAS(ES_cd_nc_expr_betw_min02_02, ES_cd_nc_expr_above_05)
TH_wilcox_03_02_ratio <- wilcoxRatioSAS(TH_cd_nc_expr_below_min03, TH_cd_nc_expr_betw_min02_02)
TH_wilcox_03_05_ratio <- wilcoxRatioSAS(TH_cd_nc_expr_below_min03, TH_cd_nc_expr_above_05)
TH_wilcox_02_05_ratio <- wilcoxRatioSAS(TH_cd_nc_expr_betw_min02_02, TH_cd_nc_expr_above_05)
MT_wilcox_03_02_ratio <- wilcoxRatioSAS(MT_cd_nc_expr_below_min03, MT_cd_nc_expr_betw_min02_02)
MT_wilcox_03_05_ratio <- wilcoxRatioSAS(MT_cd_nc_expr_below_min03, MT_cd_nc_expr_above_05)
MT_wilcox_02_05_ratio <- wilcoxRatioSAS(MT_cd_nc_expr_betw_min02_02, MT_cd_nc_expr_above_05)
BD_wilcox_03_02_ratio <- wilcoxRatioSAS(BD_cd_nc_expr_below_min03, BD_cd_nc_expr_betw_min02_02)
BD_wilcox_03_05_ratio <- wilcoxRatioSAS(BD_cd_nc_expr_below_min03, MT_cd_nc_expr_above_05)
BD_wilcox_02_05_ratio <- wilcoxRatioSAS(BD_cd_nc_expr_betw_min02_02, MT_cd_nc_expr_above_05)


# Function to scatter plot max expression versus pearson correlation
makeScrPlotExprRatio <- function(data, lim_y, p03_02, p03_05, p02_05,
	plot_title = c("ATH_all", "ATH", "CR", "ES", "TH", "MT", "BD")) {

	fname <- sprintf('%s.jpg', paste(deparse(substitute(data)), sep="_"))

	n_ac <- nrow(subset(data, class=="-0.3_nc"))
	n_nc <- nrow(subset(data, class==">-02 <02_nc"))
	n_pc <- nrow(subset(data, class==">05_nc"))

	cor03_02 = paste("ac vs. nc ", "p=", p03_02, sep="")
	cor03_05 = paste("ac vs. pc ", "p=", p03_05, sep="")
	cor02_05 = paste("nc vs. pc ", "p=", p02_05, sep="")

	blu = rgb(0, 70, 139, max = 255, alpha = 0)
	gray = rgb(131, 145, 145, max = 255, alpha = 0)
	grn = rgb(94, 200, 100, max = 255, alpha = 0)

	p <- ggplot(data, aes(x=max_expression, group=class, fill=class, colour=class, linetype=class)) +
	geom_density(adjust=1.35, size=1.25) + 
	scale_x_continuous(trans='log10', labels = prettyNum, limits = c(0.01,100), expand = c(0, 0)) +
	scale_y_continuous(limits = lim_y, expand = c(0, 0)) + 
	annotation_logticks(sides = 'b') + 
	annotate("text", x = 1.12, y = Inf, hjust = 0, vjust = 1.8, size=5.5, label = cor03_02) + 
	annotate("text", x = 1.12, y = Inf, hjust = 0, vjust = 3.425, size=5.5, label = cor03_05) + 
	annotate("text", x = 1.12, y = Inf, hjust = 0, vjust = 5.05, size=5.5, label = cor02_05)
	q <- p + ggtitle(plot_title) + theme_bw() + scale_fill_manual(values = c(gray, blu, grn)) +
		scale_color_manual(values = c("#839191", "#00468b", "#52b540")) + xlab("Expression ratio (nc:cd SAS)") + ylab("Density") + 
		scale_linetype_manual(values = c("solid","solid","solid")) + 
		theme(text=element_text(size=16), 
  		axis.ticks.length = unit(.3, "cm"),
  		plot.margin = unit(c(5.5, 13.5, 20.25, 0.5), "points"),
  		axis.text.x = element_text(colour = "black", size=16, angle=0, margin = margin(t = 5, r = 0, b = 0, l = 0)), 
  		axis.text.y = element_text(colour = "black", size=16, angle=0, margin = margin(t = 0, r = 5, b = 0, l = 0)),
  		axis.title.x = element_text(colour = "black", size=17.5, margin = margin(t = 14.5, r = 0, b = 0, l = 0)),
  		axis.title.y = element_text(colour = "black", size=17.5, margin = margin(t = 0, r = 12.5, b = 0, l = 0)),
  		plot.title = element_text(colour = "black", size=17.5, margin = margin(t = 17.25, r = 0, b = 8, l = 0), hjust = 0.5),
  		legend.position = "bottom",
  		panel.border = element_rect(colour = "black", fill=NA, size=0.5))

  	ggsave(file = file.path(out_dir, "output", "plots", fname), plot = q,
		scale = 1, width = 4.848485, height = 5.95, units = c("in"), 
		dpi = 825, limitsize = FALSE)
}

makeScrPlotExprRatio(data=ATH_all_cd_nc_max_expr_ratio, lim_y=c(0,2.105), p03_02=ATH_all_wilcox_03_02_ratio, p03_05=ATH_all_wilcox_03_05_ratio, p02_05=ATH_all_wilcox_02_05_ratio, plot_title="ATH_all")
makeScrPlotExprRatio(data=ATH_comp_cd_nc_max_expr_ratio, lim_y=c(0,1.615), p03_02=ATH_comp_wilcox_03_02_ratio, p03_05=ATH_comp_wilcox_03_05_ratio, p02_05=ATH_comp_wilcox_02_05_ratio, plot_title="ATH_comp")
makeScrPlotExprRatio(data=AL_cd_nc_max_expr_ratio, lim_y=c(0,1.42), p03_02=AL_wilcox_03_02_ratio, p03_05=AL_wilcox_03_05_ratio, p02_05=AL_wilcox_02_05_ratio, plot_title="AL_")
makeScrPlotExprRatio(data=CR_cd_nc_max_expr_ratio, lim_y=c(0,1.79), p03_02=CR_wilcox_03_02_ratio, p03_05=CR_wilcox_03_05_ratio, p02_05=CR_wilcox_02_05_ratio, plot_title="CR_")
makeScrPlotExprRatio(data=ES_cd_nc_max_expr_ratio, lim_y=c(0,1.935), p03_02=ES_wilcox_03_02_ratio, p03_05=ES_wilcox_03_05_ratio, p02_05=ES_wilcox_02_05_ratio, plot_title="ES_")
makeScrPlotExprRatio(data=TH_cd_nc_max_expr_ratio, lim_y=c(0,1.5), p03_02=TH_wilcox_03_02_ratio, p03_05=TH_wilcox_03_05_ratio, p02_05=TH_wilcox_02_05_ratio, plot_title="TH_")
makeScrPlotExprRatio(data=MT_cd_nc_max_expr_ratio, lim_y=c(0,1.685), p03_02=MT_wilcox_03_02_ratio, p03_05=MT_wilcox_03_05_ratio, p02_05=MT_wilcox_02_05_ratio, plot_title="MT_")
makeScrPlotExprRatio(data=BD_cd_nc_max_expr_ratio, lim_y=c(0,1.682), p03_02=BD_wilcox_03_02_ratio, p03_05=BD_wilcox_03_05_ratio, p02_05=BD_wilcox_02_05_ratio, plot_title="BD_")




#-------------------------------- Perform Wilcox rank sum test ---------------------------------


wilcox_pearson_cor <- sapply(SAS_pairs_list_pearson, function(x) sapply(
	SAS_pairs_list_pearson, function(y) wilcox.test(x,y)$p.value))
wilcox_spearman_cor <- sapply(SAS_pairs_list_spearman, function(x) sapply(
	SAS_pairs_list_spearman, function(y) wilcox.test(x,y)$p.value))


write.table(wilcox_pearson_cor, file=file.path(out_dir, "output", "plots", "wilcox_pearson_cor.csv"), 
	sep=";", dec=".", row.names=TRUE, col.names=NA)
write.table(wilcox_spearman_cor, file=file.path(out_dir, "output", "plots", "wilcox_spearman_cor.csv"), 
	sep=";", dec=".", row.names=TRUE, col.names=NA)




#---------------- Generate pearson plots for all species and expression ranges -----------------


# Make boxplot of result
# Pearson plot of cd-cd SAS / nc-cd SAS pairs ATH all samples and comparative samples
n_ATH_pc_all_wo_pollen <- length(ATH_coding_SAS_cor_wo_pollen_pearson)
n_ATH_nc_all_wo_pollen <- length(ATH_cd_nc_SAS_cor_wo_pollen_0.5_pearson)
n_ATH_pc_comp_wo_pollen <- length(ATH_comp_samples_coding_SAS_cor_wo_pollen_pearson)
n_ATH_nc_comp_wo_pollen <- length(ATH_comp_samples_cd_nc_SAS_cor_wo_pollen_0.5_pearson)

png(file=file.path(out_dir, "output", "plots", "cd_cd_SAS_NAT_cd_SAS_pearson_ATH_all_vs_comp.png"), 
	width = 2850, height = 4000, res = 825)
par(mar = c(4.5, 4.5, 4, 2.5))
boxplot(ATH_coding_SAS_cor_wo_pollen_pearson, ATH_cd_nc_SAS_cor_wo_pollen_0.5_pearson, 
	ATH_comp_samples_coding_SAS_cor_wo_pollen_pearson, ATH_comp_samples_cd_nc_SAS_cor_wo_pollen_0.5_pearson,
	ylim = c(-1.2, 1.2), 
	names = FALSE, 
	xaxt = 'n', 
	yaxt = 'n', 
	cex.lab = 1.1, 
	las = 2,
	cex.axis = 1.1, #adapt size of axis labels
	ylab = "Pearson ρ", 
	col = c("#a8a8a8", "#d8a900", "#a8a8a8", "#d8a900"), 
	boxwex = 0.85, 
	pars = list(outcol = "gray50"), 
	lwd = 1.35, 
	whisklty = 1, 
	at = c(1,2,3.5,4.5), 
	notch = FALSE
	)
	title("SAS pairs in ATH", adj = 0.50, line = 1.3, font.main = 1, cex.main = 1.2)
	rug(x = c(2.75), ticksize = -0.13, side = 1, lwd = 1.35, col = "gray60") #x-axis ticks
	abline(v = c(2.75), col = "gray60")
	box(lwd = 1.35)
	axis(side = 2, lwd = 1.35, las = 2)
	text(x= 1.5, y = 1.15, labels= "p<1e-100", col = "black", cex = 1) #ATH_all p-value
	text(x= 4, y = 1.15, labels= "p<1e-50", col = "black", cex = 1) #ATH_comp p-value
	text(x= 1, y= -1.175, labels= n_ATH_pc_all_wo_pollen, col= "gray40", cex= 0.97) #ATH_all no.genes
	text(x= 2, y= -1.04, labels= n_ATH_nc_all_wo_pollen, col= "gray40", cex= 0.97)
	text(x= 3.5, y= -1.175, labels= n_ATH_pc_comp_wo_pollen, col= "gray40", cex= 0.97) #AL_comp no.genes
	text(x= 4.5, y= -1.04, labels= n_ATH_nc_comp_wo_pollen, col= "gray40", cex= 0.97)
	mtext('ATH_all', side = 1, line = 0.55, at = 1.5)
	mtext('ATH_comp', side = 1, line = 0.55, at = 4)
	par(xpd=TRUE)
	legend(-0.35,-1.6,c("cd-cd SAS", "nc-cd SAS"),  
	bty='n', horiz = TRUE, fill = c("#a8a8a8", "#d8a900"), cex = 1.1, x.intersp = 0.5)
dev.off()



# Pearson plot of cd-cd SAS / nc-cd SAS pairs all species comparative samples
n_ATH_pc_comp_wo_pollen <- length(ATH_comp_samples_coding_SAS_cor_wo_pollen_pearson)
n_ATH_nc_comp_wo_pollen <- length(ATH_comp_samples_cd_nc_SAS_cor_wo_pollen_0.5_pearson)
n_AL_pc_comp_wo_pollen <- length(AL_comp_samples_coding_SAS_cor_wo_pollen_pearson)
n_AL_nc_comp_wo_pollen <- length(AL_comp_samples_cd_nc_SAS_cor_wo_pollen_0.5_pearson)
n_CR_pc_wo_pollen <- length(CR_coding_SAS_cor_wo_pollen_pearson)
n_CR_nc_wo_pollen <- length(CR_cd_nc_SAS_cor_wo_pollen_0.5_pearson)
n_ES_pc_wo_pollen <- length(ES_coding_SAS_cor_wo_pollen_pearson)
n_ES_nc_wo_pollen <- length(ES_cd_nc_SAS_cor_wo_pollen_0.5_pearson)
n_TH_pc_wo_pollen <- length(TH_coding_SAS_cor_wo_pollen_pearson)
n_TH_nc_wo_pollen <- length(TH_cd_nc_SAS_cor_wo_pollen_0.5_pearson)
n_MT_pc_wo_pollen <- length(MT_coding_SAS_cor_wo_pollen_pearson)
n_MT_nc_wo_pollen <- length(MT_cd_nc_SAS_cor_wo_pollen_0.5_pearson)
n_BD_pc_wo_pollen <- length(BD_coding_SAS_cor_wo_pollen_pearson)
n_BD_nc_wo_pollen <- length(BD_cd_nc_SAS_cor_wo_pollen_0.5_pearson)

png(file = file.path(out_dir, "output", "plots", "cd_cd_SAS_NAT_cd_SAS_pearson_wo_pollen.png"), 
	width = 7200, height = 4000, res = 825)
par(mar = c(4.5, 4.5, 4, 1.5))
boxplot(ATH_comp_samples_coding_SAS_cor_wo_pollen_pearson, ATH_comp_samples_cd_nc_SAS_cor_wo_pollen_0.5_pearson,
	AL_comp_samples_coding_SAS_cor_wo_pollen_pearson, AL_comp_samples_cd_nc_SAS_cor_wo_pollen_0.5_pearson,
	CR_coding_SAS_cor_wo_pollen_pearson, CR_cd_nc_SAS_cor_wo_pollen_0.5_pearson, 
	ES_coding_SAS_cor_wo_pollen_pearson, ES_cd_nc_SAS_cor_wo_pollen_0.5_pearson, 
	TH_coding_SAS_cor_wo_pollen_pearson, TH_cd_nc_SAS_cor_wo_pollen_0.5_pearson, 
	MT_coding_SAS_cor_wo_pollen_pearson, MT_cd_nc_SAS_cor_wo_pollen_0.5_pearson, 
	BD_coding_SAS_cor_wo_pollen_pearson, BD_cd_nc_SAS_cor_wo_pollen_0.5_pearson, 
	ylim = c(-1.2, 1.2), 
	names = FALSE, 
	xaxt='n', 
	yaxt='n', 
	cex.lab = 1.1, 
	las = 2,
	cex.axis = 1.1, #adapt size of axis labels
	ylab = "Pearson ρ", 
	col = c("#a8a8a8", "#d8a900", "#a8a8a8", "#d8a900", "#a8a8a8", "#d8a900", "#a8a8a8", 
		"#d8a900", "#a8a8a8", "#d8a900", "#a8a8a8", "#d8a900", "#a8a8a8", "#d8a900"), 
	boxwex = 0.85, 
	pars = list(outcol = "gray50"), 
	lwd = 1.35, 
	whisklty = 1, 
	at = c(1,2,4,5,7,8,10,11,13,14,16,17,19,20), 
	notch = FALSE
	)
	title("SAS pairs in all species", adj = 0.5, line = 1.3, font.main = 1, cex.main = 1.2)
	rug(x = c(3, 6, 9, 12, 15, 18), ticksize = -0.08, side = 1, lwd=1.35, col="gray60") #x-axis ticks
	abline(v = c(3, 6, 9, 12, 15, 18), col="gray60")
	box(lwd = 1.35)
	axis(side=2, lwd = 1.35, las = 2)
	text(x= 1.5, y= 1.15, labels= "p<1e-50", col= "black", cex=1) #ATH p-value
	text(x= 4.5, y= 1.15, labels= "p<1e-30", col= "black", cex=1) #AL p-value
	text(x= 7.5, y= 1.15, labels= "p<1e-50", col= "black", cex=1) #CR p-value
	text(x= 10.5, y= 1.15, labels= "p<1e-40", col= "black", cex=1) #ES p-value
	text(x= 13.5, y= 1.15, labels= "p<1e-15", col= "black", cex=1) #TH p-value
	text(x= 16.5, y= 1.15, labels= "p<1e-50", col= "black", cex=1) #MT p-value
	text(x= 19.5, y= 1.15, labels= "p<1e-40", col= "black", cex=1) #BD p-value
	text(x= 1, y= -1.175, labels= n_ATH_pc_comp_wo_pollen, col= "gray40", cex=0.97) #ATH no.genes
	text(x= 2, y= -1.04, labels= n_ATH_nc_comp_wo_pollen, col= "gray40", cex=0.97)
	text(x= 4, y= -1.175, labels= n_AL_pc_comp_wo_pollen, col= "gray40", cex=0.97) #AL no.genes
	text(x= 5, y= -1.04, labels= n_AL_nc_comp_wo_pollen, col= "gray40", cex=0.97)
	text(x= 7, y= -1.175, labels= n_CR_pc_wo_pollen, col= "gray40", cex=0.97) #CR no.genes
	text(x= 8, y= -1.04, labels= n_CR_nc_wo_pollen, col= "gray40", cex=0.97)
	text(x= 10, y= -1.175, labels= n_ES_pc_wo_pollen, col= "gray40", cex=0.97) #ES no.genes
	text(x= 11, y= -1.04, labels= n_ES_nc_wo_pollen, col= "gray40", cex=0.97)
	text(x= 13, y= -1.175, labels= n_TH_pc_wo_pollen, col= "gray40", cex=0.97) #TH no.genes
	text(x= 14, y= -1.04, labels= n_TH_nc_wo_pollen, col= "gray40", cex=0.97)
	text(x= 16, y= -1.175, labels= n_MT_pc_wo_pollen, col= "gray40", cex=0.97) #MT no.genes
	text(x= 17, y= -1.04, labels= n_MT_nc_wo_pollen, col= "gray40", cex=0.97)
	text(x= 19, y= -1.175, labels= n_BD_pc_wo_pollen, col= "gray40", cex=0.97) #BD no.genes
	text(x= 20, y= -1.04, labels= n_BD_nc_wo_pollen, col= "gray40", cex=0.97)
	mtext('ATH', side=1, line=0.5, at=1.5)
	mtext('AL', side=1, line=0.5, at=4.5)
	mtext('CR', side=1, line=0.5, at=7.5)
	mtext('ES', side=1, line=0.5, at=10.5)
	mtext('TH', side=1, line=0.5, at=13.5)
	mtext('MT', side=1, line=0.5, at=16.5)
	mtext('BD', side=1, line=0.5, at=19.5)
	par(xpd=TRUE)
	legend(6.55,-1.6,c("cd-cd SAS", "nc-cd SAS"),  
	bty='n', horiz=TRUE, fill=c("#a8a8a8", "#d8a900"), cex=1.1, x.intersp = 0.5)
dev.off()



# Pearson plot of nc-cd SAS pairs with all thresholds
make_Boxplot_All_Thresholds_Labels <- function(threshold_05_2, threshold_2_5, threshold_5_10, 
	threshold_greater10, samples=c("all","comparative")) {

	species <- sub("\\_.*", "", deparse(substitute(threshold_05_2)))
    n_values_05_2 <- length(threshold_05_2)
    n_values_2_5 <- length(threshold_2_5)
    n_values_5_10 <- length(threshold_5_10)
    n_values_greater10 <- length(threshold_greater10)

    if (is.element("all", samples))
    	title_plot = paste(species, "all", sep="_")
    if (is.element("comparative", samples))
    	title_plot = paste(species, "comp", sep="_")
    if (missing(samples))
    	title_plot = species

    fname <- sprintf('%s.png', paste(title_plot, "thresholds", sep="_")) 

	png(file = file.path(out_dir, "output", "plots", fname), 
		width = 2620, height = 4000, res = 825)
	par(mar = c(5.725, 4.5, 4, 1))
	boxplot(threshold_05_2, threshold_2_5, threshold_5_10, threshold_greater10,
		ylim = c(-1.1, 1.1), 
		yaxt='n', 
		cex.lab = 1.1, 
		las = 1,
		cex.axis = 1.1, #adapt size of axis labels
		xlab = "", 
		ylab = "Pearson ρ", 
		col = c("#d8a900", "#00bc1f", "#00c094", "#00beda"), 
		boxwex = 0.71, 
		pars = list(outcol = "gray50"), 
		lwd = 1.35, 
		whisklty = 1, 
		at = c(1,2,3,4), 
		notch = FALSE
		)
		title(title_plot, adj = 0.5, line = 1.25, font.main = 1, cex.main = 1.2)
		title(xlab = "cis-NAT expression (TPM)", line = 2.65, cex.lab = 1.1)
		box(lwd = 1.35)
		axis(side=2, lwd = 1.35, las = 2)
		text(x= 1, y= -1.05, labels= n_values_05_2, col= "gray40", cex=0.97) #threshold_>0.5
		text(x= 2, y= -1.05, labels= n_values_2_5, col= "gray40", cex=0.97) #threshold_0.5-2
		text(x= 3, y= -1.05, labels= n_values_5_10, col= "gray40", cex=0.97) #threshold_2-5
		text(x= 4, y= -1.05, labels= n_values_greater10, col= "gray40", cex=0.97) #threshold_greater5
		mtext('0.5-2', side=1, line=0.85, at=1)
		mtext('2-5', side=1, line=0.85, at=2)
		mtext('5-10', side=1, line=0.85, at=3)
		mtext('>10', side=1, line=0.85, at=4)
		par(xpd=TRUE)
	dev.off()
}


# Pearson plot of nc-cd SAS pairs with all thresholds
make_Boxplot_All_Thresholds <- function(threshold_05_2, threshold_2_5, threshold_5_10, 
	threshold_greater10, samples=c("all","comparative")) {

	species <- sub("\\_.*", "", deparse(substitute(threshold_05_2)))
    n_values_05_2 <- length(threshold_05_2)
    n_values_2_5 <- length(threshold_2_5)
    n_values_5_10 <- length(threshold_5_10)
    n_values_greater10 <- length(threshold_greater10)

    if (is.element("all", samples))
    	title_plot = paste(species, "all", sep="_")
    if (is.element("comparative", samples))
    	title_plot = paste(species, "comp", sep="_")
    if (missing(samples))
    	title_plot = species

    fname <- sprintf('%s.png', paste(title_plot, "thresholds", sep="_")) 

	png(file = file.path(out_dir, "output", "plots", fname), 
		width = 2620, height = 4000, res = 825)
	par(mar = c(5.725, 4.5, 4, 1))
	boxplot(threshold_05_2, threshold_2_5, threshold_5_10, threshold_greater10,
		ylim = c(-1.1, 1.1), 
		yaxt='n', 
		cex.lab = 1.1, 
		las = 1,
		cex.axis = 1.1, #adapt size of axis labels
		xlab = "", 
		col = c("#d8a900", "#00bc1f", "#00c094", "#00beda"), 
		boxwex = 0.71, 
		lwd = 1.35, 
		whisklty = 1, 
		at = c(1,2,3,4), 
		pars = list(outcol = "gray50"),
		notch = FALSE
		)
		title(title_plot, adj = 0.5, line = 1.25, font.main = 1, cex.main = 1.2)
		title(xlab = "cis-NAT expression (TPM)", line = 2.65, cex.lab = 1.1)
		box(lwd = 1.35)
		axis(side=2, lwd = 1.35, las = 2)
		text(x= 1, y= -1.05, labels= n_values_05_2, col= "gray40", cex=0.97) #threshold_>0.5
		text(x= 2, y= -1.05, labels= n_values_2_5, col= "gray40", cex=0.97) #threshold_0.5-2
		text(x= 3, y= -1.05, labels= n_values_5_10, col= "gray40", cex=0.97) #threshold_2-5
		text(x= 4, y= -1.05, labels= n_values_greater10, col= "gray40", cex=0.97) #threshold_greater5
		mtext('0.5-2', side=1, line=0.85, at=1)
		mtext('2-5', side=1, line=0.85, at=2)
		mtext('5-10', side=1, line=0.85, at=3)
		mtext('>10', side=1, line=0.85, at=4)
		par(xpd=TRUE)
	dev.off()
}



# ATH all samples
make_Boxplot_All_Thresholds_Labels(ATH_cd_nc_SAS_cor_wo_pollen_0.5_2_pearson, ATH_cd_nc_SAS_cor_wo_pollen_2_5_pearson, 
	ATH_cd_nc_SAS_cor_wo_pollen_5_10_pearson, ATH_cd_nc_SAS_cor_wo_pollen_10_pearson, samples = "all")

# ATH comparative samples
make_Boxplot_All_Thresholds(ATH_comp_samples_cd_nc_SAS_cor_wo_pollen_0.5_2_pearson, ATH_comp_samples_cd_nc_SAS_cor_wo_pollen_2_5_pearson, 
	ATH_comp_samples_cd_nc_SAS_cor_wo_pollen_5_10_pearson, ATH_comp_samples_cd_nc_SAS_cor_wo_pollen_10_pearson, 
	samples = "comparative")

# AL comparative samples
make_Boxplot_All_Thresholds(AL_comp_samples_cd_nc_SAS_cor_wo_pollen_0.5_2_pearson, AL_comp_samples_cd_nc_SAS_cor_wo_pollen_2_5_pearson, 
	AL_comp_samples_cd_nc_SAS_cor_wo_pollen_5_10_pearson, AL_comp_samples_cd_nc_SAS_cor_wo_pollen_10_pearson)

# CR
make_Boxplot_All_Thresholds(CR_cd_nc_SAS_cor_wo_pollen_0.5_2_pearson, CR_cd_nc_SAS_cor_wo_pollen_2_5_pearson, 
	CR_cd_nc_SAS_cor_wo_pollen_5_10_pearson, CR_cd_nc_SAS_cor_wo_pollen_10_pearson)

# ES
make_Boxplot_All_Thresholds_Labels(ES_cd_nc_SAS_cor_wo_pollen_0.5_2_pearson, ES_cd_nc_SAS_cor_wo_pollen_2_5_pearson, 
	ES_cd_nc_SAS_cor_wo_pollen_5_10_pearson, ES_cd_nc_SAS_cor_wo_pollen_10_pearson)

# TH
make_Boxplot_All_Thresholds(TH_cd_nc_SAS_cor_wo_pollen_0.5_2_pearson, TH_cd_nc_SAS_cor_wo_pollen_2_5_pearson, 
	TH_cd_nc_SAS_cor_wo_pollen_5_10_pearson, TH_cd_nc_SAS_cor_wo_pollen_10_pearson)

# MT
make_Boxplot_All_Thresholds(MT_cd_nc_SAS_cor_wo_pollen_0.5_2_pearson, MT_cd_nc_SAS_cor_wo_pollen_2_5_pearson, 
	MT_cd_nc_SAS_cor_wo_pollen_5_10_pearson, MT_cd_nc_SAS_cor_wo_pollen_10_pearson)

# BD
make_Boxplot_All_Thresholds(BD_cd_nc_SAS_cor_wo_pollen_0.5_2_pearson, BD_cd_nc_SAS_cor_wo_pollen_2_5_pearson, 
	BD_cd_nc_SAS_cor_wo_pollen_5_10_pearson, BD_cd_nc_SAS_cor_wo_pollen_10_pearson)




#--------------- ATGE/DevSeq, NAT_length and Spearman - Pearson cor scatter plots ---------------


# Correlation plots of cd-cd SAS / nc-cd SAS pairs ATH_all_samples in DevSeq and ATGE
DevSeq_pearson <- length(ATH_cd_nc_SAS_cor_wo_pollen_0.5_pearson)
ATGE_pearson <- length(ATH_cd_nc_SAS_cor_wo_pollen_0.5_in_ATGE_pearson)
DevSeq_spearman <- length(ATH_cd_nc_SAS_cor_wo_pollen_0.5_spearman)
ATGE_spearman <- length(ATH_cd_nc_SAS_cor_wo_pollen_0.5_in_ATGE_spearman)

# Pearson boxplot
jpeg(file=file.path(out_dir, "output", "plots", "nccd_SAS_pearson_ATH_all_DevSeq_all.jpeg"), 
	width = 2620, height = 4000, res = 825)
par(mar = c(5.725, 5.92, 4, 1))
boxplot(ATH_cd_nc_SAS_cor_wo_pollen_0.5_pearson, ATH_cd_nc_SAS_cor_wo_pollen_0.5_in_ATGE_pearson, 
	ylim = c(-1.0, 1.2), 
	names = FALSE, 
	xaxt = 'n', 
	yaxt = 'n', 
	cex.lab = 1.1, 
	las = 2,
	cex.axis = 1.1, #adapt size of axis labels
	xlab = "", 
	ylab = "", 
	col = c("#d8a900", "#d8a900"), 
	boxwex = 0.8, 
	lwd = 1.35, 
	whisklty = 1, 
	at = c(1,2), 
	pars = list(outcol = "gray50"), 
	notch = FALSE
	)
	title("nc-cd SAS pairs", adj = 0.50, line = 1.25, font.main = 1, cex.main = 1.2)
	title(xlab = "present in data set", line = 2.65, cex.lab = 1.1)
	title(ylab = "Pearson ρ", line = 3.0, cex.lab = 1.1)
	rug(x = c(1,2), ticksize = -0.035, side = 1, lwd = 1.35, col = "black") #x-axis ticks
	box(lwd = 1.35)
	axis(side = 2, lwd = 1.35, las = 2)
	text(x= 1.5, y = 1.135, labels= "p < 0.01", col = "black", cex = 1) #ATH_all p-value
	text(x= 1, y= -0.95, labels= DevSeq_pearson, col= "gray40", cex= 0.97) #ATH_all no.genes
	text(x= 2, y= -0.95, labels= ATGE_pearson, col= "gray40", cex= 0.97)
	mtext('DevSeq', side = 1, line = 0.85, at = 1)
	mtext('ATGE', side = 1, line = 0.85, at = 2)
	par(xpd=TRUE)
dev.off()



# Compute rsqr cor value
testRsq <- function(x, y) { 
	test <- cor(x, y) ^ 2
  	test <- round(test, digits=2)
  	return(test)
}

ATH_all_cor <- testRsq(ATH_cd_nc_SAS_cor_wo_pollen_0.5_pearson, ATH_cd_nc_SAS_cor_wo_pollen_0.5_spearman)
ATH_comp_cor <- testRsq(ATH_comp_samples_cd_nc_SAS_cor_wo_pollen_0.5_pearson, ATH_comp_samples_cd_nc_SAS_cor_wo_pollen_0.5_spearman)
AL_comp_cor <- testRsq(AL_comp_samples_cd_nc_SAS_cor_wo_pollen_0.5_pearson, AL_comp_samples_cd_nc_SAS_cor_wo_pollen_0.5_spearman)
CR_cor <- testRsq(CR_cd_nc_SAS_cor_wo_pollen_0.5_pearson, CR_cd_nc_SAS_cor_wo_pollen_0.5_spearman)
ES_cor <- testRsq(ES_cd_nc_SAS_cor_wo_pollen_0.5_pearson, ES_cd_nc_SAS_cor_wo_pollen_0.5_spearman)
TH_cor <- testRsq(TH_cd_nc_SAS_cor_wo_pollen_0.5_pearson, TH_cd_nc_SAS_cor_wo_pollen_0.5_spearman)
MT_cor <- testRsq(MT_cd_nc_SAS_cor_wo_pollen_0.5_pearson, MT_cd_nc_SAS_cor_wo_pollen_0.5_spearman)
BD_cor <- testRsq(BD_cd_nc_SAS_cor_wo_pollen_0.5_pearson, BD_cd_nc_SAS_cor_wo_pollen_0.5_spearman)



# Function to prepare data for nc-cd SAS pair overlap length in relation to pearson correlation
list_for_overlap <- list(
	ATH_cd_nc_SAS_wo_pollen_0.5_cor_length = ATH_cd_nc_SAS_cor_wo_pollen_0.5,
	ATH_comp_samples_cd_nc_SAS_wo_pollen_0.5_cor_length = ATH_comparative_samples_cd_nc_SAS_cor_wo_pollen_0.5,
	AL_cd_nc_SAS_wo_pollen_0.5_cor_length = AL_comparative_samples_cd_nc_SAS_cor_wo_pollen_0.5,
	BD_cd_nc_SAS_wo_pollen_0.5_cor_length = BD_cd_nc_SAS_cor_wo_pollen_0.5,
	CR_cd_nc_SAS_wo_pollen_0.5_cor_length = CR_cd_nc_SAS_cor_wo_pollen_0.5,
	ES_cd_nc_SAS_wo_pollen_0.5_cor_length = ES_cd_nc_SAS_cor_wo_pollen_0.5,
	MT_cd_nc_SAS_wo_pollen_0.5_cor_length = MT_cd_nc_SAS_cor_wo_pollen_0.5,
	TH_cd_nc_SAS_wo_pollen_0.5_cor_length = TH_cd_nc_SAS_cor_wo_pollen_0.5)

getPearsonPercOverlap <- function(x) {

	plus_strand_overlap <- x %>% select(id_plus_strand, width_query, biotype_query, Spearman, 
		Pearson, NAT_overlap_width)
	plus_strand_NAT_overlap <- subset(plus_strand_overlap, 
		biotype_query == "lnc_exonic_antisense" | biotype_query == "lnc_intronic_antisense")
	names(plus_strand_NAT_overlap) <- c("id", "width", "biotype", "Spearman", "Pearson", "overlap")

	minus_strand_overlap <- x %>% select(id_minus_strand, width_subject, biotype_subject, 
		Spearman, Pearson, NAT_overlap_width)
	minus_strand_NAT_overlap <- subset(minus_strand_overlap, 
		biotype_subject == "lnc_exonic_antisense" | biotype_subject == "lnc_intronic_antisense")
	names(minus_strand_NAT_overlap) <- c("id", "width", "biotype", "Spearman", "Pearson", "overlap")

	plus_minus_NAT_overlap <- rbind(plus_strand_NAT_overlap, minus_strand_NAT_overlap)

	plus_minus_NAT_overlap$percent_overlap <- (
		plus_minus_NAT_overlap$overlap / plus_minus_NAT_overlap$width) * 100

	return(plus_minus_NAT_overlap)
}

percent_overlap_pearson_list <- lapply(list_for_overlap, getPearsonPercOverlap)
list2env(percent_overlap_pearson_list, envir = .GlobalEnv)


# Function to prepare data frame and encode data density as color
scatterDensity <- function(x, y) {
	plot_data <- data.frame(x, y)
	names(plot_data) <- c("x_data", "y_data")
	
	# Use densCols() output to get density at each point
	plot_data$col <- densCols(x, y, colramp=colorRampPalette(c("#000099", "#00FEFF", "#45FE4F", 
                            "#FCFF00", "#FF9400", "#FF3100")))
	
	# Reorder "plot_data" based on "col" values - the highest density points are plotted on top
	plot_data <- plot_data[order(plot_data$col),]
	plot_data <- na.omit(plot_data)

	return(plot_data)
}


# Apply scatterDensity function
# for percent overlap vs pearson plots
perc_overlap_ATH_all <- scatterDensity(ATH_cd_nc_SAS_wo_pollen_0.5_cor_length$Pearson, ATH_cd_nc_SAS_wo_pollen_0.5_cor_length$percent_overlap)
perc_overlap_ATH_comp <- scatterDensity(ATH_comp_samples_cd_nc_SAS_wo_pollen_0.5_cor_length$Pearson, ATH_comp_samples_cd_nc_SAS_wo_pollen_0.5_cor_length$percent_overlap)
perc_overlap_AL <- scatterDensity(AL_cd_nc_SAS_wo_pollen_0.5_cor_length$Pearson, AL_cd_nc_SAS_wo_pollen_0.5_cor_length$percent_overlap)
perc_overlap_CR <- scatterDensity(CR_cd_nc_SAS_wo_pollen_0.5_cor_length$Pearson, CR_cd_nc_SAS_wo_pollen_0.5_cor_length$percent_overlap)
perc_overlap_ES <- scatterDensity(ES_cd_nc_SAS_wo_pollen_0.5_cor_length$Pearson, ES_cd_nc_SAS_wo_pollen_0.5_cor_length$percent_overlap)
perc_overlap_TH <- scatterDensity(TH_cd_nc_SAS_wo_pollen_0.5_cor_length$Pearson, TH_cd_nc_SAS_wo_pollen_0.5_cor_length$percent_overlap)
perc_overlap_MT <- scatterDensity(MT_cd_nc_SAS_wo_pollen_0.5_cor_length$Pearson, MT_cd_nc_SAS_wo_pollen_0.5_cor_length$percent_overlap)
perc_overlap_BD <- scatterDensity(BD_cd_nc_SAS_wo_pollen_0.5_cor_length$Pearson, BD_cd_nc_SAS_wo_pollen_0.5_cor_length$percent_overlap)

# for absolute overlap vs pearson plot
abs_overlap_ATH_all <- scatterDensity(ATH_cd_nc_SAS_wo_pollen_0.5_cor_length$Pearson, ATH_cd_nc_SAS_wo_pollen_0.5_cor_length$overlap)
abs_overlap_ATH_comp <- scatterDensity(ATH_comp_samples_cd_nc_SAS_wo_pollen_0.5_cor_length$Pearson, ATH_comp_samples_cd_nc_SAS_wo_pollen_0.5_cor_length$overlap)
abs_overlap_AL <- scatterDensity(AL_cd_nc_SAS_wo_pollen_0.5_cor_length$Pearson, AL_cd_nc_SAS_wo_pollen_0.5_cor_length$overlap)
abs_overlap_CR <- scatterDensity(CR_cd_nc_SAS_wo_pollen_0.5_cor_length$Pearson, CR_cd_nc_SAS_wo_pollen_0.5_cor_length$overlap)
abs_overlap_ES <- scatterDensity(ES_cd_nc_SAS_wo_pollen_0.5_cor_length$Pearson, ES_cd_nc_SAS_wo_pollen_0.5_cor_length$overlap)
abs_overlap_TH <- scatterDensity(TH_cd_nc_SAS_wo_pollen_0.5_cor_length$Pearson, TH_cd_nc_SAS_wo_pollen_0.5_cor_length$overlap)
abs_overlap_MT <- scatterDensity(MT_cd_nc_SAS_wo_pollen_0.5_cor_length$Pearson, MT_cd_nc_SAS_wo_pollen_0.5_cor_length$overlap)
abs_overlap_BD <- scatterDensity(BD_cd_nc_SAS_wo_pollen_0.5_cor_length$Pearson, BD_cd_nc_SAS_wo_pollen_0.5_cor_length$overlap)

# for pearson vs spearman plots
ATH_all_PS <- scatterDensity(ATH_cd_nc_SAS_cor_wo_pollen_0.5_pearson, ATH_cd_nc_SAS_cor_wo_pollen_0.5_spearman)
ATH_comp_PS <- scatterDensity(ATH_comp_samples_cd_nc_SAS_cor_wo_pollen_0.5_pearson, ATH_comp_samples_cd_nc_SAS_cor_wo_pollen_0.5_spearman)
AL_comp_PS <- scatterDensity(AL_comp_samples_coding_SAS_cor_wo_pollen_pearson, AL_comp_samples_coding_SAS_cor_wo_pollen_spearman)
CR_comp_PS <- scatterDensity(CR_cd_nc_SAS_cor_wo_pollen_0.5_pearson, CR_cd_nc_SAS_cor_wo_pollen_0.5_spearman)
ES_comp_PS <- scatterDensity(ES_cd_nc_SAS_cor_wo_pollen_0.5_pearson, ES_cd_nc_SAS_cor_wo_pollen_0.5_spearman)
TH_comp_PS <- scatterDensity(TH_cd_nc_SAS_cor_wo_pollen_0.5_pearson, TH_cd_nc_SAS_cor_wo_pollen_0.5_spearman)
MT_comp_PS <- scatterDensity(MT_cd_nc_SAS_cor_wo_pollen_0.5_pearson, MT_cd_nc_SAS_cor_wo_pollen_0.5_spearman)
BD_comp_PS <- scatterDensity(BD_cd_nc_SAS_cor_wo_pollen_0.5_pearson, BD_cd_nc_SAS_cor_wo_pollen_0.5_spearman)

# Calculate R squared value for relative overlap vs pearson data
rsqd_ATH_all_perc <- testRsq(ATH_cd_nc_SAS_wo_pollen_0.5_cor_length$Pearson, ATH_cd_nc_SAS_wo_pollen_0.5_cor_length$percent_overlap)
rsqd_ATH_comp_perc <- testRsq(ATH_comp_samples_cd_nc_SAS_wo_pollen_0.5_cor_length$Pearson, ATH_comp_samples_cd_nc_SAS_wo_pollen_0.5_cor_length$percent_overlap)
rsqd_AL_perc <- testRsq(AL_cd_nc_SAS_wo_pollen_0.5_cor_length$Pearson, AL_cd_nc_SAS_wo_pollen_0.5_cor_length$percent_overlap)
rsqd_CR_perc <- testRsq(CR_cd_nc_SAS_wo_pollen_0.5_cor_length$Pearson, CR_cd_nc_SAS_wo_pollen_0.5_cor_length$percent_overlap)
rsqd_ES_perc <- testRsq(ES_cd_nc_SAS_wo_pollen_0.5_cor_length$Pearson, ES_cd_nc_SAS_wo_pollen_0.5_cor_length$percent_overlap)
rsqd_TH_perc <- testRsq(TH_cd_nc_SAS_wo_pollen_0.5_cor_length$Pearson, TH_cd_nc_SAS_wo_pollen_0.5_cor_length$percent_overlap)
rsqd_MT_perc <- testRsq(MT_cd_nc_SAS_wo_pollen_0.5_cor_length$Pearson, MT_cd_nc_SAS_wo_pollen_0.5_cor_length$percent_overlap)
rsqd_BD_perc <- testRsq(BD_cd_nc_SAS_wo_pollen_0.5_cor_length$Pearson, BD_cd_nc_SAS_wo_pollen_0.5_cor_length$percent_overlap)


# Compute adjusted rsqr value
testAdjRsq <- function(x,y) { 
	test <- summary(gam(y ~ x))$r.sq
	test <- round(test, digits=2)
  	return(test)
}

# Calculate R squared value for absolute overlap vs pearson data
rsqd_ATH_all_abs <- testAdjRsq(ATH_cd_nc_SAS_wo_pollen_0.5_cor_length$Pearson, ATH_cd_nc_SAS_wo_pollen_0.5_cor_length$overlap)
rsqd_ATH_comp_abs <- testAdjRsq(ATH_comp_samples_cd_nc_SAS_wo_pollen_0.5_cor_length$Pearson, ATH_comp_samples_cd_nc_SAS_wo_pollen_0.5_cor_length$overlap)
rsqd_AL_abs <- testAdjRsq(AL_cd_nc_SAS_wo_pollen_0.5_cor_length$Pearson, AL_cd_nc_SAS_wo_pollen_0.5_cor_length$overlap)
rsqd_CR_abs <- testAdjRsq(CR_cd_nc_SAS_wo_pollen_0.5_cor_length$Pearson, CR_cd_nc_SAS_wo_pollen_0.5_cor_length$overlap)
rsqd_ES_abs <- testAdjRsq(ES_cd_nc_SAS_wo_pollen_0.5_cor_length$Pearson, ES_cd_nc_SAS_wo_pollen_0.5_cor_length$overlap)
rsqd_TH_abs <- testAdjRsq(TH_cd_nc_SAS_wo_pollen_0.5_cor_length$Pearson, TH_cd_nc_SAS_wo_pollen_0.5_cor_length$overlap)
rsqd_MT_abs <- testAdjRsq(MT_cd_nc_SAS_wo_pollen_0.5_cor_length$Pearson, MT_cd_nc_SAS_wo_pollen_0.5_cor_length$overlap)
rsqd_BD_abs <- testAdjRsq(BD_cd_nc_SAS_wo_pollen_0.5_cor_length$Pearson, BD_cd_nc_SAS_wo_pollen_0.5_cor_length$overlap)

# Calculate R squared value for pearson vs spearman data
rsqd_ATH_all_PS <- testRsq(ATH_cd_nc_SAS_cor_wo_pollen_0.5_pearson, ATH_cd_nc_SAS_cor_wo_pollen_0.5_spearman)
rsqd_ATH_comp_PS <- testRsq(ATH_comp_samples_cd_nc_SAS_cor_wo_pollen_0.5_pearson, ATH_comp_samples_cd_nc_SAS_cor_wo_pollen_0.5_spearman)
rsqd_AL_PS <- testRsq(AL_comp_samples_coding_SAS_cor_wo_pollen_pearson, AL_comp_samples_coding_SAS_cor_wo_pollen_spearman)
rsqd_CR_PS <- testRsq(CR_cd_nc_SAS_cor_wo_pollen_0.5_pearson, CR_cd_nc_SAS_cor_wo_pollen_0.5_spearman)
rsqd_ES_PS <- testRsq(ES_cd_nc_SAS_cor_wo_pollen_0.5_pearson, ES_cd_nc_SAS_cor_wo_pollen_0.5_spearman)
rsqd_TH_PS <- testRsq(TH_cd_nc_SAS_cor_wo_pollen_0.5_pearson, TH_cd_nc_SAS_cor_wo_pollen_0.5_spearman)
rsqd_MT_PS <- testRsq(MT_cd_nc_SAS_cor_wo_pollen_0.5_pearson, MT_cd_nc_SAS_cor_wo_pollen_0.5_spearman)
rsqd_BD_PS <- testRsq(BD_cd_nc_SAS_cor_wo_pollen_0.5_pearson, BD_cd_nc_SAS_cor_wo_pollen_0.5_spearman)



## Some tests to check if assumptions of linear model are met
## Examples for ATH data
# abs_lm_mod <- lm(abs_overlap_ATH_all$y_data ~ abs_overlap_ATH_all$x_data)
# summary(abs_lm_mod)
# plot(abs_lm_mod)
# perc_lm_mod <- lm(perc_overlap_ATH_all$y_data ~ perc_overlap_ATH_all$x_data)
# summary(perc_lm_mod)
# plot(perc_lm_mod)

## Check if assumptions of gam or glm gamma are met
# library(DHARMa)
# abs_gam_mod <- gam(abs_overlap_ATH_all$y_data ~ abs_overlap_ATH_all$x_data)
# summary(abs_gam_mod)
# simulationOutput <- simulateResiduals(fittedModel = abs_gam_mod)
# plot(simulationOutput)
# abs_gamma_mod <- glm(abs_overlap_ATH_all$y_data ~ abs_overlap_ATH_all$x_data, 
# 	family = Gamma(link = "log"))
# summary(abs_gamma_mod)
# simulationOutput <- simulateResiduals(fittedModel = abs_gamma_mod)
# plot(simulationOutput)



# Function to scatter plot relative overlap (%) versus pearson correlation

makeScrPlotRelOverlap <- function(data, rsqd, plot_title = c(
	"ATH_all", "ATH", "CR", "ES", "TH", "MT", "BD"), rsgd_pos, vjust_1, vjust_2) {

	fname <- sprintf('%s.jpg', paste(deparse(substitute(data)), sep="_"))
	
	rsrt_label = paste("R ^ 2", "==", ".")
	p <- ggplot(data, aes(x = x_data, y = y_data)) + 
	geom_point(size = 1.5, colour = data$col) + 
	scale_x_continuous(limits = c(-1.02,1.02), breaks=c(-1,-0.5,0,0.5,1), labels=c(-1,-0.5,0,0.5,1), expand = c(0, 0)) +
	scale_y_continuous(limits = c(0,101), expand = c(0, 0)) + 
	annotate("text", x = -Inf, y = Inf, hjust = -0.31, vjust = vjust_1, size=5.7, label = rsrt_label, parse = TRUE) + 
	annotate("text", x = -Inf, y = Inf, hjust = rsgd_pos, vjust = vjust_2, size=5.7, label = rsqd, parse = FALSE) 
	q <- p + ggtitle(plot_title) + theme_bw() + xlab("Pearson") + ylab("NAT overlap (%)") + 
		theme(text=element_text(size=16), 
		axis.ticks.length = unit(.3, "cm"),
		plot.margin = unit(c(3.0, 10.5, 45.5, 17), "points"),
		axis.text.x = element_text(colour = "black", size=16, angle=0, margin = margin(t = 5.25, r = 0, b = 0, l = 0)), 
		axis.text.y = element_text(colour = "black", size=16, angle=0, margin = margin(t = 0, r = 5.25, b = 0, l = 0)),
		axis.title.x = element_text(colour = "black", size=17.5, margin = margin(t = 14.5, r = 0, b = 1, l = 0)),
		axis.title.y = element_text(colour = "black", size=17.5, margin = margin(t = 0, r = 9, b = 0, l = 1)),
		plot.title = element_text(colour = "black", size=17.5, margin = margin(t = 19, r = 0, b = 8, l = 0), hjust = 0.5),
		legend.position = "bottom",
		panel.border = element_rect(colour = "black", fill=NA, size=0.5))

	ggsave(file = file.path(out_dir, "output", "plots", fname), plot = q,
		scale = 1, width = 5, height = 5.69697, units = c("in"), 
		dpi = 825, limitsize = FALSE)
}


makeScrPlotRelOverlap(data=perc_overlap_ATH_all, rsqd=rsqd_ATH_all_perc, plot_title="ATH_all", rsgd_pos= -0.405, vjust_1=2.9, vjust_2=5.5)
makeScrPlotRelOverlap(data=perc_overlap_ATH_comp, rsqd=rsqd_ATH_comp_perc, plot_title="ATH_comp", rsgd_pos= -0.405, vjust_1=10.5, vjust_2=16.55)
makeScrPlotRelOverlap(data=perc_overlap_AL, rsqd=rsqd_AL_perc, plot_title="AL_", rsgd_pos= -1.41, vjust_1=7.55, vjust_2=12)
makeScrPlotRelOverlap(data=perc_overlap_CR, rsqd=rsqd_CR_perc, plot_title="CR_", rsgd_pos= -0.405, vjust_1=3.8, vjust_2=6.8)
makeScrPlotRelOverlap(data=perc_overlap_ES, rsqd=rsqd_ES_perc, plot_title="ES_", rsgd_pos= -0.405, vjust_1=12.68, vjust_2=19.1)
makeScrPlotRelOverlap(data=perc_overlap_TH, rsqd=rsqd_TH_perc, plot_title="TH_", rsgd_pos= -0.405, vjust_1=7.0, vjust_2=11.25)
makeScrPlotRelOverlap(data=perc_overlap_MT, rsqd=rsqd_MT_perc, plot_title="MT_", rsgd_pos= -0.405, vjust_1=2, vjust_2=4.33)
makeScrPlotRelOverlap(data=perc_overlap_BD, rsqd=rsqd_BD_perc, plot_title="BD_", rsgd_pos= -0.405, vjust_1=1.8, vjust_2=4)




# Function to scatter plot absolute overlap (bp) versus pearson correlation

makeScrPlotAbsOverlap <- function(data, rsqd, plot_title = c(
	"ATH_all", "ATH", "CR", "ES", "TH", "MT", "BD")) {

	fname <- sprintf('%s.jpg', paste(deparse(substitute(data)), sep="_"))
	
	rsrt_label = paste("R ^ 2"," == ", rsqd)
	p <- ggplot(data, aes(x = x_data, y = y_data)) + 
	geom_point(size = 1.5, colour = data$col) + 
	scale_x_continuous(limits = c(-1.02,1.02), breaks=c(-1,-0.5,0,0.5,1), labels=c(-1,-0.5,0,0.5,1), expand = c(0, 0)) +
	scale_y_continuous(trans='log10', labels = prettyNum, breaks=c(1,10,100,1000,10000), limits=c(1, 22000), expand = c(0, 0)) + 
	geom_smooth(method="auto" , color="gray20", fill="#69b3a2", se=TRUE, size=1) +  # use gam regression model
	annotate("text", x = -Inf, y = Inf, hjust = -0.31, vjust = 1.6, size=5.7, label = rsrt_label, parse = TRUE)
	q <- p + ggtitle(plot_title) + theme_bw() + xlab("Pearson") + ylab("NAT overlap (bp)") + 
  		theme(text=element_text(size=16), 
  		axis.ticks.length = unit(.3, "cm"),
  		plot.margin = unit(c(3.0, 10.5, 45.5, 8), "points"),
		axis.text.x = element_text(colour = "black", size=16, angle=0, margin = margin(t = 5.25, r = 0, b = 0, l = 0)), 
		axis.text.y = element_text(colour = "black", size=16, angle=0, margin = margin(t = 0, r = 5.25, b = 0, l = 0)),
		axis.title.x = element_text(colour = "black", size=17.5, margin = margin(t = 14.5, r = 0, b = 1, l = 0)),
		axis.title.y = element_text(colour = "black", size=17.5, margin = margin(t = 0, r = 0, b = 0, l = 1)),
		plot.title = element_text(colour = "black", size=17.5, margin = margin(t = 19, r = 0, b = 8, l = 0), hjust = 0.5),
		legend.position = "bottom",
  		panel.border = element_rect(colour = "black", fill=NA, size=0.5))

	ggsave(file = file.path(out_dir, "output", "plots", fname), plot = q,
		scale = 1, width = 5, height = 5.69697, units = c("in"), 
		dpi = 825, limitsize = FALSE)
}


makeScrPlotAbsOverlap(data=abs_overlap_ATH_all, rsqd=rsqd_ATH_all_abs, plot_title="ATH_all")
makeScrPlotAbsOverlap(data=abs_overlap_ATH_comp, rsqd=rsqd_ATH_comp_abs, plot_title="ATH_comp")
makeScrPlotAbsOverlap(data=abs_overlap_AL, rsqd=rsqd_AL_abs, plot_title="AL_")
makeScrPlotAbsOverlap(data=abs_overlap_CR, rsqd=rsqd_CR_abs, plot_title="CR_")
makeScrPlotAbsOverlap(data=abs_overlap_ES, rsqd=rsqd_ES_abs, plot_title="ES_")
makeScrPlotAbsOverlap(data=abs_overlap_TH, rsqd=rsqd_TH_abs, plot_title="TH_")
makeScrPlotAbsOverlap(data=abs_overlap_MT, rsqd=rsqd_MT_abs, plot_title="MT_")
makeScrPlotAbsOverlap(data=abs_overlap_BD, rsqd=rsqd_BD_abs, plot_title="BD_")



# Function to scatter plot pearson versus spearman correlation

makeScrPlotPSCor <- function(data, rsqd, plot_title = c(
	"ATH_all", "ATH", "CR", "ES", "TH", "MT", "BD")) {

	fname <- sprintf('%s.jpg', paste(deparse(substitute(data)), sep="_"))
	
	rsrt_label = paste("R ^ 2"," == ", rsqd)
	p <- ggplot(data, aes(x = x_data, y = y_data)) + 
	geom_point(size = 1.25, colour = data$col) + 
	scale_x_continuous(limits = c(-1.02,1.02), expand = c(0, 0)) +
	scale_y_continuous(limits = c(-1.02,1.02), expand = c(0, 0)) + 
	geom_smooth(method="lm" , color="gray20", fill="#69b3a2", se=TRUE, size=1) +  # use linear regression model
	annotate("text", x = -Inf, y = Inf, hjust = -0.335, vjust = 1.6, size=5.35, label = rsrt_label, parse = TRUE)
	q <- p + ggtitle(plot_title) + theme_bw() + xlab("Pearson") + ylab("Spearman") + 
  		theme(text=element_text(size=16), 
  		axis.ticks.length = unit(.3, "cm"),
  		plot.margin = unit(c(3.0, 10.5, 42.5, 5.5), "points"),
  		axis.text.x = element_text(colour = "black", size=14.25, angle=0, margin = margin(t = 8.25, r = 0, b = 0, l = 0)), 
  		axis.text.y = element_text(colour = "black", size=14.25, angle=0, margin = margin(t = 0, r = 8.25, b = 0, l = 0)),
  		axis.title.x = element_text(colour = "black", margin = margin(t = 14.5, r = 0, b = 1, l = 0)),
  		axis.title.y = element_text(colour = "black", margin = margin(t = 0, r = 9, b = 0, l = 1)),
  		plot.title = element_text(colour = "black", size=17, margin = margin(t = 11.5, r = 0, b = 15.5, l = 0), hjust = 0.5),
  		legend.position = "bottom",
  		panel.border = element_rect(colour = "black", fill=NA, size=1.2))

	ggsave(file = file.path(out_dir, "output", "plots", fname), plot = q,
		scale = 1, width = 4.848485, height = 5.69697, units = c("in"), 
		dpi = 825, limitsize = FALSE)
}


makeScrPlotPSCor(data=ATH_all_PS, rsqd=rsqd_ATH_all_PS, plot_title="ATH_all")
makeScrPlotPSCor(data=ATH_comp_PS, rsqd=rsqd_ATH_comp_PS, plot_title="ATH_comp")
makeScrPlotPSCor(data=AL_comp_PS, rsqd=rsqd_AL_PS, plot_title="AL")
makeScrPlotPSCor(data=CR_comp_PS, rsqd=rsqd_CR_PS, plot_title="CR")
makeScrPlotPSCor(data=ES_comp_PS, rsqd=rsqd_ES_PS, plot_title="ES")
makeScrPlotPSCor(data=TH_comp_PS, rsqd=rsqd_TH_PS, plot_title="TH")
makeScrPlotPSCor(data=MT_comp_PS, rsqd=rsqd_MT_PS, plot_title="MT")
makeScrPlotPSCor(data=BD_comp_PS, rsqd=rsqd_BD_PS, plot_title="BD")











