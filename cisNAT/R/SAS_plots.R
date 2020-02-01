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
if (!require(ggplot2)) install.packages('dplyr')
library(ggplot2)


# Set file path and input files
in_dir_cd <- "/Volumes/User/Shared/Christoph_manuscript/DevSeq_paper/Analysis/Analysis_2019/A_thaliana_gene_exression_map/20191121_CS_coding_cisNAT_analysis/output/overlap_cd_genes"
in_dir_nc <- "/Volumes/User/Shared/Christoph_manuscript/DevSeq_paper/Analysis/Analysis_2019/A_thaliana_gene_exression_map/20191121_CS_coding_cisNAT_analysis/output/overlap_nc_genes"
in_dir_ATGE <- "/Volumes/User/Shared/Christoph_manuscript/DevSeq_paper/Analysis/Analysis_2019/A_thaliana_gene_exression_map/20191121_CS_coding_cisNAT_analysis/output/SAS_DevSeq_ATGE"
out_dir <- "/Volumes/User/Shared/Christoph_manuscript/DevSeq_paper/Analysis/Analysis_2019/A_thaliana_gene_exression_map/20191121_CS_coding_cisNAT_analysis"


# Read all csv files in input file path
readTable <- function(path, pattern = "*.csv") {
    files = list.files(path, pattern, full.names = TRUE)
    lapply(files, function(x) read.table(x, sep=";", dec=".", header=TRUE, stringsAsFactors=FALSE))
}

coding_genes_tables <- readTable(in_dir_cd)
NAT_genes_tables <- readTable(in_dir_nc)


# Get file names and save them in character vector
coding_genes_table_list <- as.character(list.files(in_dir_cd, pattern = "*.csv"))
coding_genes_table_names <- gsub('\\.csv$', '', coding_genes_table_list)

NAT_genes_table_list <- as.character(list.files(in_dir_nc, pattern = "*.csv"))
NAT_genes_table_names <- gsub('\\.csv$', '', NAT_genes_table_list)


# Change data frame names in list
names(coding_genes_tables) <- coding_genes_table_names
list2env(coding_genes_tables, envir = .GlobalEnv)

names(NAT_genes_tables) <- NAT_genes_table_names
list2env(NAT_genes_tables, envir = .GlobalEnv)


# Read ATGE_NAT ID table
ATGE_NAT_ID <- read.table(file=file.path(in_dir_ATGE, "ATH_cd_nc_SAS_cor_wo_pollen_0.5_in_ATGE.csv"), 
	sep=";", dec=".", header=TRUE, stringsAsFactors = FALSE)




#----------------------- Define TPM ranges for all species (w/o pollen)  -----------------------


# ATH all samples
ATH_cd_nc_SAS_cor_wo_pollen_05_2 <- ATH_cd_nc_SAS_cor_wo_pollen_0.5 %>% filter(
	!((id_plus_strand %in% ATH_cd_nc_SAS_cor_wo_pollen_2$id_plus_strand) &
		(id_minus_strand %in% ATH_cd_nc_SAS_cor_wo_pollen_2$id_minus_strand)))
ATH_cd_nc_SAS_cor_wo_pollen_2_5 <- ATH_cd_nc_SAS_cor_wo_pollen_2 %>% filter(
	!((id_plus_strand %in% ATH_cd_nc_SAS_cor_wo_pollen_5$id_plus_strand)
	 & (id_minus_strand %in% ATH_cd_nc_SAS_cor_wo_pollen_5$id_minus_strand)))
# ATH comparative samples
ATH_comparative_samples_cd_nc_SAS_cor_wo_pollen_0.5_2 <- ATH_comparative_samples_cd_nc_SAS_cor_wo_pollen_0.5 %>% filter(
	!((id_plus_strand %in% ATH_comparative_samples_cd_nc_SAS_cor_wo_pollen_2$id_plus_strand) &
		(id_minus_strand %in% ATH_comparative_samples_cd_nc_SAS_cor_wo_pollen_2$id_minus_strand)))
ATH_comparative_samples_cd_nc_SAS_cor_wo_pollen_2_5 <- ATH_comparative_samples_cd_nc_SAS_cor_wo_pollen_2 %>% filter(
	!((id_plus_strand %in% ATH_comparative_samples_cd_nc_SAS_cor_wo_pollen_5$id_plus_strand)
	 & (id_minus_strand %in% ATH_comparative_samples_cd_nc_SAS_cor_wo_pollen_5$id_minus_strand)))
# AL comparative samples
AL_comparative_samples_cd_nc_SAS_cor_wo_pollen_0.5_2 <- AL_comparative_samples_cd_nc_SAS_cor_wo_pollen_0.5 %>% filter(
	!((id_plus_strand %in% AL_comparative_samples_cd_nc_SAS_cor_wo_pollen_2$id_plus_strand) &
		(id_minus_strand %in% AL_comparative_samples_cd_nc_SAS_cor_wo_pollen_2$id_minus_strand)))
AL_comparative_samples_cd_nc_SAS_cor_wo_pollen_2_5 <- AL_comparative_samples_cd_nc_SAS_cor_wo_pollen_2 %>% filter(
	!((id_plus_strand %in% AL_comparative_samples_cd_nc_SAS_cor_wo_pollen_5$id_plus_strand)
	 & (id_minus_strand %in% AL_comparative_samples_cd_nc_SAS_cor_wo_pollen_5$id_minus_strand)))
# CR
CR_cd_nc_SAS_cor_wo_pollen_0.5_2 <- CR_cd_nc_SAS_cor_wo_pollen_0.5 %>% filter(
	!((id_plus_strand %in% CR_cd_nc_SAS_cor_wo_pollen_2$id_plus_strand) &
		(id_minus_strand %in% CR_cd_nc_SAS_cor_wo_pollen_2$id_minus_strand)))
CR_cd_nc_SAS_cor_wo_pollen_2_5 <- CR_cd_nc_SAS_cor_wo_pollen_2 %>% filter(
	!((id_plus_strand %in% CR_cd_nc_SAS_cor_wo_pollen_5$id_plus_strand)
	 & (id_minus_strand %in% CR_cd_nc_SAS_cor_wo_pollen_5$id_minus_strand)))
# ES
ES_cd_nc_SAS_cor_wo_pollen_0.5_2 <- ES_cd_nc_SAS_cor_wo_pollen_0.5 %>% filter(
	!((id_plus_strand %in% ES_cd_nc_SAS_cor_wo_pollen_2$id_plus_strand) &
		(id_minus_strand %in% ES_cd_nc_SAS_cor_wo_pollen_2$id_minus_strand)))
ES_cd_nc_SAS_cor_wo_pollen_2_5 <- ES_cd_nc_SAS_cor_wo_pollen_2 %>% filter(
	!((id_plus_strand %in% ES_cd_nc_SAS_cor_wo_pollen_5$id_plus_strand)
	 & (id_minus_strand %in% ES_cd_nc_SAS_cor_wo_pollen_5$id_minus_strand)))
# TH
TH_cd_nc_SAS_cor_wo_pollen_0.5_2 <- TH_cd_nc_SAS_cor_wo_pollen_0.5 %>% filter(
	!((id_plus_strand %in% TH_cd_nc_SAS_cor_wo_pollen_2$id_plus_strand) &
		(id_minus_strand %in% TH_cd_nc_SAS_cor_wo_pollen_2$id_minus_strand)))
TH_cd_nc_SAS_cor_wo_pollen_2_5 <- TH_cd_nc_SAS_cor_wo_pollen_2 %>% filter(
	!((id_plus_strand %in% TH_cd_nc_SAS_cor_wo_pollen_5$id_plus_strand)
	 & (id_minus_strand %in% TH_cd_nc_SAS_cor_wo_pollen_5$id_minus_strand)))
# MT
MT_cd_nc_SAS_cor_wo_pollen_0.5_2 <- MT_cd_nc_SAS_cor_wo_pollen_0.5 %>% filter(
	!((id_plus_strand %in% MT_cd_nc_SAS_cor_wo_pollen_2$id_plus_strand) &
		(id_minus_strand %in% MT_cd_nc_SAS_cor_wo_pollen_2$id_minus_strand)))
MT_cd_nc_SAS_cor_wo_pollen_2_5 <- MT_cd_nc_SAS_cor_wo_pollen_2 %>% filter(
	!((id_plus_strand %in% MT_cd_nc_SAS_cor_wo_pollen_5$id_plus_strand)
	 & (id_minus_strand %in% MT_cd_nc_SAS_cor_wo_pollen_5$id_minus_strand)))
# BD
BD_cd_nc_SAS_cor_wo_pollen_0.5_2 <- BD_cd_nc_SAS_cor_wo_pollen_0.5 %>% filter(
	!((id_plus_strand %in% BD_cd_nc_SAS_cor_wo_pollen_2$id_plus_strand) &
		(id_minus_strand %in% BD_cd_nc_SAS_cor_wo_pollen_2$id_minus_strand)))
BD_cd_nc_SAS_cor_wo_pollen_2_5 <- BD_cd_nc_SAS_cor_wo_pollen_2 %>% filter(
	!((id_plus_strand %in% BD_cd_nc_SAS_cor_wo_pollen_5$id_plus_strand)
	 & (id_minus_strand %in% BD_cd_nc_SAS_cor_wo_pollen_5$id_minus_strand)))




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
	ATH_cd_nc_SAS_cor_wo_pollen_05_2 = ATH_cd_nc_SAS_cor_wo_pollen_05_2,
	ATH_cd_nc_SAS_cor_wo_pollen_2_5 = ATH_cd_nc_SAS_cor_wo_pollen_2_5,
	ATH_cd_nc_SAS_cor_0.5 = ATH_cd_nc_SAS_cor_0.5, #including pollen samples for ATH
	AL_comp_samples_cd_nc_SAS_cor_wo_pollen_0.5 = AL_comparative_samples_cd_nc_SAS_cor_wo_pollen_0.5, 
	AL_comp_samples_cd_nc_SAS_cor_wo_pollen_2 = AL_comparative_samples_cd_nc_SAS_cor_wo_pollen_2, 
	AL_comp_samples_cd_nc_SAS_cor_wo_pollen_5 = AL_comparative_samples_cd_nc_SAS_cor_wo_pollen_5,
	AL_comp_samples_cd_nc_SAS_cor_wo_pollen_0.5_2 = AL_comparative_samples_cd_nc_SAS_cor_wo_pollen_0.5_2,
	AL_comp_samples_cd_nc_SAS_cor_wo_pollen_2_5 = AL_comparative_samples_cd_nc_SAS_cor_wo_pollen_2_5,
	ATH_comp_samples_cd_nc_SAS_cor_wo_pollen_0.5 = ATH_comparative_samples_cd_nc_SAS_cor_wo_pollen_0.5,
	ATH_comp_samples_cd_nc_SAS_cor_wo_pollen_2 = ATH_comparative_samples_cd_nc_SAS_cor_wo_pollen_2, 
	ATH_comp_samples_cd_nc_SAS_cor_wo_pollen_5 = ATH_comparative_samples_cd_nc_SAS_cor_wo_pollen_5,
	ATH_comp_samples_cd_nc_SAS_cor_wo_pollen_0.5_2 = ATH_comparative_samples_cd_nc_SAS_cor_wo_pollen_0.5_2,
	ATH_comp_samples_cd_nc_SAS_cor_wo_pollen_2_5 = ATH_comparative_samples_cd_nc_SAS_cor_wo_pollen_2_5,
	ATH_comp_samples_cd_nc_SAS_cor_0.5 = ATH_comparative_samples_cd_nc_SAS_cor_0.5, #including pollen samples for ATH
	BD_cd_nc_SAS_cor_wo_pollen_0.5 = BD_cd_nc_SAS_cor_wo_pollen_0.5, 
	BD_cd_nc_SAS_cor_wo_pollen_2 = BD_cd_nc_SAS_cor_wo_pollen_2,
	BD_cd_nc_SAS_cor_wo_pollen_5 = BD_cd_nc_SAS_cor_wo_pollen_5,
	BD_cd_nc_SAS_cor_wo_pollen_0.5_2 = BD_cd_nc_SAS_cor_wo_pollen_0.5_2,
	BD_cd_nc_SAS_cor_wo_pollen_2_5 = BD_cd_nc_SAS_cor_wo_pollen_2_5,
	CR_cd_nc_SAS_cor_wo_pollen_0.5 = CR_cd_nc_SAS_cor_wo_pollen_0.5,
	CR_cd_nc_SAS_cor_wo_pollen_2 = CR_cd_nc_SAS_cor_wo_pollen_2,
	CR_cd_nc_SAS_cor_wo_pollen_5 = CR_cd_nc_SAS_cor_wo_pollen_5,
	CR_cd_nc_SAS_cor_wo_pollen_0.5_2 = CR_cd_nc_SAS_cor_wo_pollen_0.5_2,
	CR_cd_nc_SAS_cor_wo_pollen_2_5 = CR_cd_nc_SAS_cor_wo_pollen_2_5,
	ES_cd_nc_SAS_cor_wo_pollen_0.5 = ES_cd_nc_SAS_cor_wo_pollen_0.5,
	ES_cd_nc_SAS_cor_wo_pollen_2 = ES_cd_nc_SAS_cor_wo_pollen_2,
	ES_cd_nc_SAS_cor_wo_pollen_5 = ES_cd_nc_SAS_cor_wo_pollen_5,
	ES_cd_nc_SAS_cor_wo_pollen_0.5_2 = ES_cd_nc_SAS_cor_wo_pollen_0.5_2,
	ES_cd_nc_SAS_cor_wo_pollen_2_5 = ES_cd_nc_SAS_cor_wo_pollen_2_5,
	MT_cd_nc_SAS_cor_wo_pollen_0.5 = MT_cd_nc_SAS_cor_wo_pollen_0.5,
	MT_cd_nc_SAS_cor_wo_pollen_2 = MT_cd_nc_SAS_cor_wo_pollen_2,
	MT_cd_nc_SAS_cor_wo_pollen_5 = MT_cd_nc_SAS_cor_wo_pollen_5,
	MT_cd_nc_SAS_cor_wo_pollen_0.5_2 = MT_cd_nc_SAS_cor_wo_pollen_0.5_2,
	MT_cd_nc_SAS_cor_wo_pollen_2_5 = MT_cd_nc_SAS_cor_wo_pollen_2_5,
	TH_cd_nc_SAS_cor_wo_pollen_0.5 = TH_cd_nc_SAS_cor_wo_pollen_0.5,
	TH_cd_nc_SAS_cor_wo_pollen_2 = TH_cd_nc_SAS_cor_wo_pollen_2,
	TH_cd_nc_SAS_cor_wo_pollen_5 = TH_cd_nc_SAS_cor_wo_pollen_5,
	TH_cd_nc_SAS_cor_wo_pollen_0.5_2 = TH_cd_nc_SAS_cor_wo_pollen_0.5_2,
	TH_cd_nc_SAS_cor_wo_pollen_2_5 = TH_cd_nc_SAS_cor_wo_pollen_2_5)


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


# Write final data tables to csv files and store them in /out_dir/output/data_tables
if (!dir.exists(file.path(out_dir, "output", "plots"))) 
	dir.create(file.path(out_dir, "output", "plots"), recursive = TRUE)




#-------------------------------- Perform Wilcox rank sum test ---------------------------------


wilcox_pearson_cor <- sapply(SAS_pairs_list_pearson, function(x) sapply(
	SAS_pairs_list_pearson, function(y) wilcox.test(x,y)$p.value))
wilcox_spearman_cor <- sapply(SAS_pairs_list_spearman, function(x) sapply(
	SAS_pairs_list_spearman, function(y) wilcox.test(x,y)$p.value))


write.table(wilcox_pearson_cor, file=file.path(out_dir, "output", "plots", "wilcox_pearson_cor.csv"), 
	sep=";", dec=".", row.names=TRUE, col.names=NA)
write.table(wilcox_spearman_cor, file=file.path(out_dir, "output", "plots", "wilcox_spearman_cor.csv"), 
	sep=";", dec=".", row.names=TRUE, col.names=NA)




#---------------------------- Generate pearson plots for figure 1 -----------------------------


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
make_Boxplot_All_Thresholds_Labels <- function(threshold_05, threshold_05_2, threshold_2_5, 
	threshold_greater5, samples=c("all","comparative")) {

	species <- sub("\\_.*", "", deparse(substitute(threshold_05)))
    n_values_05 <- length(threshold_05)
    n_values_05_2 <- length(threshold_05_2)
    n_values_2_5 <- length(threshold_2_5)
    n_values_greater5 <- length(threshold_greater5)

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
	boxplot(threshold_05, threshold_05_2, threshold_2_5, threshold_greater5,
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
		text(x= 1, y= -1.05, labels= n_values_05, col= "gray40", cex=0.97) #threshold_>0.5
		text(x= 2, y= -1.05, labels= n_values_05_2, col= "gray40", cex=0.97) #threshold_0.5-2
		text(x= 3, y= -1.05, labels= n_values_2_5, col= "gray40", cex=0.97) #threshold_2-5
		text(x= 4, y= -1.05, labels= n_values_greater5, col= "gray40", cex=0.97) #threshold_greater5
		mtext('>0.5', side=1, line=0.85, at=1)
		mtext('>0.5-2', side=1, line=0.85, at=2)
		mtext('>2-5', side=1, line=0.85, at=3)
		mtext('>5', side=1, line=0.85, at=4)
		par(xpd=TRUE)
	dev.off()
}


# Pearson plot of nc-cd SAS pairs with all thresholds
make_Boxplot_All_Thresholds <- function(threshold_05, threshold_05_2, threshold_2_5, 
	threshold_greater5, samples=c("all","comparative")) {

	species <- sub("\\_.*", "", deparse(substitute(threshold_05)))
    n_values_05 <- length(threshold_05)
    n_values_05_2 <- length(threshold_05_2)
    n_values_2_5 <- length(threshold_2_5)
    n_values_greater5 <- length(threshold_greater5)

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
	boxplot(threshold_05, threshold_05_2, threshold_2_5, threshold_greater5,
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
		text(x= 1, y= -1.05, labels= n_values_05, col= "gray40", cex=0.97) #threshold_>0.5
		text(x= 2, y= -1.05, labels= n_values_05_2, col= "gray40", cex=0.97) #threshold_0.5-2
		text(x= 3, y= -1.05, labels= n_values_2_5, col= "gray40", cex=0.97) #threshold_2-5
		text(x= 4, y= -1.05, labels= n_values_greater5, col= "gray40", cex=0.97) #threshold_greater5
		mtext('>0.5', side=1, line=0.85, at=1)
		mtext('>0.5-2', side=1, line=0.85, at=2)
		mtext('>2-5', side=1, line=0.85, at=3)
		mtext('>5', side=1, line=0.85, at=4)
		par(xpd=TRUE)
	dev.off()
}



# ATH all samples
make_Boxplot_All_Thresholds_Labels(ATH_cd_nc_SAS_cor_wo_pollen_0.5_pearson, ATH_cd_nc_SAS_cor_wo_pollen_05_2_pearson, 
	ATH_cd_nc_SAS_cor_wo_pollen_2_5_pearson, ATH_cd_nc_SAS_cor_wo_pollen_5_pearson, samples = "all")

# ATH comparative samples
make_Boxplot_All_Thresholds(ATH_comp_samples_cd_nc_SAS_cor_wo_pollen_0.5_pearson, ATH_comp_samples_cd_nc_SAS_cor_wo_pollen_0.5_2_pearson, 
	ATH_comp_samples_cd_nc_SAS_cor_wo_pollen_2_5_pearson, ATH_comp_samples_cd_nc_SAS_cor_wo_pollen_5_pearson, 
	samples = "comparative")

# AL comparative samples
make_Boxplot_All_Thresholds(AL_comp_samples_cd_nc_SAS_cor_wo_pollen_0.5_pearson, AL_comp_samples_cd_nc_SAS_cor_wo_pollen_0.5_2_pearson, 
	AL_comp_samples_cd_nc_SAS_cor_wo_pollen_2_5_pearson, AL_comp_samples_cd_nc_SAS_cor_wo_pollen_5_pearson)

# CR
make_Boxplot_All_Thresholds(CR_cd_nc_SAS_cor_wo_pollen_0.5_pearson, CR_cd_nc_SAS_cor_wo_pollen_0.5_2_pearson, 
	CR_cd_nc_SAS_cor_wo_pollen_2_5_pearson, CR_cd_nc_SAS_cor_wo_pollen_5_pearson)

# ES
make_Boxplot_All_Thresholds(ES_cd_nc_SAS_cor_wo_pollen_0.5_pearson, ES_cd_nc_SAS_cor_wo_pollen_0.5_2_pearson, 
	ES_cd_nc_SAS_cor_wo_pollen_2_5_pearson, ES_cd_nc_SAS_cor_wo_pollen_5_pearson)

# TH
make_Boxplot_All_Thresholds_Labels(TH_cd_nc_SAS_cor_wo_pollen_0.5_pearson, TH_cd_nc_SAS_cor_wo_pollen_0.5_2_pearson, 
	TH_cd_nc_SAS_cor_wo_pollen_2_5_pearson, TH_cd_nc_SAS_cor_wo_pollen_5_pearson)

# MT
make_Boxplot_All_Thresholds(MT_cd_nc_SAS_cor_wo_pollen_0.5_pearson, MT_cd_nc_SAS_cor_wo_pollen_0.5_2_pearson, 
	MT_cd_nc_SAS_cor_wo_pollen_2_5_pearson, MT_cd_nc_SAS_cor_wo_pollen_5_pearson)

# BD
make_Boxplot_All_Thresholds(BD_cd_nc_SAS_cor_wo_pollen_0.5_pearson, BD_cd_nc_SAS_cor_wo_pollen_0.5_2_pearson, 
	BD_cd_nc_SAS_cor_wo_pollen_2_5_pearson, BD_cd_nc_SAS_cor_wo_pollen_5_pearson)




#------------- ATGE/DevSeq, NAT_length and Spearman - Pearson plots for figure 2 -------------


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



# Prepare data for nc-cd SAS pair overlap length distribution in relation to pearson correlation
plus_strand_overlap = ATH_cd_nc_SAS_cor_wo_pollen_0.5 %>% select(id_plus_strand, start_plus, 
	end_plus, width_query, strand_query, biotype_query, Spearman, Pearson, NAT_overlap_width)
plus_strand_NAT_overlap <- subset(plus_strand_overlap, 
	biotype_query == "lnc_exonic_antisense" | biotype_query == "lnc_intronic_antisense")
names(plus_strand_NAT_overlap) <- c("id", "start", "end", "width", "strand", "biotype", "Spearman", "Pearson", "overlap")

minus_strand_overlap = ATH_cd_nc_SAS_cor_wo_pollen_0.5 %>% select(id_minus_strand, start_minus, 
	end_minus, width_subject, strand_subject, biotype_subject, Spearman, Pearson, NAT_overlap_width)
minus_strand_NAT_overlap <- subset(minus_strand_overlap, 
	biotype_subject == "lnc_exonic_antisense" | biotype_subject == "lnc_intronic_antisense")
names(minus_strand_NAT_overlap) <- c("id", "start", "end", "width", "strand", "biotype", "Spearman", "Pearson", "overlap")

plus_minus_NAT_overlap <- rbind(plus_strand_NAT_overlap, minus_strand_NAT_overlap)

plus_minus_NAT_overlap$percent_overlap <- (plus_minus_NAT_overlap$overlap / plus_minus_NAT_overlap$width) * 100



# Function to prepare data frame and encode data density as color
scatterPlot <- function(x, y) {
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


# Apply scatterPlot function
perc_overlap_pearson <- scatterPlot(plus_minus_NAT_overlap$Pearson, plus_minus_NAT_overlap$percent_overlap)
abs_overlap_pearson <- scatterPlot(plus_minus_NAT_overlap$Pearson, plus_minus_NAT_overlap$overlap)
ATH_all <- scatterPlot(ATH_cd_nc_SAS_cor_wo_pollen_0.5_pearson, ATH_cd_nc_SAS_cor_wo_pollen_0.5_spearman)
ATH_comp <- scatterPlot(ATH_comp_samples_cd_nc_SAS_cor_wo_pollen_0.5_pearson, ATH_comp_samples_cd_nc_SAS_cor_wo_pollen_0.5_spearman)
AL_comp <- scatterPlot(AL_comp_samples_coding_SAS_cor_wo_pollen_pearson, AL_comp_samples_coding_SAS_cor_wo_pollen_spearman)
CR_comp <- scatterPlot(CR_cd_nc_SAS_cor_wo_pollen_0.5_pearson, CR_cd_nc_SAS_cor_wo_pollen_0.5_spearman)
ES_comp <- scatterPlot(ES_cd_nc_SAS_cor_wo_pollen_0.5_pearson, ES_cd_nc_SAS_cor_wo_pollen_0.5_spearman)
TH_comp <- scatterPlot(TH_cd_nc_SAS_cor_wo_pollen_0.5_pearson, TH_cd_nc_SAS_cor_wo_pollen_0.5_spearman)
MT_comp <- scatterPlot(MT_cd_nc_SAS_cor_wo_pollen_0.5_pearson, MT_cd_nc_SAS_cor_wo_pollen_0.5_spearman)
BD_comp <- scatterPlot(BD_cd_nc_SAS_cor_wo_pollen_0.5_pearson, BD_cd_nc_SAS_cor_wo_pollen_0.5_spearman)



# Make scatter plot showing relative overlap (%) in relation to pearson correlation in ATH_all 
NAT_perc_cor <- round((cor(plus_minus_NAT_overlap$Pearson, plus_minus_NAT_overlap$percent_overlap)^2), digits=3)

jpeg(file = file.path(out_dir, "output", "plots", "perc_overlap_pearson.jpg"), 
	width = 4000, height = 4700, res = 825)
rsrt_label = paste("R ^ 2", "==", ".")
p <- ggplot(perc_overlap_pearson, aes(x = x_data, y = y_data)) + 
geom_point(size = 1.25, colour = perc_overlap_pearson$col) + 
scale_x_continuous(limits = c(-1.02,1.02), expand = c(0, 0)) +
scale_y_continuous(limits = c(0,101), expand = c(0, 0)) + 
geom_smooth(method="lm" , color="gray20", fill="#69b3a2", se=TRUE, size=1) +  # use linear regression model
annotate("text", x = -Inf, y = Inf, hjust = -0.38, vjust = 1.6, size=5.35, label = rsrt_label, parse = TRUE) + 
annotate("text", x = -Inf, y = Inf, hjust = -0.38, vjust = 3.8, size=5.35, label = NAT_perc_cor, parse = FALSE) 
p + ggtitle("ATH_all") + theme_bw() + xlab("Pearson") + ylab("NAT overlap (%)") + 
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
dev.off()


# Make scatter plot showing absolute overlap (bp) in relation to pearson correlation in ATH_all 
NAT_abs_cor <- testRsq(plus_minus_NAT_overlap$overlap, plus_minus_NAT_overlap$Pearson)

jpeg(file = file.path(out_dir, "output", "plots", "abs_overlap_pearson.jpg"), 
	width = 4000, height = 4700, res = 825)
corr_label = paste("R ^ 2"," == ", NAT_abs_cor)
p <- ggplot(abs_overlap_pearson, aes(x = x_data, y = y_data)) + 
geom_point(size = 1.25, colour = abs_overlap_pearson$col) + 
scale_x_continuous(limits = c(-1.02,1.02), expand = c(0, 0)) +
scale_y_continuous(limits = c(0,4030), expand = c(0, 0)) + 
geom_smooth(method="lm" , color="gray20", fill="#69b3a2", se=TRUE, size=1) +  # use linear regression model
annotate("text", x = -Inf, y = Inf, hjust = -0.335, vjust = 1.6, size=5.35, label = corr_label, parse = TRUE)
p + ggtitle("ATH_all") + theme_bw() + xlab("Pearson") + ylab("NAT overlap (bp)") + 
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
dev.off()



# Make Spearman-Pearson scatter plots

jpeg(file = file.path(out_dir, "output", "plots", "ATH_cor.jpg"), 
	width = 4000, height = 4700, res = 825)
corr_label = paste("R ^ 2"," == ", ATH_comp_cor)
p <- ggplot(ATH_comp, aes(x = x_data, y = y_data)) +
geom_point(size = 1.25, colour = ATH_comp$col) + 
scale_x_continuous(limits = c(-1.02,1.02), expand = c(0, 0)) +
scale_y_continuous(limits = c(-1.02,1.02), expand = c(0, 0)) + 
geom_smooth(method="lm" , color="gray20", fill="#69b3a2", se=TRUE, size=1) +  # use linear regression model
annotate("text", x = -Inf, y = Inf, hjust = -0.225, vjust = 1.6, size=5.35, label = corr_label, parse = TRUE)
p + ggtitle("ATH") + theme_bw() + xlab("Pearson") + ylab("Spearman") + 
  theme(text=element_text(size=16), 
  	axis.ticks.length = unit(.3, "cm"),
  	plot.margin = unit(c(5.5, 10.5, 42.5, 5.5), "points"),
  	axis.text.x = element_text(colour = "black", size=14.25, angle=0, margin = margin(t = 8.25, r = 0, b = 0, l = 0)), 
  	axis.text.y = element_text(colour = "black", size=14.25, angle=0, margin = margin(t = 0, r = 8.25, b = 0, l = 0)),
  	axis.title.x = element_text(colour = "black", margin = margin(t = 14.5, r = 0, b = 1, l = 0)),
  	axis.title.y = element_text(colour = "black", margin = margin(t = 0, r = 9, b = 0, l = 1)),
  	plot.title = element_text(colour = "black", size=17, margin = margin(t = 9, r = 0, b = 18, l = 0), hjust = 0.5),
  	legend.position = "bottom",
  	panel.border = element_rect(colour = "black", fill=NA, size=1.2))
dev.off()

jpeg(file = file.path(out_dir, "output", "plots", "AL_cor.jpg"), 
	width = 4000, height = 4700, res = 825)
corr_label = paste("R ^ 2"," == ", AL_comp_cor)
p <- ggplot(AL_comp, aes(x = x_data, y = y_data)) +
geom_point(size = 1.25, colour = AL_comp$col) + 
scale_x_continuous(limits = c(-1.02,1.02), expand = c(0, 0)) +
scale_y_continuous(limits = c(-1.02,1.02), expand = c(0, 0)) + 
geom_smooth(method="lm" , color="gray20", fill="#69b3a2", se=TRUE, size=1) +  # use linear regression model
annotate("text", x = -Inf, y = Inf, hjust = -0.225, vjust = 1.6, size=5.35, label = corr_label, parse = TRUE)
p + ggtitle("AL") + theme_bw() + xlab("Pearson") + ylab("Spearman") + 
  theme(text=element_text(size=16), 
  	axis.ticks.length = unit(.3, "cm"),
  	plot.margin = unit(c(5.5, 10.5, 42.5, 5.5), "points"),
  	axis.text.x = element_text(colour = "black", size=14.25, angle=0, margin = margin(t = 8.25, r = 0, b = 0, l = 0)), 
  	axis.text.y = element_text(colour = "black", size=14.25, angle=0, margin = margin(t = 0, r = 8.25, b = 0, l = 0)),
  	axis.title.x = element_text(colour = "black", margin = margin(t = 14.5, r = 0, b = 1, l = 0)),
  	axis.title.y = element_text(colour = "black", margin = margin(t = 0, r = 9, b = 0, l = 1)),
  	plot.title = element_text(colour = "black", size=17, margin = margin(t = 9, r = 0, b = 18, l = 0), hjust = 0.5),
  	legend.position = "bottom",
  	panel.border = element_rect(colour = "black", fill=NA, size=1.2))
dev.off()

jpeg(file = file.path(out_dir, "output", "plots", "CR_cor.jpg"), 
	width = 4000, height = 4700, res = 825)
corr_label = paste("R ^ 2"," == ", CR_cor)
p <- ggplot(CR_comp, aes(x = x_data, y = y_data)) +
geom_point(size = 1.25, colour = CR_comp$col) + 
scale_x_continuous(limits = c(-1.02,1.02), expand = c(0, 0)) +
scale_y_continuous(limits = c(-1.02,1.02), expand = c(0, 0)) + 
geom_smooth(method="lm" , color="gray20", fill="#69b3a2", se=TRUE, size=1) +  # use linear regression model
annotate("text", x = -Inf, y = Inf, hjust = -0.225, vjust = 1.6, size=5.35, label = corr_label, parse = TRUE)
p + ggtitle("CR") + theme_bw() + xlab("Pearson") + ylab("Spearman") + 
  theme(text=element_text(size=16), 
  	axis.ticks.length = unit(.3, "cm"),
  	plot.margin = unit(c(5.5, 10.5, 42.5, 5.5), "points"),
  	axis.text.x = element_text(colour = "black", size=14.25, angle=0, margin = margin(t = 8.25, r = 0, b = 0, l = 0)), 
  	axis.text.y = element_text(colour = "black", size=14.25, angle=0, margin = margin(t = 0, r = 8.25, b = 0, l = 0)),
  	axis.title.x = element_text(colour = "black", margin = margin(t = 14.5, r = 0, b = 1, l = 0)),
  	axis.title.y = element_text(colour = "black", margin = margin(t = 0, r = 9, b = 0, l = 1)),
  	plot.title = element_text(colour = "black", size=17, margin = margin(t = 9, r = 0, b = 18, l = 0), hjust = 0.5),
  	legend.position = "bottom",
  	panel.border = element_rect(colour = "black", fill=NA, size=1.2))
dev.off()

jpeg(file = file.path(out_dir, "output", "plots", "ES_cor.jpg"), 
	width = 4000, height = 4700, res = 825)
corr_label = paste("R ^ 2"," == ", ES_cor)
p <- ggplot(ES_comp, aes(x = x_data, y = y_data)) +
geom_point(size = 1.25, colour = ES_comp$col) + 
scale_x_continuous(limits = c(-1.02,1.02), expand = c(0, 0)) +
scale_y_continuous(limits = c(-1.02,1.02), expand = c(0, 0)) + 
geom_smooth(method="lm" , color="gray20", fill="#69b3a2", se=TRUE, size=1) +  # use linear regression model
annotate("text", x = -Inf, y = Inf, hjust = -0.225, vjust = 1.6, size=5.35, label = corr_label, parse = TRUE)
p + ggtitle("ES") + theme_bw() + xlab("Pearson") + ylab("Spearman") + 
  theme(text=element_text(size=16), 
  	axis.ticks.length = unit(.3, "cm"),
  	plot.margin = unit(c(5.5, 10.5, 42.5, 5.5), "points"),
  	axis.text.x = element_text(colour = "black", size=14.25, angle=0, margin = margin(t = 8.25, r = 0, b = 0, l = 0)), 
  	axis.text.y = element_text(colour = "black", size=14.25, angle=0, margin = margin(t = 0, r = 8.25, b = 0, l = 0)),
  	axis.title.x = element_text(colour = "black", margin = margin(t = 14.5, r = 0, b = 1, l = 0)),
  	axis.title.y = element_text(colour = "black", margin = margin(t = 0, r = 9, b = 0, l = 1)),
  	plot.title = element_text(colour = "black", size=17, margin = margin(t = 9, r = 0, b = 18, l = 0), hjust = 0.5),
  	legend.position = "bottom",
  	panel.border = element_rect(colour = "black", fill=NA, size=1.2))
dev.off()

jpeg(file = file.path(out_dir, "output", "plots", "TH_cor.jpg"), 
	width = 4000, height = 4700, res = 825)
corr_label = paste("R ^ 2"," == ", TH_cor)
p <- ggplot(TH_comp, aes(x = x_data, y = y_data)) +
geom_point(size = 1.25, colour = TH_comp$col) + 
scale_x_continuous(limits = c(-1.02,1.02), expand = c(0, 0)) +
scale_y_continuous(limits = c(-1.02,1.02), expand = c(0, 0)) + 
geom_smooth(method="lm" , color="gray20", fill="#69b3a2", se=TRUE, size=1) +  # use linear regression model
annotate("text", x = -Inf, y = Inf, hjust = -0.225, vjust = 1.6, size=5.35, label = corr_label, parse = TRUE)
p + ggtitle("TH") + theme_bw() + xlab("Pearson") + ylab("Spearman") + 
  theme(text=element_text(size=16), 
  	axis.ticks.length = unit(.3, "cm"),
  	plot.margin = unit(c(5.5, 10.5, 42.5, 5.5), "points"),
  	axis.text.x = element_text(colour = "black", size=14.25, angle=0, margin = margin(t = 8.25, r = 0, b = 0, l = 0)), 
  	axis.text.y = element_text(colour = "black", size=14.25, angle=0, margin = margin(t = 0, r = 8.25, b = 0, l = 0)),
  	axis.title.x = element_text(colour = "black", margin = margin(t = 14.5, r = 0, b = 1, l = 0)),
  	axis.title.y = element_text(colour = "black", margin = margin(t = 0, r = 9, b = 0, l = 1)),
  	plot.title = element_text(colour = "black", size=17, margin = margin(t = 9, r = 0, b = 18, l = 0), hjust = 0.5),
  	legend.position = "bottom",
  	panel.border = element_rect(colour = "black", fill=NA, size=1.2))
dev.off()

jpeg(file = file.path(out_dir, "output", "plots", "MT_cor.jpg"), 
	width = 4000, height = 4700, res = 825)
corr_label = paste("R ^ 2"," == ", MT_cor)
p <- ggplot(MT_comp, aes(x = x_data, y = y_data)) +
geom_point(size = 1.25, colour = MT_comp$col) + 
scale_x_continuous(limits = c(-1.02,1.02), expand = c(0, 0)) +
scale_y_continuous(limits = c(-1.02,1.02), expand = c(0, 0)) + 
geom_smooth(method="lm" , color="gray20", fill="#69b3a2", se=TRUE, size=1) +  # use linear regression model
annotate("text", x = -Inf, y = Inf, hjust = -0.225, vjust = 1.6, size=5.35, label = corr_label, parse = TRUE)
p + ggtitle("MT") + theme_bw() + xlab("Pearson") + ylab("Spearman") + 
  theme(text=element_text(size=16), 
  	axis.ticks.length = unit(.3, "cm"),
  	plot.margin = unit(c(5.5, 10.5, 42.5, 5.5), "points"),
  	axis.text.x = element_text(colour = "black", size=14.25, angle=0, margin = margin(t = 8.25, r = 0, b = 0, l = 0)), 
  	axis.text.y = element_text(colour = "black", size=14.25, angle=0, margin = margin(t = 0, r = 8.25, b = 0, l = 0)),
  	axis.title.x = element_text(colour = "black", margin = margin(t = 14.5, r = 0, b = 1, l = 0)),
  	axis.title.y = element_text(colour = "black", margin = margin(t = 0, r = 9, b = 0, l = 1)),
  	plot.title = element_text(colour = "black", size=17, margin = margin(t = 9, r = 0, b = 18, l = 0), hjust = 0.5),
  	legend.position = "bottom",
  	panel.border = element_rect(colour = "black", fill=NA, size=1.2))
dev.off()

jpeg(file = file.path(out_dir, "output", "plots", "BD_cor.jpg"), 
	width = 4000, height = 4700, res = 825)
corr_label = paste("R ^ 2"," == ", BD_cor)
p <- ggplot(BD_comp, aes(x = x_data, y = y_data)) +
geom_point(size = 1.25, colour = BD_comp$col) + 
scale_x_continuous(limits = c(-1.02,1.02), expand = c(0, 0)) +
scale_y_continuous(limits = c(-1.02,1.02), expand = c(0, 0)) + 
geom_smooth(method="lm" , color="gray20", fill="#69b3a2", se=TRUE, size=1) +  # use linear regression model
annotate("text", x = -Inf, y = Inf, hjust = -0.225, vjust = 1.6, size=5.35, label = corr_label, parse = TRUE)
p + ggtitle("BD") + theme_bw() + xlab("Pearson") + ylab("Spearman") + 
  theme(text=element_text(size=16), 
  	axis.ticks.length = unit(.3, "cm"),
  	plot.margin = unit(c(5.5, 10.5, 42.5, 5.5), "points"),
  	axis.text.x = element_text(colour = "black", size=14.25, angle=0, margin = margin(t = 8.25, r = 0, b = 0, l = 0)), 
  	axis.text.y = element_text(colour = "black", size=14.25, angle=0, margin = margin(t = 0, r = 8.25, b = 0, l = 0)),
  	axis.title.x = element_text(colour = "black", margin = margin(t = 14.5, r = 0, b = 1, l = 0)),
  	axis.title.y = element_text(colour = "black", margin = margin(t = 0, r = 9, b = 0, l = 1)),
  	plot.title = element_text(colour = "black", size=17, margin = margin(t = 9, r = 0, b = 18, l = 0), hjust = 0.5),
  	legend.position = "bottom",
  	panel.border = element_rect(colour = "black", fill=NA, size=1.2))
dev.off()








































































#------------------------ Older density plots used for figure V1 --------------------------


# Prepare data for ggplot2 ATH density plot cd-cd, nc-cd_spearman, nc-cd_pearson V2
combine_Species_Data <- function(coding, pearson, spearman) {

	number_values <- length(coding)+length(pearson)+length(spearman)
	species_name = as.data.frame(rep(c(sub("\\_.*", "", deparse(substitute(pearson)))),each=number_values))
	names(species_name) <- "species"

	class_0 = as.data.frame(rep(c("coding"),each=length(coding)))
	names(class_0) <- "class"
	class_1 = as.data.frame(rep(c("pearson"),each=length(pearson)))
	names(class_1) <- "class"
	class_2 = as.data.frame(rep(c("spearman"),each=length(spearman)))
	names(class_2) <- "class"

	cor_values_0 = as.data.frame(coding)
	names(cor_values_0) <- "correlation"
	cor_values_1 = as.data.frame(pearson)
	names(cor_values_1) <- "correlation"
	cor_values_2 = as.data.frame(spearman)
	names(cor_values_2) <- "correlation"

	species_df = data.frame(species_name, rbind(class_0, class_1, class_2) , 
		rbind(cor_values_0, cor_values_1, cor_values_2))
	species_df <- na.omit(species_df)

	return(species_df)
}

ATH_all_nccd_spearman_pearson <- combine_Species_Data(ATH_coding_SAS_cor_wo_pollen_pearson, 
	ATH_cd_nc_SAS_cor_wo_pollen_0.5_pearson, ATH_cd_nc_SAS_cor_wo_pollen_0.5_spearman)


# Make cd-cd, nc-cd_spearman, nc-cd_pearson density plot

jpeg(file = file.path(out_dir, "output", "plots", "ATH_all_nccd_spearman_pearson_V2.jpg"), 
	width = 4000, height = 4700, res = 825)
p <- ggplot(ATH_all_nccd_spearman_pearson, aes(x=correlation, group=class, fill=class, colour=class, linetype=class)) +
geom_density(adjust=1.5, alpha=0.35, size=1.5) + 
scale_x_continuous(limits = c(-1,1), expand = c(0, 0)) +
scale_y_continuous(limits = c(0,1.28), expand = c(0, 0))
p + ggtitle("ATH_all") + theme_bw() + scale_fill_manual(values = c("white", "#52b540", "#00468b")) +
  scale_color_manual(values = c("gray65", "#52b540", "#00468b")) + 
  scale_linetype_manual(values = c("dotted","solid","solid")) + 
  theme(text=element_text(size=16), 
  	axis.ticks.length = unit(.3, "cm"),
  	plot.margin = unit(c(5.5, 10.5, 3.0, 3.5), "points"),
  	axis.text.x = element_text(colour = "black", size=14.25, angle=0, margin = margin(t = 8.25, r = 0, b = 0, l = 0)), 
  	axis.text.y = element_text(colour = "black", size=14.25, angle=0, margin = margin(t = 0, r = 8.25, b = 0, l = 0)),
  	axis.title.x = element_text(colour = "black", margin = margin(t = 15.75, r = 0, b = 0, l = 0)),
  	axis.title.y = element_text(colour = "black", margin = margin(t = 0, r = 15.25, b = 0, l = 0)),
  	plot.title = element_text(colour = "black", size=17, margin = margin(t = 10.25, r = 0, b = 16.75, l = 0), hjust = 0.5),
  	legend.position = "bottom",
  	panel.border = element_rect(colour = "black", fill=NA, size=1.2))
dev.off()



# Prepare data for ggplot2 boxplot
combine_Species_Data <- function(coding, pearson, spearman) {

	number_values <- length(coding)+length(pearson)+length(spearman)
	species_name = as.data.frame(rep(c(sub("\\_.*", "", deparse(substitute(pearson)))),each=number_values))
	names(species_name) <- "species"

	class_0 = as.data.frame(rep(c("coding"),each=length(coding)))
	names(class_0) <- "class"
	class_1 = as.data.frame(rep(c("pearson"),each=length(pearson)))
	names(class_1) <- "class"
	class_2 = as.data.frame(rep(c("spearman"),each=length(spearman)))
	names(class_2) <- "class"

	cor_values_0 = as.data.frame(coding)
	names(cor_values_0) <- "correlation"
	cor_values_1 = as.data.frame(pearson)
	names(cor_values_1) <- "correlation"
	cor_values_2 = as.data.frame(spearman)
	names(cor_values_2) <- "correlation"

	species_df = data.frame(species_name, rbind(class_0, class_1, class_2) , 
		rbind(cor_values_0, cor_values_1, cor_values_2))
	species_df <- na.omit(species_df)

	return(species_df)
}

ATH_all_nccd_spearman_pearson <- combine_Species_Data(ATH_coding_SAS_cor_wo_pollen_pearson, 
	ATH_cd_nc_SAS_cor_wo_pollen_0.5_pearson, ATH_cd_nc_SAS_cor_wo_pollen_0.5_spearman)
ATH_comp_nccd_spearman_pearson <- combine_Species_Data(ATH_comp_samples_coding_SAS_cor_wo_pollen_pearson, 
	ATH_comp_samples_cd_nc_SAS_cor_wo_pollen_0.5_pearson, ATH_comp_samples_cd_nc_SAS_cor_wo_pollen_0.5_spearman)
AL_comp_nccd_spearman_pearson <- combine_Species_Data(AL_comp_samples_coding_SAS_cor_wo_pollen_pearson, 
	AL_comp_samples_cd_nc_SAS_cor_wo_pollen_0.5_pearson, AL_comp_samples_cd_nc_SAS_cor_wo_pollen_0.5_spearman)
CR_nccd_spearman_pearson <- combine_Species_Data(CR_coding_SAS_cor_wo_pollen_pearson, 
	CR_cd_nc_SAS_cor_wo_pollen_0.5_pearson, CR_cd_nc_SAS_cor_wo_pollen_0.5_spearman)
ES_nccd_spearman_pearson <- combine_Species_Data(ES_coding_SAS_cor_wo_pollen_pearson, 
	ES_cd_nc_SAS_cor_wo_pollen_0.5_pearson, ES_cd_nc_SAS_cor_wo_pollen_0.5_spearman)
TH_nccd_spearman_pearson <- combine_Species_Data(TH_coding_SAS_cor_wo_pollen_pearson, 
	TH_cd_nc_SAS_cor_wo_pollen_0.5_pearson, TH_cd_nc_SAS_cor_wo_pollen_0.5_spearman)
MT_nccd_spearman_pearson <- combine_Species_Data(MT_coding_SAS_cor_wo_pollen_pearson, 
	MT_cd_nc_SAS_cor_wo_pollen_0.5_pearson, MT_cd_nc_SAS_cor_wo_pollen_0.5_spearman)
BD_nccd_spearman_pearson <- combine_Species_Data(BD_coding_SAS_cor_wo_pollen_pearson, 
	BD_cd_nc_SAS_cor_wo_pollen_0.5_pearson, BD_cd_nc_SAS_cor_wo_pollen_0.5_spearman)




jpeg(file = file.path(out_dir, "output", "plots", "ATH_all_nccd_spearman_pearson_V2.jpg"), 
	width = 4000, height = 4700, res = 825)
p <- ggplot(ATH_all_nccd_spearman_pearson, aes(x=correlation, group=class, fill=class, colour=class, linetype=class)) +
geom_density(adjust=1.5, alpha=0.35, size=1.5) + 
scale_x_continuous(limits = c(-1,1), expand = c(0, 0)) +
scale_y_continuous(limits = c(0,1.28), expand = c(0, 0))
p + ggtitle("ATH_all") + theme_bw() + scale_fill_manual(values = c("white", "#52b540", "#00468b")) +
  scale_color_manual(values = c("gray65", "#52b540", "#00468b")) + 
  scale_linetype_manual(values = c("dotted","solid","solid")) + 
  theme(text=element_text(size=16), 
  	axis.ticks.length = unit(.3, "cm"),
  	plot.margin = unit(c(5.5, 10.5, 3.0, 3.5), "points"),
  	axis.text.x = element_text(colour = "black", size=14.25, angle=0, margin = margin(t = 8.25, r = 0, b = 0, l = 0)), 
  	axis.text.y = element_text(colour = "black", size=14.25, angle=0, margin = margin(t = 0, r = 8.25, b = 0, l = 0)),
  	axis.title.x = element_text(colour = "black", margin = margin(t = 15.75, r = 0, b = 0, l = 0)),
  	axis.title.y = element_text(colour = "black", margin = margin(t = 0, r = 15.25, b = 0, l = 0)),
  	plot.title = element_text(colour = "black", size=17, margin = margin(t = 10.25, r = 0, b = 16.75, l = 0), hjust = 0.5),
  	legend.position = "bottom",
  	panel.border = element_rect(colour = "black", fill=NA, size=1.2))
dev.off()



jpeg(file = file.path(out_dir, "output", "plots", "ATH_all_nccd_spearman_pearson.jpg"), 
	width = 4000, height = 4700, res = 825)
corr_label = paste("R = ", ATH_all_cor)
p <- ggplot(ATH_all_nccd_spearman_pearson, aes(x=correlation, group=class, fill=class, colour=class, linetype=class)) +
geom_density(adjust=1.5, alpha=0.35, size=1.5) + 
scale_x_continuous(limits = c(-1,1), expand = c(0, 0)) +
scale_y_continuous(limits = c(0,1.28), expand = c(0, 0)) + 
annotate("text", x = -Inf, y = Inf, hjust = -0.225, vjust = 2.1, size=5.35, label = corr_label)
p + ggtitle("ATH_all") + theme_bw() + scale_fill_manual(values = c("white", "#52b540", "#00468b")) +
  scale_color_manual(values = c("gray65", "#52b540", "#00468b")) + 
  scale_linetype_manual(values = c("dotted","solid","solid")) + 
  theme(text=element_text(size=16), 
  	axis.ticks.length = unit(.3, "cm"),
  	plot.margin = unit(c(5.5, 10.5, 5.5, 5.5), "points"),
  	axis.text.x = element_text(colour = "black", size=14.25, angle=0, margin = margin(t = 8.25, r = 0, b = 0, l = 0)), 
  	axis.text.y = element_text(colour = "black", size=14.25, angle=0, margin = margin(t = 0, r = 8.25, b = 0, l = 0)),
  	axis.title.x = element_text(colour = "black", margin = margin(t = 13.25, r = 0, b = 0, l = 0)),
  	axis.title.y = element_text(colour = "black", margin = margin(t = 0, r = 13.25, b = 0, l = 0)),
  	plot.title = element_text(colour = "black", size=17, margin = margin(t = 8, r = 0, b = 19, l = 0), hjust = 0.5),
  	legend.position = "bottom",
  	panel.border = element_rect(colour = "black", fill=NA, size=1.2))
dev.off()

jpeg(file = file.path(out_dir, "output", "plots", "ATH_comp_nccd_spearman_pearson.jpg"), 
	width = 4000, height = 4700, res = 825)
corr_label = paste("R = ", ATH_comp_cor)
p <- ggplot(ATH_comp_nccd_spearman_pearson, aes(x=correlation, group=class, fill=class, colour=class, linetype=class)) +
geom_density(adjust=1.5, alpha=0.35, size=1.5) + 
scale_x_continuous(limits = c(-1,1), expand = c(0, 0)) +
scale_y_continuous(limits = c(0,1.1315), expand = c(0, 0)) + 
annotate("text", x = -Inf, y = Inf, hjust = -0.225, vjust = 2.1, size=5.35, label = corr_label)
p + ggtitle("ATH_comp") + theme_bw() + scale_fill_manual(values = c("white", "#52b540", "#00468b")) +
  scale_color_manual(values = c("gray65", "#52b540", "#00468b")) + 
  scale_linetype_manual(values = c("dotted","solid","solid")) + 
  theme(text=element_text(size=16), 
  	axis.ticks.length = unit(.3, "cm"),
  	plot.margin = unit(c(5.5, 10.5, 5.5, 5.5), "points"),
  	axis.text.x = element_text(colour = "black", size=14.25, angle=0, margin = margin(t = 8.25, r = 0, b = 0, l = 0)), 
  	axis.text.y = element_text(colour = "black", size=14.25, angle=0, margin = margin(t = 0, r = 8.25, b = 0, l = 0)),
  	axis.title.x = element_text(colour = "black", margin = margin(t = 13.25, r = 0, b = 0, l = 0)),
  	axis.title.y = element_text(colour = "black", margin = margin(t = 0, r = 13.25, b = 0, l = 0)),
  	plot.title = element_text(colour = "black", size=17, margin = margin(t = 8, r = 0, b = 19, l = 0), hjust = 0.5),
  	legend.position = "bottom",
  	panel.border = element_rect(colour = "black", fill=NA, size=1.2))
dev.off()

jpeg(file = file.path(out_dir, "output", "plots", "AL_comp_nccd_spearman_pearson.jpg"), 
	width = 4000, height = 4700, res = 825)
corr_label = paste("R = ", AL_comp_cor)
p <- ggplot(AL_comp_nccd_spearman_pearson, aes(x=correlation, group=class, fill=class, colour=class, linetype=class)) +
geom_density(adjust=1.5, alpha=0.35, size=1.5) + 
scale_x_continuous(limits = c(-1,1), expand = c(0, 0)) +
scale_y_continuous(limits = c(0,1.002), expand = c(0, 0)) + 
annotate("text", x = -Inf, y = Inf, hjust = -0.225, vjust = 2.1, size=5.35, label = corr_label)
p + ggtitle("AL") + theme_bw() + scale_fill_manual(values = c("white", "#52b540", "#00468b")) +
  scale_color_manual(values = c("gray65", "#52b540", "#00468b")) + 
  scale_linetype_manual(values = c("dotted","solid","solid")) + 
  theme(text=element_text(size=16), 
  	axis.ticks.length = unit(.3, "cm"),
  	plot.margin = unit(c(5.5, 10.5, 5.5, 5.5), "points"),
  	axis.text.x = element_text(colour = "black", size=14.25, angle=0, margin = margin(t = 8.25, r = 0, b = 0, l = 0)), 
  	axis.text.y = element_text(colour = "black", size=14.25, angle=0, margin = margin(t = 0, r = 8.25, b = 0, l = 0)),
  	axis.title.x = element_text(colour = "black", margin = margin(t = 13.25, r = 0, b = 0, l = 0)),
  	axis.title.y = element_text(colour = "black", margin = margin(t = 0, r = 13.25, b = 0, l = 0)),
  	plot.title = element_text(colour = "black", size=17, margin = margin(t = 8, r = 0, b = 19, l = 0), hjust = 0.5),
  	legend.position = "bottom",
  	panel.border = element_rect(colour = "black", fill=NA, size=1.2))
dev.off()

jpeg(file = file.path(out_dir, "output", "plots", "CR_nccd_spearman_pearson.jpg"), 
	width = 4000, height = 4700, res = 825)
corr_label = paste("R = ", CR_cor)
p <- ggplot(CR_nccd_spearman_pearson, aes(x=correlation, group=class, fill=class, colour=class, linetype=class)) +
geom_density(adjust=1.5, alpha=0.35, size=1.5) + 
scale_x_continuous(limits = c(-1,1), expand = c(0, 0)) +
scale_y_continuous(limits = c(0,1.263), expand = c(0, 0)) + 
annotate("text", x = -Inf, y = Inf, hjust = -0.225, vjust = 2.1, size=5.35, label = corr_label)
p + ggtitle("CR") + theme_bw() + scale_fill_manual(values = c("white", "#52b540", "#00468b")) +
  scale_color_manual(values = c("gray65", "#52b540", "#00468b")) + 
  scale_linetype_manual(values = c("dotted","solid","solid")) + 
  theme(text=element_text(size=16), 
  	axis.ticks.length = unit(.3, "cm"),
  	plot.margin = unit(c(5.5, 10.5, 5.5, 5.5), "points"),
  	axis.text.x = element_text(colour = "black", size=14.25, angle=0, margin = margin(t = 8.25, r = 0, b = 0, l = 0)), 
  	axis.text.y = element_text(colour = "black", size=14.25, angle=0, margin = margin(t = 0, r = 8.25, b = 0, l = 0)),
  	axis.title.x = element_text(colour = "black", margin = margin(t = 13.25, r = 0, b = 0, l = 0)),
  	axis.title.y = element_text(colour = "black", margin = margin(t = 0, r = 13.25, b = 0, l = 0)),
  	plot.title = element_text(colour = "black", size=17, margin = margin(t = 8, r = 0, b = 19, l = 0), hjust = 0.5),
  	legend.position = "bottom",
  	panel.border = element_rect(colour = "black", fill=NA, size=1.2))
dev.off()

jpeg(file = file.path(out_dir, "output", "plots", "ES_nccd_spearman_pearson.jpg"), 
	width = 4000, height = 4700, res = 825)
corr_label = paste("R = ", ES_cor)
p <- ggplot(ES_nccd_spearman_pearson, aes(x=correlation, group=class, fill=class, colour=class, linetype=class)) +
geom_density(adjust=1.5, alpha=0.35, size=1.5) + 
scale_x_continuous(limits = c(-1,1), expand = c(0, 0)) +
scale_y_continuous(limits = c(0,1.158), expand = c(0, 0)) + 
annotate("text", x = -Inf, y = Inf, hjust = -0.225, vjust = 2.1, size=5.35, label = corr_label)
p + ggtitle("ES") + theme_bw() + scale_fill_manual(values = c("white", "#52b540", "#00468b")) +
  scale_color_manual(values = c("gray65", "#52b540", "#00468b")) + 
  scale_linetype_manual(values = c("dotted","solid","solid")) + 
  theme(text=element_text(size=16), 
  	axis.ticks.length = unit(.3, "cm"),
  	plot.margin = unit(c(5.5, 10.5, 5.5, 5.5), "points"),
  	axis.text.x = element_text(colour = "black", size=14.25, angle=0, margin = margin(t = 8.25, r = 0, b = 0, l = 0)), 
  	axis.text.y = element_text(colour = "black", size=14.25, angle=0, margin = margin(t = 0, r = 8.25, b = 0, l = 0)),
  	axis.title.x = element_text(colour = "black", margin = margin(t = 13.25, r = 0, b = 0, l = 0)),
  	axis.title.y = element_text(colour = "black", margin = margin(t = 0, r = 13.25, b = 0, l = 0)),
  	plot.title = element_text(colour = "black", size=17, margin = margin(t = 8, r = 0, b = 19, l = 0), hjust = 0.5),
  	legend.position = "bottom",
  	panel.border = element_rect(colour = "black", fill=NA, size=1.2))
dev.off()

jpeg(file = file.path(out_dir, "output", "plots", "TH_nccd_spearman_pearson.jpg"), 
	width = 4000, height = 4700, res = 825)
corr_label = paste("R = ", TH_cor)
p <- ggplot(TH_nccd_spearman_pearson, aes(x=correlation, group=class, fill=class, colour=class, linetype=class)) +
geom_density(adjust=1.5, alpha=0.35, size=1.5) + 
scale_x_continuous(limits = c(-1,1), expand = c(0, 0)) +
scale_y_continuous(limits = c(0,0.8565), expand = c(0, 0)) + 
annotate("text", x = -Inf, y = Inf, hjust = -0.225, vjust = 2.1, size=5.35, label = corr_label)
p + ggtitle("TH") + theme_bw() + scale_fill_manual(values = c("white", "#52b540", "#00468b")) +
  scale_color_manual(values = c("gray65", "#52b540", "#00468b")) + 
  scale_linetype_manual(values = c("dotted","solid","solid")) + 
  theme(text=element_text(size=16), 
  	axis.ticks.length = unit(.3, "cm"),
  	plot.margin = unit(c(5.5, 10.5, 5.5, 5.5), "points"),
  	axis.text.x = element_text(colour = "black", size=14.25, angle=0, margin = margin(t = 8.25, r = 0, b = 0, l = 0)), 
  	axis.text.y = element_text(colour = "black", size=14.25, angle=0, margin = margin(t = 0, r = 8.25, b = 0, l = 0)),
  	axis.title.x = element_text(colour = "black", margin = margin(t = 13.25, r = 0, b = 0, l = 0)),
  	axis.title.y = element_text(colour = "black", margin = margin(t = 0, r = 13.25, b = 0, l = 0)),
  	plot.title = element_text(colour = "black", size=17, margin = margin(t = 8, r = 0, b = 19, l = 0), hjust = 0.5),
  	legend.position = "bottom",
  	panel.border = element_rect(colour = "black", fill=NA, size=1.2))
dev.off()

jpeg(file = file.path(out_dir, "output", "plots", "MT_nccd_spearman_pearson.jpg"), 
	width = 4000, height = 4700, res = 825)
corr_label = paste("R = ", MT_cor)
p <- ggplot(MT_nccd_spearman_pearson, aes(x=correlation, group=class, fill=class, colour=class, linetype=class)) +
geom_density(adjust=1.5, alpha=0.35, size=1.5) + 
scale_x_continuous(limits = c(-1,1), expand = c(0, 0)) +
scale_y_continuous(limits = c(0,1.2158), expand = c(0, 0)) + 
annotate("text", x = -Inf, y = Inf, hjust = -0.225, vjust = 2.1, size=5.35, label = corr_label)
p + ggtitle("MT") + theme_bw() + scale_fill_manual(values = c("white", "#52b540", "#00468b")) +
  scale_color_manual(values = c("gray65", "#52b540", "#00468b")) + 
  scale_linetype_manual(values = c("dotted","solid","solid")) + 
  theme(text=element_text(size=16), 
  	axis.ticks.length = unit(.3, "cm"),
  	plot.margin = unit(c(5.5, 10.5, 5.5, 5.5), "points"),
  	axis.text.x = element_text(colour = "black", size=14.25, angle=0, margin = margin(t = 8.25, r = 0, b = 0, l = 0)), 
  	axis.text.y = element_text(colour = "black", size=14.25, angle=0, margin = margin(t = 0, r = 8.25, b = 0, l = 0)),
  	axis.title.x = element_text(colour = "black", margin = margin(t = 13.25, r = 0, b = 0, l = 0)),
  	axis.title.y = element_text(colour = "black", margin = margin(t = 0, r = 13.25, b = 0, l = 0)),
  	plot.title = element_text(colour = "black", size=17, margin = margin(t = 8, r = 0, b = 19, l = 0), hjust = 0.5),
  	legend.position = "bottom",
  	panel.border = element_rect(colour = "black", fill=NA, size=1.2))
dev.off()

jpeg(file = file.path(out_dir, "output", "plots", "BD_nccd_spearman_pearson.jpg"), 
	width = 4000, height = 4700, res = 825)
corr_label = paste("R = ", BD_cor)
p <- ggplot(BD_nccd_spearman_pearson, aes(x=correlation, group=class, fill=class, colour=class, linetype=class)) +
geom_density(adjust=1.5, alpha=0.35, size=1.5) + 
scale_x_continuous(limits = c(-1,1), expand = c(0, 0)) +
scale_y_continuous(limits = c(0,1.078), expand = c(0, 0)) + 
annotate("text", x = -Inf, y = Inf, hjust = -0.225, vjust = 2.1, size=5.35, label = corr_label)
p + ggtitle("BD") + theme_bw() + scale_fill_manual(values = c("white", "#52b540", "#00468b")) +
  scale_color_manual(values = c("gray65", "#52b540", "#00468b")) + 
  scale_linetype_manual(values = c("dotted","solid","solid")) + 
  theme(text=element_text(size=16), 
  	axis.ticks.length = unit(.3, "cm"),
  	plot.margin = unit(c(5.5, 10.5, 5.5, 5.5), "points"),
  	axis.text.x = element_text(colour = "black", size=14.25, angle=0, margin = margin(t = 8.25, r = 0, b = 0, l = 0)), 
  	axis.text.y = element_text(colour = "black", size=14.25, angle=0, margin = margin(t = 0, r = 8.25, b = 0, l = 0)),
  	axis.title.x = element_text(colour = "black", margin = margin(t = 13.25, r = 0, b = 0, l = 0)),
  	axis.title.y = element_text(colour = "black", margin = margin(t = 0, r = 13.25, b = 0, l = 0)),
  	plot.title = element_text(colour = "black", size=17, margin = margin(t = 8, r = 0, b = 19, l = 0), hjust = 0.5),
  	legend.position = "bottom",
  	panel.border = element_rect(colour = "black", fill=NA, size=1.2))
dev.off()









