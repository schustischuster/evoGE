
# This script loads gene expression correlation tables of sense-antisense (SAS) gene pairs 
# both coding/cisNATs SAS and coding-coding SAS

# coding_genes_tables of coding SAS have the following format:
# id_plus_strand / start_plus / end_plus / id_minus_strand / start_minus / end_minus / biotype / 
#   Spearman / Pearson


# Set file path and input files
in_dir_cd <- "./output/overlapp_cd_genes"
in_dir_nc <- "./output/overlapp_nc_genes"


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
	as.numeric(unlist(x[, 11]))})
names(NAT_genes_tables_wo_pollen_list_spearman) <- paste(names(
	NAT_genes_tables_wo_pollen_list_spearman),"_spearman", sep="")

NAT_genes_tables_wo_pollen_list_pearson <- lapply(NAT_genes_tables_wo_pollen_list, function(x) {
	as.numeric(unlist(x[, 12]))})
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




#--------------------------------------- Generate plots ----------------------------------------


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
	col = c("darkseagreen", "gold2", "darkseagreen", "gold2", "darkseagreen", "gold2", "darkseagreen", 
		"gold2", "darkseagreen", "gold2", "darkseagreen", "gold2", "darkseagreen", "gold2"), 
	boxwex = 0.85, 
	lwd = 1.35, 
	whisklty = 1, 
	at = c(1,2,3.5,4.5), 
	notch = FALSE
	)
	title("SAS pairs in ATH", adj = 0.50, line = 1.3)
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
	bty='n', horiz = TRUE, fill = c("darkseagreen", "gold2"), cex = 1.1, x.intersp = 0.5)
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
	col = c("darkseagreen", "gold2", "darkseagreen", "gold2", "darkseagreen", "gold2", "darkseagreen", 
		"gold2", "darkseagreen", "gold2", "darkseagreen", "gold2", "darkseagreen", "gold2"), 
	boxwex = 0.85, 
	lwd = 1.35, 
	whisklty = 1, 
	at = c(1,2,4,5,7,8,10,11,13,14,16,17,19,20), 
	notch = FALSE
	)
	title("SAS pairs in all species", adj = 0.5, line = 1.3)
	rug(x = c(3, 6, 9, 12, 15, 18), ticksize = -0.08, side = 1, lwd=1.35, col="gray60") #x-axis ticks
	abline(v = c(3, 6, 9, 12, 15, 18), col="gray60")
	box(lwd = 1.35)
	axis(side=2, lwd = 1.35, las = 2)
	text(x= 1.5, y= 1.15, labels= "p<1e-50", col= "black", cex=1) #ATH p-value
	text(x= 4.5, y= 1.15, labels= "p<1e-30", col= "black", cex=1) #AL p-value
	text(x= 7.5, y= 1.15, labels= "p<1e-50", col= "black", cex=1) #CR p-value
	text(x= 10.5, y= 1.15, labels= "p<1e-50", col= "black", cex=1) #ES p-value
	text(x= 13.5, y= 1.15, labels= "p<1e-15", col= "black", cex=1) #TH p-value
	text(x= 16.5, y= 1.15, labels= "p<1e-50", col= "black", cex=1) #MT p-value
	text(x= 19.5, y= 1.15, labels= "p<1e-50", col= "black", cex=1) #BD p-value
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
	bty='n', horiz=TRUE, fill=c("darkseagreen", "gold2"), cex=1.1, x.intersp = 0.5)
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
		width = 2900, height = 4000, res = 825)
	par(mar = c(5.725, 4.5, 4, 2.5))
	boxplot(threshold_05, threshold_05_2, threshold_2_5, threshold_greater5,
		ylim = c(-1.1, 1.1), 
		yaxt='n', 
		cex.lab = 1.1, 
		las = 1,
		cex.axis = 1.1, #adapt size of axis labels
		xlab = "cis-NAT expression (TPM)", 
		ylab = "Pearson ρ", 
		col = c("gold2", "orange1", "orange2", "orange3"), 
		boxwex = 0.71, 
		lwd = 1.35, 
		whisklty = 1, 
		at = c(1,2,3,4), 
		pars = list(outcol = "gray50"),
		notch = FALSE
		)
		title(title_plot, adj = 0.5, line = 1.25)
		box(lwd = 1.35)
		axis(side=2, lwd = 1.35, las = 2)
		text(x= 1, y= -1.05, labels= n_values_05, col= "gray40", cex=0.97) #threshold_>0.5
		text(x= 2, y= -1.05, labels= n_values_05_2, col= "gray40", cex=0.97) #threshold_0.5-2
		text(x= 3, y= -1.05, labels= n_values_2_5, col= "gray40", cex=0.97) #threshold_2-5
		text(x= 4, y= -1.05, labels= n_values_greater5, col= "gray40", cex=0.97) #threshold_greater5
		mtext('>0.5', side=1, line=1, at=1)
		mtext('>0.5-2', side=1, line=1, at=2)
		mtext('>2-5', side=1, line=1, at=3)
		mtext('>5', side=1, line=1, at=4)
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
		width = 2900, height = 4000, res = 825)
	par(mar = c(5.725, 4.5, 4, 2.5))
	boxplot(threshold_05, threshold_05_2, threshold_2_5, threshold_greater5,
		ylim = c(-1.1, 1.1), 
		yaxt='n', 
		cex.lab = 1.1, 
		las = 1,
		cex.axis = 1.1, #adapt size of axis labels
		xlab = "cis-NAT expression (TPM)", 
		col = c("gold2", "orange1", "orange2", "orange3"), 
		boxwex = 0.71, 
		lwd = 1.35, 
		whisklty = 1, 
		at = c(1,2,3,4), 
		pars = list(outcol = "gray50"),
		notch = FALSE
		)
		title(title_plot, adj = 0.5, line = 1.25)
		box(lwd = 1.35)
		axis(side=2, lwd = 1.35, las = 2)
		text(x= 1, y= -1.05, labels= n_values_05, col= "gray40", cex=0.97) #threshold_>0.5
		text(x= 2, y= -1.05, labels= n_values_05_2, col= "gray40", cex=0.97) #threshold_0.5-2
		text(x= 3, y= -1.05, labels= n_values_2_5, col= "gray40", cex=0.97) #threshold_2-5
		text(x= 4, y= -1.05, labels= n_values_greater5, col= "gray40", cex=0.97) #threshold_greater5
		mtext('>0.5', side=1, line=1, at=1)
		mtext('>0.5-2', side=1, line=1, at=2)
		mtext('>2-5', side=1, line=1, at=3)
		mtext('>5', side=1, line=1, at=4)
		par(xpd=TRUE)
	dev.off()
}

