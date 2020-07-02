# This script loads and analysis the data statistics and expression tables for protein-coding 
# genes, isoforms, lncRNAs and circRNAs and generates the plots for the DevSeq transcriptome  
# single-species expression figures


#------------------- Load packages, set directories and read sample tables ---------------------


# Install and load packages
if (!require(dplyr)) install.packages('dplyr')
library(dplyr)
if (!require(ggplot2)) install.packages('ggplot2')
library(ggplot2)
if (!require(data.table)) install.packages('data.table')
library(data.table)
if (!require(mgcv)) install.packages('mgcv')
library(mgcv)
if (!require(grid)) install.packages('grid')
library(grid)
if (!require(gtable)) install.packages('gtable')
library(gtable)
if (!require(scales)) install.packages('scales')
library(scales)
if (!require(factoextra)) install.packages('factoextra')
library(factoextra)
if (!require(dendextend)) install.packages('dendextend')
library(dendextend)


# Set file path and input files
in_dir_stats <- "/Volumes/User/Shared/Christoph_manuscript/DevSeq_paper/Analysis/Analysis_2019/A_thaliana_gene_exression_map/20200401_CS_exprGenes/output/mapping_statistics"
in_dir_expr_genes <- "/Volumes/User/Shared/Christoph_manuscript/DevSeq_paper/Analysis/Analysis_2019/A_thaliana_gene_exression_map/20200401_CS_exprGenes/output/expr_genes"
out_dir <- "/Volumes/User/Shared/Christoph_manuscript/DevSeq_paper/Analysis/Analysis_2019/A_thaliana_gene_exression_map/20200401_CS_exprGenes"


# Read all csv files in input file path
readTable <- function(path, pattern = "*.csv") {
    files = list.files(path, pattern, full.names = TRUE)
    lapply(files, function(x) read.table(x, sep=";", dec=".", header=TRUE, stringsAsFactors=FALSE))
}

stats_tables <- readTable(in_dir_stats)
expr_genes_tables <- readTable(in_dir_expr_genes)


# Get file names and save them in character vector
stats_table_list <- as.character(list.files(in_dir_stats, pattern = "*.csv"))
stats_table_names <- gsub('\\.csv$', '', stats_table_list)
expr_genes_table_list <- as.character(list.files(in_dir_expr_genes, pattern = "*.csv"))
expr_genes_table_names <- gsub('\\.csv$', '', expr_genes_table_list)


# Change data frame names in list
names(stats_tables) <- stats_table_names
list2env(stats_tables, envir = .GlobalEnv)
names(expr_genes_tables) <- expr_genes_table_names
list2env(expr_genes_tables, envir = .GlobalEnv)


# Create "plots" folder in /out_dir/output/plots
if (!dir.exists(file.path(out_dir, "output", "plots"))) 
	dir.create(file.path(out_dir, "output", "plots"), recursive = TRUE)




#----------------------------------- Plotting stats data ------------------------------------


# Prepare data for ggplot
prepareStats <- function(x) {

	number_values <- (nrow(x))

	classT_name = as.data.frame(rep(c("trimmed"), times = number_values))
	names(classT_name) <- "class"
	trimmed <- cbind(classT_name, select(x, Trimmed))
	names(trimmed)[2] <- "reads"

	classM_name = as.data.frame(rep(c("mapped"), times = number_values))
	names(classM_name) <- "class"
	mapped <- cbind(classM_name, select(x, Mapped))
	names(mapped)[2] <- "reads"

	classD_name = as.data.frame(rep(c("dedupl."), times = number_values))
	names(classD_name) <- "class"
	deduplicated <- cbind(classD_name, select(x, Deduplicated))
	names(deduplicated)[2] <- "reads"

	mapStats <- rbind(trimmed, mapped, deduplicated)
	return(mapStats)
}

ATH_stats_df <- prepareStats(ATH_stats)
non_ATH_stats_df <- prepareStats(non_ATH_stats)



# Make stats violin plot for ATH
makePlotStatsATH <- function(data, lim_y, medw, plot_title) {

	fname <- sprintf('%s.jpg', paste(deparse(substitute(data)), sep="_"))

	dedupl <- subset(data, class=="dedupl.")
	n_dedupl <- paste("n=", nrow(dedupl), sep="")
	total_dedupl <- paste0(round(sum(as.numeric(dedupl[,2]))/1e9,2),"B")
	total_dedupl = paste(n_dedupl, total_dedupl, sep="\n")

	ylabels = function(l) {paste0(round(l/1e6,1),"M")}

	p <- ggplot(data, aes(x=class, y=reads, fill=class)) + 
		 geom_violin(trim=TRUE, width = 1.5, size=1.25, scale="area", color="gray15") + 
		 geom_boxplot(aes(x=class, y=reads),alpha=0, color="gray15", fill="white", width=medw, 
		 	size=0.0005, fatten = 5000) + 
		 scale_y_continuous(limits = c(0,lim_y), expand = c(0, 0), 
		 	labels = function(l) { 
		 		ifelse(l==0, paste0(round(l/1e6,1)),paste0(round(l/1e6,1),"M"))
		 	}) + 
		 annotate("rect", xmin=0.25, xmax=3.85, ymin=0, ymax=lim_y, fill="white", alpha=0, 
		 	color="black", size=1.35) + 
		 annotate("text", x = 2.7, y = Inf, hjust = 0, vjust = 1.55, size=6.5, label = total_dedupl)

	q <- p + scale_fill_manual(values=c("#b2b2b2", "#d8a900", "#35bceb")) + theme_minimal() + 
	xlab("") + ylab("PE reads") + ggtitle(plot_title) + 
	geom_hline(yintercept=30e6, linetype="dashed", color = "red", size=1) + 
	theme(legend.position = "none", 
		text=element_text(size=20.75), 
  		axis.ticks.length = unit(.3, "cm"),
  		axis.ticks = element_line(colour = "gray15", size = 0.7), 
  		axis.title.y = element_text(colour = "black", size=20, 
  			margin = margin(t = 0, r = 11, b = 0, l = 0)), 
  		axis.text.x = element_text(colour = "black", size=18.5, angle=90, 
  			margin = margin(t = 4, r = 0, b = 1, l = 0), hjust = 1, vjust = 0.5),
  		axis.text.y = element_text(colour = "black", margin = margin(t = 0, r = 4, b = 0, l = 1)), 
  		plot.title = element_text(colour = "black", size=22, 
  			margin = margin(t = 18, r = 0, b = 16.5, l = 0), hjust = 0.5), 
  		plot.margin = unit(c(7.0, 30, 14.1, 5.1), "points"))
	

  	ggsave(file = file.path(out_dir, "output", "plots", fname), plot = q,
		scale = 1, width = 5.0, height = 6.95, units = c("in"), 
		dpi = 600, limitsize = FALSE)
}

makePlotStatsATH(data=ATH_stats_df, lim_y=212000000, medw = 0.277, plot_title="A.thaliana") 
# 1 data point for trimmed raw reads above lim_y



# Make stats violin plot for other species
makePlotStatsOS <- function(data, lim_y, medw, plot_title) {

	fname <- sprintf('%s.jpg', paste(deparse(substitute(data)), sep="_"))

	dedupl <- subset(data, class=="dedupl.")
	n_dedupl <- paste("n=", nrow(dedupl), sep="")
	total_dedupl <- paste0(round(sum(as.numeric(dedupl[,2]))/1e9,2),"B")
	total_dedupl = paste(n_dedupl, total_dedupl, sep="\n")

	# separate outliers with low reads from violin plot data and plot them individually as dots
	data_wo_outl <- data[c(-505:-507),]
	data_outl <- data[c(505:507),]

	p <- ggplot(data_wo_outl, aes(x=class, y=reads, fill=class)) + 
		 geom_violin(trim=TRUE, width = 1.5, size=1.25, scale="area", color="gray15") + 
		 geom_boxplot(data=data, aes(x=class, y=reads),alpha=0, color="gray15", fill="white", width=medw, 
		 	size=0.0005, fatten = 5000) +
		 geom_point(aes(x=3, y=data_outl[1,2]), shape=21, colour="gray35", size=2.25, fill="white", stroke=2) + 
		 geom_point(aes(x=3, y=data_outl[2,2]), shape=21, colour="gray35", size=2.25, fill="white", stroke=2) + 
		 geom_point(aes(x=3, y=data_outl[3,2]), shape=21, colour="gray35", size=2.25, fill="white", stroke=2) + 
		 scale_y_continuous(limits = c(0,lim_y), expand = c(0, 0), 
		 	labels = function(l) { 
		 		ifelse(l==0, paste0(round(l/1e6,1)),paste0(round(l/1e6,1),"M"))
		 	}) + 
		 annotate("rect", xmin=0.25, xmax=3.85, ymin=0, ymax=lim_y, fill="white", alpha=0, 
		 	color="black", size=1.35) + 
		 annotate("text", x = 2.7, y = Inf, hjust = 0, vjust = 1.55, size=6.5, label = total_dedupl)

	q <- p + scale_fill_manual(values=c("#b2b2b2", "#d8a900", "#35bceb")) + theme_minimal() + 
	xlab("") + ylab("PE reads") + ggtitle(plot_title) + 
	geom_hline(yintercept=30e6, linetype="dashed", color = "red", size=1) + 
	theme(legend.position = "none", 
		text=element_text(size=20.75), 
  		axis.ticks.length = unit(.3, "cm"),
  		axis.ticks = element_line(colour = "gray15", size = 0.7), 
  		axis.title.y = element_text(colour = "black", size=20, 
  			margin = margin(t = 0, r = 11, b = 0, l = 0)), 
  		axis.text.x = element_text(colour = "black", size=18.5, angle=90, 
  			margin = margin(t = 4, r = 0, b = 1, l = 0), hjust = 1, vjust = 0.5), 
  		axis.text.y = element_text(colour = "black", margin = margin(t = 0, r = 4, b = 0, l = 1)), 
  		plot.title = element_text(colour = "black", size=22, 
  			margin = margin(t = 18, r = 0, b = 14.25, l = 0), hjust = 0.5), 
  		plot.margin = unit(c(7.0, 2, 14.1, 33.1), "points"))
	

  	ggsave(file = file.path(out_dir, "output", "plots", fname), plot = q,
		scale = 1, width = 5.0, height = 6.95, units = c("in"), 
		dpi = 600, limitsize = FALSE)
}

makePlotStatsOS(data=non_ATH_stats_df, lim_y=226000000, medw = 0.415, plot_title="Other species")
# 10 data point for trimmed raw reads above lim_y



# Make Deduplicated Read Replicates
makeRepl <- function(x) {

	repl_names_ATH <- data.frame(c("root.1","root.2","root.3","hypocotyl.1","hypocotyl.2",
		"hypocotyl.3","leaf.1","leaf.2","leaf.3","apex.veg.1","apex.veg.2","apex.veg.3",
		"apex.infl.1","apex.infl.2","apex.infl.3","flower.1","flower.2","flower.3","stamen.1",
		"stamen.2","stamen.3","pollen.1","pollen.2","pollen.3","carpel.1","carpel.2","carpel.3"))
	names(repl_names_ATH) <- "Sample_repl"
	repl_names_non_ATH <- c("root.1","root.2","root.3","hypocotyl.1","hypocotyl.2","hypocotyl.3",
		"leaf.1","leaf.2","leaf.3","apex.veg.1","apex.veg.2","apex.veg.3","apex.infl.1",
		"apex.infl.2","apex.infl.3","flower.1","flower.2","flower.3","pollen.1","pollen.2",
		"pollen.3","carpel.1","carpel.2","carpel.3","stamen.1","stamen.2","stamen.3")

	repl_names_non_ATH <- as.data.frame(rep(repl_names_non_ATH, times = 6))
	names(repl_names_non_ATH) <- "Sample_repl"
	all_repl_non_ATH <- rbind(repl_names_ATH,repl_names_non_ATH)

	repl_df <- cbind(x, all_repl_non_ATH)
	return(repl_df)
}

comp_stats_df <- makeRepl(comp_stats)



# Plot number of deduplicated reads for each species
plotDedupReads <- function(data, plot_title) {

	fname <- sprintf('%s.jpg', paste(deparse(substitute(data)), sep="_"))

	level_order <- c("root.1","root.2","root.3","hypocotyl.1","hypocotyl.2","hypocotyl.3",
		"leaf.1","leaf.2","leaf.3","apex.veg.1","apex.veg.2","apex.veg.3","apex.infl.1",
		"apex.infl.2","apex.infl.3","flower.1","flower.2","flower.3","carpel.1","carpel.2",
		"carpel.3","stamen.1","stamen.2","stamen.3","pollen.1","pollen.2","pollen.3")

	species_order <- c("ATH","AL","CR","ES","TH","MT","BD")

	p <- ggplot(data, aes(x = factor(Sample_repl, level= level_order), y = Deduplicated, color = Species, group = Species)) + 
	geom_line(aes(x = factor(Sample_repl, level= level_order)), size=1.5) + 
  	geom_point(aes(x = factor(Sample_repl, level= level_order)), size=3.25) + 
  	scale_y_continuous(limits = c(0,7.4e7), expand = c(0, 0), 
		 	labels = function(l) { 
		 		ifelse(l==0, paste0(round(l/1e6,1)),paste0(round(l/1e6,1),"M"))
		 	}) + 
  	annotate("rect", xmin=0.25, xmax=27.85, ymin=0, ymax=7.4e7, fill="white", alpha=0, 
		 	color="black", size=0.7) + 
  	labs(color="Species")

	q <- p + ggtitle(plot_title) + theme_bw() + xlab("") + ylab("PE dedupl. reads") + 
	scale_color_manual(values=c("#dca207","#a8a8a8","#ea6965","#46ae12","#1fac7b","#967cee","#36a5d8"), breaks=species_order) + 
		guides(colour = guide_legend(nrow = 1)) + 
  		theme(text=element_text(size=21), 
  		axis.ticks.length = unit(.3, "cm"),
  		axis.ticks = element_line(colour = "gray15", size = 0.7), 
  		axis.title.y = element_text(colour = "black", size=21, 
  			margin = margin(t = 0, r = 11, b = 0, l = 2.2)), 
  		axis.text.x = element_text(colour = "black", size=18, angle=90, 
  			margin = margin(t = 2.5, r = 0, b = 1, l = 0), hjust = 1, vjust = 0.5), 
  		axis.text.y = element_text(colour = "black", margin = margin(t = 0, r = 5, b = 0, l = 1)), 
  		plot.title = element_text(colour = "black", size=22, 
  			margin = margin(t = 18, r = 0, b = 14.95, l = 0), hjust = 0.5), 
  		plot.margin = unit(c(5.5, 2, 3, 4.5), "points"),
		legend.position = c(0.337, 0.115),
		legend.background = element_rect(fill = NA),
		legend.key = element_rect(fill = NA),
		legend.title = element_text(colour = "black", size=19.5, face ="bold"),
		legend.text=element_text(size=19.5), 
		legend.spacing.x = unit(0.25, 'cm'),
		legend.key.size = unit(0.775, "cm"),
  		panel.border = element_rect(colour = "black", fill=NA, size=0.5))

  	png("NUL")
	r <- ggplotGrob(q)
	r$layout$clip[r$layout$name=="panel"] <- "off"

	ggsave(file = file.path(out_dir, "output", "plots", fname), plot = r,
		scale = 1, width = 9.8, height = 7.09, units = c("in"), 
		dpi = 600, limitsize = FALSE)
}

plotDedupReads(data=comp_stats_df, plot_title="Comparative samples")




#-------------------------- Plotting expressed genes data for ATH ---------------------------


# Prepare data for ggplot
prepareExprGenes <- function(biotype = c("coding", "NAT", "lincRNA"), th_0_01, th_0_05, th_0_1, 
	th_0) {

	df_names <- c("Detailed_name" , "Sample" , "Threshold", "Expressed")

	if (is.element("coding", biotype)) {
		data <- cbind(th_0_01[1,3:ncol(th_0_01)], th_0_05[1,3:ncol(th_0_05)], th_0_1[1,3:ncol(th_0_1)], th_0[1,3:ncol(th_0)])
	} else if (is.element("NAT", biotype)) {
		data <- cbind(th_0_01[2,3:ncol(th_0_01)], th_0_05[2,3:ncol(th_0_05)], th_0_1[2,3:ncol(th_0_1)], th_0[2,3:ncol(th_0)])
	} else  if (is.element("lincRNA", biotype)) {
		data <- cbind(th_0_01[3,3:ncol(th_0_01)], th_0_05[3,3:ncol(th_0_05)], th_0_1[3,3:ncol(th_0_1)], th_0[3,3:ncol(th_0)])
	}

	colnames(data) <- NULL
	data <- as.data.frame(t(data))

	detailed_sample_name <- names(th_0)[3:ncol(th_0)]
	detailed_sample_name <- as.data.frame(rep(detailed_sample_name, times=4))

	sample_names <- c("root tip 5d", "root m.zone", "whole root 5", "whole root 7", "whole rt.14d", 
		"whole rt.21d", "hypocotyl10", "3.internode", "2.internode", "1.internode", "cotyledons", 
		"leaf 1+2 7d", "leaf 1.2 10d", "leaf petiole", "leaf tip 10d", "leaf 5.6 17d", "leaf 910 27d", "leaf sen.35d", 
		"cauline leaf", "apex veg.7d", "apex veg.10", "apex veg.14", "apex inf 21d", "apex inf clv1", 
		"apex inf 28d", "flower stg.9", "flower 10.11", "flower st12", "flower st15", 
		"sepals st12", "sepals st15", "petals st12", "petals st15", "stamen st12", 
		"stamen st15", "pollen mat.", "carpel st12e", "carpel st12l", "carpels st15", 
		"fruit stg.16", "fruit stg.17a", "seeds st16", "seeds st17a", "seeds st18")
	sample_names <- as.data.frame(rep(sample_names, times=4))

	threshold <- c("0.01", "0.05", "0.1", "0.5")
	threshold <- as.data.frame(rep(threshold, each=44))

	expr_df <- cbind(detailed_sample_name,sample_names, threshold, data)
	colnames(expr_df) <- df_names

	return(expr_df)
}


expr_coding_genes_ATH <- prepareExprGenes(biotype = "coding", th_0_01 = ATH_expr_genes_0.01, 
	th_0_05 = ATH_expr_genes_0.05, th_0_1 = ATH_expr_genes_0.1, th_0 = ATH_expr_genes_0)
expr_NATs_ATH <- prepareExprGenes(biotype = "NAT", th_0_01 = ATH_expr_genes_0.01, 
	th_0_05 = ATH_expr_genes_0.05, th_0_1 = ATH_expr_genes_0.1, th_0 = ATH_expr_genes_0)
expr_lincRNAs_ATH <- prepareExprGenes(biotype = "lincRNA", th_0_01 = ATH_expr_genes_0.01, 
	th_0_05 = ATH_expr_genes_0.05, th_0_1 = ATH_expr_genes_0.1, th_0 = ATH_expr_genes_0)



# Plot number of expressed genes at different thresholds for ATH
plotExprGenes <- function(data, plot_title, biotype = c("coding","NAT","linc"), texpr) {

	fname <- sprintf('%s.jpg', paste(deparse(substitute(data)), sep="_"))

	total_expr <- paste("total:", texpr, "at 0.05" , sep=" ")

	if (is.element("coding", biotype)) {
		breaksY <- c(1.5e4,2e4,2.5e4)
		pltymin <- 1.0e4
		pltymax <- 2.78e4
		xtepos <- 20.35
		y_margin <- margin(t = 0, r = 12, b = 0, l = 0)

	} else if (is.element("NAT", biotype)) {
		breaksY <- c(1e3,2e3,3e3)
		pltymin <- 2.0e2
		pltymax <- 3.68e3
		xtepos <- 32.15
		y_margin <- margin(t = 0, r = 12, b = 0, l = 10.5)

	} else if (is.element("linc", biotype)) {
		breaksY <- c(5e2,1e3,1.5e3)
		pltymin <- 0.35e2
		pltymax <- 1.585e3
		xtepos <- 32.15
		y_margin <- margin(t = 0, r = 6.85, b = 0, l = 0)
	}

	level_order <- c("root tip 5d", "root m.zone", "whole root 5", "whole root 7", "whole rt.14d", 
		"whole rt.21d", "hypocotyl10", "3.internode", "2.internode", "1.internode", "cotyledons", 
		"leaf 1+2 7d", "leaf 1.2 10d", "leaf petiole", "leaf tip 10d", "leaf 5.6 17d", "leaf 910 27d", "leaf sen.35d", 
		"cauline leaf", "apex veg.7d", "apex veg.10", "apex veg.14", "apex inf 21d", "apex inf clv1", 
		"apex inf 28d", "flower stg.9", "flower 10.11", "flower st12", "flower st15", 
		"sepals st12", "sepals st15", "petals st12", "petals st15", "stamen st12", 
		"stamen st15", "pollen mat.", "carpel st12e", "carpel st12l", "carpels st15", 
		"fruit stg.16", "fruit stg.17a", "seeds st16", "seeds st17a", "seeds st18")

	p <- ggplot(data, aes(x = factor(Sample, level= level_order), y = Expressed, color = Threshold, group = Threshold)) + 

	geom_line(aes(x = factor(Sample, level= level_order)), size=1.7) + 
	scale_y_continuous(limits = c(pltymin,pltymax), breaks = breaksY, expand = c(0, 0), 
		 	labels = function(l) { 
		 		ifelse(l==0, paste0(round(l/1e3,1)),paste0(round(l/1e3,1),"K"))
		 	}) +
  	annotate("rect", xmin=0.25, xmax=44.75, ymin=pltymin, ymax=pltymax, fill="white", alpha=0, 
		 	color="black", size=0.7) + 
  	annotate("rect", xmin=0.25, xmax=6.5, ymin=pltymin, ymax=pltymax, fill="#747474", alpha=0.34) + 
  	annotate("rect", xmin=10.5, xmax=19.5, ymin=pltymin, ymax=pltymax, fill="#0fc94d", alpha=0.34) + 
  	annotate("rect", xmin=25.5, xmax=29.5, ymin=pltymin, ymax=pltymax, fill="#747474", alpha=0.34) + 
  	annotate("rect", xmin=38.5, xmax=44.75, ymin=pltymin, ymax=pltymax, fill="#db4a10", alpha=0.34) +
  	geom_line(aes(x = factor(Sample, level= level_order)), size=1.55) + 
  	annotate("text", x = xtepos, y = Inf, hjust = 0, vjust = 22.75, size=7.01, label = total_expr) + 
  	annotate("text", x = 1.675, y = Inf, hjust = 0, vjust = 21.075, size=7.01, label = "Threshold", fontface = 2) + 
  	annotate("text", x = 2.15, y = Inf, hjust = 0, vjust = 2.4, size=7.25, label = "root") + 
  	annotate("text", x = 7.0, y = Inf, hjust = 0, vjust = 2.4, size= 7.25, label = "stem") + 
  	annotate("text", x = 13.8, y = Inf, hjust = 0, vjust = 2.4, size= 7.25, label = "leaf") + 
  	annotate("text", x = 20.95, y = Inf, hjust = 0, vjust = 2.4, size= 7.25, label = "apex") + 
  	annotate("text", x = 26.2, y = Inf, hjust = 0, vjust = 2.4, size= 7.25, label = "flow.") + 
  	annotate("text", x = 30.55, y = Inf, hjust = 0, vjust = 2.4, size= 7.25, label = "floral organ") + 
  	annotate("text", x = 40.4, y = Inf, hjust = 0, vjust = 2.4, size= 7.25, label = "fruit") + 
  	labs(color="")

	q <- p + ggtitle(plot_title) + theme_bw() + xlab("") + ylab("Number of Genes") + 
	scale_color_manual(values=c("gray45","#ea6965","#967cee","#dca207")) + 
		guides(colour = guide_legend(nrow = 1)) + 
  		theme(text=element_text(size=23.5), 
  		panel.grid.major = element_line(colour = "white"), 
  		panel.grid.minor = element_line(colour = "white"),  
  		axis.ticks.length = unit(.3, "cm"),
  		axis.ticks = element_line(colour = "gray15", size = 0.7), 
  		axis.title.y = element_text(colour = "black", size=21, 
  			margin = y_margin), 
  		axis.text.x = element_text(colour = "black", size=16.5, angle=90, 
  			margin = margin(t = 2.5, r = 0, b = 1, l = 0), hjust = 1, vjust = 0.5), 
  		axis.text.y = element_text(colour = "black", margin = margin(t = 0, r = 4, b = 0, l = 1)), 
  		plot.title = element_text(colour = "black", size=23.5, 
  			margin = margin(t = 17, r = 0, b = 16, l = 0), hjust = 0.5), 
  		plot.margin = unit(c(0, 2, 0, 1), "points"),
		legend.position = c(0.21, 0.11275),
		legend.title = element_text(colour = "black", size=20, face ="bold"),
		legend.text = element_text(size=20), 
		legend.key.size = unit(0.775, "cm"),
		legend.key.height = unit(0.4, "cm"),
		legend.background = element_rect(fill = NA),
		legend.key = element_rect(fill = NA),
  		panel.border = element_rect(colour = "black", fill=NA, size=0.5))

  	png("NUL")
	r <- ggplotGrob(q)
	r$layout$clip[r$layout$name=="panel"] <- "off"

	ggsave(file = file.path(out_dir, "output", "plots", fname), plot = r,
		scale = 1, width = 10.25, height = 7.3, units = c("in"), 
		dpi = 600, limitsize = FALSE)
}


plotExprGenes(data=expr_coding_genes_ATH, plot_title="Expressed protein-coding genes in A.thaliana", biotype = "coding", texpr=ATH_expr_genes_0.05[1,2])
plotExprGenes(data=expr_NATs_ATH, plot_title="Expressed NATs in A.thaliana", biotype = "NAT", texpr=ATH_expr_genes_0.05[2,2])
plotExprGenes(data=expr_lincRNAs_ATH, plot_title="Expressed lincRNAs in A.thaliana", biotype = "linc", texpr=ATH_expr_genes_0.05[3,2])




#----------------------------- Plotting replicate correlations ------------------------------


# Prepare data for ggplot
prepareReplStats <- function(ATH,AL,CR,ES,TH,MT,BD) {

	number_values_ATH <- (ncol(ATH))
	number_values_AL <- (ncol(AL))
	number_values_other <- (ncol(CR))

	class_ATH_key = as.data.frame(rep(c("ATH"), times = number_values_ATH))
	names(class_ATH_key) <- "Species"
	class_ATH <- cbind(class_ATH_key, as.data.frame(t(ATH)))
	names(class_ATH)[2] <- "Correlation"
	class_AL_key = as.data.frame(rep(c("AL"), times = number_values_AL))
	names(class_AL_key) <- "Species"
	class_AL <- cbind(class_AL_key, as.data.frame(t(AL)))
	names(class_AL)[2] <- "Correlation"
	class_CR_key = as.data.frame(rep(c("CR"), times = number_values_other))
	names(class_CR_key) <- "Species"
	class_CR <- cbind(class_CR_key, as.data.frame(t(CR)))
	names(class_CR)[2] <- "Correlation"
	class_ES_key = as.data.frame(rep(c("ES"), times = number_values_other))
	names(class_ES_key) <- "Species"
	class_ES <- cbind(class_ES_key, as.data.frame(t(ES)))
	names(class_ES)[2] <- "Correlation"
	class_TH_key = as.data.frame(rep(c("TH"), times = number_values_other))
	names(class_TH_key) <- "Species"
	class_TH <- cbind(class_TH_key, as.data.frame(t(TH)))
	names(class_TH)[2] <- "Correlation"
	class_MT_key = as.data.frame(rep(c("MT"), times = number_values_other))
	names(class_MT_key) <- "Species"
	class_MT <- cbind(class_MT_key, as.data.frame(t(MT)))
	names(class_MT)[2] <- "Correlation"
	class_BD_key = as.data.frame(rep(c("BD"), times = number_values_other))
	names(class_BD_key) <- "Species"
	class_BD <- cbind(class_BD_key, as.data.frame(t(BD)))
	names(class_BD)[2] <- "Correlation"

	repStats <- rbind(class_ATH, class_AL, class_CR, class_ES, class_TH, class_MT, class_BD)
	return(repStats)
}

all_spec_repl_df <- prepareReplStats(ATH=ATH_repl_corr_tpm_0.05, AL=AL_repl_corr_tpm_0.05, 
	CR=CR_repl_corr_tpm_0.05, ES=ES_repl_corr_tpm_0.05, TH=TH_repl_corr_tpm_0.05, 
	MT=MT_repl_corr_tpm_0.05, BD=BD_repl_corr_tpm_0.05)



# Make replicate correlation plot
makePlotReplCorr <- function(data, plot_title) {

	fname <- sprintf('%s.jpg', paste(deparse(substitute(data)), sep="_"))

	p <- ggplot(data, aes(x=Species, y=Correlation, fill=Species)) + 
	     stat_boxplot(geom ='errorbar', width = 0.45, size=1.0, color="gray15") + 
		 geom_boxplot(width = 0.75, size=1.0, color="gray15", outlier.shape = 21, 
		 	outlier.size = 2.5, outlier.stroke = 1.5, outlier.fill = NA, outlier.color="gray35") + 
		 scale_y_continuous(limits = c(0.9658,1.0005), expand = c(0, 0)) + 
		 annotate("rect", xmin=0.35, xmax=7.65, ymin=0.9658, ymax=1.0005, fill="white", alpha=0, 
		 	color="black", size=1.35)

	q <- p + scale_fill_manual(values=c("#b2b2b2","#dca207","#46ae12","#1fac7b","#36a5d8","#967cee","#ea6965")) + 
	theme_minimal() + 
	xlab("Species") + ylab("Pearson's r") + ggtitle(plot_title) + 
	theme(legend.position = "none", 
		text=element_text(size=20.75), 
  		axis.ticks.length = unit(.3, "cm"),
  		axis.ticks = element_line(colour = "gray15", size = 0.7), 
  		axis.title.x = element_text(colour = "black", size=20, 
  			margin = margin(t = 17.5, r = 0, b = 0, l = 0)), 
  		axis.title.y = element_text(colour = "black", size=20, 
  			margin = margin(t = 0, r = 12.5, b = 0, l = 0.5)), 
  		axis.text.x = element_text(colour = "black", size=18.5, angle=0, 
  			margin = margin(t = 7.5, r = 0, b = 0, l = 0), hjust = 0.5, vjust = 0.5),
  		axis.text.y = element_text(colour = "black", margin = margin(t = 0, r = 5, b = 0, l = 0)), 
  		plot.title = element_text(colour = "black", size=22, 
  			margin = margin(t = 21.5, r = 0, b = 13.5, l = 0), hjust = 0.5), 
  		plot.margin = unit(c(0, 1.25, 69.75, 8.0), "points"))

  	ggsave(file = file.path(out_dir, "output", "plots", fname), plot = q,
		scale = 1, width = 9.1, height = 7.09, units = c("in"), 
		dpi = 600, limitsize = FALSE)
}

makePlotReplCorr(data=all_spec_repl_df, plot_title="Replicate correlations") 
# 1 data point for trimmed raw reads above lim_y




#-------------------- Plotting expressed genes data for non-ATH species ---------------------


# Prepare data for ggplot
prepareExprGenesOS <- function(biotype=c("coding","NAT","lincRNA"), species=c("AL","CR","ES","TH","MT","BD"), 
	th_0, th_0_01, th_0_05, th_0_1) {

	df_names <- c("Detailed_name" , "Sample" , "Threshold", "Expressed")

	if (is.element("coding", biotype)) {
		data <- cbind(th_0[1,3:ncol(th_0)], th_0_01[1,3:ncol(th_0_01)], th_0_05[1,3:ncol(th_0_05)], th_0_1[1,3:ncol(th_0_1)])
	} else if (is.element("NAT", biotype)) {
		data <- cbind(th_0[2,3:ncol(th_0)], th_0_01[2,3:ncol(th_0_01)], th_0_05[2,3:ncol(th_0_05)], th_0_1[2,3:ncol(th_0_1)])
	} else  if (is.element("lincRNA", biotype)) {
		data <- cbind(th_0[3,3:ncol(th_0)], th_0_01[3,3:ncol(th_0_01)], th_0_05[3,3:ncol(th_0_05)], th_0_1[3,3:ncol(th_0_1)])
	}

	colnames(data) <- NULL
	data <- as.data.frame(t(data))

	detailed_sample_name <- names(th_0)[3:ncol(th_0)]
	detailed_sample_name <- as.data.frame(rep(detailed_sample_name, times=4))

	if (is.element("BD", species)) {

		sample_names <- c("root", "mesocotyl", "leaf", "apex veg", "spiklet m", 
		"floret", "pollen", "carpel", "stamen")

	} else {

		sample_names <- c("root", "hypocotyl", "leaf", "apex veg", "apex inf", 
		"flower", "pollen", "carpel", "stamen")
	}

	sample_names <- as.data.frame(rep(sample_names, times=4))

	threshold <- c("0", "0.01", "0.05", "0.1")
	threshold <- as.data.frame(rep(threshold, each=9))

	expr_df <- cbind(detailed_sample_name,sample_names, threshold, data)
	colnames(expr_df) <- df_names

	expr_df <- expr_df[c(1:6,8,9,7,10:15,17,18,16,19:24,26,27,25,28:33,35,36,34),]

	return(expr_df)
}


expr_coding_genes_AL <- prepareExprGenesOS(biotype = "coding", species = "AL", th_0 = AL_expr_genes_0,
	th_0_01 = AL_expr_genes_0.01, th_0_05 = AL_expr_genes_0.05, th_0_1 = AL_expr_genes_0.1)
expr_NATs_AL <- prepareExprGenesOS(biotype = "NAT", species = "AL", th_0 = AL_expr_genes_0,
	th_0_01 = AL_expr_genes_0.01, th_0_05 = AL_expr_genes_0.05, th_0_1 = AL_expr_genes_0.1)
expr_lincRNAs_AL <- prepareExprGenesOS(biotype = "lincRNA", species = "AL", th_0 = AL_expr_genes_0,
	th_0_01 = AL_expr_genes_0.01, th_0_05 = AL_expr_genes_0.05, th_0_1 = AL_expr_genes_0.1)

expr_coding_genes_CR <- prepareExprGenesOS(biotype = "coding", species = "CR", th_0 = CR_expr_genes_0,
	th_0_01 = CR_expr_genes_0.01, th_0_05 = CR_expr_genes_0.05, th_0_1 = CR_expr_genes_0.1)
expr_NATs_CR <- prepareExprGenesOS(biotype = "NAT", species = "CR", th_0 = CR_expr_genes_0,
	th_0_01 = CR_expr_genes_0.01, th_0_05 = CR_expr_genes_0.05, th_0_1 = CR_expr_genes_0.1)
expr_lincRNAs_CR <- prepareExprGenesOS(biotype = "lincRNA", species = "CR", th_0 = CR_expr_genes_0,
	th_0_01 = CR_expr_genes_0.01, th_0_05 = CR_expr_genes_0.05, th_0_1 = CR_expr_genes_0.1)

expr_coding_genes_ES <- prepareExprGenesOS(biotype = "coding", species = "ES", th_0 = ES_expr_genes_0,
	th_0_01 = ES_expr_genes_0.01, th_0_05 = ES_expr_genes_0.05, th_0_1 = ES_expr_genes_0.1)
expr_NATs_ES <- prepareExprGenesOS(biotype = "NAT", species = "ES", th_0 = ES_expr_genes_0,
	th_0_01 = ES_expr_genes_0.01, th_0_05 = ES_expr_genes_0.05, th_0_1 = ES_expr_genes_0.1)
expr_lincRNAs_ES <- prepareExprGenesOS(biotype = "lincRNA", species = "ES", th_0 = ES_expr_genes_0,
	th_0_01 = ES_expr_genes_0.01, th_0_05 = ES_expr_genes_0.05, th_0_1 = ES_expr_genes_0.1)

expr_coding_genes_TH <- prepareExprGenesOS(biotype = "coding", species = "TH", th_0 = TH_expr_genes_0,
	th_0_01 = TH_expr_genes_0.01, th_0_05 = TH_expr_genes_0.05, th_0_1 = TH_expr_genes_0.1)
expr_NATs_TH <- prepareExprGenesOS(biotype = "NAT", species = "TH", th_0 = TH_expr_genes_0,
	th_0_01 = TH_expr_genes_0.01, th_0_05 = TH_expr_genes_0.05, th_0_1 = TH_expr_genes_0.1)
expr_lincRNAs_TH <- prepareExprGenesOS(biotype = "lincRNA", species = "TH", th_0 = TH_expr_genes_0,
	th_0_01 = TH_expr_genes_0.01, th_0_05 = TH_expr_genes_0.05, th_0_1 = TH_expr_genes_0.1)

expr_coding_genes_MT <- prepareExprGenesOS(biotype = "coding", species = "MT", th_0 = MT_expr_genes_0,
	th_0_01 = MT_expr_genes_0.01, th_0_05 = MT_expr_genes_0.05, th_0_1 = MT_expr_genes_0.1)
expr_NATs_MT <- prepareExprGenesOS(biotype = "NAT", species = "MT", th_0 = MT_expr_genes_0,
	th_0_01 = MT_expr_genes_0.01, th_0_05 = MT_expr_genes_0.05, th_0_1 = MT_expr_genes_0.1)
expr_lincRNAs_MT <- prepareExprGenesOS(biotype = "lincRNA", species = "MT", th_0 = MT_expr_genes_0,
	th_0_01 = MT_expr_genes_0.01, th_0_05 = MT_expr_genes_0.05, th_0_1 = MT_expr_genes_0.1)

expr_coding_genes_BD <- prepareExprGenesOS(biotype = "coding", species = "BD", th_0 = BD_expr_genes_0,
	th_0_01 = BD_expr_genes_0.01, th_0_05 = BD_expr_genes_0.05, th_0_1 = BD_expr_genes_0.1)
expr_NATs_BD <- prepareExprGenesOS(biotype = "NAT", species = "BD", th_0 = BD_expr_genes_0,
	th_0_01 = BD_expr_genes_0.01, th_0_05 = BD_expr_genes_0.05, th_0_1 = BD_expr_genes_0.1)
expr_lincRNAs_BD <- prepareExprGenesOS(biotype = "lincRNA", species = "BD", th_0 = BD_expr_genes_0,
	th_0_01 = BD_expr_genes_0.01, th_0_05 = BD_expr_genes_0.05, th_0_1 = BD_expr_genes_0.1)




# Plot number of deduplicated reads for each species
plotExprGenesOS <- function(data, plot_title, species = c("dicot", "BD"), texpr, 
	class = c("coding", "non-coding")) {

	fname <- sprintf('%s.jpg', paste(deparse(substitute(data)), sep="_"))

	total_expr <- paste("Total:", texpr, "at 0.05" , sep=" ")

	legend_pos <- c(0.3475,0.22)

	if (is.element("dicot", species)) {
		level_order <- c("root", "hypocotyl", "leaf", "apex veg", "apex inf", 
		"flower", "carpel", "stamen", "pollen")

	} else if (is.element("BD", species)) {
		level_order <- c("root", "mesocotyl", "leaf", "apex veg", "spiklet m", 
		"floret", "carpel", "stamen", "pollen")
	}

	if (is.element("coding", class) &&! is.element("BD", species)) {
		legend_title <- "Threshold"
		tlmargin <- margin(t = 16.5, r = 0, b = 16.0, l = 0)
		pltmargin <- c(0, 5, 15, 10)

	} else if (is.element("non-coding", class) &&! is.element("BD", species)) {
		legend_title <- ""
		tlmargin <- margin(t = 16.5, r = 0, b = 18.375, l = 0)
		pltmargin <- c(0, 5, 15, 4.75)

	} else if (is.element("coding", class) && is.element("BD", species)) {
		legend_title <- "Threshold"
		tlmargin <- margin(t = 16.5, r = 0, b = 16.0, l = 0)
		pltmargin <- c(0, 5, 9.55, 10)

	} else if (is.element("non-coding", class) && is.element("BD", species)) {
		legend_title <- ""
		tlmargin <- margin(t = 16.5, r = 0, b = 18.375, l = 0)
		pltmargin <- c(0, 5, 9.55, 4.75)
	}

	if (is.element("lincRNAs in MT", plot_title)) {
		pltmargin <- c(0, 5, 15, 20.7)
	}

	if (is.element("lincRNAs in CR", plot_title) | is.element("lincRNAs in ES", plot_title) 
		| is.element("lincRNAs in BD", plot_title)) {
		legend_pos <- "none"
	}

	p <- ggplot(data, aes(x = factor(Sample, level= level_order), y = Expressed, color = Threshold, group = Threshold)) + 

	geom_line(aes(x = factor(Sample, level= level_order)), size=1.55) + 
	scale_y_continuous(expand = c(0.15, 0.15), breaks= pretty_breaks(), 
		 	labels = function(l) { 
		 		ifelse(l==0, paste0(round(l/1e3,1)),paste0(round(l/1e3,1),"K"))
		 	}) +
  	geom_line(aes(x = factor(Sample, level= level_order)), size=1.55) + 
  	annotate("text", x = 0.825, y = Inf, hjust = 0, vjust = 19.5, size=6.85, label = total_expr) + 
  	annotate("text", x = 0.825, y = Inf, hjust = 0, vjust = 15.75, size=6.85, label = legend_title, fontface = 2) + 
  	labs(color="")

	q <- p + ggtitle(plot_title) + theme_bw() + xlab("") + ylab("Number of Genes") + 
	scale_color_manual(values=c("gray45","#ea6965","#967cee","#dca207")) + 
		guides(colour = guide_legend(nrow = 1)) + 
  		theme(   
  		axis.ticks.length = unit(.3, "cm"),
  		axis.ticks = element_line(colour = "gray15", size = 0.75), 
  		axis.title.y = element_text(colour = "black", size=22.5, 
  			margin = margin(t = 0, r = 14, b = 0, l = 1)), 
  		axis.text.x = element_text(colour = "black", size=19.5, angle=90, 
  			margin = margin(t = 3.5, r = 0, b = 0, l = 0), hjust = 1, vjust = 0.5), 
  		axis.text.y = element_text(colour = "black", size=19, margin = margin(t = 0, r = 3, b = 0, l = 0)), 
  		plot.title = element_text(colour = "black", size=22.5, 
  			margin = tlmargin, hjust = 0.5), 
  		plot.margin = unit(pltmargin, "points"),
		legend.position = legend_pos,
		legend.title = element_text(colour = "black", size=20, face ="bold"),
		legend.text = element_text(size=19.5), 
		legend.key.size = unit(0.775, "cm"),
		legend.key.height = unit(0.4, "cm"),
		legend.background = element_rect(fill = NA),
		legend.key = element_rect(fill = NA),
  		panel.border = element_rect(colour = "black", fill=NA, size=0.75))

  	png("NUL")
	r <- ggplotGrob(q)
	r$layout$clip[r$layout$name=="panel"] <- "off"

	ggsave(file = file.path(out_dir, "output", "plots", fname), plot = r,
		scale = 1, width = 6.2, height = 6.5, units = c("in"), 
		dpi = 600, limitsize = FALSE)
}


plotExprGenesOS(data=expr_coding_genes_AL, plot_title="Protein-coding genes in AL", species = "dicot", texpr=AL_expr_genes_0.05[1,2], class="coding")
plotExprGenesOS(data=expr_NATs_AL, plot_title="NATs in AL", species = "dicot", texpr=AL_expr_genes_0.05[2,2], class="non-coding")
plotExprGenesOS(data=expr_lincRNAs_AL, plot_title="lincRNAs in AL", species = "dicot", texpr=AL_expr_genes_0.05[3,2], class="non-coding")

plotExprGenesOS(data=expr_coding_genes_CR, plot_title="Protein-coding genes in CR", species = "dicot", texpr=CR_expr_genes_0.05[1,2], class="coding")
plotExprGenesOS(data=expr_NATs_CR, plot_title="NATs in CR", species = "dicot", texpr=CR_expr_genes_0.05[2,2], class="non-coding")
plotExprGenesOS(data=expr_lincRNAs_CR, plot_title="lincRNAs in CR", species = "dicot", texpr=CR_expr_genes_0.05[3,2], class="non-coding")

plotExprGenesOS(data=expr_coding_genes_ES, plot_title="Protein-coding genes in ES", species = "dicot", texpr=ES_expr_genes_0.05[1,2], class="coding")
plotExprGenesOS(data=expr_NATs_ES, plot_title="NATs in ES", species = "dicot", texpr=ES_expr_genes_0.05[2,2], class="non-coding")
plotExprGenesOS(data=expr_lincRNAs_ES, plot_title="lincRNAs in ES", species = "dicot", texpr=ES_expr_genes_0.05[3,2], class="non-coding")

plotExprGenesOS(data=expr_coding_genes_TH, plot_title="Protein-coding genes in TH", species = "dicot", texpr=TH_expr_genes_0.05[1,2], class="coding")
plotExprGenesOS(data=expr_NATs_TH, plot_title="NATs in TH", species = "dicot", texpr=TH_expr_genes_0.05[2,2], class="non-coding")
plotExprGenesOS(data=expr_lincRNAs_TH, plot_title="lincRNAs in TH", species = "dicot", texpr=TH_expr_genes_0.05[3,2], class="non-coding")

plotExprGenesOS(data=expr_coding_genes_MT, plot_title="Protein-coding genes in MT", species = "dicot", texpr=MT_expr_genes_0.05[1,2], class="coding")
plotExprGenesOS(data=expr_NATs_MT, plot_title="NATs in MT", species = "dicot", texpr=MT_expr_genes_0.05[2,2], class="non-coding")
plotExprGenesOS(data=expr_lincRNAs_MT, plot_title="lincRNAs in MT", species = "dicot", texpr=MT_expr_genes_0.05[3,2], class="non-coding")

plotExprGenesOS(data=expr_coding_genes_BD, plot_title="Protein-coding genes in BD", species = "BD", texpr=BD_expr_genes_0.05[1,2], class="coding")
plotExprGenesOS(data=expr_NATs_BD, plot_title="NATs in BD", species = "BD", texpr=BD_expr_genes_0.05[2,2], class="non-coding")
plotExprGenesOS(data=expr_lincRNAs_BD, plot_title="lincRNAs in BD", species = "BD", texpr=BD_expr_genes_0.05[3,2], class="non-coding")




#------------ Plotting expressed genes data for non-ATH species with facet wrap -------------


# Prepare data for ggplot
prepareExprGenesOS <- function(species=c("AL","CR","ES","TH","MT","BD"), 
	th_0, th_0_01, th_0_05, th_0_1) {


	th_0_01 <- th_0_01[,c(1:8,10:11,9)]
	th_0_05 <- th_0_05[,c(1:8,10:11,9)]
	th_0_1 <- th_0_1[,c(1:8,10:11,9)]
	th_0 <- th_0[,c(1:8,10:11,9)]


	threshold <- as.data.frame(rep(c("0.01", "0.05", "0.1", "0.5"), each=9))
	colnames(threshold) <- "Threshold"

	species <- as.data.frame(rep(c(species),144))
	colnames(species) <- "species"


	coding <- cbind(th_0_01[1,3:ncol(th_0_01)], th_0_05[1,3:ncol(th_0_05)], th_0_1[1,3:ncol(th_0_1)], th_0[1,3:ncol(th_0)])
	colnames(coding) <- NULL
	coding <- as.data.frame(t(coding))
	colnames(coding) <- "expressed"
	total_expressed <- c(th_0_01[1,2], th_0_05[1,2], th_0_1[1,2], th_0[1,2])
	total_expressed <- as.data.frame(rep(c(total_expressed),each=9))
	colnames(total_expressed) <- "total_expressed"
	transcript_class <- as.data.frame(rep(c("Protein-coding"),36))
	colnames(transcript_class) <- "class"
	cd_value_class <- cbind(coding, total_expressed, threshold, transcript_class)

	NAT <- cbind(th_0_01[2,3:ncol(th_0_01)], th_0_05[2,3:ncol(th_0_05)], th_0_1[2,3:ncol(th_0_1)], th_0[2,3:ncol(th_0)])
	colnames(NAT) <- NULL
	NAT <- as.data.frame(t(NAT))
	colnames(NAT) <- "expressed"
	total_expressed <- c(th_0_01[2,2], th_0_05[2,2], th_0_1[2,2], th_0[2,2])
	total_expressed <- as.data.frame(rep(c(total_expressed),each=9))
	colnames(total_expressed) <- "total_expressed"
	transcript_class <- as.data.frame(rep(c("NAT"),36))
	colnames(transcript_class) <- "class"
	NAT_value_class <- cbind(NAT, total_expressed, threshold, transcript_class)
	
	linc <- cbind(th_0_01[3,3:ncol(th_0_01)], th_0_05[3,3:ncol(th_0_05)], th_0_1[3,3:ncol(th_0_1)], th_0[3,3:ncol(th_0)])
	colnames(linc) <- NULL
	linc <- as.data.frame(t(linc))
	colnames(linc) <- "expressed"
	total_expressed <- c(th_0_01[3,2], th_0_05[3,2], th_0_1[3,2], th_0[3,2])
	total_expressed <- as.data.frame(rep(c(total_expressed),each=9))
	colnames(total_expressed) <- "total_expressed"
	transcript_class <- as.data.frame(rep(c("lincRNA"),36))
	colnames(transcript_class) <- "class"
	linc_value_class <- cbind(linc, total_expressed ,threshold, transcript_class)
	
	# Add linc data as dummy data for circRNAs
	# Change this once real data is available!
	circ <- cbind(th_0_01[3,3:ncol(th_0_01)], th_0_05[3,3:ncol(th_0_05)], th_0_1[3,3:ncol(th_0_1)], th_0[3,3:ncol(th_0)])
	colnames(circ) <- NULL
	circ <- as.data.frame(t(circ))
	colnames(circ) <- "expressed"
	total_expressed <- c(th_0_01[3,2], th_0_05[3,2], th_0_1[3,2], th_0[3,2])
	total_expressed <- as.data.frame(rep(c(total_expressed),each=9))
	colnames(total_expressed) <- "total_expressed"
	transcript_class <- as.data.frame(rep(c("circRNA"),36))
	colnames(transcript_class) <- "class"
	circ_value_class <- cbind(circ, total_expressed, threshold, transcript_class)
	

	detailed_sample_name <- names(th_0)[3:ncol(th_0)]
	detailed_sample_name <- as.data.frame(rep(detailed_sample_name, times=16))
	colnames(detailed_sample_name) <- "detailed_sample_name"

	sample_names <- c("root", "hypocotyl", "leaf", "apex veg", "apex inf", 
		"flower", "carpel", "stamen", "pollen")

	sample_names <- as.data.frame(rep(sample_names, times=16))
	colnames(sample_names) <- "sample_names"


	expr_df <- rbind(cd_value_class, NAT_value_class, linc_value_class, circ_value_class) 
	expr_df_ext <- cbind(expr_df, species, sample_names, detailed_sample_name)

	return(expr_df_ext)
}


# Get expressed genes for each species
expr_genes_AL <- prepareExprGenesOS(species = "AL", th_0 = AL_expr_genes_0,
	th_0_01 = AL_expr_genes_0.01, th_0_05 = AL_expr_genes_0.05, th_0_1 = AL_expr_genes_0.1)

expr_genes_CR <- prepareExprGenesOS(species = "CR", th_0 = CR_expr_genes_0,
	th_0_01 = CR_expr_genes_0.01, th_0_05 = CR_expr_genes_0.05, th_0_1 = CR_expr_genes_0.1)

expr_genes_ES <- prepareExprGenesOS(species = "ES", th_0 = ES_expr_genes_0,
	th_0_01 = ES_expr_genes_0.01, th_0_05 = ES_expr_genes_0.05, th_0_1 = ES_expr_genes_0.1)

expr_genes_TH <- prepareExprGenesOS(species = "TH", th_0 = TH_expr_genes_0,
	th_0_01 = TH_expr_genes_0.01, th_0_05 = TH_expr_genes_0.05, th_0_1 = TH_expr_genes_0.1)

expr_genes_MT <- prepareExprGenesOS(species = "MT", th_0 = MT_expr_genes_0,
	th_0_01 = MT_expr_genes_0.01, th_0_05 = MT_expr_genes_0.05, th_0_1 = MT_expr_genes_0.1)

expr_genes_BD <- prepareExprGenesOS(species = "BD", th_0 = BD_expr_genes_0,
	th_0_01 = BD_expr_genes_0.01, th_0_05 = BD_expr_genes_0.05, th_0_1 = BD_expr_genes_0.1)


# Combine expr_genes data from all species
expr_genes_OS <- rbind(expr_genes_AL, expr_genes_CR, expr_genes_ES, expr_genes_TH, expr_genes_MT, 
	expr_genes_BD)




# Plot number of deduplicated reads for each species
plotExprGenesOS <- function(data) {

	fname <- sprintf('%s.jpg', paste(deparse(substitute(data)), sep="_"))

	level_order <- c("root", "hypocotyl", "leaf", "apex veg", "apex inf", 
		"flower", "carpel", "stamen", "pollen")


	p <- ggplot(data, aes(x = factor(sample_names), y = expressed, color = Threshold, group = Threshold)) + 
  	geom_line(aes(x = factor(sample_names, level= level_order)), size=0.55) + 
	scale_y_continuous(expand = c(0.15, 0.15), breaks= pretty_breaks(), 
		 	labels = function(l) { 
		 		ifelse(l==0, paste0(round(l/1e3,1)),paste0(round(l/1e3,1),"K"))
		 	})

	q <- p + facet_wrap(class ~ species, scales='free_y', ncol = 6) + 
	theme_bw() + xlab("Samples") + ylab("Number of Genes") + 
	scale_color_manual(values=c("gray45","#ea6965","#967cee","#dca207")) + 
		guides(colour = guide_legend(nrow = 1)) + 
  		theme(
  		strip.text = element_blank(), 
  		strip.background = element_blank(),
  		plot.margin = unit(c(15, 0, 0, 2), "points"),
  		axis.ticks.length = unit(.075, "cm"),
  		axis.ticks = element_line(colour = "gray15", size = 0.2), 
  		panel.grid.major = element_line(size = 0.4), 
  		panel.grid.minor = element_line(size = 0.25),  
  		axis.title.x = element_text(colour = "black", size=6.75, 
  			margin = margin(t = 5, r = 0, b = -14, l = 0)), 
  		axis.title.y = element_text(colour = "black", size=6.75, 
  			margin = margin(t = 0, r = 8, b = 0, l = 0)), 
  		axis.text.x = element_text(colour = "black", size=5.5, angle=90, 
  			margin = margin(t = 0.5, r = 0, b = 0, l = 0), hjust = 1, vjust = 0.5), 
  		axis.text.y = element_text(colour = "black", size=5.0, margin = margin(t = 0, r = 0.4, b = 0, l = -3.25)),  
		legend.position = "bottom",
		legend.title = element_text(colour = "black", size=6.2, face ="bold"),
		legend.text = element_text(size=6.1), 
		legend.key.size = unit(0.5, "cm"),
		legend.key.height = unit(0.4, "cm"),
		legend.background = element_rect(fill = NA),
		legend.key = element_rect(fill = NA),
  		panel.border = element_rect(colour = "grey70", fill=NA, size=0.7))


  	pg <- ggplot(data, aes(x = factor(sample_names), y = expressed, color = Threshold, group = Threshold)) + 
  	geom_line(aes(x = factor(sample_names, level= level_order)), size=0.55) + 
	scale_y_continuous(expand = c(0.15, 0.15), breaks= pretty_breaks(), 
		 	labels = function(l) { 
		 		ifelse(l==0, paste0(round(l/1e3,1)),paste0(round(l/1e3,1),"K"))
		 	})

	qg <- pg + facet_grid(class ~ species, scales='free') + 
	theme_bw() + xlab("Samples") + ylab("Number of Genes") + 
	scale_color_manual(values=c("gray45","#ea6965","#967cee","#dca207")) + 
		guides(colour = guide_legend(nrow = 1)) + 
  		theme(
  		strip.text.x = element_text(margin = margin(0.105, 0, 0.105, 0, "cm"), size = 5.75), 
        strip.text.y = element_text(margin = margin(0, 0.105, 0, 0.105, "cm"), size = 5.75),
        strip.background = element_rect(colour="grey70", size=0.7), 
  		plot.margin = unit(c(15, 0, 0, 2), "points"),
  		axis.ticks.length = unit(.075, "cm"),
  		panel.grid.major = element_line(size = 0.4), 
  		panel.grid.minor = element_line(size = 0.25),  
  		axis.ticks = element_line(colour = "gray15", size = 0.2), 
  		axis.title.x = element_text(colour = "black", size=6.75, 
  			margin = margin(t = 5, r = 0, b = -14, l = 0)), 
  		axis.title.y = element_text(colour = "black", size=6.75, 
  			margin = margin(t = 0, r = 8, b = 0, l = 0)), 
  		axis.text.x = element_text(colour = "black", size=5.5, angle=90, 
  			margin = margin(t = 0.5, r = 0, b = 0, l = 0), hjust = 1, vjust = 0.5), 
  		axis.text.y = element_text(colour = "black", size=5.0, margin = margin(t = 0, r = 0.4, b = 0, l = -3.25)),  
		legend.position = "bottom",
		legend.title = element_text(colour = "black", size=6.2, face ="bold"),
		legend.text = element_text(size=6.1), 
		legend.key.size = unit(0.5, "cm"),
		legend.key.height = unit(0.4, "cm"),
		legend.background = element_rect(fill = NA),
		legend.key = element_rect(fill = NA),
  		panel.border = element_rect(colour = "grey70", fill=NA, size=0.7))


	gt1 = ggplot_gtable(ggplot_build(q))
	gt2 = ggplot_gtable(ggplot_build(qg))

	gt1 <- gtable_add_rows(gt1, heights = unit(0.151, 'cm'), pos = 2)
	gt1 <- gtable_add_grob(gt1, grobs = gt2$grobs[grep('strip-t', gt2$layout$name)], t = 2, l = gt1$layout[grep('strip-t.+1$', gt1$layout$name),]$l)

    gt.side1 = gtable_filter(gt2, 'strip-r-1')
    gt.side2 = gtable_filter(gt2, 'strip-r-2')
    gt.side3 = gtable_filter(gt2, 'strip-r-3')
    gt.side4 = gtable_filter(gt2, 'strip-r-4')

    gt1 = gtable_add_cols(gt1, widths = unit(0.3025, 'cm'), pos = -1)
    gt1 = gtable_add_grob(gt1, zeroGrob(), t = 1, l = ncol(gt1), b=nrow(gt1))

    panel_id <- gt1$layout[grep('panel-.+1$', gt1$layout$name),]
    gt1 = gtable_add_grob(gt1, gt.side1, t = panel_id$t[1], l = ncol(gt1))
    gt1 = gtable_add_grob(gt1, gt.side2, t = panel_id$t[2], l = ncol(gt1))
    gt1 = gtable_add_grob(gt1, gt.side3, t = panel_id$t[3], l = ncol(gt1))
    gt1 = gtable_add_grob(gt1, gt.side4, t = panel_id$t[4], l = ncol(gt1))


  	png("NUL")

	ggsave(file = file.path(out_dir, "output", "plots", fname), plot = gt1,
		scale = 1, width = 6.2, height = 4.5, units = c("in"), 
		dpi = 800, limitsize = FALSE)
}


plotExprGenesOS(data = expr_genes_OS)



# Generate blank plot and add annotations (workaround to add number of expressed genes 
# per species and transcript class for non-ATH species, since this was not possible to add 
# directly to the plot generated by facet_wrap/grid followed by gtable function...)
# This plot has same dimensions as "plotExprGenesOS", but features complete transparency
# with added text annotations

df_blank <- data.frame()

AT_total_cod <- paste("Total:", expr_genes_OS[10,2], sep=" ")
CR_total_cod <- paste("Total:", expr_genes_OS[154,2], sep=" ")
ES_total_cod <- paste("Total:", expr_genes_OS[298,2], sep=" ")
TH_total_cod <- paste("Total:", expr_genes_OS[442,2], sep=" ")
MT_total_cod <- paste("Total:", expr_genes_OS[586,2], sep=" ")
BD_total_cod <- paste("Total:", expr_genes_OS[730,2], sep=" ")


dat_text <- data.frame(
    label = c(AT_total_cod, CR_total_cod, ES_total_cod, TH_total_cod, MT_total_cod, BD_total_cod),
    x     = c(0.314, 1.376, 2.438, 3.5, 4.562, 5.625),
    y     = c(32.4, 32.4, 32.4, 32.4, 32.4, 32.4)
    )

p <- ggplot(df_blank) + geom_point() + xlim(0, 6) + ylim(0, 40) + theme_void() + 
theme(
    panel.background = element_rect(fill = "transparent"), # bg of the panel
    plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
    panel.grid.major = element_blank(), # get rid of major grid
    panel.grid.minor = element_blank(), # get rid of minor grid
    legend.background = element_rect(fill = "transparent"), # get rid of legend bg
    legend.box.background = element_rect(fill = "transparent") # get rid of legend panel bg
)

q <- p + geom_text(
    data = dat_text,
    mapping = aes(x = x, y = y, label = label), 
    size=2
    )

png("NUL")

ggsave(file = file.path(out_dir, "output", "plots", "non-ATH_annotation_layer.png"), plot = q,
	scale = 1, width = 6.2, height = 4.5, units = c("in"), dpi = 800, limitsize = FALSE, 
	bg = "transparent")




#-------------------------- Do hclust dendrograms for all species ---------------------------


# General settings

# Create shorter sample descriptions
AT_names <- rep(c("root_root tip 5d", "root mat.zone 5d", "root whole rt 5d", "root whole rt 7d",
 "root whole rt 14d", "root whole rt 21d", "hypocotyl 10d", "internode 3rd 24d", "internode 2nd 24d", 
 "internode 1st 28d", "cotyledons 7d", "leaf 1+2 7d", "leaf 1+2 10d", "leaf petiole 10d", 
 "leaf tip 10d", "leaf 5+6 17d", "leaf 9+10 27d", "leaf senesc 35d", "cauline leaf 24d", 
 "apex veg 7d", "apex veg 10d", "apex veg 14d", "apex inf 21d", "apex inf clv1 21d", "apex inf 28d", 
 "flower st9", "flower st10/11", "flower st12", "flower st15", "sepals st12", "sepals st15", 
 "petals st12", "petals st15", "stamen st12", "stamen st15", "pollen mature", "carpel early st12", 
 "carpel late st12", "fruit st15", "fruit st16", "fruit st17a", "seeds st16", "seeds st17a", 
 "seeds st18"), each=3)

replicate_tag <- rep(c(" 1"," 2"," 3"), times=44)

AT_names <- paste0(AT_names,replicate_tag)

names(ATH_th_genes_repl_tpm_0.05)[4:ncol(ATH_th_genes_repl_tpm_0.05)] <- AT_names


# Define colors based on sample name
label_col <- c(roo="#52428c", hyp="#808dc2", int="#0c703d", lea="#00994f", cot="#00994f", 
	cau="#00994f", ape="#f4dc28", flo="#de6daf", sep="#84cd6a", pet="#ead1c7", 
	sta="#f23d29", pol="#a63126", car="#f2a529", fru="#b54185", see="#e9a3b3")

# apex color for Fig1
# ape="#f4dc28"

# color setting for comparative heatamp and PCA: 
# split apex samles => apex veg="#95b73a", apex inf="#fad819" 
# hypocotyl slightly lighter hyp="#8591c7" root lighter roo="#5850a3"


# Generate hclust dendrogram using relative expression data
makeDendrogram <- function(x, coefficient = c("pearson", "spearman"), 
	biotype = c("protein_coding", "antisense" , "lnc_intergenic"), label_col, d_leaf, 
	cby_shift, cby_scale) {

	# Show error message if no scaling is chosen
	if (missing(coefficient))

		stop(
			"Please choose one of the following coefficients: 
			'pearson', 'spearman'",
			call. = TRUE
			)

	if (missing(biotype))

		stop(
			"Please choose one of the following biotypes: 
			'protein_coding', 'antisense', 'lnc_intergenic'",
			call. = TRUE
			)

	# Set filename
    dfname <- deparse(substitute(x))
    coefficient_tag <- match.arg(coefficient)
    fname <- sprintf('%s_dend.png', paste(dfname, biotype, coefficient_tag, sep="_"))
    species <- substr(dfname, start = 1, stop = 2)

	if (is.element("protein_coding", biotype)) {

		x <- subset(x, biotype=="protein_coding")

	} else if (is.element("antisense", biotype)) {

		x <- x[x$biotype %like% "antisense", ]

	} else if (is.element("lnc_intergenic", biotype)) {

		x <- subset(x, biotype=="lnc_intergenic")
	}

    df_t <- t(x[, 4:ncol(x)]) # transposes data frame so rows become columns and vice versa
    df_t[is.na(df_t)] <- 0 # replaces NAs by 0

    # Build distance matrix
    if (is.element(coefficient, c("pearson"))) {
        df_t_dist.mat <- get_dist(df_t, stand = FALSE, method = "pearson")

    } else if (is.element(coefficient, c("spearman"))) {
      df_t_dist.mat <- get_dist(df_t, stand = FALSE, method = "spearman")
    } 

    df_clust.res <- hclust(df_t_dist.mat, method = "average") # agglomerate clustering using average linkage
  
    df_dend <- dendrapply(as.dendrogram(df_clust.res), function(n){
    
    if (is.leaf(n)){
      dend_col <- label_col[substr(attr(n,"label"),1,3)]
      attr(n, "nodePar") <- list(pch = NA, lab.col = dend_col) # to define label color
      attr(n, "edgePar") <- list(col = dend_col) # to color branch
      }
    return(n)
    })

    # make branch colors extend to last common node
    brc_col <- label_col[substr(colnames(x[, 4:ncol(x)]),1,3)]
    brc_col <- brc_col[order.dendrogram(df_dend)]
    brc_col <- factor(brc_col, unique(brc_col))

    png(height = 725, width = 2600, pointsize = 10, 
    	file = file.path(out_dir, "output", "plots", fname))
    par(mar = c(10, 3.5, 0.5, 0), lwd = 6.5, cex = 2.025, cex.axis = 1)
    df_dend = color_branches(df_dend, clusters = as.numeric(brc_col), col = levels(brc_col))

    if ((species == "AT") && (biotype == "protein_coding") && (coefficient == "pearson")) { 
    	df_dend <- rotate(df_dend,c(1:15,19:33,16:18,34:39,55:75,46:54,40:45,76:78,82:87,79:81,88:132))
    }

    else if ((species == "AT") && (biotype == "antisense") && (coefficient == "pearson")) { 
    	df_dend <- rotate(df_dend,c(4:6,1:3,7:33,37:45,34:36,76:87,112:117,109:111,103:108,88:102,130:132,124:129,121:123,118:120,46:75))
    }

    else if ((species == "AT") && (biotype == "lnc_intergenic") && (coefficient == "pearson")) { 
    	df_dend <- rotate(df_dend,c(124:132,1:3,7:9,4:6,64:69,73:84,70:73,13:30,34:39,46:48,40:45,31:33,55:63,49:54,85:123))
    }

    else if ((species == "AT") && (biotype == "protein_coding") && (coefficient == "spearman")) { 
    	df_dend <- rotate(df_dend,c(1:9,13:27,10:12,28:30,34:36,31:33,37:45,49:54,46:48,91:105,112:126,106:111,127:132,55:90))
    }

    else if ((species == "AT") && (biotype == "antisense") && (coefficient == "spearman")) { 
    	df_dend <- rotate(df_dend,c(4:6,1:3,7:12,16:30,13:15,31:36,46:51,40:45,37:39,127:132,109:111,115:117,112:114,103:108,88:102,118:126,85:87,82:84,55:57,58:81,52:54))
    }

    else if ((species == "AT") && (biotype == "lnc_intergenic") && (coefficient == "spearman")) { 
    	df_dend <- rotate(df_dend,c(7:9,4:6,1:3,10:24,28:42,25:27,43:51,58:66,52:57,73:78,67:72,79:90,121:132,91:120))
    }

    # Get color vector for reordered dendrogram
    brc_col <- label_col[substr(colnames(x[, 4:ncol(x)]),1,3)]
    brc_col <- brc_col[order.dendrogram(df_dend)]
    brc_col <- factor(brc_col, unique(brc_col))

    plot(df_dend, dLeaf = d_leaf)
    df_dend = colored_bars(colors = brc_col, dend = df_dend, add=TRUE, sort_by_labels_order=FALSE, 
    	y_shift = cby_shift, y_scale = cby_scale, rowLabels = "")
    dev.off()
}

makeDendrogram(ATH_th_genes_repl_tpm_0.05, coefficient = "pearson", biotype = "protein_coding", 
	label_col = label_col, d_leaf = 0.111, cby_shift = -0.01125, cby_scale=0.094)
makeDendrogram(ATH_th_genes_repl_tpm_0.05, coefficient = "pearson", biotype = "antisense", 
	label_col = label_col, d_leaf = 0.098, cby_shift = -0.01, cby_scale=0.083)
makeDendrogram(ATH_th_genes_repl_tpm_0.05, coefficient = "pearson", biotype = "lnc_intergenic", 
	label_col = label_col,  d_leaf = 0.084, cby_shift = -0.0085, cby_scale=0.0705)
makeDendrogram(ATH_th_genes_repl_tpm_0.05, coefficient = "spearman", biotype = "protein_coding", 
	label_col = label_col, d_leaf = 0.074, cby_shift = -0.0075, cby_scale=0.0625)
makeDendrogram(ATH_th_genes_repl_tpm_0.05, coefficient = "spearman", biotype = "antisense", 
	label_col = label_col, d_leaf = 0.115, cby_shift = -0.01175, cby_scale=0.0965)
makeDendrogram(ATH_th_genes_repl_tpm_0.05, coefficient = "spearman", biotype = "lnc_intergenic", 
	label_col = label_col, d_leaf = 0.1085, cby_shift = -0.01135, cby_scale=0.091)



# -------- non-ATH species -------- 

# Create shorter sample descriptions
brass_names <- rep(c("root, whole root", "hypocotyl", "leaf 1+2", "vegetative apex", "inflorescence apex", 
	"flower stg12", "mature pollen", "carpel stg12", "stamen stg12"), each=3)
TH_names <- rep(c("root, whole root", "hypocotyl", "leaf 1+2", "vegetative apex", "inflorescence apex", 
	"flower stg12 equiv", "mature pollen", "carpel stg12 equiv", "stamen stg12 equiv"), each=3)
MT_names <- rep(c("root, whole root", "hypocotyl", "leaf 2", "vegetative apex", "inflorescence apex", 
	"flower stg8", "mature pollen", "carpel stg8", "stamen stg8"), each=3)
BD_names <- rep(c("root, whole root", "mesocotyl", "leaf 1", "vegetative apex", "spikelet meristem", "floret stg12 equiv", "mature pollen", 
	"carpel stg12 equiv", "stamen stg12 equiv"), each=3)

replicate_tag_comp_samples <- rep(c(" 1"," 2"," 3"), times=9)

brass_names <- paste0(brass_names,replicate_tag_comp_samples)
TH_names <- paste0(TH_names,replicate_tag_comp_samples)
MT_names <- paste0(MT_names,replicate_tag_comp_samples)
BD_names <- paste0(BD_names,replicate_tag_comp_samples)

names(AL_th_genes_repl_tpm_0.05)[4:ncol(AL_th_genes_repl_tpm_0.05)] <- brass_names
names(CR_th_genes_repl_tpm_0.05)[4:ncol(CR_th_genes_repl_tpm_0.05)] <- brass_names
names(ES_th_genes_repl_tpm_0.05)[4:ncol(ES_th_genes_repl_tpm_0.05)] <- brass_names
names(TH_th_genes_repl_tpm_0.05)[4:ncol(TH_th_genes_repl_tpm_0.05)] <- TH_names
names(MT_th_genes_repl_tpm_0.05)[4:ncol(MT_th_genes_repl_tpm_0.05)] <- MT_names
names(BD_th_genes_repl_tpm_0.05)[4:ncol(BD_th_genes_repl_tpm_0.05)] <- BD_names


# Define colors based on sample name
label_col_c <- c(roo="#52428c", hyp="#808dc2", mes="#808dc2", lea="#00994f", veg="#95b73a", 
	inf="#eed410", spi="#eed410", flo="#de6daf", sta="#f23d29", mat="#a63126", car="#f2a529")


# Generate hclust dendrogram for non-ATH species
makeDendrogramC <- function(x, coefficient = c("pearson", "spearman"), 
	biotype = c("protein_coding", "antisense" , "lnc_intergenic"), label_col, d_leaf, 
	cby_shift, cby_scale) {

	# Show error message if no scaling is chosen
	if (missing(coefficient))

		stop(
			"Please choose one of the following coefficients: 
			'pearson', 'spearman'",
			call. = TRUE
			)

	if (missing(biotype))

		stop(
			"Please choose one of the following biotypes: 
			'protein_coding', 'antisense', 'lnc_intergenic'",
			call. = TRUE
			)

	# Set filename
    dfname <- deparse(substitute(x))
    coefficient_tag <- match.arg(coefficient)
    fname <- sprintf('%s_dend.png', paste(dfname, biotype, coefficient_tag, sep="_"))
    species <- substr(dfname, start = 1, stop = 2)

	if (is.element("protein_coding", biotype)) {

		x <- subset(x, biotype=="protein_coding")

	} else if (is.element("antisense", biotype)) {

		x <- x[x$biotype %like% "antisense", ]

	} else if (is.element("lnc_intergenic", biotype)) {

		x <- subset(x, biotype=="lnc_intergenic")
	}

    df_t <- t(x[, 4:ncol(x)]) # transposes data frame so rows become columns and vice versa
    df_t[is.na(df_t)] <- 0 # replaces NAs by 0

    # Build distance matrix
    if (is.element(coefficient, c("pearson"))) {
        df_t_dist.mat <- get_dist(df_t, stand = FALSE, method = "pearson")

    } else if (is.element(coefficient, c("spearman"))) {
      df_t_dist.mat <- get_dist(df_t, stand = FALSE, method = "spearman")
    } 

    df_clust.res <- hclust(df_t_dist.mat, method = "average") # agglomerate clustering using average linkage
  
    df_dend <- dendrapply(as.dendrogram(df_clust.res), function(n){
    
    if (is.leaf(n)){
      dend_col <- label_col[substr(attr(n,"label"),1,3)]
      attr(n, "nodePar") <- list(pch = NA, lab.col = dend_col) # to define label color
      attr(n, "edgePar") <- list(col = dend_col) # to color branch
      }
    return(n)
    })

    # make branch colors extend to last common node
    brc_col <- label_col[substr(colnames(x[, 4:ncol(x)]),1,3)]
    brc_col <- brc_col[order.dendrogram(df_dend)]
    brc_col <- factor(brc_col, unique(brc_col))

    png(height = 1580, width = 1400, pointsize = 10, 
    	file = file.path(out_dir, "output", "plots", fname))
    par(mar = c(10, 2.5, 0.5, 0), lwd = 13, cex = 4.8, cex.axis = 0.8)
    df_dend = color_branches(df_dend, clusters = as.numeric(brc_col), col = levels(brc_col))

    if ((species == "AL") && (biotype == "protein_coding")) { 
    	df_dend <- rotate(df_dend,c(1:3,7:9,4:6,22:27,10:21))
    }
    if (((species == "CR")|(species == "ES")) && (biotype == "protein_coding") && (coefficient == "pearson")) { 
    	df_dend <- rotate(df_dend,c(1:3,7:9,4:6,10:15,25:27,22:24,19:21,16:18))
    }
    if (((species == "CR")|(species == "ES")) && (biotype == "protein_coding") && (coefficient == "spearman")) { 
    	df_dend <- rotate(df_dend,c(1:3,7:9,4:6,10:15,25:27,22:24,16:21))
    }
    if ((species == "TH") && (biotype == "protein_coding")) { 
    	df_dend <- rotate(df_dend,c(1:12,25:27,22:24,13:21))
    }
    if ((species == "MT") && (biotype == "protein_coding") && (coefficient == "pearson")) { 
    	df_dend <- rotate(df_dend,c(1:6,10:12,7:9,25:27,22:24,13:15,19:21,16:18))
    }
    if ((species == "MT") && (biotype == "protein_coding") && (coefficient == "spearman")) { 
    	df_dend <- rotate(df_dend,c(1:12,25:27,22:24,13:15,19:21,16:18))
    }
    if ((species == "BD") && (biotype == "protein_coding") && (coefficient == "pearson")) { 
    	df_dend <- rotate(df_dend,c(4:6,1:3,7:12,22:27,13:15,16:21))
    }
    if ((species == "BD") && (biotype == "protein_coding") && (coefficient == "spearman")) { 
    	df_dend <- rotate(df_dend,c(25:27,22:24,1:6,19:21,16:18,7:15))
    }

    # Get color vector for reordered dendrogram
    brc_col <- label_col[substr(colnames(x[, 4:ncol(x)]),1,3)]
    brc_col <- brc_col[order.dendrogram(df_dend)]
    brc_col <- factor(brc_col, unique(brc_col))

    plot(df_dend, dLeaf = d_leaf)
    df_dend = colored_bars(colors = brc_col, dend = df_dend, add=TRUE, sort_by_labels_order=FALSE, 
    	y_shift = cby_shift, y_scale = cby_scale, rowLabels = "")
    dev.off()
}

# Make Pearson cor dendrograms
makeDendrogramC(AL_th_genes_repl_tpm_0.05, coefficient = "pearson", biotype = "protein_coding", 
	label_col = label_col_c, d_leaf = 0.067, cby_shift = -0.01125, cby_scale=0.05)
makeDendrogramC(CR_th_genes_repl_tpm_0.05, coefficient = "pearson", biotype = "protein_coding", 
	label_col = label_col_c, d_leaf = 0.0675, cby_shift = -0.01125, cby_scale=0.0501)
makeDendrogramC(ES_th_genes_repl_tpm_0.05, coefficient = "pearson", biotype = "protein_coding", 
	label_col = label_col_c, d_leaf = 0.0678, cby_shift = -0.01125, cby_scale=0.0507)
makeDendrogramC(TH_th_genes_repl_tpm_0.05, coefficient = "pearson", biotype = "protein_coding", 
	label_col = label_col_c, d_leaf = 0.065, cby_shift = -0.011, cby_scale=0.0484)
makeDendrogramC(MT_th_genes_repl_tpm_0.05, coefficient = "pearson", biotype = "protein_coding", 
	label_col = label_col_c, d_leaf = 0.0638, cby_shift = -0.0098, cby_scale=0.0485)
makeDendrogramC(BD_th_genes_repl_tpm_0.05, coefficient = "pearson", biotype = "protein_coding", 
	label_col = label_col_c, d_leaf = 0.052625, cby_shift = -0.008, cby_scale=0.04)

# Make Spearman cor dendrograms
makeDendrogramC(AL_th_genes_repl_tpm_0.05, coefficient = "spearman", biotype = "protein_coding", 
	label_col = label_col_c, d_leaf = 0.047, cby_shift = -0.008, cby_scale=0.035)
makeDendrogramC(CR_th_genes_repl_tpm_0.05, coefficient = "spearman", biotype = "protein_coding", 
	label_col = label_col_c, d_leaf = 0.0355, cby_shift = -0.00575, cby_scale=0.0268)
makeDendrogramC(ES_th_genes_repl_tpm_0.05, coefficient = "spearman", biotype = "protein_coding", 
	label_col = label_col_c, d_leaf = 0.048, cby_shift = -0.0075, cby_scale=0.036)
makeDendrogramC(TH_th_genes_repl_tpm_0.05, coefficient = "spearman", biotype = "protein_coding", 
	label_col = label_col_c, d_leaf = 0.047, cby_shift = -0.0076, cby_scale=0.035)
makeDendrogramC(MT_th_genes_repl_tpm_0.05, coefficient = "spearman", biotype = "protein_coding", 
	label_col = label_col_c, d_leaf = 0.0415, cby_shift = -0.0068, cby_scale=0.031)
makeDendrogramC(BD_th_genes_repl_tpm_0.05, coefficient = "spearman", biotype = "protein_coding", 
	label_col = label_col_c, d_leaf = 0.0455, cby_shift = -0.0077, cby_scale=0.034)













