# This script loads and analysis the data statistics and expression tables for protein-coding 
# genes, isoforms, lncRNAs and circRNAs and generates the plots for the DevSeq transcriptome  
# single-species expression figures


#------------------- Load packages, set directories and read sample tables ---------------------


# Install and load packages
if (!require(dplyr)) install.packages('dplyr')
library(dplyr)
if (!require(ggplot2)) install.packages('ggplot2')
library(ggplot2)
if (!require(mgcv)) install.packages('mgcv')
library(mgcv)
if (!require(grid)) install.packages('grid')
library(grid)


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
		 annotate("text", x = 2.65, y = Inf, hjust = 0, vjust = 1.5, size=7, label = total_dedupl)

	q <- p + scale_fill_manual(values=c("#b2b2b2", "#d8a900", "#35bceb")) + theme_minimal() + 
	xlab("") + ylab("PE reads") + ggtitle(plot_title) + 
	geom_hline(yintercept=30e6, linetype="dashed", color = "red", size=1) + 
	theme(legend.position = "none", 
		text=element_text(size=23), 
  		axis.ticks.length = unit(.3, "cm"),
  		axis.ticks = element_line(colour = "gray15", size = 0.7), 
  		axis.title.y = element_text(colour = "black", size=22, 
  			margin = margin(t = 0, r = 15, b = 0, l = 0)), 
  		axis.text.x = element_text(colour = "black", size=22, angle=90, 
  			margin = margin(t = 5, r = 0, b = 0, l = 0), hjust = 1, vjust = 0.5),
  		axis.text.y = element_text(colour = "black", margin = margin(t = 0, r = 5, b = 0, l = 0)), 
  		plot.title = element_text(colour = "black", size=24, 
  			margin = margin(t = 16, r = 0, b = 18.5, l = 0), hjust = 0.5), 
  		plot.margin = unit(c(0, 7, 0, 14), "points"))
	

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
		 geom_point(aes(x=3, y=data_outl[1,2]), shape=21, colour="gray35", size=2.5, fill="white", stroke=2) + 
		 geom_point(aes(x=3, y=data_outl[2,2]), shape=21, colour="gray35", size=2.5, fill="white", stroke=2) + 
		 geom_point(aes(x=3, y=data_outl[3,2]), shape=21, colour="gray35", size=2.5, fill="white", stroke=2) + 
		 scale_y_continuous(limits = c(0,lim_y), expand = c(0, 0), 
		 	labels = function(l) { 
		 		ifelse(l==0, paste0(round(l/1e6,1)),paste0(round(l/1e6,1),"M"))
		 	}) + 
		 annotate("rect", xmin=0.25, xmax=3.85, ymin=0, ymax=lim_y, fill="white", alpha=0, 
		 	color="black", size=1.35) + 
		 annotate("text", x = 2.65, y = Inf, hjust = 0, vjust = 1.5, size=7, label = total_dedupl)

	q <- p + scale_fill_manual(values=c("#b2b2b2", "#d8a900", "#35bceb")) + theme_minimal() + 
	xlab("") + ylab("PE reads") + ggtitle(plot_title) + 
	geom_hline(yintercept=30e6, linetype="dashed", color = "red", size=1) + 
	theme(legend.position = "none", 
		text=element_text(size=23), 
  		axis.ticks.length = unit(.3, "cm"),
  		axis.ticks = element_line(colour = "gray15", size = 0.7), 
  		axis.title.y = element_text(colour = "black", size=22, 
  			margin = margin(t = 0, r = 15, b = 0, l = 0)), 
  		axis.text.x = element_text(colour = "black", size=22, angle=90, 
  			margin = margin(t = 5, r = 0, b = 0, l = 0), hjust = 1, vjust = 0.5), 
  		axis.text.y = element_text(colour = "black", margin = margin(t = 0, r = 5, b = 0, l = 0)), 
  		plot.title = element_text(colour = "black", size=24, 
  			margin = margin(t = 16, r = 0, b = 16.25, l = 0), hjust = 0.5), 
  		plot.margin = unit(c(0, 7, 0, 14), "points"))
	

  	ggsave(file = file.path(out_dir, "output", "plots", fname), plot = q,
		scale = 1, width = 5.0, height = 6.95, units = c("in"), 
		dpi = 600, limitsize = FALSE)
}

makePlotStatsOS(data=non_ATH_stats_df, lim_y=226000000, medw = 0.415, plot_title="Other species")
# 10 data point for trimmed raw reads above lim_y



# Make Deduplicated Read Replicates
makeRepl <- function(x) {

	repl_names_ATH <- data.frame(c("root.1","root.2","root.3","hypocotyl.1","hypocotyl.2",
		"hypocotyl.3","leaf.1","leaf.2","leaf.3","apex_veg.1","apex_veg.2","apex_veg.3",
		"apex_infl.1","apex_infl.2","apex_infl.3","flower.1","flower.2","flower.3","stamen.1",
		"stamen.2","stamen.3","pollen.1","pollen.2","pollen.3","carpel.1","carpel.2","carpel.3"))
	names(repl_names_ATH) <- "Sample_repl"
	repl_names_non_ATH <- c("root.1","root.2","root.3","hypocotyl.1","hypocotyl.2","hypocotyl.3",
		"leaf.1","leaf.2","leaf.3","apex_veg.1","apex_veg.2","apex_veg.3","apex_infl.1",
		"apex_infl.2","apex_infl.3","flower.1","flower.2","flower.3","pollen.1","pollen.2",
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
		"leaf.1","leaf.2","leaf.3","apex_veg.1","apex_veg.2","apex_veg.3","apex_infl.1",
		"apex_infl.2","apex_infl.3","flower.1","flower.2","flower.3","carpel.1","carpel.2",
		"carpel.3","stamen.1","stamen.2","stamen.3","pollen.1","pollen.2","pollen.3")

	x_labels <- rep(c("root","hypocotyl","leaf","apex.veg","apex.infl","flower","stamen",
		"pollen","carpel"),times=7) ## order of labels has to match the order in first species (ATH)

	p <- ggplot(data, aes(x = factor(Sample_repl, level= level_order), y = Deduplicated, color = Species, group = Species)) + 
	geom_line(aes(x = factor(Sample_repl, level= level_order)), size=1.5) + 
  	geom_point(aes(x = factor(Sample_repl, level= level_order)), size=3.25) + 
  	scale_y_continuous(limits = c(0,8e7), expand = c(0, 0), 
		 	labels = function(l) { 
		 		ifelse(l==0, paste0(round(l/1e6,1)),paste0(round(l/1e6,1),"M"))
		 	}) + 
  	scale_x_discrete(labels=x_labels, breaks=data$Sample_repl[seq(1, length(data$Sample_repl), by = 3)]) + 
  	annotate("rect", xmin=0.25, xmax=27.85, ymin=0, ymax=8e7, fill="white", alpha=0, 
		 	color="black", size=0.7) + 
	annotation_custom(segmentsGrob(gp=gpar(col="black", lwd=2)), xmin=1, xmax=3,ymin=-1950000, ymax=-1950000) + 
	annotation_custom(segmentsGrob(gp=gpar(col="black", lwd=2)), xmin=1, xmax=1,ymin=-1950000, ymax=-1000000) + 
	annotation_custom(segmentsGrob(gp=gpar(col="black", lwd=2)), xmin=3, xmax=3,ymin=-1950000, ymax=-1000000) + 
	annotation_custom(segmentsGrob(gp=gpar(col="black", lwd=2)), xmin=4, xmax=6,ymin=-1950000, ymax=-1950000) +
	annotation_custom(segmentsGrob(gp=gpar(col="black", lwd=2)), xmin=4, xmax=4,ymin=-1950000, ymax=-1000000) + 
	annotation_custom(segmentsGrob(gp=gpar(col="black", lwd=2)), xmin=6, xmax=6,ymin=-1950000, ymax=-1000000) + 
	annotation_custom(segmentsGrob(gp=gpar(col="black", lwd=2)), xmin=7, xmax=9,ymin=-1950000, ymax=-1950000) +
	annotation_custom(segmentsGrob(gp=gpar(col="black", lwd=2)), xmin=7, xmax=7,ymin=-1950000, ymax=-1000000) + 
	annotation_custom(segmentsGrob(gp=gpar(col="black", lwd=2)), xmin=9, xmax=9,ymin=-1950000, ymax=-1000000) + 
	annotation_custom(segmentsGrob(gp=gpar(col="black", lwd=2)), xmin=10, xmax=12,ymin=-1950000, ymax=-1950000) +
	annotation_custom(segmentsGrob(gp=gpar(col="black", lwd=2)), xmin=10, xmax=10,ymin=-1950000, ymax=-1000000) + 
	annotation_custom(segmentsGrob(gp=gpar(col="black", lwd=2)), xmin=12, xmax=12,ymin=-1950000, ymax=-1000000) + 
	annotation_custom(segmentsGrob(gp=gpar(col="black", lwd=2)), xmin=13, xmax=15,ymin=-1950000, ymax=-1950000) +
	annotation_custom(segmentsGrob(gp=gpar(col="black", lwd=2)), xmin=13, xmax=13,ymin=-1950000, ymax=-1000000) + 
	annotation_custom(segmentsGrob(gp=gpar(col="black", lwd=2)), xmin=15, xmax=15,ymin=-1950000, ymax=-1000000) + 
	annotation_custom(segmentsGrob(gp=gpar(col="black", lwd=2)), xmin=16, xmax=18,ymin=-1950000, ymax=-1950000) +
	annotation_custom(segmentsGrob(gp=gpar(col="black", lwd=2)), xmin=16, xmax=16,ymin=-1950000, ymax=-1000000) + 
	annotation_custom(segmentsGrob(gp=gpar(col="black", lwd=2)), xmin=18, xmax=18,ymin=-1950000, ymax=-1000000) + 
	annotation_custom(segmentsGrob(gp=gpar(col="black", lwd=2)), xmin=19, xmax=21,ymin=-1950000, ymax=-1950000) +
	annotation_custom(segmentsGrob(gp=gpar(col="black", lwd=2)), xmin=19, xmax=19,ymin=-1950000, ymax=-1000000) + 
	annotation_custom(segmentsGrob(gp=gpar(col="black", lwd=2)), xmin=21, xmax=21,ymin=-1950000, ymax=-1000000) + 
	annotation_custom(segmentsGrob(gp=gpar(col="black", lwd=2)), xmin=22, xmax=24,ymin=-1950000, ymax=-1950000) +
	annotation_custom(segmentsGrob(gp=gpar(col="black", lwd=2)), xmin=22, xmax=22,ymin=-1950000, ymax=-1000000) + 
	annotation_custom(segmentsGrob(gp=gpar(col="black", lwd=2)), xmin=24, xmax=24,ymin=-1950000, ymax=-1000000) + 
	annotation_custom(segmentsGrob(gp=gpar(col="black", lwd=2)), xmin=25, xmax=27,ymin=-1950000, ymax=-1950000) +
	annotation_custom(segmentsGrob(gp=gpar(col="black", lwd=2)), xmin=25, xmax=25,ymin=-1950000, ymax=-1000000) + 
	annotation_custom(segmentsGrob(gp=gpar(col="black", lwd=2)), xmin=27, xmax=27,ymin=-1950000, ymax=-1000000) + 
  	labs(color="Species")

	q <- p + ggtitle(plot_title) + theme_bw() + xlab("") + ylab("PE dedupl. reads") + 
	scale_color_manual(values=c("#ea6965","#dca207","#46ae12","#1fac7b","#36a5d8","#967cee","#e057c3")) + 
		guides(colour = guide_legend(nrow = 1)) + 
  		theme(text=element_text(size=23), 
  		axis.ticks.length = unit(.3, "cm"),
  		axis.ticks = element_line(colour = "gray15", size = 0.7), 
  		axis.title.y = element_text(colour = "black", size=22, 
  			margin = margin(t = 0, r = 15, b = 0, l = 0)), 
  		axis.text.x = element_text(colour = "black", size=22, angle=90, 
  			margin = margin(t = 5.75, r = 0, b = 0, l = 0), hjust = 1, vjust = 1.5), 
  		axis.ticks.x = element_blank(),
  		axis.text.y = element_text(colour = "black", margin = margin(t = 0, r = 5, b = 0, l = 0)), 
  		plot.title = element_text(colour = "black", size=24, 
  			margin = margin(t = 16, r = 0, b = 16.5, l = 0), hjust = 0.5), 
  		plot.margin = unit(c(0, 2, 0, 1), "points"),
		legend.position = c(0.385,0.125),
		legend.title = element_text(colour = "black", size=20.5, face ="bold"),
		legend.text=element_text(size=20.5), 
		legend.spacing.x = unit(0.25, 'cm'),
		legend.key.size = unit(0.775, "cm"),
  		panel.border = element_rect(colour = "black", fill=NA, size=0.5))

  	png("NUL")
	r <- ggplotGrob(q)
	r$layout$clip[r$layout$name=="panel"] <- "off"

	ggsave(file = file.path(out_dir, "output", "plots", fname), plot = r,
		scale = 1, width = 9.1, height = 7.09, units = c("in"), 
		dpi = 600, limitsize = FALSE)
}

plotDedupReads(data=comp_stats_df, plot_title="Comparative samples")




#------------------------------ Plotting expressed genes data -------------------------------


# Prepare data for ggplot
prepareExprGenes <- function(biotype = c("coding", "NAT", "lincRNA"), th_0, th_0_01, th_0_05, 
	th_0_1) {

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

	sample_names <- c("root tip", "root m.zone", "whole root 5d", "whole root 7d", "whole root 14", 
		"whole root 21", "hypocotyl 10d", "3rd internode", "2nd internode", "1st internode", "cotyledons", 
		"leaf 1.2 7d", "leaf 1.2 10d", "leaf petiole", "leaf tip 10d", "leaf 5.6 17d", "leaf 9.10 27d", "leaf sen.35d", 
		"cauline leaf", "apex veg 7d", "apex veg 10d", "apex veg 14d", "apex inf 21d", "apex inf clv1", 
		"apex inf 28d", "flower st9", "flower 10.11", "flower st12", "flower st15", 
		"sepals st12", "sepals st15", "petals st12", "petals st15", "stamens st12", 
		"stamens st15", "mature pollen", "carpels st12.e", "carpels st12.l", "carpels st15", 
		"fruit st16", "fruit st17a", "seeds st16", "seeds st17a", "seeds st18")
	sample_names <- as.data.frame(rep(sample_names, times=4))

	threshold <- c("0", "0.01", "0.05", "0.1")
	threshold <- as.data.frame(rep(threshold, each=44))

	expr_df <- cbind(detailed_sample_name,sample_names, threshold, data)
	colnames(expr_df) <- df_names

	return(expr_df)
}


expr_coding_genes_ATH <- prepareExprGenes(biotype = "coding", th_0 = ATH_expr_genes_0,
	th_0_01 = ATH_expr_genes_0.01, th_0_05 = ATH_expr_genes_0.05, th_0_1 = ATH_expr_genes_0.1)
expr_NATs_ATH <- prepareExprGenes(biotype = "NAT", th_0 = ATH_expr_genes_0,
	th_0_01 = ATH_expr_genes_0.01, th_0_05 = ATH_expr_genes_0.05, th_0_1 = ATH_expr_genes_0.1)
expr_lincRNAs_ATH <- prepareExprGenes(biotype = "lincRNA", th_0 = ATH_expr_genes_0,
	th_0_01 = ATH_expr_genes_0.01, th_0_05 = ATH_expr_genes_0.05, th_0_1 = ATH_expr_genes_0.1)






# Plot number of deduplicated reads for each species
plotExprGenes <- function(data, plot_title, biotype = c("coding","NAT","linc"), texpr) {

	fname <- sprintf('%s.jpg', paste(deparse(substitute(data)), sep="_"))

	total_expr <- paste("total:", texpr, "at 0.05" , sep=" ")

	if (is.element("coding", biotype)) {
		breaksY <- c(1.5e4,2e4,2.5e4)
		pltymin <- 1.375e4
		pltymax <- 2.77e4
		xtepos <- 30.25
		y_margin <- margin(t = 0, r = 15, b = 0, l = 0)

	} else if (is.element("NAT", biotype)) {
		breaksY <- c(1e3,2e3,3e3)
		pltymin <- 5.25e2
		pltymax <- 3.7e3
		xtepos <- 31.25
		y_margin <- margin(t = 0, r = 15, b = 0, l = 10)

	} else if (is.element("linc", biotype)) {
		breaksY <- c(5e2,1e3,1.5e3)
		pltymin <- 2.0e2
		pltymax <- 1.575e3
		xtepos <- 31.05
		y_margin <- margin(t = 0, r = 10, b = 0, l = 0)
	}

	level_order <- c("root tip", "root m.zone", "whole root 5d", "whole root 7d", "whole root 14", 
		"whole root 21", "hypocotyl 10d", "3rd internode", "2nd internode", "1st internode", "cotyledons", 
		"leaf 1.2 7d", "leaf 1.2 10d", "leaf petiole", "leaf tip 10d", "leaf 5.6 17d", "leaf 9.10 27d", "leaf sen.35d", 
		"cauline leaf", "apex veg 7d", "apex veg 10d", "apex veg 14d", "apex inf 21d", "apex inf clv1", 
		"apex inf 28d", "flower st9", "flower 10.11", "flower st12", "flower st15", 
		"sepals st12", "sepals st15", "petals st12", "petals st15", "stamens st12", 
		"stamens st15", "mature pollen", "carpels st12.e", "carpels st12.l", "carpels st15", 
		"fruit st16", "fruit st17a", "seeds st16", "seeds st17a", "seeds st18")

	p <- ggplot(data, aes(x = factor(Sample, level= level_order), y = Expressed, color = Threshold, group = Threshold)) + 

	geom_line(aes(x = factor(Sample, level= level_order)), size=1.55) + 
	scale_y_continuous(limits = c(pltymin,pltymax), breaks = breaksY, expand = c(0, 0), 
		 	labels = function(l) { 
		 		ifelse(l==0, paste0(round(l/1e3,1)),paste0(round(l/1e3,1),"K"))
		 	}) +
  	annotate("rect", xmin=0.25, xmax=44.75, ymin=pltymin, ymax=pltymax, fill="white", alpha=0, 
		 	color="black", size=0.7) + 
  	annotate("rect", xmin=0.25, xmax=6.5, ymin=pltymin, ymax=pltymax, fill="#747474", alpha=0.175) + 
  	annotate("rect", xmin=10.5, xmax=19.5, ymin=pltymin, ymax=pltymax, fill="#0fc941", alpha=0.175) + 
  	annotate("rect", xmin=25.5, xmax=29.5, ymin=pltymin, ymax=pltymax, fill="#747474", alpha=0.175) + 
  	annotate("rect", xmin=38.5, xmax=44.75, ymin=pltymin, ymax=pltymax, fill="#db1010", alpha=0.175) +
  	geom_line(aes(x = factor(Sample, level= level_order)), size=1.55) + 
  	annotate("text", x = xtepos, y = Inf, hjust = 0, vjust = 21.51, size=7.01, label = total_expr) + 
  	annotate("text", x = 1.675, y = Inf, hjust = 0, vjust = 19.75, size=7.01, label = "Threshold", fontface = 2) + 
  	annotate("text", x = 2, y = Inf, hjust = 0, vjust = 2.4, size=7.25, label = "root") + 
  	annotate("text", x = 6.75, y = Inf, hjust = 0, vjust = 2.4, size= 7.25, label = "stem") + 
  	annotate("text", x = 13.65, y = Inf, hjust = 0, vjust = 2.4, size= 7.25, label = "leaf") + 
  	annotate("text", x = 20.75, y = Inf, hjust = 0, vjust = 2.4, size= 7.25, label = "apex") + 
  	annotate("text", x = 26.1, y = Inf, hjust = 0, vjust = 2.4, size= 7.25, label = "flow") + 
  	annotate("text", x = 30.05, y = Inf, hjust = 0, vjust = 2.4, size= 7.25, label = "floral organ") + 
  	annotate("text", x = 40.15, y = Inf, hjust = 0, vjust = 2.4, size= 7.25, label = "fruit") + 
  	labs(color="")

	q <- p + ggtitle(plot_title) + theme_bw() + xlab("") + ylab("# expressed genes") + 
	scale_color_manual(values=c("gray45","#ea6965","#967cee","#e5a907")) + 
		guides(colour = guide_legend(nrow = 1)) + 
  		theme(text=element_text(size=23), 
  		panel.grid.major = element_line(colour = "white"), 
  		panel.grid.minor = element_line(colour = "white"),  
  		axis.ticks.length = unit(.3, "cm"),
  		axis.ticks = element_line(colour = "gray15", size = 0.7), 
  		axis.title.y = element_text(colour = "black", size=22, 
  			margin = y_margin), 
  		axis.text.x = element_text(colour = "black", size=13.25, angle=90, 
  			margin = margin(t = 3.5, r = 0, b = 0, l = 0), hjust = 1, vjust = 0.5), 
  		axis.text.y = element_text(colour = "black", margin = margin(t = 0, r = 5, b = 0, l = 0)), 
  		plot.title = element_text(colour = "black", size=24, 
  			margin = margin(t = 16, r = 0, b = 16.5, l = 0), hjust = 0.5), 
  		plot.margin = unit(c(0, 2, 18, 4), "points"),
		legend.position = c(0.225,0.115),
		legend.title = element_text(colour = "black", size=20, face ="bold"),
		legend.text = element_text(size=20), 
		legend.key.size = unit(0.775, "cm"),
		legend.key.height = unit(0.4, "cm"),
		legend.background = element_rect(fill = NA),
		legend.key = element_rect(fill = "white"),
  		panel.border = element_rect(colour = "black", fill=NA, size=0.5))

  	png("NUL")
	r <- ggplotGrob(q)
	r$layout$clip[r$layout$name=="panel"] <- "off"

	ggsave(file = file.path(out_dir, "output", "plots", fname), plot = r,
		scale = 1, width = 9.1, height = 7.09, units = c("in"), 
		dpi = 600, limitsize = FALSE)
}


plotExprGenes(data=expr_coding_genes_ATH, plot_title="Expressed protein-coding genes", biotype = "coding", texpr=ATH_expr_genes_0.05[1,2])
plotExprGenes(data=expr_NATs_ATH, plot_title="Expressed NATs", biotype = "NAT", texpr=ATH_expr_genes_0.05[2,2])
plotExprGenes(data=expr_lincRNAs_ATH, plot_title="Expressed lincRNAs", biotype = "linc", texpr=ATH_expr_genes_0.05[3,2])




























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
	z <- wilcox.test(x[,21], y[,21], paired = FALSE)$p.value
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
	plot_title = c("ATH_all", "ATH", "CR", "ES", "TH", "MT", "BD"), yadj) {

	fname <- sprintf('%s.jpg', paste(deparse(substitute(data)), sep="_"))

	cor03_02 = paste("vs.    ", " p=", p03_02, sep="")
	cor03_05 = paste("vs.    ", " p=", p03_05, sep="")
	cor02_05 = paste("vs.    ", " p=", p02_05, sep="")

	blu = rgb(0, 70, 139, max = 255, alpha = 0)
	gray = rgb(131, 145, 145, max = 255, alpha = 0)
	grn = rgb(94, 200, 100, max = 255, alpha = 0)

	p <- ggplot(data, aes(x=max_expression, group=class, fill=class, colour=class, linetype=class)) +
	geom_density(adjust=1.35, size=1.6) +
	scale_x_continuous(trans='log10', labels = prettyNum, limits = c(0.01,100), expand = c(0, 0)) +
	scale_y_continuous(limits = lim_y, expand = c(0, 0)) + 
	annotation_logticks(sides = 'b') + 
	annotate("rect", xmin=6.0, xmax=77, ymin=1.457, ymax=2.014, color="black", fill="white", size=0.55) + 
	annotate("rect", xmin=c(7.4,15.95,7.4,15.95,7.4,15.95), 
		xmax=c(9.6,20.72,9.6,20.72,9.6,20.72), 
		ymin=c(1.7837*yadj,1.7837*yadj,1.6641*yadj,1.6641*yadj,1.5445*yadj,1.5445*yadj), 
		ymax=c(1.792*yadj,1.792*yadj,1.6724*yadj,1.6724*yadj,1.5528*yadj,1.5528*yadj), 
		color=c("#49b43c","#839191","#49b43c","#00468b","#839191","#00468b"), 
		size=1.2, fill=c(grn,gray,grn,blu,gray,blu)) + 
	annotate("text", x = 7.1, y = Inf, hjust = 0, vjust = 2.98, size=5.5, label = "Wilcoxon test", fontface=2) +
	annotate("text", x = 10.6, y = Inf, hjust = 0, vjust = 4.805, size=5.5, label = cor03_02) +
	annotate("text", x = 10.6, y = Inf, hjust = 0, vjust = 6.435, size=5.5, label = cor03_05) + 
	annotate("text", x = 10.6, y = Inf, hjust = 0, vjust = 8.060, size=5.5, label = cor02_05) + 
	annotate("rect", xmin=6.0, xmax=77, ymin=0.817, ymax=1.371, color="black", fill="white", size=0.55) + 
	annotate("rect", xmin=c(7.4,7.4,7.4), 
		xmax=c(9.6,9.6,9.6), 
		ymin=c(1.144*yadj,1.0247*yadj,0.9055*yadj), 
		ymax=c(1.1523*yadj,1.033*yadj,0.9138*yadj), 
		color=c("#49b43c","#839191","#00468b"), 
		size=1.2, fill=c(grn,gray,blu)) + 
	annotate("text", x = 7.1, y = Inf, hjust = 0, vjust = 11.735, size=5.5, label = "Pearson's r", fontface=2) + 
	annotate("text", x = 10.6, y = Inf, hjust = 0, vjust = 13.505, size=5.5, label = "r < -0.3") + 
	annotate("text", x = 10.6, y = Inf, hjust = 0, vjust = 15.13, size=5.5, label = "-0.2 < r < 0.2") + 
	annotate("text", x = 10.6, y = Inf, hjust = 0, vjust = 16.75, size=5.5, label = "r > 0.5") 
	q <- p + ggtitle(plot_title) + theme_bw() + scale_fill_manual(values = c(blu, gray, grn)) +
		scale_color_manual(values = c("#00468b", "#839191", "#52b540")) + xlab("Expression ratio (nc:cd SAS)") + ylab("Density") + 
		scale_linetype_manual(values = c("solid","solid","solid")) + 
		theme(text=element_text(size=16), 
  		axis.ticks.length = unit(.3, "cm"),
  		plot.margin = unit(c(5.5, 13.5, 20.25, 0.5), "points"),
  		axis.text.x = element_text(colour = "black", size=16, angle=0, margin = margin(t = 5, r = 0, b = 0, l = 0)), 
  		axis.text.y = element_text(colour = "black", size=16, angle=0, margin = margin(t = 0, r = 5, b = 0, l = 0)),
  		axis.title.x = element_text(colour = "black", size=18, margin = margin(t = 14.1, r = 0, b = 0, l = 0)),
  		axis.title.y = element_text(colour = "black", size=18, margin = margin(t = 0, r = 12.1, b = 0, l = 0)),
  		plot.title = element_text(colour = "black", size=17.5, margin = margin(t = 17.25, r = 0, b = 8, l = 0), hjust = 0.5),
  		legend.position = "bottom",
  		panel.border = element_rect(colour = "black", fill=NA, size=1.25))

  	ggsave(file = file.path(out_dir, "output", "plots", fname), plot = q,
		scale = 1, width = 8, height = 6.68, units = c("in"),
		dpi = 825, limitsize = FALSE)
}

makeScrPlotExprRatio(data=ATH_all_cd_nc_max_expr_ratio, lim_y=c(0,2.105), p03_02=ATH_all_wilcox_03_02_ratio, p03_05=ATH_all_wilcox_03_05_ratio, p02_05=ATH_all_wilcox_02_05_ratio, plot_title="SAS expression in ATH", yadj=1)

# makeScrPlotExprRatio(data=ATH_comp_cd_nc_max_expr_ratio, lim_y=c(0,1.615), p03_02=ATH_comp_wilcox_03_02_ratio, p03_05=ATH_comp_wilcox_03_05_ratio, p02_05=ATH_comp_wilcox_02_05_ratio, plot_title="ATH_comp", yadj=0.7673)
# makeScrPlotExprRatio(data=AL_cd_nc_max_expr_ratio, lim_y=c(0,1.42), p03_02=AL_wilcox_03_02_ratio, p03_05=AL_wilcox_03_05_ratio, p02_05=AL_wilcox_02_05_ratio, plot_title="AL_", yadj=0.6746)
# makeScrPlotExprRatio(data=CR_cd_nc_max_expr_ratio, lim_y=c(0,1.79), p03_02=CR_wilcox_03_02_ratio, p03_05=CR_wilcox_03_05_ratio, p02_05=CR_wilcox_02_05_ratio, plot_title="CR_", yadj=0.8504)
# makeScrPlotExprRatio(data=ES_cd_nc_max_expr_ratio, lim_y=c(0,1.935), p03_02=ES_wilcox_03_02_ratio, p03_05=ES_wilcox_03_05_ratio, p02_05=ES_wilcox_02_05_ratio, plot_title="ES_", yadj=0.9193)
# makeScrPlotExprRatio(data=TH_cd_nc_max_expr_ratio, lim_y=c(0,1.5), p03_02=TH_wilcox_03_02_ratio, p03_05=TH_wilcox_03_05_ratio, p02_05=TH_wilcox_02_05_ratio, plot_title="TH_", yadj=0.7126)
# makeScrPlotExprRatio(data=MT_cd_nc_max_expr_ratio, lim_y=c(0,1.685), p03_02=MT_wilcox_03_02_ratio, p03_05=MT_wilcox_03_05_ratio, p02_05=MT_wilcox_02_05_ratio, plot_title="MT_", yadj=0.8005)
# makeScrPlotExprRatio(data=BD_cd_nc_max_expr_ratio, lim_y=c(0,1.682), p03_02=BD_wilcox_03_02_ratio, p03_05=BD_wilcox_03_05_ratio, p02_05=BD_wilcox_02_05_ratio, plot_title="BD_", yadj=0.7991)



# Combine all pearson correlation SAS data to create data table for stacked bar chart
species <- c(rep("ATH_all", 3), rep("ATH", 3), rep("AL", 3), rep("CR", 3), rep("ES", 3), rep("TH", 3), rep("MT", 3), rep("BD", 3))
condition <- rep(c("-0.3" , ">-02 <02" , ">05") , 8)
value <- c(
	nrow(subset(ATH_all_cd_nc_max_expr_ratio, class=="-0.3")),
	nrow(subset(ATH_all_cd_nc_max_expr_ratio, class==">-02 <02")),
	nrow(subset(ATH_all_cd_nc_max_expr_ratio, class==">05")),
	nrow(subset(ATH_comp_cd_nc_max_expr_ratio, class=="-0.3")),
	nrow(subset(ATH_comp_cd_nc_max_expr_ratio, class==">-02 <02")),
	nrow(subset(ATH_comp_cd_nc_max_expr_ratio, class==">05")),
	nrow(subset(AL_cd_nc_max_expr_ratio, class=="-0.3")),
	nrow(subset(AL_cd_nc_max_expr_ratio, class==">-02 <02")),
	nrow(subset(AL_cd_nc_max_expr_ratio, class==">05")),
	nrow(subset(CR_cd_nc_max_expr_ratio, class=="-0.3")),
	nrow(subset(CR_cd_nc_max_expr_ratio, class==">-02 <02")),
	nrow(subset(CR_cd_nc_max_expr_ratio, class==">05")),
	nrow(subset(ES_cd_nc_max_expr_ratio, class=="-0.3")),
	nrow(subset(ES_cd_nc_max_expr_ratio, class==">-02 <02")),
	nrow(subset(ES_cd_nc_max_expr_ratio, class==">05")),
	nrow(subset(TH_cd_nc_max_expr_ratio, class=="-0.3")),
	nrow(subset(TH_cd_nc_max_expr_ratio, class==">-02 <02")),
	nrow(subset(TH_cd_nc_max_expr_ratio, class==">05")),
	nrow(subset(MT_cd_nc_max_expr_ratio, class=="-0.3")),
	nrow(subset(MT_cd_nc_max_expr_ratio, class==">-02 <02")),
	nrow(subset(MT_cd_nc_max_expr_ratio, class==">05")),
	nrow(subset(BD_cd_nc_max_expr_ratio, class=="-0.3")),
	nrow(subset(BD_cd_nc_max_expr_ratio, class==">-02 <02")),
	nrow(subset(BD_cd_nc_max_expr_ratio, class==">05"))
	)
SAS_class_abundances <- data.frame(species,condition,value)


# Function for stacked bar chart of pearson correlation SAS groups
makeAbndPlot <- function(data) {

	fname <- 'SAS_class_abundances.jpg'
	data$species <- factor(data$species, levels = unique(data$species))

	blu = "#0c5093"
	gray = "#839191"
	grn = "#5ec864"

	p <- ggplot(data, aes(fill=condition, x=species, y=value)) + 
	scale_y_continuous(limits = c(0,2550), expand = c(0, 0)) + 
    geom_bar(width=0.5, position="stack", stat="identity")
    q <- p + ggtitle("SAS class abundances") +
    theme_bw() + scale_fill_manual(values = c(grn,gray,blu), labels=c("r < -0.3 ","-0.2 < r < 0.2 ","r > 0.5 ")) + 
    xlab("Species") + ylab("Number of cis-NAT pairs") + 
    labs(fill="Pearson's r:", fontface=2) + 
    theme(text=element_text(size=16), 
  		axis.ticks.length = unit(.3, "cm"),
  		plot.margin = unit(c(7.9, 13.5, 60.5, 0.5), "points"),
  		axis.text.x = element_text(colour = "black", size=16, angle=0, margin = margin(t = 5, r = 0, b = 0, l = 0)), 
  		axis.text.y = element_text(colour = "black", size=16, angle=0, margin = margin(t = 0, r = 5, b = 0, l = 0)),
  		axis.title.x = element_text(colour = "black", size=18, margin = margin(t = 14.1, r = 0, b = 0, l = 0)),
  		axis.title.y = element_text(colour = "black", size=18, margin = margin(t = 0, r = 12.1, b = 0, l = 0)),
  		plot.title = element_text(colour = "black", size=17.5, margin = margin(t = 17.25, r = 0, b = 8, l = 0), hjust = 0.5),
  		legend.position=c(0.5,0.9),
  		legend.direction = "horizontal", 
  		legend.box = "horizontal",
  		legend.key.size = unit(1.35, 'lines'),
  		legend.title=element_text(size=17.5), 
    	legend.text=element_text(size=16.5),
  		panel.border = element_rect(colour = "black", fill=NA, size=1.25))
    
  	ggsave(file = file.path(out_dir, "output", "plots", fname), plot = q,
		scale = 1, width = 9.52, height = 6.68, units = c("in"), #smaller img settings: width = 4.848485, height = 5.95
		dpi = 400, limitsize = FALSE)
}

makeAbndPlot(data=SAS_class_abundances)


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
n_ATH_same_strand_PCT <- length(ATH_same_strand_PCT[,16])
n_ATH_SAS_PCT <- length(ATH_SAS_PCT[,16])


png(file=file.path(out_dir, "output", "plots", "cd_cd_SAS_NAT_cd_SAS_pearson_ATH_all_vs_comp.png"), 
	width = 2850, height = 4000, res = 825)
par(mar = c(4.5, 4.5, 4, 2.4))
boxplot(ATH_same_strand_PCT[,16], ATH_SAS_PCT[,16], 
	ATH_coding_SAS_cor_wo_pollen_pearson, ATH_cd_nc_SAS_cor_wo_pollen_0.5_pearson, 
	ylim = c(-1.2, 1.35), 
	names = FALSE, 
	xaxt = 'n', 
	yaxt = 'n', 
	cex.lab = 1.1, 
	las = 2,
	cex.axis = 1.1, #adapt size of axis labels
	ylab = "Pearson ρ", 
	col = c("#e0e0e0", "#c4c4c4", "#a8a8a8", "#d8a900"), 
	boxwex = 0.75, 
	pars = list(outcol = "gray50"), 
	lwd = 1.35, 
	whisklty = 1, 
	at = c(1,2,3,4), 
	notch = FALSE
	)
	title("SAS pairs in ATH", adj = 0.50, line = 1.3, font.main = 1, cex.main = 1.2)
	box(lwd = 1.35)
	axis(side = 2, lwd = 1.35, las = 2)
	text(x= 2.5, y = 1.3, labels= "p<1e-100", col = "black", cex = 1) #ATH_all p-value
	segments(x0=1,x1=3,y0=1.1,y1=1.1, col="gray10", lwd = 1.35)
	segments(x0=2,x1=4,y0=1.2,y1=1.2, col="gray10", lwd = 1.35)
	segments(x0=2,x1=2,y0=1.1,y1=1.2, col="gray10", lwd = 1.35)
	segments(x0=4,x1=4,y0=1.1,y1=1.2, col="gray10", lwd = 1.35)
	text(x= 1, y= -1.175, labels= n_ATH_same_strand_PCT, col= "gray40", cex= 0.97) #ATH_all no.genes
	text(x= 2, y= -1.04, labels= n_ATH_SAS_PCT, col= "gray40", cex= 0.97)
	text(x= 3, y= -1.175, labels= n_ATH_pc_all_wo_pollen, col= "gray40", cex= 0.97) #AL_comp no.genes
	text(x= 4, y= -1.04, labels= n_ATH_nc_all_wo_pollen, col= "gray40", cex= 0.97)
	par(xpd=TRUE)
	legend(-0.35,-1.385,c("cd-cd SSN", "cd-cd OSN"),  
	bty='n', horiz = TRUE, fill = c("#e0e0e0", "#c4c4c4"), cex = 1.1, x.intersp = 0.5, 
	text.width=c(2,1.81))
	legend(-0.35,-1.625,c("cd-cd SAS", "nc-cd SAS"),  
	bty='n', horiz = TRUE, fill = c("#a8a8a8", "#d8a900"), cex = 1.1, x.intersp = 0.5, 
	text.width=c(2,1.81))
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
	ATH_comp_samples_cd_nc_SAS_cor_wo_pollen_5_10_pearson, ATH_comp_samples_cd_nc_SAS_cor_wo_pollen_10_pearson)

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


# Correlation plots of cd-cd SAS / nc-cd SAS pairs ATH_all_samples in DevSeq, Araport and ATGE
DevSeq_pearson <- length(ATH_cd_nc_SAS_cor_wo_pollen_0.5_pearson)
DevSeq_DevSeq <- rbind(
	subset(ATH_cd_nc_SAS_cor_wo_pollen_0.5, biotype_query=="lnc_exonic_antisense" & gene_source_query == "DevSeq"), 
	subset(ATH_cd_nc_SAS_cor_wo_pollen_0.5, biotype_query=="lnc_intronic_antisense" & gene_source_query == "DevSeq"), 
	subset(ATH_cd_nc_SAS_cor_wo_pollen_0.5, biotype_subject=="lnc_exonic_antisense" & gene_source_subject == "DevSeq"), 
	subset(ATH_cd_nc_SAS_cor_wo_pollen_0.5, biotype_subject=="lnc_intronic_antisense" & gene_source_subject == "DevSeq") 
	)
DevSeq_DevSeq_wo_pollen_0.5_pearson <- unlist(select(DevSeq_DevSeq, Pearson))
DevSeq_DevSeq_pearson <- length(DevSeq_DevSeq_wo_pollen_0.5_pearson)
DevSeq_Araport <- rbind(
	subset(ATH_cd_nc_SAS_cor_wo_pollen_0.5, biotype_query=="lnc_exonic_antisense" & gene_source_query == "araport11"), 
	subset(ATH_cd_nc_SAS_cor_wo_pollen_0.5, biotype_query=="lnc_intronic_antisense" & gene_source_query == "araport11"), 
	subset(ATH_cd_nc_SAS_cor_wo_pollen_0.5, biotype_subject=="lnc_exonic_antisense" & gene_source_subject == "araport11"), 
	subset(ATH_cd_nc_SAS_cor_wo_pollen_0.5, biotype_subject=="lnc_intronic_antisense" & gene_source_subject == "araport11") 
	)
DevSeq_Araport_wo_pollen_0.5_pearson <- unlist(select(DevSeq_Araport, Pearson))
DevSeq_Araport_pearson <- length(DevSeq_Araport_wo_pollen_0.5_pearson)
ATGE_pearson <- length(ATH_cd_nc_SAS_cor_wo_pollen_0.5_in_ATGE_pearson)
DevSeq_spearman <- length(ATH_cd_nc_SAS_cor_wo_pollen_0.5_spearman)
ATGE_spearman <- length(ATH_cd_nc_SAS_cor_wo_pollen_0.5_in_ATGE_spearman)


# Pearson boxplot
jpeg(file=file.path(out_dir, "output", "plots", "nccd_SAS_pearson_ATH_DevSeq_Araport.jpeg"), 
	width = 4000, height = 4700, res = 825) 
par(mar = c(6.73, 4.95, 3.14, 0.75))
boxplot(ATH_cd_nc_SAS_cor_wo_pollen_0.5_pearson, DevSeq_DevSeq_wo_pollen_0.5_pearson, 
	DevSeq_Araport_wo_pollen_0.5_pearson, ATH_cd_nc_SAS_cor_wo_pollen_0.5_in_ATGE_pearson, 
	ylim = c(-1.0, 1.0), 
	names = FALSE, 
	xaxt = 'n', 
	yaxt = 'n', 
	cex.lab = 1.1, 
	las = 2,
	cex.axis = 1.1, #adapt size of axis labels
	xlab = "", 
	ylab = "", 
	col = c("#d8a900", "#00bc1f", "#00c094", "#00beda"), 
	boxwex = 0.75, 
	lwd = 1.7, 
	whisklty = 1, 
	at = c(1,2,3,4), 
	pars = list(outcol = "gray50"), 
	notch = FALSE
	)
	title("nc-cd SAS pairs", adj = 0.50, line = 1.3, font.main = 1, cex.main = 1.41)
	title(xlab = "Data set", line = 2.92, cex.lab = 1.34)
	title(ylab = "Pearson ρ", line = 3.5, cex.lab = 1.34)
	rug(x = c(1,2,3,4), ticksize = -0.032, side = 1, lwd = 1.5, col = "black") #x-axis ticks
	box(lwd = 1.7)
	axis(side = 2, lwd = 1.5, las = 2, cex.axis = 1.2, tck = -0.032, mgp=c(3,1.2,0))
	text(x= 1, y= -0.95, labels= DevSeq_pearson, col= "gray40", cex= 1.15) #ATH_all no.genes
	text(x= 2, y= -0.95, labels= DevSeq_DevSeq_pearson, col= "gray40", cex= 1.15)
	text(x= 3, y= -0.95, labels= DevSeq_Araport_pearson, col= "gray40", cex= 1.15)
	text(x= 4, y= -0.95, labels= ATGE_pearson, col= "gray40", cex= 1.15)
	mtext('DevSeq', side = 1, line = 1.12, at = 1, cex = 1.2)
	mtext('DS-DS', side = 1, line = 1.12, at = 2, cex = 1.2)
	mtext('DS-AP', side = 1, line = 1.12, at = 3, cex = 1.2)
	mtext('DS-ATGE', side = 1, line = 1.12, at = 4, cex = 1.2)
	par(xpd=TRUE)
dev.off()



# Compute rsqr cor value
testRsq <- function(x, y) { 
	test <- cor(x, y, use = "complete.obs") ^ 2
  	test <- round(test, digits=2)
  	return(test)
}

# Calculate R squared value for pearson vs spearman data
rsqd_ATH_all_PS <- testRsq(ATH_cd_nc_SAS_cor_wo_pollen_0.5_pearson, ATH_cd_nc_SAS_cor_wo_pollen_0.5_spearman)
rsqd_ATH_comp_PS <- testRsq(ATH_comp_samples_cd_nc_SAS_cor_wo_pollen_0.5_pearson, ATH_comp_samples_cd_nc_SAS_cor_wo_pollen_0.5_spearman)
rsqd_AL_PS <- testRsq(AL_comp_samples_coding_SAS_cor_wo_pollen_pearson, AL_comp_samples_coding_SAS_cor_wo_pollen_spearman)
rsqd_CR_PS <- testRsq(CR_cd_nc_SAS_cor_wo_pollen_0.5_pearson, CR_cd_nc_SAS_cor_wo_pollen_0.5_spearman)
rsqd_ES_PS <- testRsq(ES_cd_nc_SAS_cor_wo_pollen_0.5_pearson, ES_cd_nc_SAS_cor_wo_pollen_0.5_spearman)
rsqd_TH_PS <- testRsq(TH_cd_nc_SAS_cor_wo_pollen_0.5_pearson, TH_cd_nc_SAS_cor_wo_pollen_0.5_spearman)
rsqd_MT_PS <- testRsq(MT_cd_nc_SAS_cor_wo_pollen_0.5_pearson, MT_cd_nc_SAS_cor_wo_pollen_0.5_spearman)
rsqd_BD_PS <- testRsq(BD_cd_nc_SAS_cor_wo_pollen_0.5_pearson, BD_cd_nc_SAS_cor_wo_pollen_0.5_spearman)



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
	"ATH_all", "ATH", "AL", "CR", "ES", "TH", "MT", "BD")) {

	fname <- sprintf('%s.jpg', paste(deparse(substitute(data)), sep="_"))

	ifelse (is.element("ATH_all", plot_title), b_adj<-15.5, b_adj<-17.9)
	
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
  		plot.title = element_text(colour = "black", size=17, margin = margin(t = 11.5, r = 0, b = b_adj, l = 0), hjust = 0.5),
  		legend.position = "bottom",
  		panel.border = element_rect(colour = "black", fill=NA, size=1.2))

	ggsave(file = file.path(out_dir, "output", "plots", fname), plot = q,
		scale = 1, width = 4.848485, height = 5.69697, units = c("in"), 
		dpi = 825, limitsize = FALSE)
}


makeScrPlotPSCor(data=ATH_all_PS, rsqd=rsqd_ATH_all_PS, plot_title="ATH_all")
makeScrPlotPSCor(data=ATH_comp_PS, rsqd=rsqd_ATH_comp_PS, plot_title="ATH")
makeScrPlotPSCor(data=AL_comp_PS, rsqd=rsqd_AL_PS, plot_title="AL")
makeScrPlotPSCor(data=CR_comp_PS, rsqd=rsqd_CR_PS, plot_title="CR")
makeScrPlotPSCor(data=ES_comp_PS, rsqd=rsqd_ES_PS, plot_title="ES")
makeScrPlotPSCor(data=TH_comp_PS, rsqd=rsqd_TH_PS, plot_title="TH")
makeScrPlotPSCor(data=MT_comp_PS, rsqd=rsqd_MT_PS, plot_title="MT")
makeScrPlotPSCor(data=BD_comp_PS, rsqd=rsqd_BD_PS, plot_title="BD")




#--------------- Generate count plots for neighboring protein-coding gene pairs ----------------


# Create intergenic ranges lists
getDistanceRange <- function(x, min_dist, max_dist) {
	data <- subset(x, distance >= min_dist & distance < max_dist)
	return(data)
}

# for cd-cd gene pairs that are on same strand
range_1_50_SSN <- getDistanceRange(ATH_same_strand_PCT, min_dist = 1, max_dist = 50)
range_50_100_SSN <- getDistanceRange(ATH_same_strand_PCT, min_dist = 50, max_dist = 100)
range_100_200_SSN <- getDistanceRange(ATH_same_strand_PCT, min_dist = 100, max_dist = 200)
range_200_500_SSN <- getDistanceRange(ATH_same_strand_PCT, min_dist = 200, max_dist = 500)
range_500_1000_SSN <- getDistanceRange(ATH_same_strand_PCT, min_dist = 500, max_dist = 1000)
range_1000_2000_SSN <- getDistanceRange(ATH_same_strand_PCT, min_dist = 1000, max_dist = 2000)
range_2000_5000_SSN <- getDistanceRange(ATH_same_strand_PCT, min_dist = 2000, max_dist = 5000)
range_5000_SSN <- getDistanceRange(ATH_same_strand_PCT, min_dist = 5000, max_dist = 50000)

# for cd-cd gene pairs that are on opposite strand
range_1_50_OSN <- getDistanceRange(ATH_SAS_PCT, min_dist = 1, max_dist = 50)
range_50_100_OSN <- getDistanceRange(ATH_SAS_PCT, min_dist = 50, max_dist = 100)
range_100_200_OSN <- getDistanceRange(ATH_SAS_PCT, min_dist = 100, max_dist = 200)
range_200_500_OSN <- getDistanceRange(ATH_SAS_PCT, min_dist = 200, max_dist = 500)
range_500_1000_OSN <- getDistanceRange(ATH_SAS_PCT, min_dist = 500, max_dist = 1000)
range_1000_2000_OSN <- getDistanceRange(ATH_SAS_PCT, min_dist = 1000, max_dist = 2000)
range_2000_5000_OSN <- getDistanceRange(ATH_SAS_PCT, min_dist = 2000, max_dist = 5000)
range_5000_OSN <- getDistanceRange(ATH_SAS_PCT, min_dist = 5000, max_dist = 50000)


SSN_list <- list(range_1_50_SSN=range_1_50_SSN,range_50_100_SSN=range_50_100_SSN,
	range_100_200_SSN=range_100_200_SSN,range_200_500_SSN=range_200_500_SSN,
	range_500_1000_SSN=range_500_1000_SSN,range_1000_2000_SSN=range_1000_2000_SSN,
	range_2000_5000_SSN=range_2000_5000_SSN,range_5000_SSN=range_5000_SSN)

OSN_list <- list(range_1_50_OSN=range_1_50_OSN,range_50_100_OSN=range_50_100_OSN,
	range_100_200_OSN=range_100_200_OSN,range_200_500_OSN=range_200_500_OSN,
	range_500_1000_OSN=range_500_1000_OSN,range_1000_2000_OSN=range_1000_2000_OSN,
	range_2000_5000_OSN=range_2000_5000_OSN,range_5000_OSN=range_5000_OSN)


# Get median of Pearson values
med_SSN_pearson <- lapply(SSN_list, function(x) {
	median(x[, 16])
})
names(med_SSN_pearson) <- paste(names(
	med_SSN_pearson),"_median_pearson", sep="")
list2env(med_SSN_pearson, envir = .GlobalEnv)
med_SSN_pearson <- unlist(med_SSN_pearson)
med_SSN_pearson <- as.data.frame(med_SSN_pearson)

med_OSN_pearson <- lapply(OSN_list, function(x) {
	median(x[, 16])
})
names(med_OSN_pearson) <- paste(names(
	med_OSN_pearson),"_median_pearson", sep="")
list2env(med_OSN_pearson, envir = .GlobalEnv)
med_OSN_pearson <- unlist(med_OSN_pearson)
med_OSN_pearson <- as.data.frame(med_OSN_pearson)

dist_range <- c("1-50","50-100","100-200","200-500","500-1K","1K-2K","2K-5K",">5K")
dist_range <- as.data.frame(dist_range)

med_dist_cor <- cbind(dist_range, med_SSN_pearson, med_OSN_pearson)
rownames(med_dist_cor) <- c()


# Get number of genes per intergenic range
n_range_1_50_SSN <- nrow(range_1_50_SSN)
n_range_50_100_SSN <- nrow(range_50_100_SSN)
n_range_100_200_SSN <- nrow(range_100_200_SSN)
n_range_200_500_SSN <- nrow(range_200_500_SSN)
n_range_500_1000_SSN <- nrow(range_500_1000_SSN)
n_range_1000_2000_SSN <- nrow(range_1000_2000_SSN)
n_range_2000_5000_SSN <- nrow(range_2000_5000_SSN)
n_range_5000_SSN <- nrow(range_5000_SSN)

n_ranges_SSN <- c(n_range_1_50_SSN,n_range_50_100_SSN,n_range_100_200_SSN,n_range_200_500_SSN,
	n_range_500_1000_SSN,n_range_1000_2000_SSN,n_range_2000_5000_SSN,n_range_5000_SSN)
n_ranges_SSN <- t(as.data.frame(n_ranges_SSN))
colnames(n_ranges_SSN) <- c("1-50","50-100","100-200","200-500","500-1K","1K-2K","2K-5K",">5K")
rownames(n_ranges_SSN) <- c()

n_range_1_50_OSN <- nrow(range_1_50_OSN)
n_range_50_100_OSN <- nrow(range_50_100_OSN)
n_range_100_200_OSN <- nrow(range_100_200_OSN)
n_range_200_500_OSN <- nrow(range_200_500_OSN)
n_range_500_1000_OSN <- nrow(range_500_1000_OSN)
n_range_1000_2000_OSN <- nrow(range_1000_2000_OSN)
n_range_2000_5000_OSN <- nrow(range_2000_5000_OSN)
n_range_5000_OSN <- nrow(range_5000_OSN)

n_ranges_OSN <- c(n_range_1_50_OSN,n_range_50_100_OSN,n_range_100_200_OSN,n_range_200_500_OSN,
	n_range_500_1000_OSN,n_range_1000_2000_OSN,n_range_2000_5000_OSN,n_range_5000_OSN)
n_ranges_OSN <- t(as.data.frame(n_ranges_OSN))
colnames(n_ranges_OSN) <- c("1-50","50-100","100-200","200-500","500-1K","1K-2K","2K-5K",">5K")
rownames(n_ranges_OSN) <- c()



# Make connected scatter plot
makeScrPlotDistCor <- function(data, plot_title, n_ranges_SSN, n_ranges_OSN) {

	fname <- sprintf('%s.jpg', paste(deparse(substitute(data)), sep="_"))

	# Get number of gene pairs from same strand
    range_1_50_SSN <- n_ranges_SSN[,1]
    range_50_100_SSN <- n_ranges_SSN[,2]
    range_100_200_SSN <- n_ranges_SSN[,3]
    range_200_500_SSN <- n_ranges_SSN[,4]
    range_500_1000_SSN <- n_ranges_SSN[,5]
    range_1000_2000_SSN <- n_ranges_SSN[,6]
    range_2000_5000_SSN <- n_ranges_SSN[,7]
    range_5000_SSN <- n_ranges_SSN[,8]

    # Get number of gene pairs from opposite strands
    range_1_50_OSN <- n_ranges_OSN[,1]
    range_50_100_OSN <- n_ranges_OSN[,2]
    range_100_200_OSN <- n_ranges_OSN[,3]
    range_200_500_OSN <- n_ranges_OSN[,4]
    range_500_1000_OSN <- n_ranges_OSN[,5]
    range_1000_2000_OSN <- n_ranges_OSN[,6]
    range_2000_5000_OSN <- n_ranges_OSN[,7]
    range_5000_OSN <- n_ranges_OSN[,8]
	
	p <- ggplot(data, aes(dist_range,group = 1)) + 
	geom_line(aes(y = med_SSN_pearson, color = "med_SSN_pearson",group = 1), size=1.125) + 
	geom_line(aes(y = med_OSN_pearson, color = "med_OSN_pearson",group = 1), size=1.125) + 
  	geom_point(aes(y = med_SSN_pearson, color = "med_SSN_pearson",group = 1), size=3.25) + 
  	geom_point(aes(y = med_OSN_pearson, color = "med_OSN_pearson",group = 1), size=3.25) + 
  	labs(color="Adjecent gene pair") + 
  	scale_color_manual(labels = c("same strand", "opposite strands"), values = c("#18b3b7", "#f35e5a")) + 
	scale_y_continuous(limits = c(0,0.25)) + 
	scale_x_discrete(breaks=data$dist_range,labels=data$dist_range, 
		limits=c("1-50","50-100","100-200","200-500","500-1K","1K-2K","2K-5K",">5K")) + 
	annotate("text", x=0.76, y=Inf, hjust=0, vjust=10.5, size=5.5, label=range_1_50_SSN, col="gray40") + 
	annotate("text", x=1.82, y=Inf, hjust=0, vjust=13.05, size=5.5, label=range_50_100_SSN, col="gray40") + 
	annotate("text", x=2.76, y=Inf, hjust=0, vjust=14.36, size=5.5, label=range_100_200_SSN, col="gray40") + 
	annotate("text", x=3.74, y=Inf, hjust=0, vjust=16.25, size=5.5, label=range_200_500_SSN, col="gray40") + 
	annotate("text", x=4.76, y=Inf, hjust=0, vjust=21.58, size=5.5, label=range_500_1000_SSN, col="gray40") + 
	annotate("text", x=5.76, y=Inf, hjust=0, vjust=20.57, size=5.5, label=range_1000_2000_SSN, col="gray40") + 
	annotate("text", x=6.80, y=Inf, hjust=0, vjust=20.9, size=5.5, label=range_2000_5000_SSN, col="gray40") + 
	annotate("text", x=7.82, y=Inf, hjust=0, vjust=27.68, size=5.5, label=range_5000_SSN, col="gray40") + 
	annotate("text", x=0.825, y=Inf, hjust=0, vjust=2.975, size=5.5, label=range_1_50_OSN, col="gray40") + 
	annotate("text", x=1.825, y=Inf, hjust=0, vjust=5.4, size=5.5, label=range_50_100_OSN, col="gray40") + 
	annotate("text", x=2.8275, y=Inf, hjust=0, vjust=6.9, size=5.5, label=range_100_200_OSN, col="gray40") + 
	annotate("text", x=3.74, y=Inf, hjust=0, vjust=10.89, size=5.5, label=range_200_500_OSN, col="gray40") + 
	annotate("text", x=4.8275, y=Inf, hjust=0, vjust=10.05, size=5.5, label=range_500_1000_OSN, col="gray40") + 
	annotate("text", x=5.8275, y=Inf, hjust=0, vjust=10.37, size=5.5, label=range_1000_2000_OSN, col="gray40") + 
	annotate("text", x=6.8275, y=Inf, hjust=0, vjust=16.02, size=5.5, label=range_2000_5000_OSN, col="gray40") + 
	annotate("text", x=7.8275, y=Inf, hjust=0, vjust=18.8, size=5.5, label=range_5000_OSN, col="gray40")

	q <- p + ggtitle(plot_title) + theme_bw() + xlab("Intergenic distance (bp)") + ylab("Pearson's r") + 
  		theme(text=element_text(size=16), 
  		axis.ticks.length = unit(.25, "cm"),
  		plot.margin = unit(c(3.0, 10.5, 20, 8), "points"),
		axis.text.x = element_text(colour = "black", size=16, angle=0, margin = margin(t = 9, r = 0, b = 0, l = 0)), 
		axis.text.y = element_text(colour = "black", size=16, angle=0, margin = margin(t = 0, r = 9, b = 0, l = 0)),
		axis.title.x = element_text(colour = "black", size=18.25, margin = margin(t = 28, r = 0, b = 1, l = 0), face ="bold"),
		axis.title.y = element_text(colour = "black", size=18.25, margin = margin(t = 0, r = 27, b = 0, l = 10), face ="bold"),
		plot.title = element_text(colour = "black", size=20.25, margin = margin(t = 25, r = 0, b = 28, l = 0), hjust = 0.5),
		legend.position = "right",
		legend.title = element_text(colour = "black", size=17, face ="bold"),
		legend.text=element_text(size=17), 
		legend.spacing.x = unit(0.25, 'cm'),
		legend.key.size = unit(0.775, "cm"),
  		panel.border = element_rect(colour = "black", fill=NA, size=0.5))

	ggsave(file = file.path(out_dir, "output", "plots", fname), plot = q,
		scale = 1, width = 12.5, height = 8, units = c("in"), 
		dpi = 450, limitsize = FALSE)
}

makeScrPlotDistCor(data=med_dist_cor, plot_title="Correlation between coexpression and intergenic distance", n_ranges_SSN=n_ranges_SSN, n_ranges_OSN=n_ranges_OSN)








