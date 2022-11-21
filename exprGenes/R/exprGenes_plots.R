# This script loads and analysis the data statistics and expression tables for protein-coding 
# genes, coding isoforms, lncRNAs and LTRs and generates the plots for the DevSeq transcriptome  
# single-species expression figures

# NOTES: Add lines 84+85, 102 (colors), 108-109, 115-116, 162 (colors), 168+169, 175+176, 139+140, 229, 245 (colors), 248, 252+253, 260-262, update lines 364-366 + 400-402 (ylim), 425-430 (domain color, add geom_point), lines 420+429 (line size), 541+542, 548, 549, 553 (colors) 710-798, delete 802-853, 804ff (ATH hclust dendrogram), delete 976ff (hclust dendrogram non-ATH species)
#------------------------------------ Read sample tables ------------------------------------


# Read all csv files in input file path
readTable <- function(path, pattern = "*.csv") {
    files = list.files(path, pattern, full.names = TRUE)
    lapply(files, function(x) read.table(x, sep=";", dec=".", header=TRUE, stringsAsFactors=FALSE))
}

stats_tables <- readTable(file.path(out_dir, "output", "mapping_statistics"))
expr_genes_tables <- readTable(file.path(out_dir, "output", "expr_genes"))
ATH_th_genes_repl_inter_counts <- read.table(file.path(in_dir, "Expression_data", "AT_genes_inter_norm_count_mat_vsd_sample_names.csv"), sep=";", dec=".", header=TRUE, stringsAsFactors=FALSE)
ATH_th_genes_repl_intra_counts <- read.table(file.path(in_dir, "Expression_data", "AT_genes_intra_norm_count_mat_vsd_sample_names.csv"), sep=";", dec=".", header=TRUE, stringsAsFactors=FALSE)
ATH_genes_complete_tpm <- read.table(file.path(in_dir, "Expression_data", "AT_genes_complete_table_tpm_sample_names.csv"), sep=";", dec=".", header=TRUE, stringsAsFactors=FALSE)


# Get file names and save them in character vector
stats_table_list <- as.character(list.files(file.path(out_dir, "output", "mapping_statistics"), pattern = "*.csv"))
stats_table_names <- gsub('\\.csv$', '', stats_table_list)
expr_genes_table_list <- as.character(list.files(file.path(out_dir, "output", "expr_genes"), pattern = "*.csv"))
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

	classT_name = as.data.frame(rep(c("Trimmed"), times = number_values))
	names(classT_name) <- "class"
	trimmed <- cbind(classT_name, select(x, Trimmed))
	names(trimmed)[2] <- "reads"

	classM_name = as.data.frame(rep(c("Mapped"), times = number_values))
	names(classM_name) <- "class"
	mapped <- cbind(classM_name, select(x, Mapped))
	names(mapped)[2] <- "reads"

	classD_name = as.data.frame(rep(c("Dedupl."), times = number_values))
	names(classD_name) <- "class"
	deduplicated <- cbind(classD_name, select(x, Deduplicated))
	names(deduplicated)[2] <- "reads"

	mapStats <- rbind(trimmed, mapped, deduplicated)
	return(mapStats)
}

ATH_stats_df <- prepareStats(ATH_stats)
non_ATH_stats_df <- prepareStats(non_ATH_stats)



# Make stats violin plot for ATH
makePlotStatsATH <- function(data, lim_y, plot_title) {

	fname <- sprintf('%s.pdf', paste(deparse(substitute(data)), sep = "_"))

	dedupl <- subset(data, class=="Dedupl.")
	n_dedupl <- paste("n=", nrow(dedupl), sep = "")
	total_dedupl <- paste0(round(sum(as.numeric(dedupl[,2]))/1e9,2),"B")
	total_dedupl = paste(n_dedupl, total_dedupl, sep = "\n")

	data$class <- factor(data$class, levels = unique(data$class))

	ylabels = function(l) {paste0(round(l/1e6,1),"M")}

	stat_sum_single <- function(fun, geom = "crossbar", ...) {

		stat_summary(fun.y = fun, colour = "gray15", geom = geom, size = 1, 
			mapping = aes(ymin = ..y.., ymax = ..y..), width = 0.26, ...)
	}

	p <- ggplot(data, aes(x = class, y = reads, fill = class)) + 
		 geom_violin(trim = TRUE, width = 1.5, size = 1.25, scale = "area", color = "gray15") +
		 stat_sum_single(fun = median) + 
		 scale_y_continuous(limits = c(0,lim_y), expand = c(0, 0), 
		 	labels = function(l) { 
		 		ifelse(l==0, paste0(round(l/1e6,1)),paste0(round(l/1e6,1),"M"))
		 	}) + 
		 scale_x_discrete(labels = c("Dedupl." = "Mapped\n+dedupl")) + 
		 annotate("rect", xmin = 0.25, xmax = 3.85, ymin = 0, ymax = lim_y, fill = "white", alpha = 0, 
		 	color = "black", size = 1.35) + 
		 annotate("text", x = 2.7, y = Inf, hjust = 0, vjust = 1.55, size = 6.5, label = total_dedupl)

	q <- p + scale_fill_manual(values = c("#b2b2b2", "#eeb722", "#5bb1e2")) + theme_minimal() + 
	xlab("") + ylab("Number of PE reads") + ggtitle(plot_title) + 
	geom_hline(yintercept = 30e6, linetype = "dashed", color = "red", size = 1) + 
	theme(legend.position = "none", 
		text = element_text(size = 20.75), 
  		axis.ticks.length = unit(.2, "cm"),
  		axis.ticks = element_line(colour = "gray10", size = 0.8), 
  		axis.line = element_line(colour = "gray10", size = 0.38), 
  		axis.title.y = element_text(colour = "black", size = 20, 
  			margin = margin(t = 0, r = 8.5, b = 0, l = 1.5)), 
  		axis.text.x = element_text(colour = "black", size=18.0, angle = 45, 
  			margin = margin(t = 2.0, r = 0, b = 2.5, l = 0), hjust = 1, vjust = 1),
  		axis.text.y = element_text(colour = "grey50", margin = margin(t = 0, r = 3, b = 0, l = 2)), 
  		plot.title = element_text(colour = "black", size = 21.5, face = "italic", 
  			margin = margin(t = 21.5, r = 0, b = 11.0, l = 0), hjust = 0.5), 
  		plot.margin = unit(c(7.0, 15, 13.55, 16.1), "points"))
	

  	ggsave(file = file.path(out_dir, "output", "plots", fname), plot = q,
		scale = 1, width = 5.0, height = 6.95, units = c("in"), limitsize = FALSE)
}

makePlotStatsATH(data = ATH_stats_df, lim_y = 226000000, plot_title = "A.thaliana") 
# 1 data point for trimmed raw reads above lim_y



# Make stats violin plot for other species
makePlotStatsOS <- function(data, lim_y, plot_title) {

	fname <- sprintf('%s.pdf', paste(deparse(substitute(data)), sep = "_"))

	dedupl <- subset(data, class=="Dedupl.")
	n_dedupl <- paste("n=", nrow(dedupl), sep = "")
	total_dedupl <- paste0(round(sum(as.numeric(dedupl[,2]))/1e9,2),"B")
	total_dedupl = paste(n_dedupl, total_dedupl, sep = "\n")

	data$class <- factor(data$class, levels = unique(data$class))

	stat_sum_single <- function(fun, geom = "crossbar", ...) {

		stat_summary(fun.y = fun, colour = "gray15", geom = geom, size = 1, 
			mapping = aes(ymin = ..y.., ymax = ..y..), width = 0.4075, ...)
	}

	# separate outliers with low reads from violin plot data and plot them individually as dots
	data_wo_outl <- data[c(-505:-507),]
	data_outl <- data[c(505:507),]

	p <- ggplot(data_wo_outl, aes(x = class, y = reads, fill = class)) + 
		 geom_violin(trim = TRUE, width = 1.5, size = 1.25, scale = "area", color = "gray15") + 
		 stat_sum_single(fun = median) + 
		 geom_point(aes(x=3, y=data_outl[1,2]), shape = 21, colour = "gray35", size = 2.25, fill = "white", stroke = 2) + 
		 geom_point(aes(x=3, y=data_outl[2,2]), shape = 21, colour = "gray35", size = 2.25, fill = "white", stroke = 2) + 
		 geom_point(aes(x=3, y=data_outl[3,2]), shape = 21, colour = "gray35", size = 2.25, fill = "white", stroke = 2) + 
		 scale_y_continuous(limits = c(0,lim_y), expand = c(0, 0), 
		 	labels = function(l) { 
		 		ifelse(l==0, paste0(round(l/1e6,1)),paste0(round(l/1e6,1),"M"))
		 	}) + 
		 scale_x_discrete(labels = c("Dedupl." = "Mapped\n+dedupl")) + 
		 annotate("rect", xmin = 0.25, xmax = 3.85, ymin = 0, ymax = lim_y, fill = "white", alpha = 0, 
		 	color = "black", size = 1.35) + 
		 annotate("text", x = 2.7, y = Inf, hjust = 0, vjust = 1.55, size = 6.5, label = total_dedupl)

	q <- p + scale_fill_manual(values = c("#b2b2b2", "#eeb722", "#5bb1e2")) + theme_minimal() + 
	xlab("") + ylab("Number of PE reads") + ggtitle(plot_title) + 
	geom_hline(yintercept = 30e6, linetype = "dashed", color = "red", size = 1) + 
	theme(legend.position = "none", 
		text = element_text(size = 20.75), 
  		axis.ticks.length = unit(.2, "cm"),
  		axis.ticks = element_line(colour = "gray10", size = 0.8), 
  		axis.line = element_line(colour = "gray10", size = 0.38), 
  		axis.title.y = element_text(colour = "black", size = 20, 
  			margin = margin(t = 0, r = 8.5, b = 0, l = 1.5)), 
  		axis.text.x = element_text(colour = "black", size = 18.0, angle = 45, 
  			margin = margin(t = 2.0, r = 0, b = 2.5, l = 0), hjust = 1, vjust = 1), 
  		axis.text.y = element_text(colour = "grey50", margin = margin(t = 0, r = 3, b = 0, l = 2)),  
  		plot.title = element_text(colour = "black", size = 21.5, 
  			margin = margin(t = 21.5, r = 0, b = 8.75, l = 0), hjust = 0.5), 
  		plot.margin = unit(c(7.0, 5.0, 13.55, 26.1), "points"))
	

  	ggsave(file = file.path(out_dir, "output", "plots", fname), plot = q,
		scale = 1, width = 5.0, height = 6.95, units = c("in"), 
		dpi = 600, limitsize = FALSE)
}

makePlotStatsOS(data = non_ATH_stats_df, lim_y = 226000000, plot_title = "Other species")
# 10 data point for trimmed raw reads above lim_y



# Make Deduplicated Read Replicates
makeRepl <- function(x) {

	repl_names_ATH <- data.frame(c("Root.1","Root.2","Root.3","Hypocotyl.1","Hypocotyl.2",
		"Hypocotyl.3","Leaf.1","Leaf.2","Leaf.3","Apex.veg.1","Apex.veg.2","Apex.veg.3",
		"Apex.infl.1","Apex.infl.2","Apex.infl.3","Flower.1","Flower.2","Flower.3","Stamen.1",
		"Stamen.2","Stamen.3","Pollen.1","Pollen.2","Pollen.3","Carpel.1","Carpel.2","Carpel.3"))
	names(repl_names_ATH) <- "Sample_repl"
	repl_names_non_ATH <- c("Root.1","Root.2","Root.3","Hypocotyl.1","Hypocotyl.2","Hypocotyl.3",
		"Leaf.1","Leaf.2","Leaf.3","Apex.veg.1","Apex.veg.2","Apex.veg.3","Apex.infl.1",
		"Apex.infl.2","Apex.infl.3","Flower.1","Flower.2","Flower.3","Pollen.1","Pollen.2",
		"Pollen.3","Carpel.1","Carpel.2","Carpel.3","Stamen.1","Stamen.2","Stamen.3")

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

	level_order <- c("Root.1","Root.2","Root.3","Hypocotyl.1","Hypocotyl.2","Hypocotyl.3",
		"Leaf.1","Leaf.2","Leaf.3","Apex.veg.1","Apex.veg.2","Apex.veg.3","Apex.infl.1",
		"Apex.infl.2","Apex.infl.3","Flower.1","Flower.2","Flower.3","Carpel.1","Carpel.2",
		"Carpel.3","Stamen.1","Stamen.2","Stamen.3","Pollen.1","Pollen.2","Pollen.3")

	samplelabs <- c("Root.1",".2",".3","Hypocot.1",".2",".3",
		"Leaf.1",".2",".3","Apex.v.1",".2",".3","Apex.i.1",
		".2",".3","Flower.1",".2",".3","Carpel.1",".2",
		".3","Stamen.1",".2",".3","Pollen.1",".2",".3")
	data$Species <- factor(data$Species, levels = unique(data$Species))
	species_order <- c("AT","AL","CR","ES","TH","MT","BD")

	p <- ggplot(data, aes(x = factor(Sample_repl, level= level_order), y = Deduplicated, color = Species, group = Species)) + 
	geom_line(aes(x = factor(Sample_repl, level= level_order)), size=1.825) + 
  	geom_point(aes(x = factor(Sample_repl, level= level_order)), size=3.75) + 
  	scale_y_continuous(limits = c(0,7.07e7), expand = c(0, 0), 
		 	labels = function(l) { 
		 		ifelse(l==0, paste0(round(l/1e6,1)),paste0(round(l/1e6,1),"M"))
		 	}) + 
  	scale_x_discrete(labels = samplelabs) + 
  	annotate("rect", xmin=0.25, xmax=27.85, ymin=0, ymax=7.07e7, fill="white", alpha=0, 
		 	color="black", size=0.7) + 
  	labs(color="Species")

	q <- p + ggtitle(plot_title) + theme_bw() + xlab("") + ylab("Number of PE dedupl. reads") + 
	scale_color_manual(values = c("#b2b2b2","#e8a215","#f0d737","#069870","#0770ab","#4fb6f0","#ea6965") 
		# Order of color vector is in alphabetical order of species (AL/AT/BD/CR/ES/MT/TH)
		# It uses a slightly modified colorblind-friendly palette from Wong (Nature Methods, 2011)
		) + 
		guides(colour = guide_legend(nrow = 1)) +  
  		theme(text=element_text(size=21), 
  		axis.ticks.length = unit(.2, "cm"),
  		axis.ticks = element_line(colour = "gray10", size = 0.9), 
  		axis.line = element_line(colour = "gray10", size = 0.9), 
  		panel.grid = element_blank(), 
  		axis.title.y = element_text(colour = "black", size=20, 
  			margin = margin(t = 0, r = 5.85, b = 0, l = 28.5)), 
  		axis.text.x = element_text(colour = "black", size=16.5, angle=45, 
  			margin = margin(t = 0.25, r = 0, b = 0.75, l = 0), hjust = 1, vjust = 1), 
  		axis.text.y = element_text(colour = "black", margin = margin(t = 0, r = 3.95, b = 0, l = 2)), 
  		plot.title = element_text(colour = "black", size = 21.5, 
  			margin = margin(t = 21.5, r = 0, b = 9.45, l = 0), hjust = 0.5), 
  		plot.margin = unit(c(5.5, -3.5, 37.25, 0), "points"),
		legend.position = c(0.331, 0.115),
		legend.background = element_rect(fill = NA),
		legend.key = element_rect(fill = NA),
		legend.title = element_text(colour = "black", size=19.5, face ="bold"),
		legend.text=element_text(size=19.5), 
		legend.spacing.x = unit(0.25, 'cm'),
		legend.key.size = unit(0.775, "cm"),
  		panel.border = element_rect(colour = "black", fill=NA, size=0.5))

	r <- ggplotGrob(q)
	r$layout$clip[r$layout$name=="panel"] <- "off"

	ggsave(file = file.path(out_dir, "output", "plots", fname), plot = r,
		scale = 1, width = 9.8, height = 7.09, units = c("in"), 
		dpi = 600, limitsize = FALSE)
}

plotDedupReads(data=comp_stats_df, plot_title="Comparative samples")




#-------------------------- Plotting expressed genes data for ATH ---------------------------


# Prepare data for ggplot
prepareExprGenes <- function(biotype = c("coding", "NAT", "lincRNA", "LTR", "transcripts"), th_0_01, th_0_05, 
	th_0_1, th_0) {

	df_names <- c("Detailed_name" , "Sample" , "Threshold", "Expressed")

	if (is.element("coding", biotype)) {
		data <- cbind(th_0_01[1,3:ncol(th_0_01)], th_0_05[1,3:ncol(th_0_05)], th_0_1[1,3:ncol(th_0_1)], th_0[1,3:ncol(th_0)])
	} else if (is.element("NAT", biotype)) {
		data <- cbind(th_0_01[2,3:ncol(th_0_01)], th_0_05[2,3:ncol(th_0_05)], th_0_1[2,3:ncol(th_0_1)], th_0[2,3:ncol(th_0)])
	} else  if (is.element("lincRNA", biotype)) {
		data <- cbind(th_0_01[3,3:ncol(th_0_01)], th_0_05[3,3:ncol(th_0_05)], th_0_1[3,3:ncol(th_0_1)], th_0[3,3:ncol(th_0)])
    } else  if (is.element("LTR", biotype)) {
		data <- cbind(th_0_01[4,3:ncol(th_0_01)], th_0_05[4,3:ncol(th_0_05)], th_0_1[4,3:ncol(th_0_1)], th_0[4,3:ncol(th_0)])
	} else  if (is.element("transcripts", biotype)) {
		data <- cbind(th_0_01[1,3:ncol(th_0_01)], th_0_05[1,3:ncol(th_0_05)], th_0_1[1,3:ncol(th_0_1)], th_0[1,3:ncol(th_0)])
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
expr_LTRs_ATH <- prepareExprGenes(biotype = "LTR", th_0_01 = ATH_expr_genes_0.01, 
	th_0_05 = ATH_expr_genes_0.05, th_0_1 = ATH_expr_genes_0.1, th_0 = ATH_expr_genes_0)
expr_transcripts_ATH <- prepareExprGenes(biotype = "transcripts", th_0_01 = ATH_expr_coding_transcripts_0.01, 
	th_0_05 = ATH_expr_coding_transcripts_0.05, th_0_1 = ATH_expr_coding_transcripts_0.1, th_0 = ATH_expr_coding_transcripts_0)




# Plot number of expressed genes at different thresholds for ATH
# This plotting function generates plots with organ instead of individual sample labels
plotExprGenes <- function(data, plot_title, biotype = c("coding","NAT","linc","LTR","iso"), texpr) {

	fname <- sprintf('%s.jpg', paste(deparse(substitute(data)), "domain", sep="_"))

	total_expr <- paste("Total:", texpr, "at 0.05" , sep=" ")

	yLabelsK = function(l) { 
		 		ifelse(l==0, 0, paste0(round(l/1e3,1),"K"))
		 	}

	y_tick_pos <- data.frame(x=c(6.315,10.315,19.315,25.315,29.315,38.315))
	y_depl_pos <- data.frame(x=c(5, 9, 18, 24, 28, 37))

	if (is.element("coding", biotype)) {
		breaksY <- c(1.2e4,1.5e4,1.8e4,2.1e4,2.4e4)
		pltymin <- 1.09e4
		pltymax <- 2.55e4
		xtepos <- 21.75
		x_axt_mar <- 25.55
		y_margin <- margin(t = 0, r = 9, b = 0, l = 1.5)
		y_axs_title <- "Number of Genes"
		y_labels <- yLabelsK
		x_v_adj <- 0

	} else if (is.element("NAT", biotype)) {
		breaksY <- c(0,0.5e3,1e3,1.5e3,2e3,2.5e3)
		pltymin <- -50.75
		pltymax <- 2.825e3
		xtepos <- 31.42
		x_axt_mar <- 25.575
		y_margin <- margin(t = 0, r = 3.8, b = 0, l = 1.55)
		y_axs_title <- "Number of NATs"
		y_labels <- yLabelsK
		x_v_adj <- 0.008

	} else if (is.element("linc", biotype)) {
		breaksY <- c(0,5e2,1e3,1.5e3)
		pltymin <- -20
		pltymax <- 1.102e3
		xtepos <- 31.45
		x_axt_mar <- 25.575
		y_margin <- margin(t = 0, r = 3.8, b = 0, l = 1.55)
		y_axs_title <- "Number of lincRNAs"
		y_labels <- yLabelsK
		x_v_adj <- 0.008

	} else if (is.element("iso", biotype)) {
		breaksY <- c(1.5e4,2e4,2.5e4,3.0e4,3.5e4,4.0e4)
		pltymin <- 1.47e4
		pltymax <- 4.35e4
		xtepos <- 21.75
		x_axt_mar <- 25.55
		y_margin <- margin(t = 0, r = 6.65, b = 0, l = 1.55)
		y_axs_title <- "Number of Transcripts"
		y_labels <- yLabelsK
		x_v_adj <- 0
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

	geom_line(aes(x = factor(Sample, level= level_order)), size=2.15) + 
	scale_y_continuous(limits = c(pltymin,pltymax), breaks = breaksY, expand = c(0, 0), 
		 	labels = y_labels) +
  	annotate("rect", xmin=0.25, xmax=44.75, ymin=pltymin, ymax=pltymax, fill="white", alpha=1, 
		 	color="black", size=0.7) + 
  	annotate("rect", xmin=0.25, xmax=6.5, ymin=pltymin, ymax=pltymax, fill="#d7d7d7", alpha=1) + 
  	annotate("rect", xmin=10.5, xmax=19.5, ymin=pltymin, ymax=pltymax, fill="#d2ebc7", alpha=1) + 
  	annotate("rect", xmin=25.5, xmax=29.5, ymin=pltymin, ymax=pltymax, fill="#d7d7d7", alpha=1) + 
  	annotate("rect", xmin=38.5, xmax=44.75, ymin=pltymin, ymax=pltymax, fill="#fad0c8", alpha=1) +
  	geom_line(aes(x = factor(Sample, level= level_order)), size=2.15) + 
  	geom_point(aes(x = factor(Sample, level= level_order)), size=2.75) + 
  	annotate("text", x = xtepos, y = Inf, hjust = 0, vjust = 22.9, size=7.01, label = total_expr) + 
  	annotate("text", x = 1.675, y = Inf, hjust = 0, vjust = 21.075, size=7.01, label = "Threshold", fontface = 2) + 
  	annotate("text", x = y_tick_pos$x, y = Inf, hjust = 0, vjust = 25.515, size=7.0, label = "I", col="gray10") + 
  	annotate("text", x = y_depl_pos$x, y = Inf, hjust = 0, vjust = 8.44+x_v_adj, size=20.5, label = "_", col="white", fontface = 2) + 
  	annotate("text", x = 1.85, y = Inf, hjust = 0, vjust = x_axt_mar, size=7.2, label = "Root") + 
  	annotate("text", x = 6.85, y = Inf, hjust = 0, vjust = x_axt_mar, size= 7.2, label = "Stem") + 
  	annotate("text", x = 13.57, y = Inf, hjust = 0, vjust = x_axt_mar, size= 7.2, label = "Leaf") + 
  	annotate("text", x = 20.88, y = Inf, hjust = 0, vjust = x_axt_mar, size= 7.2, label = "Apex") + 
  	annotate("text", x = 25.95, y = Inf, hjust = 0, vjust = x_axt_mar, size= 7.2, label = "Flow.") + 
  	annotate("text", x = 30.29, y = Inf, hjust = 0, vjust = x_axt_mar, size= 7.2, label = "Floral organ") + 
  	annotate("text", x = 40.2, y = Inf, hjust = 0, vjust = x_axt_mar, size= 7.2, label = "Fruit") + 
  	labs(color="")

	q <- p + theme_bw() + labs(title = plot_title,  x = "Organ", y = y_axs_title) + 
	scale_color_manual(values=c("gray35","#fe5651","#967cee","#dea80c")) + 
		guides(colour = guide_legend(nrow = 1)) + 
  		theme(text = element_text(size=23.5), 
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
  		panel.grid.major = element_line(colour = "white"), 
  		panel.grid.minor = element_line(colour = "white"),  
  		axis.ticks.length = unit(.3, "cm"),
  		axis.ticks = element_line(colour = "gray10", size = 0.925), 
  		axis.line = element_line(colour = "gray10", size = 0.8),
  		axis.title.x = element_text(colour = "black", size=21.5, 
  			margin = margin(t = 30.25, r = 0, b = 50, l = 0)),  
  		axis.title.y = element_text(colour = "black", size=21.5, 
  			margin = y_margin), 
  		axis.text.y = element_text(colour = "black", size=18.5, margin = margin(t = 0, r = 3, b = 0, l = 2)), 
  		plot.title = element_text(colour = "black", size=22.88, 
  			margin = margin(t = 37.1, r = 0, b = 8.5, l = 0), hjust = 0.5), 
  		plot.margin = unit(c(0, 15, 0, 1), "points"),
		legend.position = c(0.22, 0.108),
		legend.title = element_text(colour = "black", size=20, face ="bold"),
		legend.text = element_text(size=20), 
		legend.key.size = unit(0.775, "cm"),
		legend.key.height = unit(0.4, "cm"),
		legend.background = element_rect(fill = NA),
		legend.key = element_rect(fill = NA),
  		panel.border = element_rect(colour = "gray10", fill=NA, size=0.9))

	r <- ggplotGrob(q)
	r$layout$clip[r$layout$name=="panel"] <- "off"

	ggsave(file = file.path(out_dir, "output", "plots", fname), plot = r,
		scale = 1, width = 10.25, height = 7.3, units = c("in"), 
		dpi = 600, limitsize = FALSE)
}
plotExprGenes(data=expr_coding_genes_ATH, plot_title=expression(paste("Expressed protein-coding genes in ", italic("A.thaliana"))), biotype = "coding", texpr=ATH_expr_genes_0.05[1,2])
plotExprGenes(data=expr_NATs_ATH, plot_title=expression(paste("Expressed NATs in ", italic("A.thaliana"))), biotype = "NAT", texpr=ATH_expr_genes_0.05[2,2])
plotExprGenes(data=expr_lincRNAs_ATH, plot_title=expression(paste("Expressed lincRNAs in ", italic("A.thaliana"))), biotype = "linc", texpr=ATH_expr_genes_0.05[3,2])
plotExprGenes(data=expr_transcripts_ATH, plot_title=expression(paste("Expressed protein-coding transcripts in ", italic("A.thaliana"))), biotype = "iso", texpr=ATH_expr_coding_transcripts_0.05[1,2])




#----------------------------- Plotting replicate correlations ------------------------------


# Prepare data for ggplot
prepareReplStats <- function(ATH,AL,CR,ES,TH,MT,BD) {

	number_values_ATH <- ncol(ATH)
	number_values_AL <- ncol(AL)
	number_values_other <- ncol(CR)

	class_ATH_key = as.data.frame(rep(c("AT"), times = number_values_ATH))
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

all_spec_repl_df <- prepareReplStats(ATH=ATH_repl_corr_counts_0.05, AL=AL_repl_corr_counts_0.05, 
	CR=CR_repl_corr_counts_0.05, ES=ES_repl_corr_counts_0.05, TH=TH_repl_corr_counts_0.05, 
	MT=MT_repl_corr_counts_0.05, BD=BD_repl_corr_counts_0.05)



# Make replicate correlation plot
makePlotReplCorr <- function(data, plot_title) {

	fname <- sprintf('%s.jpg', paste(deparse(substitute(data)), sep="_"))

	data$Species <- factor(data$Species, levels = unique(data$Species))

	p <- ggplot(data, aes(x=Species, y=Correlation, fill=Species)) + 
	     stat_boxplot(geom ='errorbar', width = 0.45, size=1.0, color="gray15") + 
		 geom_boxplot(width = 0.75, size=1.0, color="gray15", outlier.shape = 21, 
		 	outlier.size = 2.5, outlier.stroke = 1.5, outlier.fill = NA, outlier.color="gray35") + 
		 scale_y_continuous(limits = c(0.9692, 1.0007), expand = c(0, 0)) + 
		 annotate("rect", xmin=0.35, xmax=7.65, ymin=0.9692, ymax=1.0007, fill="white", alpha=0,  
		 	color="black", size=1.35)

	q <- p + scale_fill_manual(values=c("#b2b2b2","#e8a215","#f0d737","#069870","#0770ab","#4fb6f0","#ea6965")) + 
	# Uses a slightly modified colorblind-friendly palette from Wong (Nature Methods, 2011)
	theme_minimal() + 
	xlab("Species") + ylab("Pearson's r") + ggtitle(plot_title) + 
	theme(text = element_text(size=23.5), 
  		panel.grid.major = element_line(colour = "white"), 
  		panel.grid.minor = element_line(colour = "white"),  
  		axis.ticks.length = unit(.215, "cm"),
  		axis.ticks = element_line(colour = "gray10", size = 0.9),
  		axis.line = element_line(colour = "gray10", size = 0.55),
  		axis.title.x = element_text(colour = "black", size=21.55, 
  			margin = margin(t = 12.5, r = 0, b = 50.2, l = 0)),  
  		axis.title.y = element_text(colour = "black", size=21.5, 
  			margin = margin(t = 0, r = 5.8, b = 0, l = 1.5)), 
  		axis.text.x = element_text(colour = "black", margin = margin(t = 3.5, r = 0, b = 1.6, l = 0), size=20.5), 
  		axis.text.y = element_text(colour = "black", margin = margin(t = 0, r = 3.25, b = 0, l = 4), size=18.55), 
  		plot.title = element_text(colour = "black", size = 22.85, 
  			margin = margin(t = 36.25, r = 0, b = 11.05, l = 0), hjust = 0.5), 
  		plot.margin = unit(c(0, 0.5, 0, 1), "points"),
		legend.position = "none",
		legend.title = element_text(colour = "black", size=20, face ="bold"),
		legend.text = element_text(size=20), 
		legend.background = element_rect(fill = NA),
  		panel.border = element_rect(colour = "gray10", fill = NA, size = 0.5))

  	ggsave(file = file.path(out_dir, "output", "plots", fname), plot = q,
		scale = 1, width = 10.25, height = 7.3, units = c("in"), 
		dpi = 600, limitsize = FALSE)
}

makePlotReplCorr(data=all_spec_repl_df, plot_title="Replicate correlations")




#------------ Plotting expressed genes data for non-ATH species with facet wrap -------------


# Prepare data for ggplot
prepareExprGenesOS <- function(species=c("AL","CR","ES","TH","MT","BD"), 
	th_0, th_0_01, th_0_05, th_0_1) {


	th_0_01 <- th_0_01[,c(1:8,10:11,9)]
	th_0_05 <- th_0_05[,c(1:8,10:11,9)]
	th_0_1 <- th_0_1[,c(1:8,10:11,9)]
	th_0 <- th_0[,c(1:8,10:11,9)]


	species_df <- as.data.frame(rep(c(species), 180))
	colnames(species_df) <- "species"


	th_list <- list(th_0.01=th_0_01, th_0.05=th_0_05, th_0.1=th_0_1, th_0.5=th_0)

	getNumExpr <- function(df, trans_type) {

		if (trans_type == "Genes") {
			sel_row <- 1
		} else if (trans_type == "NATs") {
			sel_row <- 2
		} else if (trans_type == "lincRNAs") {
			sel_row <- 3
		} else if (trans_type == "LTR TEs") {
			sel_row <- 4
		} else if (trans_type == "Transcripts") {
			sel_row <- 5
		}

		n_expr <- df[sel_row, 3:ncol(df)]
		colnames(n_expr) <- NULL
		n_expr <- as.data.frame(t(n_expr))
		colnames(n_expr) <- "expressed"
		total_expressed <- df[sel_row, 2]
		total_expressed <- as.data.frame(rep(c(total_expressed), 9))
		colnames(total_expressed) <- "total_expressed"
		transcript_class <- as.data.frame(rep(c(trans_type), 9))
		colnames(transcript_class) <- "class"
		biotype <- df[sel_row, 1]
		biotype <- as.data.frame(rep(c(biotype), 9))
		colnames(biotype) <- "biotype"
		value_class <- cbind(n_expr, total_expressed, biotype, transcript_class)

		return(value_class)

	}

	cd_value_class <- do.call(rbind, lapply(th_list, getNumExpr, trans_type = "Genes"))
	threshold <- data.frame(Threshold=gsub(".*_","", sub(".[^.]+$", "", rownames(cd_value_class))))
	cd_value_class <- data.frame(cd_value_class[,1:2], threshold, cd_value_class[,3:4])

	NAT_value_class <- do.call(rbind, lapply(th_list, getNumExpr, trans_type = "NATs"))
	threshold <- data.frame(Threshold=gsub(".*_","", sub(".[^.]+$", "", rownames(NAT_value_class))))
	NAT_value_class <- data.frame(NAT_value_class[,1:2], threshold, NAT_value_class[,3:4])

	linc_value_class <- do.call(rbind, lapply(th_list, getNumExpr, trans_type = "lincRNAs"))
	threshold <- data.frame(Threshold=gsub(".*_","", sub(".[^.]+$", "", rownames(linc_value_class))))
	linc_value_class <- data.frame(linc_value_class[,1:2], threshold, linc_value_class[,3:4])

	circ_value_class <- do.call(rbind, lapply(th_list, getNumExpr, trans_type = "LTR TEs"))
	threshold <- data.frame(Threshold=gsub(".*_","", sub(".[^.]+$", "", rownames(circ_value_class))))
	circ_value_class <- data.frame(circ_value_class[,1:2], threshold, circ_value_class[,3:4])

	iso_value_class <- do.call(rbind, lapply(th_list, getNumExpr, trans_type = "Transcripts"))
	threshold <- data.frame(Threshold=gsub(".*_","", sub(".[^.]+$", "", rownames(iso_value_class))))
	iso_value_class <- data.frame(iso_value_class[,1:2], threshold, iso_value_class[,3:4])
	

	detailed_sample_name <- names(th_0)[3:ncol(th_0)]
	detailed_sample_name <- as.data.frame(rep(detailed_sample_name, times = 20))
	colnames(detailed_sample_name) <- "detailed_sample_name"

	sample_names <- c("Root", "Hypocotyl", "Leaf", "Apex veg", "Apex inf", 
		"Flower", "Carpel", "Stamen", "Pollen")

	if (species == "BD") {

		sample_names <- c("Root_b", "Mesocotyl_b", "Leaf_b", "Apex veg_b", "Spikelet m_b", 
		"Floret_b", "Carpel_b", "Stamen_b", "Pollen_b")
	}

	sample_names <- as.data.frame(rep(sample_names, times = 20))
	colnames(sample_names) <- "sample_names"


	expr_df <- rbind(cd_value_class, iso_value_class, NAT_value_class, linc_value_class, circ_value_class) 
	expr_df_ext <- cbind(expr_df, species_df, sample_names, detailed_sample_name)

	return(expr_df_ext)
}


# Get expressed genes for each species
expr_genes_AL <- prepareExprGenesOS(species = "AL", th_0 = rbind(AL_expr_genes_0, AL_expr_coding_transcripts_0),
	th_0_01 = rbind(AL_expr_genes_0.01, AL_expr_coding_transcripts_0.01), th_0_05 = rbind(AL_expr_genes_0.05, AL_expr_coding_transcripts_0.05), 
	th_0_1 = rbind(AL_expr_genes_0.1, AL_expr_coding_transcripts_0.1))

expr_genes_CR <- prepareExprGenesOS(species = "CR", th_0 = rbind(CR_expr_genes_0, CR_expr_coding_transcripts_0),
	th_0_01 = rbind(CR_expr_genes_0.01, CR_expr_coding_transcripts_0.01), th_0_05 = rbind(CR_expr_genes_0.05, CR_expr_coding_transcripts_0.05), 
	th_0_1 = rbind(CR_expr_genes_0.1, CR_expr_coding_transcripts_0.1))

expr_genes_ES <- prepareExprGenesOS(species = "ES", th_0 = rbind(ES_expr_genes_0, ES_expr_coding_transcripts_0),
	th_0_01 = rbind(ES_expr_genes_0.01, ES_expr_coding_transcripts_0.01), th_0_05 = rbind(ES_expr_genes_0.05, ES_expr_coding_transcripts_0.05), 
	th_0_1 = rbind(ES_expr_genes_0.1, ES_expr_coding_transcripts_0.1))

expr_genes_TH <- prepareExprGenesOS(species = "TH", th_0 = rbind(TH_expr_genes_0, TH_expr_coding_transcripts_0),
	th_0_01 = rbind(TH_expr_genes_0.01, TH_expr_coding_transcripts_0.01), th_0_05 = rbind(TH_expr_genes_0.05, TH_expr_coding_transcripts_0.05), 
	th_0_1 = rbind(TH_expr_genes_0.1, TH_expr_coding_transcripts_0.1))

expr_genes_MT <- prepareExprGenesOS(species = "MT", th_0 = rbind(MT_expr_genes_0, MT_expr_coding_transcripts_0),
	th_0_01 = rbind(MT_expr_genes_0.01, MT_expr_coding_transcripts_0.01), th_0_05 = rbind(MT_expr_genes_0.05, MT_expr_coding_transcripts_0.05), 
	th_0_1 = rbind(MT_expr_genes_0.1, MT_expr_coding_transcripts_0.1))

expr_genes_BD <- prepareExprGenesOS(species = "BD", th_0 = rbind(BD_expr_genes_0, BD_expr_coding_transcripts_0),
	th_0_01 = rbind(BD_expr_genes_0.01, BD_expr_coding_transcripts_0.01), th_0_05 = rbind(BD_expr_genes_0.05, BD_expr_coding_transcripts_0.05), 
	th_0_1 = rbind(BD_expr_genes_0.1, BD_expr_coding_transcripts_0.1))


# Combine expr_genes data from all species
expr_genes_OS <- rbind(expr_genes_AL, expr_genes_CR, expr_genes_ES, expr_genes_TH, expr_genes_MT, 
	expr_genes_BD)




# Plot number of expressed genes for each species
plotExprGenesOS <- function(data) {

    trans_class <- unique(data$class)

    if (unique(data$class) == "Genes") {
        y_scale_factor <- 1
        th_label <- "Total expressed:"
    } else if (unique(data$class) == "Transcripts") {
        y_scale_factor <- 1
        th_label <- ""
    } else if (unique(data$class) == "NATs") {
        y_scale_factor <- 0.85
        th_label <- ""
    } else {
        y_scale_factor <- 0.74
        th_label <- ""
    }

    expr_genes_df <- subset(data, data$Threshold == 0.05)
    expr_genes_df <- expr_genes_df[seq(1, nrow(expr_genes_df), 
        nrow(expr_genes_df)/length(unique(expr_genes_df$species))), ]
    expr_genes_df <- subset(expr_genes_df[c("species", "class", "Threshold", "total_expressed")])
    expr_genes_df$x <- rep(1.25)
    expr_genes_df$y <- rep(y_scale_factor*sapply(
        split(data$expressed, rep(1:length(unique(data$species)), each = nrow(data)/length(unique(data$species)))), min), 
    each=length(unique(expr_genes_df$Threshold)))
    expr_genes_df$value <- paste0("n = ", expr_genes_df$total_expressed, " (0.05)")
    expr_genes_df$th_value <- c(th_label, rep( "", nrow(expr_genes_df)-1))

    fname <- paste0("other_species_", trans_class, ".jpg")

    y_label_form <- function(l) { 
        ifelse(l<100, l, paste0(round(l/1e3,1),"K"))
    }

    data$sample_names <- case_when(data$sample_names == "Root" ~ "Rt", data$sample_names == "Root_b" ~ "Rt ", data$sample_names == "Hypocotyl" ~ "Hy", data$sample_names == "Mesocotyl_b" ~ "Me ", 
        data$sample_names == "Leaf" ~ "Lf", data$sample_names == "Leaf_b" ~ "Lf ", data$sample_names == "Apex veg" ~ "Av", data$sample_names == "Apex veg_b" ~ "Av ", data$sample_names == "Apex inf" ~ "Ai", data$sample_names == "Spikelet m_b" ~ "Sm ", 
        data$sample_names == "Flower" ~ "Fl", data$sample_names == "Floret_b" ~ "Fl ", data$sample_names == "Carpel" ~ "Ca", data$sample_names == "Carpel_b" ~ "Ca ", data$sample_names == "Stamen" ~ "St", data$sample_names == "Stamen_b" ~ "St ", 
        data$sample_names == "Pollen" ~ "Pl", data$sample_names == "Pollen_b" ~ "Pl ")

    data$species <- factor(data$species, levels = unique(data$species))
    data$sample_names <- factor(data$sample_names, levels = unique(data$sample_names))
    p <- ggplot(data = data, color = Threshold, aes(x = sample_names, y = expressed)) + 
            geom_line(size = 2.45, data = data, aes(x = sample_names, y = expressed, group = Threshold, color = Threshold)) + 
            geom_point(size = 3.25, data = data, aes(x = sample_names, y = expressed, group = Threshold, color = Threshold)) + 
            scale_y_continuous(expand = c(0.1, 0), breaks = pretty_breaks(n = 4), labels= y_label_form) + 
            scale_color_manual(values = c("gray35", "#fe5651", "#967cee", "#dea80c")) + 
            scale_x_discrete(expand = c(0.05, 0)) + 
            guides(shape = guide_legend(override.aes = list(stroke = 7.75)))

            q <- p + theme_classic() + xlab("") + ylab(paste(gsub('.{1}$', '', trans_class), "count", sep=" ")) + 
            geom_text(data = expr_genes_df, mapping = aes(x = x, y = y, label = value), 
                size = 7.9, vjust = 0.28, hjust = 0, color = "grey35") + 
            geom_text(data = expr_genes_df, mapping = aes(x = x, y = y, label = th_value), 
                size = 7.9, vjust = -1.45, hjust = 0.02, color = "grey35") + 
            theme(text=element_text(size = 16), 
                strip.text = element_text(size = 22.75), 
                strip.text.x = element_text(margin = margin(0.4, 0, 0.4, 0, "cm")), 
                strip.background = element_rect(colour = 'black', fill = NA, size = 2.5), 
                axis.ticks.length = unit(0.25, "cm"), 
                axis.ticks = element_line(colour = "black", size = 1.1), 
                axis.line = element_line(colour = 'black', size = 1.1), 
                plot.margin = unit(c(0.75, 4.5, 1.07, 4.05),"cm"), 
                axis.title.y = element_text(size = 24.25, margin = margin(t = 0, r = 7, b = 0, l = 12.5), 
                    colour = "black", face = "plain"), 
                axis.title.x = element_text(size = 24.25, margin = margin(t = 6.5, r = 0, b = 5.75, l = 0), 
                    colour = "black", face = "plain"), 
                axis.text.x = element_text(size=20.9, margin = margin(t = 4, b = 7.75), colour = "grey35", 
                    angle = 0, vjust = 1, hjust = 0.5), 
                axis.text.y = element_text(size = 20.9, angle = 0, margin = margin(l = 0.75, r = 1.5), colour = "grey35"), 
                plot.title = element_text(size = 27.35, margin = margin(t = 0, b = 15), face = "plain"), 
                panel.spacing = unit(0.2, "cm"), 
                panel.grid.major = element_blank(),
                panel.grid.minor.x = element_blank(), 
                panel.grid.minor.y = element_blank(), 
                legend.margin = margin(t = -1.0, b = 2.0, unit = "cm"), 
                legend.text = element_text(size = 25.0), 
                legend.title = element_text(size = 25.0), 
                legend.key.size = unit(2, "line"), 
                legend.position = "bottom") 

            q <- q + facet_wrap(~ factor(species, levels = c("AL", "CR", "ES", "TH", "MT", "BD")) , nrow = 1, scales = "free")

            ggsave(file = file.path(out_dir, "output", "plots", fname), plot = q, 
                width = 28.5, height = 6.5, dpi = 300, units = c("in"), limitsize = FALSE)

}
trans_class <- c("Genes", "Transcripts", "NATs", "lincRNAs") # split expr_genes_OS df into transcript class list
trans_class_ls <- setNames(as.list(trans_class), trans_class)
trans_class_ls <- lapply(trans_class_ls, function(x){dplyr::filter(expr_genes_OS, grepl(x, class))})

lapply(trans_class_ls, plotExprGenesOS)




#-------------------------- Do hclust dendrograms for all species ---------------------------


# General settings
# Get pollen samples from intra count table
AT_intra_count_pollen <- select(ATH_th_genes_repl_intra_counts, 
	"flowers_mature_pollen_28d_.1.", "flowers_mature_pollen_28d_.2.", "flowers_mature_pollen_28d_.3.")
AT_intra_count_pollen <- tibble::rownames_to_column(AT_intra_count_pollen, "gene_id")
ATH_th_genes_repl_inter_counts <- tibble::rownames_to_column(ATH_th_genes_repl_inter_counts, "gene_id")
ATH_th_genes_repl_counts <- merge(ATH_th_genes_repl_inter_counts, AT_intra_count_pollen, by = "gene_id")

# Get biotypes
ATH_genes_complete_tpm <- select(ATH_genes_complete_tpm, "id", "biotype", "source")
colnames(ATH_genes_complete_tpm)[which(names(ATH_genes_complete_tpm) == "id")] <- "gene_id"
ATH_th_genes_repl_counts <- merge(ATH_genes_complete_tpm, ATH_th_genes_repl_counts, by = "gene_id")

# Create shorter sample descriptions
AT_names <- rep(c("root_root tip 5d", "root mat.zone 5d", "root whole rt 5d", "root whole rt 7d",
 "root whole rt 14d", "root whole rt 21d", "hypocotyl 10d", "internode 3rd 24d", "internode 2nd 24d", 
 "internode 1st 28d", "cotyledons 7d", "leaf 1+2 7d", "leaf 1+2 10d", "leaf petiole 10d", 
 "leaf tip 10d", "leaf 5+6 17d", "leaf 9+10 27d", "leaf senesc 35d", "cauline leaf 24d", 
 "apex veg 7d", "apex veg 10d", "apex veg 14d", "apex inf 21d", "apex inf clv1 21d", "apex inf 28d", 
 "flower st9", "flower st10/11", "flower st12", "flower st15", "sepals st12", "sepals st15", 
 "petals st12", "petals st15", "stamen st12", "stamen st15", "carpel early st12", 
 "carpel late st12", "fruit st15", "fruit st16", "fruit st17a", "seeds st16", "seeds st17a", 
 "seeds st18", "pollen mature"), each=3)

replicate_tag <- rep(c(" 1"," 2"," 3"), times=44)

AT_names <- paste0(AT_names,replicate_tag)

names(ATH_th_genes_repl_counts)[4:ncol(ATH_th_genes_repl_counts)] <- AT_names


# Define colors based on sample name
label_col <- c(roo="#5d4a95", hyp="#5bb1e2", int="#0c703d", lea="#00994f", cot="#00994f", 
	cau="#00994f", ape="#f0d737", flo="#de6daf", sep="#84cd6a", pet="#ead1c7", 
	sta="#f23d29", pol="#a63126", car="#e8a215", fru="#b54185", see="#e9a3b3")


# apex color for hclust dendrogram
# ape="#f4dc28"

# color setting for comparative heatmap and PCA: 
# split apex samles => apex veg="#95b73a", apex inf="#fad819" 
# hypocotyl slightly lighter hyp="#8591c7" root lighter roo="#5850a3"


# Generate hclust dendrogram using relative expression data
makeDendrogram <- function(x, coefficient = c("pearson", "spearman"), label_col) {

	# Show error message if no scaling is chosen
	if (missing(coefficient))

		stop(
			"Please choose one of the following coefficients: 
			'pearson', 'spearman'",
			call. = TRUE
			)

	# Set filename
    dfname <- deparse(substitute(x))
    coefficient_tag <- match.arg(coefficient)
    fname <- sprintf('%s_dend.jpg', paste(dfname, coefficient_tag, sep="_"))


	x_cd <- subset(x, biotype == "protein_coding")
	x_as <- x[x$biotype %like% "antisense", ]
	x_li <- subset(x, biotype == "lnc_intergenic")

    df_t_cd <- t(x_cd[, 4:ncol(x_cd)]) # transposes data frame so rows become columns and vice versa
    df_t_cd[is.na(df_t_cd)] <- 0 # replaces NAs by 0
    df_t_as <- t(x_as[, 4:ncol(x_as)]) # transposes data frame so rows become columns and vice versa
    df_t_as[is.na(df_t_as)] <- 0 # replaces NAs by 0
    df_t_li <- t(x_li[, 4:ncol(x_li)]) # transposes data frame so rows become columns and vice versa
    df_t_li[is.na(df_t_li)] <- 0 # replaces NAs by 0

    # Build distance matrix
    if (is.element(coefficient, c("pearson"))) {
        df_t_cd_dist.mat <- get_dist(df_t_cd, stand = FALSE, method = "pearson")
        df_t_as_dist.mat <- get_dist(df_t_as, stand = FALSE, method = "pearson")
        df_t_li_dist.mat <- get_dist(df_t_li, stand = FALSE, method = "pearson")

    } else if (is.element(coefficient, c("spearman"))) {
      df_t_cd_dist.mat <- get_dist(df_t_cd, stand = FALSE, method = "spearman")
      df_t_as_dist.mat <- get_dist(df_t_as, stand = FALSE, method = "spearman")
      df_t_li_dist.mat <- get_dist(df_t_li, stand = FALSE, method = "spearman")
    } 

    df_cd_clust.res <- hclust(df_t_cd_dist.mat, method = "average") # agglomerate clustering using average linkage
    df_as_clust.res <- hclust(df_t_as_dist.mat, method = "average") # agglomerate clustering using average linkage
    df_li_clust.res <- hclust(df_t_li_dist.mat, method = "average") # agglomerate clustering using average linkage
  
    df_dend_cd <- dendrapply(as.dendrogram(df_cd_clust.res), function(n){
    
    if (is.leaf(n)){
      dend_col <- label_col[substr(attr(n,"label"),1,3)]
      attr(n, "nodePar") <- list(pch = NA, lab.col = dend_col) # to define label color
      attr(n, "edgePar") <- list(col = dend_col) # to color branch
      }
    return(n)
    })

    df_dend_as <- dendrapply(as.dendrogram(df_as_clust.res), function(n){
    
    if (is.leaf(n)){
      dend_col <- label_col[substr(attr(n,"label"),1,3)]
      attr(n, "nodePar") <- list(pch = NA, lab.col = dend_col) # to define label color
      attr(n, "edgePar") <- list(col = dend_col) # to color branch
      }
    return(n)
    })

    df_dend_li <- dendrapply(as.dendrogram(df_li_clust.res), function(n){
    
    if (is.leaf(n)){
      dend_col <- label_col[substr(attr(n,"label"),1,3)]
      attr(n, "nodePar") <- list(pch = NA, lab.col = dend_col) # to define label color
      attr(n, "edgePar") <- list(col = dend_col) # to color branch
      }
    return(n)
    })

    # make branch colors extend to last common node
    brc_col_cd <- label_col[substr(colnames(x_cd[, 4:ncol(x_cd)]),1,3)]
    brc_col_cd <- brc_col_cd[order.dendrogram(df_dend_cd)]
    brc_col_cd <- factor(brc_col_cd, unique(brc_col_cd))
    brc_col_as <- label_col[substr(colnames(x_as[, 4:ncol(x_as)]),1,3)]
    brc_col_as <- brc_col_as[order.dendrogram(df_dend_as)]
    brc_col_as <- factor(brc_col_as, unique(brc_col_as))
    brc_col_li <- label_col[substr(colnames(x_li[, 4:ncol(x_li)]),1,3)]
    brc_col_li <- brc_col_li[order.dendrogram(df_dend_li)]
    brc_col_li <- factor(brc_col_li, unique(brc_col_li))

    jpeg(height = 7.3, width = 10.25, pointsize = 10, units = c("in"), res = 300,  
    	file = file.path(out_dir, "output", "plots", fname))
    par(mar = c(5.0, 3.5, 5.3, 8.8), lwd = 5, cex = 0.1, cex.axis = 3.5, cex.main = 3.5)
    df_dend_cd = color_branches(df_dend_cd, clusters = as.numeric(brc_col_cd), col = levels(brc_col_cd))
    df_dend_as = color_branches(df_dend_as, clusters = as.numeric(brc_col_as), col = levels(brc_col_as))
    df_dend_li = color_branches(df_dend_li, clusters = as.numeric(brc_col_li), col = levels(brc_col_li))

    df_dend_cd <- rotate(df_dend_cd, c(1:15,19:33,16:18,34:93,106:114,100:105,94:99,115:120,124:132,121:123))
    df_dend_as <- rotate(df_dend_as, c(130:132,127:129,1:24,28:36,25:27,37:90,99:107,96:98,91:95,108:120,124:126,121:123))
    df_dend_li <- rotate(df_dend_li, c(124:132,1:3,7:12,4:6,13:30,85:123,31:45,67:84,58:66,46:57))


    # Get color vector for reordered dendrogram
    brc_col_cd <- label_col[substr(colnames(x_cd[, 4:ncol(x_cd)]),1,3)]
    brc_col_cd <- brc_col_cd[order.dendrogram(df_dend_cd)]
    brc_col_cd <- factor(brc_col_cd, unique(brc_col_cd))
    brc_col_as <- label_col[substr(colnames(x_as[, 4:ncol(x_as)]),1,3)]
    brc_col_as <- brc_col_as[order.dendrogram(df_dend_as)]
    brc_col_as <- factor(brc_col_as, unique(brc_col_as))
    brc_col_li <- label_col[substr(colnames(x_li[, 4:ncol(x_li)]),1,3)]
    brc_col_li <- brc_col_li[order.dendrogram(df_dend_li)]
    brc_col_li <- factor(brc_col_li, unique(brc_col_li))


    par(mfrow = c(3,1))

    df_dend_cd %>% set("labels_cex", 0.1) %>% plot(main = "Protein-coding")
    colored_bars(colors = brc_col_cd, dend = df_dend_cd, sort_by_labels_order = FALSE, y_shift = -0.035, y_scale = 0.340, rowLabels = "")
    par(mfrow = c(3,1), new = TRUE, mfg = c(2, 1))
    df_dend_as %>% set("labels_cex", 0.1) %>% plot(main = "NATs")
    colored_bars(colors = brc_col_as, dend = df_dend_as, sort_by_labels_order = FALSE, y_shift = -0.0225, y_scale = 0.253, rowLabels = "")
    par(mfrow = c(3,1), new = TRUE)
    df_dend_li %>% set("labels_cex", 0.1) %>% plot(main = "lincRNAs")
    colored_bars(colors = brc_col_li, dend = df_dend_li, sort_by_labels_order = FALSE, y_shift = -0.0175, y_scale = 0.2025, rowLabels = "")

    dev.off()
}

makeDendrogram(ATH_th_genes_repl_counts, coefficient = "pearson", label_col = label_col)







