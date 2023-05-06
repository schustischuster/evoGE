# Plot results of cisNAT analysis
# Input tables: (1) cisNAT/PC pairwise correlation tables ("output", "NAT_expr_cor"), 
# (2) overlapping PC/PC pairwise cor and feature tables ("output", "overlap_cd_genes"),
# (3) cisNAT/PC gene pair feature tables (total length/overlap length), 
# (4) neigbouring PC/PC gene pair tables (pairwise correlation)


#------------------- Load packages, set directories and read sample tables ---------------------

# Set file path and input files
in_dir_NAT_cor <- file.path(out_dir, "output", "NAT_expr_cor")
in_dir_PC_pairs <- file.path(out_dir, "output", "overlap_pc_genes")
in_dir_NCPC_pairs <- file.path(out_dir, "output", "overlap_nc_genes")


# Read all csv files in input file path
readTable <- function(path, pattern = "*.csv") {
    files = list.files(path, pattern, full.names = TRUE)
    lapply(files, function(x) read.table(x, sep=";", dec=".", header=TRUE, stringsAsFactors=FALSE))
}

NAT_expr_cor_ls <- readTable(in_dir_NAT_cor)
overlap_cd_genes_ls <- readTable(in_dir_PC_pairs)
overlap_nccd_genes_ls <- readTable(in_dir_NCPC_pairs)


# Get file names and save them in character vector
NAT_gene_tables_list <- as.character(list.files(in_dir_NAT_cor, pattern = "*.csv"))
NAT_gene_tables_names <- gsub('\\.csv$', '', NAT_gene_tables_list)

coding_gene_tables_list <- as.character(list.files(in_dir_PC_pairs, pattern = "*.csv"))
coding_gene_tables_names <- gsub('\\.csv$', '', coding_gene_tables_list)

nc_cd_gene_tables_list <- as.character(list.files(in_dir_NCPC_pairs, pattern = "*.csv"))
nc_cd_gene_tables_names <- gsub('\\.csv$', '', nc_cd_gene_tables_list)


# Change data frame names in list
names(NAT_expr_cor_ls) <- NAT_gene_tables_names
list2env(NAT_expr_cor_ls, envir = .GlobalEnv)

names(overlap_cd_genes_ls) <- coding_gene_tables_names
list2env(overlap_cd_genes_ls, envir = .GlobalEnv)

names(overlap_nccd_genes_ls) <- nc_cd_gene_tables_names
list2env(overlap_nccd_genes_ls, envir = .GlobalEnv)



# Show message
message("Generating plots...")


### Prepare data for ggplot2 ###


# Use VST pairwise cor values for accuracy along TPM expression values for easy visualization
prepareCisNAT <- function(x, y) {

  df <- merge(x, y, by = c("gene_id", "prt_id", "biotype", "source", "info"))

  df <- subset(df, select = -c(Spearman.x, Pearson.x, 18:25))

  names(df) <- gsub(pattern = "\\..*", replacement = "", x = names(df))

  spec_name <- gsub(pattern = "\\_.*", replacement = "", x = deparse(substitute(x)))

  df$Species <- rep(spec_name)

  return(df)

}

AT_cd <- prepareCisNAT(AT_comparative_samples_cd_nc_cor_tpm, AT_comparative_samples_cd_nc_cor_count)
AL_cd <- prepareCisNAT(AL_comparative_samples_cd_nc_cor_tpm, AL_comparative_samples_cd_nc_cor_count)
CR_cd <- prepareCisNAT(CR_cd_nc_cor_tpm, CR_cd_nc_cor_count)
ES_cd <- prepareCisNAT(ES_cd_nc_cor_tpm, ES_cd_nc_cor_count)
TH_cd <- prepareCisNAT(TH_cd_nc_cor_tpm, TH_cd_nc_cor_count)
MT_cd <- prepareCisNAT(MT_cd_nc_cor_tpm, MT_cd_nc_cor_count) # check data!
BD_cd <- prepareCisNAT(BD_cd_nc_cor_tpm, BD_cd_nc_cor_count) # check data!


# Clean up data
MT_cd <- MT_cd[- grep("ALlnc", MT_cd$gene_id),] # rm one wrong id from list
BD_cd <- BD_cd[- grep("ALlnc", BD_cd$gene_id),] # rm one wrong id from list


# Add overlap length info to data
addLength <- function(x, y) {

  y1 <- y[y$gene_biotype1 != "protein_coding",]
  y1 <- data.frame(gene_id = y1$gene_id1, prt_id = y1$gene_id2, start = y1$start1, end = y1$end1, 
    strand = y1$strand1, width = y1$width1, overlap = y1$overlap)

  y2 <- y[y$gene_biotype2 != "protein_coding",]
  y2 <- data.frame(gene_id = y2$gene_id2, prt_id = y2$gene_id1, start = y2$start2, end = y2$end2, 
    strand = y2$strand2, width = y2$width2, overlap = y2$overlap)

  nc_overlap <- rbind(y1, y2)

  out_df <- merge(x, nc_overlap, by = c("gene_id", "prt_id"))

}

AT_cd <- addLength(AT_cd, AT_nd_cd_overlap)
AL_cd <- addLength(AL_cd, AL_nd_cd_overlap)
CR_cd <- addLength(CR_cd, CR_nd_cd_overlap)
ES_cd <- addLength(ES_cd, ES_nd_cd_overlap)
TH_cd <- addLength(TH_cd, TH_nd_cd_overlap)
MT_cd <- addLength(MT_cd, MT_nd_cd_overlap)
BD_cd <- addLength(BD_cd, BD_nd_cd_overlap)



all_spec_ls <- list(AT_cd, AL_cd, CR_cd, ES_cd, TH_cd, MT_cd, BD_cd)



  # Do correlation test
  testCor <- function(t) {

    p_val_NC <- cor.test(t$maxNC, t$Pearson, method = "pearson")$estimate
    t$maxNCPea <- rep(p_val_NC)

    p_val_Ratio <- cor.test(t$maxRatio, t$Pearson, method = "pearson")$estimate
    t$maxRatioPea <- rep(p_val_Ratio)

    p_val_Overlap <- cor.test(t$overlap, t$Pearson, method = "pearson")$estimate
    t$OverlapPea <- rep(p_val_Overlap)

    return(t)

  }

  all_spec_ls <- lapply(all_spec_ls, testCor)



# Define density plot coloursd
dcols <- colorRampPalette(c(
    "#3b4086", "#3b458e", "#3b4a95", "#3b4f9d",  "#3a54a5", "#316cb7", "#3083c5", 
    "#3c9ad1", "#53b0db","#5bbed8","#70cbd2","#75c9b5","#82ca97","#93cc79","#a7d059",
    "#bfd735","#bfd735", "#dae11e","#e7e71c","#f4ed1a", "#f8e410", "#f0d737" ,"#f3bf2f",
    "#f2a72f" ,"#f2871a","#f16214","#ed311c","#de2c1e","#cf2820","#c02420","#b12020",
    "#a21d20"))(256)


# Function to prepare data frame and encode data density as color for maxNC expression
scatterDensMAX <- function(x) {

  x$col250 <- densCols(x$Pearson, x$maxNC, nbin = 250, colramp = colorRampPalette(c("black", "white")))
  x$dens <- col2rgb(x$col250)[1,] + 1L

  x$col <- dcols[x$dens]
  
  # Reorder "plot_data" based on "dens" values - the highest density points are plotted on top
  x <- x[order(x$dens),]
  x <- na.omit(x)

  return(x)
}


all_cd_nc_cor_max <- do.call("rbind", lapply(all_spec_ls, scatterDensMAX))


# Function to prepare data frame and encode data density as color for NAT/PC expression ratio
scatterDensRATIO <- function(x) {

  x$col250 <- densCols(x$Pearson, x$maxRatio, nbin = 2000, colramp = colorRampPalette(c("black", "white")))
  x$dens <- col2rgb(x$col250)[1,] + 1L

  x$col <- dcols[x$dens]
  
  x <- x[order(x$dens),]
  x <- na.omit(x)

  return(x)
}

all_cd_nc_cor_ratio <- do.call("rbind", lapply(all_spec_ls, scatterDensRATIO))


# Function to prepare data frame and encode data density as color for NAT/PC expression ratio
scatterDensOverlap <- function(x) {

  x$col250 <- densCols(x$Pearson, x$overlap, nbin = 1000, colramp = colorRampPalette(c("black", "white")))
  x$dens <- col2rgb(x$col250)[1,] + 1L

  x$col <- dcols[x$dens]
  
  x <- x[order(x$dens),]
  x <- na.omit(x)

  return(x)
}

all_cd_nc_cor_overlap <- do.call("rbind", lapply(all_spec_ls, scatterDensOverlap))



# Generate plots
plotPC.NAT.feat <- function(data, feat_type) {

  yLabelsK = function(l) { 

    ifelse(l==0, 0, paste0(round(l/1e3,1),"K"))
  }

  if (feat_type == "maxNC") {

    p_title <- "Maximum NAT expression in relation to pairwise NAT/PC gene correlation"

    data$feat <- data$maxNC

    fname <- "cd_nc_cor_maxNC.pdf"

    y_lab <- "    Expression level (log2 TPM)"

    plt_mar <- c(0.5, 1.75, 1.5, 1.75)

    p_df1 <- data.frame(Species = unique(data$Species), 
      x = rep(-0.8, length(unique(data$Species))),
      y = rep(11.8, length(unique(data$Species))))

    p_df2 <- data.frame(Species = unique(data$Species), 
      x = rep(-0.53, length(unique(data$Species))),
      y = rep(12.1, length(unique(data$Species))), 
      label = c(round(unique(data$maxNCPea), digits = 2)))

    r_df <- data.frame(Species = unique(data$Species), 
      start = c(-1, -1, -1, -1, -1, -1, -1),
      end = c(-1, -1, -1, -1, -1, -1, -1))

    ytlmar <- margin(t = 0, r = 11.75, b = 0, l = 1)

    yrmin = -1
    yrmax = -1

  } else if (feat_type == "maxRatio") {

    p_title <- "NAT/PC maximum expression ratio in relation to pairwise NAT/PC gene correlation" 

    data$feat <- data$maxRatio

    fname <- "cd_nc_cor_maxRatio.pdf"

    y_lab <- " Expression ratio"

    plt_mar <- c(0.5, 1.75, 1.5, 1.75)

    p_df1 <- data.frame(Species = unique(data$Species), 
      x = rep(-0.8, length(unique(data$Species))),
      y = rep(44.0, length(unique(data$Species))))

    p_df2 <- data.frame(Species = unique(data$Species), 
      x = rep(-0.53, length(unique(data$Species))),
      y = rep(55.0, length(unique(data$Species))), 
      label = c(round(unique(data$maxRatioPea), digits = 2)))

    r_df <- data.frame(Species = unique(data$Species), 
      start = c(-1, -1, -1, -1, -1, -1, -1),
      end = c(-1, -1, -1, -1, -1, 0.2, -1))

    ytlmar <- margin(t = 0, r = 8, b = 0, l = 1)

    yrmin = 1.5
    yrmax = 1000

  } else if (feat_type == "overlap") {

    p_title <- "NAT/PC overlap length in relation to pairwise NAT/PC gene correlation" 

    data$feat <- data$overlap # Adjust column name!

    fname <- "cd_nc_cor_overlap_length.pdf"

    y_lab <- " Overlap length (bp)"

    plt_mar <- c(0.5, 1.75, 1.5, 1.75)

    p_df1 <- data.frame(Species = unique(data$Species), 
      x = rep(-0.8, length(unique(data$Species))),
      y = rep(3800, length(unique(data$Species))))

    p_df2 <- data.frame(Species = unique(data$Species), 
      x = rep(-0.53, length(unique(data$Species))),
      y = rep(3900, length(unique(data$Species))), 
      label = c(round(unique(data$OverlapPea), digits = 2)))

    r_df <- data.frame(Species = unique(data$Species), 
      start = c(-1, -1, -1, -1, -0.4, -1, -1),
      end = c(-0.055, -0.4, -0.055, -0.055, 0.2, -0.37, -0.05))

    ytlmar <- margin(t = 0, r = 9.25, b = 0, l = 1)

    yrmin = 3780
    yrmax = 4300 

  }

  data$Species <- factor(data$Species, levels = unique(data$Species))

  p <- ggplot(data, aes(x = Pearson, y = feat)) + 
  geom_point(size = 2.5, colour = data$col) + 
  geom_smooth(method = 'lm', formula = y ~ x, size = 2.5, col = "white") + 
  geom_rect(data = r_df, aes(NULL, NULL, xmin = start, xmax = end), 
    ymin = yrmin, ymax = yrmax , colour = "white", fill = "white", alpha = 1) + 
  geom_text(data = p_df1, mapping = aes(x = x, y = y, 
    label = as.character(expression(paste(rho, " = ")))
    ), size = 9.275, colour = "black", parse = TRUE, hjust = 0.325, vjust = 0) + 
  geom_text(data = p_df2, mapping = aes(x = x, y = y, label = label), size = 9.275, colour = "black", 
    parse = TRUE, hjust = 0, vjust = 0) + 
  scale_x_continuous(expand = c(0.05, 0), limits = c(-1, 1), 
    labels = function(x) sub(".0+$", "", x)) + 

  if (feat_type == "maxNC") {

    scale_y_continuous(expand = c(0, 0), limits = c(-0.77, 13.7)) 

  } else if (feat_type == "maxRatio") {

    scale_y_log10(expand = c(0, 0), limits = c(0.014, 150), breaks = c(0.1, 1, 10, 100), 
      labels = scales::trans_format("log10", scales::math_format(10^.x)))

  } else if (feat_type == "overlap") {

    scale_y_continuous(expand = c(0.05, 0), limits = c(-70, 4200), labels = yLabelsK)
  }

  q <- p + theme_classic() + xlab("Pearson's r") + ylab(y_lab) + ggtitle(p_title) + 
  theme(text = element_text(size = 23.5), 
    strip.text = element_text(size = 24.1, face = "plain"), 
    strip.text.x = element_text(margin = margin(0.4457, 0, 0.4457, 0, "cm")), 
    strip.background = element_rect(colour = 'black', fill = NA, size = 2.75), 
    axis.ticks.length = unit(0.25, "cm"), 
    axis.ticks = element_line(colour = "black", size = 1.4), 
    axis.line = element_line(colour = 'black', size = 1.4), 
    plot.margin = unit(plt_mar, "cm"), 
    axis.title.y = element_text(size = 25, margin = ytlmar, 
      colour = "black", face = "plain"), 
    axis.title.x = element_text(size = 25, margin = margin(t = 3.5, r = 0, b = 8.75, l = 0), 
      colour = "black", face = "plain"), 
    axis.text.x = element_text(size = 21.5, margin = margin(t = 4, b = 7.75), colour = "grey35", 
      angle = 0, vjust = 1, hjust = 0.5), 
    axis.text.y = element_text(size = 21.5, angle = 0, margin = margin(l = 0.75, r = 1.5), colour = "grey35"), 
    plot.title = element_text(size = 25.25, margin = margin(t = 5.0, b = 15.8), face = "plain"), 
    panel.spacing = unit(0.55, "cm"), 
    panel.grid.major = element_blank(),
    panel.grid.minor.x = element_blank(), 
    panel.grid.minor.y = element_blank(),  
    legend.position = "none")

  q <- q + facet_wrap(~ factor(Species, levels = c("AT", "AL", "CR", "ES", "TH", "MT", "BD")) , nrow = 1, scales = "free_x")

  ggsave(file = file.path(out_dir, "output", "plots", fname), plot = q, width = 28.5, 
    height = 6.5, units = c("in"))
}

suppressWarnings(plotPC.NAT.feat(data = all_cd_nc_cor_max, feat_type = "maxNC"))
suppressWarnings(plotPC.NAT.feat(data = all_cd_nc_cor_ratio, feat_type = "maxRatio"))
suppressWarnings(plotPC.NAT.feat(data = all_cd_nc_cor_overlap, feat_type = "overlap"))




# Generate density info scalebar
pdf(file = file.path(out_dir, "output", "plots", "density_scalebar.pdf"))

my.colors = colorRampPalette(dcols)
z = matrix(1:256, nrow = 1)
x = 1
y = seq(0, 256, len = 256) # Range of colour data values
image(x, y, z, col = my.colors(100), axes = FALSE, xlab = "", ylab = "")
axis(2)

dev.off()



