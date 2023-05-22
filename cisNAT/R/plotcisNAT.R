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
in_dir_NOPC_pairs <- file.path(out_dir, "output", "cd_gene_pairs")


# Read table containing tandem duplicate genes in AT
# Data from Liu et al., GBE (2011)
tan_dupl <- read.table(file = file.path(in_dir, "AT_tandem_dupl", "AT_tandem_dupl.csv"), 
  sep = ";", dec = ".", header = TRUE, stringsAsFactors = FALSE)


# Read all csv files in input file path
readTable <- function(path, pattern = "*.csv") {
    files = list.files(path, pattern, full.names = TRUE)
    lapply(files, function(x) read.table(x, sep=";", dec=".", header=TRUE, stringsAsFactors=FALSE))
}

NAT_expr_cor_ls <- readTable(in_dir_NAT_cor)
overlap_cd_genes_ls <- readTable(in_dir_PC_pairs)
overlap_nccd_genes_ls <- readTable(in_dir_NCPC_pairs)
non_overlap_cdcd_genes_ls <- readTable(in_dir_NOPC_pairs)


# Get file names and save them in character vector
NAT_gene_tables_list <- as.character(list.files(in_dir_NAT_cor, pattern = "*.csv"))
NAT_gene_tables_names <- gsub('\\.csv$', '', NAT_gene_tables_list)

coding_gene_tables_list <- as.character(list.files(in_dir_PC_pairs, pattern = "*.csv"))
coding_gene_tables_names <- gsub('\\.csv$', '', coding_gene_tables_list)

nc_cd_gene_tables_list <- as.character(list.files(in_dir_NCPC_pairs, pattern = "*.csv"))
nc_cd_gene_tables_names <- gsub('\\.csv$', '', nc_cd_gene_tables_list)

nopc_gene_tables_list <- as.character(list.files(in_dir_NOPC_pairs, pattern = "*.csv"))
nopc_gene_tables_names <- gsub('\\.csv$', '', nopc_gene_tables_list)


# Change data frame names in list
names(NAT_expr_cor_ls) <- NAT_gene_tables_names
list2env(NAT_expr_cor_ls, envir = .GlobalEnv)

names(overlap_cd_genes_ls) <- coding_gene_tables_names
list2env(overlap_cd_genes_ls, envir = .GlobalEnv)

names(overlap_nccd_genes_ls) <- nc_cd_gene_tables_names
list2env(overlap_nccd_genes_ls, envir = .GlobalEnv)

names(non_overlap_cdcd_genes_ls) <- nopc_gene_tables_names
list2env(non_overlap_cdcd_genes_ls, envir = .GlobalEnv)



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





# Make correlation plots for PC/PC and NAT/PC gene pairs for all species

nc_pc_gene_ls <- list(AT_nc_pc = AT_comparative_samples_cd_nc_cor_count, AL_nc_pc = AL_comparative_samples_cd_nc_cor_count, 
  CR_nc_pc = CR_cd_nc_cor_count, ES_nc_pc = ES_cd_nc_cor_count, TH_nc_pc = TH_cd_nc_cor_count, 
  MT_nc_pc = MT_cd_nc_cor_count, BD_nc_pc = BD_cd_nc_cor_count)

pc_pc_gene_ls <- list(AT_pc_pc = AT_comp_cd_cd_SAS_cor_0.5, AL_pc_pc = AL_comp_cd_cd_SAS_cor_0.5, 
  CR_pc_pc = CR_cd_cd_SAS_cor_0.5, ES_pc_pc = ES_cd_cd_SAS_cor_0.5, TH_pc_pc = TH_cd_cd_SAS_cor_0.5, 
  MT_pc_pc = MT_cd_cd_SAS_cor_0.5, BD_pc_pc = BD_cd_cd_SAS_cor_0.5)


formatNcPc <- function(x) {

  df <- data.frame(species = rep(substr(x[1, 1], start = 1, stop = 2), nrow(x)), 
    gene_pair = rep("NAT/PC", nrow(x)), 
    gene_biotype1 = x$biotype, gene_biotype2 = rep("protein_coding", nrow(x)), 
    Spearman = x$Spearman, Pearson = x$Pearson)

  return(df)
}

nc_pc_cor_df <- do.call(rbind, lapply(nc_pc_gene_ls, formatNcPc))


formatPcPc <- function(y) {

  df <- data.frame(species = rep(substr(y[1, 2], start = 1, stop = 2), nrow(y)), 
    gene_pair = rep("PC/PC", nrow(y)), 
    gene_biotype1 = y$gene_biotype1, gene_biotype2 = y$gene_biotype2, 
    Spearman = y$Spearman, Pearson = y$Pearson)

  # Change species names to fit nomenclature in NAT/PC table
  df$species <- gsub("Ca", "CR", df$species)
  df$species <- gsub("Th", "ES", df$species)
  df$species <- gsub("10", "TH", df$species)
  df$species <- gsub("Me", "MT", df$species)
  df$species <- gsub("Br", "BD", df$species)

  return(df)
}

pc_pc_cor_df <- do.call(rbind, lapply(pc_pc_gene_ls, formatPcPc))

pc_pc_nc_pc_cor <- rbind(pc_pc_cor_df, nc_pc_cor_df)



# Do Wilcoxon rank sum test between pc/pc and nc/pc of same species
getMWU <- function(p) {

  nat_pc <- p[p$gene_pair == "NAT/PC",]
  pc_pc <- p[p$gene_pair == "PC/PC",]

  nat_pc_ls <- split(nat_pc, nat_pc$species)
  pc_pc_ls <- split(pc_pc, pc_pc$species)

  getPVal <- function(x, y){ 

    p_val <- wilcox.test(x$Pearson, y$Pearson)$p.value

    return(p_val)
  }

  p_df <- mapply(getPVal, nat_pc_ls, pc_pc_ls)

  mwu_df <- data.frame(species = names(p_df), p_value = p_df)

  return(mwu_df)

  }

p_val_df <- getMWU(pc_pc_nc_pc_cor)



# Define specific notation
set_scientific <- function(l) {
  # turn in to character string in scientific notation
  l <- formatC(l, format = "e", digits = 0)
  # quote the part before the exponent to keep all the digits
  l <- gsub("^(.*)e", "'\\1'e", l)
  # turn the 'e+' into plotmath format
  l <- gsub("e", "%*%10^", l)
  # return this as an expression
  parse(text=l)
}



# Generate plots
plotCdNcCor <- function(data) {

   # Create df for FDR p-value mapping
   mwu_df <- data.frame(
    species = p_val_df$species,
    p_val = p_val_df$p_value, 
    y = rep(c(1.155), 7),
    label = ifelse(p_val_df$p_value < 1e-07, "**** ", 

      c(paste("italic('P =')~", set_scientific(p_val_df$p_value))))
    )
   # Create df for gem_segments
   h_seg_df <- data.frame(
    x = c(0.75, 1.75, 2.75, 3.75, 4.75, 5.75, 6.75), 
    xend = c(1.25, 2.25, 3.25, 4.25, 5.25, 6.25, 7.25), 
    y = rep(1.19, 7), 
    yend = rep(1.19, 7) 
    )

   v_seg_df <- data.frame(
    x = c(0.75, 1.25, 1.75, 2.25, 2.75, 3.25, 3.75, 4.25, 4.75, 5.25, 5.75, 6.25, 6.75, 7.25), 
    xend = c(0.75, 1.25, 1.75, 2.25, 2.75, 3.25, 3.75, 4.25, 4.75, 5.25, 5.75, 6.25, 6.75, 7.25), 
    y = c(1.12, 1.12, 1.12, 1.12, 1.12, 1.12, 1.12, 1.12, 1.12, 1.12, 1.12, 1.12, 1.12, 1.12), 
    yend = c(1.19, 1.19, 1.19, 1.19, 1.19, 1.19, 1.19, 1.19, 1.19, 1.19, 1.19, 1.19, 1.19, 1.19)
    )

   # Adjust position of p-value labels
   mwu_df$label <- paste0(mwu_df$label, c("", "              "))

   fname <- sprintf('%s.pdf', paste(deparse(substitute(data)), sep="_"))
   data$gene_pair <- factor(data$gene_pair, levels = unique(data$gene_pair))
   data$species <- factor(data$species, levels = unique(data$species))

   p <- ggplot(data, aes(x = species, y = Pearson, color = gene_pair)) + 
   geom_boxplot(aes(fill = gene_pair), colour = "black", width = 0.65, outlier.shape = NA, 
    size = 0.8, fatten = 2.8, notch = TRUE, position = position_dodge(width = 0.83), show.legend = FALSE) + 
   geom_point(size = -1, ) + 
   scale_x_discrete(expand = c(0.025, 0)) + 
   scale_y_continuous(limits = c(-1.05, 1.5), expand = c(0, 0), breaks = c(-1, -0.5, 0, 0.5, 1))

   q <- p + 
   scale_fill_manual(values = c("PC/PC" = "#f7ddb0", "NAT/PC" = "#cdbee5")) + 
   scale_colour_manual(values = c("PC/PC" = "#f7ddb0", "NAT/PC" = "#cdbee5")) + 
   geom_text(data = mwu_df, mapping = aes(x = c(0.775, 1.925, 2.775, 3.925, 4.775, 5.925, 6.775), y = y, label = label), 
    size = 8.9, colour = "black", parse = FALSE, hjust = 0.1, vjust = 0) + 
   geom_segment(data = h_seg_df, mapping = aes(x = x, xend = xend, y = y, yend = yend), 
    size = 0.8, colour = "black") + 
   geom_segment(data = v_seg_df, mapping = aes(x = x, xend = xend, y = y, yend = yend), 
    size = 0.8, colour = "black") + guides(colour = guide_legend(override.aes = list(size = 7, shape = 15))) + 
   theme_classic() + 
   xlab("Species") + ylab("Pearson's r     ") + ggtitle("") + labs(colour = 'Gene pair') + 
   theme(text = element_text(size = 23.5),  
        axis.ticks.length = unit(0.2, "cm"), 
        axis.ticks = element_line(colour = "black", size = 0.95), 
        axis.line = element_line(colour = 'black', size = 0.95), 
        plot.margin = unit(c(0.25, 32.14, 1.7, 0.1), "cm"), 
        axis.title.y = element_text(size = 18.4, margin = margin(t = 0, r = 4, b = 0, l = 1), 
          colour = "black", face = "plain"), 
        axis.title.x = element_text(size = 18.4, margin = margin(t = 2.8, r = 0, b = 23.525, l = 0), 
          colour = "black", face = "plain"), 
        axis.text.x = element_text(size = 16.25, margin = margin(t = 3.5, b = 2.0), colour = "black", 
          angle = 0, vjust = 1, hjust = 0.5), 
        axis.text.y = element_text(size = 16.5, angle = 0, margin = margin(l = 0, r = 1.5), colour = "black"), 
        panel.grid.major = element_blank(),
        panel.grid.minor.x = element_blank(), 
        panel.grid.minor.y = element_blank(),  
        legend.position = "top", 
        legend.title = element_text(size = 18.4, face = "bold"), 
        legend.text = element_text(size = 18.4), 
        legend.margin = margin(t = 4, b = -7))

   ggsave(file = file.path(out_dir, "output", "plots", fname), plot = q, 
    width = 20, height = 5.75, units = c("in"))
 }

 plotCdNcCor(data = pc_pc_nc_pc_cor)




#-------------------------- Plots AT coding and non-coding gene data --------------------------


# Prepare data
# AT_cd_cd_NO_cor_0.5 -> PC/PC non-overlapping gene pair correlation all AT samples (10932)
# AT_cd_cd_SAS_cor_0.5 -> PC/PC overlapping gene pair correlation all AT samples (4183)
# AT_cd_nc_cor_count -> cisNAT/PC overlapping gene pair correlation all AT samples (3175)

# Remove 466 tandem duplicate gene pairs from neighbouring protein-coding gene list in AT
AT_cd_cd_NO_cor_0.5 <- filter(AT_cd_cd_NO_cor_0.5, !(gene_id1 %in% unlist(tan_dupl)))
AT_cd_cd_NO_cor_0.5 <- filter(AT_cd_cd_NO_cor_0.5, !(gene_id2 %in% unlist(tan_dupl)))


# Get gene pairs where genes are located on same strand
sstr1 <- AT_cd_cd_NO_cor_0.5[AT_cd_cd_NO_cor_0.5$strand1 == "+",]
sstr1 <- sstr1[sstr1$strand2 == "+",]
sstr2 <- AT_cd_cd_NO_cor_0.5[AT_cd_cd_NO_cor_0.5$strand1 == "-",]
sstr2 <- sstr2[sstr2$strand2 == "-",]
sstr <- rbind(sstr1, sstr2)

# Get gene pairs where genes are located on opposite strand
osst1 <- AT_cd_cd_NO_cor_0.5[AT_cd_cd_NO_cor_0.5$strand1 == "+",]
osst1 <- osst1[osst1$strand2 == "-",]
osst2 <- AT_cd_cd_NO_cor_0.5[AT_cd_cd_NO_cor_0.5$strand1 == "-",]
osst2 <- osst2[osst2$strand2 == "+",]
ostr <- rbind(osst1, osst2)



# Prepare data for ggpot2
row_num <- nrow(AT_cd_cd_NO_cor_0.5) + nrow(AT_cd_cd_SAS_cor_0.5) + nrow(AT_cd_nc_cor_count)

at_all_samples_cd_nc_df <- data.frame(

  Species = rep("AT", row_num),

  Feature = c(rep("PCSS", nrow(sstr)), 
              rep("PCOS", nrow(ostr)), 
              rep("PC/PC", nrow(AT_cd_cd_SAS_cor_0.5)), 
              rep("NAT/PC", nrow(AT_cd_nc_cor_count))),

  Pearson = c(sstr$Pearson, 
              ostr$Pearson, 
              AT_cd_cd_SAS_cor_0.5$Pearson, 
              AT_cd_nc_cor_count$Pearson), 

  Spearman = c(sstr$Spearman, 
               ostr$Spearman, 
               AT_cd_cd_SAS_cor_0.5$Spearman, 
               AT_cd_nc_cor_count$Spearman)
  )



# Do Wilcoxon rank sum test between pc/pc and nc/pc of same species
getWRS <- function(z) {

  pc_ss <- z[z$Feature == "PCSS",]
  pc_os <- z[z$Feature == "PCOS",]
  pc_pc <- z[z$Feature == "PC/PC",]
  nat_pc <- z[z$Feature == "NAT/PC",]

  getPVal <- function(x, y){ 

    p_val <- wilcox.test(x$Pearson, y$Pearson)$p.value

    return(p_val)
  }

  p_1 <- getPVal(pc_ss, nat_pc)
  p_2 <- getPVal(pc_os, nat_pc)
  p_3 <- getPVal(pc_pc, nat_pc)

  mwu_df <- data.frame(comp = c("PCSS_NAT", "PCOS_NAT", "PCPC_NAT"), 
                       p_value = c(p_1, p_2, p_3))

  return(mwu_df)

  }

p_val_table <- getWRS(at_all_samples_cd_nc_df)



# Generate plots
plotATCor <- function(data) {


  # Create df for FDR p-value mapping
  mwu_df <- data.frame(
    p_val = p_val_table$p_value, 
    x = c(1, 2, 3),
    y = rep(c(1.155), 3),
    label = ifelse(p_val_table$p_value < 1e-07, "**** ", 

      c(paste("italic('P =')~", set_scientific(p_val_table$p_value))))
  )

  # Create df for gem_segments
  h_seg_df <- data.frame(
    x = c(0.75, 1.75, 2.75), 
    xend = c(1.25, 2.25, 3.25), 
    y = rep(1.19, 3), 
    yend = rep(1.19, 3)
  )

  v_seg_df <- data.frame(
    x = c(0.75, 1.25, 1.75, 2.25, 2.75, 3.25), 
    xend = c(0.75, 1.25, 1.75, 2.25, 2.75, 3.25), 
    y = c(1.12, 1.12, 1.12, 1.12, 1.12, 1.12), 
    yend = c(1.19, 1.19, 1.19, 1.19, 1.19, 1.19)
  )

   # Adjust position of p-value labels
   mwu_df$label <- paste0(mwu_df$label, c("", "              ", ""))

   fname <- sprintf('%s.pdf', paste(deparse(substitute(data)), sep="_"))
   data$Feature <- factor(data$Feature, levels = unique(data$Feature))
   data$Species <- factor(data$Species, levels = unique(data$Species))

   p <- ggplot(data, aes(x = Feature, y = Pearson, color = Feature)) + 
   geom_boxplot(aes(fill = Feature), colour = "black", width = 0.65, outlier.shape = NA, 
    size = 0.8, fatten = 2.8, notch = TRUE, position = position_dodge(width = 0.83), show.legend = FALSE) + 
   geom_point(size = -1, ) + 
   scale_x_discrete(expand = c(0.025, 0)) + 
   scale_y_continuous(limits = c(-1.05, 1.5), expand = c(0, 0), breaks = c(-1, -0.5, 0, 0.5, 1))

   q <- p + 
   scale_fill_manual(values = c("PCSS" = "white", "PCOS" = "white", "PC/PC" = "#f7ddb0", "NAT/PC" = "#cdbee5")) + 
   scale_colour_manual(values = c("PCSS" = "white", "PCOS" = "white", "PC/PC" = "#f7ddb0", "NAT/PC" = "#cdbee5")) + 
   geom_text(data = mwu_df, mapping = aes(x = x, y = y, label = label), 
    size = 8.9, colour = "black", parse = FALSE, hjust = 0.1, vjust = 0) + 
   geom_segment(data = h_seg_df, mapping = aes(x = x, xend = xend, y = y, yend = yend), 
    size = 0.8, colour = "black") + 
   geom_segment(data = v_seg_df, mapping = aes(x = x, xend = xend, y = y, yend = yend), 
    size = 0.8, colour = "black") + guides(colour = guide_legend(override.aes = list(size = 7, shape = 15))) + 
   theme_classic() + 
   xlab("Species") + ylab("Pearson's r     ") + ggtitle("") + labs(colour = 'Gene pair') + 
   theme(text = element_text(size = 23.5),  
        axis.ticks.length = unit(0.2, "cm"), 
        axis.ticks = element_line(colour = "black", size = 0.95), 
        axis.line = element_line(colour = 'black', size = 0.95), 
        plot.margin = unit(c(0.25, 32.14, 1.7, 0.1), "cm"), 
        axis.title.y = element_text(size = 18.4, margin = margin(t = 0, r = 4, b = 0, l = 1), 
          colour = "black", face = "plain"), 
        axis.title.x = element_text(size = 18.4, margin = margin(t = 2.8, r = 0, b = 23.525, l = 0), 
          colour = "black", face = "plain"), 
        axis.text.x = element_text(size = 16.25, margin = margin(t = 3.5, b = 2.0), colour = "black", 
          angle = 0, vjust = 1, hjust = 0.5), 
        axis.text.y = element_text(size = 16.5, angle = 0, margin = margin(l = 0, r = 1.5), colour = "black"), 
        panel.grid.major = element_blank(),
        panel.grid.minor.x = element_blank(), 
        panel.grid.minor.y = element_blank(),  
        legend.position = "top", 
        legend.title = element_text(size = 18.4, face = "bold"), 
        legend.text = element_text(size = 18.4), 
        legend.margin = margin(t = 4, b = -7))

   ggsave(file = file.path(out_dir, "output", "plots", fname), plot = q, 
    width = 20, height = 5.75, units = c("in"))
 }

 plotATCor(data = at_all_samples_cd_nc_df)



