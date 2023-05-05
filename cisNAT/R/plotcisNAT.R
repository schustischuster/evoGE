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
    "#4a3191","#483b97","#47459d","#4650a4","#445aaa","#4365b1","#426fb7","#417abe",
    "#429dd6","#32c6f4","#61cbe6","#70cbd2","#75c9b5","#82ca97","#93cc79","#a7d059",
    "#bfd735","#bfd735", "#dae11e","#e7e71c","#f4ed1a", "#f8e410","#fcdb05","#fdb713","#f7951e","#f47321","#f05323",
    "#ee3523","#ed2224","#e61d25","#c52026","#c52026","#a21d20"))(256)


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
pdf(file = file.path(out_dir, "output", "plots", "density_scalebar"))

my.colors = colorRampPalette(dcols)
z = matrix(1:256, nrow = 1)
x = 1
y = seq(0, 256, len = 256) # Range of colour data values
image(x, y, z, col = my.colors(100), axes = FALSE, xlab = "", ylab = "")
axis(2)

dev.off()


































      #-- Analyse distribution of max expression for coding+non-coding genes across species --


      expr_table_ls_br <- expr_table_ls[c(1:4)]

      expr_table_repl_ls <- lapply(expr_table_ls_br, calculateAvgExpr)


      # Select comparative organs for AT and AL
      expr_table_repl_ls$AT_expr <- dplyr::select(expr_table_repl_ls$AT_expr, c("root_whole_root_5d", 
         "hypocotyl_10d", "leaf_12_7d", "apex_vegetative_7d", "apex_inflorescence_21d", 
         "flower_stg12_21d", "flower_stg12_stamens_21d", "flower_early_stg12_carpels_21d"))


      # Add dataset identifier to list elements
      spec_exp_names <- lapply(seq_along(expr_table_repl_ls), function(i) { 
         paste(names(expr_table_repl_ls)[[i]])
      })

      for (i in seq_along(expr_table_repl_ls)) {
         expr_table_repl_ls[[i]]$dataset <- rep(spec_exp_names[i], nrow(expr_table_repl_ls[[i]]))
      }


      
      # Get maximum expression values for coding and non-coding genes for Brassicaceae species

      getMaxExprDist <- function(df, scripttype = c("coding", "lncRNA"), c_level = c("all", "non-core" , "core")) {

         spec_id <- unique(sub("\\_.*", "", df$dataset))
         df <- within(df, rm(dataset))

         Core_expr <- Core_expr[!grepl("ERCC", Core_expr[,1]),]

         `%nin%` = Negate(`%in%`)


         # Get protein-coding and non-coding core IDs for non-AT species

         if (spec_id == "AT") {

            compl_table <- AT_expr_compl
            core_ids <- as.data.frame(sapply(Core_expr[,1], function(x) unlist(strsplit(x, "\\:"))[1]))
            core_lnc_ids <- as.data.frame(sapply(Brass_nc_expr[,1], function(x) unlist(strsplit(x, "\\:"))[1]))

         } else if (spec_id == "AL") {

            compl_table <- AL_expr_compl
            core_ids <- as.data.frame(sapply(Core_expr[,1], function(x) unlist(strsplit(x, "\\:"))[2]))
            core_lnc_ids <- as.data.frame(sapply(Brass_nc_expr[,1], function(x) unlist(strsplit(x, "\\:"))[2]))

         } else if (spec_id == "CR") {

            compl_table <- CR_expr_compl
            core_ids <- as.data.frame(sapply(Core_expr[,1], function(x) unlist(strsplit(x, "\\:"))[3]))
            core_lnc_ids <- as.data.frame(sapply(Brass_nc_expr[,1], function(x) unlist(strsplit(x, "\\:"))[3]))

         } else if (spec_id == "ES") {

            compl_table <- ES_expr_compl
            core_ids <- as.data.frame(sapply(Core_expr[,1], function(x) unlist(strsplit(x, "\\:"))[4]))
            core_lnc_ids <- as.data.frame(sapply(Brass_nc_expr[,1], function(x) unlist(strsplit(x, "\\:"))[4]))

         }

         # -------------------------- Process protein-coding data --------------------------

         if ((c_level == "all") && (scripttype == "coding")) {

            df <- df[rownames(df) %in% compl_table[compl_table$biotype == "protein_coding",]$id, ]


         # ------------------------------ Process lncRNA data ------------------------------

         } else if ((c_level == "all") && (scripttype == "lncRNA")) {

            df <- df[rownames(df) %in% subset(compl_table, subset = biotype %in% c(
               "lnc_exonic_antisense", "lnc_intronic_antisense", "lnc_intergenic"))$id, ]


         # ----------------- Process non-core ortholog protein-coding data -----------------

         } else if ((c_level == "non-core") && (scripttype == "coding")) {

            df <- df[rownames(df) %in% compl_table[compl_table$biotype == "protein_coding",]$id, ]
            df <- df[rownames(df) %nin% as.character(core_ids[,1]),]


         # --------------------- Process non-core ortholog lncRNA data ---------------------
         
         } else if ((c_level == "non-core") && (scripttype == "lncRNA")) {

            df <- df[rownames(df) %in% subset(compl_table, subset = biotype %in% c(
               "lnc_exonic_antisense", "lnc_intronic_antisense", "lnc_intergenic"))$id, ]
            df <- df[rownames(df) %nin% as.character(core_lnc_ids[,1]), ]


         # ------------------- Process core ortholog protein-coding data -------------------

         } else if ((c_level == "core") && (scripttype == "coding")) {

            df <- df[rownames(df) %in% as.character(core_ids[,1]),]


         # ----------------------- Process core ortholog lncRNA data -----------------------
         
         } else if ((c_level == "core") && (scripttype == "lncRNA")) {

            df <- df[rownames(df) %in% as.character(core_lnc_ids[,1]), ]
         }


         # Get max expression value per gene
         df$max <- apply(df, 1, max)

         df_out <- data.frame(
            species = rep(spec_id), 
            biotype = rep(scripttype), 
            conservation = rep(c_level), 
            class = rep(paste(scripttype, c_level, sep = "_")),
            max_expr = df$max
            )

         return(df_out)

      }

      cd_expr_dist_all <- do.call(rbind, lapply(expr_table_repl_ls, getMaxExprDist, scripttype = "coding", 
         c_level = "all"))
      cd_expr_dist_n_core <- do.call(rbind, lapply(expr_table_repl_ls, getMaxExprDist, scripttype = "coding", 
         c_level = "non-core"))
      cd_expr_dist_core <- do.call(rbind, lapply(expr_table_repl_ls, getMaxExprDist, scripttype = "coding", 
         c_level = "core"))
      nc_expr_dist_all <- do.call(rbind, lapply(expr_table_repl_ls, getMaxExprDist, scripttype = "lncRNA", 
         c_level = "all"))
      nc_expr_dist_n_core <- do.call(rbind, lapply(expr_table_repl_ls, getMaxExprDist, scripttype = "lncRNA", 
         c_level = "non-core"))
      nc_expr_dist_brass <- do.call(rbind, lapply(expr_table_repl_ls[1:4], getMaxExprDist, scripttype = "lncRNA", 
         c_level = "core")) # Ortholog lncRNA dataset is limited to Brassicaceae


      # Combine data
      max_expr_dist <- rbind(cd_expr_dist_all, cd_expr_dist_n_core, cd_expr_dist_core, nc_expr_dist_all, 
         nc_expr_dist_n_core, nc_expr_dist_brass)



      # ---------------------------- ggplot2 helper functions ----------------------------

      # horizontal nudge position adjustment
      position_hnudge <- function(x = 0) {
         ggproto(NULL, PositionHNudge, x = x)
      }

      PositionHNudge <- ggproto("PositionHNudge", Position,
         x = 0,
         required_aes = "x",
         setup_params = function(self, data) {
            list(x = self$x)
         },
         compute_layer = function(data, params, panel) {
            transform_position(data, function(x) x + params$x)
         }
      )


      # Function to create split violin plot
      "%||%" <- function(a, b) {
         if (!is.null(a)) a else b
      }

      geom_flat_violin <- function(mapping = NULL, data = NULL, stat = "ydensity",
                        position = "dodge", trim = TRUE, scale = "area",
                        show.legend = NA, inherit.aes = TRUE, ...) {

         layer(
            data = data,
            mapping = mapping,
            stat = stat,
            geom = GeomFlatViolin,
            position = position,
            show.legend = show.legend,
            inherit.aes = inherit.aes,
            params = list(
               trim = trim,
               scale = scale,
               ...
               )
            )
      }

      GeomFlatViolin <-
      ggproto("GeomFlatViolin", Geom,
         setup_data = function(data, params) {
            data$width <- data$width %||%
            params$width %||% (resolution(data$x, FALSE) * 0.9)

            # ymin, ymax, xmin, and xmax define the bounding rectangle for each group
            data %>%
            group_by(group) %>%
            mutate(ymin = min(y),
               ymax = max(y),
               xmin = x - width / 2,
               xmax = x)
         },

         draw_group = function(data, panel_scales, coord) {
         # Find the points for the line to go all the way around
            data <- transform(data, 
               xmaxv = x,
               xminv = x + violinwidth * (xmin - x))

            # Make sure it's sorted properly to draw the outline
            newdata <- rbind(plyr::arrange(transform(data, x = xminv), y),
               plyr::arrange(transform(data, x = xmaxv), -y))

            # Close the polygon: set first and last point the same
            # Needed for coord_polar and such
            newdata <- rbind(newdata, newdata[1,])

            ggplot2:::ggname("geom_flat_violin", GeomPolygon$draw_panel(newdata, panel_scales, coord))
         },

         draw_key = draw_key_polygon,

         default_aes = aes(weight = 1, colour = "grey20", fill = "white", size = 0.5,
            alpha = NA, linetype = "solid"),

         required_aes = c("x", "y")
      )



      # Wilcoxon rank sum test all genes vs all-ortho// all genes vs ortho
      getPMWU <- function(z, spec = c("AT", "AL", "CR", "ES")) {

         sp_cd_all <- z[z$species == spec & z$class == "coding_all",]
         sp_cd_ncore <- z[z$species == spec & z$class == "coding_non-core",]
         sp_cd_core <- z[z$species == spec & z$class == "coding_core",]

         sp_nc_all <- z[z$species == spec & z$class == "lncRNA_all",]
         sp_nc_ncore <- z[z$species == spec & z$class == "lncRNA_non-core",]
         sp_nc_core <- z[z$species == spec & z$class == "lncRNA_core",]

         cd_all_vs_core <- wilcox.test(sp_cd_all$max_expr, sp_cd_core$max_expr)$p.value
         cd_ncore_vs_core <- wilcox.test(sp_cd_ncore$max_expr, sp_cd_core$max_expr)$p.value

         nc_all_vs_core <- wilcox.test(sp_nc_all$max_expr, sp_nc_core$max_expr)$p.value
         nc_ncore_vs_core <- wilcox.test(sp_nc_ncore$max_expr, sp_nc_core$max_expr)$p.value

         pmwu <- data.frame(species = rep(spec), 
            comparison = c("cd_all_vs_core", "cd_ncore_vs_core", "nc_all_vs_core", "nc_ncore_vs_core"), 
            p_value = c(cd_all_vs_core, cd_ncore_vs_core, nc_all_vs_core, nc_ncore_vs_core))

         return(pmwu)
      }

      pmwu_at <- getPMWU(z = max_expr_dist, spec = "AT")
      pmwu_al <- getPMWU(z = max_expr_dist, spec = "AL")
      pmwu_cr <- getPMWU(z = max_expr_dist, spec = "CR")
      pmwu_es <- getPMWU(z = max_expr_dist, spec = "ES")

      p_mwu <- rbind(pmwu_at, pmwu_al, pmwu_cr, pmwu_es)


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



      # Split data AT/non-AT
      max_expr_dist_non_AT <- max_expr_dist[!grepl("AT", max_expr_dist$species),]
      max_expr_dist_AT <- max_expr_dist[max_expr_dist$species == "AT", ]



      # Generate plots
      plotMaxExprDist <- function(data, species) {

         if (species == "ACE") {

            p_mwu <- p_mwu[!grepl("AT", p_mwu$species),]

            # Create df for FDR p-value mapping
            mwu_df <- data.frame(
                class = rep(c("coding_non-core", "coding_core", "lncRNA_non-core", "lncRNA_core"), 
                    times = 3), 
                y = rep(c(19.68, 18.25), times = 6),
                label = ifelse(p_mwu$p_value < 1e-07, "****", 

                    c(paste("italic('P =')~", set_scientific(p_mwu$p_value)))), 

                species = rep(c("A.lyrata", "C.rubella", "E.salsugineum"), each = 4)
            )

            # Create df for gem_segments
            h_seg_df <- data.frame(
                x = rep(c(1.105, 2.105, 4.107, 5.107), times = 3), 
                xend = rep(c(3.105, 3.105, 6.107, 6.107), times = 3), 
                y = rep(c(20.08, 18.65, 20.08, 18.65), times = 3), 
                yend = rep(c(20.08, 18.65, 20.08, 18.65), times = 3), 
                species = rep(c("A.lyrata", "C.rubella", "E.salsugineum"), each = 4)
            )

            v_seg_df <- data.frame(
                x = rep(c(1.107, 3.107, 2.107, 3.107, 4.107, 6.107, 5.107, 6.107), times = 3), 
                xend = rep(c(1.107, 3.107, 2.107, 3.107, 4.107, 6.107, 5.107, 6.107), times = 3), 
                y = rep(c(19.64, 19.64, 18.2, 18.2, 19.64, 19.64, 18.2, 18.2), times = 3), 
                yend = rep(c(20.08, 20.08, 18.65, 18.65, 20.08, 20.08, 18.65, 18.65), times = 3), 
                species = rep(c("A.lyrata", "C.rubella", "E.salsugineum"), each = 4)
            )

            # Adjust position of p-value labels
            mwu_df$label <- paste0(mwu_df$label, c("", "              "))

            y_scale <- c(2.9, 21.125)

            plt_mar <- c(0.1, 1.55, 1.7, 0.55)

            stp_mar <- margin(0.25, 0, 0.25, 0, "cm")

         } else if (species == "AT") { 

          p_mwu <- p_mwu[!grepl("AT", p_mwu$species),]

            # Create df for FDR p-value mapping
            mwu_df <- data.frame(
                class = rep(c("coding_non-core", "coding_core", "lncRNA_non-core", "lncRNA_core"), 
                    times = 1), 
                y = rep(c(18.77, 17.6), times = 2),
                label = ifelse(p_mwu$p_value < 1e-07, "****", 

                    c(paste("italic('P =')~", set_scientific(p_mwu$p_value)))), 

                species = rep(c("A.thaliana"), each = 4)
            )

            # Create df for gem_segments
            h_seg_df <- data.frame(
                x = rep(c(1.105, 2.105, 4.107, 5.107), times = 1), 
                xend = rep(c(3.105, 3.105, 6.107, 6.107), times = 1), 
                y = rep(c(19.08, 17.925, 19.08, 17.925), times = 1), 
                yend = rep(c(19.08, 17.925, 19.08, 17.925), times = 1), 
                species = rep(c("A.thaliana"), each = 4)
            )

            v_seg_df <- data.frame(
                x = rep(c(1.107, 3.107, 2.107, 3.107, 4.107, 6.107, 5.107, 6.107), times = 1), 
                xend = rep(c(1.107, 3.107, 2.107, 3.107, 4.107, 6.107, 5.107, 6.107), times = 1), 
                y = rep(c(18.64, 18.64, 17.45, 17.45, 18.64, 18.64, 17.45, 17.45), times = 1), 
                yend = rep(c(19.08, 19.08, 17.9, 17.9, 19.08, 19.08, 17.9, 17.9), times = 1), 
                species = rep(c("A.thaliana"), each = 4)
            )

            # Adjust position of p-value labels
            mwu_df$label <- paste0(mwu_df$label, c("", "              "))

            y_scale <- c(5.5, 19.89)

            plt_mar <- c(0.1, 32.475, 1.7, 0.55)

            stp_mar <- margin(0.24, 0, 0.26, 0, "cm")

         }

         fname <- sprintf('%s.pdf', paste(deparse(substitute(data)), sep="_"))

         x_lab <- c(Root = "Rt", Hypocotyl = "Hc", Leaf = "Lf", Apex_veg = "Av", 
            Apex_inf = "Ai", Flower = "Fl", Stamen = "St", Carpel = "Ca")

         x_labels = c("coding_all" = expression(atop(NA, atop(textstyle('All'), textstyle('PC')))), 
            "coding_non-core" = expression(atop(NA, atop(textstyle('PC w/o'), textstyle('Ortho')))), 
            "coding_core" = expression(atop(NA, atop(textstyle('Ortho'), textstyle('PC')))), 
            "lncRNA_all" = expression(atop(NA, atop(textstyle('All'), textstyle('lnc')))), 
            "lncRNA_non-core" = expression(atop(NA, atop(textstyle('lnc w/o'), textstyle('Ortho')))), 
            "lncRNA_core" = expression(atop(NA, atop(textstyle('Ortho'), textstyle('lnc')))))

         data$species <- gsub("AT", "A.thaliana", data$species)
         data$species <- gsub("AL", "A.lyrata", data$species)
         data$species <- gsub("CR", "C.rubella", data$species)
         data$species <- gsub("ES", "E.salsugineum", data$species)

         data$class <- factor(data$class, levels = unique(data$class))
         data$conservation <- factor(data$conservation, levels = unique(data$conservation))
         data$species <- factor(data$species, levels = unique(data$species))

         p <- ggplot(data, aes(x = class, y = max_expr, color = class)) + geom_flat_violin(aes(fill = class), colour = "black", position = position_nudge(x = -0.037, y = 0), alpha = 1, size = 0.8) + 
         geom_boxplot(aes(fill = class), colour = "black", width = 0.44, outlier.shape = NA, position = position_hnudge(x = 0.25), size = 0.8, fatten = 2.8, notch = TRUE) +
         scale_x_discrete(expand = c(0.005, 0), labels = x_labels) + 
         scale_y_continuous(limits = y_scale, expand = c(0, 0), breaks = c(5,7.5,10,12.5,15,17.5))
         q <- p + 
         scale_fill_manual(values = c("coding_all" = "#f7ddb0", "coding_non-core" = "#edbb5c", 
            "coding_core" = "#e7a007", "lncRNA_all" = "#cdbee5", "lncRNA_non-core" = "#A689CE", 
            "lncRNA_core" = "#8055b8")) + 
         geom_text(data = mwu_df, mapping = aes(x = class, y = y, label = label), size = 9.275, colour = "black", 
            parse = FALSE, hjust = 0.325, vjust = 0) + 
         geom_segment(data = h_seg_df, mapping = aes(x = x, xend = xend, y = y, yend = yend), size = 0.8, colour = "black") + 
         geom_segment(data = v_seg_df, mapping = aes(x = x, xend = xend, y = y, yend = yend), size = 0.8, colour = "black") + 
         theme_classic() + 
         xlab("") + ylab("Maximum expression    \n (VST-normalized counts)     ") + ggtitle("") + 
         theme(text = element_text(size = 23.5), 
            strip.text = element_text(size = 19.5, face = "italic"), 
                strip.text.x = element_text(margin = stp_mar), 
                strip.background = element_rect(colour = 'white', fill = NA, size = 0.1), 
                axis.ticks.length = unit(0.2, "cm"), 
                axis.ticks = element_line(colour = "black", size = 0.95), 
                axis.line = element_line(colour = 'black', size = 0.95), 
                plot.margin = unit(plt_mar, "cm"), 
                axis.title.y = element_text(size = 18.4, margin = margin(t = 0, r = 6.4, b = 0, l = 3.38), 
                    colour = "black", face = "plain"), 
                axis.title.x = element_text(size = 18.75, margin = margin(t = 6.5, r = 0, b = 5.75, l = 0), 
                    colour = "black", face = "bold"), 
                axis.text.x = element_text(size = 16.25, margin = margin(t = -7, b = 2), colour = "black", 
                    angle = 0, vjust = 1, hjust = 0.5), 
                axis.text.y = element_text(size = 16.5, angle = 0, margin = margin(l = 0, r = 1.5), colour = "black"), 
                panel.spacing = unit(0.7, "cm"), 
                panel.grid.major = element_blank(),
                panel.grid.minor.x = element_blank(), 
                panel.grid.minor.y = element_blank(),  
                legend.position ="none")

         q <- q + facet_wrap(~ factor(species, levels = c("A.thaliana", "A.lyrata", "C.rubella", "E.salsugineum")) , nrow = 1, scales = "free_x")

            ggsave(file = file.path(out_dir, "output", "plots", fname), plot = q, 
                width = 20, height = 5.75, units = c("in"))
      }

      plotMaxExprDist(data = max_expr_dist_non_AT, species = "ACE")
      plotMaxExprDist(data = max_expr_dist_AT, species = "AT")


      # Wilcoxon rank sum test with continuity correction (all genes vs core genes)
      # W = 1729000, p-value < 2.2e-16 for comparison cd/lnc for all species




   }



}





