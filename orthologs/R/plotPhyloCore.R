# Prepare Brawand and DevSeq phylogenetic tree and ortholog plot
# Data input: Pairwise divergence times and orthologr ortholog gene output
# Pairwise TimeTree divergence were obtained from TimeTree (http://www.timetree.org/)


plotPhyloCore <- function(div_times = c("Median", "Estimated")) {

    
#--------------------------- Make phylogram of DevSeq marker species ---------------------------


  if (is.element("Median", div_times)) {

    devseq_phylo <- "((((((Arabidopsis_thaliana:5.8, Arabidopsis_lyrata:5.8):2.7,Capsella_rubella:8.5):18.9,Eutrema_salsugineum:27.4):23.6,Tarenaya_hassleriana:51):57,Medicago_truncatula:108):44,Brachypodium_distachyon:152);"


  } else if (is.element("Estimated", div_times)) {

    devseq_phylo <- "((((((Arabidopsis_thaliana:7.1, Arabidopsis_lyrata:7.1):2.3,Capsella_rubella:9.4):16.2,Eutrema_salsugineum:25.6):20.4,Tarenaya_hassleriana:46):60,Medicago_truncatula:106):54,Brachypodium_distachyon:160);"
  }


  vert_tree <- read.tree(text = devseq_phylo)



  # Create "plots" folder in /out_dir/output/plots
  if (!dir.exists(file.path(out_dir, "output", "plots"))) 
  dir.create(file.path(out_dir, "output", "plots"), recursive = TRUE)


  if (is.element("Median", div_times)) {

    treename <- sprintf('%s.png', paste("Marker_species_phylogeny", "median_times" , sep="_"))

  } else if (is.element("Estimated", div_times)) {

    treename <- sprintf('%s.png', paste("Marker_species_phylogeny", "estimated_times" , sep="_"))
  }

     
  # Make marker species phylogeny plot
  png(file = file.path(out_dir, "output", "plots", treename), width = 3250, height = 3700, 
    res = 800)
  par(mar = c(0.5, 0.25, 0, 0.25), bg=NA)

  plot(vert_tree, type = "phylogram", use.edge.length = TRUE, show.node.label = FALSE, 
    edge.width = 1.55, edge.lty = 1, font = 3, root.edge = FALSE, label.offset = 2, 
    direction = "upwards")

  add.scale.bar(x = 1, y = 1, cex = 0.8, font = 5, col = "red")

  dev.off()




#---------------------- Plot pairwise orthologs for DevSeq marker species ----------------------


  # Read coding and non-coding pairwise_orthologs ATH vs species X
  pairwise_ortholog <- read.table(file=file.path(in_dir, "summarized_pairwise_orthologs.csv"), 
    sep=";", dec=".", header=TRUE, stringsAsFactors = FALSE)


  # Reshape data for ggplot2
  species <- rep(c(pairwise_ortholog$species[1:6]), times = 3)

  class <- rep(c("coding", "lncRNA", "circRNA"), each = 6)

  # divergence times are estimated taxon pair times from TimeTree
  # http://www.timetree.org/
  if (is.element("Median", div_times)) {

    div_time <- rep(c(5.8, 8.5, 27.4, 51, 108, 152), times=3)

  } else if (is.element("Estimated", div_times)) {

    div_time <- rep(c(7.1, 9.4, 25.6, 46, 106, 160), times=3)
  }

  div_time <- as.data.frame(div_time)
  colnames(div_time) <- "div_time"

  expressed <- data.frame(rbind(data.frame(expressed = pairwise_ortholog$protein_coding[1:6]), 
  data.frame(expressed = pairwise_ortholog$lncRNA[1:6]), data.frame(expressed = pairwise_ortholog$circRNA[1:6])))

  ortholo_data <- data.frame(cbind(species, div_time, class, expressed))

  ortholo_data$species <- factor(ortholo_data$species)
  ortholo_data$class <- factor(ortholo_data$class)
  ortholo_data$div_time <- as.numeric(ortholo_data$div_time)
  ortholo_data$expressed <- as.numeric(ortholo_data$expressed)

  coding_lnc <- ortholo_data[1:12,]
  circ_lnc <- ortholo_data[7:18,]
  core_orthologs <- as.numeric(pairwise_ortholog[7,2:4])



  # Plot number of deduplicated reads for each species
  plotOrthologs <- function(data) {

    fname <- sprintf('%s.png', paste(deparse(substitute(data)), div_times, sep="_"))

    if (is.element("Median", div_times)) {

      species_time <- c(6, 9, 27, 51, 108, 152)
      AL <- 6
      CR <- 9
      ES <- 27
      TH <- 51
      MT <- 108

    } else if (is.element("Estimated", div_times)) {

      species_time <- c(7, 9, 26, 46, 106, 160)
      AL <- 7
      CR <- 9
      ES <- 26
      TH <- 46
      MT <- 106
    }

    if (deparse(substitute(data)) == "coding_lnc") { 

      dat_circle <- data.frame(
        class   = c("coding", "lncRNA"),
        x     = c(156, 156),
        y     = c(core_orthologs[1], core_orthologs[2])
      )

      dat_text <- data.frame(
          label = c(core_orthologs[1], core_orthologs[2]),
          class = c("coding", "lncRNA"),
          x     = c(138.1, 146.7),
          y     = c(core_orthologs[1], core_orthologs[2] + 20)
      )

      y_labels <- function(l) { 
        ifelse(l==0, paste0(round(l/1e3,1)),paste0(round(l/1e3,1),"K"))
      }

      title_y_rmg <- 1

    } else if (deparse(substitute(data)) == "circ_lnc") {

      dat_circle <- data.frame(
          class   = c("circRNA", "lncRNA"),
          x       = c(155.7, 155.7),
          y       = c(core_orthologs[3], core_orthologs[2])
      )

      dat_text <- data.frame(
          label = c(core_orthologs[3], core_orthologs[2]),
          class = c("circRNA", "lncRNA"),
          x     = c(146.4, 146.4),
          y     = c(core_orthologs[3] + 5, core_orthologs[2] + 20)
      )

      y_labels <- waiver()
      title_y_rmg <- 0.4

    }

    p <- ggplot(data, aes(x = div_time, y = expressed)) + 
        geom_line(aes(x = div_time), size = 1.0) + expand_limits(y=0) + 
        scale_x_continuous(breaks = species_time) + 
        scale_y_continuous(expand = c(0.05, 0.15), breaks= pretty_breaks(), labels = y_labels) + 
        geom_text(aes(y = expressed * 1.077, label = "")) + 
        geom_segment(x = AL, xend = AL, y = -1000, yend = 0, color = "gray15", size = 0.45) + 
        geom_segment(x = CR, xend = CR, y = -1000, yend = 0, color = "gray15", size = 0.45) + 
        geom_segment(x = ES, xend = ES, y = -1000, yend = 0, color = "gray15", size = 0.45) + 
        geom_segment(x = TH, xend = TH, y = -1000, yend = 0, color = "gray15", size = 0.45) + 
        geom_segment(x = MT, xend = MT, y = -1000, yend = 0, color = "gray15", size = 0.45) + 
        geom_hline(yintercept = 0, colour = "grey95", size = 0.5) + 
        geom_point(data = dat_circle, mapping = aes(x = x, y = y), shape = 21, colour = "black", 
          fill = "red", size = 3.25, stroke = 0.75) + 
        geom_text(data = dat_text, mapping = aes(x = x, y = y, label=label), color = "red", size = 4)

    q <- p + facet_wrap( ~ class, scales='free_y', ncol = 2) + 
        theme_bw() + xlab("Divergence time (Myr)") + ylab("Number of orthologs AT vs X") + 
        scale_color_manual(values = "gray35") + 
        theme(
          panel.background = element_rect(fill = "transparent"), 
          plot.background = element_rect(fill = "transparent", color = NA), 
          strip.text = element_blank(), 
          strip.background = element_blank(),
          plot.margin = unit(c(15, 0, 0, 2), "points"),
          axis.ticks.length = unit(0.15, "cm"),
          axis.ticks = element_line(colour = "gray15", size = 0.45), 
          panel.grid.major = element_line(size = 0.5, colour = "grey95"), 
          panel.grid.minor = element_blank(), 
          panel.grid.minor.y = element_line(size = 0.3, colour = "grey95"), 
          axis.title.x = element_text(colour = "black", size=12, 
            margin = margin(t = 4.75, r = 0, b = -14, l = 0)), 
          axis.title.y = element_text(colour = "black", size=12, 
            margin = margin(t = 0, r = title_y_rmg, b = 0, l = 0.125)), 
          axis.text.x = element_text(colour = "black", size=11.5, 
            margin = margin(t = 2, r = 0, b = 12.5, l = 0), vjust = 0.5), 
          axis.text.y = element_text(colour = "black", size=11.5, margin = margin(t = 0, r = 0.75, b = 0, l = 0)), 
          panel.border = element_rect(colour = "grey70", fill=NA, size=1))


    pg <- ggplot(data, aes(x = div_time, y = expressed)) + 
        geom_line(aes(x = div_time), size=0.55) + expand_limits(y=0) + 
        scale_x_continuous(breaks = species_time) + 
        scale_y_continuous(expand = c(0.05, 0.15), breaks= pretty_breaks(), labels = y_labels) + 
        geom_text(aes(y = expressed * 1.077, label = ""))

        qg <- pg + facet_grid( ~ class, scales='free') + 
        theme_bw() + xlab("Divergence time (Myr)") + ylab("Number of orthologs AT vs X") + 
        scale_color_manual(values = "gray35") + 
        theme(
          panel.background = element_rect(fill = "transparent"), 
          plot.background = element_rect(fill = "transparent", color = NA), 
          strip.text.x = element_text(margin = margin(0.2125, 0, 0.2125, 0, "cm"), size = 12), 
          strip.text.y = element_text(margin = margin(0, 0.2125, 0, 0.2125, "cm"), size = 12),
          strip.background = element_rect(colour="grey70", size=1), 
          plot.margin = unit(c(15, 0, 0, 2), "points"),
          axis.ticks.length = unit(0.15, "cm"),
          panel.grid.major = element_line(size = 0.5, colour = "grey95"), 
          panel.grid.minor = element_blank(), 
          panel.grid.minor.y = element_line(size = 0.3, colour = "grey95"), 
          axis.ticks = element_line(colour = "gray15", size = 0.4), 
          axis.title.x = element_text(colour = "black", size=12, 
            margin = margin(t = 4.75, r = 0, b = -14, l = 0)), 
          axis.title.y = element_text(colour = "black", size=12, 
            margin = margin(t = 0, r = title_y_rmg, b = 0, l = 0.125)), 
          axis.text.x = element_text(colour = "black", size=11.5, 
            margin = margin(t = 2, r = 0, b = 12.5, l = 0), vjust = 0.5), 
          axis.text.y = element_text(colour = "black", size=11.5, margin = margin(t = 0, r = 0.75, b = 0, l = 0)),  
          panel.border = element_rect(colour = "grey70", fill=NA, size=1))


    gt1 = ggplot_gtable(ggplot_build(q))
    gt2 = ggplot_gtable(ggplot_build(qg))

    gt1 <- gtable_add_rows(gt1, heights = unit(0.325, 'cm'), pos = 2)
    gt1 <- gtable_add_grob(gt1, grobs = gt2$grobs[grep('strip-t', gt2$layout$name)], t = 2, l = gt1$layout[grep('strip-t.+1$', gt1$layout$name),]$l)

    gt.side1 = gtable_filter(gt2, 'strip-t-1')
    gt.side2 = gtable_filter(gt2, 'strip-t-2')

    gt1 = gtable_add_cols(gt1, widths = unit(0.01, 'cm'), pos = -1)
    gt1 = gtable_add_grob(gt1, zeroGrob(), t = 1, l = ncol(gt1), b=nrow(gt1))

    panel_id <- gt1$layout[grep('panel-.+1$', gt1$layout$name),]
    gt1 = gtable_add_grob(gt1, gt.side1, t = panel_id$t[1], l = ncol(gt1))
    gt1 = gtable_add_grob(gt1, gt.side2, t = panel_id$t[2], l = ncol(gt1))


    ggsave(file = file.path(out_dir, "output", "plots", fname), plot = gt1, 
      scale = 1, width = 6.2, height = 2.5, units = c("in"), 
      dpi = 800, limitsize = FALSE, bg = "transparent")
  }

  plotOrthologs(data = coding_lnc)
  plotOrthologs(data = circ_lnc)

}


plotPhyloCore("Estimated")
plotPhyloCore("Median")

