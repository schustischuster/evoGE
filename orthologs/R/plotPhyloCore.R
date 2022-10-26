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

    treename <- sprintf('%s.pdf', paste("Marker_species_phylogeny", "median_times" , sep="_"))

  } else if (is.element("Estimated", div_times)) {

    treename <- sprintf('%s.pdf', paste("Marker_species_phylogeny", "estimated_times" , sep="_"))
  }

     
  # Make marker species phylogeny plot
  pdf(file = file.path(out_dir, "output", "plots", treename), 
    width = 4.01, height = 4.5825)
  par(mar = c(0.88, 0.225, 0, 0.25), bg = NA)

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
  species <- rep(c(pairwise_ortholog$species[1:6]), times = 2)

  class <- rep(c("Protein-coding", "lncRNA"), each = 6)

  # divergence times are estimated taxon pair times from TimeTree
  # http://www.timetree.org/
  if(is.element("Median", div_times)) {

    div_time <- rep(c(5.8, 8.5, 27.4, 51, 108, 152), times=2)

  } else if (is.element("Estimated", div_times)) {

    div_time <- rep(c(7.1, 9.4, 25.6, 46, 106, 160), times=2)
  }

  div_time <- as.data.frame(div_time)
  colnames(div_time) <- "div_time"

  expressed <- data.frame(rbind(data.frame(expressed = pairwise_ortholog$protein_coding[1:6]), 
    data.frame(expressed = pairwise_ortholog$lncRNA[1:6])))

  ortholo_data <- data.frame(cbind(species, div_time, class, expressed))

  ortholo_data$species <- factor(ortholo_data$species)
  ortholo_data$class <- factor(ortholo_data$class)
  ortholo_data$div_time <- as.numeric(ortholo_data$div_time)
  ortholo_data$expressed <- as.numeric(ortholo_data$expressed)

  coding_lnc <- ortholo_data[1:12,]
  core_orthologs <- as.numeric(pairwise_ortholog[7,2:3])
  brass_orthologs <- as.numeric(pairwise_ortholog[8,2:3])



  # Plot number of deduplicated reads for each species
  plotOrthologs <- function(data) {

    fname <- sprintf('%s.pdf', paste(deparse(substitute(data)), div_times, sep="_"))

    if (is.element("Median", div_times)) {

      species_time <- c(6, 9, 27, 51, 108, 152)
      AL <- 6
      CR <- 9
      ES <- 27
      TH <- 51
      MT <- 108
      time_labels <- c(6, "", 27, 51, 108, 152)

    } else if (is.element("Estimated", div_times)) {

      species_time <- c(7, 9, 26, 46, 106, 160)
      AL <- 7
      CR <- 9
      ES <- 26
      TH <- 46
      MT <- 106
      time_labels <- c(7, "", 26, 46, 106, 160)
    }

    dat_circle <- data.frame(
      class   = c("Protein-coding", "lncRNA"),
      x       = c(156, 156),
      y       = c(core_orthologs[1], core_orthologs[2])
      )

    dat_text <- data.frame(
      label = c(core_orthologs[1], core_orthologs[2]),
      class = c("Protein-coding", "lncRNA"),
      x     = c(138.5, 147.25),
      y     = c(core_orthologs[1], core_orthologs[2] + 20)
      )

    y_labels <- function(l) { 
      ifelse(l==0, paste0(round(l/1e3,1)),paste0(round(l/1e3,1),"K"))
    }

    spec_label <- data.frame(
      label = c("Brass.", "TH", "MT", "BD",  "Brass.", "TH", "MT"),
      class   = c("Protein-coding", "Protein-coding", "Protein-coding", "Protein-coding", "lncRNA", "lncRNA", "lncRNA"),
      x     = c(18.5, 46, 106, 157, 18.5, 46, 106),
      y     = c(1200, 1200, 1200, 1200, 85.4, 85.4, 85.4)
      )
    p <- ggplot(transform(data, class=factor(class, levels=c("Protein-coding", "lncRNA"))), 
      aes(x = div_time, y = expressed)) + 
    geom_line(aes(x = div_time), size = 1.0) + expand_limits(y=0) + 
    scale_x_continuous(breaks = species_time, labels = time_labels) + 
    scale_y_continuous(expand = c(0.05, 0.15), breaks= pretty_breaks(), labels = y_labels) + 
    geom_text(aes(y = expressed * 1.077, label = "")) + 
        geom_hline(yintercept = 0, colour = "grey95", size = 0.5) + 
    geom_point(data = dat_circle, mapping = aes(x = x, y = y), shape = 21, colour = "black", 
      fill = "red", size = 3.0, stroke = 0.5) + 
    geom_text(data = dat_text, mapping = aes(x = x, y = y, label=label), color = "red", size = 3.8)

    q <- p + facet_wrap( ~ class, scales = 'free', ncol = 1) + 
    theme_bw() + xlab("Divergence time (Myr)") + ylab("Number of pairwise 1-1 orthologs w/ A.thaliana") + 
    scale_color_manual(values = "gray35") + 
    geom_text(data = spec_label, mapping = aes(x = x, y = y, label = label), size=3.75) + 
    theme(
      panel.background = element_rect(fill = "transparent"), 
      plot.background = element_rect(fill = "transparent", color = NA), 
      strip.background = element_rect(fill = "grey81",colour = "grey47", size=1.0), 
      strip.text.x = element_text(size=12, margin = margin(t = 0.1875,r = 0,b = 0.1525,l = 0,"cm")),
      plot.margin = unit(c(1, 2.085, 0.25, 0.295), "points"),
      axis.ticks.length = unit(0.15, "cm"),
      axis.ticks = element_line(colour = "gray15", size = 0.5), 
      panel.grid.major = element_line(size = 0.5, colour = "grey95"), 
      panel.grid.minor = element_blank(), 
      panel.spacing = unit(0.125, "lines"), 
      panel.grid.minor.y = element_line(size = 0.3, colour = "grey95"), 
      axis.title.x = element_text(colour = "black", size = 12, 
        margin = margin(t = 1.0, r = 0, b = 10.75, l = 0)), 
      axis.title.y = element_text(colour = "black", size = 12, 
        margin = margin(t = 0, r = 4.5, b = 0, l = 0)), 
      axis.text.x = element_text(colour = "black", size = 10.45, 
        margin = margin(t = 2, r = 0, b = 5.0, l = 0), vjust = 0.5), 
      axis.text.y = element_text(colour = "black", size = 10.45, margin = margin(t = 0, r = 0.75, b = 0, l = 0)), 
      panel.border = element_rect(colour = "grey47", fill = NA, size = 1.0))
    ggsave(file = file.path(out_dir, "output", "plots", fname), plot = q, 
      scale = 1, width = 3, height = 4.475, units = c("in"), 
      dpi = 800, limitsize = FALSE, bg = "transparent")
  }

  plotOrthologs(data = coding_lnc)

}



