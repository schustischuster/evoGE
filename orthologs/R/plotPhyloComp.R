# Prepare Brawand and DevSeq simplified phylogenetic trees
# Data input: Pairwise divergence times that were obtained from TimeTree (http://www.timetree.org/)


plotPhyloComp <- function(div_times = c("Median", "Estimated")) {

    
#--------------------------- Make phylogram of DevSeq marker species ---------------------------


  if (is.element("Median", div_times)) {

    devseq_phylo <- "((((((A.thaliana:5.8, A.lyrata:5.8):2.7,C.rubella:8.5):18.9,E.salsugineum:27.4):23.6,T.hassleriana:51):57,M.truncatula:108):44,B.distachyon:152);"
    tetra_phylo <- "((((((P.troglodytes:6.4, H.sapiens:6.4):2.2,G.gorilla:8.6):6.6,P.pygmaeus:15.2):13.6,M.mulatta:28.8):60.2,M.musculus:89):71,M.domestica:160);"


  } else if (is.element("Estimated", div_times)) {

    devseq_phylo <- "((((((A.thaliana:7.1, A.lyrata:7.1):2.3,C.rubella:9.4):16.2,E.salsugineum:25.6):20.4,T.hassleriana:46):60,M.truncatula:106):54,B.distachyon:160);"
    tetra_phylo <- "((((((P.troglodytes:6.7,H.sapiens:6.7):2.4,G.gorilla:9.1):6.7,P.pygmaeus:15.8):13.6,M.mulatta:29.4):60.6,M.musculus:90):69,M.domestica:159);"
  }


  devseq_tree <- read.tree(text = devseq_phylo)
  tetra_tree <- read.tree(text = tetra_phylo)



  # Create "plots" folder in /out_dir/output/plots
  if (!dir.exists(file.path(out_dir, "output", "plots"))) 
  dir.create(file.path(out_dir, "output", "plots"), recursive = TRUE)



  plotTree <- function(t) {

    # Make marker species phylogeny plot
    fname <- sprintf('%s.png', paste(deparse(substitute(t)), div_times , sep = "_"))

    # Set label col
    if (deparse(substitute(t)) == "devseq_tree") {
      cols <- c(rep("red", 1), rep("black", 6))
    } else {
      cols <- c("black", rep("red", 1), rep("black", 5))
    }


    # Convert ape phylogeny object to dendextend dendrogram object
    png(file = file.path(out_dir, "output", "plots", fname), 
    width = 2200, height = 3700, res = 800)
    par(mar = c(0.5, 5, 0, 0.25), bg = NA)

    p <- ggtree(t, size = 0.8) + theme_tree2(plot.margin = margin(40, 2, 5, -2.85), text = element_text(size = 17.7), 
        line = element_line(size = 0.8), axis.ticks.length = unit(0.2, "cm")) + geom_tiplab(color = cols, size = 5.35) + 
    ggplot2::scale_x_continuous(limits = c(0, 325), breaks = c(0,50,100,150)) + 
    ggplot2::scale_y_discrete(expand = c(0.1, 0)) + labs(caption = "Divergence time (Myr)        ") + 
    ggplot2::theme(
      axis.line = element_line(colour = "black", size = 0.8), 
      axis.ticks = element_line(colour = "black", size = 0.8), 
      axis.text.x = element_text(size = 14.25, margin = margin(3, 0, 2.5, 0), color = "black"))

    if (deparse(substitute(t)) == "devseq_tree") {
      p <- rotate(p, 13)
    }

    plot(p)

    dev.off()

  }

  suppressMessages(plotTree(devseq_tree))
  suppressMessages(plotTree(tetra_tree))


}


