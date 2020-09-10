# Prepare Brawand and DevSeq phylogenetic tree and ortholog plot
# Data input: Pairwise divergence times and orthologr ortholog gene output
# Pairwise TimeTree divergence were obtained from TimeTree (http://www.timetree.org/)


#-------------------------------------- Read data tables ---------------------------------------


plotPhyloCore <- function(div_times = c("Median", "Estimated")) {


	if(is.element("Median", div_times)) {

		devseq_phylo <- "((((((Arabidopsis_thaliana:5.8, Arabidopsis_lyrata:5.8):2.7,Capsella_rubella:8.5):18.9,Eutrema_salsugineum:27.4):23.6,Tarenaya_hassleriana:51):57,Medicago_truncatula:108):44,Brachypodium_distachyon:152);"


	 } else if(is.element("Estimated", div_times)) {

	 	devseq_phylo <- "((((((Arabidopsis_thaliana:7.1, Arabidopsis_lyrata:7.1):2.3,Capsella_rubella:9.4):16.2,Eutrema_salsugineum:25.6):20.4,Tarenaya_hassleriana:46):60,Medicago_truncatula:106):54,Brachypodium_distachyon:160);"
	 }


	 vert_tree <- read.tree(text = devseq_phylo)



	 # Create "plots" folder in /out_dir/output/plots
     if (!dir.exists(file.path(out_dir, "output", "plots"))) 
	 dir.create(file.path(out_dir, "output", "plots"), recursive = TRUE)


     if(is.element("Median", div_times)) {

     	treename <- sprintf('%s.png', paste("Marker_species_phylogeny", "median_times" , sep="_"))

      } else if(is.element("Estimated", div_times)) {

      	treename <- sprintf('%s.png', paste("Marker_species_phylogeny", "estimated_times" , sep="_"))
      }

     
     # Make marker species phylogeny plot
	 png(file = file.path(out_dir, "output", "plots", treename), 
		width = 3250, height = 3700, res = 800)
	 par(mar = c(0.25, 0.25, 0, 0.25), bg=NA)

	 plot(vert_tree, type = "phylogram", use.edge.length = TRUE, show.node.label = FALSE, 
	 	edge.width = 1.5, edge.lty = 1, font = 3, root.edge = FALSE, label.offset = 2, 
	 	direction = "upwards")

	 add.scale.bar(x = 1, y = 1, cex = 0.8, font = 5, col = "red")

	 dev.off()

}

plotPhyloCore("Estimated")
plotPhyloCore("Median")

