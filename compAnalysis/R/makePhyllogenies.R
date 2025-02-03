# Construct expression phylogenies for orthologous protein coding genes
# Prepare Brawand and DevSeq comparative ortholog gene expression data
# Thresholds: DevSeq 0.05 ERCC; Brawand 0.5 TPM (no ERCC spike-ins available)
# Data input: Brawand and DevSeq TPM expression tables of all samples



#-------------------------------------- Read data tables ---------------------------------------


makePhyllogenies <- function(expr_estimation = c("TPM", "counts"), 
	coefficient = c("pearson", "spearman"), devseq_spec = c("Brassicaceae", "all")) {
	

   	# Show error message if expression estimation or unknown expression estimation is chosen
    if ((missing(expr_estimation)) || (!is.element(expr_estimation, c("TPM", "counts"))))
   
       stop(
       "Please choose one of the available expression estimations: 
	   'TPM', 'counts'",
	   call. = TRUE
       )

    # Show error message if no correlation or unknown correlation coefficient is chosen
    if ((missing(coefficient)) || (!is.element(coefficient, c("pearson", "spearman"))))
   
       stop(
       "Please choose one of the available correlation coefficients: 
	   'pearson', 'spearman'",
	   call. = TRUE
       )

    # Show error message if no devseq_spec or unknown devseq_spec is chosen
    if ((missing(devseq_spec)) || (!is.element(devseq_spec, c("Brassicaceae", "all"))))
   
       stop(
       "Please choose one of the available DevSeq species sets: 
	   'Brassicaceae', 'all'",
	   call. = TRUE
       )


    # Show startup message
    message("Reading data...")


	# Set expression input file
    if ((is.element("TPM", expr_estimation)) && (is.element("Brassicaceae", devseq_spec))) {
        genesExpr = file.path(in_dir, "Expression_data", "AT_brass_inter_tpm_mat_deseq_sample_names.csv")

    } else if ((is.element("TPM", expr_estimation)) && (is.element("all", devseq_spec))) {
        genesExpr = file.path(in_dir, "Expression_data", "AT_core_inter_tpm_mat_deseq_sample_names.csv")

    } else if ((is.element("counts", expr_estimation)) && (is.element("Brassicaceae", devseq_spec))) {
        genesExpr = file.path(in_dir, "Expression_data", "AT_brass_inter_count_mat_vsd_sample_names.csv")

    } else if ((is.element("counts", expr_estimation)) && (is.element("all", devseq_spec))) {
        genesExpr = file.path(in_dir, "Expression_data", "AT_core_inter_count_mat_vsd_sample_names.csv")
    }


	x <- read.table(genesExpr, sep=";", dec=".", header=TRUE, stringsAsFactors=FALSE)


    # Stop function here to allow specific analysis of a single data set
    # For DevSeq
    # return_list <- list("expr_estimation" = expr_estimation, "x" = x, "coefficient" = coefficient, "devseq_spec" = devseq_spec)
    # return(return_list)
    # }
    # return_objects <- makePhyllogenies(expr_estimation="counts", coefficient="pearson", devseq_spec="all") # read in DevSeq expression data
    # list2env(return_objects, envir = .GlobalEnv)

    # set column names
    colnames(x)[1] <- "gene_id"


    # Create "plots" folder in /out_dir/output/plots
    if (!dir.exists(file.path(out_dir, "output", "plots"))) 
        dir.create(file.path(out_dir, "output", "plots"), recursive = TRUE)

    # Show message
    message("Starting analysis and generate plots...")


    x[is.na(x)] <- 0 # replaces NAs by 0

    # Remove ERCC spike-ins from data
    x <- x[!grepl("ERCC", x$gene_id),]

    x_df <- x




#------------------------ Prepare DevSeq and Brawand data for phylotree ------------------------


    if(expr_estimation == "TPM") { x_df[,2:ncol(x_df)] <- log2(x_df[,2:ncol(x_df)] + 1) }


    # Replace pearson distance function from TreeExp2 package with metric pearson distance
    dist.pea = function (expMat = NULL) {

      object_n <- ncol(expMat)
      gene_n <- nrow(expMat)

      dis.mat <- matrix(0, nrow = object_n, ncol = object_n)


      for (i in 1:(object_n-1)) {

        for (j in (i+1):object_n) {

          dis.mat[j,i] <- sqrt(1/2*(1 - cor(expMat[,i],expMat[,j])))

        }

      }

      #browser()
      colnames(dis.mat) <- colnames(expMat)
      rownames(dis.mat) <- colnames(dis.mat)
      dis.mat  + t(dis.mat)

    }

    assignInNamespace("dist.pea", dist.pea, ns="TreeExp")



    # Construct distance matrix
    x_dist <- dist.pea(x_df[, 2:ncol(x_df)])


    # Perform analysis for all species and ATH as reference
    if (devseq_spec == "all") {

        # Construct organ gene expression phylogenies
        getPhyloTree <- function(org_dist) {

            spec_names <- rep(c("A.thaliana", "A.lyrata", "C.rubella", "E.salsugineum", "T.hassleriana",
                "M.truncatula", "B.distachyon"), each=3)
            repl_name <- rep(c("1", "2", "3"), 7)
            sample_name <- paste(spec_names, repl_name, sep=".")

            colnames(org_dist) <- sample_name
            rownames(org_dist) <- sample_name

            tr <- NJ(org_dist)

            tr <- root(tr, "B.distachyon.1", resolve.root = TRUE)

            return(tr)
        }

        root_tr <- getPhyloTree(org_dist = x_dist[1:21, 1:21])
        hypo_tr <- getPhyloTree(org_dist = x_dist[22:42, 22:42])
        leaf_tr <- getPhyloTree(org_dist = x_dist[43:63, 43:63])
        apex_v_tr <- getPhyloTree(org_dist = x_dist[64:84, 64:84])
        apex_i_tr <- getPhyloTree(org_dist = x_dist[85:105, 85:105])
        flower_tr <- getPhyloTree(org_dist = x_dist[106:126, 106:126])
        stamen_tr <- getPhyloTree(org_dist = x_dist[127:147, 127:147])
        carpel_tr <- getPhyloTree(org_dist = x_dist[148:168, 148:168])
        pollen_tr <- getPhyloTree(org_dist = x_dist[169:189, 169:189])


        # Perform boostrap analysis
        getPhyloBS <- function(tree, df) {

            spec_names <- rep(c("A.thaliana", "A.lyrata", "C.rubella", "E.salsugineum", "T.hassleriana",
                "M.truncatula", "B.distachyon"), each=3)
            repl_name <- rep(c("1", "2", "3"), 7)
            sample_name <- paste(spec_names, repl_name, sep=".")

            colnames(df) <- sample_name

            bs <- boot.exphy(phy = tree, x = df, method = 'pea', outgroup = "B.distachyon.1", 
                B = 1000, trees = TRUE)

            return(bs)
        }

        root_bs <- getPhyloBS(tree = root_tr, df = x_df[,2:22])
        hypo_bs <- getPhyloBS(tree = hypo_tr, df = x_df[,23:43])
        leaf_bs <- getPhyloBS(tree = leaf_tr, df = x_df[,44:64])
        apex_v_bs <- getPhyloBS(tree = apex_v_tr, df = x_df[,65:85])
        apex_i_bs <- getPhyloBS(tree = apex_i_tr, df = x_df[,86:106])
        flower_bs <- getPhyloBS(tree = flower_tr, df = x_df[,107:127])
        stamen_bs <- getPhyloBS(tree = stamen_tr, df = x_df[,128:148])
        carpel_bs <- getPhyloBS(tree = carpel_tr, df = x_df[,149:169])
        pollen_bs <- getPhyloBS(tree = pollen_tr, df = x_df[,170:190])


        # Get total tree length from bootstraps
        getBSValues <- function(bs) {

            tree_length <- sapply(bs, function(x){sum(x[]$edge.length)})

            return(tree_length)
        }

        root_bsl <- getBSValues(root_bs$trees)
        hypo_bsl <- getBSValues(hypo_bs$trees)
        leaf_bsl <- getBSValues(leaf_bs$trees)
        apex_v_bsl <- getBSValues(apex_v_bs$trees)
        apex_i_bsl <- getBSValues(apex_i_bs$trees)
        flower_bsl <- getBSValues(flower_bs$trees)
        stamen_bsl <- getBSValues(stamen_bs$trees)
        carpel_bsl <- getBSValues(carpel_bs$trees)
        pollen_bsl <- getBSValues(pollen_bs$trees)


        # Check bootstraps e.g.
        # root_tr$node.label = root_bs$BP 
        # plot(root_tr, show.node.label = TRUE)
        # Generate enhanced phylo trees using ggtree
        # Set ladderize = TRUE to display ladderized tree
        # ladderize = FALSE is the default in plot.phylo()
        # bootstrap colors: < 0.95 = yellow; < 0.9 = #ffb100 (orange); 
        # < 0.7 = #eb0000 (red)

        # Make phyloplots
        # Root
        png(height = 1750, width = 1450, pointsize = 100, res = 325, file = file.path(out_dir, "output", "plots", "root_tr.png"))
        p <- ggtree(root_tr, ladderize = FALSE, size=0.55) + 
        geom_tiplab() + 
        geom_nodepoint(shape=21, size=2.7, fill='white', alpha=1, stroke=0.55) + 
        geom_point2(aes(subset=(node==37)), shape=21, size=2.75, fill='#ffb100', alpha=1, stroke=0.1) +
        geom_point2(aes(subset=(node==24)), shape=21, size=2.75, fill='#ffb100', alpha=1, stroke=0.1) + 
        geom_point2(aes(subset=(node==25)), shape=21, size=2.75, fill='#ffb100', alpha=1, stroke=0.1) + 
        geom_point2(aes(subset=(node==40)), shape=21, size=2.75, fill='#eb0000', alpha=1, stroke=0.1) + 
        geom_treescale(color="white", y = 0.75) + 
        # geom_text(aes(label=node), hjust=-.3) + 
        xlim(0, 0.671) + 
        theme_tree2(plot.margin=margin(5, 5, 35, 5), text = element_text(size = 12.5), 
            line = element_line(size = 0.4), axis.ticks.length=unit(.125, "cm"), legend.position = "none")
        # Rotate nodes
        p2 <- rotate(p, 23) %>% rotate(24) %>% rotate(28) %>% rotate(34)
        plot(p2)
        dev.off()

        # Hypocotyl
        hypo_tr2 <- groupClade(hypo_tr, c(30, 34))
        png(height = 1750, width = 1450, pointsize = 100, res = 325, file = file.path(out_dir, "output", "plots", "hypo_tr.png"))
        p <- ggtree(hypo_tr2, ladderize = FALSE, size=0.55, aes(linetype=group)) + 
        geom_tiplab() + 
        geom_nodepoint(shape=21, size=2.7, fill='white', alpha=1, stroke=0.55) + 
        geom_point2(aes(subset=(node==37)), shape=21, size=2.75, fill='#ffb100', alpha=1, stroke=0.1) + 
        geom_point2(aes(subset=(node==35)), shape=21, size=2.75, fill='#ffb100', alpha=1, stroke=0.1) +
        geom_point2(aes(subset=(node==31)), shape=21, size=2.75, fill='#ffb100', alpha=1, stroke=0.1) + 
        geom_point2(aes(subset=(node==29)), shape=21, size=2.75, fill='#eb0000', alpha=1, stroke=0.1) + 
        geom_treescale(color="white", y = 0.75) + 
        # geom_text(aes(label=node), hjust=-.3) + 
        xlim(0, 0.671) + 
        theme_tree2(plot.margin=margin(5, 5, 35, 5), text = element_text(size = 12.5), 
            line = element_line(size = 0.4), axis.ticks.length=unit(.125, "cm"), legend.position = "none")
        # Rotate nodes
        p2 <- flip(p, 25, 30)
        plot(p2)
        dev.off()

        # Leaf
        png(height = 1750, width = 1450, pointsize = 100, res = 325, file = file.path(out_dir, "output", "plots", "leaf_tr.png"))
        p <- ggtree(leaf_tr, ladderize = FALSE, size=0.55) + 
        geom_tiplab() + 
        geom_nodepoint(shape=21, size=2.7, fill='white', alpha=1, stroke=0.55) + 
        geom_point2(aes(subset=(node==39)), shape=21, size=2.75, fill='#ffb100', alpha=1, stroke=0.1) + 
        geom_point2(aes(subset=(node==29)), shape=21, size=2.75, fill='#ffb100', alpha=1, stroke=0.1) +
        geom_treescale(color="white", y = 0.75) + 
        # geom_text(aes(label=node), hjust=-.3) + 
        xlim(0, 0.671) + 
        theme_tree2(plot.margin=margin(5, 5, 35, 5), text = element_text(size = 12.5), 
            line = element_line(size = 0.4), axis.ticks.length=unit(.125, "cm"), legend.position = "none")
        # Rotate nodes
        p2 <- flip(p, 25, 30)
        plot(p2)
        dev.off()

        # Apex veg
        png(height = 1750, width = 1450, pointsize = 100, res = 325, file = file.path(out_dir, "output", "plots", "apex_v_tr.png"))
        p <- ggtree(apex_v_tr, ladderize = FALSE, size=0.55) + 
        geom_tiplab() + 
        geom_nodepoint(shape=21, size=2.7, fill='white', alpha=1, stroke=0.55) + 
        geom_point2(aes(subset=(node==37)), shape=21, size=2.75, fill='#ffb100', alpha=1, stroke=0.1) + 
        geom_point2(aes(subset=(node==35)), shape=21, size=2.75, fill='yellow', alpha=1, stroke=0.1) + 
        geom_point2(aes(subset=(node==33)), shape=21, size=2.75, fill='#eb0000', alpha=1, stroke=0.1) + 
        geom_treescale(color="white", y = 0.75) + 
        # geom_text(aes(label=node), hjust=-.3) + 
        xlim(0, 0.671) + 
        theme_tree2(plot.margin=margin(5, 5, 35, 5), text = element_text(size = 12.5), 
            line = element_line(size = 0.4), axis.ticks.length=unit(.125, "cm"), legend.position = "none")
        # Rotate nodes
        p2 <- flip(p, 28,30)
        plot(p2)
        dev.off()

        # Apex inf
        png(height = 1750, width = 1450, pointsize = 100, res = 325, file = file.path(out_dir, "output", "plots", "apex_i_tr.png"))
        p <- ggtree(apex_i_tr, ladderize = FALSE, size=0.55) + 
        geom_tiplab() + 
        geom_nodepoint(shape=21, size=2.7, fill='white', alpha=1, stroke=0.55) + 
        geom_point2(aes(subset=(node==39)), shape=21, size=2.75, fill='#eb0000', alpha=1, stroke=0.1) + 
        geom_point2(aes(subset=(node==33)), shape=21, size=2.75, fill='yellow', alpha=1, stroke=0.1) + 
        geom_point2(aes(subset=(node==31)), shape=21, size=2.75, fill='#eb0000', alpha=1, stroke=0.1) + 
        geom_treescale(color="white", y = 0.75) + 
        # geom_text(aes(label=node), hjust=-.3) + 
        xlim(0, 0.671) + 
        theme_tree2(plot.margin=margin(5, 5, 35, 5), text = element_text(size = 12.5), 
            line = element_line(size = 0.4), axis.ticks.length=unit(.125, "cm"), legend.position = "none")
        # Rotate nodes
        p2 <- rotate(p, 23) %>% rotate(24) %>% rotate(34)
        plot(p2)
        dev.off()

        # Flower
        flower_tr2 <- groupClade(flower_tr, c(32, 34))
        png(height = 1750, width = 1450, pointsize = 100, res = 325, file = file.path(out_dir, "output", "plots", "flower_tr.png"))
        p <- ggtree(flower_tr2, ladderize = FALSE, size=0.55, aes(linetype=group)) + 
        geom_tiplab() + 
        geom_nodepoint(shape=21, size=2.7, fill='white', alpha=1, stroke=0.55) + 
        geom_point2(aes(subset=(node==39)), shape=21, size=2.75, fill='#eb0000', alpha=1, stroke=0.1) + 
        geom_point2(aes(subset=(node==33)), shape=21, size=2.75, fill='#eb0000', alpha=1, stroke=0.1) + 
        geom_point2(aes(subset=(node==23)), shape=21, size=2.75, fill='#ffb100', alpha=1, stroke=0.1) + 
        geom_treescale(color="white", y = 0.75) + 
        # geom_text(aes(label=node), hjust=-.3) + 
        xlim(0, 0.671) + 
        theme_tree2(plot.margin=margin(5, 5, 35, 5), text = element_text(size = 12.5), 
            line = element_line(size = 0.4), axis.ticks.length=unit(.125, "cm"), legend.position = "none")
        # Rotate nodes
        p2 <- rotate(p, 24) %>% rotate(25) %>% rotate(36)
        plot(p2)
        dev.off()

        # Stamen
        stamen_tr2 <- groupClade(stamen_tr, c(36, 38))
        png(height = 1750, width = 1450, pointsize = 100, res = 325, file = file.path(out_dir, "output", "plots", "stamen_tr.png"))
        p <- ggtree(stamen_tr2, ladderize = FALSE, size=0.55, aes(linetype=group)) + 
        geom_tiplab() + 
        scale_linetype_manual(values=c("solid", 42, 22)) + 
        geom_nodepoint(shape=21, size=2.7, fill='white', alpha=1, stroke=0.55) + 
        geom_point2(aes(subset=(node==31)), shape=21, size=2.75, fill='#eb0000', alpha=1, stroke=0.1) + 
        geom_point2(aes(subset=(node==24)), shape=21, size=2.75, fill='#ffb100', alpha=1, stroke=0.1) + 
        geom_point2(aes(subset=(node==23)), shape=21, size=2.75, fill='#ffb100', alpha=1, stroke=0.1) + 
        geom_treescale(color="white", y = 0.75) + 
        # geom_text(aes(label=node), hjust=-.3) + 
        xlim(0, 0.671) + 
        theme_tree2(plot.margin=margin(5, 5, 35, 5), text = element_text(size = 12.5), 
            line = element_line(size = 0.4), axis.ticks.length=unit(.125, "cm"), legend.position = "none")
        # Rotate nodes
        p2 <- rotate(p, 23) %>% rotate(27) %>% rotate(24)
        plot(p2)
        dev.off()

        # Carpel
        png(height = 1750, width = 1450, pointsize = 100, res = 325, file = file.path(out_dir, "output", "plots", "carpel_tr.png"))
        p <- ggtree(carpel_tr, ladderize = FALSE, size=0.55) + 
        geom_tiplab() + 
        geom_nodepoint(shape=21, size=2.7, fill='white', alpha=1, stroke=0.55) + 
        geom_point2(aes(subset=(node==39)), shape=21, size=2.75, fill='yellow', alpha=1, stroke=0.1) + 
        geom_point2(aes(subset=(node==37)), shape=21, size=2.75, fill='#eb0000', alpha=1, stroke=0.1) + 
        geom_point2(aes(subset=(node==33)), shape=21, size=2.75, fill='#eb0000', alpha=1, stroke=0.1) + 
        geom_point2(aes(subset=(node==31)), shape=21, size=2.75, fill='#ffb100', alpha=1, stroke=0.1) + 
        geom_treescale(color="white", y = 0.75) + 
        # geom_text(aes(label=node), hjust=-.3) + 
        xlim(0, 0.671) + 
        theme_tree2(plot.margin=margin(5, 5, 35, 5), text = element_text(size = 12.5), 
            line = element_line(size = 0.4), axis.ticks.length=unit(.125, "cm"), legend.position = "none")
        # Rotate nodes
        p2 <- rotate(p, 23) %>% rotate(25)
        plot(p2)
        dev.off()

        # Pollen
        png(height = 1750, width = 1450, pointsize = 100, res = 325, file = file.path(out_dir, "output", "plots", "pollen_tr.png"))
        p <- ggtree(pollen_tr, ladderize = FALSE, size=0.55) + 
        geom_tiplab() + 
        geom_nodepoint(shape=21, size=2.7, fill='white', alpha=1, stroke=0.55) + 
        geom_point2(aes(subset=(node==39)), shape=21, size=2.75, fill='#ffb100', alpha=1, stroke=0.1) + 
        geom_point2(aes(subset=(node==35)), shape=21, size=2.75, fill='#ffb100', alpha=1, stroke=0.1) + 
        geom_point2(aes(subset=(node==33)), shape=21, size=2.75, fill='#eb0000', alpha=1, stroke=0.1) + 
        geom_point2(aes(subset=(node==31)), shape=21, size=2.75, fill='#ffb100', alpha=1, stroke=0.1) + 
        geom_point2(aes(subset=(node==29)), shape=21, size=2.75, fill='#ffb100', alpha=1, stroke=0.1) + 
        geom_treescale(color="white", y = 0.75) + 
        # geom_text(aes(label=node), hjust=-.3) + 
        xlim(0, 0.671) + 
        theme_tree2(plot.margin=margin(5, 5, 35, 5), text = element_text(size = 12.5), 
            line = element_line(size = 0.4), axis.ticks.length=unit(.125, "cm"), legend.position = "none")
        # Rotate nodes
        p2 <- flip(p, 25, 32)
        plot(p2)
        dev.off()



        # Generate df with combined boostrap lengths and reshape data for ggplot2
        ttree_length <- data.frame(tree_length = c(root_bsl, hypo_bsl, leaf_bsl, apex_v_bsl, apex_i_bsl, 
            flower_bsl, carpel_bsl, stamen_bsl, pollen_bsl))
        organs_names <- data.frame(organ = rep(c("Root", "Hypocotyl", "Leaf", "veg_apex", "inf_apex", 
            "Flower", "Carpel", "Stamen", "Pollen"), each = length(root_bsl)))

        original_tree <- c(sum(root_tr$edge.length), sum(hypo_tr$edge.length), sum(leaf_tr$edge.length), 
            sum(apex_v_tr$edge.length), sum(apex_i_tr$edge.length), sum(flower_tr$edge.length), 
            sum(carpel_tr$edge.length), sum(stamen_tr$edge.length), sum(pollen_tr$edge.length))
        original_tree_df <- data.frame(orig_tree_length = rep(original_tree, each = length(root_bsl)))

        phyloBSTreeL <- cbind(organs_names, ttree_length, original_tree_df)
        phyloBSTreeL$tree_length <- as.numeric(phyloBSTreeL$tree_length)
        phyloBSTreeL$orig_tree_length <- as.numeric(phyloBSTreeL$orig_tree_length)
        phyloBSTreeL$organ <- factor(phyloBSTreeL$organ)

        write.table(phyloBSTreeL, file=file.path(out_dir, "output", "data", "phyloBSTreeL_all_ATH.txt"), 
            sep="\t", col.names=TRUE, row.names=FALSE, dec=".", quote = FALSE)


    } else if (devseq_spec == "Brassicaceae") {

        # Construct organ gene expression phylogenies
        getPhyloTree <- function(org_dist) {

            spec_names <- rep(c("A.thaliana", "A.lyrata", "C.rubella", "E.salsugineum"), each=3)
            repl_name <- rep(c("1", "2", "3"), 4)
            sample_name <- paste(spec_names, repl_name, sep=".")

            colnames(org_dist) <- sample_name
            rownames(org_dist) <- sample_name

            tr <- NJ(org_dist)

            tr <- root(tr, "E.salsugineum.1", resolve.root = TRUE)

            return(tr)
        }

        root_tr <- getPhyloTree(org_dist = x_dist[1:12, 1:12])
        hypo_tr <- getPhyloTree(org_dist = x_dist[13:24, 13:24])
        leaf_tr <- getPhyloTree(org_dist = x_dist[25:36, 25:36])
        apex_v_tr <- getPhyloTree(org_dist = x_dist[37:48, 37:48])
        apex_i_tr <- getPhyloTree(org_dist = x_dist[49:60, 49:60])
        flower_tr <- getPhyloTree(org_dist = x_dist[61:72, 61:72])
        stamen_tr <- getPhyloTree(org_dist = x_dist[73:84, 73:84])
        carpel_tr <- getPhyloTree(org_dist = x_dist[85:96, 85:96])
        pollen_tr <- getPhyloTree(org_dist = x_dist[97:108, 97:108])


        # Perform boostrap analysis
        getPhyloBS <- function(tree, df) {

            spec_names <- rep(c("A.thaliana", "A.lyrata", "C.rubella", "E.salsugineum"), each=3)
            repl_name <- rep(c("1", "2", "3"), 4)
            sample_name <- paste(spec_names, repl_name, sep=".")

            colnames(df) <- sample_name

            bs <- boot.exphy(phy = tree, x = df, method = 'pea', outgroup = "E.salsugineum.1", 
                B = 1000, trees = TRUE)

            return(bs)
        }

        root_bs <- getPhyloBS(tree = root_tr, df = x_df[,2:13])
        hypo_bs <- getPhyloBS(tree = hypo_tr, df = x_df[,14:25])
        leaf_bs <- getPhyloBS(tree = leaf_tr, df = x_df[,26:37])
        apex_v_bs <- getPhyloBS(tree = apex_v_tr, df = x_df[,38:49])
        apex_i_bs <- getPhyloBS(tree = apex_i_tr, df = x_df[,50:61])
        flower_bs <- getPhyloBS(tree = flower_tr, df = x_df[,62:73])
        stamen_bs <- getPhyloBS(tree = stamen_tr, df = x_df[,74:85])
        carpel_bs <- getPhyloBS(tree = carpel_tr, df = x_df[,86:97])
        pollen_bs <- getPhyloBS(tree = pollen_tr, df = x_df[,98:109])


        # Get total tree length from bootstraps
        getBSValues <- function(bs) {

            tree_length <- sapply(bs, function(x){sum(x[]$edge.length)})

            return(tree_length)
        }

        root_bsl <- getBSValues(root_bs$trees)
        hypo_bsl <- getBSValues(hypo_bs$trees)
        leaf_bsl <- getBSValues(leaf_bs$trees)
        apex_v_bsl <- getBSValues(apex_v_bs$trees)
        apex_i_bsl <- getBSValues(apex_i_bs$trees)
        flower_bsl <- getBSValues(flower_bs$trees)
        stamen_bsl <- getBSValues(stamen_bs$trees)
        carpel_bsl <- getBSValues(carpel_bs$trees)
        pollen_bsl <- getBSValues(pollen_bs$trees)



        # Generate df with combined boostrap lengths and reshape data for ggplot2
        ttree_length <- data.frame(tree_length = c(root_bsl, hypo_bsl, leaf_bsl, apex_v_bsl, apex_i_bsl, 
            flower_bsl, carpel_bsl, stamen_bsl, pollen_bsl))
        organs_names <- data.frame(organ = rep(c("Root", "Hypocotyl", "Leaf", "veg_apex", "inf_apex", 
            "Flower", "Carpel", "Stamen", "Pollen"), each = length(root_bsl)))

        original_tree <- c(sum(root_tr$edge.length), sum(hypo_tr$edge.length), sum(leaf_tr$edge.length), 
            sum(apex_v_tr$edge.length), sum(apex_i_tr$edge.length), sum(flower_tr$edge.length), 
            sum(carpel_tr$edge.length), sum(stamen_tr$edge.length), sum(pollen_tr$edge.length))
        original_tree_df <- data.frame(orig_tree_length = rep(original_tree, each = length(root_bsl)))

        phyloBSTreeL <- cbind(organs_names, ttree_length, original_tree_df)
        phyloBSTreeL$tree_length <- as.numeric(phyloBSTreeL$tree_length)
        phyloBSTreeL$orig_tree_length <- as.numeric(phyloBSTreeL$orig_tree_length)
        phyloBSTreeL$organ <- factor(phyloBSTreeL$organ)

        write.table(phyloBSTreeL, file=file.path(out_dir, "output", "data", "phyloBSTreeL_Brassicaceae_ATH.txt"), 
            sep="\t", col.names=TRUE, row.names=FALSE, dec=".", quote = FALSE)

    }


}








