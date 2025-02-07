# Construct expression phylogenies for orthologous lncRNAs
# Thresholds: DevSeq 0.05 ERCC
# Data input: DevSeq VST count and TPM expression tables of Brassicaceae samples



#-------------------------------------- Read data tables ---------------------------------------


getNCPhyllogenies <- function(expr_estimation = c("TPM", "counts"), 
    coefficient = c("pearson", "spearman"), transcripttype = c("coding", "non-coding")) {
	

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

    # Show error message if no transcript type or unknown type is chosen
    if ((missing(transcripttype)) || (!is.element(transcripttype, c("coding", "non-coding"))))
   
       stop(
       "Please choose one of the available transcript types: 
       'coding', 'non-coding'",
       call. = TRUE
       )


    # Show startup message
    message("Reading data...")


    # Set expression input file
    if ((is.element("TPM", expr_estimation)) && (is.element("non-coding", transcripttype))) {
        genesExpr = file.path(in_dir, "Expression_data", "lnc_AT_brass_inter_combined_tpm_mat_deseq_sample_names.csv")

    } else if ((is.element("counts", expr_estimation)) && (is.element("non-coding", transcripttype))) {
        genesExpr = file.path(in_dir, "Expression_data", "lnc_AT_brass_inter_combined_count_mat_vsd_sample_names.csv")
    
    } else if ((is.element("TPM", expr_estimation)) && (is.element("coding", transcripttype))) {
        genesExpr = file.path(in_dir, "Expression_data", "AT_brass_inter_tpm_mat_deseq_sample_names.csv")

    } else if ((is.element("counts", expr_estimation)) && (is.element("coding", transcripttype))) {
        genesExpr = file.path(in_dir, "Expression_data", "AT_brass_inter_count_mat_vsd_sample_names.csv")
    }


    x <- read.table(genesExpr, sep=";", dec=".", header=TRUE, stringsAsFactors=FALSE)


    # Stop function here to allow specific analysis of a single data set
    # return_list <- list("expr_estimation" = expr_estimation, "x" = x, "coefficient" = coefficient, "transcripttype" = transcripttype)
    # return(return_list)
    # }
    # return_objects <- getNCPhyllogenies(expr_estimation="counts", coefficient="spearman", transcripttype="non-coding")
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

    dist.spe = function (expMat = NULL) {

      object_n <- ncol(expMat)
      gene_n <- nrow(expMat)

      dis.mat <- matrix(0, nrow = object_n, ncol = object_n)


      for (i in 1:(object_n-1)) {

        for (j in (i+1):object_n) {

          dis.mat[j,i] <- sqrt(1/2*(1 - cor(expMat[,i],expMat[,j],
                                  method = "spearman")))

        }

      }

      #browser()
      colnames(dis.mat) <- colnames(expMat)
      rownames(dis.mat) <- colnames(dis.mat)
      dis.mat  + t(dis.mat)

    }

    assignInNamespace("dist.spe", dist.spe, ns="TreeExp")


    # Construct distance matrix
    if (coefficient == "pearson") {

        x_dist <- dist.pea(x_df[, 2:ncol(x_df)])

    } else if (coefficient == "spearman") {

        x_dist <- dist.spe(x_df[, 2:ncol(x_df)])
    }


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

        if (coefficient == "pearson") {

            dist_method <- 'pea'

        } else if (coefficient == "spearman") {

            dist_method <- 'spe'
        }

        spec_names <- rep(c("A.thaliana", "A.lyrata", "C.rubella", "E.salsugineum"), each=3)
        repl_name <- rep(c("1", "2", "3"), 4)
        sample_name <- paste(spec_names, repl_name, sep=".")

        colnames(df) <- sample_name

        bs <- boot.exphy(phy = tree, x = df, method = dist_method, outgroup = "E.salsugineum.1", 
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


    # Check bootstraps e.g.
    # root_tr$node.label = root_bs$BP 
    # plot(root_tr, show.node.label = TRUE)
    # nodelabels(cex=0.5, width=0.1)
    # Generate enhanced phylo trees using ggtree
    # Set ladderize = TRUE to display ladderized tree
    # ladderize = FALSE is the default in plot.phylo()
    # bootstrap colors: < 0.95 = yellow; < 0.9 = #ffbc00 (orange); 
    # < 0.7 = red; < 0.5 = #7e0000

    # Make phyloplots
    if ((coefficient == "spearman") && (expr_estimation ==  "counts") && (transcripttype == "non-coding")) {

    # Root
    png(height = 1150, width = 1100, pointsize = 100, res = 325, file = file.path(out_dir, "output", "plots", "lnc_root_tr.png"))
    p <- ggtree(root_tr, ladderize = FALSE, size=0.55) + 
    geom_tiplab() + 
    geom_nodepoint(shape=21, size=2.8, fill='white', alpha=1, stroke=0.55) + 
    geom_point2(aes(subset=(node==14)), shape=21, size=2.85, fill='#ffbc00', alpha=1, stroke=0.1) +
    geom_point2(aes(subset=(node==17)), shape=21, size=2.85, fill='red', alpha=1, stroke=0.1) + 
    geom_point2(aes(subset=(node==18)), shape=21, size=2.85, fill='red', alpha=1, stroke=0.1) + 
    geom_point2(aes(subset=(node==20)), shape=21, size=2.85, fill='#7e0000', alpha=1, stroke=0.1) + 
    geom_point2(aes(subset=(node==22)), shape=21, size=2.85, fill='red', alpha=1, stroke=0.1) + 
    geom_treescale(color="white", y = 0.75) + 
    # geom_text(aes(label=node), hjust=-.3) + 
    xlim(-0.005, 0.565) + 
    theme_tree2(plot.margin=margin(5, 5, 20, 5), text = element_text(size = 11.75), 
        line = element_line(size = 0.4), axis.ticks.length=unit(.125, "cm"))
    p2 <- flip(p, 19, 21)
    plot(p2)
    dev.off()

    # Hypocotyl
    png(height = 1150, width = 1100, pointsize = 100, res = 325, file = file.path(out_dir, "output", "plots", "lnc_hypo_tr.png"))
    p <- ggtree(hypo_tr, ladderize = FALSE, size=0.55) + 
    geom_tiplab() + 
    geom_nodepoint(shape=21, size=2.8, fill='white', alpha=1, stroke=0.55) + 
    geom_point2(aes(subset=(node==14)), shape=21, size=2.85, fill='#ffbc00', alpha=1, stroke=0.1) + 
    geom_point2(aes(subset=(node==19)), shape=21, size=2.85, fill='#ffbc00', alpha=1, stroke=0.1) +
    geom_point2(aes(subset=(node==16)), shape=21, size=2.85, fill='red', alpha=1, stroke=0.1) + 
    geom_point2(aes(subset=(node==20)), shape=21, size=2.85, fill='#7e0000', alpha=1, stroke=0.1) + 
    geom_point2(aes(subset=(node==22)), shape=21, size=2.85, fill='#7e0000', alpha=1, stroke=0.1) + 
    geom_treescale(color="white", y = 0.75) + 
    # geom_text(aes(label=node), hjust=-.3) + 
    xlim(-0.005, 0.565) + 
    theme_tree2(plot.margin=margin(5, 5, 20, 5), text = element_text(size = 11.75), 
        line = element_line(size = 0.4), axis.ticks.length=unit(.125, "cm"))
    # Rotate nodes
    plot(p)
    dev.off()

    # Leaf
    png(height = 1150, width = 1100, pointsize = 100, res = 325, file = file.path(out_dir, "output", "plots", "lnc_leaf_tr.png"))
    p <- ggtree(leaf_tr, ladderize = FALSE, size=0.55) + 
    geom_tiplab() + 
    geom_nodepoint(shape=21, size=2.8, fill='white', alpha=1, stroke=0.55) + 
    geom_point2(aes(subset=(node==18)), shape=21, size=2.85, fill='yellow', alpha=1, stroke=0.1) + 
    geom_point2(aes(subset=(node==22)), shape=21, size=2.85, fill='yellow', alpha=1, stroke=0.1) + 
    geom_point2(aes(subset=(node==14)), shape=21, size=2.85, fill='#7e0000', alpha=1, stroke=0.1) + 
    geom_point2(aes(subset=(node==17)), shape=21, size=2.85, fill='red', alpha=1, stroke=0.1) + 
    geom_point2(aes(subset=(node==20)), shape=21, size=2.85, fill='red', alpha=1, stroke=0.1) + 
    geom_treescale(color="white", y = 0.75) + 
    # geom_text(aes(label=node), hjust=-.3) + 
    xlim(-0.005, 0.565) + 
    theme_tree2(plot.margin=margin(5, 5, 20, 5), text = element_text(size = 11.75), 
        line = element_line(size = 0.4), axis.ticks.length=unit(.125, "cm"))
    # Rotate nodes
    p2 <- rotate(p, 15) %>% rotate(14)
    p3 <- flip(p2, 19, 21)
    plot(p3)
    dev.off()

    # Apex veg
    png(height = 1150, width = 1100, pointsize = 100, res = 325, file = file.path(out_dir, "output", "plots", "lnc_apex_v_tr.png"))
    p <- ggtree(apex_v_tr, ladderize = FALSE, size=0.55) + 
    geom_tiplab() + 
    geom_nodepoint(shape=21, size=2.8, fill='white', alpha=1, stroke=0.55) + 
    geom_point2(aes(subset=(node==19)), shape=21, size=2.85, fill='#ffbc00', alpha=1, stroke=0.1) + 
    geom_point2(aes(subset=(node==14)), shape=21, size=2.85, fill='red', alpha=1, stroke=0.1) + 
    geom_point2(aes(subset=(node==16)), shape=21, size=2.85, fill='red', alpha=1, stroke=0.1) + 
    geom_point2(aes(subset=(node==20)), shape=21, size=2.85, fill='red', alpha=1, stroke=0.1) + 
    geom_point2(aes(subset=(node==22)), shape=21, size=2.85, fill='red', alpha=1, stroke=0.1) + 
    geom_treescale(color="white", y = 0.75) + 
    # geom_text(aes(label=node), hjust=-.3) + 
    xlim(-0.005, 0.565) + 
    theme_tree2(plot.margin=margin(5, 5, 20, 5), text = element_text(size = 11.75), 
        line = element_line(size = 0.4), axis.ticks.length=unit(.125, "cm"))
    # Rotate nodes
    plot(p)
    dev.off()

    # Apex inf
    png(height = 1150, width = 1100, pointsize = 100, res = 325, file = file.path(out_dir, "output", "plots", "lnc_apex_i_tr.png"))
    p <- ggtree(apex_i_tr, ladderize = FALSE, size=0.55) + 
    geom_tiplab() + 
    geom_nodepoint(shape=21, size=2.8, fill='white', alpha=1, stroke=0.55) + 
    geom_point2(aes(subset=(node==14)), shape=21, size=2.85, fill='#ffbc00', alpha=1, stroke=0.1) + 
    geom_point2(aes(subset=(node==16)), shape=21, size=2.85, fill='#ffbc00', alpha=1, stroke=0.1) + 
    geom_point2(aes(subset=(node==18)), shape=21, size=2.85, fill='red', alpha=1, stroke=0.1) + 
    geom_point2(aes(subset=(node==20)), shape=21, size=2.85, fill='red', alpha=1, stroke=0.1) + 
    geom_point2(aes(subset=(node==22)), shape=21, size=2.85, fill='#7e0000', alpha=1, stroke=0.1) + 
    geom_treescale(color="white", y = 0.75) + 
    # geom_text(aes(label=node), hjust=-.3) + 
    xlim(-0.005, 0.565) + 
    theme_tree2(plot.margin=margin(5, 5, 20, 5), text = element_text(size = 11.75), 
        line = element_line(size = 0.4), axis.ticks.length=unit(.125, "cm"))
    plot(p)
    dev.off()

    # Flower
    png(height = 1150, width = 1100, pointsize = 100, res = 325, file = file.path(out_dir, "output", "plots", "lnc_flower_tr.png"))
    p <- ggtree(flower_tr, ladderize = FALSE, size=0.55) + 
    geom_tiplab() + 
    geom_nodepoint(shape=21, size=2.8, fill='white', alpha=1, stroke=0.55) + 
    geom_point2(aes(subset=(node==14)), shape=21, size=2.85, fill='red', alpha=1, stroke=0.1) + 
    geom_point2(aes(subset=(node==16)), shape=21, size=2.85, fill='red', alpha=1, stroke=0.1) + 
    geom_point2(aes(subset=(node==18)), shape=21, size=2.85, fill='red', alpha=1, stroke=0.1) + 
    geom_point2(aes(subset=(node==20)), shape=21, size=2.85, fill='red', alpha=1, stroke=0.1) + 
    geom_point2(aes(subset=(node==22)), shape=21, size=2.85, fill='red', alpha=1, stroke=0.1) + 
    geom_treescale(color="white", y = 0.75) + 
    # geom_text(aes(label=node), hjust=-.3) + 
    xlim(-0.005, 0.565) + 
    theme_tree2(plot.margin=margin(5, 5, 20, 5), text = element_text(size = 11.75), 
        line = element_line(size = 0.4), axis.ticks.length=unit(.125, "cm"))
    plot(p)
    dev.off()

    # Stamen
    png(height = 1150, width = 1100, pointsize = 100, res = 325, file = file.path(out_dir, "output", "plots", "lnc_stamen_tr.png"))
    p <- ggtree(stamen_tr, ladderize = FALSE, size=0.55) + 
    geom_tiplab() + 
    geom_nodepoint(shape=21, size=2.8, fill='white', alpha=1, stroke=0.55) + 
    geom_point2(aes(subset=(node==14)), shape=21, size=2.85, fill='#ffbc00', alpha=1, stroke=0.1) + 
    geom_point2(aes(subset=(node==19)), shape=21, size=2.85, fill='#ffbc00', alpha=1, stroke=0.1) + 
    geom_point2(aes(subset=(node==16)), shape=21, size=2.85, fill='red', alpha=1, stroke=0.1) + 
    geom_point2(aes(subset=(node==21)), shape=21, size=2.85, fill='red', alpha=1, stroke=0.1) + 
    geom_treescale(color="white", y = 0.75) + 
    # geom_text(aes(label=node), hjust=-.3) + 
    xlim(-0.005, 0.565) + 
    theme_tree2(plot.margin=margin(5, 5, 20, 5), text = element_text(size = 11.75), 
        line = element_line(size = 0.4), axis.ticks.length=unit(.125, "cm"))
    plot(p)
    dev.off()

    # Carpel
    png(height = 1150, width = 1100, pointsize = 100, res = 325, file = file.path(out_dir, "output", "plots", "lnc_carpel_tr.png"))
    p <- ggtree(carpel_tr, ladderize = FALSE, size=0.55) + 
    geom_tiplab() + 
    geom_nodepoint(shape=21, size=2.8, fill='white', alpha=1, stroke=0.55) + 
    geom_point2(aes(subset=(node==22)), shape=21, size=2.85, fill='#ffbc00', alpha=1, stroke=0.1) + 
    geom_point2(aes(subset=(node==14)), shape=21, size=2.85, fill='#7e0000', alpha=1, stroke=0.1) + 
    geom_point2(aes(subset=(node==16)), shape=21, size=2.85, fill='red', alpha=1, stroke=0.1) + 
    geom_point2(aes(subset=(node==19)), shape=21, size=2.85, fill='red', alpha=1, stroke=0.1) + 
    geom_point2(aes(subset=(node==20)), shape=21, size=2.85, fill='#7e0000', alpha=1, stroke=0.1) + 
    geom_treescale(color="white", y = 0.75) + 
    # geom_text(aes(label=node), hjust=-.3) + 
    xlim(-0.005, 0.565) + 
    theme_tree2(plot.margin=margin(5, 5, 20, 5), text = element_text(size = 11.75), 
        line = element_line(size = 0.4), axis.ticks.length=unit(.125, "cm"))
    plot(p)
    dev.off()

    # Pollen
    png(height = 1150, width = 1100, pointsize = 100, res = 325, file = file.path(out_dir, "output", "plots", "lnc_pollen_tr.png"))
    p <- ggtree(pollen_tr, ladderize = FALSE, size=0.55) + 
    geom_tiplab() + 
    geom_nodepoint(shape=21, size=2.8, fill='white', alpha=1, stroke=0.55) + 
    geom_point2(aes(subset=(node==15)), shape=21, size=2.85, fill='red', alpha=1, stroke=0.1) + 
    geom_point2(aes(subset=(node==17)), shape=21, size=2.85, fill='red', alpha=1, stroke=0.1) + 
    geom_point2(aes(subset=(node==18)), shape=21, size=2.85, fill='red', alpha=1, stroke=0.1) + 
    geom_point2(aes(subset=(node==20)), shape=21, size=2.85, fill='#ffbc00', alpha=1, stroke=0.1) + 
    geom_point2(aes(subset=(node==22)), shape=21, size=2.85, fill='red', alpha=1, stroke=0.1) + 
    geom_treescale(color="white", y = 0.75) + 
    # geom_text(aes(label=node), hjust=-.3) + 
    xlim(-0.005, 0.565) + 
    theme_tree2(plot.margin=margin(5, 5, 20, 5), text = element_text(size = 11.75), 
        line = element_line(size = 0.4), axis.ticks.length=unit(.125, "cm"))
    p2 <- rotate(p, 15) %>% rotate(14)
    plot(p2)
    dev.off()

    }


    
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
    phyloBSTreeL$organ <- factor(phyloBSTreeL$organ, levels = unique(phyloBSTreeL$organ))

    sfname <- sprintf('%s.txt', paste("Brass", transcripttype, "total_tree_length", coefficient, expr_estimation, sep="_"))

    write.table(phyloBSTreeL, file=file.path(out_dir, "output", "data", sfname), 
            sep="\t", col.names=TRUE, row.names=FALSE, dec=".", quote = FALSE)


}








