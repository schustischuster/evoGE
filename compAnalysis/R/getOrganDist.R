# Compute intra-species inter-organ distances based on angiospern and mammalian ortholog gene expression data
# Thresholds: DevSeq 0.05 ERCC; Brawand 0.5 TPM (no ERCC spike-ins available)
# Data input: Brawand and DevSeq TPM and VST count expression tables of all samples


#-------------------------------------- Read data tables ---------------------------------------


getOrganDist <- function(expr_estimation = c("TPM", "counts"), 
    coefficient = c("pearson", "spearman")) {

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


    # Show startup message
    message("Reading data...")


	# Set expression input file
    if (is.element("TPM", expr_estimation)) {
        brExpr = file.path(in_dir, "Expression_data", "Brawand_inter_tpm_mat_deseq_sample_names_0_5_threshold.csv")
        dsExpr = file.path(in_dir, "Expression_data", "AT_core_inter_tpm_mat_deseq_sample_names.csv")

    } else if (is.element("counts", expr_estimation)) {
        brExpr = file.path(in_dir, "Expression_data", "Brawand_inter_count_mat_vsd_sample_names_0_5_threshold.csv")
        dsExpr = file.path(in_dir, "Expression_data", "AT_core_inter_count_mat_vsd_sample_names.csv")
    }


	# Read expression data
	brExpr <- read.table(brExpr, sep=";", dec=".", header=TRUE, stringsAsFactors=FALSE)
	dsExpr <- read.table(dsExpr, sep=";", dec=".", header=TRUE, stringsAsFactors=FALSE)


    # Stop function here to allow specific analysis of a single data set
    # return_list <- list("expr_estimation" = expr_estimation, "brExpr" = brExpr, "dsExpr" = dsExpr, "coefficient" = coefficient)
    # return(return_list)
    # }
    # return_objects <- getOrganDist(expr_estimation="TPM", coefficient="pearson") # read in expression data
    # list2env(return_objects, envir = .GlobalEnv)


    # Create "plots" folder in /out_dir/output/plots
    if (!dir.exists(file.path(out_dir, "output", "plots"))) 
        dir.create(file.path(out_dir, "output", "plots"), recursive = TRUE)

    # Show message
    message("Starting analysis...")

    
    # Update column names
    ds_col_names <- rep(c("Root", "Hypocotyl", "Leaf", "veg_apex", "inf_apex", "Flower", 
        "Stamen", "Carpel", "Pollen"), each=21)
    ds_replicate_tag_samples <- rep(c(".1",".2",".3"), times=9)
    ds_col_names <- paste0(ds_col_names, ds_replicate_tag_samples)
    ds_spec_names <- rep(c("_AT", "_AL", "_CR", "_ES", "_TH", "_MT", "_BD"), each=3)
    ds_spec_names <- rep(ds_spec_names, times=9)
    ds_col_names <- paste0(ds_col_names, ds_spec_names)
    ds_col_names <- c("gene_id", ds_col_names)

    # set column names
    colnames(dsExpr) <- ds_col_names
    

    # Generate a sequence to replace missing gene_id column
    # Remove this in case input table will have gene_id column again
    ID_repl_Br <- as.data.frame(seq(1:nrow(brExpr)))
    colnames(ID_repl_Br) <- "gene_id"
    brExpr <- cbind(ID_repl_Br, brExpr)


    # Log-transform data if TPM and Pearson correlation is chosen
    if (is.element("pearson", coefficient) && is.element("TPM", expr_estimation)) {

        brExpr[,2:ncol(brExpr)] <- log2(brExpr[,2:ncol(brExpr)] + 1)
        dsExpr[,2:ncol(dsExpr)] <- log2(dsExpr[,2:ncol(dsExpr)] + 1)
    }


    # Merge Brawand replicates
    # Need to do this manually because different number of replicates across organs and species
    brExpr <- data.frame(cbind("gene_id"=brExpr[,1], rowMeans(brExpr[,2:5]), rowMeans(brExpr[,6:8]), 
        rowMeans(brExpr[,9:14]), rowMeans(brExpr[,15:16]), rowMeans(brExpr[,17:18]), rowMeans(brExpr[,19:21]), 
        rowMeans(brExpr[,22:24]), rowMeans(brExpr[,25:26]), rowMeans(brExpr[,27:28]), rowMeans(brExpr[,29:30]), 
        rowMeans(brExpr[,31:32]), rowMeans(brExpr[,33:34]), brExpr[,35], rowMeans(brExpr[,36:37]), 
        rowMeans(brExpr[,38:40]), rowMeans(brExpr[,41:42]), rowMeans(brExpr[,43:44]), rowMeans(brExpr[,45:46]), 
        rowMeans(brExpr[,47:48]), rowMeans(brExpr[,49:50]), rowMeans(brExpr[,51:52]), rowMeans(brExpr[,53:54]), 
        rowMeans(brExpr[,55:57]), rowMeans(brExpr[,58:59]), rowMeans(brExpr[,60:61]), rowMeans(brExpr[,62:63]), 
        rowMeans(brExpr[,64:65]), rowMeans(brExpr[,66:67]), rowMeans(brExpr[,68:69]), rowMeans(brExpr[,70:71]), 
        rowMeans(brExpr[,72:74]), brExpr[,75], rowMeans(brExpr[,76:77]), rowMeans(brExpr[,78:79]), 
        rowMeans(brExpr[,80:81]), rowMeans(brExpr[,82:83]), rowMeans(brExpr[,84:85]), rowMeans(brExpr[,86:87]), 
        rowMeans(brExpr[,88:90]), rowMeans(brExpr[,91:92]), rowMeans(brExpr[,93:94]), brExpr[,95:97],  
        rowMeans(brExpr[,98:99]), rowMeans(brExpr[,100:101]), rowMeans(brExpr[,102:103])))

    tetra_organs <- rep(c("br", "cb", "ht", "kd", "lv", "ts"), each=8)
    tetra_species <- rep(c("Hsa", "Ppa", "Ptr", "Ggo", "Ppy", "Mml", "Mmu", "Mdo"), times=6)
    tetra_colnames <- paste(tetra_organs, tetra_species, sep="_")
    tetra_colnames <- tetra_colnames[-45]

    colnames(brExpr)[2:ncol(brExpr)] <- tetra_colnames


    brExpr[is.na(brExpr)] <- 0 # replaces NAs by 0
    dsExpr[is.na(dsExpr)] <- 0 # replaces NAs by 0

    # Remove ERCC spike-ins from data
    dsExpr <- dsExpr[!grepl("ERCC", dsExpr$gene_id),]



    # Calculate average expression for DevSeq replicates
    calculateAvgExpr <- function(df) {

        # Split data frame by sample replicates into a list
        # then get rowMeans for each subset and bind averaged data to gene_id column

        averaged_replicates <- do.call(cbind, lapply(split.default(df[2:ncol(df)], 
                rep(seq_along(df), 
                each = 3, 
                length.out=ncol(df)-1)
                ), rowMeans)
              )

        averaged_replicates <- cbind(df[1], averaged_replicates)
        
        return(averaged_replicates)
    }


    dsExpr <- calculateAvgExpr(dsExpr)

    DevSeq_col_names <- rep(c("Root", "Hypocotyl", "Leaf", "veg_apex", "inf_apex", "Flower", 
            "Stamen", "Carpel", "Pollen"), each=7)
    DevSeq_spec_names <- rep(c("_AT", "_AL", "_CR", "_ES", "_TH", "_MT", "_BD"), times=9)
    repl_names <- paste0(DevSeq_col_names, DevSeq_spec_names)

    colnames(dsExpr)[2:ncol(dsExpr)] <- repl_names

    # Remove pollen samples from angiosperm data
    dsExpr <- dsExpr %>% select(-c(Pollen_AT, Pollen_AL, Pollen_CR, Pollen_ES, Pollen_TH, 
                 Pollen_MT, Pollen_BD))


    # Compute metric inter-organ distances
    getOrganDist <- function(df, Species, OrganSet) {

        x_df <- select(df, contains(Species))

        if (OrganSet == "selBr") {

            x_df.1 <- select(x_df, contains("br"))
            x_df.2 <- select(x_df, contains("ht"))
            x_df.3 <- select(x_df, contains("kd"))
            x_df.4 <- select(x_df, contains("lv"))
            x_df.5 <- select(x_df, contains("ts"))

            x_df <- cbind(x_df.1, x_df.2, x_df.3, x_df.4, x_df.5)

        } else if (OrganSet == "selDS") {

            x_df.1 <- select(x_df, contains("Root"))
            x_df.2 <- select(x_df, contains("Hypocotyl"))
            x_df.3 <- select(x_df, contains("Leaf"))
            x_df.4 <- select(x_df, contains("Carpel"))
            x_df.5 <- select(x_df, contains("Stamen"))

            x_df <- cbind(x_df.1, x_df.2, x_df.3, x_df.4, x_df.5)

        }

        if (coefficient == "pearson") {
            x_cor_avg <- sqrt(1/2*(1-cor(x_df, method = "pearson")))

        } else if (coefficient == "spearman") {
            x_cor_avg <- sqrt(1/2*(1-cor(x_df, method = "spearman")))
        }

        x_cor_avg <- as.data.frame(unique(as.vector(x_cor_avg)))[-1,]
        x_cor_avg <- as.data.frame(x_cor_avg)
        colnames(x_cor_avg) <- "Distance"

        spec_tag <- rep(Species, nrow(x_cor_avg))
        spec_tag <- as.data.frame(spec_tag)
        names(spec_tag) <- "Species"

        x_cor_avg <- cbind(spec_tag, x_cor_avg)

        return(x_cor_avg)
    }


    br_list <- list("Hsa", "Ppa", "Ptr", "Ggo", "Ppy", "Mml", "Mmu", "Mdo")

    br_dist <- as.data.frame(do.call(rbind, lapply(br_list, getOrganDist, df = brExpr, OrganSet = "all")))
    comp_study_br <- data.frame("Class" = rep("Mammals", nrow(br_dist)))
    br_dist <- cbind(br_dist, comp_study_br)


    ds_list <- list("AT", "AL", "CR", "ES", "TH", "MT", "BD")

    ds_dist <- as.data.frame(do.call(rbind, lapply(ds_list, getOrganDist, df = dsExpr, OrganSet = "all")))
    comp_study_ds <- data.frame("Class" = rep("Angiosperms", nrow(ds_dist)))
    ds_dist <- cbind(ds_dist, comp_study_ds)

    # Combine angiosperm and mammalian data
    dist_df <- rbind(br_dist, ds_dist)


    # Computete distances for selected organ data
    br_dist_sel <- as.data.frame(do.call(rbind, lapply(br_list, getOrganDist, df = brExpr, OrganSet = "selBr")))
    comp_study_br_sel <- data.frame("Class" = rep("Mammals", nrow(br_dist_sel)))
    br_dist_sel <- cbind(br_dist_sel, comp_study_br_sel)

    ds_dist_sel <- as.data.frame(do.call(rbind, lapply(ds_list, getOrganDist, df = dsExpr, OrganSet = "selDS")))
    comp_study_ds_sel <- data.frame("Class" = rep("Angiosperms", nrow(ds_dist_sel)))
    ds_dist_sel <- cbind(ds_dist_sel, comp_study_ds_sel)

    # Combine data
    dist_df_sel <- rbind(br_dist_sel, ds_dist_sel)



    # Get pairwise distances organ "outliers" (e.g.distance for brain-cerebellum)
    getOutlDist <- function(df, Species, OrganSet) {

        x_df <- select(df, contains(Species))

        if (OrganSet == "selBr") {

            x_df.1 <- select(x_df, contains("br"))
            x_df.2 <- select(x_df, contains("cb"))

            x_df <- cbind(x_df.1, x_df.2)

        } else if (OrganSet == "selDS") {

            x_df.1 <- select(x_df, contains("veg_apex"))
            x_df.2 <- select(x_df, contains("inf_apex"))
            x_df.3 <- select(x_df, contains("Flower"))
            x_df.4 <- select(x_df, contains("Carpel"))
            x_df.5 <- select(x_df, contains("Stamen"))

            x_df <- cbind(x_df.1, x_df.2, x_df.3, x_df.4, x_df.5)
        }

        if (coefficient == "pearson") {
            x_cor_avg <- sqrt(1/2*(1-cor(x_df, method = "pearson")))

        } else if (coefficient == "spearman") {
            x_cor_avg <- sqrt(1/2*(1-cor(x_df, method = "spearman")))
        }

        # select pairwise organ distances
        if (OrganSet == "selBr") {

            x_cor_avg <- x_cor_avg[1,2] # brain-cerebellum

        } else if (OrganSet == "selDS") {

            x_cor_avg.1 <- x_cor_avg[1,2] # apex_v-apex_i
            x_cor_avg.2 <- x_cor_avg[1,4] # apex_v-carpel
            x_cor_avg.3 <- x_cor_avg[2,4] # apex_i-carpel
            x_cor_avg.4 <- x_cor_avg[3,4] # flower-carpel
            x_cor_avg.5 <- x_cor_avg[3,5] # flower-stamen

            x_cor_avg <- cbind(x_cor_avg.1, x_cor_avg.2, x_cor_avg.3, x_cor_avg.4, x_cor_avg.5)
        }

        x_cor_avg <- data.frame("Distance" = as.vector(x_cor_avg))

        spec_tag <- rep(Species, nrow(x_cor_avg))
        spec_tag <- as.data.frame(spec_tag)
        names(spec_tag) <- "Species"

        x_cor_avg <- cbind(spec_tag, x_cor_avg)

        return(x_cor_avg)
    }


    br_dist_outl <- as.data.frame(do.call(rbind, lapply(br_list, getOutlDist, df = brExpr, OrganSet = "selBr")))
    comp_study_br_outl <- data.frame("Class" = rep("Mammals", nrow(br_dist_outl)))
    comp_organ_br_outl <- data.frame("Organ pair" = rep(c("Brain-Cereb"), 8))
    br_dist_outl <- cbind(br_dist_outl, comp_organ_br_outl, comp_study_br_outl)

    ds_dist_outl <- as.data.frame(do.call(rbind, lapply(ds_list, getOutlDist, df = dsExpr, OrganSet = "selDS")))
    comp_study_ds_outl <- data.frame("Class" = rep("Angiosperms", nrow(ds_dist_outl)))
    comp_organ_ds_outl <- data.frame("Organ pair" = rep(c("Apex_v-Apex_i", "Apex_v-Carpel", "Apex_i-Carpel", 
        "Flower-Carpel", "Flower-Stamen"), 7))
    ds_dist_outl <- cbind(ds_dist_outl, comp_organ_ds_outl, comp_study_ds_outl)

    # Combine data
    dist_df_outl <- rbind(br_dist_outl, ds_dist_outl)



    # Make distance plot of individual species
    plotDist <- function(data, data_outl) {

        fname <- sprintf('%s.jpg', paste("Intra-species_inter-organ_distances"))

        cor_colors <- c(rep(c("#4c74b0"), each=115), rep(c("#4ca130"), each=196))
        box_colors <- c(rep(c("#2082dd"), each=8), rep(c("#32af18"), each=7))
        outl_shape <- c(rep(c(8),8), rep(c(16,16,16,16,16),7))
        shape_col <- c(rep(c("red"),8), rep(c("#4ca130"),35))
        shape_size <- c(rep(c(7),8), rep(c(4.65,4.65,4.65,4.65,4.65),7))
        outl_shape2 <- c(rep(c(16),8), rep(c(16,16,16,16,16),7))
        shape_col2 <- c(rep(c("#4c74b0"),8), rep(c("#4ca130"),35))
        shape_size2 <- c(rep(c(4.55),8), rep(c(4.65,4.65,4.65,4.65,4.65),7))
        
        p <- ggplot(data=data, aes(x = Species, y = Distance)) + 
        geom_boxplot(width = 0.75, size=1.5, fatten=2, color="black", fill=box_colors, outlier.shape = NA, alpha = 0.44) + 
        geom_beeswarm(data = data, groupOnX = TRUE, priority = c("ascending"), colour=cor_colors, cex=2, size=6) + 
        geom_beeswarm(data = data_outl, groupOnX = TRUE, priority = c("ascending"), colour=shape_col, cex=2, size=shape_size, shape=outl_shape, stroke=3) + 
        geom_beeswarm(data = data_outl, groupOnX = TRUE, priority = c("ascending"), colour=shape_col2, cex=2, size=shape_size2, shape=outl_shape2, stroke=3) + 
        scale_y_continuous(expand = c(0.05, 0), labels = comma) + 
        scale_x_discrete(labels=c("Hsa" = "Hsa", "Ppa" = "Ppa", 
            "Ptr" = "Ptr", "Ggo" = "Ggo", "Ppy" = "Ppy", 
            "Mml" = "Mml", "Mmu" = "Mm", "Mdo" = "Mdo", 
            "AT" = "AT", "AL" = "AL", "CR" = "CR", 
            "ES" = "ES", "TH" = "TH", "MT" = "MT", 
            "BD" = "BD")) + 
        guides(shape = guide_legend(override.aes = list(stroke=1.5)))

        q <- p + theme_classic() + xlab("Species") + ylab("Pearson distance") + 
        theme(text=element_text(size = 16), 
            strip.text = element_text(size = 37), 
            strip.text.x = element_text(margin = margin(0.685, 0, 0.685, 0, "cm")), 
            strip.background = element_rect(colour = 'black', fill = NA, size = 2.5), 
            axis.ticks.length = unit(0.5825, "cm"), 
            axis.ticks = element_line(colour = "black", size = 1.25), 
            axis.line = element_line(colour = 'black', size = 1.25), 
            plot.margin = unit(c(1, 0.85, 1, 1),"cm"), 
            axis.title.y = element_text(size=37, margin = margin(t = 0, r = 20, b = 0, l = 10), colour="black"), 
            axis.title.x = element_text(size=37, margin = margin(t = 10, r = 0, b = 0, l = 0), colour="black"), 
            axis.text.x = element_text(size=32.5, margin = margin(t = 7.5, b = 10), colour="black"), 
            axis.text.y = element_text(size=32, angle=0, margin = margin(r = 7.5), colour="black"), 
            panel.spacing = unit(1, "cm"), 
            panel.grid.major = element_blank(),
            panel.grid.minor.x = element_blank(), 
            panel.grid.minor.y = element_blank(), 
            legend.position = "none", 
            legend.title = element_blank(), 
            legend.text = element_text(size = 32.5), 
            legend.spacing.x = unit(0, 'cm'), 
            legend.key.size = unit(1.7, "cm"), 
            legend.background=element_blank()) 

        q <- q + facet_grid(~ Class, scales = "free_x")

        ggsave(file = file.path(out_dir, "output", "plots", fname), plot = q, 
            width = 23, height = 12, dpi = 300, units = c("in"), limitsize = FALSE) 
    }

    plotDist(data = dist_df, data_outl = dist_df_outl)



    # Plot averaged distances for ammals and angiosperms
    # Prepare data
    dist_data_all <- data.frame("Data" = rep(c("All data"), nrow(dist_df)))
    dist_df_data <- cbind(dist_df, dist_data_all)
    dist_data_sel <- data.frame("Data" = rep(c("Subset"), nrow(dist_df_sel)))
    dist_df_sel_data <- cbind(dist_df_sel, dist_data_sel)
    dist_data_avg <- rbind(dist_df_data, dist_df_sel_data)


    plotAvgDist <- function(data) {

        fname <- sprintf('%s.jpg', paste("Averaged_intra-species_inter-organ_distances"))

        cor_colors <- c(rep(c("#4c74b0"), each=115), rep(c("#4ca130"), each=196), 
            rep(c("#4c74b0"), each=76), rep(c("#4ca130"), each=70))
        box_colors <- c(rep(c("#2082dd", "#32af18"), 2))
        
        p <- ggplot(data=data, aes(x = Class, y = Distance)) + 
        geom_boxplot(width = 0.75, size=1.5, fatten=2, color="black", fill=box_colors, outlier.shape = NA, alpha = 0.44) + 
        geom_beeswarm(data = data, groupOnX = TRUE, priority = c("ascending"), colour=cor_colors, cex=2, size=5.25) + 
        scale_y_continuous(expand = c(0.05, 0), labels = comma) + 
        scale_x_discrete(labels=c("Mammals" = "Mammals", "Angiosperms" = "Angiosperms", 
            "Mammals" = "Mammals", "Angiosperms" = "Angiosperms")) + 
        guides(shape = guide_legend(override.aes = list(stroke=1.5)))

        q <- p + theme_classic() + xlab("Data set") + ylab("Pearson distance") + 
        theme(text=element_text(size = 16), 
            strip.text = element_text(size = 37), 
            strip.text.x = element_text(margin = margin(0.685, 0, 0.685, 0, "cm")), 
            strip.background = element_rect(colour = 'black', fill = NA, size = 2.5), 
            axis.ticks.length = unit(0.5825, "cm"), 
            axis.ticks = element_line(colour = "black", size = 1.25), 
            axis.line = element_line(colour = 'black', size = 1.25), 
            plot.margin = unit(c(1, 1, 1, 21),"cm"), 
            axis.title.y = element_text(size=37, margin = margin(t = 0, r = 20, b = 0, l = 10), colour="black"), 
            axis.title.x = element_text(size=37, margin = margin(t = 10, r = 0, b = 0, l = 0), colour="black"), 
            axis.text.x = element_text(size=32.5, margin = margin(t = 7.5, b = 10), colour="black"), 
            axis.text.y = element_text(size=32, angle=0, margin = margin(r = 7.5), colour="black"), 
            panel.spacing = unit(1, "cm"), 
            panel.grid.major = element_blank(),
            panel.grid.minor.x = element_blank(), 
            panel.grid.minor.y = element_blank(), 
            legend.position = "none", 
            legend.title = element_blank(), 
            legend.text = element_text(size = 32.5), 
            legend.spacing.x = unit(0, 'cm'), 
            legend.key.size = unit(1.7, "cm"), 
            legend.background=element_blank()) 

        q <- q + facet_grid(~ Data, scales = "free_x")

        ggsave(file = file.path(out_dir, "output", "plots", fname), plot = q, 
            width = 23, height = 12, dpi = 300, units = c("in"), limitsize = FALSE) 
    }

    plotAvgDist(data = dist_data_avg)








}


getOrganDist(expr_estimation="TPM", coefficient="pearson")
getOrganDist(expr_estimation="TPM", coefficient="spearman")
getOrganDist(expr_estimation="counts", coefficient="pearson")
getOrganDist(expr_estimation="counts", coefficient="spearman")





