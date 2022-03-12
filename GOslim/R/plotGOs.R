# Visualize results of GO analysis 
# Input files: GO annotation files from TAIR, quantile expression gene lists of interest, 
# GOslim terms containing core orthologs that evolve at different rate than background,
# Arabidopsis thaliana DevSeq gene expression tables


plotGOs <- function(...) {


    # Set file path for input files
    GOSLIM = file.path(in_dir, "ATH_GO_GOSLIM.txt")
    GOCAT = file.path(in_dir, "TAIR_GO_slim_categories.txt")
    orthoTPM = file.path(in_dir, "AT_core_inter_tpm_mat_deseq_sample_names.csv")
    atTPM = file.path(in_dir, "AT_genes_inter_norm_tpm_mat_deseq_sample_names.csv")

    glob_q1_gop = file.path(in_dir, "global_q1_express_enrichment.csv")
    glob_q14_gop = file.path(in_dir, "global_q14_express_enrichment.csv")
    root_q1_gop = file.path(in_dir, "root_q1_express_enrichment.csv")
    emb_dev_gop = file.path(in_dir, "embryo_post_embryonic_dev_genes.csv")

    filenames <- list.files(file.path(out_dir, "output", "data"), pattern="*.txt", full.names=TRUE)
    filenames <- filenames[!filenames %in% "./GOslim/output/data/cor_bsv_traject_1000.txt"]
    names_txt <- gsub('.*/', '', filenames)
    names <- gsub("\\..*","", names_txt)


    ldf <- lapply(filenames, function(x) {
        df <- read.table(x, sep="\t", dec=".", quote = "", header=TRUE, stringsAsFactors=FALSE)
        return(df)
    })

    names(ldf) <- names
    res <- lapply(ldf, summary) # check if file input is OK

    GOSLIM <- read.table(GOSLIM, sep="\t", dec=".", quote = "", header=FALSE, skip=4, fill = TRUE, stringsAsFactors=FALSE)
    GOCAT <- read.table(GOCAT, sep="\t", dec=".", header=TRUE, skip=7, fill = TRUE, stringsAsFactors=FALSE)
    orthoTPM <- read.table(orthoTPM, sep=";", dec=".", header=TRUE, stringsAsFactors=FALSE)
    atTPM <- read.table(atTPM, sep=";", dec=".", header=TRUE, stringsAsFactors=FALSE)

    glob_q1_go <- read.table(glob_q1_gop, sep=",", dec=".", header=TRUE, fill = TRUE, stringsAsFactors=FALSE)
    glob_q14_go <- read.table(glob_q14_gop, sep=",", dec=".", header=TRUE, fill = TRUE, stringsAsFactors=FALSE)
    root_q1_go <- read.table(root_q1_gop, sep=",", dec=".", header=TRUE, fill = TRUE, stringsAsFactors=FALSE)
    emb_dev_go <- read.table(emb_dev_gop, sep=",", dec=".", header=TRUE, fill = TRUE, stringsAsFactors=FALSE)


    # return_list <- list("ldf" = ldf, "glob_q1_go" = glob_q1_go, "glob_q14_go" = glob_q14_go, "root_q1_go" = root_q1_go, "emb_dev_go" = emb_dev_go, "GOSLIM" = GOSLIM, "GOCAT" = GOCAT, "orthoTPM" = orthoTPM, "atTPM" = atTPM)
    # return(return_list)
    # }
    # return_objects <- plotGOs()
    # list2env(return_objects, envir = .GlobalEnv)

    # Show message
    message("Starting analysis...")

    list2env(ldf, envir = .GlobalEnv)



    #--------- Extract GOslim data for DevSeq core orthologs and check GO enrichment ----------


    selCAT <- dplyr::filter(GOCAT, !grepl("cellular component", ONTOLOGY.ASPECT))


    slim_names <- selCAT$SLIM_NAME
    slim_name_ls <- setNames(as.list(c(slim_names)), c(slim_names)) # generate named list with all GOSLIM terms

    # Extract all genes for each GOslim term of slim_name_ls list
    # V1 in GOSLIM table contains the gene ID, V9 contains the GOslim terms
    slim_genes_ls <- lapply(slim_name_ls, function(x){dplyr::filter(GOSLIM, grepl(x,V9))})

    # Get list with items reduced to unique gene id's
    slim_genes_uls <- lapply(slim_genes_ls, function(x){x[!duplicated(x[,1]),]})

    coreOrthologs <- sub("\\:.*", "", orthoTPM[,1]) # get AT orthologs
    coreOrthologs <- coreOrthologs[!grepl("ERCC", coreOrthologs)] # Rm spike-ins from ortholog list
    coreOrthologs <- as.data.frame(coreOrthologs)

    # Prepare A.thaliana expression data
    atExpr <- tibble::rownames_to_column(atTPM, "gene_id")

    # Remove "other" categories because not informative
    # Remove "response to light stimulus" because plants were grown at constant light
    # Remove "cell cycle" because those genes have strong, cell-type specific expression that makes it impossible to match controls
    slim_genes_uls <- slim_genes_uls %>% purrr::list_modify("other metabolic processes" = NULL)
    slim_genes_uls <- slim_genes_uls %>% purrr::list_modify("other cellular processes" = NULL)
    slim_genes_uls <- slim_genes_uls %>% purrr::list_modify("other biological processes" = NULL)
    slim_genes_uls <- slim_genes_uls %>% purrr::list_modify("response to light stimulus" = NULL)
    slim_genes_uls <- slim_genes_uls %>% purrr::list_modify("cell cycle" = NULL)
    # Delete "Cell Growth" GOslim term from "Growth" category
    slim_genes_uls[["growth"]] <- dplyr::filter(slim_genes_uls[["growth"]], !grepl("cell growth", V9))
    # Correct GOslim term for class "RNA binding" (it contains another term called "translation factor activity, RNA binding")
    slim_genes_uls[["RNA binding"]]$V9 <- rep("RNA binding", nrow(slim_genes_uls[["RNA binding"]]))


    # Get list of orthologous genes that are associated with a GOSLIM term
    slim_ortho_ls <- lapply(slim_genes_uls, function(x){dplyr::filter(x, (V1 %in% coreOrthologs[,1]))})

    # Create stats table
    stats_ortho_genes <- lapply(slim_ortho_ls, function(x){nrow(x)})

    stats_ortho_genes_df <- data.frame(ortho_genes = matrix(unlist(stats_ortho_genes), byrow=TRUE), stringsAsFactors=FALSE)
    stats_ortho_genes_df$goslim_term <- names(slim_ortho_ls)


    # Add number of genes to each GOSLIM category of stats tables
    perm_stats_biological_process <- merge(subset(perm_stats_biological_process, 
        perm_stats_biological_process$p_value_FDR < 0.005), stats_ortho_genes_df)
    perm_stats_molecular_function <- merge(subset(perm_stats_molecular_function, 
        perm_stats_molecular_function$p_value_FDR < 0.005), stats_ortho_genes_df)
    wilcox_stats_biological_process <- merge(subset(wilcox_stats_biological_process, 
        wilcox_stats_biological_process$p_value_FDR < 0.005), stats_ortho_genes_df)
    wilcox_stats_molecular_function <- merge(subset(wilcox_stats_molecular_function, 
        wilcox_stats_molecular_function$p_value_FDR < 0.005), stats_ortho_genes_df)


    # Check if GE of functional group shows stronger or weaker conservation compared to control
    getGoslimStats_biological_process$diff <- getGoslimStats_biological_process$nlm_slope - 
                                              getGoslimStats_biological_process$nlm_slope_control

    getGoslimStats_molecular_function$diff <- getGoslimStats_molecular_function$nlm_slope - 
                                              getGoslimStats_molecular_function$nlm_slope_control

    bp_slp <- split(getGoslimStats_biological_process, getGoslimStats_biological_process$goslim_term)
    mf_slp <- split(getGoslimStats_molecular_function, getGoslimStats_molecular_function$goslim_term)

    n_cont_bp <- unlist(lapply(bp_slp, function(x){length(unique(x$nlm_slope_control))/length(unique(x$organ))}))
    n_cont_mf <- unlist(lapply(mf_slp, function(x){length(unique(x$nlm_slope_control))/length(unique(x$organ))}))

    bp_slp <- unlist(lapply(bp_slp, function(x){mean(x$diff)}))
    mf_slp <- unlist(lapply(mf_slp, function(x){mean(x$diff)}))

    bp_slp <- data.frame(goslim_term=names(bp_slp), delta_slope=bp_slp, n_cont=n_cont_bp)
    perm_stats_biological_process <- merge(perm_stats_biological_process, bp_slp)
    wilcox_stats_biological_process <- merge(wilcox_stats_biological_process, bp_slp)

    mf_slp <- data.frame(goslim_term=names(mf_slp), delta_slope=mf_slp, n_cont=n_cont_mf)
    perm_stats_molecular_function <- merge(perm_stats_molecular_function, mf_slp)
    wilcox_stats_molecular_function <- merge(wilcox_stats_molecular_function, mf_slp)

    perm_stats <- rbind(perm_stats_biological_process, perm_stats_molecular_function)
    wilcox_stats <- rbind(wilcox_stats_biological_process, wilcox_stats_molecular_function)

    # Remove GO categories that have only one background control set from wilcox_stats table
    wilcox_stats <- wilcox_stats[!(wilcox_stats$goslim_term %in% c("protein binding" ,"RNA binding")), ]



    plotGOCat <- function(df) {

        fname <- sprintf('%s.jpg', paste(deparse(substitute(df)), sep="_"))

        if (deparse(substitute(df)) == "glob_q1_go"){

            df <- subset(df, Enrichment.FDR < 0.01 & nGenes >= 8)
            plot_mar <- c(0.2, -0.025, 9.677, 0.825)
            enr_br <- c(3, 6, 9, 12)
            enr_l <- c("3", "6", "9", "12")
            fdr_br <- c(5, 15, 25)
            fdr_l <- c("5", "15", "25")
            leg_b <- "horizontal"
            plt_title <- "GO enrichment global q1"
            

        } else if (deparse(substitute(df)) == "glob_q14_go"){
            
            df <- subset(df, Enrichment.FDR < 0.000000001 & nGenes >= 8 & Fold.Enrichment >= 5)
            plot_mar <- c(1.125, -0.05, 0, -0.025)
            enr_br <- c(5.01, 6, 7, 8, 9, 10)
            enr_l <- c("5", "6", "7", "8", "9", "10")
            fdr_br <- c(10, 20, 30)
            fdr_l <- c("10", "20", "30")
            leg_b <- "vertical"
            plt_title <- "GO enrichment global q14"
            # Remove redundant GO categories
            go_rm_ls <- c("Peptide biosynthetic process ", "Ribonucleotide complex subunit organization ", 
            "Ribonucleoprotein complex subunit organization ")
            df <- df[ ! df$Pathway %in% go_rm_ls, ]
            df$Pathway <- gsub("compounds", "comp", df$Pathway)

        } else if (deparse(substitute(df)) == "root_q1_go"){

            df <- subset(df, Enrichment.FDR < 0.01 & nGenes >= 6 & Fold.Enrichment >= 3)
            plot_mar <- c(0.2, -0.05, 3.5575, 0.6925)
            enr_br <- c(3.06, 6, 9, 12)
            enr_l <- c("3", "6", "9", "12")
            fdr_br <- c(5, 10, 15)
            fdr_l <- c("5", "10", "15")
            leg_b <- "vertical"
            plt_title <- "GO enrichment root q1"
            # Remove redundant GO categories
            go_rm_ls <- c("Anatomical structure formation involved in morphogenesis ", "Cellular component assembly involved in morphogenesis ", 
            "Mitochondrial mRNA modification ", "Pollen exine formation ", "Anatomical structure morphogenesis ", 
            "Cellular component morphogenesis ")
            df <- df[ ! df$Pathway %in% go_rm_ls, ]
            df$Pathway <- gsub("transport in photosystem I", "transport in photos I", df$Pathway)
        }

        # Reorder df by number of genes
        df$Pathway <- reorder(df$Pathway, df$nGenes)
        df$FDR <- -log(df$Enrichment.FDR, 10)

        p <- ggplot(df, aes(x = nGenes, y = Pathway, colour = FDR)) +
        geom_point(mapping=aes(size = Fold.Enrichment)) + 
        scale_x_continuous(expand = c(0.1, 0)) + 
        scale_colour_continuous(low = "#e7e700", high = "red", breaks = fdr_br, labels = fdr_l, 
            name = expression(-log[10]*"(FDR)")) + 
        scale_size_continuous(range = c(5, 10.55), breaks = enr_br, labels = enr_l, 
            name = "Enrichment") + 
        guides(size = guide_legend(order = 2, keywidth = 1.9, keyheight = 1.9), colour=guide_colourbar(order = 1)) + 
        labs(x = "Number of Genes", y = NULL) + 
        ggtitle(plt_title) + 
        theme(panel.background = element_blank(), 
            axis.ticks.length = unit(0.25, "cm"), 
            axis.ticks = element_line(colour = "black", size = 1.25), 
            axis.line = element_line(colour = 'black', size = 1.25), 
            plot.margin = unit(plot_mar, "cm"), 
            plot.title = element_text(size=22.75, margin = margin(t = 0, r = 0, b = 9, l = 0), hjust = 0.5),
            legend.text=element_text(size=17.5), 
            legend.title=element_text(size=18.0),
            legend.direction = "vertical", 
            legend.box = leg_b,
            legend.key = element_rect(colour = "transparent", fill = "white"),
            axis.title.y = element_text(size=22.75, margin = margin(t = 0, r = 7.0, b = 0, l = 10), 
                colour="black", face = "bold"), 
            axis.title.x = element_text(size=22.75, margin = margin(t = 0.5, r = 0, b = 8.15, l = 0), 
                colour="black", face = "bold"), 
            axis.text.x = element_text(size=18.8, margin = margin(t = 3.5, b = 7), colour="grey20"), 
            axis.text.y = element_text(size=19.0, angle=0, margin = margin(l = 10, r = -2), colour="grey5")
        )

        ggsave(file = file.path(out_dir, "output", "plots", fname), plot = p, 
               width = 11.5, height = 6.5, dpi = 300, units = c("in"), limitsize = FALSE)

    }

    plotGOCat(glob_q1_go)
    plotGOCat(glob_q14_go)
    plotGOCat(root_q1_go)



    # Visualize GOslim categories that show different rates of GE divergence than background
    # Plot results of permutation test for SI
    plotGOslimSts <- function(df) {

        fname <- sprintf('%s.jpg', paste(deparse(substitute(df)), sep="_"))

        # Capitalize first string of each GO term category name
        CapStr <- function(y) {
            c <- strsplit(y, " ")[[1]]
            paste(toupper(substring(c, 1,1)), substring(c, 2),
                sep="", collapse=" ")
        }

        df$goslim_term <- paste(sapply(gsub("([A-Za-z]+).*", "\\1", df$goslim_term), CapStr), 
            sub(".*? ", "", df$goslim_term))

        df$goslim_term <- gsub("Growth growth", "Growth", df$goslim_term)

        # Reorder df by number of genes
        df$goslim_term <- gsub("nucleobase-containing", "nucleobase", df$goslim_term)
        df$goslim_term <- gsub("Embryo development", "Embryo and post-embryonic development", df$goslim_term)
        df$goslim_term <- reorder(df$goslim_term, desc(df$p_value_FDR))
        df$p_value_FDR <- -log(df$p_value_FDR, 10)
        df$color = ifelse(df$delta_slope < 0, "#4b71ae", "#dc580c")

        p <- ggplot(df, aes(x = p_value_FDR, y = goslim_term, colour = color)) +
        geom_point(mapping=aes(size = ortho_genes, colour = color)) + 
        scale_x_continuous(expand = c(0, 0), limits = c(1.924, 3.59), breaks = c(2, 2.5, 3, 3.5), 
            labels = c("2", "", "3", "")) + 
        scale_color_identity(breaks = c("#dc580c", "#4b71ae"), labels = c("High", "Low"), 
            guide = "legend", name = "Divergence") + 
        scale_size_continuous(range = c(5, 10.55), breaks = c(415, 900, 1600), 
            labels = c("415", "900", "1600"), name = "Gene count") + 
        guides(size = guide_legend(order = 2), colour = guide_legend(order = 1, override.aes=list(size = 8))) + 
        labs(x = expression(-log[10]*"(FDR)"), y = NULL) + 
        ggtitle("Expression divergence across functional groups") + 
        theme(panel.background = element_blank(), 
            axis.ticks.length = unit(0.25, "cm"), 
            axis.ticks = element_line(colour = "black", size = 1.25), 
            axis.line = element_line(colour = 'black', size = 1.25), 
            plot.margin = unit(c(7.25, 0.15, 0, 1.37), "cm"), 
            plot.title = element_text(size=22.75, margin = margin(t = 0, r = 0, b = 9, l = 0), hjust = 0.62),
            legend.text=element_text(size=17.5), 
            legend.title=element_text(size=18.0),
            legend.direction = "vertical", 
            legend.box = "vertical",
            legend.key.size = unit(1.75,"line"),
            legend.key = element_rect(colour = "transparent", fill = "white"),
            legend.margin=margin(t = 0.228, l = 0.41, b = -0.25, unit='cm'),
            axis.title.y = element_text(size=22.75, margin = margin(t = 0, r = 7.0, b = 0, l = 10), 
                colour="black", face = "bold"), 
            axis.title.x = element_text(size=22.75, margin = margin(t = 0, r = 0, b = 1, l = 0), 
                colour="black", face = "bold"), 
            axis.text.x = element_text(size=18.8, margin = margin(t = 3.5, b = 7.15), colour="grey20"), 
            axis.text.y = element_text(size=19.0, angle=0, margin = margin(l = 8.8, r = 3.25), colour="grey20")
        )

        ggsave(file = file.path(out_dir, "output", "plots", fname), plot = p, 
               width = 11.5, height = 6.5, dpi = 300, units = c("in"), limitsize = FALSE)

    }

    plotGOslimSts(perm_stats)



    # Plot results of GOslim rank sum test for Fig4
    plotGOslimStsM <- function(df) {

        fname <- sprintf('%s.jpg', paste(deparse(substitute(df)), sep="_"))

        # Capitalize first string of each GO term category name
        CapStr <- function(y) {
            c <- strsplit(y, " ")[[1]]
            paste(toupper(substring(c, 1,1)), substring(c, 2),
                sep="", collapse=" ")
        }

        df$goslim_term <- paste(sapply(gsub("([A-Za-z]+).*", "\\1", df$goslim_term), CapStr), 
            sub(".*? ", "", df$goslim_term))

        df$goslim_term <- gsub("Growth growth", "Growth", df$goslim_term)

        # Reorder df by number of genes
        df$goslim_term <- gsub("nucleobase-containing", "nucleobase", df$goslim_term)
        df$goslim_term <- gsub("Embryo development", "Embryo and post-embryonic development", df$goslim_term)
        df$goslim_term <- reorder(df$goslim_term, desc(df$p_value_FDR))
        df$p_value_FDR <- -log(df$p_value_FDR, 10)
        df$color = ifelse(df$delta_slope < 0, "#4b71ae", "#dc580c")

        # Get number of control gene sets
        if (nrow(df) == 9) {
        cont_df <- data.frame(y = df$goslim_term, n_cont = df$n_cont, 
            x = df$p_value_FDR + 0.35, col = rep("#888888", nrow(df)))
        }

        p <- ggplot(df, aes(x = p_value_FDR, y = goslim_term, colour = color)) +
        geom_point(mapping=aes(size = ortho_genes, colour = color)) + 
        scale_x_continuous(expand = c(0.05, 0), limits = c(2, 5.025), breaks = c(2, 2.5, 3, 3.5, 4, 4.5, 5), 
            labels = c("2", "", "3", "", "4", "", "5")) + 
        scale_y_discrete(expand = c(0.075, 0)) + 
        scale_color_identity(breaks = c("#dc580c", "#4b71ae"), labels = c("High", "Low"), 
            guide = "legend", name = "Divergence") + 
        scale_size_continuous(range = c(5.9, 12.0), breaks = c(415, 700, 1100, 1600), 
            labels = c("415", "700", "1100", "1600"), name = "Gene count") + 
        guides(size = guide_legend(order = 2), colour = guide_legend(order = 1, override.aes=list(size = 8))) + 
        labs(x = expression(-log[10]*"(FDR)"), y = NULL) + 
        ggtitle("Gene expression divergence \nacross functional groups") + 
        theme(panel.background = element_blank(), 
            axis.ticks.length = unit(0.26, "cm"), 
            axis.ticks = element_line(colour = "black", size = 1.025), 
            axis.line = element_line(colour = 'black', size = 1.025), 
            plot.margin = unit(c(1.1225, 1.75, 0, 2.25), "cm"), 
            plot.title = element_text(size=21.25, margin = margin(t = 8, r = 0, b = 11, l = 0), 
                hjust = 0.5, lineheight = 1),
            legend.position = c(0.815, 0.42), 
            legend.text = element_text(size=19.25), 
            legend.title = element_text(size=19.75),
            legend.direction = "vertical", 
            legend.box = "vertical",
            legend.key.size = unit(1.75,"line"),
            legend.key = element_rect(colour = "transparent", fill = "white"),
            legend.margin=margin(t = 0.228, l = 0.41, b = 0, unit='cm'),
            axis.title.y = element_text(size=21.25, margin = margin(t = 0, r = 7.0, b = 0, l = 10), 
                colour="black", face = "bold"), 
            axis.title.x = element_text(size=21.25, margin = margin(t = 0, r = 0, b = 1, l = 0), 
                colour="black", face = "bold"), 
            axis.text.x = element_text(size=18.95, margin = margin(t = 4.25, b = 7.15), colour="grey20"), 
            axis.text.y = element_text(size=20.28, angle=0, margin = margin(l = 8.8, r = 4), colour="grey20")
        )

        ggsave(file = file.path(out_dir, "output", "plots", fname), plot = p, 
               width = 11.5, height = 6.5, dpi = 300, units = c("in"), limitsize = FALSE)

    }

    plotGOslimStsM(wilcox_stats)




    # Heatmap showing expression of embryo development genes in Arabidopsis thaliana
    at_embr_dev <- atExpr[atExpr$gene_id %in% slim_ortho_ls[["embryo development"]]$V1,]


    calculateAvgEx <- function(df) {

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

    at_embr_dev <- calculateAvgEx(at_embr_dev)


    # Scale data to the unit interval
    scaleTPM <- function(x){(x-min(x))/(max(x)-min(x))}
    at_embr_dev_scd <- as.data.frame(t(apply(at_embr_dev[,2:ncol(at_embr_dev)], 1, scaleTPM)))


    color.palette <- function(steps, n.steps.between=NULL, ...) {

        if (is.null(n.steps.between)) 
            n.steps.between <- rep(0, (length(steps)-1))

        if (length(n.steps.between) != length(steps)-1)
            stop("Must have one less n.steps.between value than steps")

        fill.steps <- cumsum(rep(1, length(steps)) + c(0,n.steps.between))
        RGB <- matrix(NA, nrow = 3, ncol = fill.steps[length(fill.steps)])
        RGB[,fill.steps] <- col2rgb(steps)

        for (i in which(n.steps.between > 0)) {
            col.start = RGB[,fill.steps[i]]
            col.end = RGB[,fill.steps[i + 1]]

            for (j in seq(3)) {
                vals <-seq(col.start[j], col.end[j], length.out=n.steps.between[i]+2)[2:(2+n.steps.between[i]-1)]  
                RGB[j,(fill.steps[i] + 1):(fill.steps[i + 1] - 1)] <- vals
            }
        }

        new.steps <- rgb(RGB[1, ], RGB[2, ], RGB[3, ], maxColorValue = 255)
        pal <- colorRampPalette(new.steps, ...)

        return(pal)
    }

    # Define colors and number of steps for the plot
    steps <- c("#fae85a", "#f7ea40", "#fdc91c", "#ffa700", "#fe8300", "#f85b17", 
        "#ea2828", "#ea285a")

    pal <- color.palette(steps, c(2, 10, 11, 12, 13, 14, 5), space = "rgb")


    # Create heatmap with reversed RowSideColors
    png(height = 880, width = 1600, pointsize = 10, file = file.path(out_dir, "output", "plots", "at_embr_dev_scaled.png"))
    cc <- c(rep("#6a54a9",6), rep("#53b0db",4), rep("#2c8654",9), rep("#96ba37",3), rep("#fad819",3), rep("#e075af",4), rep("red3",8), rep("grey70",6))

    heatmap.2(as.matrix(at_embr_dev_scd), 
        density.info = "none",
        labRow = FALSE, 
        labCol = FALSE,
        dendrogram = "none", 
        col = pal(100), 
        scale = "none",
        trace = "none",
        lmat = rbind(c(0,0,0,0,0), c(0,5,0,4,0), c(0,3,0,2,0), c(0,0,0,1,0), c(0,0,0,0,0)), 
        lhei = c(0,2.5,5,0.28,0.1),
        lwid = c(0.1,2.4,0.25,5,0.5),
        key.par = list(cex = 2.8), 
        ColSideColors = cc, 
        margins = c(2, 2),
        key = TRUE,
        key.xlab = "",
        key.title = "",
        distfun = function(x) as.dist(sqrt(1/2*(1-cor(t(x))))),
        hclustfun = function(x) hclust(x, method = "average"),
        Rowv = TRUE, 
        Colv = FALSE
        )

    dev.off()



    # Plot GO categories present in embryo and post-embryonic development GOslim parent term
    plotEDCat <- function(df) {

        fname <- sprintf('%s.jpg', paste(deparse(substitute(df)), "child_terms" , sep="_"))

        # Remove some GO categories
        go_rm_ls <- c("Developmental process ", "Anatomical structure development ", 
            "Multicellular organismal process ", "Multicellular organism development ", 
            "System development ", "Post-embryonic development ", "Embryo development ", 
            "Reproductive process ", "Reproduction ", "Developmental process involved in reproduction ", 
            "Reproductive system development ", "Cellular component organization or biogenesis ", 
            "Cellular component organization ", "Organelle organization ", "RNA processing ", 
            "NcRNA metabolic process ", "Cell division ")
        df <- df[ ! df$Pathway %in% go_rm_ls, ]
        df$Pathway <- gsub("Embryo development ending in seed dormancy ", 
            "Embryo dev ending in seed dormancy ", df$Pathway)


        p <- ggplot(df, aes(x = reorder(Pathway, nGenes), y = nGenes)) +
        geom_bar(stat = "identity", width = 0.75, fill = "#37764f") + 
        coord_flip() + 
        scale_x_discrete(expand = c(0.025,0)) + 
        scale_y_continuous(expand = c(0.03,0)) + 
        labs(x = NULL, y = "Number of Genes") + 
        ggtitle("Embryo and post-embryonic development") + 
        theme(panel.background = element_blank(), 
            axis.ticks.length = unit(0.24, "cm"), 
            axis.ticks = element_line(colour = "black", size = 1.1), 
            axis.line = element_line(colour = 'black', size = 1.1), 
            plot.margin = unit(c(0.5, 0.5, 1.0, 0.5), "cm"), 
            plot.title = element_text(size = 19.85, margin = margin(t = 1, r = 0, b = 4, l = 0), hjust = 0.5),
            axis.title.y = element_text(size = 22.0, margin = margin(t = 0, r = 7.0, b = 0, l = 10), 
                colour="black", face = "bold"), 
            axis.title.x = element_text(size = 22.0, margin = margin(t = 0.5, r = 0, b = 8.15, l = 0), 
                colour="black", face = "bold"), 
            axis.text.x = element_text(size = 18.8, margin = margin(t = 3.5, b = 7), colour = "grey5"), 
            axis.text.y = element_text(size = 19.0, angle = 0, margin = margin(l = 10, r = -2), colour = "grey5")
        )

        ggsave(file = file.path(out_dir, "output", "plots", fname), plot = p, 
               width = 11.5, height = 6.5, dpi = 300, units = c("in"), limitsize = FALSE)

    }

    plotEDCat(emb_dev_go)



}


