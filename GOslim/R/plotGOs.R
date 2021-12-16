# Visualize results of GO analysis 
# Input files: GO annotation files from TAIR, quantile expression gene lists of interest, 
# GOslim terms cntaining core orthologs that evolve at different rate than background


plotGOs <- function(...) {


    # Set file path for input files
    GOSLIM = file.path(in_dir, "ATH_GO_GOSLIM.txt")
    GOCAT = file.path(in_dir, "TAIR_GO_slim_categories.txt")
    orthoTPM = file.path(in_dir, "AT_core_inter_tpm_mat_deseq_sample_names.csv")
    atTPM = file.path(in_dir, "AT_genes_inter_norm_tpm_mat_deseq_sample_names.csv")

    glob_q1_gop = file.path(in_dir, "global_q1_express_ids.csv")
    glob_q14_gop = file.path(in_dir, "global_q14_express_ids.csv")
    root_q1_gop = file.path(in_dir, "root_q1_express_ids.csv")

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


    # return_list <- list("ldf" = ldf, "glob_q1_go" = glob_q1_go, "glob_q14_go" = glob_q14_go, "root_q1_go" = root_q1_go, "GOSLIM" = GOSLIM, "GOCAT" = GOCAT, "orthoTPM" = orthoTPM, "atTPM" = atTPM)
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

    bp_slp <- unlist(lapply(bp_slp, function(x){mean(x$diff)}))
    mf_slp <- unlist(lapply(mf_slp, function(x){mean(x$diff)}))

    bp_slp <- data.frame(goslim_term=names(bp_slp), delta_slope=bp_slp)
    perm_stats_biological_process <- merge(perm_stats_biological_process, bp_slp)
    wilcox_stats_biological_process <- merge(wilcox_stats_biological_process, bp_slp)

    mf_slp <- data.frame(goslim_term=names(mf_slp), delta_slope=mf_slp)
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
            enr_br <- c(4, 8, 12)
            enr_l <- c("4", "8", "12")
            fdr_br <- c(5, 15, 25)
            fdr_l <- c("5", "15", "25")
            leg_b <- "horizontal"
            plt_title <- "GO enrichment global q1"
            

        } else if (deparse(substitute(df)) == "glob_q14_go"){
            
            df <- subset(df, Enrichment.FDR < 0.000000001 & nGenes >= 8 & Fold.Enrichment >= 5)
            plot_mar <- c(1.125, -0.05, 0, -0.025)
            enr_br <- c(6, 8, 10)
            enr_l <- c("6", "8", "10")
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
            enr_br <- c(4, 6, 8, 10)
            enr_l <- c("4", "6", "8", "10")
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
        scale_colour_continuous(low = "#eaea00", high = "red", breaks = fdr_br, labels = fdr_l, 
            name = expression(-log[10]*"(FDR)")) + 
        scale_size_continuous(range = c(5, 10.55), breaks = enr_br, labels = enr_l, 
            name = "Enrichment") + 
        guides(size = guide_legend(order = 2), colour=guide_colourbar(order = 1)) + 
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
            legend.key = element_blank(),
            axis.title.y = element_text(size=22.75, margin = margin(t = 0, r = 7.0, b = 0, l = 10), 
                colour="black", face = "bold"), 
            axis.title.x = element_text(size=22.75, margin = margin(t = 0.5, r = 0, b = 8.15, l = 0), 
                colour="black", face = "bold"), 
            axis.text.x = element_text(size=18.8, margin = margin(t = 3.5, b = 7), colour="grey20"), 
            axis.text.y = element_text(size=19.0, angle=0, margin = margin(l = 10, r = -2), colour="grey20")
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
        df$color = ifelse(df$delta_slope < 0, "blueviolet", "#e79600")

        p <- ggplot(df, aes(x = p_value_FDR, y = goslim_term, colour = color)) +
        geom_point(mapping=aes(size = ortho_genes, colour = color)) + 
        scale_x_continuous(expand = c(0.05, 0), limits = c(2, 3.4), breaks = c(2, 2.5, 3), 
            labels = c("2", "", "3")) + 
        scale_color_identity(breaks = c("#e79600", "blueviolet"), labels = c("High", "Low"), 
            guide = "legend", name = "Expression \ndivergence") + 
        scale_size_continuous(range = c(5, 10.55), breaks = c(425, 850, 1500), 
            labels = c("425", "850", "1500"), name = "Gene count") + 
        guides(size = guide_legend(order = 2), colour = guide_legend(order = 1, override.aes=list(size = 8))) + 
        labs(x = expression(-log[10]*"(FDR)"), y = NULL) + 
        ggtitle("Expression divergence across functional groups") + 
        theme(panel.background = element_blank(), 
            axis.ticks.length = unit(0.25, "cm"), 
            axis.ticks = element_line(colour = "black", size = 1.25), 
            axis.line = element_line(colour = 'black', size = 1.25), 
            plot.margin = unit(c(7.25, 0.12, 0, 1.37), "cm"), 
            plot.title = element_text(size=22.75, margin = margin(t = 0, r = 0, b = 9, l = 0), hjust = 0.62),
            legend.text=element_text(size=17.5), 
            legend.title=element_text(size=18.0),
            legend.direction = "vertical", 
            legend.box = "vertical",
            legend.key = element_blank(),
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
    # plotGOslimSts(wilcox_stats)



    # Heatmap showing expression of embryo development genes in Arabidopsis thaliana
    at_embr_dev <-atExpr[atExpr$gene_id %in% slim_ortho_ls[["embryo development"]]$V1,]

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
        "#ee2727", "#ea2828", "#ea285a")

    pal <- color.palette(steps, c(2, 10, 11, 12, 13, 14, 15, 16), space = "rgb")


    # Create heatmap with reversed RowSideColors
    png(height = 880, width = 1600, pointsize = 10, file = file.path(out_dir, "output", "plots", "at_embr_dev_scaled.png"))
    cc <- c(rep("#6a54a9",18), rep("#53b0db",12), rep("#2c8654",27), rep("#96ba37",9), rep("#fad819",9), rep("#e075af",12), rep("red3",24), rep("grey70",18))

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
        hclustfun = function(x) hclust(x, method = "complete"),
        Rowv = TRUE, 
        Colv = FALSE
        )

    dev.off()

























    # Compute GO enrichment p value and FDR corrected p value
    pwdata <- split(goslim_northo_stats, rep(1:(nrow(goslim_northo_stats)), each = 1))
    lst <- setNames(vector('list', length(pwdata)), 1:length(pwdata))

    for (i in 1:length(pwdata)) {
       allg <- as.numeric(as.character(unlist(pwdata[[i]][,2])))
       ortg <- as.numeric(as.character(unlist(pwdata[[i]][,3])))
       probabilities <- dhyper(c(0:allg), allg, (28553-allg), 7003, log = FALSE)
       pvalue <- 2*(sum(probabilities[(ortg+1):(allg+1)]))
       lst[[i]] <- pvalue
    }

    p_value <- c(do.call(rbind, lst))
    padj<- data.frame(p_value=p.adjust(p_value, method = "fdr", n = length(p_value))) # FDR correction
    goslim_northo_enrich_stats <- cbind(goslim_northo_stats, padj)


    # Remove all ortholog GOslim lists wth fewer entries than defines sample_size
    slim_ortho_ls <- Filter(function(dt) nrow(dt) >= sample_size, slim_ortho_ls)



    #------------ Combine DevSeq core ortholog expression tables with GOslim data -------------


    # Prepare angiosperm ortholog data
    orthoExpr <- data.frame(gene_id=sub("\\:.*", "", orthoTPM[,1]),orthoTPM[,2:ncol(orthoTPM)])
    orthoExpr[,2:ncol(orthoExpr)] <- log2(orthoExpr[,2:ncol(orthoExpr)] + 1)
    orthoExpr <- orthoExpr[!grepl("ERCC", orthoExpr$gene_id),]


    # Negate dplyr %in%
    `%!in%` = Negate(`%in%`)

    # Remove pollen samples
    orthoExpr <- orthoExpr %>% select (-c(
        A.thaliana_flowers_mature_pollen_28d_.2.:B.distachyon_flowers_mature_pollen_32d_.1.))


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

        x_avg <- calculateAvgExpr(orthoExpr)

        DevSeq_col_names <- rep(c("Root", "Hypocotyl", "Leaf", "veg_apex", "inf_apex", "Flower", 
            "Stamen", "Carpel"), each=7)
        DevSeq_spec_names <- rep(c("_AT", "_AL", "_CR", "_ES", "_TH", "_MT", "_BD"), times=8)
        repl_names <- paste0(DevSeq_col_names, DevSeq_spec_names)

        colnames(x_avg)[2:ncol(x_avg)] <- repl_names


        # Compute average expression and sd for each organ
        calculateSpAvg <- function(df) {

            averaged_spec <- do.call(cbind, lapply(split.default(df[2:ncol(df)], 
                rep(seq_along(df), 
                each = 7, 
                length.out=ncol(df)-1)
                ), rowMeans)
              )

            base_averaged <- rowMeans(df[2:ncol(df)])
            # base_min <- apply(df[2:ncol(df)], 1, FUN = min)
            # base_max <- apply(df[2:ncol(df)], 1, FUN = max)
            # quartiles <- as.data.frame(t(apply(df[2:ncol(df)], 1, quantile, c(0.25, 0.75))))
            # names(quartiles) <- c("q25", "q75")

            RowSD <- function(x) {
                sqrt(rowSums((x - rowMeans(x))^2)/(dim(x)[2] - 1))
            }

            averaged_sd <- do.call(cbind, lapply(split.default(df[2:ncol(df)], 
                rep(seq_along(df), 
                each = 7, 
                length.out=ncol(df)-1)
                ), RowSD)
              )

            names_averaged_spec <- unique(sub("\\_.*", "", colnames(df)[2:ncol(df)]))
            avg_names <- paste("avg", names_averaged_spec, sep="_")
            colnames(averaged_spec) <- avg_names
            # sd_names <- paste("sd", names_averaged_spec, sep="_")
            # colnames(averaged_sd) <- sd_names

            averaged_spec <- cbind(df[1], averaged_spec, base_averaged)
        
            return(averaged_spec)
        }

        spec_avg <- calculateSpAvg(x_avg)

        orthoExDf <- merge(spec_avg, x_avg)



    #--------- Preprocess control data for GOslim term analysis by 1:1/k:1 matching -----------


    # Match control genes to each ortholog GOslim class
    getGoslimStats <- as.data.frame(do.call(rbind, lapply(slim_ortho_ls, function(x){

        x_df <- as.data.frame(x)
        x_df <- data.frame(gene_id=x_df$V1, goslim=x_df$V9)
        sign <- as.numeric(rep(c(1), nrow(x_df)))
        x_df <- cbind(x_df, sign)
        control <- data.frame(gene_id=subset(coreOrthologs[,1], coreOrthologs[,1] %!in% x_df[,1]))
        control_df <- data.frame(goslim=rep(unique(x_df$goslim), nrow(control)), 
            sign=as.numeric(rep(c(0), nrow(control))))
        control_df <- cbind(control, control_df)
        comb_df <- rbind(x_df, control_df)
        comb_exdf <- merge(comb_df, orthoExDf)

        # Set ratio for control:treatment
        ntreat <- nrow(x_df)

        cratio <- round((nrow(control_df)*0.70)/ntreat)

        # Create "plots" folder in /out_dir/output/plots
        if (!dir.exists(file.path(out_dir, "output", "plots", "MatchIt"))) 
            dir.create(file.path(out_dir, "output", "plots", "MatchIt"), recursive = TRUE)

        # Create background gene set
        matchSample <- function(x) {

            success <- FALSE
            while (!success) {

                # Create background gene set
                match_res <- matchit(sign ~ base_averaged, x, method="nearest", 
                    distance="mahalanobis", replace=FALSE, m.order="data", ratio=cratio)
                match_res_m <- match_res$match.matrix

                # Extract standard mean difference from matchIt summary data
                comp <- as.data.frame(summary(match_res, standardize = TRUE)["sum.matched"])
                stmdif <- abs(comp[1,3])
                varR <- abs(comp[1,4])

                cratio <- cratio-1

                # check for success
                success <- ((stmdif <= 0.01) && (varR >= 1))
            }

            return(match_res_m)
        }

        match_res_df <- matchSample(comb_exdf)

        
        control_out <- apply(match_res_df, 2, function(x) {

            control_out <- comb_exdf[x,]
            control_out <- control_out %>% select (-c(avg_Root, avg_Hypocotyl, avg_Leaf, avg_veg, 
            avg_inf, avg_Flower, avg_Stamen, avg_Carpel, base_averaged))
            return(control_out)

        })

        goslim_out <- merge(x_df, orthoExDf)

        goslim_out <- goslim_out %>% select (-c(avg_Root, avg_Hypocotyl, avg_Leaf, avg_veg, 
            avg_inf, avg_Flower, avg_Stamen, avg_Carpel, base_averaged))

        
        # Plot results of k:1 matching
        mplot <- data.frame(do.call(cbind, lapply(control_out, function(c) { 
            tp <- unlist(c[-1:-3])
            return(tp)
        })))

        colnames(mplot) <- paste0(rep("m", ncol(mplot)), 1:ncol(mplot))

        mbplot <- data.frame(t=unlist(goslim_out[-1:-3]), mplot)
        mbplotcat <- rep(colnames(mbplot), each=nrow(mbplot))
        ggmbplot <- data.frame(exp=unlist(mbplot), class=mbplotcat)
        allc <- data.frame(exp=unlist(subset(comb_exdf, sign==0)[-1:-12]))
        allcc <- cbind(allc, data.frame(class=rep("c", nrow(allc))))
        gg2mbplot <- rbind(ggmbplot, allcc)


        plotMatchIt <- function(data) {

            fname <- sprintf('%s.png', paste(unique(x_df$goslim), "matchIt", sep="_"))

            data$class <- factor(data$class, levels = unique(data$class))

            plt_title <- paste(unique(x_df$goslim), " (n=", ntreat, ")", sep="")

            p <- ggplot(data=data, aes(x = class, y = exp)) + 
            geom_boxplot(data = data, aes(x = class, y = exp)) + 
            ggtitle(plt_title) + 
            xlab("GOslim(t) + Matched_control(m) + All_control(c)") + ylab("Expression (log2[TPM+1])")

            ggsave(file = file.path(out_dir, "output", "plots", "MatchIt", fname), plot = p, 
                width = 5.9, height = 5.9, dpi = 300, units = c("in"), limitsize = FALSE) 
        }

        plotMatchIt(data = gg2mbplot)

        
        # Make plot of "biosynthetic process" matching for supplement
        if (unique(x_df$goslim) == "biosynthetic process") {

            # Get expression of all unmatched control genes
            cids <- data.frame(do.call(rbind, lapply(control_out, function(k) k[1])))

            cgogeneexpress <- data.frame(subset(comb_exdf[,c(3,13:ncol(comb_exdf))], comb_exdf[,1] %!in% cids[,1]))
            cgeneexpress <- subset(cgogeneexpress, sign == 0)[-1]
            cgeneexpress_df <- data.frame(exp=unlist(cgeneexpress), class=rep("unmatched controls", length(unlist(cgeneexpress))))
            tcdata <- ggmbplot
            tcdata$class <- gsub("t", unique(x_df$goslim), tcdata$class)
            tcdata$class <- gsub("m", "matched controls #" , tcdata$class)
            biosplotdf <- rbind(tcdata, cgeneexpress_df)

            plotMatchItSI <- function(data) {

                fname <- sprintf('%s.png', paste(unique(x_df$goslim), "matchIt_SI", sep="_"))

                data$class <- factor(data$class, levels = unique(data$class))

                ntreatmc <- nrow(x_df)
                nuc <- nrow(cgeneexpress)

                n_text <- data.frame(x = c(1, 2, 3, 4, 5), y = c(13.15, 13.15, 13.15, 13.15, 13.15), 
                    label = c(rep(ntreatmc, 4), nuc))

                p <- ggplot(data=data, aes(x = class, y = exp)) + 
                geom_boxplot(data = data, aes(x = class, y = exp), fill=c("#28a100","grey50","grey50","grey50","#ed0000"), 
                    size=1.1, fatten=1.5, outlier.size = 2, alpha=0.35) + 
                scale_y_continuous(limits = c(0, 14.75))

                q <- p + theme_classic() + xlab("") + ylab("log2(TPM+1)") +
                geom_text(data = n_text, mapping = aes(x = x, y = y, label = label), 
                    size=6.1, hjust = 0, angle=90, col="grey55") + 
                theme(text=element_text(size = 16), 
                axis.ticks.length = unit(0.29, "cm"), 
                axis.ticks = element_line(colour = "black", size = 1.25), 
                axis.line = element_line(colour = 'black', size = 1.25), 
                plot.margin = unit(c(0.2, 0.1, 0, 0),"cm"), 
                axis.title.y = element_text(size=22.5, margin = margin(t = 0, r = 1.0, b = 0, l = 10), 
                    colour="black", face = "bold"), 
                axis.title.x = element_text(size=8, margin = margin(t = 0, r = 0, b = 0, l = 0), 
                    colour="black", face = "bold"), 
                axis.text.x = element_text(size=18.8, angle=90, margin = margin(t = 2.5, b = 0), 
                    colour="grey20", hjust=1, vjust=0.5), 
                axis.text.y = element_text(size=18.8, angle=0, margin = margin(l = 2.5, r = 1.5), 
                    colour="grey20"), 
                panel.spacing = unit(0.5, "cm"), 
                panel.grid.major = element_blank(),
                panel.grid.minor.x = element_blank(), 
                panel.grid.minor.y = element_blank(), 
                legend.position = "none")

                ggsave(file = file.path(out_dir, "output", "plots", "MatchIt", fname), plot = q, 
                    width = 3.05, height = 6.5, dpi = 300, units = c("in"), limitsize = FALSE) 
            }

            plotMatchItSI(data = biosplotdf)

        }



        #---------------- Get gene expression divergence rates for ATH/AL vs species X -----------------


        # Show message
        message("Calculating expression distances...")


        ortho_control <- lapply(control_out, function(x) {

            df <- merge(x, orthoExpr)
            df <- df[-2:-59]
            return(df)
        })

        ortho_go <- merge(goslim_out, orthoExpr)
        ortho_go <- ortho_go[-2:-59]


        # Use pearson correlation, inter-organ normalization and TPM for ms

        getDSOrganCor <- function(df, organ) {

            # Select rows for each organ
            if (organ == "Root"){
                df <- df[,2:22]
            } else if (organ == "Hypocotyl"){
                df <- df[,23:43]
            } else if (organ == "Leaf"){
                df <- df[,44:64]
            } else if (organ == "Apex_veg"){
                df <- df[,65:85]
            } else if (organ == "Apex_inf"){
                df <- df[,86:106]
            } else if (organ == "Flower"){
                df <- df[,107:127]
            } else if (organ == "Stamen"){
                df <- df[,128:148]
            } else if (organ == "Carpel"){
                df <- df[,149:169]
            }

            df_cor <- sqrt(1/2*(1 - cor(df, method="pearson")))

            replcor <- function(xrepl){

                xrepl <- c(xrepl)[xrepl>0][!(duplicated(c(xrepl)[xrepl>0]))]
                return(xrepl)
            }

            sp0_repl <- c(mean(replcor(df_cor[1:3,1:3])), mean(replcor(df_cor[4:6,4:6])), 
                mean(replcor(df_cor[7:9,7:9])), mean(replcor(df_cor[10:12,10:12])), 
                mean(replcor(df_cor[13:15,13:15])), mean(replcor(df_cor[16:18,16:18])), 
                mean(replcor(df_cor[19:21,19:21]))) # AT-AT

            # Get mean
            df_cor <- as.data.frame(df_cor, stringsAsFactors=FALSE)


            avgRepl <- function(x_df) {

                getRepl <- function(x) {

                    split.default(x, 
                        rep(seq_along(x), 
                            each = 3, 
                            length.out=ncol(x)
                            )
                        )
                }

                repl_lst <- getRepl(x_df)
                repl_sum <- lapply(repl_lst, sum)
                repl_mean <- as.numeric(unlist(repl_sum))/9

                return(repl_mean)
            }

            sp1_repl <- avgRepl(df_cor[4:6,1:3]) # AT-AL
            sp2_repl <- avgRepl(df_cor[7:9,1:6]) # AT-CR AL-CR
            sp3_repl <- avgRepl(df_cor[10:12,1:9]) # AT-ES AL-ES CR-ES
            sp4_repl <- avgRepl(df_cor[13:15,1:12]) # AT-TH AL-TH CR-TH ES-TH
            sp5_repl <- avgRepl(df_cor[16:18,1:15]) # AT-MT AL-MT CR-MT ES-MT TH-MT
            sp6_repl <- avgRepl(df_cor[19:21,1:18]) # AT-BD AL-BD CR-BD ES-BD TH-BD MT-BD

            getError <- function(cor_data) {
                std <- sd(cor_data, na.rm=TRUE)
                num <- length(cor_data)
                error <- std/sqrt(num)
                return(error)
            } # Use this function to replace sd if reqired

            df_cor_sd <- data.frame(error = c(rep(as.numeric(c(sd(sp0_repl))),length(sp0_repl)),
                    as.numeric(c(sd(sp1_repl))), rep(as.numeric(c(sd(sp2_repl))),length(sp2_repl)), 
                    rep(as.numeric(c(sd(sp3_repl))),length(sp3_repl)), rep(as.numeric(c(sd(sp4_repl))),length(sp4_repl)), 
                    rep(as.numeric(c(sd(sp5_repl))),length(sp5_repl)), rep(as.numeric(c(sd(sp6_repl))),length(sp6_repl))))

            df_cor_avg <- data.frame(correlation = c(sp0_repl, sp1_repl, sp2_repl, sp3_repl, sp4_repl, sp5_repl, sp6_repl))
            div_tag <- data.frame(clade = c(rep("T0", length(sp0_repl)), "T1", rep("T2", length(sp2_repl)), rep("T3", length(sp3_repl)), 
                rep("T4", length(sp4_repl)), rep("T5", length(sp5_repl)), rep("T6", length(sp6_repl))))
            organ_id <- data.frame(comp_organ = rep(organ, nrow(df_cor_avg)))
            div_times <- data.frame(div_times = c(rep(0, length(sp0_repl)), 7.1, rep(9.4, length(sp2_repl)), rep(25.6, length(sp3_repl)), 
                rep(46, length(sp4_repl)), rep(106, length(sp5_repl)), rep(160, length(sp6_repl))))
            dataset <- data.frame(dataset = rep("Angiosperms ", nrow(df_cor_avg)))
            df_cor_avg <- cbind(div_tag, organ_id, div_times, df_cor_avg, df_cor_sd, dataset)

            return(df_cor_avg)

        }

        root_divc <- lapply(ortho_control, getDSOrganCor, organ="Root")
        hypocotyl_divc <- lapply(ortho_control, getDSOrganCor, organ="Hypocotyl")
        leaf_divc <- lapply(ortho_control, getDSOrganCor, organ="Leaf")
        veg_apex_divc <- lapply(ortho_control, getDSOrganCor, organ="Apex_veg")
        inf_apex_divc <- lapply(ortho_control, getDSOrganCor, organ="Apex_inf")
        flower_divc <- lapply(ortho_control, getDSOrganCor, organ="Flower")
        stamen_divc <- lapply(ortho_control, getDSOrganCor, organ="Stamen")
        carpel_divc <- lapply(ortho_control, getDSOrganCor, organ="Carpel")
    

        root_divg <- getDSOrganCor(df=ortho_go, organ="Root")
        hypocotyl_divg <- getDSOrganCor(df=ortho_go, organ="Hypocotyl")
        leaf_divg <- getDSOrganCor(df=ortho_go, organ="Leaf")
        veg_apex_divg <- getDSOrganCor(df=ortho_go, organ="Apex_veg")
        inf_apex_divg <- getDSOrganCor(df=ortho_go, organ="Apex_inf")
        flower_divg <- getDSOrganCor(df=ortho_go, organ="Flower")
        stamen_divg <- getDSOrganCor(df=ortho_go, organ="Stamen")
        carpel_divg <- getDSOrganCor(df=ortho_go, organ="Carpel")

        ortho_div_rates <- rbind(root_divg, hypocotyl_divg, leaf_divg, veg_apex_divg, inf_apex_divg, 
            flower_divg, stamen_divg, carpel_divg)

        # Set up lists containing metric pearson expression distances
        ortho_organ_lst <- list(ortho_div_rates[1:28,], ortho_div_rates[29:56,], 
            ortho_div_rates[57:84,], ortho_div_rates[85:112,], ortho_div_rates[113:140,], 
            ortho_div_rates[141:168,], ortho_div_rates[169:196,], ortho_div_rates[197:224,])




        #---- Apply non-linear regression to sOU and pearson dist expression data and compare slopes -----

        # Non-linear regression using negative exponential law fit: pairwise expression differences
        # between species saturate with evolutionary time in a power law relationship
        # Fits assumption of OU model underlying stabilizing GE selection as a decelarated process

        getNLEstimates <- function(corrdata) {

            comp_organ <- unique(corrdata$comp_organ)


            nl_model <- function(a, b, c, x){

                y = a * exp(c * x) + b * (1 - exp(c * x))
                return(y)
            }
            # a + b defines maximum y value
            # a defines intercept


            x_DS_grid <- seq(0, 160, length = 200)  ## prediction grid

            cor_0 <- corrdata$correlation[corrdata$clade=="T0"]

            weights <- c(rep(3.5,7), 0.5, rep(1,2), rep(1.5,3), rep(2,4), rep(2.5,5), rep(3,6))

            # Compute data points for DevSeq_AL_pearson_dist based on model
            # First try to manually find rough parameters, then use nls to fine tune
            mcoeff <- nls(correlation ~ a * exp(div_times * c) + b * (1-(exp(div_times * c))), 
                start = list(a = 0.01, b = 0.5, c = -0.01), data = corrdata, control = list(maxiter = 500), 
                weights=weights)
            coeff <- as.data.frame(summary(mcoeff)["coefficients"])

            model_expr_dist <- data.frame(y = do.call(rbind, lapply(x_DS_grid, nl_model, 
                a = coeff["a",1], b = coeff["b",1], c = coeff["c",1])))

            model_coord <- data.frame(x = x_DS_grid, model_expr_dist)

            # Get slope values
            slopes = diff(model_coord$y)/diff(model_coord$x)
            slopes_avg <- mean(slopes)

            nlm_coord <- data.frame(model_coord, organ=rep(comp_organ, nrow(model_coord)), 
                nlm_slope=rep(slopes_avg, nrow(model_coord)))

            return(nlm_coord)

        }


        # Get mean nlm regression slopes for goslim category
        go_nlm_slopes <- data.frame(do.call(rbind, lapply(ortho_organ_lst, getNLEstimates)))
        go_nlm_mean_slp <- data.frame(nlm_slope=unique(go_nlm_slopes$nlm_slope), organ=unique(go_nlm_slopes$organ))


        # Get mean nlm regression slopes for all goslim control sets
        rtc_nlm_slopes <- data.frame(do.call(rbind, lapply(root_divc, getNLEstimates)))
        hcc_nlm_slopes <- data.frame(do.call(rbind, lapply(hypocotyl_divc, getNLEstimates)))
        lfc_nlm_slopes <- data.frame(do.call(rbind, lapply(leaf_divc, getNLEstimates)))
        avc_nlm_slopes <- data.frame(do.call(rbind, lapply(veg_apex_divc, getNLEstimates)))
        aic_nlm_slopes <- data.frame(do.call(rbind, lapply(inf_apex_divc, getNLEstimates)))
        flc_nlm_slopes <- data.frame(do.call(rbind, lapply(flower_divc, getNLEstimates)))
        stc_nlm_slopes <- data.frame(do.call(rbind, lapply(stamen_divc, getNLEstimates)))
        clc_nlm_slopes <- data.frame(do.call(rbind, lapply(carpel_divc, getNLEstimates)))


        # Set up list containing all organ goslim control sets
        control_nlm_slp_lst <- list(rtc_nlm_slopes, hcc_nlm_slopes, lfc_nlm_slopes, 
            avc_nlm_slopes, aic_nlm_slopes, flc_nlm_slopes, stc_nlm_slopes, clc_nlm_slopes)

        # Get loess mean and CI for control sets
        getCNLMStats <- function(cset) {

            cslope <- data.frame(nlm_slope_control=unique(cset$nlm_slope))

            gset <- data.frame(cslope, organ=unique(cset$organ, nrow(cslope)))

            return(gset)
        }

        control_nlm_slopes <- data.frame(do.call(rbind, lapply(control_nlm_slp_lst, getCNLMStats)))



        # Merge GOslim and control organ slope value tables 
        nlm_slope_df <- merge(go_nlm_mean_slp, control_nlm_slopes, by="organ", sort=FALSE)

        nlm_slope_df$goslim_term <- rep(unique(x_df$goslim), nrow(nlm_slope_df))




        #----------------- Prepare data and define color palette for corrplot -----------------

        # Create "plots" folder in /out_dir/output/plots
        if (!dir.exists(file.path(out_dir, "output", "plots"))) 
            dir.create(file.path(out_dir, "output", "plots"), recursive = TRUE)

        # Show message
        message("Starting analysis and generate plots...")


        # Combine all goslim organ list elements to one data frame
        ortho_organ_df <- do.call("rbind", ortho_organ_lst)
        ortho_organ_df$group <- rep(unique(x_df$goslim))


        # Compute mean corr distances if number of control sets is greater than 1
        if (ntreat > 1) {

            getMeanCorrSD <- function(lsel) {

                comb_cor <- rowMeans(as.data.frame(cbind(sapply(lsel, `[[`, "correlation"))))
                comb_sd <- rowMeans(as.data.frame(cbind(sapply(lsel, `[[`, "error"))))
                clade_name <- lsel[[1]]$clade
                organ_name <- lsel[[1]]$comp_organ
                div_times_c <- lsel[[1]]$div_times
                dataset_name <- lsel[[1]]$dataset
                comb_data <- data.frame(clade = clade_name, comp_organ = organ_name, 
                    div_times = div_times_c, correlation = comb_cor, error = comb_sd, 
                    dataset = dataset_name, group = rep("control"))

                return(comb_data)
            }

            root_divc_avg <- getMeanCorrSD(root_divc)
            hypocotyl_divc_avg <- getMeanCorrSD(hypocotyl_divc)
            leaf_divc_avg <- getMeanCorrSD(leaf_divc)
            veg_apex_divc_avg <- getMeanCorrSD(veg_apex_divc)
            inf_apex_divc_avg <- getMeanCorrSD(inf_apex_divc)
            flower_divc_avg <- getMeanCorrSD(flower_divc)
            stamen_divc_avg <- getMeanCorrSD(stamen_divc)
            carpel_divc_avg <- getMeanCorrSD(carpel_divc)

            control_organ_df <- rbind(root_divc_avg, hypocotyl_divc_avg, leaf_divc_avg, 
                veg_apex_divc_avg, inf_apex_divc_avg, flower_divc_avg, stamen_divc_avg, 
                carpel_divc_avg)

        } else if (ntreat == 1) {
            # Add data processing here
        }


        # Combine ortho and control corr data into final table for plotting
        ortho_control_dist <- rbind(ortho_organ_df, control_organ_df)


        # Change organ names for facet strip
        formOrganNames <- function(dist_df) {

            dist_df$comp_organ <- dist_df$comp_organ %<>% 
            gsub("Apex_veg", "Apex veg", .) %>% 
            gsub("Apex_inf", "Apex inf", .)

            dist_df$comp_organ <- factor(dist_df$comp_organ, 
                levels=c("Root","Hypocotyl","Leaf","Apex veg","Apex inf","Flower","Stamen","Carpel"))

            return(dist_df)
        }

        ortho_control_dist <- formOrganNames(ortho_control_dist)


        # Get organ slopes for averaged control sets
        control_organ_lst <- list(root_divc_avg, hypocotyl_divc_avg, leaf_divc_avg, 
                veg_apex_divc_avg, inf_apex_divc_avg, flower_divc_avg, stamen_divc_avg, 
                carpel_divc_avg)

        control_avg_nlm_slopes <- data.frame(do.call(rbind, lapply(control_organ_lst, getNLEstimates)))


        # Add group label to goslim and control slope tables
        go_nlm_slopes$group <- rep(unique(x_df$goslim))
        control_avg_nlm_slopes$group <- rep("control")

        # Combine ortho and control corr data into final table for plotting
        ortho_control_slopes <- rbind(go_nlm_slopes, control_avg_nlm_slopes)
        colnames(ortho_control_slopes)[which(names(ortho_control_slopes)=="x")] <- "div_times"
        colnames(ortho_control_slopes)[which(names(ortho_control_slopes)=="y")] <- "correlation"
        colnames(ortho_control_slopes)[which(names(ortho_control_slopes)=="organ")] <- "comp_organ"
        ortho_control_slopes <- formOrganNames(ortho_control_slopes)


        # Define specific notation
        set_scientific <- function(l) {
            # turn in to character string in scientific notation
            l <- format(l, scientific = TRUE)
            # quote the part before the exponent to keep all the digits
            l <- gsub("^(.*)e", "'\\1'e", l)
            # turn the 'e+' into plotmath format
            l <- gsub("e", "%*%10^", l)
            # return this as an expression
            parse(text=l)
        }



        # Plot pea distances and slopes of goslim and control data
        plotGOSLIM.pea.NLM <- function(data, data2) {

            fname <- sprintf('%s.jpg', paste(unique(x_df$goslim), "nlm_regression_slopes", sep="_"))

            # Define goslim colors for selected categories
            gocat <- as.character(unique(x_df$goslim))

            # Set plot title
            tname <- paste(unique(x_df$goslim))
            tltname <- paste(toupper(substr(tname, 1, 1)), substr(tname, 2, nchar(tname)), sep="")
            tltname <- paste0(tltname, " (n = ", ntreat, ")")

            # Get number of genes in GOslim category and number of control groups
            if (ntreat >= 900) {
                xpos <- 80
            } else xpos <- 73

            rt_data <- data2[data2$comp_organ == "Root",]
            y1pos <- (max(rt_data$correlation)*1.05)/3.44
            y2pos <- (max(rt_data$correlation)*1.05)/7
            hc_data <- data2[data2$comp_organ == "Hypocotyl",]
            y3pos <- (max(hc_data$correlation)*1.05)/7
            lf_data <- data2[data2$comp_organ == "Leaf",]
            y4pos <- (max(lf_data$correlation)*1.12)/7

            if (gocat == "response to chemical") {

                colscale <- c("#adadad", "#1e9ac7")
                p.value <- c(paste("italic('P =')~", set_scientific(0.002)))
                y3pos <- y3pos-0.001
                y4pos <- y4pos-0.002

            } else if (gocat == "embryo development") {

                colscale <- c("#adadad", "#cb0000")
                p.value <- c(paste("italic('P =')~", set_scientific(0.0004)))

            } else if (gocat == "nucleobase-containing compound metabolic process") {

                colscale <- c("#adadad", "#ee7500")
                p.value <- c(paste("italic('P =')~", set_scientific(0.00004)))
                y3pos <- y3pos-0.0015
                y4pos <- y4pos-0.0004

            } else if (gocat == "DNA binding") {

                colscale <- c("#adadad", "#08ac39")
                p.value <- c(paste("italic('P =')~", set_scientific(0.005)))
                y3pos <- y3pos-0.0025
                y4pos <- y4pos-0.0055

            } else if (gocat == "cellular component organization") {

                colscale <- c("#adadad", "#835bba")
                p.value <- c(paste("italic('P =')~", set_scientific(0.001)))
                y3pos <- y3pos-0.0025
                y4pos <- y4pos-0.0045

            } else {
                colscale <- c("#adadad", "black")
                p.value <- c("")
            }

            corg <- c("Root", "Hypocotyl", "Leaf", "Apex veg", "Apex inf", "Flower", "Stamen", "Carpel")

            goslim_lb <- data.frame(x = 82, y = y1pos, label = c("GOterm","","","","","","",""), 
                comp_organ = corg)
            control_lb <- data.frame(x = 82, y = y2pos, label = c("Control","","","","","","",""), 
                comp_organ = corg)
            c_text <- data.frame(x = xpos-53, y = y3pos, label = c("",paste("control sets:", length(root_divc)),"","","","","",""), 
                comp_organ = corg)

            p_text <- data.frame(x = 48, y = y4pos, label = c("","",p.value,"","","","",""), 
                comp_organ = corg)

            corgcat <- factor("Root", levels = c("Root", "Hypocotyl", "Leaf", 
                    "Apex veg", "Apex inf", "Flower", "Stamen", "Carpel"))

            go_line <- data.frame(x = 47, xend = 77, y = y1pos, yend = y1pos, 
                comp_organ = corgcat)
            cont_line <- data.frame(x = 47, xend = 77, y = y2pos, yend = y2pos, 
                comp_organ = corgcat)

            go_circ <- data.frame(x = 62, y = y1pos, comp_organ = corgcat)
            cont_circ <- data.frame(x = 62, y = y2pos, comp_organ = corgcat)

            data$group <- factor(data$group, c("control", paste(unique(x_df$goslim))))
            data2$group <- factor(data2$group, c("control", paste(unique(x_df$goslim))))

            p <- ggplot(data=data, color = group, aes(x=div_times, y=correlation)) + 
            geom_point(data=data2, alpha = 0.5, aes(stroke = 0.5, size = 1.5, color = group, shape = group, fill = group)) + 
            geom_line(size = 2.5, data = data, aes(x = div_times, y = correlation, group = group, color = group)) + 
            scale_y_continuous(expand = c(0.1, 0), breaks = pretty_breaks()) + 
            scale_x_continuous(expand = c(0.075, 0), breaks=c(0, 50, 100, 150)) + 
            scale_shape_manual(values = c(21,21)) + 
            scale_color_manual(values = colscale) + 
            scale_fill_manual(values = colscale) + 
            guides(shape = guide_legend(override.aes = list(stroke = 7.75)))

            q <- p + theme_classic() + xlab("Divergence time (Myr)") + ylab("Pearson distance") + 
            labs(title = tltname) + 
            geom_text(data = goslim_lb, mapping = aes(x = x, y = y, label = label), size=8, hjust = 0) + 
            geom_text(data = control_lb, mapping = aes(x = x, y = y, label = label), size=8, hjust = 0) + 
            geom_text(data = c_text, mapping = aes(x = x, y = y, label = label), size=8, hjust = 0) + 
            geom_text(data = p_text, mapping = aes(x = x, y = y, label = label), 
                parse=TRUE, size=8, hjust = 0) + 
            geom_segment(data = go_line, mapping = aes(x = x, xend = xend, y = y, yend = yend), 
                colour = colscale[2], show.legend = FALSE, size = 2.5) + 
            geom_point(data = go_circ, mapping = aes(x = x, y = y), size = 5, shape = 16, color = colscale[2]) + 
            geom_segment(data = cont_line, mapping = aes(x = x, xend = xend, y = y, yend = yend), 
                colour = colscale[1], show.legend = FALSE, size = 2.5) + 
            geom_point(data = cont_circ, mapping = aes(x = x, y = y), size = 5, shape = 16, color = colscale[1]) + 
            theme(text=element_text(size = 16), 
                strip.text = element_text(size = 23.75), 
                strip.text.x = element_text(margin = margin(0.46, 0, 0.46, 0, "cm")), 
                strip.background = element_rect(colour = 'black', fill = NA, size = 2.5), 
                axis.ticks.length = unit(0.29, "cm"), 
                axis.ticks = element_line(colour = "black", size = 1.25), 
                axis.line = element_line(colour = 'black', size = 1.25), 
                plot.margin = unit(c(1, 0.25, 3.0, 0),"cm"), 
                axis.title.y = element_text(size=27.35, margin = margin(t = 0, r = 12.75, b = 0, l = 12.0), 
                    colour="black", face = "bold"), 
                axis.title.x = element_text(size=27.35, margin = margin(t = 6.5, r = 0, b = 5.75, l = 0), 
                    colour="black", face = "bold"), 
                axis.text.x = element_text(size=22.5, margin = margin(t = 2.5, b = 8), colour="grey20"), 
                axis.text.y = element_text(size=22.5, angle=0, margin = margin(l = 2.5, r = 1.5), colour="grey20"), 
                plot.title = element_text(size=27.35, colour=colscale[2], margin = margin(t = 0, b = 15), face = "plain"), 
                panel.spacing = unit(0.2, "cm"), 
                panel.grid.major = element_blank(),
                panel.grid.minor.x = element_blank(), 
                panel.grid.minor.y = element_blank(), 
                legend.position = "none") 

            q <- q + facet_wrap(~ comp_organ, nrow = 1, scales = "free")

            ggsave(file = file.path(out_dir, "output", "plots", fname), plot = q, 
                width = 28.5, height = 6.5, dpi = 300, units = c("in"), limitsize = FALSE) 
        }

        plotGOSLIM.pea.NLM(data = ortho_control_slopes, data2 = ortho_control_dist)



        return(nlm_slope_df)




    })))






  getGoslimStats_lst <- split(getGoslimStats, getGoslimStats$goslim_term)
  wilcox_stats <- do.call(rbind, lapply(getGoslimStats_lst, function(i) {

    goslim_slope <- unique(i$nlm_slope)
    control_slope <- i$nlm_slope_control
    p_value <- wilcox.test(goslim_slope, control_slope)$p.value
    teststat <- data.frame(goslim_term=unique(i$goslim_term), p_value=p_value)
    return(teststat)

  }))


  wilcox_stats$p_value_FDR <- p.adjust(wilcox_stats$p_value, method = "fdr")


  # Perform permutation test as alternative to Wilcox rank sum test

  perm_stats <- do.call(rbind, lapply(getGoslimStats_lst, function(z) {

  	message("Performing permutation test...")
    treat <- unique(z$nlm_slope)
  	control <- z$nlm_slope_control
  	median_treat <- median(treat)
  	median_control <- median(control)
  	test_stat <- abs(median_treat-median_control)

  	set.seed(1357) # Set seed for reproducibility
  	all_slope_data <- data.frame(slope_value = c(treat, control), 
  		goslim_term = c(rep("GO", length(treat)),rep("Control", length(control))))
  	n_perm <- 100000
  	n_data <- nrow(all_slope_data)

  	perm_samples <- matrix(0, nrow = n_data, ncol = n_perm)
  	for (i in 1:n_perm) { 
  		perm_samples[,i] <- sample(all_slope_data$slope_value, size = n_data, replace = FALSE)
  	}

  	perm_test <- rep(0, n_perm)
  	for (i in 1:n_perm) { 
  		perm_test[i] <-  abs(
  			median(perm_samples[all_slope_data$goslim_term == "GO", i]) - 
  			median(perm_samples[all_slope_data$goslim_term == "Control", i])) 
  	}

  	perm_p <- sum(perm_test >= test_stat)/n_perm
  	permstat <- data.frame(goslim_term=unique(z$goslim_term), p_value=perm_p)
    return(permstat)

  }))


  perm_stats$p_value_FDR <- p.adjust(perm_stats$p_value, method = "fdr")


  
  # Make simple barplot of the results from Wilcoxon rank sum and permutation test
  term <- c(rownames(perm_stats), rownames(wilcox_stats))
  stat_summary_p <- c(round(perm_stats$p_value_FDR,4), round(wilcox_stats$p_value_FDR,4))
  testclass <- rep(c(" Permutation test  ", " Wilcoxon rank-sum test  "), each=nrow(perm_stats))
  statdata <- data.frame(term=term, data=stat_summary_p, test=testclass)


  plotTestStats <- function(pdata) {

    fname <- sprintf('%s.jpg', paste(aspect, "permutation_wilcox_p_values", sep="_"))

    if (aspect == "biological_process") {

        plt_w <- 25.5
        mb <- 11
        legpos <- c(0.5,0.95)
        btmm <- 0.25

    } else if (aspect == "molecular_function") {

        plt_w <- 12
        mb <- 17
        legpos <- "none"
        btmm <- 5.35
    }
    
    p <- ggplot(pdata, aes(factor(term), data, fill = test)) + 
    geom_bar(stat="identity", position = "dodge") + 
    scale_fill_brewer(palette = "Set1") + 
    scale_y_continuous(limits = c(0, 0.95)) + 
    xlab("") + 
    ylab("p value (FDR adjusted)") +  
    ggtitle(paste("GO", aspect, sep=" ")) +
    geom_text(aes(label=data), position=position_dodge(width=0.9), vjust=0.5, hjust=-0.25, angle=90, size=8.5) + 
    theme(axis.title = element_text(size=26),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=26), 
        axis.text.y = element_text(size=26, margin=margin(0,2,0,7)),
        legend.title = element_text(size=0),
        legend.text = element_text(size=26),
        legend.position = legpos,
        legend.direction = "horizontal",
        legend.key.size = unit(0.7, "cm"),
        axis.ticks.length = unit(0.2, "cm"),
        axis.ticks = element_line(colour = "black", size = 1),
        plot.title = element_text(size=28, margin=margin(20,20,mb,20)),
        plot.margin=unit(c(0.25,0.75,btmm,0.5), "cm"))

    ggsave(file = file.path(out_dir, "output", "plots", fname), plot = p, 
                width = plt_w, height = 23.25, dpi = 300, units = c("in"), limitsize = FALSE)
    }

    plotTestStats(pdata = statdata)




  # Write goslim and control slope data and test statistics to file
  # Show message
  message("Writing data tables...")

  # Create "data" folder in /out_dir/output
  if (!dir.exists(file.path(out_dir, "output", "data"))) 
    dir.create(file.path(out_dir, "output", "data"), recursive = TRUE)

  goslim_out_list <- list(getGoslimStats = getGoslimStats, wilcox_stats = wilcox_stats, 
  	perm_stats = perm_stats)

  for(i in names(goslim_out_list)){
    write.table(goslim_out_list[[i]], file=file.path(out_dir, "output", "data", paste0(i, "_", aspect, ".txt")), 
        sep="\t", col.names=TRUE, row.names=FALSE, dec=".", quote = FALSE)
  }





}
