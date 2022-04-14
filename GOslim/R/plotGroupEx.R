# Plot mean expression of orthologs from all GOslim categories that show a significant
# different proportion of evolutionarily stable and variable genes
# Make plots for both matched and all genes for each GOslim category
# Also show distribution of mean coefficient of variation (CV) before and after matching
# GOslim categories retrieved from TAIR version 20211101


plotGroupEx <- function(sample_size, ...) {


    # Show error message if no sample_size for GO term size is chosen
    if ((missing(sample_size)) || (sample_size < 1))

        stop("Please choose GOslim sample size cutoff",
            call. = TRUE
            )

    # Set file path for input files
    GOSLIM = file.path(in_dir, "ATH_GO_GOSLIM.txt")
    GOCAT = file.path(in_dir, "TAIR_GO_slim_categories.txt")


    orthoEst = file.path(in_dir, "AT_core_inter_count_mat_vsd_sample_names.csv")
    orthoTPM = file.path(in_dir, "AT_core_inter_TPM_mat_deseq_sample_names.csv")


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
    orthoEst <- read.table(orthoEst, sep=";", dec=".", header=TRUE, stringsAsFactors=FALSE)
    orthoTPM <- read.table(orthoTPM, sep=";", dec=".", header=TRUE, stringsAsFactors=FALSE)


    # return_list <- list("ldf" = ldf, "orthoEst" = orthoEst, "orthoTPM" = orthoTPM, "GOSLIM" = GOSLIM, "GOCAT" = GOCAT, "sample_size" = sample_size, "res" = res)
    # return(return_list)
    # }
    # return_objects <- plotGroupEx(sample_size = 412)
    # list2env(return_objects, envir = .GlobalEnv)

    # Show message
    message("Starting analysis...")

    list2env(ldf, envir = .GlobalEnv)



    #------------------- Calculate average core ortholog expression and CV --------------------


    # Negate dplyr %in%
    `%!in%` = Negate(`%in%`)


    # Prepare angiosperm ortholog data
    procExprTable <- function(z) {

        zExpr <- data.frame(gene_id=sub("\\:.*", "", z[,1]), z[,2:ncol(z)])
        zExpr <- zExpr[!grepl("ERCC", zExpr$gene_id),]

        # Remove pollen samples
        zExpr <- zExpr %>% select (-c(
        A.thaliana_flowers_mature_pollen_28d_.2.:B.distachyon_flowers_mature_pollen_32d_.1.))

        return(zExpr)
    }

    orthoExpr <- procExprTable(orthoEst)
    orthoTPM <- procExprTable(orthoTPM)


    orthoTPM[,2:ncol(orthoTPM)] <- log2(orthoTPM[,2:ncol(orthoTPM)] + 1)


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
    x_tavg <- calculateAvgExpr(orthoTPM)

    DevSeq_col_names <- rep(c("Root", "Hypocotyl", "Leaf", "veg_apex", "inf_apex", "Flower", 
    	"Stamen", "Carpel"), each=7)
    DevSeq_spec_names <- rep(c("_AT", "_AL", "_CR", "_ES", "_TH", "_MT", "_BD"), times=8)
    repl_names <- paste0(DevSeq_col_names, DevSeq_spec_names)

    colnames(x_avg)[2:ncol(x_avg)] <- repl_names
    colnames(x_tavg)[2:ncol(x_tavg)] <- repl_names


    # Compute average expression and CV for each organ
    calculateAvgCV <- function(df) {

            averaged_spec <- do.call(cbind, lapply(split.default(df[2:ncol(df)], 
                rep(seq_along(df), 
                each = 7, 
                length.out=ncol(df)-1)
                ), rowMeans)
              )

            base_averaged <- rowMeans(df[2:ncol(df)])
   

            RowSD <- function(x) {
                sqrt(rowSums((x - rowMeans(x))^2)/(dim(x)[2] - 1))
            }

            averaged_sd <- do.call(cbind, lapply(split.default(df[2:ncol(df)], 
                rep(seq_along(df), 
                each = 7, 
                length.out=ncol(df)-1)
                ), RowSD)
              )

            organCV <- averaged_sd/averaged_spec

            CV_averaged <- rowMeans(organCV)

            names_organ <- unique(sub("\\_.*", "", colnames(df)[2:ncol(df)]))
            avg_names <- paste("CV", names_organ, sep="_")
            colnames(organCV) <- avg_names

            averaged_spec <- cbind(df[1], organCV, CV_averaged, base_averaged)
        
            return(averaged_spec)
    }

    spec_CV <- calculateAvgCV(x_avg)


    # Order coeff var table by base expression value
    spec_CV <- spec_CV[order(spec_CV$CV_averaged),]

    spec_CV$sign <- c(rep(1, 3501), rep(0, 3502))


    # Perform nearest distance matching with caliper option
    # Check balance statistics to find optimal caliper
    message("Matching genes...")

    calpr <- 0.5

        matchSample <- function(x) {

            success <- FALSE
            while (!success) {

                # Create background gene set
                match_res <- matchit(sign ~ base_averaged, x, method = "nearest", caliper = calpr, 
                    std.caliper = TRUE, distance = "logit", replace = FALSE, m.order = "data", 
                    ratio = 1)
                match_res_m <- match_res$match.matrix

                # Extract standard mean difference from matchIt summary data
                comp <- as.data.frame(summary(match_res, standardize = TRUE)["sum.matched"])
                stmdif <- abs(comp[1,3])
                varR <- abs(comp[1,4])

                calpr <- calpr-0.01

                # check for success
                success <- ((stmdif <= 0.01) && (varR >= 0.99))
            }

            return(match_res_m)
        }

    match_res_df <- matchSample(spec_CV)

    match_res_sep <- data.frame(treat = rownames(match_res_df), control = match_res_df)
    match_res_sep <- match_res_sep[!is.na(match_res_sep$control),]

    stable_genes <- subset(x_avg, rownames(x_avg) %in% match_res_sep$treat)
    variable_genes <- subset(x_avg, rownames(x_avg) %in% match_res_sep$control)

    stable_genes <- merge(stable_genes, spec_CV)
    variable_genes <- merge(variable_genes, spec_CV)


    # Generate table of all matched genes (4896) containing average CV and expression values
    # spec_CV table = table containing average CV and expression values for all genes (7003)
    x_avg_match <- rbind(stable_genes, variable_genes)
    spec_CV_match <- x_avg_match %>% select (-c(Root_AT:Carpel_BD))
    spec_CV_all <- spec_CV

    spec_CV_all$group <- rep("All", nrow(spec_CV_all))
    spec_CV_match$group <- rep("Matched", nrow(spec_CV_match))
    spec_CV_df <- rbind(spec_CV_all, spec_CV_match)

    spec_CV_df$sign[spec_CV_df$sign == 1] <- "stable"
    spec_CV_df$sign[spec_CV_df$sign == 0] <- "variable"



    # Plot pea distances and slopes of goslim and control data
    plotCVExpr <- function(data, plotclass) {

            if (plotclass == "mCV") {

                data <- data.frame(value = data$CV_averaged, class = data$sign, group = data$group)
                fname <- sprintf('%s.jpg', paste("mean_coeff_var_VST", sample_size, sep="_"))
                x_title <- "Orthologous genes"
                y_title <- "Mean coefficient of variation"
                y_lim <- c(-0.0225, max(data$value))

            } else if (plotclass == "express") {

                data <- data.frame(value = data$base_averaged, class = data$sign, group = data$group)
                fname <- sprintf('%s.jpg', paste("coeff_var_expr_VST", sample_size, sep="_"))
                x_title <- "Orthologous genes"
                y_title <- "Expression (VST counts)"
                y_lim <- c(min(data$value)-(0.14*min(data$value)), max(data$value))
            }

            n_all <- sum(data['class'] == "stable" & data['group'] == "All")
            n_sub <- sum(data['class'] == "stable" & data['group'] == "Matched")

            dat_text <- data.frame(
                label = c(n_all, n_sub, n_all, n_sub),
                group   = c( "All", "Matched"),
                x     = c(1, 1, 2, 2),
                y     = c(rep(y_lim[1], 4))
                )

            p <- ggplot(data = data, color = class, aes(x=class, y=value)) + 
            geom_boxplot(colour = "black", size = 1.2, fatten = 2.5, notch = TRUE, 
                outlier.shape = 21, outlier.size = 3.5, fill = rep(c("orange", "blueviolet"), 2), 
                outlier.color = "grey25", outlier.alpha = 0.28, outlier.fill = "grey25") + 
            scale_y_continuous(expand = c(0.035, 0), limits = y_lim) + 
            guides(shape = guide_legend(override.aes = list(stroke = 7.75)))

            q <- p + theme_classic() + xlab(x_title) + ylab(y_title) + 
            geom_text(data = dat_text, mapping = aes(x = x, y = y, label = label), 
                size = 6.5, vjust = 0, colour = "grey25") + 
            theme(text=element_text(size = 16), 
                strip.text = element_text(size = 19.85), 
                strip.text.x = element_text(margin = margin(0.375, 0, 0.375, 0, "cm")), 
                strip.background = element_rect(colour = 'black', fill = NA, size = 2.5), 
                axis.ticks.length = unit(0.24, "cm"), 
                axis.ticks = element_line(colour = "black", size = 1.2), 
                axis.line = element_line(colour = 'black', size = 1.2), 
                plot.margin = unit(c(0.5, 15.95, 1.25, 0), "cm"), 
                axis.title.y = element_text(size = 22.25, margin = margin(t = 0, r = 6.4, b = 0, l = 10), 
                    colour="black", face = "bold"), 
                axis.title.x = element_text(size = 22.25, margin = margin(t = 4.0, r = 0, b = 7.0, l = 0), 
                    colour="black", face = "bold"), 
                axis.text.x = element_text(size = 19.0, margin = margin(t = 3.5, b = 8), colour="grey20"), 
                axis.text.y = element_text(size = 18.8, angle = 0, margin = margin(l = 2.5, r = 2.5), colour="grey20"), 
                panel.spacing = unit(0.5, "cm"), 
                panel.grid.major = element_blank(),
                panel.grid.minor.x = element_blank(), 
                panel.grid.minor.y = element_blank(), 
                legend.position = "none") 

            q <- q + facet_wrap(~ group, ncol = 2)

            ggsave(file = file.path(out_dir, "output", "plots", fname), plot = q, 
                width = 11.5, height = 6.5, dpi = 300, units = c("in"), limitsize = FALSE) 
        }

        plotCVExpr(data = spec_CV_df, plotclass = "mCV")
        plotCVExpr(data = spec_CV_df, plotclass = "express")




    #--- Extract GOslim data for DevSeq core orthologs and classify genes (stable/variable) ---

    message("Processing GOslim terms...")
    

    selCAT <- dplyr::filter(GOCAT, grepl("biological process", ONTOLOGY.ASPECT))

    slim_names <- selCAT$SLIM_NAME
    slim_name_ls <- setNames(as.list(c(slim_names)), c(slim_names)) # generate named list with all GOSLIM terms

    # Extract all genes for each GOslim term of slim_name_ls list
    # V1 in GOSLIM table contains the gene ID, V9 contains the GOslim terms
    slim_genes_ls <- lapply(slim_name_ls, function(x){dplyr::filter(GOSLIM, grepl(x,V9))})

    # Get list with items reduced to unique gene id's
    slim_genes_uls <- lapply(slim_genes_ls, function(x){x[!duplicated(x[,1]),]})

    coreOrthologs <- sub("\\:.*", "", orthoEst[,1]) # get AT orthologs
    coreOrthologs <- coreOrthologs[!grepl("ERCC", coreOrthologs)] # Rm spike-ins from ortholog list
    coreOrthologs <- as.data.frame(coreOrthologs)


    # Clean up GOslim term list
    
    # Remove "other" categories because not informative
    # Remove "response to light stimulus" because plants were grown at constant light
    # Remove "cell cycle" because those genes have strong, cell-type specific expression that makes it impossible to match controls
    slim_genes_uls <- slim_genes_uls %>% purrr::list_modify("other metabolic processes" = NULL)
    slim_genes_uls <- slim_genes_uls %>% purrr::list_modify("other cellular processes" = NULL)
    slim_genes_uls <- slim_genes_uls %>% purrr::list_modify("other biological processes" = NULL)
    slim_genes_uls <- slim_genes_uls %>% purrr::list_modify("response to light stimulus" = NULL)
    slim_genes_uls <- slim_genes_uls %>% purrr::list_modify("cell cycle" = NULL)
    # Remove some remaining functional ontologies from list
    slim_genes_uls <- lapply(slim_genes_uls, function(x){dplyr::filter(x, !grepl("F", V8))})
    # Delete "Cell Growth" GOslim term from "Growth" category
    slim_genes_uls[["growth"]] <- dplyr::filter(slim_genes_uls[["growth"]], !grepl("cell growth", V9))


    # Get list of orthologous genes that are associated with a GOSLIM term
    slim_ortho_ls <- lapply(slim_genes_uls, function(x){dplyr::filter(x, (V1 %in% coreOrthologs[,1]))})

    # Remove all ortholog GOslim lists wth fewer entries than defines sample_size
    slim_ortho_ls <- Filter(function(dt) nrow(dt) >= sample_size, slim_ortho_ls)


    # Does the proportion of stable/variable genes in GO category differ from the proportion
    # of stable/variable genes across all orthologs => Chi-Square Test/Fisher's Exact Test
    checkStability <- function(df) {

        go_term <- unique(df$V9)

        names(df)[1] <- "gene_id" 

        stable_GO <- merge(df, stable_genes, by = "gene_id")
        variable_GO <- merge(df, variable_genes, by = "gene_id")

        # Create contingency table
        tab <- matrix(c(nrow(stable_GO), nrow(variable_GO), nrow(stable_genes)-nrow(stable_GO), 
            nrow(variable_genes)-nrow(variable_GO)), ncol=2, byrow=TRUE)

        colnames(tab) <- c('stable', 'variable')
        rownames(tab) <- c('GO', 'other')
        tab <- as.table(tab)

        chisq <- chisq.test(tab, correct=T)$p.value

        fishtest <- fisher.test(tab)$p.value

        output <- data.frame(GO_term = go_term,
            stable_genes = nrow(stable_GO), 
            variable_genes = nrow(variable_GO),
            chisq_test = chisq, 
            fisher_test = fishtest
            )

        return(output)
    }

    cv_stats <- do.call(rbind, lapply(slim_ortho_ls, checkStability))
    rownames(cv_stats) <- NULL


    # FDR correction of Chi-Square p-value
    cv_stats$chisq_FDR <- p.adjust(cv_stats$chisq_test, method = "fdr")
    cv_stats$fisher_FDR <- p.adjust(cv_stats$fisher_test, method = "fdr")

    # Remove all GOslim categories with a p-adjust value greater than 0.001
    cv_stats <- dplyr::filter(cv_stats, cv_stats$chisq_FDR <= 0.001, cv_stats$fisher_FDR <= 0.001)


    # Select GOslim categories that show significant chisq_test p-value
    slim_ortho_all_ls <- slim_ortho_ls[c(as.character(cv_stats$GO_term))] # all core orthologos

    getGOCat <- function(r) {
        gene_id <- r$V1
        category <- r$V9
        df <- cbind(gene_id = gene_id, category = category)
    }

    slim_ortho_all_df <- data.frame(do.call(rbind, lapply(slim_ortho_all_ls, getGOCat)))

    x_tavg$base_averaged <- rowMeans(x_tavg[2:ncol(x_tavg)])
    x_tavg <- x_tavg %>% select(gene_id, base_averaged)

    # Replace mean VST values with mean TPM expression values in all/matched groups
    spec_CV_all <- merge(subset(spec_CV_all, select=-c(base_averaged)), x_tavg)
    spec_CV_match <- merge(subset(spec_CV_match, select=-c(base_averaged)), x_tavg)

    slim_ortho_all_expr <- merge(slim_ortho_all_df, spec_CV_all, by = "gene_id") %>% 
                           select(-c(CV_Root:CV_averaged, sign)) %>% 
                           arrange(category)

    slim_ortho_match_expr <- merge(slim_ortho_all_df, spec_CV_match, by = "gene_id") %>% 
                             select(-c(CV_Root:CV_averaged, sign)) %>% 
                             arrange(category)

    slim_ortho_expr <- rbind(slim_ortho_all_expr, slim_ortho_match_expr)



    # Plot mean expression of GOslim groups across samples
    plotGOEx <- function(data) {

        fname <- sprintf('%s.jpg', paste("mean_expr_across_samples", sample_size, sep="_"))

        # Capitalize first string of each GO term category name
        CapStr <- function(y) {
            c <- strsplit(y, " ")[[1]]
            paste(toupper(substring(c, 1,1)), substring(c, 2),
                sep="", collapse=" ")
        }

        data$category <- paste(sapply(gsub("([A-Za-z]+).*", "\\1", data$category), CapStr), 
            sub(".*? ", "", data$category))

        # Correct some formating errors
        data$category <- gsub("Transport transport", "Transport", data$category)
        data$category <- gsub("Reproduction reproduction", "Reproduction", data$category)

        data$category <- gsub("Post development", "Post-embryonic development", data$category)
        data$category <- gsub("Embryo development", "Embryo and post-embryonic development", data$category)

        y_title <- "log2(TPM+1)"
        y_lim <- c(min(data$base_averaged), max(data$base_averaged)*1.1)


        # Extract number of genes for each GO category
        # Do this for both all genes and matched genes
        get_all <- data[data$group == "All",]
        get_match <- data[data$group == "Matched",]
        get_cat <- unique(get_all$category)

        extr_all_num <- function(t) {
            num <- sum(get_all$category == t)
            out <- data.frame(category = t, num = num, group = rep("All"))
            return(out)
        }

        all_num <- data.frame(do.call(rbind, lapply(get_cat, extr_all_num)))

        extr_match_num <- function(u) {
            num <- sum(get_match$category == u)
            out <- data.frame(category = u, num = num, group = rep("Matched"))
            return(out)
        }
        
        match_num <- data.frame(do.call(rbind, lapply(get_cat, extr_match_num)))

        dat_text <- rbind(all_num, match_num)
        dat_text$y <- rep(10.2)

        # Label facet strips
        strip_names <- c(
                    `All` = paste0("All ","(", "n=", nrow(spec_CV_all), ")"),
                    `Matched` = paste0("Matched ","(", "n=", nrow(spec_CV_match), ")")
                    )

        p <- ggplot(data = data, color = group, aes(x=reorder(category, base_averaged, .fun='median'), y=base_averaged)) + 
        geom_boxplot(colour = "black", size = 1.2, fatten = 2.5, notch = TRUE, 
            outlier.shape = NA, fill = rep(c("orange", "blueviolet"), each = 13)) + 
        coord_flip(ylim = quantile(data$base_averaged, c(0.0025, 0.9975))) + 
        scale_x_discrete(expand = c(0.025, 0)) + 
        guides(shape = guide_legend(override.aes = list(stroke = 7.75)))

        q <- p + theme_classic() + xlab("") + ylab(y_title) + 
        # geom_text(data = dat_text, mapping = aes(x = category, y = y, label = num, group = group), 
            # size = 6.65, vjust = 0.5, hjust = 1, colour = "grey50") + 
        theme(text=element_text(size = 16), 
            strip.text = element_text(size = 19.5), 
            strip.text.x = element_text(margin = margin(0.37, 0, 0.37, 0, "cm")), 
            strip.background = element_rect(colour = 'black', fill = NA, size = 2.5), 
            axis.ticks.length = unit(0.24, "cm"), 
            axis.ticks = element_line(colour = "black", size = 1.2), 
            axis.line = element_line(colour = 'black', size = 1.2), 
            plot.margin = unit(c(0.5, 0.2, 1.25, -0.4), "cm"), 
            axis.title.y = element_text(size = 22.0, margin = margin(t = 0, r = 0, b = 0, l = 0), 
                colour="black", face = "bold"), 
            axis.title.x = element_text(size = 22.0, margin = margin(t = 4.0, r = 0, b = 7.0, l = 0), 
                colour="black", face = "bold"), 
            axis.text.x = element_text(size = 18.8, margin = margin(t = 3.5, b = 8), colour="grey20"), 
            axis.text.y = element_text(size = 19.0, angle = 0, margin = margin(l = 0, r = 2.5), colour="grey20"), 
            panel.spacing = unit(0.5, "cm"), 
            panel.grid.major = element_blank(),
            panel.grid.minor.x = element_blank(), 
            panel.grid.minor.y = element_blank(), 
            legend.position = "none") 

        q <- q + facet_wrap(~ group, ncol = 2, labeller = as_labeller(strip_names))

        ggsave(file = file.path(out_dir, "output", "plots", fname), plot = q, 
            width = 11.5, height = 6.5, dpi = 300, units = c("in"), limitsize = FALSE) 
    }

    plotGOEx(data = slim_ortho_expr)





}






