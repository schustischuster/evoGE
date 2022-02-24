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


    # return_list <- list("ldf" = ldf, "orthoTPM" = orthoTPM, "GOSLIM" = GOSLIM, "GOCAT" = GOCAT, "sample_size" = sample_size, "res" = res)
    # return(return_list)
    # }
    # return_objects <- plotGroupEx(sample_size = 412)
    # list2env(return_objects, envir = .GlobalEnv)

    # Show message
    message("Starting analysis...")

    list2env(ldf, envir = .GlobalEnv)



    #------------------- Calculate average core ortholog expression and CV --------------------


    # Prepare angiosperm ortholog data
    orthoExpr <- data.frame(gene_id=sub("\\:.*", "", orthoEst[,1]),orthoEst[,2:ncol(orthoEst)])
    
    if (estimate == "TPM") { 
        orthoExpr[,2:ncol(orthoExpr)] <- log2(orthoExpr[,2:ncol(orthoExpr)] + 1)
    }

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
                fname <- sprintf('%s.jpg', paste("mean_coeff_var", estimate, sample_size, sep="_"))
                x_title <- "Orthologous genes"
                y_title <- "Mean coefficient of variation"
                y_lim <- c(-0.025, max(data$value))

            } else if (plotclass == "express") {

                data <- data.frame(value = data$base_averaged, class = data$sign, group = data$group)
                fname <- sprintf('%s.jpg', paste("coeff_var_expr", estimate, sample_size, sep="_"))
                x_title <- "Orthologous genes"

                if (estimate == "VST") {
                    y_title <- "Expression (VST counts)"
                } else if (estimate == "TPM") {
                    y_title <- "log2(TPM+1)"
                }
                
                y_lim <- c(min(data$value)-(0.15*min(data$value)), max(data$value))
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
            scale_y_continuous(expand = c(0.028, 0), limits = y_lim) + 
            guides(shape = guide_legend(override.aes = list(stroke = 7.75)))

            q <- p + theme_classic() + xlab(x_title) + ylab(y_title) + 
            geom_text(data = dat_text, mapping = aes(x = x, y = y, label = label), 
                size = 6.5, vjust = 0, colour = "grey25") + 
            theme(text=element_text(size = 16), 
                strip.text = element_text(size = 19.5), 
                strip.text.x = element_text(margin = margin(0.37, 0, 0.37, 0, "cm")), 
                strip.background = element_rect(colour = 'black', fill = NA, size = 2.5), 
                axis.ticks.length = unit(0.24, "cm"), 
                axis.ticks = element_line(colour = "black", size = 1.2), 
                axis.line = element_line(colour = 'black', size = 1.2), 
                plot.margin = unit(c(0.5, 15.5, 1.25, 0), "cm"), 
                axis.title.y = element_text(size = 22.0, margin = margin(t = 0, r = 6.4, b = 0, l = 10), 
                    colour="black", face = "bold"), 
                axis.title.x = element_text(size = 22.0, margin = margin(t = 4.0, r = 0, b = 7.0, l = 0), 
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
    

    if (aspect == "biological_process") {

        selCAT <- dplyr::filter(GOCAT, grepl("biological process", ONTOLOGY.ASPECT))

    } else if (aspect == "molecular_function") {

        selCAT <- dplyr::filter(GOCAT, grepl("molecular function", ONTOLOGY.ASPECT))
    }

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


    # Clean up GOslim term list: Remove specific categories
    if (aspect == "biological_process") {
        
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
    
    } else if (aspect == "molecular_function") {

        # Correct GOslim term for class "RNA binding" (it contains another term called "translation factor activity, RNA binding")
        slim_genes_uls[["RNA binding"]]$V9 <- rep("RNA binding", nrow(slim_genes_uls[["RNA binding"]]))

    }


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


    # Reshape data for ggplot2
    cv_stats_rs <- data.frame(
        GO_term = rep(cv_stats$GO_term,2), CV_cat = rep(c("stable genes", "variable genes"), each = nrow(cv_stats)), 
        n_genes = c(cv_stats$stable_genes, cv_stats$variable_genes), 
        chisq_test = rep(cv_stats$chisq_test,2), fisher_test = rep(cv_stats$fisher_test,2), 
        chisq_FDR = rep(cv_stats$chisq_FDR,2), fisher_FDR = rep(cv_stats$fisher_FDR,2))
    cv_stats_rs <- cv_stats_rs[order(cv_stats_rs$GO_term),]


    # Define specific notation
        set_scientific <- function(l) {
            # turn in to character string in scientific notation
            l <- formatC(l, format = "e", digits = 0)
            # quote the part before the exponent to keep all the digits
            l <- gsub("^(.*)e", "'\\1'e", l)
            # turn the 'e+' into plotmath format
            l <- gsub("e", "%*%10^", l)
            # return this as an expression
            parse(text=l)
        }



    # Plot mean expression of orthologs for GOslim categories that show significant different 
    # proportion of stable/variable genes
    # Table of core orthologs with average expression value = x_avg
    # List of all GOslim categories above size threshold = slim_ortho_ls
    # Table with all GOslim categories that show significant different proportion of table/
    # variable genes = cv_stats_rs

}






