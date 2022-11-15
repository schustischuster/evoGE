# Compute coefficient of variation (CV) for 7003 angiosperm core orthologous genes
# Classifiy genes into evolutionarily stable and variable genes based on average CV across organs
# Match stable and variable genes based on their average gene expression level, using caliper matching
# Similar procedure as described in Berthelot et al., Nat Ecol Evol (2018)
# Check proportion of stable/variable genes for GOslim terms of interest  
# GOslim categories retrieved from TAIR version 20211101


getCV <- function(aspect = c("biological_process", "molecular_function"), estimate = c("VST", "TPM"), 
    sample_size, ...) {

    # Show error message if no/unknown GO aspect is chosen
    if ((missing(aspect)) || (!is.element(aspect, c("biological_process", "molecular_function"))))

        stop("Please choose one of the available aspects: 
            'biological_process', 'molecular_function'",
            call. = TRUE
            )

    # Show error message if no expression estimate is chosen
    if ((missing(estimate)) || (!is.element(estimate, c("VST", "TPM"))))

        stop("Please choose one of the available expression estimates: 
            'VST', 'TPM'",
            call. = TRUE
            )

    # Show error message if no sample_size for GO term size is chosen
    if ((missing(sample_size)) || (sample_size < 1))

        stop("Please choose GOslim sample size cutoff",
            call. = TRUE
            )

    # Set file path for input files
	GOSLIM = file.path(in_dir, "ATH_GO_GOSLIM.txt")
	GOCAT = file.path(in_dir, "TAIR_GO_slim_categories.txt")

    if (estimate == "VST" ) {

        orthoEst = file.path(in_dir, "AT_core_inter_count_mat_vsd_sample_names.csv")

    } else { 

        orthoEst = file.path(in_dir, "AT_core_inter_TPM_mat_deseq_sample_names.csv")
    }

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


    # return_list <- list("ldf" = ldf, "orthoEst" = orthoEst, "GOSLIM" = GOSLIM, "GOCAT" = GOCAT, "aspect" = aspect, "estimate" = estimate, "sample_size" = sample_size, "res" = res)
    # return(return_list)
    # }
    # return_objects <- getCV(aspect = "biological_process", estimate = "VST", sample_size = 412)
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


    # Plot results of CV analysis
    plotCV <- function(df, cat) {

        fname <- sprintf('%s.pdf', paste(deparse(substitute(df)), cat, sep="_"))

        # Capitalize first string of each GO term category name
        CapStr <- function(y) {
            c <- strsplit(y, " ")[[1]]
            paste(toupper(substring(c, 1,1)), substring(c, 2),
                sep="", collapse=" ")
        }

        df$GO_term <- paste(sapply(gsub("([A-Za-z]+).*", "\\1", df$GO_term), CapStr), 
            sub(".*? ", "", df$GO_term))

        # Correct some formating errors
        df$GO_term <- gsub("Transport transport", "Transport", df$GO_term)
        df$GO_term <- gsub("Reproduction reproduction", "Reproduction", df$GO_term)

        df$GO_term <- gsub("Post development", "Post-embryonic development", df$GO_term)
        df$GO_term <- gsub("Embryo development", "Embryo and post-embryonic development", df$GO_term)

        # Subset and reorder data
        if (cat == "stable") {

            gene_list <- c("Post-embryonic development", "Reproduction", 
                "Cellular component organization", "Transport", 
                "Nucleobase compound metabolic process", "Embryo and post-embryonic development", 
                "Protein metabolic process")

            df <- df[df$GO_term %in% gene_list,]
            df$GO_term <- factor(df$GO_term, levels = gene_list)

            plt_title <- "GO categories with larger fraction of stable genes   "
            plot_mar = unit(c(0.9, -0.25, 0.47, -0.325), "cm")
            legend_pos <- "none"
            x_pand <- c(0.0435, 0)
            y_dat_pos <- 567
            y_lim <- c(0, 755)
            y_breaks <- c(0, 100, 200, 300, 400, 500)

        } else if (cat == "variable") {

            gene_list <- c("Response to abiotic stimulus", "Response to chemical", 
                "Response to external stimulus", "Response to biotic stimulus", 
                "Response to endogenous stimulus", "Response to stress")

            df <- df[df$GO_term %in% gene_list,]
            df$GO_term <- factor(df$GO_term, levels = gene_list)

            plt_title <- "GO categories with larger fraction of variable genes"
            plot_mar = unit(c(0.9, 1.9, 0.47, -0.11), "cm")
            legend_pos <- "top"
            x_pand <- c(0.048, 0)
            y_dat_pos <- 675
            y_lim <- c(0, 900)
            y_breaks <- c(0, 200, 400, 600)

        }

        # Create df for FDR p-value mapping
        FDR_df <- data.frame(x = df[seq(1, nrow(df), 2), "GO_term"], 
            y = rep(y_dat_pos, nrow(df)/2), 
            p_val = c(paste("italic('P =')~", set_scientific(df[seq(1, nrow(df), 2), "chisq_FDR"]))), 
            CV_cat = df[seq(1, nrow(df), 2), "CV_cat"]
            )

        df$CV_cat <- factor(df$CV_cat, levels = c('variable genes', 'stable genes'))

        p <- ggplot(df, aes(x = reorder(GO_term, n_genes), y = n_genes, fill = CV_cat)) + 
        geom_bar(position = "dodge", stat = "identity", width = 0.725, size = 0.7, col = "black") + 
        coord_flip() + 
        scale_y_continuous(breaks = y_breaks, limits = y_lim, expand = c(0, 0)) + 
        scale_x_discrete(expand = x_pand) + 
        labs(y = "Number of Genes                    " , x = NULL) + 
        scale_fill_manual(values = c("#1e77d1", "#fbc11a"), limits = c("stable genes", "variable genes"), labels=c(" Stable  ", " Variable  ")) + 
        guides(fill = guide_legend(keywidth = 0.305, keyheight = 0.305, default.unit="inch")) + 
        ggtitle(plt_title) +  
        geom_text(data = FDR_df, aes(x = x, y = y, label = p_val), size = 7, parse=TRUE, hjust = 0, vjust = 0.375) + 
        theme(panel.background = element_blank(), 
            axis.ticks.length = unit(0.26, "cm"), 
            axis.ticks = element_line(colour = "black", size = 0.91), 
            axis.line = element_line(colour = 'black', size = 0.91), 
            plot.margin = plot_mar, 
            plot.title = element_text(size = 21.5, margin = margin(t = 0, r = 0, b = 2, l = 0), hjust = 2.0), 
            legend.box.margin = margin(2.95, 82, -17.05, 0), 
            legend.text = element_text(size = 20.5), 
            legend.title = element_blank(), 
            legend.direction = "horizontal", 
            legend.position = legend_pos, 
            legend.key = element_rect(colour = "transparent", fill = "white"),
            axis.title.y = element_text(size = 20.5, margin = margin(t = 0, r = 7.0, b = 0, l = 10), 
                colour="black", face = "plain"), 
            axis.title.x = element_text(size = 20.5, margin = margin(t = 2.5, r = 0, b = 8.15, l = 0), 
                colour="black", face = "plain"), 
            axis.text.x = element_text(size = 19.35, margin = margin(t = 4.5, b = 8), colour = "black"), 
            axis.text.y = element_text(size = 20.28, angle=0, margin = margin(l = 10, r = 2.5), colour = "black")
        )

        ggsave(file = file.path(out_dir, "output", "plots", fname), plot = p, 
               width = 11.5, height = 6.5, dpi = 300, units = c("in"), limitsize = FALSE)

    }

    plotCV(cv_stats_rs, cat = "stable")
    plotCV(cv_stats_rs, cat = "variable")

}






