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

        stable_G <- merge(df, stable_genes, by = "gene_id")
        variable_G <- merge(df, variable_genes, by = "gene_id")

        # Create contingency table
        tab <- matrix(c(nrow(stable_G), nrow(variable_G), nrow(stable_genes)-nrow(stable_G), 
            nrow(variable_genes)-nrow(variable_G)), ncol=2, byrow=TRUE)

        colnames(tab) <- c('stable', 'variable')
        rownames(tab) <- c('GO', 'other')
        tab <- as.table(tab)

        chisq <- chisq.test(tab, correct=T)$p.value

        fishtest <- fisher.test(tab)$p.value

        output <- data.frame(GO_term = go_term,
            stable_genes = nrow(stable_G), 
            variable_genes = nrow(variable_G),
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








