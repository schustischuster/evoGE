# Extract GOslim terms with aspect of interest (biological process, molecular function) 
# from TAIR list downloaded on 19th April 2021
# GOSLIM table contains the GOslim terms for 28553 A.thaliana genes with "AT" identifier
# Create GOSLIM lists of 7003 angiosperm orthologous genes

library(dplyr)
library(MatchIt)

getGOSLIM <- function(aspect = c("biological_process", "molecular_function"), sample_size) {

    # Show error message if no/unknown GO aspect is chosen
	if ((missing(aspect)) || (!is.element(aspect, c("biological_process", "molecular_function"))))

		stop("Please choose one of the available aspects: 
			'biological_process', 'molecular_function'",
			call. = TRUE
			)

	# Show error message if no sample_size for GO term size is chosen
	if ((missing(sample_size)) || (sample_size < 1))

		stop("Please choose one of the available aspects",
			call. = TRUE
			)


	# Set file path for input files
	GOSLIM = file.path(in_dir, "Functional_groups", "ATH_GO_GOSLIM.txt")
	GOCAT = file.path(in_dir, "Functional_groups", "TAIR_GO_slim_categories.txt")
	orthoTPM = file.path(in_dir, "Expression_data", "AT_core_inter_tpm_mat_deseq_sample_names.csv")

	GOSLIM <- read.table(GOSLIM, sep="\t", dec=".", quote = "", header=FALSE, skip=4, fill = TRUE, stringsAsFactors=FALSE)
	GOCAT <- read.table(GOCAT, sep="\t", dec=".", header=TRUE, skip=7, fill = TRUE, stringsAsFactors=FALSE)
	orthoTPM <- read.table(orthoTPM, sep=";", dec=".", header=TRUE, stringsAsFactors=FALSE)


    # return_list <- list("orthoTPM" = orthoTPM, "GOSLIM" = GOSLIM, "GOCAT" = GOCAT, "aspect" = aspect, "sample_size" = sample_size)
    # return(return_list)
    # }
    # return_objects <- getGOSLIM(aspect = "biological_process", sample_size = 150)
    # list2env(return_objects, envir = .GlobalEnv)

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

    coreOrthologs <- sub("\\:.*", "", orthoTPM[,1]) # get AT orthologs
    coreOrthologs <- coreOrthologs[!grepl("ERCC", coreOrthologs)] # Rm spike-ins from ortholog list
    coreOrthologs <- as.data.frame(coreOrthologs)

    # Get list of orthologous genes that are associated with a GOSLIM term
    slim_ortho_ls <- lapply(slim_genes_uls, function(x){dplyr::filter(x, (V1 %in% coreOrthologs[,1]))})

    # Create stats table
    stats_all_genes <- lapply(slim_genes_uls, function(x){nrow(x)})
    stats_ortho_genes <- lapply(slim_ortho_ls, function(x){nrow(x)})

    stats_all_genes_df <- data.frame(all_genes = matrix(unlist(stats_all_genes), byrow=TRUE), stringsAsFactors=FALSE)
    stats_ortho_genes_df <- data.frame(ortho_genes = matrix(unlist(stats_ortho_genes), byrow=TRUE), stringsAsFactors=FALSE)
    goslim_northo_stats <- cbind(stats_all_genes_df, stats_ortho_genes_df)
    goslim_northo_stats$goslim_term <- names(slim_ortho_ls)
    goslim_northo_stats$expected <- (7003/28553)*goslim_northo_stats$all_genes
    goslim_northo_stats$fold_enrichment <- goslim_northo_stats$ortho_genes/goslim_northo_stats$expected
    goslim_northo_stats <- goslim_northo_stats[c("goslim_term", "all_genes", "ortho_genes", "expected", "fold_enrichment")]


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


    # Prepare angiosperm ortholog data
    orthoExpr <- data.frame(gene_id=sub("\\:.*", "", orthoTPM[,1]),orthoTPM[,2:ncol(orthoTPM)])
    orthoExpr[,2:ncol(orthoExpr)] <- log2(orthoExpr[,2:ncol(orthoExpr)] + 1)
    orthoExpr <- orthoExpr[!grepl("ERCC", orthoExpr$gene_id),]

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
            "Stamen", "Carpel", "Pollen"), each=7)
        DevSeq_spec_names <- rep(c("_AT", "_AL", "_CR", "_ES", "_TH", "_MT", "_BD"), times=9)
        repl_names <- paste0(DevSeq_col_names, DevSeq_spec_names)

        colnames(x_avg)[2:ncol(x_avg)] <- repl_names

        x_avg <- x_avg %>% select (-c(Pollen_AT, Pollen_AL, Pollen_CR, Pollen_ES, Pollen_TH, 
            Pollen_MT, Pollen_BD))


        # Compute average expression for each species
        calculateSpAvg <- function(xdf) {

            df <- xdf

            # Re-oder columns to allow easier splitting of species
            avg_names <- colnames(df[2:ncol(df)])
            spec_names <- gsub(".*_","", avg_names)
            avg_names <- paste(spec_names, avg_names, sep="_")
            avg_names <- gsub('.{3}$', '', avg_names)
            colnames(df)[2:ncol(df)] <- avg_names
            df <- df[,order(colnames(df))]
            df <- select(df, "gene_id", everything())

            # Split data frame by sample replicates into a list
            # then get rowMeans for each subset and bind averaged data to gene_id column

            averaged_spec <- do.call(cbind, lapply(split.default(df[2:ncol(df)], 
                rep(seq_along(df), 
                each = 8, 
                length.out=ncol(df)-1)
                ), rowMeans)
              )

            names_averaged_spec <- unique(sub("\\_.*", "", colnames(df)[2:ncol(df)]))
            avg_names <- paste("avg", names_averaged_spec, sep="_")
            colnames(averaged_spec) <- avg_names

            averaged_spec <- cbind(xdf[1], averaged_spec)
        
            return(averaged_spec)
        }

        spec_avg <- calculateSpAvg(x_avg)

        orthoExDf <- merge(spec_avg, x_avg)




    # Match control genes to each ortholog GOslim class
    getControls <- lapply(slim_ortho_ls, function(x){

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

        # Create background gene set
        background <- c()
        match_res <- matchit(sign ~ avg_AT + avg_AL + avg_CR + avg_ES + avg_TH + avg_MT + avg_BD, 
            comb_exdf, method="genetic", distance="mahalanobis", pop.size=1, replace=FALSE, 
            ratio=1)
        background <- c(background, match_res$match.matrix[,1]) 
        out <- comb_exdf[background,]



    })










}
