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

    # Show message
    message("Starting analysis...")



    #--------- Extract GOslim data for DevSeq core orthologs and check GO enrichment ----------


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


    # Clean up GOslim term list: Remove specific categories
    if (aspect == "biological_process") {
        
        # Remove "other" categories because not informative
        # Remove "Response to light stimulus" because plants were grown at constant light
        slim_genes_uls <- slim_genes_uls %>% purrr::list_modify("other metabolic processes" = NULL)
        slim_genes_uls <- slim_genes_uls %>% purrr::list_modify("other cellular processes" = NULL)
        slim_genes_uls <- slim_genes_uls %>% purrr::list_modify("other biological processes" = NULL)
        slim_genes_uls <- slim_genes_uls %>% purrr::list_modify("response to light stimulus" = NULL)
        # Remove some remaining functional ontologies from list
        slim_genes_uls <- lapply(slim_genes_uls, function(x){dplyr::filter(x, !grepl("F", V8))})
        # Delete "Cell Growth" GOslim term from "Growth" category
        slim_genes_uls[["growth"]] <- dplyr::filter(slim_genes_uls[["growth"]], !grepl("cell growth", V9))
    
    } else if (aspect == "molecular_function") {

    }


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



    #------------ Combine DevSeq core ortholog expression tables with GOslim data -------------


    # Prepare angiosperm ortholog data
    orthoExpr <- data.frame(gene_id=sub("\\:.*", "", orthoTPM[,1]),orthoTPM[,2:ncol(orthoTPM)])
    orthoExpr[,2:ncol(orthoExpr)] <- log2(orthoExpr[,2:ncol(orthoExpr)] + 1)
    orthoExpr <- orthoExpr[!grepl("ERCC", orthoExpr$gene_id),]


    # Negate dplyr %in%
    `%!in%` = Negate(`%in%`)


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


        # Compute average expression and sd for each organ
        calculateSpAvg <- function(df) {

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

        if (ntreat <= 1000 & ntreat >= 501) {
            cratio <- 2
        } else if (ntreat <= 500 & ntreat >= 351) {
            cratio <- 3
        } else if (ntreat <= 350 & ntreat >= 251) {
            cratio <- 4
        } else if (ntreat <= 250 & ntreat >= 176) {
            cratio <- 5
        } else if (ntreat <= 175) {
            cratio <- 6
        } else cratio <- 1

        # Create background gene set
        background <- c()
        match_res <- matchit(sign ~ base_averaged, comb_exdf, 
            method="nearest", distance="mahalanobis", replace=FALSE, ratio=cratio)
        background <- c(background, match_res$match.matrix) 
        control_out <- comb_exdf[background,]

        control_out <- control_out %>% select (-c(avg_Root, avg_Hypocotyl, avg_Leaf, avg_veg, 
            avg_inf, avg_Flower, avg_Stamen, avg_Carpel, base_averaged))

        goslim_out <- merge(x_df, orthoExDf)

        goslim_out <- goslim_out %>% select (-c(avg_Root, avg_Hypocotyl, avg_Leaf, avg_veg, 
            avg_inf, avg_Flower, avg_Stamen, avg_Carpel, base_averaged))




        #---------------- Get gene expression divergence rates for ATH/AL vs species X -----------------


        # Show message
        message("Calculating expression distances...")


        ortho_control <- merge(control_out, orthoExpr)
        ortho_go <- merge(goslim_out, orthoExpr)

        ortho_control <- ortho_control[-2:-59]
        ortho_go <- ortho_go[-2:-59]


        # Use pearson correlation, inter-organ normalization and TPM for ms

        getDSOrganCor <- function(df, organ) {

            df_cor <- sqrt(1/2*(1 - cor(df, method="pearson")))
            df_cor <- df_cor[4:nrow(df_cor),]

            sp1 <- mean(df_cor[1:3,1:3])
            sp2 <- mean(df_cor[4:6,1:6])
            sp3 <- mean(df_cor[7:9,1:9])
            sp4 <- mean(df_cor[10:12,1:12])
            sp5 <- mean(df_cor[13:15,1:15])
            sp6 <- mean(df_cor[16:18,1:18])

            # Get errors
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

            sp1_repl <- avgRepl(df_cor[1:3,1:3]) # AT-AL
            sp2_repl <- avgRepl(df_cor[4:6,1:6]) # AT-CR AL-CR
            sp3_repl <- avgRepl(df_cor[7:9,1:9]) # AT-ES AL-ES CR-ES
            sp4_repl <- avgRepl(df_cor[10:12,1:12]) # AT-TH AL-TH CR-TH ES-TH
            sp5_repl <- avgRepl(df_cor[13:15,1:15]) # AT-MT AL-MT CR-MT ES-MT TH-MT
            sp6_repl <- avgRepl(df_cor[16:18,1:18]) # AT-BD AL-BD CR-BD ES-BD TH-BD MT-BD

            getError <- function(cor_data) {
                std <- sd(cor_data, na.rm=TRUE)
                num <- length(cor_data)
                error <- std/sqrt(num)
                return(error)
            } # Use this function to replace sd if reqired

            df_cor_error <- data.frame(error = c(as.numeric(c(sd(sp1_repl))),
                    as.numeric(c(sd(sp2_repl))), as.numeric(c(sd(sp3_repl))), 
                    as.numeric(c(sd(sp4_repl))), as.numeric(c(sd(sp5_repl))), 
                    as.numeric(c(sd(sp6_repl)))))

            df_cor_avg <- data.frame(correlation = c(sp1, sp2, sp3, sp4, sp5, sp6))
            div_tag <- data.frame(clade = c("T1", "T2", "T3", "T4", "T5", "T6"))
            organ_id <- data.frame(comp_organ = rep(organ, 6))
            div_times <- data.frame(div_times = c(7.1, 9.4, 25.6, 46, 106, 160))
            dataset <- data.frame(dataset = rep("Angiosperms ", 6))
            df_cor_avg <- cbind(div_tag, organ_id, div_times, df_cor_avg, df_cor_error, dataset)

            return(df_cor_avg)

        }

        root_divc <- getDSOrganCor(df=ortho_control[,2:22], organ="Root")
        hypocotyl_divc <- getDSOrganCor(df=ortho_control[,23:43], organ="Hypocotyl")
        leaf_divc <- getDSOrganCor(df=ortho_control[,44:64], organ="Leaf")
        veg_apex_divc <- getDSOrganCor(df=ortho_control[,65:85], organ="Apex veg")
        inf_apex_divc <- getDSOrganCor(df=ortho_control[,86:106], organ="Apex inf")
        flower_divc <- getDSOrganCor(df=ortho_control[,107:127], organ="Flower")
        stamen_divc <- getDSOrganCor(df=ortho_control[,128:148], organ="Stamen")
        carpel_divc <- getDSOrganCor(df=ortho_control[,149:169], organ="Carpel")

        control_div_rates <- rbind(root_divc, hypocotyl_divc, leaf_divc, veg_apex_divc, inf_apex_divc, 
            flower_divc, stamen_divc, carpel_divc)
    

        root_divg <- getDSOrganCor(df=ortho_go[,2:22], organ="Root")
        hypocotyl_divg <- getDSOrganCor(df=ortho_go[,23:43], organ="Hypocotyl")
        leaf_divg <- getDSOrganCor(df=ortho_go[,44:64], organ="Leaf")
        veg_apex_divg <- getDSOrganCor(df=ortho_go[,65:85], organ="Apex veg")
        inf_apex_divg <- getDSOrganCor(df=ortho_go[,86:106], organ="Apex inf")
        flower_divg <- getDSOrganCor(df=ortho_go[,107:127], organ="Flower")
        stamen_divg <- getDSOrganCor(df=ortho_go[,128:148], organ="Stamen")
        carpel_divg <- getDSOrganCor(df=ortho_go[,149:169], organ="Carpel")

        ortho_div_rates <- rbind(root_divg, hypocotyl_divg, leaf_divg, veg_apex_divg, inf_apex_divg, 
            flower_divg, stamen_divg, carpel_divg)

        # Set up lists containing metric pearson expression distances
        control_organ_lst <- list(control_div_rates[1:6,], control_div_rates[7:12,], 
            control_div_rates[13:18,], control_div_rates[19:24,], control_div_rates[25:30,], 
            control_div_rates[31:36,], control_div_rates[37:42,], control_div_rates[43:48,])

        ortho_organ_lst <- list(ortho_div_rates[1:6,], ortho_div_rates[7:12,], 
            ortho_div_rates[13:18,], ortho_div_rates[19:24,], ortho_div_rates[25:30,], 
            ortho_div_rates[31:36,], ortho_div_rates[37:42,], ortho_div_rates[43:48,])


        getLOESS.Slopes <- function(organ_data) {

            comp_organ <- unique(organ_data$comp_organ)

            # Use quadratic polynomes for all organs
            temp <- loess.smooth(organ_data$div_times, organ_data$correlation, span = 1, 
                degree = 2, family="gaussian", evaluation = 200)

            # Get slope values
            slopes = diff(temp$y)/diff(temp$x)
            slopes_avg <- mean(slopes)
            slopes_avg <- as.numeric(as.data.frame(slopes_avg))

            return(slopes_avg)
        }

        control_loess_slopes <- as.data.frame(do.call(rbind, lapply(control_organ_lst, getLOESS.Slopes)))
        go_loess_slopes <- as.data.frame(do.call(rbind, lapply(ortho_organ_lst, getLOESS.Slopes)))

        loess_out <- data.frame(

            goslim = rep(unique(x_df$goslim), nrow(go_loess_slopes)*2),
            group = c(rep("control", nrow(control_loess_slopes)), rep("functional", nrow(go_loess_slopes))),
            slopes = c(control_loess_slopes[,1], go_loess_slopes[,1]),
            p_value = rep(wilcox.test(as.numeric(control_loess_slopes[,1]), 
                as.numeric(go_loess_slopes[,1]))$p.value, nrow(go_loess_slopes)*2)
        )

        return(loess_out)



    })))










}
