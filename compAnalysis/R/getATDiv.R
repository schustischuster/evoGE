# Prepare Brawand and DevSeq comparative expression data
# Thresholds: 0.5 TPM (since there are no ERCC spike-ins in Brawand data)
# Data input: Brawand and DevSeq TPM expression tables of all samples



#-------------------------------------- Read data tables ---------------------------------------


getATDiv <- function(expr_estimation = c("TPM", "counts"), coefficient = c("pearson", "spearman")) {
	

    # Show error message if no correlation or unknown correlation coefficient is chosen
    if ((missing(coefficient)) || (!is.element(coefficient, c("pearson", "spearman"))))
   
       stop(
       "Please choose one of the available correlation coefficients: 
	   'pearson', 'spearman'",
	   call. = TRUE
       )

    # Show error message if expression estimation or unknown expression estimation is chosen
    if ((missing(expr_estimation)) || (!is.element(expr_estimation, c("TPM", "counts"))))
   
       stop(
       "Please choose one of the available expression estimations: 
       'TPM', 'counts'",
       call. = TRUE
       )


    # Show startup message
    message("Reading data...")


    # Set file path for input files
    if (is.element("TPM", expr_estimation)) {
        
        genesExprDS = file.path(in_dir, "Expression_data", "AT_core_inter_tpm_mat_deseq_sample_names.csv")
        genesExprBr = file.path(in_dir, "Expression_data", "Brawand_inter_tpm_mat_deseq_sample_names_0_5_threshold.csv")

    } else if (is.element("counts", expr_estimation)) {
        
        genesExprDS = file.path(in_dir, "Expression_data", "AT_core_inter_count_mat_vsd_sample_names.csv")
        genesExprBr = file.path(in_dir, "Expression_data", "Brawand_inter_count_mat_vsd_sample_names_0_5_threshold.csv")

    }

    # Read original Brawand ortholog expression data (RPKM) from 2011 Nature publication
    genesExprBr2011 = file.path(in_dir, "Expression_data", "Brawand_Supplementary_Data1", "NormalizedRPKM_ConstitutiveAlignedExons_Amniote1to1Orthologues.txt")

    
    # Set colnames
    col_namesDS <- rep(c("Root", "Hypocotyl", "Leaf", "veg_apex", "inf_apex", "Flower", "Stamen", 
        "Carpel", "Pollen"), each=21)
    replicate_tag_samples <- rep(c(".1",".2",".3"), times=9)
    col_namesDS <- paste0(col_namesDS,replicate_tag_samples)
    spec_namesDS <- rep(c("_AT", "_AL", "_CR", "_ES", "_TH", "_MT", "_BD"), each=3)
    spec_namesDS <- rep(spec_namesDS, times=9)
    col_namesDS <- paste0(col_namesDS, spec_namesDS)
    col_namesDS <- c("gene_id", col_namesDS)


    # Read DevSeq table
	x_DS <- read.table(genesExprDS, sep=";", dec=".", header=TRUE, stringsAsFactors=FALSE)
    
    # set column names
    colnames(x_DS) <- col_namesDS


	# Read Brawand table and set colnames
	x_Br <- read.table(genesExprBr, sep=";", dec=".", header=TRUE, stringsAsFactors=FALSE)

    # Remove later on once expression table has gene_id column
    ID_repl <- as.data.frame(seq(1:nrow(x_Br)))
    colnames(ID_repl) <- "gene_id"
    x_Br <- cbind(ID_repl, x_Br)


    # Read original Brawand expression table from Nature 2011
    x_Br2011 <- read.table(genesExprBr2011, sep="\t", dec=".", header=TRUE, stringsAsFactors=FALSE)

    # Remove biological replicates that show log2 RPMK Pearson's r < 0.85
    # remove ortholog names, two hsa samples with low cor and platypus and chicken data 
    x_Br2011 <- x_Br2011 %>% select (-c(hsa:gga, hsa.br.M.4, hsa.br.M.5, oan.br.M.1:gga.ts.M.2))
    ID_repl <- as.data.frame(seq(1:nrow(x_Br2011)))
    colnames(ID_repl) <- "gene_id"
    x_Br2011 <- cbind(ID_repl, x_Br2011)	


    # Stop function here to allow specific analysis of a single data set
    # return_list <- list("x_DS" = x_DS, "x_Br" = x_Br, "expr_estimation" = expr_estimation, "coefficient" = coefficient, "x_Br2011" = x_Br2011)
    # return(return_list)
    # }
    # return_objects <- getATDiv(expr_estimation = "TPM", coefficient = "pearson") # read in Comparative expression data
    # list2env(return_objects, envir = .GlobalEnv)

    


#--------------------- Calculate correlation and prepare data for corrplot ---------------------


    # Create "plots" folder in /out_dir/output/plots
    if (!dir.exists(file.path(out_dir, "output", "plots"))) 
        dir.create(file.path(out_dir, "output", "plots"), recursive = TRUE)

    # Show message
    message("Starting analysis and generate plots...")



    x_Br[is.na(x_Br)] <- 0 # replaces NAs by 0
    x_DS[is.na(x_DS)] <- 0 # replaces NAs by 0
    # Remove ERCC spike-ins from data
    x_DS <- x_DS[!grepl("ERCC", x_DS$gene_id),]
    x_Br2011[is.na(x_Br2011)] <- 0 # replaces NAs by 0


    if (is.element("pearson", coefficient) && is.element("TPM", expr_estimation)) {
            
        x_Br[,2:ncol(x_Br)] <- log2(x_Br[,2:ncol(x_Br)] + 1)
        x_DS[,2:ncol(x_DS)] <- log2(x_DS[,2:ncol(x_DS)] + 1)
        x_Br2011[,2:ncol(x_Br2011)] <- log2(x_Br2011[,2:ncol(x_Br2011)] + 1)

    }




#----------------- Read taxa objects for DevSeq and Brawand with replicates  -------------------


    if (is.element("TPM", expr_estimation)) {

   
        # Construc taxa object
        x_Br_taxa_objects = TEconstruct(ExpValueFP = file.path(out_dir, 
            "output", "data", 'x_Br_taxobj_input.txt'), taxa = "all", subtaxa = 'all')

        x_Br2011_taxa_objects = TEconstruct(ExpValueFP = file.path(out_dir, 
            "output", "data", 'x_Br2011_taxobj_input.txt'), taxa = "all", subtaxa = 'all')

        x_DS_taxa_objects = TEconstruct(ExpValueFP = file.path(out_dir, 
            "output", "data", 'x_DS_taxobj_input.txt'), taxa = "all", subtaxa = 'all')



        Brawand_organ_list <- list("brain", "cerebellum", "heart", "kidney", "liver", "testis")

        Brawand2011_organ_list <- list("br", "cb", "ht", "kd", "lv", "ts")

        DevSeq_organ_list <- list("Root", "Hypocotyl", "Leaf", "vegApex", "infApex", "Flower", 
            "Stamen", "Carpel")



        # Apply extended OU model with dynamic expression optimum ("variable-µ method")
        getExtOUBr <- function(organ, taxa_obj) {

            sou_v_out <- expdist(taxa_obj, taxa = "all",
                subtaxa = organ,
                method = "sou_v")

            # sou_v_pi <- sou_v_out$pi ##### To retrieve pi #####
            sou_v_distance <- as.data.frame(sou_v_out$distance)

            getError <- function(cor_data) {
                std <- sd(cor_data, na.rm=TRUE)
                num <- length(cor_data)
                error <- std/sqrt(num)
                return(error)
            } # Use this function to replace sd if reqired

            if (organ == "ts") {

                sou_v_distance_div <- data.frame(correlation = c(sou_v_distance[3, 1], 
                    mean(as.numeric(c(sou_v_distance[4, c(1:3)]))), 
                    NA, # ppy (orangutan) data missing
                    mean(as.numeric(c(sou_v_distance[5, c(1:4)]))), 
                    mean(as.numeric(c(sou_v_distance[6, c(1:5)]))), 
                    mean(as.numeric(c(sou_v_distance[7, c(1:6)])))))

                sou_v_distance_error <- data.frame(error = c(NA, 
                    as.numeric(c(sd(sou_v_distance[4, c(1:3)]))), 
                    NA, # ppy (orangutan) data missing
                    as.numeric(c(sd(sou_v_distance[5, c(1:4)]))), 
                    as.numeric(c(sd(sou_v_distance[6, c(1:5)]))), 
                    as.numeric(c(sd(sou_v_distance[7, c(1:6)])))))

            } else if (organ == "testis") {

                sou_v_distance_div <- data.frame(correlation = c(sou_v_distance[2, 1], 
                    mean(as.numeric(c(sou_v_distance[4, c(1:3)]))), 
                    NA, # ppy (orangutan) data missing
                    mean(as.numeric(c(sou_v_distance[5, c(1:4)]))), 
                    mean(as.numeric(c(sou_v_distance[6, c(1:5)]))), 
                    mean(as.numeric(c(sou_v_distance[7, c(1:6)])))))

                sou_v_distance_error <- data.frame(error = c(NA, 
                    as.numeric(c(sd(sou_v_distance[4, c(1:3)]))), 
                    NA, # ppy (orangutan) data missing
                    as.numeric(c(sd(sou_v_distance[5, c(1:4)]))), 
                    as.numeric(c(sd(sou_v_distance[6, c(1:5)]))), 
                    as.numeric(c(sd(sou_v_distance[7, c(1:6)])))))

            } else {

                sou_v_distance_div <- data.frame(correlation = c(
                    mean(as.numeric(c(sou_v_distance[2:3, 1]))), 
                    mean(as.numeric(c(sou_v_distance[4, c(1:3)]))), 
                    mean(as.numeric(c(sou_v_distance[5, c(1:4)]))), 
                    mean(as.numeric(c(sou_v_distance[6, c(1:5)]))), 
                    mean(as.numeric(c(sou_v_distance[7, c(1:6)]))), 
                    mean(as.numeric(c(sou_v_distance[8, c(1:7)])))))

                sou_v_distance_error <- data.frame(error = c(
                    as.numeric(c(sd(sou_v_distance[2:3, 1]))),
                    as.numeric(c(sd(sou_v_distance[4, c(1:3)]))), 
                    as.numeric(c(sd(sou_v_distance[5, c(1:4)]))), 
                    as.numeric(c(sd(sou_v_distance[6, c(1:5)]))), 
                    as.numeric(c(sd(sou_v_distance[7, c(1:6)]))), 
                    as.numeric(c(sd(sou_v_distance[8, c(1:7)])))))
            }

            div_tag <- data.frame(clade = c("T1", "T2", "T3", "T4", "T5", "T6"))
            organ_id <- data.frame(comp_organ = rep(unique(sub(".*_", "", rownames(sou_v_distance))),6))
            sou_v_distance_div <- cbind(div_tag, organ_id, sou_v_distance_div, sou_v_distance_error)

            return(sou_v_distance_div)
        }


        Br_sou_v <- as.data.frame(do.call(rbind, lapply(Brawand_organ_list, getExtOUBr,
            taxa_obj = x_Br_taxa_objects)))

        Br2011_sou_v <- as.data.frame(do.call(rbind, lapply(Brawand2011_organ_list, getExtOUBr,
            taxa_obj = x_Br2011_taxa_objects)))


        getExtOUDS <- function(organ, taxa_obj) {

            sou_v_out <- expdist(taxa_obj, taxa = "all",
                subtaxa = organ,
                method = "sou_v")

            # sou_v_pi <- sou_v_out$pi ##### To retrieve pi #####
            sou_v_distance <- as.data.frame(sou_v_out$distance)

            getError <- function(cor_data) {
                std <- sd(cor_data, na.rm=TRUE)
                num <- length(cor_data)
                error <- std/sqrt(num)
                return(error)
            } # Use this function to replace sd if reqired


            if (organ == "Stamen") {

                sou_v_distance_div <- data.frame(correlation = c(sou_v_distance[2,1], 
                    mean(as.numeric(c(sou_v_distance[3, c(1:2)]))), 
                    mean(as.numeric(c(sou_v_distance[4, c(1:3)]))), 
                    mean(as.numeric(c(sou_v_distance[5, c(1:4)]))), 
                    mean(as.numeric(c(sou_v_distance[6, c(1:5)]))), 
                    mean(as.numeric(c(sou_v_distance[7, c(1:4,6)])))))

                sou_v_distance_error <- data.frame(error = c(NA, 
                    as.numeric(c(sd(sou_v_distance[3, c(1:2)]))), 
                    as.numeric(c(sd(sou_v_distance[4, c(1:3)]))), 
                    as.numeric(c(sd(sou_v_distance[5, c(1:4)]))), 
                    as.numeric(c(sd(sou_v_distance[6, c(1:5)]))), 
                    as.numeric(c(sd(sou_v_distance[7, c(1:4,6)])))))

            } else {

                sou_v_distance_div <- data.frame(correlation = c(sou_v_distance[2,1], 
                    mean(as.numeric(c(sou_v_distance[3, c(1:2)]))), 
                    mean(as.numeric(c(sou_v_distance[4, c(1:3)]))), 
                    mean(as.numeric(c(sou_v_distance[5, c(1:4)]))), 
                    mean(as.numeric(c(sou_v_distance[6, c(1:5)]))), 
                    mean(as.numeric(c(sou_v_distance[7, c(1:6)])))))

                sou_v_distance_error <- data.frame(error = c(NA, 
                    as.numeric(c(sd(sou_v_distance[3, c(1:2)]))), 
                    as.numeric(c(sd(sou_v_distance[4, c(1:3)]))), 
                    as.numeric(c(sd(sou_v_distance[5, c(1:4)]))), 
                    as.numeric(c(sd(sou_v_distance[6, c(1:5)]))), 
                    as.numeric(c(sd(sou_v_distance[7, c(1:6)])))))
            }

            div_tag <- data.frame(clade = c("T1", "T2", "T3", "T4", "T5", "T6"))
            organ_id <- data.frame(comp_organ = rep(unique(sub(".*_", "", rownames(sou_v_distance))),6))
            sou_v_distance_div <- cbind(div_tag, organ_id, sou_v_distance_div, sou_v_distance_error)

            return(sou_v_distance_div)
        }

        DS_sou_v <- as.data.frame(do.call(rbind, lapply(DevSeq_organ_list, getExtOUDS,
            taxa_obj = x_DS_taxa_objects)))

    }




#---------------- Get gene expression divergence rates for ATH/AL vs species X -----------------


    # Use pearson correlation, inter-organ normalization and TPM for ms

    getDSOrganCor <- function(df, organ, coefficient) {

        df_cor <- sqrt(1/2*(1 - cor(df, method=coefficient)))
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
    

    root_div <- getDSOrganCor(df=x_DS[,2:22], organ="Root", coefficient=coefficient)
    hypocotyl_div <- getDSOrganCor(df=x_DS[,23:43], organ="Hypocotyl", coefficient=coefficient)
    leaf_div <- getDSOrganCor(df=x_DS[,44:64], organ="Leaf", coefficient=coefficient)
    veg_apex_div <- getDSOrganCor(df=x_DS[,65:85], organ="Apex veg", coefficient=coefficient)
    inf_apex_div <- getDSOrganCor(df=x_DS[,86:106], organ="Apex inf", coefficient=coefficient)
    flower_div <- getDSOrganCor(df=x_DS[,107:127], organ="Flower", coefficient=coefficient)
    stamen_div <- getDSOrganCor(df=x_DS[,128:148], organ="Stamen", coefficient=coefficient)
    carpel_div <- getDSOrganCor(df=x_DS[,149:169], organ="Carpel", coefficient=coefficient)

    DevSeq_div_rates <- rbind(root_div, hypocotyl_div, leaf_div, veg_apex_div, inf_apex_div, 
        flower_div, stamen_div, carpel_div)



    if (is.element("TPM", expr_estimation)) {

        ds_dataset <- data.frame(dataset = rep("Angiosperms ", 48))
        div_times <- data.frame(div_times = rep(c(7.1, 9.4, 25.6, 46, 106, 160),8))

        DevSeq_sou_v_div_rates <- data.frame(cbind(DS_sou_v, ds_dataset, div_times), 
            stringsAsFactors=FALSE)

        DevSeq_sou_v_div_rates <- DevSeq_sou_v_div_rates %>% select(clade, comp_organ, div_times, everything())

    }




#---------------- Get gene expression divergence rates for Human vs species X ------------------


    # Use pearson correlation, inter-organ normalization and TPM for ms
    getError <- function(cor_data) {
                std <- sd(cor_data, na.rm=TRUE)
                num <- length(cor_data)
                error <- std/sqrt(num)
                return(error)
            } # Use this function to replace sd if reqired

    getBrBrainCor <- function(df, organ, coefficient) {

        df_cor <- sqrt(1/2*(1 - cor(df, method=coefficient)))
        df_cor <- df_cor[5:nrow(df_cor), ]

        sp1_repl <- c(mean(df_cor[1:3,1:4]), mean(df_cor[4:9,1:4])) # Hsa-Ppa Hsa-Ptr
        sp2_repl <- c(mean(df_cor[10:11,1:4]), mean(df_cor[10:11,5:7]), mean(df_cor[10:11,8:13])) # Hsa-Ggo Ppa-Ggo Ptr-Ggo
        sp3_repl <- c(mean(df_cor[12:13,1:4]), mean(df_cor[12:13,5:7]), mean(df_cor[12:13,8:13]), mean(df_cor[12:13,14:15])) # Hsa-Ppy Ppa-Ppy Ptr-Ppy Ggo-Ppy
        sp4_repl <- c(mean(df_cor[14:16,1:4]), mean(df_cor[14:16,5:7]), mean(df_cor[14:16,8:13]), mean(df_cor[14:16,14:15]), 
            mean(df_cor[14:16,16:17])) # Hsa-Mml Ppa-Mml Ptr-Mml Ggo-Mml Ppy-Mml
        sp5_repl <- c(mean(df_cor[17:19,1:4]), mean(df_cor[17:19,5:7]), mean(df_cor[17:19,8:13]), mean(df_cor[17:19,14:15]), 
            mean(df_cor[17:19,16:17]), mean(df_cor[17:19,18:20])) # Hsa-Mmu Ppa-Mmu Ptr-Mmu Ggo-Mmu Ppy-Mmu Mml-Mmu
        sp6_repl <- c(mean(df_cor[20:21,1:4]), mean(df_cor[20:21,5:7]), mean(df_cor[20:21,8:13]), mean(df_cor[20:21,14:15]), 
            mean(df_cor[20:21,16:17]), mean(df_cor[20:21,18:20]), mean(df_cor[20:21,21:23])) # Hsa-Mdo Ppa-Mdo Ptr-Mdo Ggo-Mdo Ppy-Mdo Mml-Mdo Mmu-Mdo

        # Get mean and SE
        df_cor_avg <- data.frame(correlation = c(mean(sp1_repl), mean(sp2_repl), mean(sp3_repl), 
            mean(sp4_repl), mean(sp5_repl), mean(sp6_repl)))

        df_cor_error <- data.frame(error = c(as.numeric(c(sd(sp1_repl))),
            as.numeric(c(sd(sp2_repl))), as.numeric(c(sd(sp3_repl))), 
            as.numeric(c(sd(sp4_repl))), as.numeric(c(sd(sp5_repl))), 
            as.numeric(c(sd(sp6_repl)))))

        div_tag <- data.frame(clade = c("T1", "T2", "T3", "T4", "T5", "T6"))
        organ_id <- data.frame(comp_organ = rep(organ, 6))
        df_cor_avg <- cbind(div_tag, organ_id, df_cor_avg, df_cor_error)

        return(df_cor_avg)

    }

    brain_div <- getBrBrainCor(df=x_Br[,2:26], organ="Brain", coefficient=coefficient)

    brain_div_11 <- getBrBrainCor(df=x_Br2011[,c(2:5,33:35,18:23,45:46,56:57,65:67,78:80,95:96)], 
        organ="Brain", coefficient=coefficient)


    getBrCerebCor <- function(df, organ, coefficient) {

        df_cor <- sqrt(1/2*(1 - cor(df, method=coefficient)))
        df_cor <- df_cor[3:nrow(df_cor), ]

        sp1_repl <- c(mean(df_cor[1:2,1:2]), mean(df_cor[3:4,1:2])) # Hsa-Ppa Hsa-Ptr
        sp2_repl <- c(mean(df_cor[5:6,1:2]), mean(df_cor[5:6,3:4]), mean(df_cor[5:6,5:6])) # Hsa-Ggo Ppa-Ggo Ptr-Ggo
        sp3_repl <- c(mean(df_cor[7,1:2]), mean(df_cor[7,3:4]), mean(df_cor[7,5:6]), mean(df_cor[7,7:8])) # Hsa-Ppy Ppa-Ppy Ptr-Ppy Ggo-Ppy
        sp4_repl <- c(mean(df_cor[8:9,1:2]), mean(df_cor[8:9,3:4]), mean(df_cor[8:9,5:6]), mean(df_cor[8:9,7:8]), 
            mean(df_cor[8:9,9])) # Hsa-Mml Ppa-Mml Ptr-Mml Ggo-Mml Ppy-Mml
        sp5_repl <- c(mean(df_cor[10:12,1:2]), mean(df_cor[10:12,3:4]), mean(df_cor[10:12,5:6]), mean(df_cor[10:12,7:8]), 
            mean(df_cor[10:12,9]), mean(df_cor[10:12,10:11])) # Hsa-Mmu Ppa-Mmu Ptr-Mmu Ggo-Mmu Ppy-Mmu Mml-Mmu
        sp6_repl <- c(mean(df_cor[13:14,1:2]), mean(df_cor[13:14,3:4]), mean(df_cor[13:14,5:6]), mean(df_cor[13:14,7:8]), 
            mean(df_cor[13:14,9]), mean(df_cor[13:14,10:11]), mean(df_cor[13:14,12:14])) # Hsa-Mdo Ppa-Mdo Ptr-Mdo Ggo-Mdo Ppy-Mdo Mml-Mdo Mmu-Mdo

        # Get mean and SE
        df_cor_avg <- data.frame(correlation = c(mean(sp1_repl), mean(sp2_repl), mean(sp3_repl), 
            mean(sp4_repl), mean(sp5_repl), mean(sp6_repl)))

        df_cor_error <- data.frame(error = c(as.numeric(c(sd(sp1_repl))),
            as.numeric(c(sd(sp2_repl))), as.numeric(c(sd(sp3_repl))), 
            as.numeric(c(sd(sp4_repl))), as.numeric(c(sd(sp5_repl))), 
            as.numeric(c(sd(sp6_repl)))))

        div_tag <- data.frame(clade = c("T1", "T2", "T3", "T4", "T5", "T6"))
        organ_id <- data.frame(comp_organ = rep(organ, 6))
        df_cor_avg <- cbind(div_tag, organ_id, df_cor_avg, df_cor_error)

        return(df_cor_avg)

    }

    cereb_div <- getBrCerebCor(df=x_Br[,27:42], organ="Cerebellum", coefficient=coefficient)

    cereb_div_11 <- getBrCerebCor(df=x_Br2011[,c(6:7,36:37,24:25,47:48,58,68:69,81:83,97:98)], 
        organ="Cerebellum", coefficient=coefficient)


    getBrHtKdLvCor <- function(df, organ, coefficient) {

        df_cor <- sqrt(1/2*(1 - cor(df, method=coefficient)))
        df_cor <- df_cor[3:nrow(df_cor), ]

        sp1_repl <- c(mean(df_cor[1:2,1:2]), mean(df_cor[3:4,1:2])) # Hsa-Ppa Hsa-Ptr
        sp2_repl <- c(mean(df_cor[5:6,1:2]), mean(df_cor[5:6,3:4]), mean(df_cor[5:6,5:6])) # Hsa-Ggo Ppa-Ggo Ptr-Ggo
        sp3_repl <- c(mean(df_cor[7:8,1:2]), mean(df_cor[7:8,3:4]), mean(df_cor[7:8,5:6]), mean(df_cor[7:8,7:8])) # Hsa-Ppy Ppa-Ppy Ptr-Ppy Ggo-Ppy
        sp4_repl <- c(mean(df_cor[9:10,1:2]), mean(df_cor[9:10,3:4]), mean(df_cor[9:10,5:6]), mean(df_cor[9:10,7:8]), 
            mean(df_cor[9:10,9:10])) # Hsa-Mml Ppa-Mml Ptr-Mml Ggo-Mml Ppy-Mml
        sp5_repl <- c(mean(df_cor[11:13,1:2]), mean(df_cor[11:13,3:4]), mean(df_cor[11:13,5:6]), mean(df_cor[11:13,7:8]), 
            mean(df_cor[11:13,9:10]), mean(df_cor[11:13,11:12])) # Hsa-Mmu Ppa-Mmu Ptr-Mmu Ggo-Mmu Ppy-Mmu Mml-Mmu
        if (organ == "Kidney") {
            sp6_repl <- c(mean(df_cor[14,1:2]), mean(df_cor[14,3:4]), mean(df_cor[14,5:6]), mean(df_cor[14,7:8]), 
            mean(df_cor[14,9:10]), mean(df_cor[14,11:12]), mean(df_cor[14,13:15])) # Hsa-Mdo Ppa-Mdo Ptr-Mdo Ggo-Mdo Ppy-Mdo Mml-Mdo Mmu-Mdo
        } else {
            sp6_repl <- c(mean(df_cor[14:15,1:2]), mean(df_cor[14:15,3:4]), mean(df_cor[14:15,5:6]), mean(df_cor[14:15,7:8]), 
            mean(df_cor[14:15,9:10]), mean(df_cor[14:15,11:12]), mean(df_cor[14:15,13:15])) # Hsa-Mdo Ppa-Mdo Ptr-Mdo Ggo-Mdo Ppy-Mdo Mml-Mdo Mmu-Mdo
        }

        # Get mean and SE
        df_cor_avg <- data.frame(correlation = c(mean(sp1_repl), mean(sp2_repl), mean(sp3_repl), 
            mean(sp4_repl), mean(sp5_repl), mean(sp6_repl)))

        df_cor_error <- data.frame(error = c(as.numeric(c(sd(sp1_repl))),
            as.numeric(c(sd(sp2_repl))), as.numeric(c(sd(sp3_repl))), 
            as.numeric(c(sd(sp4_repl))), as.numeric(c(sd(sp5_repl))), 
            as.numeric(c(sd(sp6_repl)))))

        div_tag <- data.frame(clade = c("T1", "T2", "T3", "T4", "T5", "T6"))
        organ_id <- data.frame(comp_organ = rep(organ, 6))
        df_cor_avg <- cbind(div_tag, organ_id, df_cor_avg, df_cor_error)

        return(df_cor_avg)

    }

    heart_div <- getBrHtKdLvCor(df=x_Br[,43:59], organ="Heart", coefficient=coefficient)
    kidney_div <- getBrHtKdLvCor(df=x_Br[,60:75], organ="Kidney", coefficient=coefficient)
    liver_div <- getBrHtKdLvCor(df=x_Br[,76:92], organ="Liver", coefficient=coefficient)

    liver_div_11 <- getBrHtKdLvCor(df=x_Br2011[,c(14:15,42:43,30:31,53:54,63:64,74:75,90:92,103:104)], 
        organ="Liver", coefficient=coefficient)


    getBrHtKdCor <- function(df, organ, coefficient) {

        df_cor <- sqrt(1/2*(1 - cor(df, method=coefficient)))
        df_cor <- df_cor[4:nrow(df_cor), ]

        sp1_repl <- c(mean(df_cor[1:2,1:3]), mean(df_cor[3:4,1:3])) # Hsa-Ppa Hsa-Ptr
        sp2_repl <- c(mean(df_cor[5:6,1:3]), mean(df_cor[5:6,4:5]), mean(df_cor[5:6,6:7])) # Hsa-Ggo Ppa-Ggo Ptr-Ggo
        sp3_repl <- c(mean(df_cor[7:8,1:3]), mean(df_cor[7:8,4:5]), mean(df_cor[7:8,6:7]), mean(df_cor[7:8,8:9])) # Hsa-Ppy Ppa-Ppy Ptr-Ppy Ggo-Ppy
        sp4_repl <- c(mean(df_cor[9:10,1:3]), mean(df_cor[9:10,4:5]), mean(df_cor[9:10,6:7]), mean(df_cor[9:10,8:9]), 
            mean(df_cor[9:10,10:11])) # Hsa-Mml Ppa-Mml Ptr-Mml Ggo-Mml Ppy-Mml
        sp5_repl <- c(mean(df_cor[11:13,1:3]), mean(df_cor[11:13,4:5]), mean(df_cor[11:13,6:7]), mean(df_cor[11:13,8:9]), 
            mean(df_cor[11:13,10:11]), mean(df_cor[11:13,12:13])) # Hsa-Mmu Ppa-Mmu Ptr-Mmu Ggo-Mmu Ppy-Mmu Mml-Mmu
        sp6_repl <- c(mean(df_cor[14:15,1:3]), mean(df_cor[14:15,4:5]), mean(df_cor[14:15,6:7]), mean(df_cor[14:15,8:9]), 
            mean(df_cor[14:15,10:11]), mean(df_cor[14:15,12:13]), mean(df_cor[14:15,14:16])) # Hsa-Mdo Ppa-Mdo Ptr-Mdo Ggo-Mdo Ppy-Mdo Mml-Mdo Mmu-Mdo

        # Get mean and SE
        df_cor_avg <- data.frame(correlation = c(mean(sp1_repl), mean(sp2_repl), mean(sp3_repl), 
            mean(sp4_repl), mean(sp5_repl), mean(sp6_repl)))

        df_cor_error <- data.frame(error = c(as.numeric(c(sd(sp1_repl))),
            as.numeric(c(sd(sp2_repl))), as.numeric(c(sd(sp3_repl))), 
            as.numeric(c(sd(sp4_repl))), as.numeric(c(sd(sp5_repl))), 
            as.numeric(c(sd(sp6_repl)))))

        div_tag <- data.frame(clade = c("T1", "T2", "T3", "T4", "T5", "T6"))
        organ_id <- data.frame(comp_organ = rep(organ, 6))
        df_cor_avg <- cbind(div_tag, organ_id, df_cor_avg, df_cor_error)

        return(df_cor_avg)

    }

    heart_div_11 <- getBrHtKdCor(df=x_Br2011[,c(8:10,38:39,26:27,49:50,59:60,70:71,84:86,99:100)], 
        organ="Heart", coefficient=coefficient)
    kidney_div_11 <- getBrHtKdCor(df=x_Br2011[,c(11:13,40:41,28:29,51:52,61:62,72:73,87:89,101:102)], 
        organ="Kidney", coefficient=coefficient)


    getBrTestisCor <- function(df, organ, coefficient) {

        df_cor <- sqrt(1/2*(1 - cor(df, method=coefficient)))
        df_cor <- df_cor[3:nrow(df_cor), ]

        sp1_repl <- mean(df_cor[1,1:2]) # Hsa-Ppa
        sp2_repl <- c(mean(df_cor[3,1:2]), mean(df_cor[3,3]), mean(df_cor[3,4])) # Hsa-Ggo Ppa-Ggo Ptr-Ggo
        sp3_repl <- NA # orangutan: no data available
        sp4_repl <- c(mean(df_cor[4:5,1:2]), mean(df_cor[4:5,3]), mean(df_cor[4:5,4]), mean(df_cor[4:5,5])) # Hsa-Mml Ppa-Mml Ptr-Mml Ggo-Mml
        sp5_repl <- c(mean(df_cor[6:7,1:2]), mean(df_cor[6:7,3]), mean(df_cor[6:7,4]), mean(df_cor[6:7,5]), 
            mean(df_cor[6:7,6:7])) # Hsa-Mmu Ppa-Mmu Ptr-Mmu Ggo-Mmu Mml-Mmu
        sp6_repl <- c(mean(df_cor[8:9,1:2]), mean(df_cor[8:9,3]), mean(df_cor[8:9,4]), mean(df_cor[8:9,5]), 
            mean(df_cor[8:9,6:7]), mean(df_cor[8:9,8:9])) # Hsa-Mdo Ppa-Mdo Ptr-Mdo Ggo-Mdo Mml-Mdo Mmu-Mdo

        # Get mean and SE
        df_cor_avg <- data.frame(correlation = c(mean(sp1_repl), mean(sp2_repl), mean(sp3_repl), 
            mean(sp4_repl), mean(sp5_repl), mean(sp6_repl)))

        df_cor_error <- data.frame(error = c(as.numeric(c(sd(sp1_repl))),
            as.numeric(c(sd(sp2_repl))), as.numeric(c(sd(sp3_repl))), 
            as.numeric(c(sd(sp4_repl))), as.numeric(c(sd(sp5_repl))), 
            as.numeric(c(sd(sp6_repl)))))

        div_tag <- data.frame(clade = c("T1", "T2", "T3", "T4", "T5", "T6"))
        organ_id <- data.frame(comp_organ = rep(organ, 6))
        df_cor_avg <- cbind(div_tag, organ_id, df_cor_avg, df_cor_error)

        return(df_cor_avg)

    }

    testis_div <- getBrTestisCor(df=x_Br[,93:103], organ="Testis", coefficient=coefficient)

    testis_div_11 <- getBrTestisCor(df=x_Br2011[,c(16:17,44,32,55,76:77,93:94,105:106)], 
        organ="Testis", coefficient=coefficient)


    Brawand_organ_cor <- rbind(brain_div, cereb_div, heart_div, kidney_div, liver_div, testis_div)

    Brawand11_organ_cor <- rbind(brain_div_11, cereb_div_11, heart_div_11, kidney_div_11, 
        liver_div_11, testis_div_11)

    # divergence times are estimated taxon pair times from TimeTree
    # http://www.timetree.org/
    br_dataset <- data.frame(dataset = rep("Mammals(re-analyzed)", 36))
    br_div_times <- data.frame(div_times = rep(c(6.7, 9.1, 15.8, 29.4, 90, 159), times=6))

    Brawand_div_rates <- data.frame(cbind(Brawand_organ_cor, br_dataset, br_div_times), 
        stringsAsFactors=FALSE)

    Brawand_div_rates <- Brawand_div_rates %>% select(clade, comp_organ, div_times, everything())

    # Remove Orangutan testis (missing data)
    Brawand_div_rates <- Brawand_div_rates[c(-33),]

    # Combine DevSeq and Brawand GE divergence data
    compDivRates <- rbind(DevSeq_div_rates, Brawand_div_rates)


    if (is.element("TPM", expr_estimation)) {

        Brawand_sou_v_div_rates <- data.frame(cbind(Br_sou_v, br_dataset, br_div_times), 
            stringsAsFactors=FALSE)

        Brawand_sou_v_div_rates <- Brawand_sou_v_div_rates %>% select(clade, comp_organ, div_times, everything())

        # Remove Orangutan testis (missing data)
        Brawand_sou_v_div_rates <- Brawand_sou_v_div_rates[c(-33),]

        # Combine DevSeq and Brawand GE divergence data
        compSouVDivRates <- rbind(DevSeq_sou_v_div_rates, Brawand_sou_v_div_rates)

    }


    # Reshape original Brawand data table
    br11_dataset <- data.frame(dataset = rep("Mammals", 36))

    Brawand11_div_rates <- data.frame(cbind(Brawand11_organ_cor, br11_dataset, br_div_times), 
        stringsAsFactors=FALSE)

    Brawand11_div_rates <- Brawand11_div_rates %>% select(clade, comp_organ, div_times, everything())

    # Remove Orangutan testis (missing data)
    Brawand11_div_rates <- Brawand11_div_rates[c(-33),]

    # Combine DevSeq and Brawand 2011 GE divergence data
    compDivRates11 <- rbind(DevSeq_div_rates, Brawand11_div_rates)


    if (is.element("TPM", expr_estimation)) {

        Brawand11_sou_v_div_rates <- data.frame(cbind(Br2011_sou_v, br11_dataset, br_div_times), 
            stringsAsFactors=FALSE)

        Brawand11_sou_v_div_rates <- Brawand11_sou_v_div_rates %>% select(clade, comp_organ, div_times, everything())

        # Remove Orangutan testis (missing data)
        Brawand11_sou_v_div_rates <- Brawand11_sou_v_div_rates[c(-33),]

        # Combine DevSeq and Brawand 2011 GE divergence data
        compSouVDivRates11 <- rbind(DevSeq_sou_v_div_rates, Brawand11_sou_v_div_rates)

    }




#---- Apply LOESS regression to sOU and pearson dist expression data and compare slopes -----

  if(expr_estimation == "TPM") {


    getLOESS.Coord <- function(organ_data) {

      comp_organ <- unique(organ_data$comp_organ)

      if (comp_organ == "Root") {

        poly_deg <- 1

      } else { poly_deg <- 2 } # Use quadratic polynomes for all organs except root

        temp <- loess.smooth(organ_data$div_times, organ_data$correlation, span = 1, 
          degree = poly_deg, family="gaussian", evaluation = 200)

      # Obtain coordinates of the smooth curve
      div_times <- as.numeric(unlist(temp$x))
      correlation <- as.numeric(unlist(temp$y))
      comp_organ <- rep(comp_organ, 200)

      experiment <- unique(organ_data$dataset)
      experiment <- rep(experiment, 200)

      loess_out <- data.frame(comp_organ = comp_organ, div_times = div_times, correlation = correlation, 
        dataset = experiment)

      return(loess_out)
    }


    # Set up lists containing sOU expression distances
    devseqSouV_organ_lst <- list(compSouVDivRates11[1:6,], compSouVDivRates11[7:12,], compSouVDivRates11[13:18,], 
      compSouVDivRates11[19:24,], compSouVDivRates11[25:30,], compSouVDivRates11[31:36,], compSouVDivRates11[37:42,], 
      compSouVDivRates11[43:48,])

    devseqSouV_organ_lst_sel <- list(compSouVDivRates11[1:6,], compSouVDivRates11[7:11,], compSouVDivRates11[13:18,], 
      compSouVDivRates11[19:24,], compSouVDivRates11[25:30,], compSouVDivRates11[31:36,], compSouVDivRates11[37:42,], 
      compSouVDivRates11[43:48,])

    brawandSouV11_organ_lst <- list(compSouVDivRates11[49:54,], compSouVDivRates11[55:60,], compSouVDivRates11[61:66,], 
      compSouVDivRates11[67:72,],compSouVDivRates11[73:78,], compSouVDivRates11[79:83,])

    brawandSouV_organ_lst <- list(compSouVDivRates[49:54,], compSouVDivRates[55:60,], compSouVDivRates[61:66,], 
      compSouVDivRates[67:72,], compSouVDivRates[73:78,], compSouVDivRates[79:83,])


    # Set up lists containing metric pearson expression distances
    devseq_organ_lst <- list(compDivRates11[1:6,], compDivRates11[7:12,], compDivRates11[13:18,], 
      compDivRates11[19:24,], compDivRates11[25:30,], compDivRates11[31:36,], compDivRates11[37:42,], 
      compDivRates11[43:48,])

    devseq_organ_lst_sel <- list(compDivRates11[1:6,], compDivRates11[7:11,], compDivRates11[13:18,], 
      compDivRates11[19:24,], compDivRates11[25:30,], compDivRates11[31:36,], compDivRates11[37:42,], 
      compDivRates11[43:48,])

    brawand11_organ_lst <- list(compDivRates11[49:54,], compDivRates11[55:60,], compDivRates11[61:66,], 
      compDivRates11[67:72,],compDivRates11[73:78,], compDivRates11[79:83,])

    brawand_organ_lst <- list(compDivRates[49:54,], compDivRates[55:60,], compDivRates[61:66,], 
      compDivRates[67:72,], compDivRates[73:78,], compDivRates[79:83,])


    # Get LOESS coordinates for DevSeq and Brawand data
    DevSeqSouV_loess_coord <- as.data.frame(do.call(rbind, lapply(devseqSouV_organ_lst, getLOESS.Coord)))
    Brawand11SouV_loess_coord <- as.data.frame(do.call(rbind, lapply(brawandSouV11_organ_lst, getLOESS.Coord)))
    loessSouV_coor11_AT <- rbind(DevSeqSouV_loess_coord, Brawand11SouV_loess_coord)


    # Format loess df for ggplot2
    formatLOESS.DF <- function(df) {

      df$div_times <- as.numeric(df$div_times)
      df$correlation <- as.numeric(df$correlation)
      df$dataset <- as.factor(df$dataset)

      return(df)
    }

    loessSouV_coor11_AT <- formatLOESS.DF(loessSouV_coor11_AT)



    getLOESS.Slopes <- function(organ_data) {

      comp_organ <- unique(organ_data$comp_organ)

      if (comp_organ == "Root") {

        poly_deg <- 1
        alpha <- 1

      } else {

        poly_deg <- 2
        alpha <- 1

      } # Use quadratic polynomes for all organs except root

      temp <- loess.smooth(organ_data$div_times, organ_data$correlation, span = alpha, 
        degree = poly_deg, family="gaussian", evaluation = 200)

      # Get slope values
      slopes = diff(temp$y)/diff(temp$x)
      slopes_avg <- mean(slopes)
      slopes_avg <- as.numeric(as.data.frame(slopes_avg))

      return(slopes_avg)
    }


    # Get LOESS slopes for DevSeq and Brawand data
    # For sOU expression distances
    DevSeqSouV_AT_loess_slopes <- as.data.frame(do.call(rbind, lapply(devseqSouV_organ_lst, getLOESS.Slopes)))
    DevSeqSouV_sel_AT_loess_slopes <- as.data.frame(do.call(rbind, lapply(devseqSouV_organ_lst_sel, getLOESS.Slopes))) ## hypocotyl slope is 0.00797 instead 0.00934 if BD is left out
    Brawand11SouV_loess_slopes <- as.data.frame(do.call(rbind, lapply(brawandSouV11_organ_lst, getLOESS.Slopes)))
    sOU_loess_DevSeq_AT_Br11_wilcox <- wilcox.test(as.numeric(unlist(DevSeqSouV_sel_AT_loess_slopes)), as.numeric(unlist(Brawand11SouV_loess_slopes)))$p.value

    BrawandSouV_loess_slopes <- as.data.frame(do.call(rbind, lapply(brawandSouV_organ_lst, getLOESS.Slopes)))
    sOU_loess_Br_Br11_wilcox <- wilcox.test(as.numeric(unlist(Brawand11SouV_loess_slopes)), as.numeric(unlist(BrawandSouV_loess_slopes)))$p.value

    # For metric pearson expression distances
    DevSeq_AT_loess_slopes <- as.data.frame(do.call(rbind, lapply(devseq_organ_lst, getLOESS.Slopes)))
    DevSeq_AT_sel_loess_slopes <- as.data.frame(do.call(rbind, lapply(devseq_organ_lst_sel, getLOESS.Slopes)))
    Brawand11_loess_slopes <- as.data.frame(do.call(rbind, lapply(brawand11_organ_lst, getLOESS.Slopes)))

    Brawand_loess_slopes <- as.data.frame(do.call(rbind, lapply(brawand_organ_lst, getLOESS.Slopes)))


    # Write slope values to csv file
    DevSeq_slopes <- cbind(DevSeqSouV_sel_AT_loess_slopes, DevSeqSouV_sel_AT_loess_slopes, DevSeq_AT_loess_slopes, DevSeq_AT_loess_slopes)
    colnames(DevSeq_slopes) <- c("DS_AT_Br11_sOU_loess", "DS_AT_Br_sOU_loess", "DS_AT_Br11_pea_loess", "DS_AT_Br_pea_loess")
    ds_organs <- data.frame(sample = c("Root", "Hypocotyl", "Leaf", "Apex_veg", "Apex_inf", "Flower", "Stamen", "Carpel"))
    DevSeq_slopes <- cbind(ds_organs, DevSeq_slopes)

    Brawand_slopes <- cbind(Brawand11SouV_loess_slopes, BrawandSouV_loess_slopes, Brawand11_loess_slopes, Brawand_loess_slopes)
    colnames(Brawand_slopes) <- c("DS_AT_Br11_sOU_loess", "DS_AT_Br_sOU_loess", "DS_AT_Br11_pea_loess", "DS_AT_Br_pea_loess")
    br_organs <- data.frame(sample = c("Brain", "Cerebellum", "Heart", "Kidney", "Liver", "Testis"))
    Brawand_slopes <- cbind(br_organs, Brawand_slopes)

    DS_AT_Br_loess_slopes <- rbind(DevSeq_slopes, Brawand_slopes)


    # Get Welch test and Wilcox test p values
    getLoessStats <- function(x) {

      Welch_test <- t.test(x[1:8], x[9:14])$p.value
      Wilcox_test <- wilcox.test(x[1:8], x[9:14])$p.value

      stats_out <- rbind(Welch_test, Wilcox_test)
      return(stats_out)
    }


    DS_AT_Br_loess_lst <- list(DS_AT_Br_loess_slopes[,2], DS_AT_Br_loess_slopes[,3], 
      DS_AT_Br_loess_slopes[,4], DS_AT_Br_loess_slopes[,5])

    DS_AT_Br_loess_stats <- as.data.frame(do.call(cbind, lapply(DS_AT_Br_loess_lst, getLoessStats)))
    statsTests <- c("Welch_test", "Wilcox_test")
    DS_AT_Br_loess_stats <- cbind(statsTests, DS_AT_Br_loess_stats)
    colnames(DS_AT_Br_loess_stats) <- colnames(DS_AT_Br_loess_slopes)
    rownames(DS_AT_Br_loess_stats) <- NULL
    DS_AT_Br_loess_slopes <- rbind(DS_AT_Br_loess_slopes, DS_AT_Br_loess_stats)



    # Show message
    message("Writing data tables...")

    write.table(DS_AT_Br_loess_slopes, 
      file=file.path(out_dir, "output", "data", "DS_AT_Br_loess_slopes.txt"), sep="\t", 
      col.names=TRUE, row.names=FALSE, dec=".", quote = FALSE)


    # Create p-value containing test strings for plots
    p_dat_text <- data.frame(
        label = c(paste("italic('P =')~", formatC(sOU_loess_DevSeq_AT_Br11_wilcox, format = "e", digits = 0))),
        x = c(16.25),
        y = c(1.477)
    )



   # Make sOU GE divergence plot showing individual organ regressions for SI
   makeOrgRegPlot <- function(data1, data2, coefficient, expr_estimation, pos) {

      fname <- sprintf('%s.jpg', paste("compSouVDivRates11_loess", expr_estimation, pos, sep="_"))

      ymin <- 0.05
      ymax <- 1.45

      if (pos == "main") {

            plot_wdt <- 12.535
            plot_hdt <- 8
            legend_x_pos <- 0.2427
            legend_y_pos <- 0.914
            linewd <- 3.2
            point_size <- 5.85

      } else {

            plot_wdt <- 9.5 # condenced plot width for suppl
            plot_hdt <- 6.75 # condenced plot width for suppl
            legend_x_pos <- 0.317
            legend_y_pos <- 0.9
            linewd <- 2.5
            point_size <- 5
      }

      ds_col <- rep(c("#677ebc"), 48)
      bw_col <- rep(c("red"), 35)
      col_scale <- c("#677ebc", "red")
      fill_scale <- c("#677ebc", "red")
      col_breaks <- c("Angiosperms ", "Mammals")
      fill_breaks <- c("Angiosperms ", "Mammals")
      shape_scale <- c(16, 17)

      fill_col <- c(as.character(ds_col), as.character(bw_col))

      p <- ggplot() 
      p <- p + geom_line(size = linewd, data = data1, aes(x = div_times, y = correlation, group = comp_organ, 
        colour = dataset, linetype = dataset))
      p <- p + geom_point(size = point_size, data = data2, aes(x = div_times, y = correlation, group = comp_organ, 
        colour = dataset, shape = dataset))
      p <- p + geom_line(size = linewd, data = data1[2601:2800,], aes(x = div_times, y = correlation, group = comp_organ, 
        colour = dataset, linetype = dataset)) # plot testis data over angiosperm data for better visibility
      p <- p + geom_point(size = point_size, data = data2[17,], aes(x = div_times, y = correlation, group = comp_organ, 
        colour = dataset, shape = dataset)) # plot leaf MT data point over line for better visibility
      p <- p + scale_x_continuous(limits = c(0,161), expand = c(0.02,0), breaks = c(0,20,40,60,80,100,120,140,160)) + 
      scale_y_continuous(limits = c(0.055, 1.763), expand = c(0.02, 0), breaks = c(0.2,0.4,0.6,0.8,1,1.2,1.4,1.6)) + 
      scale_color_manual(values = col_scale, breaks = col_breaks) + 
      scale_fill_manual(values = fill_scale, breaks = fill_breaks) + 
      scale_shape_manual(values = shape_scale) + 
      scale_size(range = c(0.5, 12)) + 
      geom_text(data = p_dat_text, mapping = aes(x = x, y = y, label = format(label, scientific=TRUE)), 
        size = 8, parse = TRUE) + 
      guides(color = guide_legend(ncol=2, keywidth = 0.4, keyheight = 0.4, default.unit = "inch"))

      q <- p + theme_bw() + xlab("Divergence time (Myr)") + ylab("Expression distance") + 
      theme(text=element_text(size=16), 
        axis.ticks.length=unit(0.35, "cm"), 
        axis.ticks = element_line(colour = "black", size = 0.7),  
        plot.margin = unit(c(0.55, 1.175, 0.5, 0.4),"cm"), 
        axis.title.y = element_text(size=25, margin = margin(t = 0, r = 14.5, b = 0, l = 11.5), colour="black"), 
        axis.title.x = element_text(size=25, margin = margin(t = 13.25, r = 0, b = 3.5, l = 0), colour="black"), 
        axis.text.x = element_text(size=21.25, angle=0, margin = margin(t = 5.5), colour="black"), 
        axis.text.y = element_text(size=21.25, angle=0, margin = margin(r = 5.5), colour="black"), 
        legend.box.background = element_rect(colour = "#d5d5d5", fill=NA, size=1.0), 
        panel.border = element_rect(colour = "black", fill=NA, size=1.75), 
        panel.grid = element_blank(), 
        legend.position = c(legend_x_pos, legend_y_pos), 
        legend.title = element_blank(), 
        legend.text = element_text(size=22), 
        legend.spacing.x = unit(0.5, 'cm'), 
        legend.key.size = unit(0.95, "cm"), 
        legend.margin = margin(c(5.75,8,5.75,5.75)), 
        legend.background=element_blank()) 

        ggsave(file = file.path(out_dir, "output", "plots", fname), plot = q, 
            width = plot_wdt, height = plot_hdt, dpi = 300, units = c("in"), limitsize = FALSE) 
  }

  # makeOrgRegPlot(data1 = loessSouV_coor11_AT, data2 = compSouVDivRates11, coefficient = coefficient, 
  #  expr_estimation = expr_estimation, pos = "ext")

  makeOrgRegPlot(data1 = loessSouV_coor11_AT, data2 = compSouVDivRates11, coefficient = coefficient, 
    expr_estimation = expr_estimation, pos = "main")


  }




#---- Make gene expression divergence rates plot for Brawand data (original and re-analyzed) ---


   makeGEDivPlotBr <- function(data) {

      fname <- sprintf('%s.jpg', paste(deparse(substitute(data)), coefficient, expr_estimation, sep="_"))

      if (deparse(substitute(data)) == "Brawand_div_rates") {
        y_min <- 0.24
        y_max <- 0.48
        yseg_min <- 0.14
        yseg_max <- 0.2415
        y_text <- 0.25

      } else if (deparse(substitute(data)) == "Brawand11_div_rates") {
        y_min <- 0.17
        y_max <- 0.461
        yseg_min <- 0.12
        yseg_max <- 0.1715
        y_text <- 0.1822
      }

      p <- ggplot(data=data, aes(x=div_times, y=correlation, group=comp_organ, colour=comp_organ)) + 
      geom_line(size = 3) +  
      scale_x_continuous(limits = c(6,160), expand = c(0.02,0), breaks = c(7,9,16,29,90,159)) + 
      scale_y_continuous(limits = c(y_min, y_max), expand = c(0.02, 0)) + 
      annotate("text", x=15.5, y=y_text, label= "Primates", size=8) + 
      annotate("text", x=90, y=y_text, label= "Mouse", size=8) + 
      annotate("text", x=151, y=y_text, label= "Opossum", size=8) + 
      geom_segment(x=7, xend=7, y=yseg_min, yend=yseg_max, color="black", size=0.7) + 
      geom_segment(x=9, xend=9, y=yseg_min, yend=yseg_max, color="black", size=0.7) + 
      geom_segment(x=16, xend=16, y=yseg_min, yend=yseg_max, color="black", size=0.7) + 
      geom_segment(x=29, xend=29, y=yseg_min, yend=yseg_max, color="black", size=0.7) + 
      geom_segment(x=90, xend=90, y=yseg_min, yend=yseg_max, color="black", size=0.7) + 
      geom_segment(x=159, xend=159, y=yseg_min, yend=yseg_max, color="black", size=0.7) + 
      guides(color = guide_legend(ncol = 3))

      q <- p + theme_bw() + xlab("Divergence time from HSA (Myr)") + ylab("Pearson distance") + 
      theme(text=element_text(size=16), 
        axis.ticks.length=unit(0.35, "cm"), 
        axis.ticks = element_line(colour = "black", size = 0.7),  
        plot.margin = unit(c(0.55, 1.175, 0.5, 0.4),"cm"), 
        axis.title.y = element_text(size=25, margin = margin(t = 0, r = 17, b = 0, l = 9), colour="black"), 
        axis.title.x = element_text(size=25, margin = margin(t = 14.75, r = 0, b = 2, l = 0), colour="black"), 
        axis.text.x = element_text(size=21.25, angle=0, margin = margin(t = 5.5), colour="black"), 
        axis.text.y = element_text(size=21.25, angle=0, margin = margin(r = 5.5), colour="black"), 
        legend.box.background = element_rect(colour = "#d5d5d5", fill=NA, size=1.0), 
        panel.border = element_rect(colour = "black", fill=NA, size=1.75), 
        panel.grid.major = element_line(color="#d5d5d5"),
        panel.grid.minor.x = element_blank(), 
        panel.grid.minor.y = element_blank(), 
        legend.position = c(0.244, 0.889), 
        legend.title = element_blank(), 
        legend.text = element_text(size=21.5), 
        legend.spacing.x = unit(0.5, 'cm'), 
        legend.key.size = unit(0.95, "cm"), 
        legend.background=element_blank()) 

      ggsave(file = file.path(out_dir, "output", "plots", fname), plot = q, 
        width = 12.535, height = 8, dpi = 300, units = c("in"), limitsize = FALSE) 
  }

  makeGEDivPlotBr(data = Brawand_div_rates)
  makeGEDivPlotBr(data = Brawand11_div_rates)


  # Write gene expression divergence rate to file
  # For both metric Pearson distance and sOU-v expression distance

  # Show message
  message("Writing data tables...")

  div_rates_list <- list(compDivRates = compDivRates, compDivRates11 = compDivRates11, 
    compSouVDivRates = compSouVDivRates, compSouVDivRates11 = compSouVDivRates11)

  for(i in names(div_rates_list)){
    write.table(div_rates_list[[i]], file=file.path(out_dir, "output", "data", paste0(i,".txt")), 
        sep="\t", col.names=TRUE, row.names=FALSE, dec=".", quote = FALSE)
  }

}






