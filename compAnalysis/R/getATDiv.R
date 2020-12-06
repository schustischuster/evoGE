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

        x_Br2011_all_taxa_objects = TEconstruct(ExpValueFP = file.path(out_dir, 
            "output", "data", 'x_Br2011_all_taxobj_input.txt'), taxa = "all", subtaxa = 'all')

        x_DS_taxa_objects = TEconstruct(ExpValueFP = file.path(out_dir, 
            "output", "data", 'x_DS_taxobj_input.txt'), taxa = "all", subtaxa = 'all')



        Brawand_organ_list <- list("brain", "cerebellum", "heart", "kidney", "liver", "testis")

        Brawand2011_organ_list <- list("br", "cb", "ht", "kd", "lv", "ts")

        DevSeq_organ_list <- list("Root", "Hypocotyl", "Leaf", "vegApex", "infApex", "Flower", 
            "Stamen", "Carpel", "Pollen")



        # Function to apply extended OU model with dynamic expression optimum ("variable-µ method")
        getExtOU <- function(organ, taxa_obj, samples) {

            sou_v_out <- expdist(taxa_obj, taxa = "all",
                subtaxa = organ,
                method = "sou_v")

            sou_v_pi <- sou_v_out$pi
            sou_v_distance <- as.data.frame(sou_v_out$distance)
            sou_v_distance_div <- as.data.frame(sou_v_distance[,1])
            sou_v_distance_div <- rbind(sou_v_distance_div, sou_v_pi)

            spec_id <- c(sub("\\_.*", "", rownames(sou_v_distance)), "pi")
            organ_id <- unique(sub(".*_", "", rownames(sou_v_distance)))

            if ((organ_id == "ts") && (samples == "sel")) {

                ppy_testis <- "NA"
                sou_v_distance_div <- as.data.frame(c(sou_v_distance_div[1:3,], ppy_testis, 
                    sou_v_distance_div[4:7,]), stringsAsFactors = FALSE)
                spec_id <- c("hsa", "ppa", "ggo", "ppy", "mml", "mmu", "mdo", "pi")
            
            } 

            if ((organ_id == "ts") && (samples == "all")) {

                ppy_testis <- "NA"
                sou_v_distance_div <- as.data.frame(c(sou_v_distance_div[1:3,], ppy_testis, 
                    sou_v_distance_div[4:9,]), stringsAsFactors = FALSE)
                spec_id <- c("hsa", "ppa", "ggo", "ppy", "mml", "mmu", "mdo", "oan", "gga", "pi")

            }

            if (organ_id == "testis") {

                Orangutan_testis <- "NA"
                sou_v_distance_div <- as.data.frame(c(sou_v_distance_div[1:3,], Orangutan_testis, 
                    sou_v_distance_div[4:7,]), stringsAsFactors = FALSE)
                spec_id <- c("Human", "Bonobo", "Gorilla", "Orangutan", "Macaque", "Mouse", 
                    "Opossum", "pi")
            }

            rownames(sou_v_distance_div) <- spec_id
            colnames(sou_v_distance_div) <- organ_id

            return(sou_v_distance_div)
        }


        Br_sou_v <- as.data.frame(do.call(cbind, lapply(Brawand_organ_list, getExtOU,
            taxa_obj = x_Br_taxa_objects)))

        rows_to_remove_Br <- "Human"
        Br_sou_v <- Br_sou_v[!(row.names(Br_sou_v) %in% rows_to_remove_Br), ]
        Br_sou_v$testis <- suppressWarnings(as.numeric(Br_sou_v$testis))


        Br2011_sou_v <- as.data.frame(do.call(cbind, lapply(Brawand2011_organ_list, getExtOU,
            taxa_obj = x_Br2011_taxa_objects, samples = "sel")))
        rows_to_remove_Br2011 <- "hsa"
        Br2011_sou_v <- Br2011_sou_v[!(row.names(Br2011_sou_v) %in% rows_to_remove_Br2011), ]
        Br2011_sou_v$ts <- suppressWarnings(as.numeric(Br2011_sou_v$ts))


        Br2011_all_sou_v <- as.data.frame(do.call(cbind, lapply(Brawand2011_organ_list, getExtOU,
            taxa_obj = x_Br2011_all_taxa_objects, samples = "all")))
        rows_to_remove_Br2011_all <- "hsa"
        Br2011_all_sou_v <- Br2011_all_sou_v[!(row.names(Br2011_all_sou_v) %in% rows_to_remove_Br2011_all), ]
        Br2011_all_sou_v$ts <- suppressWarnings(as.numeric(Br2011_all_sou_v$ts))


        DS_sou_v <- as.data.frame(do.call(cbind, lapply(DevSeq_organ_list, getExtOU,
            taxa_obj = x_DS_taxa_objects)))
        rows_to_remove_DS <- "ATH"
        DS_sou_v <- DS_sou_v[!(row.names(DS_sou_v) %in% rows_to_remove_DS), ]
        # Note: The computation for THA and BDY Pollen distance will produce NaN due to large
        # number of genes having only "0" TPM values in this combination; We can ignore this warning
        # since we carry out downstream analyses only for ATH vs.X and ALY vs.X combinations



        # Reshape data for ggplot2
        Brawand_sou_v_div <- as.data.frame(c(Br_sou_v[1:6, 1], Br_sou_v[1:6, 2], Br_sou_v[1:6, 3], 
            Br_sou_v[1:6, 4], Br_sou_v[1:6, 5], Br_sou_v[1:6, 6]))
        colnames(Brawand_sou_v_div) <- "correlation"

        Brawand2011_sou_v_div <- as.data.frame(c(Br2011_sou_v[1:6, 1], Br2011_sou_v[1:6, 2], 
            Br2011_sou_v[1:6, 3], Br2011_sou_v[1:6, 4], Br2011_sou_v[1:6, 5], Br2011_sou_v[1:6, 6]))
        colnames(Brawand2011_sou_v_div) <- "correlation"

        Brawand2011_all_sou_v_div <- as.data.frame(c(Br2011_all_sou_v[1:8, 1], Br2011_all_sou_v[1:8, 2], 
            Br2011_all_sou_v[1:8, 3], Br2011_all_sou_v[1:8, 4], Br2011_all_sou_v[1:8, 5], Br2011_all_sou_v[1:8, 6]))
        colnames(Brawand2011_all_sou_v_div) <- "correlation"

        DevSeq_sou_v_div <- as.data.frame(c(DS_sou_v[1:6, 1], DS_sou_v[1:6, 2], DS_sou_v[1:6, 3], 
            DS_sou_v[1:6, 4], DS_sou_v[1:6, 5], DS_sou_v[1:6, 6], DS_sou_v[1:6, 7], DS_sou_v[1:6, 8]))
        colnames(DevSeq_sou_v_div) <- "correlation"

    }




#---------------- Get gene expression divergence rates for ATH/AL vs species X -----------------


    # Use pearson correlation, inter-organ normalization and TPM for ms

    getDSOrganCor <- function(df, organ, coefficient) {

        df_cor <- sqrt(1 - cor(df, method=coefficient))
        df_cor <- df_cor[4:nrow(df_cor), 1:3]

        # Reshape cor data frame to one column
        df_cor_rs <- data.frame(newcol = c(t(df_cor)), stringsAsFactors=FALSE)

        sp1 <- mean(df_cor_rs[1:9,])
        sp2 <- mean(df_cor_rs[10:18,])
        sp3 <- mean(df_cor_rs[19:27,])
        sp4 <- mean(df_cor_rs[28:36,])
        sp5 <- mean(df_cor_rs[37:45,])
        sp6 <- mean(df_cor_rs[46:54,])

        df_cor_avg <- rbind(sp1, sp2, sp3, sp4, sp5, sp6)
        colnames(df_cor_avg) <- organ

        getRowNames = function(x,n){ substring(x,nchar(x)-n+1) }
        row_names_repl <- getRowNames(rownames(df_cor),2)
        rnames_div_rates <- unique(row_names_repl)
        rownames(df_cor_avg) <- rnames_div_rates

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

    DevSeq_organ_cor <- cbind(root_div, hypocotyl_div, leaf_div, veg_apex_div, inf_apex_div, 
        flower_div, stamen_div, carpel_div)


    # Reshape data table for ggplot
    # divergence times are estimated taxon pair times from TimeTree
    # http://www.timetree.org/
    div_times <- rep(c(7.1, 9.4, 25.6, 46, 106, 160), times=8)
    comp_organ <- rep(colnames(DevSeq_organ_cor), each=6)
    comp_spec <- rep(rownames(DevSeq_organ_cor), times=8)
    dataset <- rep("Angiosperms ", 48)

    DevSeq_GE_div <- rbind(root_div, hypocotyl_div, leaf_div, veg_apex_div, inf_apex_div, 
        flower_div, stamen_div, carpel_div)
    rownames(DevSeq_GE_div) <- NULL
    colnames(DevSeq_GE_div) <- "correlation"

    DevSeq_div_rates <- data.frame(cbind(comp_spec, comp_organ, div_times, DevSeq_GE_div, dataset), 
        stringsAsFactors=FALSE)

    DevSeq_div_rates$div_times <- as.numeric(DevSeq_div_rates$div_times)
    DevSeq_div_rates$correlation <- as.numeric(DevSeq_div_rates$correlation)
      
    DevSeq_div_rates$comp_organ <- factor(DevSeq_div_rates$comp_organ, 
        levels = unique(DevSeq_div_rates$comp_organ))


    if (is.element("TPM", expr_estimation)) {

        DevSeq_sou_v_div_rates <- data.frame(cbind(comp_spec, comp_organ, div_times, DevSeq_sou_v_div, dataset), 
            stringsAsFactors=FALSE)

        DevSeq_sou_v_div_rates$div_times <- as.numeric(DevSeq_sou_v_div_rates$div_times)
        DevSeq_sou_v_div_rates$correlation <- as.numeric(DevSeq_sou_v_div_rates$correlation)

        DevSeq_sou_v_div_rates$comp_organ <- factor(DevSeq_sou_v_div_rates$comp_organ, 
            levels = unique(DevSeq_sou_v_div_rates$comp_organ))

    }




#---------------- Get gene expression divergence rates for Human vs species X ------------------


    # Use pearson correlation, inter-organ normalization and TPM for ms

    getBrBrainCor <- function(df, organ, coefficient) {

        df_cor <- sqrt(1 - cor(df, method=coefficient))
        df_cor <- df_cor[5:nrow(df_cor), 1:4]

        # Reshape cor data frame to one column
        df_cor_rs <- data.frame(newcol = c(t(df_cor)), stringsAsFactors=FALSE)

        Ppa <- mean(df_cor_rs[1:12,]) # bonobo
        Pan <- mean(df_cor_rs[13:36,]) # chimp
        Ggo <- mean(df_cor_rs[37:44,]) # gorilla
        Ppy <- mean(df_cor_rs[45:52,]) # orangutan
        Mml <- mean(df_cor_rs[53:64,]) # macaque
        Mmu <- mean(df_cor_rs[65:76,]) # mouse
        Mdo <- mean(df_cor_rs[77:84,]) # opossum

        df_cor_avg <- rbind(Ppa, Ggo, Ppy, Mml, Mmu, Mdo)
        colnames(df_cor_avg) <- organ

        return(df_cor_avg)

    }

    brain_div <- getBrBrainCor(df=x_Br[,2:26], organ="Brain", coefficient=coefficient)


    getBrCerebCor <- function(df, organ, coefficient) {

        df_cor <- sqrt(1 - cor(df, method=coefficient))
        df_cor <- df_cor[3:nrow(df_cor), 1:2]

        # Reshape cor data frame to one column
        df_cor_rs <- data.frame(newcol = c(t(df_cor)), stringsAsFactors=FALSE)

        Ppa <- mean(df_cor_rs[1:4,]) # bonobo
        Pan <- mean(df_cor_rs[5:8,]) # chimp
        Ggo <- mean(df_cor_rs[9:12,]) # gorilla
        Ppy <- mean(df_cor_rs[13:14,]) # orangutan
        Mml <- mean(df_cor_rs[15:18,]) # macaque
        Mmu <- mean(df_cor_rs[19:24,]) # mouse
        Mdo <- mean(df_cor_rs[25:28,]) # opossum

        df_cor_avg <- rbind(Ppa, Ggo, Ppy, Mml, Mmu, Mdo)
        colnames(df_cor_avg) <- organ

        return(df_cor_avg)

    }

    cereb_div <- getBrCerebCor(df=x_Br[,27:42], organ="Cerebellum", coefficient=coefficient)


    getBrHtKdLvCor <- function(df, organ, coefficient) {

        df_cor <- sqrt(1 - cor(df, method=coefficient))
        df_cor <- df_cor[3:nrow(df_cor), 1:2]

        # Reshape cor data frame to one column
        df_cor_rs <- data.frame(newcol = c(t(df_cor)), stringsAsFactors=FALSE)

        Ppa <- mean(df_cor_rs[1:4,]) # bonobo
        Pan <- mean(df_cor_rs[5:8,]) # chimp
        Ggo <- mean(df_cor_rs[9:12,]) # gorilla
        Ppy <- mean(df_cor_rs[13:16,]) # orangutan
        Mml <- mean(df_cor_rs[17:20,]) # macaque
        Mmu <- mean(df_cor_rs[21:26,]) # mouse
        if(organ == "Kidney") {
        	Mdo <- mean(df_cor_rs[27:28,]) # opossum
        } else {
        	Mdo <- mean(df_cor_rs[27:30,]) # opossum
        }

        df_cor_avg <- rbind(Ppa, Ggo, Ppy, Mml, Mmu, Mdo)
        colnames(df_cor_avg) <- organ

        return(df_cor_avg)

    }

    heart_div <- getBrHtKdLvCor(df=x_Br[,43:59], organ="Heart", coefficient=coefficient)
    kidney_div <- getBrHtKdLvCor(df=x_Br[,60:75], organ="Kidney", coefficient=coefficient)
    liver_div <- getBrHtKdLvCor(df=x_Br[,76:92], organ="Liver", coefficient=coefficient)


    getBrTestisCor <- function(df, organ, coefficient) {

        df_cor <- sqrt(1 - cor(df, method=coefficient))
        df_cor <- df_cor[3:nrow(df_cor), 1:2]

        # Reshape cor data frame to one column
        df_cor_rs <- data.frame(newcol = c(t(df_cor)), stringsAsFactors=FALSE)

        Ppa <- mean(df_cor_rs[1:2,]) # bonobo
        Pan <- mean(df_cor_rs[3:4,]) # chimp
        Ggo <- mean(df_cor_rs[5:6,]) # gorilla
        Ppy <- NA # orangutan: no data available
        Mml <- mean(df_cor_rs[7:10,]) # macaque
        Mmu <- mean(df_cor_rs[11:14,]) # mouse
        Mdo <- mean(df_cor_rs[15:18,]) # opossum

        df_cor_avg <- rbind(Ppa, Ggo, Ppy, Mml, Mmu, Mdo)
        colnames(df_cor_avg) <- organ

        return(df_cor_avg)

    }

    testis_div <- getBrTestisCor(df=x_Br[,93:103], organ="Testis", coefficient=coefficient)

    Brawand_organ_cor <- cbind(brain_div, cereb_div, heart_div, kidney_div, liver_div, testis_div)


    # Reshape data table for ggplot
    # divergence times are estimated taxon pair times from TimeTree
    # http://www.timetree.org/
    div_times <- rep(c(6.7, 9.1, 15.8, 29.4, 90, 159), times=6)
    comp_organ <- rep(colnames(Brawand_organ_cor), each=6)
    comp_spec <- rep(rownames(Brawand_organ_cor), times=6)
    dataset <- rep("Mammals(re-analyzed)", 36)

    Brawand_GE_div <- rbind(brain_div, cereb_div, heart_div, kidney_div, liver_div, testis_div)
    rownames(Brawand_GE_div) <- NULL
    colnames(Brawand_GE_div) <- "correlation"

    Brawand_div_rates <- data.frame(cbind(comp_spec, comp_organ, div_times, Brawand_GE_div, dataset), 
        stringsAsFactors=FALSE)

    Brawand_div_rates$div_times <- as.numeric(Brawand_div_rates$div_times)
    Brawand_div_rates$correlation <- as.numeric(Brawand_div_rates$correlation)

    # Remove Orangutan testis (missing data)
    Brawand_div_rates <- Brawand_div_rates[c(-33),]

    Brawand_div_rates$comp_organ <- factor(Brawand_div_rates$comp_organ, 
        levels = unique(Brawand_div_rates$comp_organ))


    # Combine DevSeq and Brawand GE divergence data
    compDivRates <- rbind(DevSeq_div_rates, Brawand_div_rates)


    if (is.element("TPM", expr_estimation)) {

        Brawand_sou_v_div_rates <- data.frame(cbind(comp_spec, comp_organ, div_times, 
        Brawand_sou_v_div, dataset), stringsAsFactors=FALSE)

        Brawand_sou_v_div_rates$div_times <- as.numeric(Brawand_sou_v_div_rates$div_times)
        Brawand_sou_v_div_rates$correlation <- as.numeric(Brawand_sou_v_div_rates$correlation)

        # Remove Orangutan testis (missing data)
        Brawand_sou_v_div_rates <- Brawand_sou_v_div_rates[c(-33),]

        Brawand_sou_v_div_rates$comp_organ <- factor(Brawand_sou_v_div_rates$comp_organ, 
            levels = unique(Brawand_sou_v_div_rates$comp_organ))

        # Combine DevSeq and Brawand GE divergence data
        compSouVDivRates <- rbind(DevSeq_sou_v_div_rates, Brawand_sou_v_div_rates)

    }




#---- Get GE divergence rates for Human vs species X (original data from Brawand 2011 paper)----


    # Use pearson correlation, inter-organ normalization and TPM for ms

    getBrBrainCor11 <- function(df, organ, coefficient) {

        df_cor <- sqrt(1 - cor(df, method=coefficient))
        df_cor <- df_cor[5:nrow(df_cor), 1:4]

        # Reshape cor data frame to one column
        df_cor_rs <- data.frame(newcol = c(t(df_cor)), stringsAsFactors=FALSE)

        Ptr <- mean(df_cor_rs[1:24,]) # chimp
        Ppa <- mean(df_cor_rs[25:36,]) # bonobo
        Ggo <- mean(df_cor_rs[37:44,]) # gorilla
        Ppy <- mean(df_cor_rs[45:52,]) # orangutan
        Mml <- mean(df_cor_rs[53:64,]) # macaque
        Mmu <- mean(df_cor_rs[65:76,]) # mouse
        Mdo <- mean(df_cor_rs[77:84,]) # opossum

        df_cor_avg <- rbind(Ppa, Ggo, Ppy, Mml, Mmu, Mdo)
        colnames(df_cor_avg) <- organ

        return(df_cor_avg)

    }

    brain_div_11 <- getBrBrainCor11(df=x_Br2011[,c(2:5,18:23,33:35,45:46,56:57,65:67,78:80,95:96)], 
        organ="Brain", coefficient=coefficient)


    getBrCerebCor11 <- function(df, organ, coefficient) {

        df_cor <- sqrt(1 - cor(df, method=coefficient))
        df_cor <- df_cor[3:nrow(df_cor), 1:2]

        # Reshape cor data frame to one column
        df_cor_rs <- data.frame(newcol = c(t(df_cor)), stringsAsFactors=FALSE)

        Ptr <- mean(df_cor_rs[1:4,]) # chimp
        Ppa <- mean(df_cor_rs[5:8,]) # bonobo
        Ggo <- mean(df_cor_rs[9:12,]) # gorilla
        Ppy <- mean(df_cor_rs[13:14,]) # orangutan
        Mml <- mean(df_cor_rs[15:18,]) # macaque
        Mmu <- mean(df_cor_rs[19:24,]) # mouse
        Mdo <- mean(df_cor_rs[25:28,]) # opossum

        df_cor_avg <- rbind(Ppa, Ggo, Ppy, Mml, Mmu, Mdo)
        colnames(df_cor_avg) <- organ

        return(df_cor_avg)

    }

    cereb_div_11 <- getBrCerebCor11(df=x_Br2011[,c(6:7,24:25,36:37,47:48,58,68:69,81:83,97:98)], 
        organ="Cerebellum", coefficient=coefficient)


    getBrHtKdCor11 <- function(df, organ, coefficient) {

        df_cor <- sqrt(1 - cor(df, method=coefficient))
        df_cor <- df_cor[4:nrow(df_cor), 1:3]

        # Reshape cor data frame to one column
        df_cor_rs <- data.frame(newcol = c(t(df_cor)), stringsAsFactors=FALSE)

        Ptr <- mean(df_cor_rs[1:6,]) # chimp
        Ppa <- mean(df_cor_rs[7:12,]) # bonobo
        Ggo <- mean(df_cor_rs[13:18,]) # gorilla
        Ppy <- mean(df_cor_rs[19:24,]) # orangutan
        Mml <- mean(df_cor_rs[25:30,]) # macaque
        Mmu <- mean(df_cor_rs[31:39,]) # mouse
        Mdo <- mean(df_cor_rs[40:45,]) # opossum

        df_cor_avg <- rbind(Ppa, Ggo, Ppy, Mml, Mmu, Mdo)
        colnames(df_cor_avg) <- organ

        return(df_cor_avg)

    }

    heart_div_11 <- getBrHtKdCor11(df=x_Br2011[,c(8:10,26:27,38:39,49:50,59:60,70:71,84:86,99:100)], 
        organ="Heart", coefficient=coefficient)
    kidney_div_11 <- getBrHtKdCor11(df=x_Br2011[,c(11:13,28:29,40:41,51:52,61:62,72:73,87:89,101:102)], 
        organ="Kidney", coefficient=coefficient)


    getBrLvCor11 <- function(df, organ, coefficient) {

        df_cor <- sqrt(1 - cor(df, method=coefficient))
        df_cor <- df_cor[3:nrow(df_cor), 1:2]

        # Reshape cor data frame to one column
        df_cor_rs <- data.frame(newcol = c(t(df_cor)), stringsAsFactors=FALSE)

        Ptr <- mean(df_cor_rs[1:4,]) # chimp
        Ppa <- mean(df_cor_rs[5:8,]) # bonobo
        Ggo <- mean(df_cor_rs[9:12,]) # gorilla
        Ppy <- mean(df_cor_rs[13:16,]) # orangutan
        Mml <- mean(df_cor_rs[17:20,]) # macaque
        Mmu <- mean(df_cor_rs[21:26,]) # mouse
        Mdo <- mean(df_cor_rs[27:30,]) # opossum

        df_cor_avg <- rbind(Ppa, Ggo, Ppy, Mml, Mmu, Mdo)
        colnames(df_cor_avg) <- organ

        return(df_cor_avg)

    }

    liver_div_11 <- getBrLvCor11(df=x_Br2011[,c(14:15,30:31,42:43,53:54,63:64,74:75,90:92,103:104)], 
        organ="Liver", coefficient=coefficient)


    getBrTestisCor11 <- function(df, organ, coefficient) {

        df_cor <- sqrt(1 - cor(df, method=coefficient))
        df_cor <- df_cor[3:nrow(df_cor), 1:2]

        # Reshape cor data frame to one column
        df_cor_rs <- data.frame(newcol = c(t(df_cor)), stringsAsFactors=FALSE)

        Ptr <- mean(df_cor_rs[1:2,]) # chimp
        Ppa <- mean(df_cor_rs[3:4,]) # bonobo
        Ggo <- mean(df_cor_rs[5:6,]) # gorilla
        Ppy <- NA # orangutan: no data available
        Mml <- mean(df_cor_rs[7:10,]) # macaque
        Mmu <- mean(df_cor_rs[11:14,]) # mouse
        Mdo <- mean(df_cor_rs[15:18,]) # opossum

        df_cor_avg <- rbind(Ppa, Ggo, Ppy, Mml, Mmu, Mdo)
        colnames(df_cor_avg) <- organ

        return(df_cor_avg)

    }

    testis_div_11 <- getBrTestisCor11(df=x_Br2011[,c(16:17,32,44,55,76:77,93:94,105:106)], 
        organ="Testis", coefficient=coefficient)

    Brawand11_organ_cor <- cbind(brain_div_11, cereb_div_11, heart_div_11, kidney_div_11, 
        liver_div_11, testis_div_11)


    # Reshape data table for ggplot
    # divergence times are estimated taxon pair times from TimeTree
    # http://www.timetree.org/
    div_times <- rep(c(6.7, 9.1, 15.8, 29.4, 90, 159), times=6)
    comp_organ <- rep(colnames(Brawand11_organ_cor), each=6)
    comp_spec <- rep(rownames(Brawand11_organ_cor), times=6)
    dataset <- rep("Mammals", 36)

    Brawand11_GE_div <- rbind(brain_div_11, cereb_div_11, heart_div_11, kidney_div_11, 
        liver_div_11, testis_div_11)
    rownames(Brawand11_GE_div) <- NULL
    colnames(Brawand11_GE_div) <- "correlation"

    Brawand11_div_rates <- data.frame(cbind(comp_spec, comp_organ, div_times, Brawand11_GE_div, dataset), 
        stringsAsFactors=FALSE)

    Brawand11_div_rates$div_times <- as.numeric(Brawand11_div_rates$div_times)
    Brawand11_div_rates$correlation <- as.numeric(Brawand11_div_rates$correlation)

    # Remove Orangutan testis (missing data) and Opossum kidney (replicate corr < 0.85) data
    Brawand11_div_rates <- Brawand11_div_rates[c(-33),]

    Brawand11_div_rates$comp_organ <- factor(Brawand11_div_rates$comp_organ, 
        levels = unique(Brawand11_div_rates$comp_organ))


    # Combine DevSeq and Brawand GE divergence data
    compDivRates11 <- rbind(DevSeq_div_rates, Brawand11_div_rates)


    if (is.element("TPM", expr_estimation)) {

        Brawand11_sou_v_div_rates <- data.frame(cbind(comp_spec, comp_organ, div_times, 
        Brawand2011_sou_v_div, dataset), stringsAsFactors=FALSE)

        Brawand11_sou_v_div_rates$div_times <- as.numeric(Brawand11_sou_v_div_rates$div_times)
        Brawand11_sou_v_div_rates$correlation <- as.numeric(Brawand11_sou_v_div_rates$correlation)

        # Remove ppy testis sample (has NA value)
        Brawand11_sou_v_div_rates <- Brawand11_sou_v_div_rates[c(-33),]

        Brawand11_sou_v_div_rates$comp_organ <- factor(Brawand11_sou_v_div_rates$comp_organ, 
            levels = unique(Brawand11_sou_v_div_rates$comp_organ))


        # Combine DevSeq and Brawand 2011 GE divergence data
        compSouVDivRates11 <- rbind(DevSeq_sou_v_div_rates, Brawand11_sou_v_div_rates)


        div_times_Br_all <- rep(c(6.7, 9.1, 15.8, 29.4, 90, 159, 177, 312), times=6)
        comp_organ_all <- rep(colnames(Brawand11_organ_cor), each=8)
        comp_spec_all <- rep(rownames(Br2011_all_sou_v[1:8,]), times=6)
        dataset_all <- rep("Mammals", 48)

        Brawand11_all_sou_v_div_rates <- data.frame(cbind(comp_spec_all, comp_organ_all, div_times_Br_all, 
            Brawand2011_all_sou_v_div, dataset_all), stringsAsFactors=FALSE)

        Brawand11_all_sou_v_div_rates$div_times_Br_all <- as.numeric(Brawand11_all_sou_v_div_rates$div_times_Br_all)
        Brawand11_all_sou_v_div_rates$correlation <- as.numeric(Brawand11_all_sou_v_div_rates$correlation)

        # Remove ppy testis sample (has NA value)
        Brawand11_all_sou_v_div_rates <- Brawand11_all_sou_v_div_rates[c(-43),]

        Brawand11_all_sou_v_div_rates$comp_organ_all <- factor(Brawand11_all_sou_v_div_rates$comp_organ_all, 
            levels = unique(Brawand11_all_sou_v_div_rates$comp_organ_all))

    }


    # Generate data set with both Brawand data (re-analyzed = "Mammals_DevSeq"; original = "Mammals_Brawand")
    Brawand_div_rates_comp <- Brawand_div_rates
    Brawand11_div_rates_comp <- Brawand11_div_rates
    Brawand_div_rates_comp$dataset[Brawand_div_rates_comp$dataset == 'Mammals'] <- 'Mammals(re-analyzed)'
    Brawand11_div_rates_comp$dataset[Brawand11_div_rates_comp$dataset == 'Mammals'] <- 'Mammals '
    compDivRatesBr <- rbind(Brawand_div_rates_comp, Brawand11_div_rates_comp)

    Brawand_sou_v_div_rates_comp <- Brawand_sou_v_div_rates
    Brawand11_sou_v_div_rates_comp <- Brawand11_sou_v_div_rates
    Brawand_sou_v_div_rates_comp$dataset <- 'Mammals(re-analyzed)'
    Brawand11_sou_v_div_rates_comp$dataset <- 'Mammals '
    compSouVDivRatesBr <- rbind(Brawand_sou_v_div_rates_comp, Brawand11_sou_v_div_rates_comp)




#---------------------------- Test the logarithmic regression model ----------------------------


    getModelFormula <- function(model) {

        formula <- as.formula(paste0("y ~ ", round(coefficients(model)[1], 2), "", 
            paste(sprintf(" %+.2f*%s ", 
                coefficients(model)[-1], 
                names(coefficients(model)[-1])), 
            collapse="")
            )
        )

        return(formula)
    }


    DevSeq_log_pea_lm <- lm(correlation ~ log(div_times), data = compDivRates[1:48,])
    Br_log_pea_lm <- lm(correlation ~ log(div_times), data = compDivRates11[49:nrow(compDivRates),])

    DevSeq_log_pea_lm_form <- getModelFormula(model = DevSeq_log_pea_lm)
    # y ~ 1.05 - 0.08 * log(div_times)
    Br_log_pea_lm_form <- getModelFormula(model = Br_log_pea_lm)
    # y ~ 0.99 - 0.04 * log(div_times)


    # log regression model for DevSeq pearson
    DevSeq_log_pea_lm_form_eq <- function(x) {
        pea_cor <- 1.04705 - 0.07686 * log(x)
        return(pea_cor)
    }
    # At 8.2e+05 Myr (8.2e+11 years), the DevSeq log regression model would reach a pearson cor of 0
    # This is 500 times longer than ATH-red algea divergence time

    # log regression model for DevSeq pearson
    Br_log_pea_lm_form_eq <- function(x) {
        pea_cor <-  0.99001 - 0.04469 * log(x)
        return(pea_cor)
    }
    # At 4.1e+09 Myr (4.1e+15 years), the Brawand log regression model would reach a pearson cor of 0


    # A log model will go to infinite -> check if pearson cor would reach 0 at emergence of plants
    # divergence times are estimated taxon pair times from TimeTree
    # http://www.timetree.org/ 
    # Angiosperms (Basal angiosperm, Amborella; 181)
    # Spermatophyta (seed plants: angiosperms, gymnosperms; 313)
    # Euphyllophyta (ferns and seed plants; 402), 
    # Tracheophyta (vascular plants; 431)
    # Bryophyta (non-vascular land plants: liverworts, hornworts and mosses; 471)
    # Charophyta (green algae, e.g. Coleochaetophyceae - one of the closest relatives to land plants; 778)
    # Chlorophyta (green algae: Chlorophyta and Charophyta/Streptophyta, e.g. volvox sp.; 1150)
    # Viridiplantae (green plants: 1150)
    # Rhodophyta (red algea; 1660)
    div_times_for_DevSeq <- c(7.1, 9.4, 25.6, 46, 106, 160, 181, 313, 402, 431, 471, 778, 1150, 1660)
    div_times_list_for_DevSeq <- list(div_times_for_DevSeq)
    species_for_DevSeq <- data.frame(Taxon=c("A.lyrata", "C.rubella", "E.salsugineum", "T.hassleriana", 
        "M.truncatula", "B.distachyon", "Angiosperms", "Spermatophyta", "Euphyllophyta", "Tracheophyta", 
        "Bryophyta", "Charophyta", "Chlorophyta", "Rhodophyta"))


    estimated_cor_for_DevSeq <- t(as.data.frame(do.call(cbind, lapply(div_times_for_DevSeq, 
        DevSeq_log_pea_lm_form_eq))))
    colnames(estimated_cor_for_DevSeq) <- "Estimated_cor"
    div_times_for_DevSeq <- data.frame(Divergence_time_Myr=div_times_for_DevSeq)
    DevSeq_log_pea_est <- cbind(div_times_for_DevSeq, species_for_DevSeq, estimated_cor_for_DevSeq)


    # Check Brawand sOU model based on GE data for HS to MDO
    # Compare model predictions with complete Brawand data
    Br_log_sOU_lm <- lm(correlation ~ log(div_times), data = Brawand11_sou_v_div_rates)
    Br_log_sOU_lm_form <- getModelFormula(model = Br_log_sOU_lm)
    # y ~ -0.15 + 0.14 * log(div_times)

    Br_all_log_sOU_lm <- lm(correlation ~ log(div_times_Br_all), data = Brawand11_all_sou_v_div_rates)
    Br_all_log_sOU_lm_form <- getModelFormula(model = Br_all_log_sOU_lm)
    # y ~ -0.2 + 0.16 * log(div_times_Br_all)


    # log regression model for DevSeq sOU_v
    Br_log_sOU_lm_form_eq <- function(x) {
        sOU <- -0.1456 + 0.1429 * log(x)
        return(sOU)
    }

    # log regression model for DevSeq sOU_v
    Br_all_log_sOU_lm_form_eq <- function(x) {
        sOU <-  -0.2014 + 0.1594 * log(x)
        return(sOU)
    }


# Predict fit of log model and confidence intervals
    predictBr.log.sOU <- function(x) {

        log_lm <- predict(Br_log_sOU_lm, data.frame(div_times = x), se.fit = TRUE, 
            interval = "confidence", level = 0.95)

        log_fit <- log_lm$fit
        log_fit <- as.data.frame(log_fit)
        x_data <- data.frame(div_time = x)
        model_fit <- cbind(x_data, log_fit)
        names(model_fit)[names(model_fit) == "fit"] <- "sOU_value"

        return(model_fit)
    }



    # Generate data based on log regression model that was built with Brawand11
    # sOU model with variable-µ distance and expression data from 7-150Myr
    # Use this model to predict transcriptome diostance at 350 Myr, and compare this
    # with real data from Brawand chicken (gga) ortholog genes

    x_Br_grid <- seq(6.7, 350, length = 250)  ## prediction grid
    x_Br_grid_list <- list(x_Br_grid)

    # Compute data points based on model
    Br_regr_pred <- as.data.frame(do.call(rbind, lapply(x_Br_grid, predictBr.log.sOU)))

    gga_ge_div <- rbind(Brawand11_all_sou_v_div_rates[8,], Brawand11_all_sou_v_div_rates[16,], 
        Brawand11_all_sou_v_div_rates[24,], Brawand11_all_sou_v_div_rates[32,], Brawand11_all_sou_v_div_rates[40,], 
        Brawand11_all_sou_v_div_rates[47,])

    Br_log_sOU_lm_form_txt <- paste(format(Br_log_sOU_lm_form))

    Br_log_sOU_lm_form_txt <- Br_log_sOU_lm_form_txt %<>% 
                                   gsub("~", "=", .) %>% 
                                   gsub("log", "ln", .) %>% 
                                   gsub("div_times", "x", .)



    # Make plot with predictive log regression
    plotDivPredict <- function(data, data2, regr_form) {

      fname <- sprintf('%s.jpg', paste("Br_prediction", sep="_"))

      p <- ggplot(data=data, aes(x=div_time, y=sOU_value)) + 
      geom_line(size = 3) + 
      geom_point(data=data2, aes(x=div_times_Br_all, y=correlation), size=5, colour="red") + 
      geom_ribbon(data=data, aes(ymin = lwr, ymax = upr), alpha=0.14) + 
      scale_x_continuous(limits = c(2.0, 350), expand = c(0.02,0)) + 
      scale_y_continuous(limits = c(0, 1.31), expand = c(0.02, 0)) + 
      geom_text(label = regr_form, x = 101, y = 1.193, color = "black", size=7.4) + 
      geom_rect(xmin = 12, xmax = 26, ymin = 1.1845, ymax = 1.2, color="black", fill="black", size=0.7) + 
      geom_rect(xmin = 0, xmax = 159, ymin = 0, ymax = 1.09, color="blue3", fill=NA, size=0.7) + 
      annotate("text", x = 312, y = 0.137, label= "Mouse", size=8, angle = 90) + 
      geom_segment(x = 312, xend = 312, y = -0.05, yend = 0.0025, color="black", size=0.7) + 
      guides(color = guide_legend(ncol = 3))

      q <- p + theme_bw() + xlab("Divergence time from HSA (Myr)") + ylab("Expression distance") + 
      theme(text=element_text(size=16), 
        axis.ticks.length=unit(0.35, "cm"), 
        axis.ticks = element_line(colour = "black", size = 0.95),  
        plot.margin = unit(c(0.55, 1.175, 0.5, 0.4),"cm"), 
        axis.title.y = element_text(size=24.5, margin = margin(t = 0, r = 15, b = 0, l = 11), colour="black", 
            face = "bold"), 
        axis.title.x = element_text(size=24.5, margin = margin(t = 15.75, r = 0, b = 1, l = 0), colour="black", 
            face = "bold"), 
        axis.text.x = element_text(size=21.5, angle=0, margin = margin(t = 5.5), colour="black"), 
        axis.text.y = element_text(size=21.75, angle=0, margin = margin(r = 5.5), colour="black"), 
        legend.box.background = element_rect(colour = "#d5d5d5", fill=NA, size=1.2), 
        panel.border = element_rect(colour = "black", fill=NA, size=1.8), 
        panel.grid.major = element_blank(),
        panel.grid.minor.x = element_blank(), 
        panel.grid.minor.y = element_blank(), 
        legend.position = c(0.757, 0.888), 
        legend.title = element_blank(), 
        legend.text = element_text(size=22.5), 
        legend.spacing.x = unit(0.5, 'cm'), 
        legend.key.size = unit(0.95, "cm"), 
        legend.background=element_blank())  

      ggsave(file = file.path(out_dir, "output", "plots", fname), plot = q, 
        width = 9.5, height = 6.75, dpi = 300, units = c("in"), limitsize = FALSE) 
    }

    plotDivPredict(data = Br_regr_pred, data2 = gga_ge_div, regr_form = Br_log_sOU_lm_form_txt)




#------ Compare slopes of regression models and perform analysis of covariance (ANCOVA) ------


    # Generate df that contains angiosperm, mammalian (Brawand 2011) and re-analyzed mammalian (Brawand-DevSeq) data
    compDivRates_sOU_all <- rbind(DevSeq_sou_v_div_rates, Brawand11_sou_v_div_rates, Brawand_sou_v_div_rates)

    model_sOU_lm <- lm(correlation ~ div_times * dataset, data = compDivRates_sOU_all) # model that assumes an interaction between the slopes of the regression lines of the two data sets and the grouping variable
    # this is same as
    model_sOU_lm = lm (correlation ~ div_times + dataset + div_times:dataset, data = compDivRates_sOU_all)
    model_sOU_lst_test <- lstrends(model_sOU_lm, "dataset", var="div_times") # get slopes
    sOU_pairs_test <- pairs(model_sOU_lst_test) # perform tukey test on all estimates

    sOU_pairs_test_df <- as.data.frame(summary(sOU_pairs_test))

    sOU_tukey_p_value <- paste("p =", formatC(sOU_pairs_test_df$p.value[1], format = "e",  digits = 1))
    sOU11_tukey_p_value <- paste("p =", formatC(sOU_pairs_test_df$p.value[2], format = "e",  digits = 1))
    sOUBr_tukey_p_value <- paste("p =", round(sOU_pairs_test_df$p.value[3], 2))



    # Compare angiosperm and mammalian divergence rates based on individual organ regressions
    # Prepare data for original and re-analyzed Brawand data
    organ_names_Br <- rep(c("Brain","Cerebellum","Heart","Kidney","Liver","Testis"), each=6)
    organ_names_Br11 <- rep(c("Brain.11","Cerebellum.11","Heart.11","Kidney.11","Liver.11","Testis.11"), each=6)
    organ_names_Br <- as.data.frame(organ_names_Br[-36])
    organ_names_Br11 <- as.data.frame(organ_names_Br11[-36])
    names(organ_names_Br) <- "comp_organ"
    names(organ_names_Br11) <- "comp_organ"
    organ_names <- rbind(organ_names_Br, organ_names_Br11)
    compSouVDivRatesBr_io <- compSouVDivRatesBr
    compSouVDivRatesBr_io[,2] <- organ_names
    compSouVDivRatesBr_io$comp_organ <- factor(compSouVDivRatesBr_io$comp_organ)


    # Retrieve slope value from individual organ regressions and compute p value
    getTrendsP <- function(corrdata, reg_model = c("lm_reg", "poly2_reg", "log_reg", "sqrt_reg")) {

        if ((missing(reg_model)) || (!is.element(reg_model, c("lm_reg", "poly2_reg", "log_reg", "sqrt_reg")))) 

            stop(
                "Regression model either missing or not one of: 
                'lm_reg', 'poly2_reg', 'log_reg', 'sqrt_reg'",
                call. = TRUE
                )

        if (reg_model == "lm_reg") {

            model_var <- correlation ~ div_times * comp_organ #linear regression

        } else if (reg_model == "poly2_reg") {

            model_var <- correlation ~ poly(div_times, 2, raw = TRUE) * comp_organ #polynomial regression w/quadratic term

        } else if (reg_model == "log_reg") {

            model_var <- correlation ~ log(div_times) * comp_organ #linear regression with log-x-transform

        } else if (reg_model == "sqrt_reg") {

            model_var <- correlation ~ sqrt(div_times) * comp_organ #linear regression with sqrt-x-transform
        }

        model_lmp = lm(model_var, data = corrdata)

        div_trend <- lstrends(model_lmp, "comp_organ", var = "div_times") # get slopes for regressions
        div_trend_df <- summary(div_trend)

        if (deparse(substitute(corrdata)) == "compDivRates" | deparse(substitute(corrdata)) == "compDivRates11"
            | deparse(substitute(corrdata)) == "compSouVDivRates" | deparse(substitute(corrdata)) == "compSouVDivRates11") {

            trend_stat_welch <- t.test(div_trend_df[1:8,2], div_trend_df[9:14,2], paired = FALSE)
            trend_stat_wilcox <- wilcox.test(div_trend_df[1:8,2], div_trend_df[9:14,2], paired = FALSE)

        } else if (deparse(substitute(corrdata)) == "compSouVDivRatesBr_io") {

            trend_stat_welch <- t.test(div_trend_df[1:6,2], div_trend_df[7:12,2], paired = FALSE)
            trend_stat_wilcox <- wilcox.test(div_trend_df[1:6,2], div_trend_df[7:12,2], paired = FALSE)
        }

        welch_p_value <- as.data.frame(trend_stat_welch$p.value)
        wilcox_p_value <- as.data.frame(trend_stat_wilcox$p.value)

        colnames(welch_p_value) <- reg_model
        rownames(welch_p_value) <- "Welch_test"

        colnames(wilcox_p_value) <- reg_model
        rownames(wilcox_p_value) <- "Wilcox_test"

        slopes <- as.data.frame(div_trend_df[,2])
        colnames(slopes) <- reg_model
        rownames(slopes) <- div_trend_df[,1]

        trends_p_value <- rbind(slopes, welch_p_value, wilcox_p_value)

        return(trends_p_value)
    }


    # Create list of models
    model_list <- list("lm_reg", "poly2_reg", "log_reg", "sqrt_reg")

    # Slopes and p-values for DevSeq angiosperm vs. re-analyzed Brawand mammalian data
    p_values_compDivRates_io <- as.data.frame(do.call(cbind, lapply(model_list, getTrendsP, corrdata = compDivRates)))
    # Slopes and p-values for DevSeq angiosperm vs. original 2011 Brawand mammalian data
    p_values_compDivRates11_io <- as.data.frame(do.call(cbind, lapply(model_list, getTrendsP, corrdata = compDivRates11)))
    # Slopes and p-values for DevSeq angiosperm vs. re-analyzed Brawand mammalian data
    p_values_compSouVDivRates_io <- as.data.frame(do.call(cbind, lapply(model_list, getTrendsP, corrdata = compSouVDivRates)))
    # Slopes and p-values for DevSeq angiosperm vs. original 2011 Brawand mammalian data
    p_values_compSouVDivRates11_io <- as.data.frame(do.call(cbind, lapply(model_list, getTrendsP, corrdata = compSouVDivRates11)))
    

    # Get text string of p-values (rank-sum test) for divergence plots with log models
    p_value_io <- paste("p =", round(p_values_compDivRates_io["Wilcox_test", "log_reg"], 4))
    p_value11_io <- paste("p =", round(p_values_compDivRates11_io["Wilcox_test", "log_reg"], 4))
    p_value_SouVio <- paste("p =", round(p_values_compSouVDivRates_io["Wilcox_test", "log_reg"], 4))
    p_value11_SouVio <- paste("p =", round(p_values_compSouVDivRates11_io["Wilcox_test", "log_reg"], 3))




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



    getLOESS.Slopes <- function(organ_data, data_set) {

      comp_organ <- unique(organ_data$comp_organ)

      if ((data_set == "selected") && (comp_organ == "Root")) {

        poly_deg <- 1
        alpha <- 1

      } else if ((data_set == "selected") && (comp_organ == "Hypocotyl")) {

        poly_deg <- 1
        alpha <- 1

      } else if ((data_set == "complete") && (comp_organ == "Root")) {

        poly_deg <- 1
        alpha <- 1

      } else {

        poly_deg <- 2
        alpha <- 1

      } # Use quadratic poynomes for all organs except root

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
    DevSeqSouV_AT_loess_slopes <- as.data.frame(do.call(rbind, lapply(devseqSouV_organ_lst, getLOESS.Slopes, data_set="complete")))
    DevSeqSouV_sel_AT_loess_slopes <- as.data.frame(do.call(rbind, lapply(devseqSouV_organ_lst_sel, getLOESS.Slopes, data_set="selected"))) ## hypocotyl slope is 0.00779 instead 0.00758 if BD is left out
    Brawand11SouV_loess_slopes <- as.data.frame(do.call(rbind, lapply(brawandSouV11_organ_lst, getLOESS.Slopes, data_set="complete")))
    sOU_loess_DevSeq_AT_Br11_wilcox <- wilcox.test(as.numeric(unlist(DevSeqSouV_AT_loess_slopes)), as.numeric(unlist(Brawand11SouV_loess_slopes)))$p.value

    BrawandSouV_loess_slopes <- as.data.frame(do.call(rbind, lapply(brawandSouV_organ_lst, getLOESS.Slopes, data_set="complete")))
    sOU_loess_Br_Br11_wilcox <- wilcox.test(as.numeric(unlist(Brawand11SouV_loess_slopes)), as.numeric(unlist(BrawandSouV_loess_slopes)))$p.value

    # For metric pearson expression distances
    DevSeq_AT_loess_slopes <- as.data.frame(do.call(rbind, lapply(devseq_organ_lst, getLOESS.Slopes, data_set="complete")))
    Brawand11_loess_slopes <- as.data.frame(do.call(rbind, lapply(brawand11_organ_lst, getLOESS.Slopes, data_set="complete")))

    Brawand_loess_slopes <- as.data.frame(do.call(rbind, lapply(brawand_organ_lst, getLOESS.Slopes, data_set="complete")))


    # Write slope values to csv file
    DevSeq_slopes <- cbind(DevSeqSouV_AT_loess_slopes, DevSeqSouV_AT_loess_slopes, DevSeq_AT_loess_slopes, DevSeq_AT_loess_slopes)
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


    # Get table with slope values and test statistics for log regression
    DS_AT_Br_log_slopes <- cbind(p_values_compSouVDivRates11_io$log_reg, p_values_compSouVDivRates_io$log_reg, 
      p_values_compDivRates11_io$log_reg, p_values_compDivRates_io$log_reg)
    DS_AT_Br_log_slopes <- as.data.frame(DS_AT_Br_log_slopes)
    DS_AT_Br_log_slopes <- cbind(DS_AT_Br_loess_slopes$sample, DS_AT_Br_log_slopes)
    colnames(DS_AT_Br_log_slopes) <- c("sample", "DS_AT_Br11_sOU_log", "DS_AT_Br_sOU_log", "DS_AT_Br11_pea_log", "DS_AT_Br_pea_log")


    # Show message
    message("Writing data tables...")

    write.table(DS_AT_Br_loess_slopes, 
      file=file.path(out_dir, "output", "data", "DS_AT_Br_loess_slopes.txt"), sep="\t", 
      col.names=TRUE, row.names=FALSE, dec=".", quote = FALSE)

    write.table(DS_AT_Br_log_slopes, 
      file=file.path(out_dir, "output", "data", "DS_AT_Br_log_slopes.txt"), sep="\t", 
      col.names=TRUE, row.names=FALSE, dec=".", quote = FALSE)


    # Create p-value containing test strings for plots
    p_dat_text <- data.frame(
        label = c(paste("italic('P =')~", formatC(sOU_loess_DevSeq_AT_Br11_wilcox, format = "e", digits = 0))),
        x = c(16.25),
        y = c(1.3328)
    )



   # Make sOU GE divergence plot showing individual organ regressions for SI
   makeOrgRegPlot <- function(data1, data2, coefficient, expr_estimation, p_value, pos) {

      fname <- sprintf('%s.jpg', paste("compSouVDivRates11_loess", expr_estimation, pos, sep="_"))

      ymin <- 0.05
      ymax <- 1.45

      if (pos == "main") {

            plot_wdt <- 12.535
            plot_hdt <- 8
            legend_x_pos <- 0.241
            legend_y_pos <- 0.914
            linewd <- 3
            point_size <- 5.75

      } else {

            plot_wdt <- 9.5 # condenced plot width for suppl
            plot_hdt <- 6.75 # condenced plot width for suppl
            legend_x_pos <- 0.317
            legend_y_pos <- 0.9
            linewd <- 2.5
            point_size <- 5
      }

      ds_col <- rep(c("#798dc4"), 48)
      bw_col <- rep(c("red"), 35)
      col_scale <- c("#798dc4", "red")
      fill_scale <- c("#798dc4", "red")
      col_breaks <- c("Angiosperms ", "Mammals")
      fill_breaks <- c("Angiosperms ", "Mammals")
      shape_scale <- c(16, 17)

      fill_col <- c(as.character(ds_col), as.character(bw_col))

      p <- ggplot() 
      p <- p + geom_line(size = linewd, data = data1, aes(x = div_times, y = correlation, group = comp_organ, 
        colour = dataset, linetype = dataset))
      p <- p + geom_point(size = point_size, data = data2, aes(x = div_times, y = correlation, group = comp_organ, 
        colour = dataset, shape = dataset))
      p <- p + scale_x_continuous(limits = c(0,161), expand = c(0.02,0), breaks = c(0,20,40,60,80,100,120,140,160)) + 
      scale_y_continuous(limits = c(0.055, 1.588), expand = c(0.02, 0), breaks = c(0.2,0.4,0.6,0.8,1,1.2,1.4)) + 
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
        legend.background=element_blank()) 

        ggsave(file = file.path(out_dir, "output", "plots", fname), plot = q, 
            width = plot_wdt, height = plot_hdt, dpi = 300, units = c("in"), limitsize = FALSE) 
  }

  makeOrgRegPlot(data1 = loessSouV_coor11_AT, data2 = compSouVDivRates11, coefficient = coefficient, 
    expr_estimation = expr_estimation, pos = "ext", p_value = "")

  makeOrgRegPlot(data1 = loessSouV_coor11_AT, data2 = compSouVDivRates11, coefficient = coefficient, 
    expr_estimation = expr_estimation, pos = "main", p_value = sOU_loess_DevSeq_AT_Br11_slope_p)




   # Make sOU plot of LOESS regression from Brawand11 and Brawand_re-analyzed data for SI
   makeGlobalBrRegPlot <- function(data, coefficient, expr_estimation, p_value) {

      if (deparse(substitute(data)) == "compSouVDivRatesBr") {

        fname <- sprintf('%s.jpg', paste("Brawand_vs_Brawand_2011_sOU_loess", expr_estimation, sep="_"))

        ymin <- 0.05
        ymax <- 1.45

        ds_col <- rep(c("red3"), 35)
        bw_col <- rep(c("red"), 35)
        col_scale <- c("red", "red3")
        fill_scale <- c("red", "red3")
        col_breaks <- c("Mammals ", "Mammals(re-analyzed)")
        fill_breaks <- c("Mammals ", "Mammals(re-analyzed)")
        shape_scale <- c(17, 15)
        linetype_a <- rep("dashed",35)
        linetype_b <- rep("solid",35)
        leg_pos_x <- 0.4051
        leg_pos_y <- 0.9

      } else if (deparse(substitute(data)) == "compSouVDivRates11") {

        fname <- sprintf('%s.jpg', paste("compSouVDivRates11_loess_avg", expr_estimation, sep="_"))

        ymin <- 0.05
        ymax <- 1.45

        ds_col <- rep(c("#798dc4"), 48)
        bw_col <- rep(c("red"), 35)
        col_scale <- c("#798dc4", "red")
        fill_scale <- c("#798dc4", "red")
        col_breaks <- c("Angiosperms ", "Mammals")
        fill_breaks <- c("Angiosperms ", "Mammals")
        shape_scale <- c(16, 17)
        linetype_a <- rep("dashed",48)
        linetype_b <- rep("solid",35)
        leg_pos_x <- 0.317
        leg_pos_y <- 0.9
      }

      fill_col <- c(as.character(ds_col), as.character(bw_col))

      p <- ggplot() 
      p <- p + geom_smooth(size = 2.5, data = data, span=1, alpha = 0.17, aes(x = div_times, y = correlation, 
        group = dataset, colour = dataset, linetype = c(linetype_a, linetype_b), fill = fill_col))
      p <- p + geom_point(size = 5, data = data, aes(x = div_times, y = correlation, group = dataset, 
        colour = dataset, shape = dataset))
      p <- p + scale_x_continuous(limits = c(0,160), expand = c(0.02,0), breaks = c(0,20,40,60,80,100,120,140,160)) + 
      scale_y_continuous(limits = c(ymin, ymax), expand = c(0.02, 0)) + 
      scale_color_manual(values = col_scale, breaks = col_breaks) + 
      scale_fill_manual(values = fill_scale, breaks = fill_breaks) + 
      scale_shape_manual(values = shape_scale) + 
      scale_size(range = c(0.5, 12)) + 
      guides(color = guide_legend(ncol=2, keywidth = 0.4, keyheight = 0.4, default.unit = "inch"), 
        linetype="none")

      q <- p + theme_bw() + xlab("Divergence time (Myr)") + ylab("Expression distance") + 
      theme(text=element_text(size=16), 
        axis.ticks.length=unit(0.35, "cm"), 
        axis.ticks = element_line(colour = "black", size = 0.7),  
        plot.margin = unit(c(0.55, 1.175, 0.5, 0.4),"cm"), 
        axis.title.y = element_text(size=25, margin = margin(t = 0, r = 15, b = 0, l = 11), colour="black"), 
        axis.title.x = element_text(size=25, margin = margin(t = 14.75, r = 0, b = 2, l = 0), colour="black"), 
        axis.text.x = element_text(size=21.25, angle=0, margin = margin(t = 5.5), colour="black"), 
        axis.text.y = element_text(size=21.25, angle=0, margin = margin(r = 5.5), colour="black"), 
        legend.box.background = element_rect(colour = "#d5d5d5", fill=NA, size=1.0), 
        panel.border = element_rect(colour = "black", fill=NA, size=1.75), 
        panel.grid.major = element_line(color="#d5d5d5"),
        panel.grid.minor.x = element_blank(), 
        panel.grid.minor.y = element_blank(), 
        legend.position = c(leg_pos_x, leg_pos_y), 
        legend.title = element_blank(), 
        legend.text = element_text(size=22), 
        legend.spacing.x = unit(0.5, 'cm'), 
        legend.key.size = unit(0.95, "cm"), 
        legend.background=element_blank()) 

        ggsave(file = file.path(out_dir, "output", "plots", fname), plot = q, 
            width = 9.5, height = 6.75, dpi = 300, units = c("in"), limitsize = FALSE) 
      
    }

    makeGlobalBrRegPlot(data = compSouVDivRates11, coefficient = coefficient, 
      expr_estimation = expr_estimation)
    makeGlobalBrRegPlot(data = compSouVDivRatesBr, coefficient = coefficient, 
      expr_estimation = expr_estimation)

  }




#---- Make gene expression divergence rates plot for Brawand data (original and re-analyzed) ---


   makeGEDivPlotBr <- function(data) {

      fname <- sprintf('%s.jpg', paste(deparse(substitute(data)), coefficient, expr_estimation, sep="_"))

      if (deparse(substitute(data)) == "Brawand_div_rates") {
        y_min <- 0.34
        y_max <- 0.692
        yseg_min <- 0.24
        yseg_max <- 0.3418
        y_text <- 0.3525

      } else if (deparse(substitute(data)) == "Brawand11_div_rates") {
        y_min <- 0.20
        y_max <- 0.71
        yseg_min <- 0.15
        yseg_max <- 0.2027
        y_text <- 0.218
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


getATDiv(expr_estimation = "TPM", coefficient = "pearson")



