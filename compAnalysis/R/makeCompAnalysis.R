# Prepare Brawand and DevSeq comparative ortholog gene expression data
# Thresholds: DevSeq 0.05 ERCC; Brawand 0.5 TPM (no ERCC spike-ins available)
# Data input: Brawand and DevSeq TPM expression tables of all samples



#-------------------------------------- Read data tables ---------------------------------------


makeCompAnylsis <- function(dataset = c("Brawand", "DevSeq"), expr_estimation = c("TPM", "counts"), 
	coefficient = c("pearson", "spearman"), devseq_spec = c("Brassicaceae", "all"), 
    data_norm = c("intra-organ", "inter-organ"), devseq_organs = c("all", "subset")) {
	
	# Show error message if no species or unknown data set is chosen
    if ((missing(dataset)) || (!is.element(dataset, c("Brawand", "DevSeq"))))
   
       stop(
       "Please choose one of the available data sets: 
	   'Brawand', 'DevSeq'",
	   call. = TRUE
       )

   	# Show error message if expression estimation or unknown expression estimation is chosen
    if ((missing(expr_estimation)) || (!is.element(expr_estimation, c("TPM", "counts"))))
   
       stop(
       "Please choose one of the available expression estimations: 
	   'TPM', 'counts'",
	   call. = TRUE
       )

    # Show error message if no correlation or unknown correlation coefficient is chosen
    if ((missing(coefficient)) || (!is.element(coefficient, c("pearson", "spearman"))))
   
       stop(
       "Please choose one of the available correlation coefficients: 
	   'pearson', 'spearman'",
	   call. = TRUE
       )

    # Show error message if no devseq_spec or unknown devseq_spec is chosen
    if ((is.element("DevSeq", dataset)) && ((missing(devseq_spec)) || (!is.element(devseq_spec, c("Brassicaceae", "all")))))
   
       stop(
       "Please choose one of the available DevSeq species sets: 
	   'Brassicaceae', 'all'",
	   call. = TRUE
       )

    # Show error message if no devseq_spec or unknown devseq_spec is chosen
    if ((is.element("DevSeq", dataset) && is.element("all", devseq_spec)) && ((missing(devseq_organs)) || (!is.element(devseq_organs, c("all", "subset")))))
   
       stop(
       "Please choose one of the available DevSeq organ sets: 
       'all', 'subset'",
       call. = TRUE
       )

    # Show error message if expression estimation or unknown expression estimation is chosen
    if ((missing(data_norm)) || (!is.element(data_norm, c("intra-organ", "inter-organ"))))
   
       stop(
       "Please choose one of the available data_norm data normalizations: 
       'intra-organ', 'inter-organ'",
       call. = TRUE
       )


    # Show startup message
    message("Reading data...")


	# Set expression input file
    if ((is.element("Brawand", dataset)) && (is.element("TPM", expr_estimation)) 
        && (is.element("intra-organ", data_norm))) {
        genesExpr = file.path(in_dir, "Expression_data", "Brawand_intra_tpm_mat_deseq_sample_names_0_5_threshold.csv")

    } else if ((is.element("Brawand", dataset)) && (is.element("TPM", expr_estimation)) 
        && (is.element("inter-organ", data_norm))) {
        genesExpr = file.path(in_dir, "Expression_data", "Brawand_inter_tpm_mat_deseq_sample_names_0_5_threshold.csv")

    } else if ((is.element("Brawand", dataset)) && (is.element("counts", expr_estimation))
        && (is.element("intra-organ", data_norm))) {
        genesExpr = file.path(in_dir, "Expression_data", "Brawand_intra_count_mat_vsd_sample_names_0_5_threshold.csv")

    } else if ((is.element("Brawand", dataset)) && (is.element("counts", expr_estimation))
        && (is.element("inter-organ", data_norm))) {
        genesExpr = file.path(in_dir, "Expression_data", "Brawand_inter_count_mat_vsd_sample_names_0_5_threshold.csv")

    } else if ((is.element("DevSeq", dataset)) && (is.element("TPM", expr_estimation)) 
        && (is.element("Brassicaceae", devseq_spec)) && (is.element("intra-organ", data_norm))) {
        genesExpr = file.path(in_dir, "Expression_data", "AT_brass_intra_tpm_mat_deseq_sample_names.csv")

    } else if ((is.element("DevSeq", dataset)) && (is.element("TPM", expr_estimation)) 
        && (is.element("Brassicaceae", devseq_spec)) && (is.element("inter-organ", data_norm))) {
        genesExpr = file.path(in_dir, "Expression_data", "AT_brass_inter_tpm_mat_deseq_sample_names.csv")

    } else if ((is.element("DevSeq", dataset)) && (is.element("TPM", expr_estimation)) 
        && (is.element("all", devseq_spec)) && (is.element("intra-organ", data_norm))) {
        genesExpr = file.path(in_dir, "Expression_data", "AT_core_intra_tpm_mat_deseq_sample_names.csv")

    } else if ((is.element("DevSeq", dataset)) && (is.element("TPM", expr_estimation)) 
        && (is.element("all", devseq_spec)) && (is.element("inter-organ", data_norm))) {
        genesExpr = file.path(in_dir, "Expression_data", "AT_core_inter_tpm_mat_deseq_sample_names.csv")

    } else if ((is.element("DevSeq", dataset)) && (is.element("counts", expr_estimation)) 
        && (is.element("Brassicaceae", devseq_spec)) && (is.element("intra-organ", data_norm))) {
        genesExpr = file.path(in_dir, "Expression_data", "AT_brass_intra_count_mat_vsd_sample_names.csv")

    } else if ((is.element("DevSeq", dataset)) && (is.element("counts", expr_estimation)) 
        && (is.element("Brassicaceae", devseq_spec)) && (is.element("inter-organ", data_norm))) {
        genesExpr = file.path(in_dir, "Expression_data", "AT_brass_inter_count_mat_vsd_sample_names.csv")

    } else if ((is.element("DevSeq", dataset)) && (is.element("counts", expr_estimation)) 
        && (is.element("all", devseq_spec)) && (is.element("intra-organ", data_norm))) {
        genesExpr = file.path(in_dir, "Expression_data", "AT_core_intra_count_mat_vsd_sample_names.csv")

    } else if ((is.element("DevSeq", dataset)) && (is.element("counts", expr_estimation)) 
        && (is.element("all", devseq_spec)) && (is.element("inter-organ", data_norm))) {
        genesExpr = file.path(in_dir, "Expression_data", "AT_core_inter_count_mat_vsd_sample_names.csv")
    }


    # Get data set is
    if (is.element("Brawand", dataset)) {
        dataset_id <- "Brawand"

    } else if (is.element("DevSeq", dataset)) {
        dataset_id <- "DevSeq"
    }


    # Get data normalization method
    if (is.element("intra-organ", data_norm)) {
        data_norm <- "intra-organ"

    } else if (is.element("inter-organ", data_norm)) {
        data_norm <- "inter-organ"
    }


	# Define simplified Brawand and DevSeq column names
	if (is.element("DevSeq", dataset) && (is.element("all", devseq_spec))) {
        col_names <- rep(c("Root", "Hypocotyl", "Leaf", "veg_apex", "inf_apex", 
            "Flower", "Stamen", "Carpel", "Pollen"), each=21)
        replicate_tag_samples <- rep(c(".1",".2",".3"), times=9)
        col_names <- paste0(col_names,replicate_tag_samples)
        spec_names <- rep(c("_AT", "_AL", "_CR", "_ES", "_TH", "_MT", "_BD"), each=3)
        spec_names <- rep(spec_names, times=9)
        col_names <- paste0(col_names, spec_names)
        col_names <- c("gene_id", col_names)

    } else if (is.element("DevSeq", dataset) && (is.element("Brassicaceae", devseq_spec))) {
        col_names <- rep(c("Root", "Hypocotyl", "Leaf", "veg_apex", "inf_apex", 
            "Flower", "Stamen", "Carpel", "Pollen"), each=12)
        replicate_tag_samples <- rep(c(".1",".2",".3"), times=4)
        col_names <- paste0(col_names,replicate_tag_samples)
        spec_names <- rep(c("_AT", "_AL", "_CR", "_ES"), each=3, times=9)
        col_names <- paste0(col_names, spec_names)
        col_names <- c("gene_id", col_names) 
    }


	# Read expression data
	if (is.element("DevSeq", dataset)) {
		x <- read.table(genesExpr, sep=";", dec=".", header=TRUE, stringsAsFactors=FALSE)

	} else if (is.element("Brawand", dataset)) {
		x <- read.table(genesExpr, sep=";", dec=".", header=TRUE, stringsAsFactors=FALSE)
	}


    # Stop function here to allow specific analysis of a single data set
    # For DevSeq
    # return_list <- list("dataset_id" = dataset_id, "expr_estimation" = expr_estimation, "x" = x, "coefficient" = coefficient, "col_names" = col_names, "devseq_spec" = devseq_spec, "data_norm" = data_norm, "devseq_organs" = devseq_organs)
    # For Brawand
    # return_list <- list("dataset_id" = dataset_id, "expr_estimation" = expr_estimation, "x" = x, "coefficient" = coefficient, "data_norm" = data_norm)
    # return(return_list)
    # }
    # return_objects <- makeCompAnylsis(dataset="DevSeq", expr_estimation="counts", coefficient="pearson", devseq_spec="all", data_norm="inter-organ", devseq_organs="all") # read in DevSeq expression data
    # return_objects <- makeCompAnylsis(dataset="Brawand", expr_estimation="counts", coefficient="pearson", data_norm="inter-organ") # read in Brawand expression data
    # list2env(return_objects, envir = .GlobalEnv)

    # Update column names
    if (dataset_id == "DevSeq") {

    	# set column names
    	colnames(x) <- col_names
    
    } else if (dataset_id == "Brawand") {

        # Generate a sequence to replace missing gene_id column
        # Remove this in case input table will have gene_id column again
        ID_repl <- as.data.frame(seq(1:nrow(x)))
        colnames(ID_repl) <- "gene_id"
        x <- cbind(ID_repl, x)

        # Merge replicates
        # Need to do this manually because different number of replicates across organs and species
        x <- data.frame(cbind("gene_id"=x[,1], rowMeans(x[,2:5]), rowMeans(x[,6:8]), 
            rowMeans(x[,9:14]), rowMeans(x[,15:16]), rowMeans(x[,17:18]), rowMeans(x[,19:21]), 
            rowMeans(x[,22:24]), rowMeans(x[,25:26]), rowMeans(x[,27:28]), rowMeans(x[,29:30]), 
            rowMeans(x[,31:32]), rowMeans(x[,33:34]), x[,35], rowMeans(x[,36:37]), 
            rowMeans(x[,38:40]), rowMeans(x[,41:42]), rowMeans(x[,43:44]), rowMeans(x[,45:46]), 
            rowMeans(x[,47:48]), rowMeans(x[,49:50]), rowMeans(x[,51:52]), rowMeans(x[,53:54]), 
            rowMeans(x[,55:57]), rowMeans(x[,58:59]), rowMeans(x[,60:61]), rowMeans(x[,62:63]), 
            rowMeans(x[,64:65]), rowMeans(x[,66:67]), rowMeans(x[,68:69]), rowMeans(x[,70:71]), 
            rowMeans(x[,72:74]), x[,75], rowMeans(x[,76:77]), rowMeans(x[,78:79]), 
            rowMeans(x[,80:81]), rowMeans(x[,82:83]), rowMeans(x[,84:85]), rowMeans(x[,86:87]), 
            rowMeans(x[,88:90]), rowMeans(x[,91:92]), rowMeans(x[,93:94]), x[,95:97],  
            rowMeans(x[,98:99]), rowMeans(x[,100:101]), rowMeans(x[,102:103])))

        tetra_organs <- rep(c("br", "cb", "ht", "kd", "lv", "ts"), each=8)
        tetra_species <- rep(c("Hsa", "Ppa", "Ptr", "Ggo", "Ppy", "Mml", "Mmu", "Mdo"), times=6)
        tetra_colnames <- paste(tetra_organs, tetra_species, sep="_")
        tetra_colnames <- tetra_colnames[-45]

        colnames(x)[2:ncol(x)] <- tetra_colnames

    }




#--------------------- Prepare data and define color palette for corrplot ---------------------


    # Create "plots" folder in /out_dir/output/plots
    if (!dir.exists(file.path(out_dir, "output", "plots"))) 
        dir.create(file.path(out_dir, "output", "plots"), recursive = TRUE)

    # Show message
    message("Starting analysis and generate plots...")


    # This is a wrapper function for the colorRampPalette. It allows to define intermediate colors
    # between the main colors, and enables to stretch out colors that should predominate the spectrum
    # color.palette is based on Marc's post on stackoverflow (https://stackoverflow.com/a/13327326)

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
    steps <- c("#d93527", "#d93527", "#faa92e", "#fff421", "#fefef4")

    if (dataset_id == "DevSeq") {
        pal <- color.palette(steps, c(30, 29, 26, 10), space = "rgb")
    } else if (dataset_id == "Brawand") {
        pal <- color.palette(steps, c(25, 29, 26, 15), space = "rgb")
    }

    # Set filename
    fname <- sprintf('%s.png', paste(dataset_id, expr_estimation, coefficient, sep="_"))

    if (devseq_organs == "subset") {

        fname <- sprintf('%s.png', paste(dataset_id, expr_estimation, coefficient, devseq_organs, sep="_"))
    }


    # Define column and row colors for color bars based on sample names and experiment
    # Build distance matrix & dendrogram then get dendrogram leaf colors to create color vector
    if (dataset_id == "DevSeq") {

        exp_col <- c(AT="gray1", AL="cornsilk3", CR="#fff7e8", ES="#95928c", TH="#bd0000", 
            MT="lightgoldenrod2", BD="#a63ea2") # last two letters of sample name

        species_col <- c(Roo="#6a54a9", Hyp="#53b0db", Lea="#2c8654", veg="#96ba37", 
            inf="#fad819", Flo="#e075af", Sta="#ed311c", Mat="#a63126", Car="#f2a72f")

    } else if (dataset_id == "Brawand") {

        exp_col <- c(sa="gray1", tr="cornsilk3", pa="#fff7e8", go="#998c74", py="lightgoldenrod2", 
            ml="#bd0000", mu="#8f428c", do="green") # last two letters of sample name

        species_col <- c(br_="#42448c", cb_="#6a76ad", ht_="gold1", kd_="olivedrab3", lv_="palegreen4", 
            ts_="#e17a21")
    }


    x[is.na(x)] <- 0 # replaces NAs by 0

    # Remove ERCC spike-ins from data
    x <- x[!grepl("ERCC", x$gene_id),]

    x_df <- x

    # Remove pollen samples for hclust heatmap
    if ((dataset_id == "DevSeq") && (devseq_spec == "all") && (devseq_organs == "all")) {
        x_df <- x_df %>% select (-c(Pollen.1_AT, Pollen.2_AT, Pollen.3_AT, Pollen.1_AL, Pollen.2_AL, 
            Pollen.3_AL, Pollen.1_CR, Pollen.2_CR, Pollen.3_CR, Pollen.1_ES, Pollen.2_ES, Pollen.3_ES, 
            Pollen.1_TH, Pollen.2_TH, Pollen.3_TH, Pollen.1_MT, Pollen.2_MT, Pollen.3_MT, Pollen.1_BD, 
            Pollen.2_BD, Pollen.3_BD))

    } else if ((dataset_id == "DevSeq") && (devseq_spec == "all") && (devseq_organs == "subset")) {

        x_df <- x_df %>% select (-c(Pollen.1_AT, Pollen.2_AT, Pollen.3_AT, Pollen.1_AL, Pollen.2_AL, 
            Pollen.3_AL, Pollen.1_CR, Pollen.2_CR, Pollen.3_CR, Pollen.1_ES, Pollen.2_ES, Pollen.3_ES, 
            Pollen.1_TH, Pollen.2_TH, Pollen.3_TH, Pollen.1_MT, Pollen.2_MT, Pollen.3_MT, Pollen.1_BD, 
            Pollen.2_BD, Pollen.3_BD, 
            veg_apex.1_AT, veg_apex.2_AT, veg_apex.3_AT, veg_apex.1_AL, veg_apex.2_AL, 
            veg_apex.3_AL, veg_apex.1_CR, veg_apex.2_CR, veg_apex.3_CR, veg_apex.1_ES, veg_apex.2_ES, veg_apex.3_ES, 
            veg_apex.1_TH, veg_apex.2_TH, veg_apex.3_TH, veg_apex.1_MT, veg_apex.2_MT, veg_apex.3_MT, veg_apex.1_BD, 
            veg_apex.2_BD, veg_apex.3_BD, 
            inf_apex.1_AT, inf_apex.2_AT, inf_apex.3_AT, inf_apex.1_AL, inf_apex.2_AL, 
            inf_apex.3_AL, inf_apex.1_CR, inf_apex.2_CR, inf_apex.3_CR, inf_apex.1_ES, inf_apex.2_ES, inf_apex.3_ES, 
            inf_apex.1_TH, inf_apex.2_TH, inf_apex.3_TH, inf_apex.1_MT, inf_apex.2_MT, inf_apex.3_MT, inf_apex.1_BD, 
            inf_apex.2_BD, inf_apex.3_BD, 
            Flower.1_AT, Flower.2_AT, Flower.3_AT, Flower.1_AL, Flower.2_AL, 
            Flower.3_AL, Flower.1_CR, Flower.2_CR, Flower.3_CR, Flower.1_ES, Flower.2_ES, Flower.3_ES, 
            Flower.1_TH, Flower.2_TH, Flower.3_TH, Flower.1_MT, Flower.2_MT, Flower.3_MT, Flower.1_BD, 
            Flower.2_BD, Flower.3_BD))

    } else if ((dataset_id == "DevSeq") && (devseq_spec == "Brassicaceae")) {
        x_df <- x_df %>% select (-c(Pollen.1_AT, Pollen.2_AT, Pollen.3_AT, Pollen.1_AL, Pollen.2_AL, 
            Pollen.3_AL, Pollen.1_CR, Pollen.2_CR, Pollen.3_CR, Pollen.1_ES, Pollen.2_ES, Pollen.3_ES))
    }




#------------------------ Prepare DevSeq and Brawand data for corrplot -------------------------


    if (((dataset_id == "DevSeq") && (devseq_spec == "all")) || dataset_id == "Brawand") {

        # Log2 transform data if expr_estimation=TPM is chosen
        # Compute correlation and build distance matrix
        # get_dist is an wrapper implemented in the factoextra package around
        # as.dist(1-cor(scale(df))) if scale = TRUE; or
        # as.dist(1-cor(df)) if scale = FALSE
        if (is.element("pearson", coefficient) && is.element("counts", expr_estimation)) {
            x_cor <- cor(x_df[, 2:ncol(x_df)], method = "pearson")
            x_dist <- as.dist(sqrt(1/2*(1-cor(x_df[, 2:ncol(x_df)], method = "pearson"))))
            scale_data <- TRUE 
            # correlation matrix does not need to be transposed for get_dist method (since ncol=nrow)

        } else if (is.element("pearson", coefficient) && is.element("TPM", expr_estimation)) {
            x_df[,2:ncol(x_df)] <- log2(x_df[,2:ncol(x_df)] + 1)
            x_cor <- cor(x_df[, 2:ncol(x_df)], method = "pearson")
            x_dist <- as.dist(sqrt(1/2*(1-cor(x_df[, 2:ncol(x_df)], method = "pearson"))))
            scale_data <- TRUE

        } else if (is.element("spearman", coefficient)) {
            x_cor <- cor(x_df[, 2:ncol(x_df)], method = "spearman")
            x_dist <- as.dist(sqrt(1/2*(1-cor(x_df[, 2:ncol(x_df)], method = "spearman"))))
            scale_data <- FALSE 
        }     



        df_clust.res <- hclust(x_dist, method = "average") # agglomerate clustering

        # Get dendrogram w/ species colors
        col_dend <- dendrapply(as.dendrogram(df_clust.res), function(y){
        
            if (is.leaf(y)){
                dend_col <- species_col[substr(attr(y, "label"), 1, 3)]
                attr(y, "edgePar") <- list(col = dend_col) # color branch
            }   
            return(y)
        })

        col_labels <- get_leaves_branches_col(col_dend) # get branch colors
        col_cols <- col_labels[order(order.dendrogram(col_dend))] # order color vector

        # Get dendrogram w/ experiment colors
        row_dend <- dendrapply(as.dendrogram(df_clust.res), function(z){
        
            if (is.leaf(z)){
                name_string <- attr(z, "label")
                dend_col <- exp_col[substr(name_string, (nchar(name_string))-1, nchar(name_string))]
                attr(z, "edgePar") <- list(col = dend_col) # color branch
            }  
            return(z)
        })

        row_labels <- get_leaves_branches_col(row_dend) # get branch colors
        row_cols <- row_labels[order(order.dendrogram(row_dend))] # order color vector

    
        # Rotate leaves of dendrogram so that clusters appear in evolutionary order wherevever possible
        # This also re-orders the colors of the colorbar where possible to make them more distinguishable 
        if (is.element("pearson", coefficient) && is.element("counts", expr_estimation) && 
            (dataset_id == "DevSeq") && is.element("inter-organ", data_norm) && is.element("all", devseq_organs)) {

            dend_order=dendextend::rotate(as.dendrogram(df_clust.res),c(145:147,151:153,148:150,154:168,1:18,22:24,19:21,25:30,37:42,34:36,31:33,46:57,43:45,58:72,76:78,73:75,79:84,88:90,85:87,91:144))

        } else if (is.element("pearson", coefficient) && is.element("counts", expr_estimation) && 
            (dataset_id == "DevSeq") && is.element("inter-organ", data_norm) && is.element("subset", devseq_organs)) {

            dend_order=dendextend::rotate(as.dendrogram(df_clust.res),c(1:24,28:30,25:27,43:87,91:93,88:90,94:99,103:105,100:102,31:33,40:42,37:39,34:36))

        } else if (is.element("spearman", coefficient) && is.element("counts", expr_estimation) && 
            (dataset_id == "DevSeq") && is.element("inter-organ", data_norm) && is.element("all", devseq_organs)) {

            dend_order = dendextend::rotate(as.dendrogram(df_clust.res),c(1:3,7:9,4:6,10:72,76:84,73:75,88:90,85:87,91:99,103:105,100:102,106:129,130:156,160:168,157:159))

        } else if (is.element("spearman", coefficient) && is.element("counts", expr_estimation) && 
            (dataset_id == "DevSeq") && is.element("inter-organ", data_norm) && is.element("subset", devseq_organs)) {

            dend_order=dendextend::rotate(as.dendrogram(df_clust.res),c(91:105,1:9,13:15,10:12,28:36,40:42,37:39,46:48,52:54,49:51,43:45,58:60,64:66,61:63,55:57,67:90,16:27))

        } else if (is.element("pearson", coefficient) && is.element("counts", expr_estimation) && 
            (dataset_id == "Brawand") && is.element("inter-organ", data_norm)) {

            dend_order = dendextend::rotate(as.dendrogram(df_clust.res),c(1:13,20:23,16:17,14:15,18:19,45:47,42:44,24:27,29,28,30:33,35,34,36,40:41,37:39))

        } else if (is.element("spearman", coefficient) && is.element("counts", expr_estimation) && 
            (dataset_id == "Brawand") && is.element("inter-organ", data_norm)) {

            dend_order = dendextend::rotate(as.dendrogram(df_clust.res),c(1:13,22:23,20:21,14:19,45:47,42:44,24,29,28,25,26,27,30:36,41,40,37:39))

        } else dend_order = TRUE 

    } 


#---------------------------- Make corrplot for DevSeq and Brawand -----------------------------



    if ((dataset_id == "DevSeq") && (is.element("all", devseq_spec)) && is.element("inter-organ", data_norm)
        && is.element("counts", expr_estimation)) {

        # Make corrplots
        png(height = 3500, width = 3500, pointsize = 20, file = file.path(out_dir, "output", "plots", fname))
        par(lwd = 18.5) # dendrogram line width
        getRowOrder = heatmap.2(x_cor,
            revC = F,
            ColSideColors = col_cols, 
            RowSideColors = row_cols,
            # Order dendrogram in a way that it starts with distant species BD, MT, TH
            Rowv = dend_order, 
            Colv = "Rowv"
            )

        # Get order of rows and rearrange "row_cols" vector
        # fixes gplots heatmap.2 RowSideColors bug (colorbar does not reverse when revC=T)
        ordinary_order = getRowOrder$rowInd
        reversal = cbind(ordinary_order, rev(ordinary_order))
        rev_col = row_cols[reversal[,2]]; rev_col = rev_col[order(reversal[,1])];

        # Create heatmap with reversed RowSideColors
        heatmap.2(x_cor, 
            revC = T,
            ColSideColors = col_cols, 
            RowSideColors = rev_col, 
            density.info = "none",
            trace = "none",
            col = pal(800),
            cexRow = 2,
            cexCol = 2,
            margins = c(12, 12),
            key = FALSE,
            lwid = c(0.2,2.3,28.5), # column width
            lhei = c(0.2,2.3,28.5), # column height
            offsetRow = 1,
            offsetCol = 1,
            key.xlab = NA,
            key.title = NULL,
            # Order dendrogram in a way that it starts with distant species BD, MT, TH
            Rowv = dend_order, 
            Colv = "Rowv"
            )

        dev.off()

        # Save colorbar
        png(height = 1000, width = 850, pointsize = 20, file = file.path(out_dir, "output", "plots", "DevSeq_comp_colorbar.png"))
        par(cex = 7, mar = c(1,2.75,1,1), cex.lab = 1.25)
        x_c = 1
        y_c = seq(0,1,len = 100)
        z_c = matrix(1:100, nrow = 1)
        ylabel <- seq(0, 1, by = 0.5)
        image(x_c,y_c,z_c,col = pal(800), axes = FALSE, xlab = "", ylab = "", cex=10)
        axis(2, at = ylabel, las = 1, lwd = 15)
        dev.off()

    } else if ((dataset_id == "Brawand") && is.element("inter-organ", data_norm)) {

        # Make corrplots
        png(height = 3500, width = 3500, pointsize = 20, file = file.path(out_dir, "output", "plots", fname))
        par(lwd = 18.5) # dendrogram line width
        getRowOrder = heatmap.2(x_cor,
            revC = F,
            ColSideColors = col_cols, 
            RowSideColors = row_cols,
            Rowv = dend_order, 
            Colv = "Rowv")

        # Get order of rows and rearrange "row_cols" vector
        # fixes gplots heatmap.2 RowSideColors bug (colorbar does not reverse when revC=T)
        ordinary_order = getRowOrder$rowInd
        reversal = cbind(ordinary_order, rev(ordinary_order))
        rev_col = row_cols[reversal[,2]]; rev_col = rev_col[order(reversal[,1])];

        # Create heatmap with reversed RowSideColors
        heatmap.2(x_cor, 
            revC = T,
            ColSideColors = col_cols, 
            RowSideColors = rev_col, 
            density.info = "none",
            trace = "none",
            col = pal(800),
            Rowv = dend_order, 
            Colv = "Rowv", 
            cexRow = 2,
            cexCol = 2,
            margins = c(12, 12),
            key = FALSE,
            lwid = c(0.2,2.3,28.5), # column width
            lhei = c(0.2,2.3,28.5), # column height
            offsetRow = 1,
            offsetCol = 1,
            key.xlab = NA,
            key.title = NULL)
        dev.off()

        # Save colorbar
        png(height = 1000, width = 850, pointsize = 20, file = file.path(out_dir, "output", "plots", "Brawand_comp_colorbar.png"))
        par(cex = 7, mar = c(1,2.75,1,1), cex.lab = 1.25)
        x_c = 1
        y_c = seq(0,1,len = 100)
        z_c = matrix(1:100, nrow = 1)
        ylabel <- seq(0, 1, by = 0.5)
        image(x_c,y_c,z_c,col = pal(800), axes = FALSE, xlab = "", ylab = "", cex=10)
        axis(2, at = ylabel, las = 1, lwd = 15)
        dev.off()

    }
    



#-------------------------------------- Merge replicates ---------------------------------------

    # Prepare data for inter-organ distance analysis (DevSeq data set all species)
    # Perform this part of analysis with TPM data, not with VST counts

    if ((dataset_id == "DevSeq") && (devseq_spec == "all") && (expr_estimation == "TPM") && (devseq_organs == "all")) {


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


        # log-transform data if TPM and Pearson are chosen
        if (is.element("pearson", coefficient) && is.element("TPM", expr_estimation)) {
            x_repl <- x
            x_repl[,2:ncol(x_repl)] <- log2(x_repl[,2:ncol(x_repl)] + 1)

        } else {
            x_repl <- x
        }


        x_avg <- calculateAvgExpr(x_repl)


        DevSeq_col_names <- rep(c("Root", "Hypocotyl", "Leaf", "veg_apex", "inf_apex", "Flower", 
            "Stamen", "Carpel", "Pollen"), each=7)
        DevSeq_spec_names <- rep(c("_AT", "_AL", "_CR", "_ES", "_TH", "_MT", "_BD"), times=9)
        repl_names <- paste0(DevSeq_col_names, DevSeq_spec_names)

        colnames(x_avg)[2:ncol(x_avg)] <- repl_names




        #------ Get replicate correlations and generate DevSeq inter-organ distance plot -------


        x_avg <- x_avg %>% select(-c(Pollen_AT, Pollen_AL, Pollen_CR, Pollen_ES, Pollen_TH, 
                 Pollen_MT, Pollen_BD))


        # Reorder data frame columns
        x_avg <- x_avg[,c(1, #gene_id column
                          2,9,16,23,30,37,44,51, #AT
                          3,10,17,24,31,38,45,52, #AL
                          4,11,18,25,32,39,46,53, #CR
                          5,12,19,26,33,40,47,54, #ES
                          6,13,20,27,34,41,48,55, #TH
                          7,14,21,28,35,42,49,56, #MT
                          8,15,22,29,36,43,50,57)] #BD


        if (coefficient == "pearson") {
            x_cor_avg <- cor(x_avg[2:ncol(x_avg)], method = "pearson") 
        } else if (coefficient == "spearman") {
            x_cor_avg <- cor(x_avg[2:ncol(x_avg)], method = "spearman") 
        }


        AT <- as.data.frame(unique(as.vector(x_cor_avg[1:8,1:8])))[-1,]
        AL <- as.data.frame(unique(as.vector(x_cor_avg[9:16,9:16])))[-1,]
        CR <- as.data.frame(unique(as.vector(x_cor_avg[17:24,17:24])))[-1,]
        ES <- as.data.frame(unique(as.vector(x_cor_avg[25:32,25:32])))[-1,]
        TH <- as.data.frame(unique(as.vector(x_cor_avg[33:40,33:40])))[-1,]
        MT <- as.data.frame(unique(as.vector(x_cor_avg[41:48,41:48])))[-1,]
        BD <- as.data.frame(unique(as.vector(x_cor_avg[49:56,49:56])))[-1,]

        df_names <- c("Species" , "Correlation")

        species_DevSeq <- c("AT", "AL", "CR", "ES", "TH", "MT", "BD")
        species_names <- as.data.frame(rep(species_DevSeq, each=28))
        organ_cor <- as.data.frame(c(AT,AL,CR,ES,TH,MT,BD))

        cor_df <- cbind(species_names, organ_cor)
        colnames(cor_df) <- df_names

        # Cor values of meristematic tssues apex_inf apex_veg and carpel
        AT_m <- as.data.frame(unique(as.vector(x_cor_avg[c(4:5,8),c(4:5,8)])))[-1,]
        AL_m <- as.data.frame(unique(as.vector(x_cor_avg[c(12:13,16),c(12:13,16)])))[-1,]
        CR_m <- as.data.frame(unique(as.vector(x_cor_avg[c(20:21,24),c(20:21,24)])))[-1,]
        ES_m <- as.data.frame(unique(as.vector(x_cor_avg[c(28:29,32),c(28:29,32)])))[-1,]
        TH_m <- as.data.frame(unique(as.vector(x_cor_avg[c(36:37,40),c(36:37,40)])))[-1,]
        MT_m <- as.data.frame(unique(as.vector(x_cor_avg[c(44:45,48),c(44:45,48)])))[-1,]
        BD_m <- as.data.frame(unique(as.vector(x_cor_avg[c(52:53,56),c(52:53,56)])))[-1,]

        species_names_m <- as.data.frame(rep(species_DevSeq, each=3))
        organ_cor_m <- as.data.frame(c(AT_m,AL_m,CR_m,ES_m,TH_m,MT_m,BD_m))

        cor_df_m <- cbind(species_names_m, organ_cor_m)
        colnames(cor_df_m) <- df_names


        cor_df$Correlation <- as.numeric(cor_df$Correlation)
        cor_df_m$Correlation <- as.numeric(cor_df_m$Correlation)


        # Make inter-organ distance plot
        plotOrganDist <- function(data, data_m, data_norm) {

            fname <- paste('DevSeq_inter_organ_dist_', coefficient, '.jpg', sep="")

            coef_lab <- ifelse(coefficient == "pearson", "Pearson's r", "Spearman's rho")

            y_max = 1.035

            if((coefficient == "pearson") && (data_norm == "intra-organ")) {
                y_min = 0.37
            } else if((coefficient == "pearson") && (data_norm == "inter-organ")) {
                y_min = 0.364
            } else if((coefficient == "spearman") && (data_norm == "intra-organ")) {
                y_min = 0.3575
            } else if((coefficient == "spearman") && (data_norm == "inter-organ")) {
                y_min = 0.357
            }

            cor_colors <- rep(c("#808080","#b1971e","#419730","#2e9d76","#3889af","#7d67c7","#b72f2d"), each=28)

            p <- ggplot(data, aes(x=factor(Species, levels=unique(Species)), y=Correlation, fill=Species)) + 
            stat_boxplot(geom ='errorbar', width = 0.45, size=1.0, color="gray15") + 
            geom_boxplot(width = 0.75, size=1.0, color="gray15", outlier.shape = 21, 
                outlier.size = 1.5, outlier.stroke = 1.5, outlier.fill = NA, outlier.color="gray35") + 
            geom_beeswarm(data = data, priority = c("ascending"), colour=cor_colors, cex=2, size=3.25) +
            geom_beeswarm(data = data_m, priority = c("ascending"), colour="red", cex=2, size=2.5, shape=8, stroke=2) + 
            scale_y_continuous(limits = c(y_min, y_max), expand = c(0, 0)) + 
            annotate("rect", xmin=0.35, xmax=7.65, ymin=y_min, ymax=y_max, fill="white", alpha=0, 
                color="gray15", size=1.35)

            q <- p + scale_fill_manual(values=c("#f0dd91","#cacaca","#ebb4b1","#b3daad","#ace1ce","#ccc3f0","#b2d7e7")) + 
            theme_minimal() + 
            xlab("Species") + 
            ylab("") + 
            theme(legend.position = "none", 
                text=element_text(size=20), 
                panel.grid.major.y = element_line(size = 0.75, color = c("gray85")),
                panel.grid.minor.y = element_line(size = 0.75, color = c("gray85")), 
                panel.grid.major.x = element_line(size = 0.75, color = c("gray85")), 
                panel.border = element_rect(colour = "black", fill=NA, size=1.7),
                axis.ticks.length = unit(.35, "cm"),
                axis.ticks = element_line(colour = "black", size = 0.725), 
                axis.title.x = element_text(colour = "black", size=24, 
                    margin = margin(t = 12, r = 0, b = 0, l = 0)), 
                axis.title.y = element_text(colour = "black", size=24, 
                    margin = margin(t = 0, r = 0, b = 0, l = -11.5)), 
                axis.text.x = element_text(colour = "black", size=21, angle=0, 
                    margin = margin(t = 5, r = 0, b = 0, l = 0)),
                axis.text.y = element_text(colour = "black", size=21, 
                    margin = margin(t = 0, r = 5, b = 0, l = 0)),  
                plot.margin = unit(c(0.55, 0.1, 1.74, 0), "cm"))

            ggsave(file = file.path(out_dir, "output", "plots", fname), plot = q,
                width = 7, height = 8.0, units = c("in"), dpi = 300, limitsize = FALSE)
        }

        plotOrganDist(data=cor_df, data_m=cor_df_m, data_norm=data_norm)

    }




#---------------- Get gene expression divergence rates for ATH/AL vs species X -----------------


   # Use pearson correlation, intra-organ normalization and TPM
   # Use previously merged replicates of DevSeq data including pollen sampless

   if ((dataset_id == "DevSeq") && (devseq_spec == "all") && (expr_estimation == "TPM") && (devseq_organs == "all")) {

      getOrganCor <- function(df, organ, coefficient, expr_estimation) {

         # log-transform data if TPM and Pearson are chosen
         if ((coefficient == "pearson") && (expr_estimation == "TPM")) {
            df <- log2(df + 1)
         }

         df_cor <- cor(df, method=coefficient)
         df_cor <- df_cor[4:nrow(df_cor), 1:3]

         getError <- function(cor_data) {
            std <- sd(cor_data)
            n_value <- length(cor_data)
            error <- qt(0.995, df = n_value-1) * std/sqrt(n_value)
            return(error)
         }

         sp1 <- mean(c(df_cor[1:3,]))
         sp1_li <- sp1 - getError(c(df_cor[1:3,]))
         sp1_ri <- sp1 + getError(c(df_cor[1:3,]))

         sp2 <- mean(c(df_cor[4:6,]))
         sp2_li <- sp2 - getError(c(df_cor[4:6,]))
         sp2_ri <- sp2 + getError(c(df_cor[4:6,]))

         sp3 <- mean(c(df_cor[7:9,]))
         sp3_li <- sp3 - getError(c(df_cor[7:9,]))
         sp3_ri <- sp3 + getError(c(df_cor[7:9,]))

         sp4 <- mean(c(df_cor[10:12,]))
         sp4_li <- sp4 - getError(c(df_cor[10:12,]))
         sp4_ri <- sp4 + getError(c(df_cor[10:12,]))

         sp5 <- mean(c(df_cor[13:15,]))
         sp5_li <- sp5 - getError(c(df_cor[13:15,]))
         sp5_ri <- sp5 + getError(c(df_cor[13:15,]))

         sp6 <- mean(c(df_cor[16:18,]))
         sp6_li <- sp6 - getError(c(df_cor[16:18,]))
         sp6_ri <- sp6 + getError(c(df_cor[16:18,]))

         df_cor_avg <- rbind(sp1, sp2, sp3, sp4, sp5, sp6)
         colnames(df_cor_avg) <- organ
         lower <- rbind(sp1_li, sp2_li, sp3_li, sp4_li, sp5_li, sp6_li)
         colnames(lower) <- "lower"
         upper <- rbind(sp1_ri, sp2_ri, sp3_ri, sp4_ri, sp5_ri, sp6_ri)
         colnames(upper) <- "upper"

         getRowNames = function(x,n){ substring(x,nchar(x)-n+1) }
         row_names_repl <- getRowNames(rownames(df_cor),2)
         rnames_div_rates <- unique(row_names_repl)
         rownames(df_cor_avg) <- rnames_div_rates
         df_cor_avg <- cbind(df_cor_avg, lower, upper)

         return(df_cor_avg)

      }

      root_div <- getOrganCor(df=x[,2:22], organ="Root  ", coefficient=coefficient, expr_estimation=expr_estimation)
      hypocotyl_div <- getOrganCor(df=x[,23:43], organ="Hypocotyl  ", coefficient=coefficient, expr_estimation=expr_estimation)
      leaf_div <- getOrganCor(df=x[,44:64], organ="Leaf  ", coefficient=coefficient, expr_estimation=expr_estimation)
      veg_apex_div <- getOrganCor(df=x[,65:85], organ="Apex veg  ", coefficient=coefficient, expr_estimation=expr_estimation)
      inf_apex_div <- getOrganCor(df=x[,86:106], organ="Apex inf  ", coefficient=coefficient, expr_estimation=expr_estimation)
      flower_div <- getOrganCor(df=x[,107:127], organ="Flower  ", coefficient=coefficient, expr_estimation=expr_estimation)
      stamen_div <- getOrganCor(df=x[,128:148], organ="Stamen  ", coefficient=coefficient, expr_estimation=expr_estimation)
      carpel_div <- getOrganCor(df=x[,149:169], organ="Carpel  ", coefficient=coefficient, expr_estimation=expr_estimation)
      pollen_div <- getOrganCor(df=x[,170:190], organ="Pollen  ", coefficient=coefficient, expr_estimation=expr_estimation)


      # Reshape data table for ggplot
      # divergence times are estimated taxon pair times from TimeTree
      # http://www.timetree.org/
      div_times <- rep(c(7.1, 9.4, 25.6, 46, 106, 160), times=9)
      comp_organ <- rep(c(colnames(root_div)[1], colnames(hypocotyl_div)[1], colnames(leaf_div)[1], 
        colnames(veg_apex_div)[1], colnames(inf_apex_div)[1], colnames(flower_div)[1], 
        colnames(stamen_div)[1], colnames(carpel_div)[1], colnames(pollen_div)[1]), each=6)
      comp_spec <- c(rownames(root_div), rownames(hypocotyl_div), rownames(leaf_div), rownames(veg_apex_div), 
        rownames(inf_apex_div), rownames(flower_div), rownames(stamen_div), rownames(carpel_div), 
        rownames(pollen_div))

      DevSeq_GE_div <- rbind(root_div, hypocotyl_div, leaf_div, veg_apex_div, inf_apex_div, 
        flower_div, stamen_div, carpel_div, pollen_div)
      rownames(DevSeq_GE_div) <- NULL
      colnames(DevSeq_GE_div) <- c("correlation", "lower", "upper")

      DevSeq_div_rates <- data.frame(cbind(comp_spec, comp_organ, div_times, DevSeq_GE_div), 
        stringsAsFactors=FALSE)

      DevSeq_div_rates$div_times <- as.numeric(DevSeq_div_rates$div_times)
      DevSeq_div_rates$correlation <- as.numeric(DevSeq_div_rates$correlation)
      DevSeq_div_rates$lower <- as.numeric(DevSeq_div_rates$lower)
      DevSeq_div_rates$upper <- as.numeric(DevSeq_div_rates$upper)
      
      # Change order of organs in df
      DevSeq_div_rates <- DevSeq_div_rates[c(7:12,37:42,31:36,1:6,19:30,43:48,13:18,49:54),]
      DevSeq_div_rates_pollen <- DevSeq_div_rates[49:54,]
      DevSeq_div_rates$comp_organ <- factor(DevSeq_div_rates$comp_organ, 
        levels = unique(DevSeq_div_rates$comp_organ))
      DevSeq_div_rates_pollen$comp_organ <- factor(DevSeq_div_rates_pollen$comp_organ, 
        levels = unique(DevSeq_div_rates_pollen$comp_organ))

      

      # Make GE divergence plot
      makeGEDivPlot <- function(data1, data2, plot_title, coefficient) {

        fname <- sprintf('%s.jpg', paste("GE_divergence_rates", coefficient, sep="_"))

        p <- ggplot(data=data1, aes(x=div_times, y=correlation, group=comp_organ, colour=comp_organ)) + 
        geom_ribbon(aes(ymin = data1$lower, ymax = data1$upper, fill= comp_organ), alpha = 0.25, 
            linetype = 0, show.legend = FALSE) + 
        scale_fill_manual(values = c("Hypocotyl  "="#53b0db", "Stamen  "="#ee412e", "Flower  "="#e075af", 
                "Root  "="#6a54a9", "Apex veg  "="#96ba37", "Apex inf  "="#fad819", "Carpel  "="#f2a72f", 
                "Leaf  "="#2c8654", "Pollen  "="#a63126")) + 
        geom_line(size = 3.1) +  
        scale_x_continuous(limits = c(7,160), expand = c(0.02,0), breaks = c(7,9,25,46,106,160)) + 
        scale_y_continuous(limits = c(0.4375, 0.9075), expand = c(0.02, 0)) + 
        scale_color_manual(values = c("#53b0db", "#ee412e", "#e075af", "#6a54a9", "#96ba37", "#fad819", 
            "#f2a72f", "#2c8654", "#a63126"), 
            # organ order: hypocotyl/stamen/flower/root/veg_apex/inf_apex/carpel/leaf
            breaks=c("Root  ", "Hypocotyl  ", "Leaf  ", "Apex veg  ", "Apex inf  ", "Flower  ", 
                "Stamen  ", "Carpel  ", "Pollen  ")) + 
        geom_line(aes(x=div_times, y=correlation), data=data2, color = "white", lty = "solid", 
            lwd = 3.1) + # pollen
        geom_line(aes(x=div_times, y=correlation), data=data2, color = "#a63126", lty = "22", 
            lwd = 3.1) + # pollen
        geom_segment(x=157.5, xend=157.5, y=0.421, yend=0.461, color="white", size=12.5) + 
        annotate("text", x=19.5, y=0.4545, label= "Brassicaceae", size=8) + 
        annotate("text", x=46, y=0.4545, label= "TH", size=8) + 
        annotate("text", x=106, y=0.4545, label= "MT", size=8) + 
        annotate("text", x=158.6, y=0.4545, label= "BD", size=8) + 
        geom_segment(x=7, xend=7, y=0.405, yend=0.44, color="black", size=0.7) + 
        geom_segment(x=9, xend=9, y=0.405, yend=0.44, color="black", size=0.7) + 
        geom_segment(x=25, xend=25, y=0.405, yend=0.44, color="black", size=0.7) + 
        geom_segment(x=46, xend=46, y=0.405, yend=0.44, color="black", size=0.7) + 
        geom_segment(x=106, xend=106, y=0.405, yend=0.44, color="black", size=0.7) + 
        geom_segment(x=160, xend=160, y=0.405, yend=0.44, color="black", size=0.7) + 
        guides(color = guide_legend(ncol = 3))

        q <- p + theme_bw() + xlab("Divergence time from A.thaliana (Myr)") + ylab("Pearson's r w/ A.thaliana") + 
        theme(text=element_text(size=16), 
            axis.ticks.length=unit(0.35, "cm"), 
            axis.ticks = element_line(colour = "black", size = 0.775),  
            plot.margin = unit(c(0.55, 1.175, 0.5, 0.4),"cm"), 
            axis.title.y = element_text(size=25, margin = margin(t = 0, r = 14.5, b = 0, l = 11.5), colour="black"), 
            axis.title.x = element_text(size=25, margin = margin(t = 13.25, r = 0, b = 3.5, l = 0), colour="black"), 
            axis.text.x = element_text(size=21.25, angle=0, margin = margin(t = 5.5), colour="black"), 
            axis.text.y = element_text(size=21.25, angle=0, margin = margin(r = 5.5), colour="black"), 
            legend.box.background = element_rect(colour = "#d5d5d5", fill=NA, size=1.0), 
            panel.border = element_rect(colour = "black", fill=NA, size=1.75), 
            panel.grid.major.x = element_line(color="#d5d5d5"),
            panel.grid.minor.x = element_blank(), 
            panel.grid.minor.y = element_blank(), 
            panel.grid.major.y = element_blank(), 
            legend.position = c(0.723, 0.88), 
            legend.title = element_blank(), 
            legend.text = element_text(size=21.5), 
            legend.spacing.x = unit(0.5, 'cm'), 
            legend.key.size = unit(0.95, "cm"), 
            legend.background=element_blank()) 

        ggsave(file = file.path(out_dir, "output", "plots", fname), plot = q, 
            width = 12.535, height = 8, dpi = 300, units = c("in"), limitsize = FALSE) 
      }

      makeGEDivPlot(data1 = DevSeq_div_rates, data2 = DevSeq_div_rates_pollen, 
        coefficient = coefficient)


   }
   



#------------------------------------------ Make PCAs ------------------------------------------


  if ((dataset_id == "DevSeq") && (devseq_organs == "all")) { 
    performPCA <- function(data, data_scale = c("raw", "scaled"), ntop) {

    makePCA <- function(df, pca_dim = c("1_2", "2_3", "1_3"), data_scale, ntop) {

        # Get highly variable genes
        x_mat <- data.matrix(x_df[,2:ncol(x_df)])
        rv <- rowVars(x_mat)
        rv <- as.data.frame(rv)
        colnames(rv) <- "rowVars"

        x_df <- cbind(x_df, rv)
        x_df <- x_df[order(x_df$rowVars, decreasing = TRUE), ]
        x_df <- subset(x_df, select = -rowVars)
        x_df <- x_df[1:ntop,]


        # Transpose data
        x_df <- t(x_df[,-1])

        if (data_scale == "raw") {
            scale_v <- FALSE
        } else scale_v <- TRUE


        # Perform PCA analysis on log-transformed data (if TPM chosen) or rlog counts
        comp_pca <- prcomp(x_df, center = TRUE, scale = scale_v) 

        # Get eigenvalues, explained variance (%) and cumulative variance (%) 
        eig_val <- get_eigenvalue(comp_pca)
        pc1_var <- round((eig_val$variance.percent[1]), digits = 1) # variance explained by PC1
        pc2_var <- round((eig_val$variance.percent[2]), digits = 1) # variance explained by PC2
        pc3_var <- round((eig_val$variance.percent[3]), digits = 1) # variance explained by PC3

        # Get PCA coordinates of individuals
        pc_scores <- get_pca_ind(comp_pca)          # Show output results
        pca_coord <- as.data.frame(pc_scores$coord) # Get coordinates as dataframe

        if (is.element("1_2", pca_dim)) {
            pca_coord <- pca_coord[,1:2]

        } else if (is.element("2_3", pca_dim)) {
            pca_coord <- pca_coord[,2:3]
        
        } else if (is.element("1_3", pca_dim)) {
            pca_coord <- pca_coord[,c(1,3)]
        }

        pca_sample_names <- rownames(pca_coord)

        # Get last two characters of each string from pca_sample_names
        substrRight <- function(x, n){ 
            substr(x, nchar(x)-n+1, nchar(x))
        }

        pca_species <- substrRight(pca_sample_names, 2)

        # Get organ names from pca_sample_names
        pca_organs <- gsub('.{5}$', '', pca_sample_names)

        pca_df <- cbind(pca_coord, pca_species, pca_organs)
        colnames(pca_df) <- c("PC1", "PC2", "Species", "Organ")

        return_list <- list("pca_df"=pca_df, "pc1_var"=pc1_var, "pc2_var"=pc2_var, "pc3_var"=pc3_var, "eig_val"=eig_val)
        return(return_list)
    }



    if (is.element("Brassicaceae", devseq_spec)) {

        # Make PCA for DevSeq w/ stamen data PC1/2
        pca_return_objects <- makePCA(df = x_df, pca_dim = "1_2", data_scale = data_scale, ntop = ntop)
        list2env(pca_return_objects, envir = .GlobalEnv)
        DevSeq_pca_1_2_w_stamen <- pca_df
        DevSeq_pc1_var_w_stamen <- pc1_var
        DevSeq_pc2_var_w_stamen <- pc2_var
        DevSeq_pc3_var_w_stamen <- pc3_var

        # Make PCA for DevSeq w/ stamen data PC2/3
        pca_return_objects <- makePCA(df = x_df, pca_dim = "2_3", data_scale = data_scale, ntop = ntop)
        list2env(pca_return_objects, envir = .GlobalEnv)
        DevSeq_pca_2_3_w_stamen <- pca_df

        pca_return_objects <- makePCA(df = x_df, pca_dim = "1_3", data_scale = data_scale, ntop = ntop)
        list2env(pca_return_objects, envir = .GlobalEnv)
        DevSeq_pca_1_3_w_stamen <- pca_df

        # Replace species abbreviations by more detailed names
        DevSeq_pca_1_2_w_stamen$Species <- DevSeq_pca_1_2_w_stamen$Species %>% gsub('AT', 'A.thaliana', .) %>% 
        gsub('AL', 'A.lyrata', .) %>% gsub('CR', 'C.rubella', .) %>% gsub('ES', 'E.salsug.', .)

        DevSeq_pca_2_3_w_stamen$Species <- DevSeq_pca_2_3_w_stamen$Species %>% gsub('AT', 'A.thaliana', .) %>% 
        gsub('AL', 'A.lyrata', .) %>% gsub('CR', 'C.rubella', .) %>% gsub('ES', 'E.salsug.', .)

        DevSeq_pca_1_3_w_stamen$Species <- DevSeq_pca_1_3_w_stamen$Species %>% gsub('AT', 'A.thaliana', .) %>% 
        gsub('AL', 'A.lyrata', .) %>% gsub('CR', 'C.rubella', .) %>% gsub('ES', 'E.salsug.', .)

        # Modify organ names
        DevSeq_pca_1_2_w_stamen$Organ <- DevSeq_pca_1_2_w_stamen$Organ %>% gsub('veg_apex', 'Apex veg', .) %>% 
        gsub('inf_apex', 'Apex inf', .)

        DevSeq_pca_2_3_w_stamen$Organ <- DevSeq_pca_2_3_w_stamen$Organ %>% gsub('veg_apex', 'Apex veg', .) %>% 
        gsub('inf_apex', 'Apex inf', .)

        DevSeq_pca_1_3_w_stamen$Organ <- DevSeq_pca_1_3_w_stamen$Organ %>% gsub('veg_apex', 'Apex veg', .) %>% 
        gsub('inf_apex', 'Apex inf', .)

    } else if (is.element("all", devseq_spec)) { 
        
        # Make PCA for DevSeq w/ stamen data PC1/2
        pca_return_objects <- makePCA(df = x_df, pca_dim = "1_2", data_scale = data_scale, ntop = ntop)
        list2env(pca_return_objects, envir = .GlobalEnv)
        DevSeq_pca_1_2_w_stamen <- pca_df
        DevSeq_pc1_var_w_stamen <- pc1_var
        DevSeq_pc2_var_w_stamen <- pc2_var
        DevSeq_pc3_var_w_stamen <- pc3_var

        # Make PCA for DevSeq w/ stamen data PC2/3
        pca_return_objects <- makePCA(df = x_df, pca_dim = "2_3", data_scale = data_scale, ntop = ntop)
        list2env(pca_return_objects, envir = .GlobalEnv)
        DevSeq_pca_2_3_w_stamen <- pca_df

        pca_return_objects <- makePCA(df = x_df, pca_dim = "1_3", data_scale = data_scale, ntop = ntop)
        list2env(pca_return_objects, envir = .GlobalEnv)
        DevSeq_pca_1_3_w_stamen <- pca_df

        # Replace species abbreviations by more detailed names
        DevSeq_pca_1_2_w_stamen$Species <- DevSeq_pca_1_2_w_stamen$Species %>% gsub('AT', 'A.thaliana', .) %>% 
        gsub('AL', 'A.lyrata', .) %>% gsub('CR', 'C.rubella', .) %>% gsub('ES', 'E.salsug.', .)  %>% 
        gsub('TH', 'T.hassler.', .) %>% gsub('MT', 'M.truncat.', .) %>% gsub('BD', 'B.distach.', .)

        DevSeq_pca_2_3_w_stamen$Species <- DevSeq_pca_2_3_w_stamen$Species %>% gsub('AT', 'A.thaliana', .) %>% 
        gsub('AL', 'A.lyrata', .) %>% gsub('CR', 'C.rubella', .) %>% gsub('ES', 'E.salsugineum', .)  %>% 
        gsub('TH', 'T.hassleriana', .) %>% gsub('MT', 'M.truncatula', .) %>% gsub('BD', 'B.distachyon', .)

        DevSeq_pca_1_3_w_stamen$Species <- DevSeq_pca_1_3_w_stamen$Species %>% gsub('AT', 'A.thaliana', .) %>% 
        gsub('AL', 'A.lyrata', .) %>% gsub('CR', 'C.rubella', .) %>% gsub('ES', 'E.salsugineum', .)  %>% 
        gsub('TH', 'T.hassleriana', .) %>% gsub('MT', 'M.truncatula', .) %>% gsub('BD', 'B.distachyon', .)

        # Modify organ names
        DevSeq_pca_1_2_w_stamen$Organ <- DevSeq_pca_1_2_w_stamen$Organ %>% gsub('veg_apex', 'Apex veg', .) %>% 
        gsub('inf_apex', 'Apex inf', .)

        DevSeq_pca_2_3_w_stamen$Organ <- DevSeq_pca_2_3_w_stamen$Organ %>% gsub('veg_apex', 'Apex veg', .) %>% 
        gsub('inf_apex', 'Apex inf', .)

        DevSeq_pca_1_3_w_stamen$Organ <- DevSeq_pca_1_3_w_stamen$Organ %>% gsub('veg_apex', 'Apex veg', .) %>% 
        gsub('inf_apex', 'Apex inf', .)
    }



    StatBag <- ggproto("Statbag", Stat,
                   compute_group = function(data, scales, prop = 0.5) {

                     #################################
                     #################################
                     # originally from aplpack package, plotting functions removed
                     plothulls_ <- function(x, y, fraction, n.hull = 1,
                                            col.hull, lty.hull, lwd.hull, density=0, ...){
                       # function for data peeling:
                       # x,y : data
                       # fraction.in.inner.hull : max percentage of points within the hull to be drawn
                       # n.hull : number of hulls to be plotted (if there is no fractiion argument)
                       # col.hull, lty.hull, lwd.hull : style of hull line
                       # plotting bits have been removed, BM 160321
                       # pw 130524
                       if(ncol(x) == 2){ y <- x[,2]; x <- x[,1] }
                       n <- length(x)
                       if(!missing(fraction)) { # find special hull
                         n.hull <- 1
                         if(missing(col.hull)) col.hull <- 1
                         if(missing(lty.hull)) lty.hull <- 1
                         if(missing(lwd.hull)) lwd.hull <- 1
                         x.old <- x; y.old <- y
                         idx <- chull(x,y); x.hull <- x[idx]; y.hull <- y[idx]
                         for( i in 1:(length(x)/3)){
                           x <- x[-idx]; y <- y[-idx]
                           if( (length(x)/n) < fraction ){
                             return(cbind(x.hull,y.hull))
                           }
                           idx <- chull(x,y); x.hull <- x[idx]; y.hull <- y[idx];
                         }
                       }
                       if(missing(col.hull)) col.hull <- 1:n.hull
                       if(length(col.hull)) col.hull <- rep(col.hull,n.hull)
                       if(missing(lty.hull)) lty.hull <- 1:n.hull
                       if(length(lty.hull)) lty.hull <- rep(lty.hull,n.hull)
                       if(missing(lwd.hull)) lwd.hull <- 1
                       if(length(lwd.hull)) lwd.hull <- rep(lwd.hull,n.hull)
                       result <- NULL
                       for( i in 1:n.hull){
                         idx <- chull(x,y); x.hull <- x[idx]; y.hull <- y[idx]
                         result <- c(result, list( cbind(x.hull,y.hull) ))
                         x <- x[-idx]; y <- y[-idx]
                         if(0 == length(x)) return(result)
                       }
                       result
                     } # end of definition of plothulls
                     #################################


                     # prepare data to go into function below
                     the_matrix <- matrix(data = c(data$x, data$y), ncol = 2)

                     # get data out of function as df with names
                     setNames(data.frame(plothulls_(the_matrix, fraction = prop)), nm = c("x", "y"))
                     # how can we get the hull and loop vertices passed on also?
                   },

                   required_aes = c("x", "y")
    )

    # inheritParams: ggplot2::stat_identity
    # param: prop Proportion of all the points to be included in the bag (default is 0.5)

    stat_bag <- function(mapping = NULL, data = NULL, geom = "polygon",
                     position = "identity", na.rm = FALSE, show.legend = NA, 
                     inherit.aes = TRUE, prop = 0.5, alpha = 0.3, ...) {
        layer(
            stat = StatBag, data = data, mapping = mapping, geom = geom, 
            position = position, show.legend = show.legend, inherit.aes = inherit.aes,
            params = list(na.rm = na.rm, prop = prop, alpha = alpha, ...)
        )
    }


    geom_bag <- function(mapping = NULL, data = NULL,
                     stat = "identity", position = "identity",
                     prop = 0.5, 
                     alpha = 0.3,
                     ...,
                     na.rm = FALSE,
                     show.legend = NA,
                     inherit.aes = TRUE) {
        layer(
            data = data,
            mapping = mapping,
            stat = StatBag,
            geom = GeomBag,
            position = position,
            show.legend = show.legend,
            inherit.aes = inherit.aes,
            params = list(
                na.rm = na.rm,
                alpha = alpha,
                prop = prop,
                ...
                )
            )
    }

    # rdname: ggplot2-ggproto
    # format: NULL
    # usage: NULL

    GeomBag <- ggproto("GeomBag", Geom,
                draw_group = function(data, panel_scales, coord) {
                    n <- nrow(data)
                    if (n == 1) return(zeroGrob())

                    munched <- coord_munch(coord, data, panel_scales)
                    # Sort by group to make sure that colors, fill, etc. come in same order
                    munched <- munched[order(munched$group), ]

                    # For gpar(), there is one entry per polygon (not one entry per point).
                    # We'll pull the first value from each group, and assume all these values
                    # are the same within each group.
                    first_idx <- !duplicated(munched$group)
                    first_rows <- munched[first_idx, ]

                    ggplot2:::ggname("geom_bag",
                        grid:::polygonGrob(munched$x, munched$y, default.units = "native",
                            id = munched$group,
                            gp = grid::gpar(
                                col = first_rows$colour,
                                fill = alpha(first_rows$fill, first_rows$alpha),
                                lwd = first_rows$size * .pt,
                                lty = first_rows$linetype
                            )
                        )
                    )
        },

        default_aes = aes(colour = "NA", fill = "grey20", size = 0.5, linetype = 1,
                          alpha = NA, prop = 0.5),
                          handle_na = function(data, params) {data},
                          required_aes = c("x", "y"),
                          draw_key = draw_key_polygon
    )



    # Set order for DevSeq organs
    DevSeq_pca_1_2_w_stamen$Organ <- factor(DevSeq_pca_1_2_w_stamen$Organ, c("Root", "Hypocotyl", 
            "Leaf", "Apex veg", "Apex inf", "Carpel", "Stamen", "Flower"))
    DevSeq_pca_2_3_w_stamen$Organ <- factor(DevSeq_pca_2_3_w_stamen$Organ, c("Root", "Hypocotyl", 
            "Leaf", "Apex veg", "Apex inf", "Carpel", "Stamen", "Flower"))
    DevSeq_pca_1_3_w_stamen$Organ <- factor(DevSeq_pca_1_3_w_stamen$Organ, c("Root", "Hypocotyl", 
            "Leaf", "Apex veg", "Apex inf", "Carpel", "Stamen", "Flower"))


    # Make PCA plots for main figure
    plotPCA <- function(data, pc_var1, pc_var2, set=c("pc1_2", "pc2_3", "pc1_3"), 
        spec=c("Brassicaceae", "all"), data_norm) {

        fname <- sprintf('%s.png', paste(deparse(substitute(data)), "pca", expr_estimation, spec, data_scale, sep="_"))
        
        if(set == "pc1_2") {
            x_coord <- "PC1 ("
            y_coord <- "PC2 ("
        } else if(set == "pc2_3") {
            x_coord <- "PC2 ("
            y_coord <- "PC3 ("
        } else if(set == "pc1_3") {
            x_coord <- "PC1 ("
            y_coord <- "PC3 ("
        }

        if(spec == "all") {
            spec_shape <- c(16, 17, 0, 18, 15, 2, 6)
            spec_shape_size <- c(5.0, 4.5, 3.5, 6.75, 4.5, 3.15, 3.15)
            legend_title <- element_blank()
            col_guide <- FALSE
            order_guide <- 0
            legend_spacing <- 15
        } else if(spec == "Brassicaceae") {
            spec_shape <- c(16, 17, 18, 15)
            spec_shape_size <- c(5.0, 4.5, 6.75, 4.5)
            legend_title <- element_blank()
            col_guide <- "legend"
        }

        if((spec == "all") && (data_norm == "inter-organ") && (data_scale == "scaled")) {
            legend_col <- 1
            legend_pos <- c(0.154, 0.8728)
            x_lim <- NULL
            y_title_mrg <- margin(t = 0, r = 8.0, b = 0, l = 7.0)
            plot_margin <- unit(c(0.55, 0.55, 0.5, 0.9),"cm")
        } else if((spec == "all") && (data_norm == "inter-organ") && (data_scale == "raw")) {
            legend_col <- 1
            legend_pos <- c(0.154, 0.8728)
            x_lim <- NULL
            y_title_mrg <- margin(t = 0, r = 3.5, b = 0, l = 6.75)
            plot_margin <- unit(c(0.55, 0.55, 0.5, 0.9),"cm")
        } else if((spec == "Brassicaceae") && (data_norm == "inter-organ") && (data_scale == "scaled")) {
            legend_col <- 1
            legend_pos <- "none"
            legend_spacing <- 0
            order_guide <- 0
            x_lim <- NULL
            y_title_mrg <- margin(t = 0, r = -3.5, b = 0, l = 18.0)
            plot_margin <- unit(c(0.55, 0.55, 0.5, 0.5),"cm")
        } else if((spec == "Brassicaceae") && (data_norm == "inter-organ") && (data_scale == "raw")) {
            legend_col <- 1
            legend_pos <- c(0.849, 0.505)
            legend_spacing <- 5.2
            order_guide <- 1
            x_lim <- NULL
            y_title_mrg <- margin(t = 0, r = -3.5, b = 0, l = 18.0)
            plot_margin <- unit(c(0.55, 0.55, 0.5, 0.5),"cm")
        }

        if (((devseq_spec == "Brassicaceae") && (set == "pc2_3")) || 
            (devseq_spec == "all")) {

            legend_pos <- "none"
        }

        x_lab <- paste(x_coord, pc_var1, "%)", sep="")
        y_lab <- paste(y_coord, pc_var2, "%)", sep="")

        plot <- ggplot(data, aes(x = PC1, y = PC2, colour=Organ, fill = Organ)) + 
        stat_bag(prop = 0.95, size=1.5) + 
        geom_point(aes(shape=Species, color=Organ, size=Species, stroke=2.25)) + 
        scale_shape_manual(values=spec_shape)  + 
        scale_x_continuous(expand = c(0.05, 0), limits = x_lim) + 
        scale_y_continuous(expand = c(0.05, 0)) +
        guides(shape = guide_legend(override.aes = list(size = spec_shape_size, stroke=2.75), 
               order = order_guide, ncol = legend_col), size = FALSE) + 
        guides(colour = guide_legend(override.aes = list(size=5, linetype = "blank", alpha=1))) + 
        scale_size_manual(values=spec_shape_size) + 
        # shapes = filled round, filled rect, empty square, filled square_rot, filled square, empty rect, inverted empty rect
        # colors = dark moderate violet, soft blue, dark green, moderate green, vivid yellow, orange, vivid red, soft pink
        scale_color_manual(values=c('#6a54a9','#53b0db', '#2c8654', '#96ba37','#fad819', '#f2a72f', '#ee412e', '#e075af'), 
            guide = col_guide) + 
        scale_fill_manual(values=c('#6a54a9','#53b0db', '#2c8654', '#96ba37','#fad819', '#f2a72f', '#ee412e', '#e075af'), 
            guide = col_guide) + 
        labs(x = x_lab, y = y_lab) + 
        theme(panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(), 
            panel.border = element_rect(colour = "black", fill=NA, size=1.75), 
            panel.background = element_blank(), 
            axis.title.y = element_text(size=24, margin = y_title_mrg, colour="black"), 
            axis.title.x = element_text(size=24, margin = margin(t = 12.75, r = 0, b = 4.5, l = 0), colour="black"), 
            axis.text.x = element_text(size=21.25, angle=0, margin = margin(t = 5), colour="black"), 
            axis.text.y = element_text(size=21.25, angle=0, margin = margin(r = 5), colour="black"), 
            axis.ticks.length=unit(0.35, "cm"), 
            axis.ticks = element_line(colour = "black", size = 0.7), 
            legend.key = element_rect(colour = "transparent", fill = "white"), 
            legend.key.size = unit(1.65,"line"), # default is 1.2
            legend.text = element_text(size=21.25), 
            legend.background = element_rect(fill = "transparent"),
            legend.title = legend_title,
            legend.spacing.y = unit(legend_spacing,'cm'), 
            legend.position = legend_pos,
            plot.margin = plot_margin)

        ggsave(file = file.path(out_dir, "output", "plots", fname), plot = plot,
            width = 8.35, height = 8, dpi = 300, units = c("in"), 
            limitsize = FALSE)
    }


    if (is.element("Brassicaceae", devseq_spec) && (data_scale == "raw")) {
        plotPCA(data=DevSeq_pca_1_2_w_stamen, pc_var1=DevSeq_pc1_var_w_stamen, pc_var2=DevSeq_pc2_var_w_stamen, 
            set="pc1_2", spec="Brassicaceae", data_norm=data_norm) 

        plotPCA(data=DevSeq_pca_2_3_w_stamen, pc_var1=DevSeq_pc2_var_w_stamen, pc_var2=DevSeq_pc3_var_w_stamen, 
            set="pc2_3", spec="Brassicaceae", data_norm=data_norm) 

        plotPCA(data=DevSeq_pca_1_3_w_stamen, pc_var1=DevSeq_pc1_var_w_stamen, pc_var2=DevSeq_pc3_var_w_stamen, 
            set="pc1_3", spec="Brassicaceae", data_norm=data_norm) 

    } else if (is.element("Brassicaceae", devseq_spec) && (data_scale == "scaled")) {
        plotPCA(data=DevSeq_pca_1_2_w_stamen, pc_var1=DevSeq_pc1_var_w_stamen, pc_var2=DevSeq_pc2_var_w_stamen, 
            set="pc1_2", spec="Brassicaceae", data_norm=data_norm) 

    } else if (is.element("all", devseq_spec) && (data_scale == "raw")) { 
        plotPCA(data=DevSeq_pca_1_2_w_stamen, pc_var1=DevSeq_pc1_var_w_stamen, pc_var2=DevSeq_pc2_var_w_stamen, 
            set="pc1_2", spec="all", data_norm=data_norm) 

        plotPCA(data=DevSeq_pca_2_3_w_stamen, pc_var1=DevSeq_pc2_var_w_stamen, pc_var2=DevSeq_pc3_var_w_stamen, 
            set="pc2_3", spec="all", data_norm=data_norm) 

        plotPCA(data=DevSeq_pca_1_3_w_stamen, pc_var1=DevSeq_pc1_var_w_stamen, pc_var2=DevSeq_pc3_var_w_stamen, 
            set="pc1_3", spec="all", data_norm=data_norm) 

    } else if (is.element("all", devseq_spec) && (data_scale == "scaled")) {
        plotPCA(data=DevSeq_pca_1_2_w_stamen, pc_var1=DevSeq_pc1_var_w_stamen, pc_var2=DevSeq_pc2_var_w_stamen, 
            set="pc1_2", spec="all", data_norm=data_norm) 
    }

  }


  if (devseq_spec == "Brassicaceae") {

    performPCA(data = x_df, data_scale = "raw", ntop = 10000)
    performPCA(data = x_df, data_scale = "scaled", ntop = nrow(x_df))
  
  } else {

    performPCA(data = x_df, data_scale = "raw", ntop = nrow(x_df))
    performPCA(data = x_df, data_scale = "scaled", ntop = nrow(x_df))

  }

 }

}


makeCompAnylsis(dataset="DevSeq", expr_estimation="TPM", coefficient="pearson", devseq_spec="all", data_norm="inter-organ", devseq_organs="all")
makeCompAnylsis(dataset="DevSeq", expr_estimation="counts", coefficient="pearson", devseq_spec="Brassicaceae", data_norm="inter-organ", devseq_organs="all")
makeCompAnylsis(dataset="DevSeq", expr_estimation="counts", coefficient="pearson", devseq_spec="all", data_norm="inter-organ", devseq_organs="all")
makeCompAnylsis(dataset="DevSeq", expr_estimation="counts", coefficient="pearson", devseq_spec="all", data_norm="inter-organ", devseq_organs="subset")
makeCompAnylsis(dataset="Brawand", expr_estimation="counts", coefficient="pearson", data_norm="inter-organ")





