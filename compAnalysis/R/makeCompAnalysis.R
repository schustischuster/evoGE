# Prepare Brawand and DevSeq comparative expression data
# Thresholds: 0.5 TPM (since there are no ERCC spike-ins in Brawand data)
# Data input: Brawand and DevSeq TPM expression tables of all samples



#-------------------------------------- Read data tables ---------------------------------------


makeCompAnylsis <- function(dataset = c("Brawand", "DevSeq"), expr_estimation = c("TPM", "counts"), 
	coefficient = c("pearson", "spearman"), devseq_spec = c("Brassicaceae", "all")) {
	
	# Show error message if no species or unknown data set is chosen
    if ((missing(dataset)) | (!is.element(dataset, c("Brawand", "DevSeq"))))
   
       stop(
       "Please choose one of the available data sets: 
	   'Brawand', 'DevSeq'",
	   call. = TRUE
       )

   	# Show error message if expression estimation or unknown expression estimation is chosen
    if ((missing(expr_estimation)) | (!is.element(expr_estimation, c("TPM", "counts"))))
   
       stop(
       "Please choose one of the available expression estimations: 
	   'TPM', 'counts'",
	   call. = TRUE
       )

    # Show error message if no correlation or unknown correlation coefficient is chosen
    if ((missing(coefficient)) | (!is.element(coefficient, c("pearson", "spearman"))))
   
       stop(
       "Please choose one of the available correlation coefficients: 
	   'pearson', 'spearman'",
	   call. = TRUE
       )

    # Show error message if no devseq_spec or unknown devseq_spec is chosen
    if ((is.element("DevSeq", dataset)) && ((missing(devseq_spec)) | (!is.element(devseq_spec, c("Brassicaceae", "all")))))
   
       stop(
       "Please choose one of the available DevSeq species sets: 
	   'Brassicaceae', 'all'",
	   call. = TRUE
       )


    # Show startup message
    message("Reading data...")


	# Set expression input file
    if ((is.element("Brawand", dataset)) && (is.element("TPM", expr_estimation))) {
        genesExpr = file.path(in_dir, "Expression_data", "TPM_Brawand_norm.csv")
        dataset_id <- "Brawand"

    } else if ((is.element("Brawand", dataset)) && (is.element("counts", expr_estimation))) {
        genesExpr = file.path(in_dir, "Expression_data", "rlog_Brawand_norm.csv")
        dataset_id <- "Brawand"

    } else if ((is.element("DevSeq", dataset)) && (is.element("TPM", expr_estimation)) 
    	&& (is.element("Brassicaceae", devseq_spec))) {
		genesExpr = file.path(in_dir, "Expression_data", "protein_gene_tpm_DESeq2_norm_Brassicaceae.csv")
		dataset_id <- "DevSeq"

	} else if ((is.element("DevSeq", dataset)) && (is.element("TPM", expr_estimation)) 
    	&& (is.element("all", devseq_spec))) {
		genesExpr = file.path(in_dir, "Expression_data", "protein_gene_tpm_DESeq2_norm.csv")
		dataset_id <- "DevSeq"

    } else if ((is.element("DevSeq", dataset)) && (is.element("counts", expr_estimation)) 
    	&& (is.element("Brassicaceae", devseq_spec))) {
		genesExpr = file.path(in_dir, "Expression_data", "rlog_counts_DESeq_norm_Brassicaceae.csv")
		dataset_id <- "DevSeq"

    } else if ((is.element("DevSeq", dataset)) && (is.element("counts", expr_estimation)) 
    	&& (is.element("all", devseq_spec))) {
		genesExpr = file.path(in_dir, "Expression_data", "rlog_counts_DESeq_norm.csv")
		dataset_id <- "DevSeq"
    }


	# Define simplified Brawand and DevSeq column names
	if (is.element("DevSeq", dataset) && (is.element("all", devseq_spec))) {
		col_names <- rep(c("root", "hypocotyl", "leaf", "veg_apex", "inf_apex", 
			"flower", "carpel", "stamen"), each=3)
		replicate_tag_samples <- rep(c(".1",".2",".3"), times=8)
		col_names <- paste0(col_names,replicate_tag_samples)
		col_names <- rep(col_names, times=7)
		spec_names <- rep(c("_AT", "_AL", "_CR", "_ES", "_TH", "_MT", "_BD"), each=24)
		col_names <- paste0(col_names, spec_names)

	} else if (is.element("DevSeq", dataset) && (is.element("Brassicaceae", devseq_spec))) {
		col_names <- rep(c("root", "hypocotyl", "leaf", "veg_apex", "inf_apex", 
			"flower", "carpel", "stamen"), each=3)
		replicate_tag_samples <- rep(c(".1",".2",".3"), times=8)
		col_names <- paste0(col_names,replicate_tag_samples)
		col_names <- rep(col_names, times=4)
		spec_names <- rep(c("_AT", "_AL", "_CR", "_ES"), each=24)
		col_names <- paste0(col_names, spec_names)
	}


	# Read expression data
	if (is.element("DevSeq", dataset)) {
		x <- read.table(genesExpr, sep=",", dec=".", skip = 1, header=FALSE, stringsAsFactors=FALSE)
		colnames(x)[1] <- "gene_id"

	} else if (is.element("Brawand", dataset)) {
		x <- read.table(genesExpr, sep=",", dec=".", header=TRUE, stringsAsFactors=FALSE)
		colnames(x)[1] <- "gene_id"
		Br_spec <- substring(names(x)[-1], 1, 3)
		Br_org <- substring(names(x)[-1], 5)
		colnames(x)[-1] <- paste0(Br_org, "_", Br_spec)
	}


    # Stop function here to allow specific analysis of a single data set
    # For DevSeq
    # return_list <- list("dataset_id" = dataset_id, "expr_estimation" = expr_estimation, "x" = x, "coefficient" = coefficient, "col_names" = col_names, "devseq_spec" = devseq_spec)
    # For Brawand
    # return_list <- list("dataset_id" = dataset_id, "expr_estimation" = expr_estimation, "x" = x, "coefficient" = coefficient)
    # return(return_list)
    # }
    # return_objects <- makeCompAnylsis(dataset="DevSeq", expr_estimation="TPM", coefficient="pearson", devseq_spec="all") # read in DevSeq expression data
    # return_objects <- makeCompAnylsis(dataset="Brawand", expr_estimation="counts", coefficient="spearman") # read in Brawand expression data
    # list2env(return_objects, envir = .GlobalEnv)

    # remove additional A.lyrata stamens samples from DevSeq data and update column names
    if ((dataset_id == "DevSeq") && is.element("counts", expr_estimation)) {

    	x <- x[,-(50:61),drop=FALSE]
    	x <- x[c(1:19,23:25,20:22,26:ncol(x))] # ATH stamen and carpel samples in inverted compared to other species

    	# set column names
    	colnames(x)[2:ncol(x)] <- col_names

    } else if ((dataset_id == "DevSeq") && is.element("TPM", expr_estimation)) {

        x <- x[c(1:19,23:25,20:22,26:ncol(x))] # ATH stamen and carpel samples in inverted compared to other species

        # set column names
        colnames(x)[2:ncol(x)] <- col_names
    }




#--------------------- Calculate correlation and prepare data for corrplot ---------------------


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
    steps <- c("#c42d2d", "#cf3c1f", "#faa11b", "#fff415", "#fcfce2")

    if (dataset_id == "DevSeq") {
        pal <- color.palette(steps, c(25, 25, 25, 5), space = "rgb")
    } else if (dataset_id == "Brawand") {
        pal <- color.palette(steps, c(25, 25, 25, 15), space = "rgb")
    }

    # Set filename
    fname <- sprintf('%s.png', paste(dataset_id, expr_estimation, coefficient, sep="_"))


    # Define column and row colors for color bars based on sample names and experiment
    # Build distance matrix & dendrogram then get dendrogram leaf colors to create color vector
    if (dataset_id == "DevSeq") {

        exp_col <- c(AT="gray1", AL="cornsilk3", CR="floralwhite", ES="wheat4", TH="lightgoldenrod2", 
            MT="#9f0000", BD="#7c3979") # last two letters of sample name

        species_col <- c(roo="#52428c", hyp="#8591c7", mes="#8591c7", lea="#00994f", veg="#95b73a", 
            inf="#fad819", spi="#fad819", flo="#de6daf", sta="#f23d29", mat="#a63126", car="#f2a529")

    } else if (dataset_id == "Brawand") {

        exp_col <- c(sa="gray1", tr="cornsilk3", pa="floralwhite", go="wheat4", py="lightgoldenrod2", 
            ml="#9f0000", mu="#7c3979", do="green") # last two letters of sample name

        species_col <- c(br_="#42448c", cb_="#6a76ad", ht_="gold1", kd_="olivedrab3", lv_="palegreen4", 
            ts_="#e17a21")
    }


    x[is.na(x)] <- 0 # replaces NAs by 0

    # Log2 transform data if expr_estimation=TPM is chosen
    # Compute correlation and build distance matrix
    if (is.element("pearson", coefficient) && is.element("counts", expr_estimation)) {
        x_cor <- cor(x[, 2:ncol(x)], method = "pearson")
        x_dist <- get_dist(x_cor, stand = FALSE, method = "pearson")
        # correlation matrix does not need to be transposed for get_dist method (since ncol=nrow)

    } else if (is.element("pearson", coefficient) && is.element("TPM", expr_estimation)) {
        x[,2:ncol(x)] <- log2(x[,2:ncol(x)] + 1)
        x_cor <- cor(x[, 2:ncol(x)], method = "pearson")
        x_dist <- get_dist(x_cor, stand = FALSE, method = "pearson")

    } else if (is.element("spearman", coefficient)) {
        x_cor <- cor(x[, 2:ncol(x)], method = "spearman")
        x_dist <- get_dist(x_cor, stand = FALSE, method = "spearman")
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




#---------------------------- Make corrplot for DevSeq and Brawand -----------------------------



	if ((dataset_id == "DevSeq") && (is.element("all", devseq_spec))) {

        # Make corrplots
        png(height = 3500, width = 3500, pointsize = 20, file = file.path(out_dir, "output", "plots", fname))
        par(lwd = 17) # dendrogram line width
        getRowOrder = heatmap.2(x_cor,
            revC = F,
            ColSideColors = col_cols, 
            RowSideColors = row_cols,
            distfun = function(c) get_dist(x_cor, stand = FALSE, method = "pearson"), 
            hclustfun = function(x) hclust(x_dist, method = "average"))

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
            Colv=TRUE, 
            cexRow = 2,
            cexCol = 2,
            margins = c(35, 35),
            key.par = list(cex = 2.75),
            lwid = c(0.3,2.5,17.5), # column width
            lhei = c(0.3,2.5,17.5), # column height
            offsetRow = 1,
            offsetCol = 1,
            key.xlab = NA,
            key.title = NULL,
            distfun = function(c) get_dist(x_cor, stand = FALSE, method = "pearson"), 
            hclustfun = function(x) hclust(x_dist, method = "average"))
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

    } else if (dataset_id == "Brawand") {

        # Make corrplots
        png(height = 4000, width = 4000, pointsize = 20, file = file.path(out_dir, "output", "plots", fname))
        par(lwd = 19) # dendrogram line width
        getRowOrder = heatmap.2(x_cor,
            revC = F,
            ColSideColors = col_cols, 
            RowSideColors = row_cols,
            distfun = function(c) get_dist(x_cor, stand = FALSE, method = "pearson"), 
            hclustfun = function(x) hclust(x_dist, method = "average"))

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
            Colv=TRUE, 
            cexRow = 2,
            cexCol = 2,
            margins = c(30, 30),
            key.par = list(cex = 2.75),
            lwid = c(0.35,2.75,17.5), # column width
            lhei = c(0.35,2.75,17.5), # column height
            offsetRow = 1,
            offsetCol = 1,
            key.xlab = NA,
            key.title = NULL,
            distfun = function(c) get_dist(x_cor, stand = FALSE, method = "pearson"), 
            hclustfun = function(x) hclust(x_dist, method = "average"))
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


    x_avg <- calculateAvgExpr(x)

    if ((dataset_id == "DevSeq") && (is.element("all", devseq_spec))) {

        DevSeq_col_names <- rep(c("root", "hypocotyl", "leaf", "veg_apex", "inf_apex", 
            "flower", "carpel", "stamen"), times=7)
        DevSeq_spec_names <- rep(c("_AT", "_AL", "_CR", "_ES", "_TH", 
            "_MT", "_BD"), each=8)
        repl_names <- paste0(DevSeq_col_names, DevSeq_spec_names)

        colnames(x_avg)[2:ncol(x_avg)] <- repl_names

    } else if ((dataset_id == "DevSeq") && (is.element("Brassicaceae", devseq_spec))) {

        DevSeq_col_names <- rep(c("root", "hypocotyl", "leaf", "veg_apex", "inf_apex", 
            "flower", "carpel", "stamen"), times=4)
        DevSeq_spec_names <- rep(c("_AT", "_AL", "_CR", "_ES"), each=8)
        repl_names <- paste0(DevSeq_col_names, DevSeq_spec_names)

        colnames(x_avg)[2:ncol(x_avg)] <- repl_names
    }




#---------- Get replicate correlations and generate DevSeq inter-organ distance plot ----------


    if ((dataset_id == "DevSeq") && (is.element("all", devseq_spec))) {

        x_cor_avg <- cor(x_avg[2:ncol(x_avg)]) 

        AT <- as.data.frame(unique(as.vector(x_cor_avg[1:8,1:8])))[-1,]
        AL <- as.data.frame(unique(as.vector(x_cor_avg[9:16,9:16])))[-1,]
        CR <- as.data.frame(unique(as.vector(x_cor_avg[17:24,17:24])))[-1,]
        ES <- as.data.frame(unique(as.vector(x_cor_avg[25:32,25:32])))[-1,]
        TH <- as.data.frame(unique(as.vector(x_cor_avg[33:40,33:40])))[-1,]
        MT <- as.data.frame(unique(as.vector(x_cor_avg[41:48,41:48])))[-1,]
        BD <- as.data.frame(unique(as.vector(x_cor_avg[49:56,49:56])))[-1,]

        df_names <- c("Species" , "Correlation")

        species_DevSeq <- c("A.thaliana", "A.lyrata", "C.rubella", "E.salsugineum", 
            "T.hassleriana", "M.truncatula", "B.distachyon")
        species_names <- as.data.frame(rep(species_DevSeq, each=28))
        organ_cor <- as.data.frame(c(AT,AL,CR,ES,TH,MT,BD))

        cor_df <- cbind(species_names, organ_cor)
        colnames(cor_df) <- df_names

        # Cor values of meristematic tssues apex_inf apex_veg and carpel
        AT_m <- as.data.frame(unique(as.vector(x_cor_avg[c(4:5,7),c(4:5,7)])))[-1,]
        AL_m <- as.data.frame(unique(as.vector(x_cor_avg[c(12:13,15),c(12:13,15)])))[-1,]
        CR_m <- as.data.frame(unique(as.vector(x_cor_avg[c(20:21,23),c(20:21,23)])))[-1,]
        ES_m <- as.data.frame(unique(as.vector(x_cor_avg[c(28:29,31),c(28:29,31)])))[-1,]
        TH_m <- as.data.frame(unique(as.vector(x_cor_avg[c(36:37,39),c(36:37,39)])))[-1,]
        MT_m <- as.data.frame(unique(as.vector(x_cor_avg[c(44:45,47),c(44:45,47)])))[-1,]
        BD_m <- as.data.frame(unique(as.vector(x_cor_avg[c(52:53,55),c(52:53,55)])))[-1,]

        species_names_m <- as.data.frame(rep(species_DevSeq, each=3))
        organ_cor_m <- as.data.frame(c(AT_m,AL_m,CR_m,ES_m,TH_m,MT_m,BD_m))

        cor_df_m <- cbind(species_names_m, organ_cor_m)
        colnames(cor_df_m) <- df_names




        # Make inter-organ distance plot
        plotOrganDist <- function(data, data_m) {

            fname <- paste('DevSeq_inter_organ_dist_', coefficient, '.jpg', sep="")

            coef_lab <- ifelse(coefficient == "pearson", "Pearson's r", "Spearman's rho")

            if(coefficient == "pearson") {
            	y_min = 0.435
            	y_max = 1.025
            } else if(coefficient == "spearman") {
            	y_min = 0
            	y_max = 1.025
            }

            cor_colors <- rep(c("#808080","#c99407","#4fa722","#1c9c6f","#2690c1","#795fcf","#c91e1b"), each=28)

            p <- ggplot(data, aes(x=factor(Species, levels = Species), y=Correlation, fill=Species)) + 
            stat_boxplot(geom ='errorbar', width = 0.45, size=1.0, color="gray15") + 
            geom_boxplot(width = 0.75, size=1.0, color="gray15", outlier.shape = 21, 
                outlier.size = 2.5, outlier.stroke = 1.5, outlier.fill = NA, outlier.color="gray35") + 
            geom_beeswarm(data = data, priority = c("ascending"), colour=cor_colors, cex=2, size=3) +
            geom_beeswarm(data = data_m, priority = c("descending"), colour="black", cex=2, size=5, shape=18, fill="black") + 
            scale_y_continuous(limits = c(y_min, y_max), expand = c(0, 0)) + 
            annotate("rect", xmin=0.35, xmax=7.65, ymin=y_min, ymax=y_max, fill="white", alpha=0, 
                color="gray15", size=1.35)

            q <- p + scale_fill_manual(values=c("#f5d88c","#cacaca","#efb0ae","#add89a","#8adcc0","#d0c7f1","#9ed0e7")) + 
            theme_minimal() + 
            xlab("") + 
            ylab(coef_lab) + 
            theme(legend.position = "none", 
                text=element_text(size=20.5), 
                panel.grid.major.y = element_line(size = 0.7, color = c("gray85")),
                panel.grid.minor.y = element_line(size = 0.5, color = c("gray85")), 
                panel.grid.major.x = element_line(size = 0.5, color = c("gray85")), 
                axis.ticks.length = unit(.3, "cm"),
                axis.ticks = element_line(colour = "gray15", size = 0.7), 
                axis.title.x = element_text(colour = "black", size=21, 
                    margin = margin(t = 17.5, r = 0, b = 0, l = 0)), 
                axis.title.y = element_text(colour = "black", size=21, 
                    margin = margin(t = 0, r = 12.5, b = 0, l = 0.5)), 
                axis.text.x = element_text(colour = "black", size=18.5, angle=90, 
                    margin = margin(t = 7.0, r = 0, b = 0, l = 0), hjust=1.0, vjust=0.35),
                axis.text.y = element_text(colour = "black", margin = margin(t = 0, r = 5, b = 0, l = 0)),  
                plot.margin = unit(c(35, 35, 0, 35), "points"))

            ggsave(file = file.path(out_dir, "output", "plots", fname), plot = q,
                scale = 1, width = 9.1, height = 8.0, units = c("in"), 
                dpi = 600, limitsize = FALSE)
        }

        plotOrganDist(data=cor_df, data_m=cor_df_m) 

    }




#------------------------------------------ Make PCAs ------------------------------------------


    makePCA <- function(df, pca_dim = c("1_2", "2_3")) {

        # Perform PCA analysis on log-transformed data (if TPM chosen) or rlog counts
        comp_pca <- prcomp(x[,2:ncol(x)], center = TRUE, scale = TRUE) 

        # Get eigenvalues, explained variance (%) and cumulative variance (%) 
	    eig_val <- get_eigenvalue(comp_pca)
	    pc1_var <- round((eig_val$variance.percent[1]), digits = 1) # variance explained by PC1
	    pc2_var <- round((eig_val$variance.percent[2]), digits = 1) # variance explained by PC2
	    pc3_var <- round((eig_val$variance.percent[3]), digits = 1) # variance explained by PC3

	    # Results for Variables
	    res_var <- get_pca_var(comp_pca)          # Show output results
	    pca_coord <- as.data.frame(res_var$coord) # Get coordinates as dataframe

	    if (is.element("1_2", pca_dim)) {
	        pca_coord <- pca_coord[,1:2]

	    } else if (is.element("2_3", pca_dim)) {
	        pca_coord <- pca_coord[,2:3]
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
        colnames(pca_df) <- c("PC1", "PC2", "Species", "Tissue")

        return_list <- list("pca_df"=pca_df, "pc1_var"=pc1_var, "pc2_var"=pc2_var, "pc3_var"=pc3_var, "eig_val"=eig_val)
        return(return_list)
    }


    # Make PCA for DevSeq w/ stamen data PC1/2
    pca_return_objects <- makePCA(df = x, pca_dim = "1_2")
    list2env(pca_return_objects, envir = .GlobalEnv)
    DevSeq_pca_1_2_w_stamen <- pca_df
    DevSeq_pc1_var_w_stamen <- pc1_var
    DevSeq_pc2_var_w_stamen <- pc2_var
    DevSeq_pc3_var_w_stamen <- pc3_var
    # Make PCA for DevSeq w/ stamen data PC2/3
    pca_return_objects <- makePCA(df = x, pca_dim = "2_3")
    list2env(pca_return_objects, envir = .GlobalEnv)
    DevSeq_pca_2_3_w_stamen <- pca_df



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

#' @inheritParams ggplot2::stat_identity
#' @param prop Proportion of all the points to be included in the bag (default is 0.5)

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

#' @rdname ggplot2-ggproto
#' @format NULL
#' @usage NULL
#' @export

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
        DevSeq_pca_1_2_w_stamen$Tissue <- factor(DevSeq_pca_1_2_w_stamen$Tissue, c("root", "hypocotyl", 
        	"leaf", "veg_apex", "inf_apex", "carpel", "stamen", "flower"))
        DevSeq_pca_2_3_w_stamen$Tissue <- factor(DevSeq_pca_2_3_w_stamen$Tissue, c("root", "hypocotyl", 
        	"leaf", "veg_apex", "inf_apex", "carpel", "stamen", "flower"))


    # Make PCA plots for main figure
    plotPCA <- function(data, pc_var1, pc_var2, set=c("pc1_2","pc2_3"), spec=c("Brassicaceae", "all")) {

        fname <- sprintf('%s.png', paste(deparse(substitute(data)), "pca", expr_estimation, sep="_"))
        
        if(set == "pc1_2") {
        	x_coord <- "PC1 ("
        	y_coord <- "PC2 ("
        } else if(set == "pc2_3") {
        	x_coord <- "PC2 ("
        	y_coord <- "PC3 ("
        }

        x_lab <- paste(x_coord, pc_var1, "%)", sep="")
        y_lab <- paste(y_coord, pc_var2, "%)", sep="")

        plot <- ggplot(data, aes(x = PC1, y = PC2, colour=Tissue, fill = Tissue)) + 
        stat_bag(prop = 0.95, size=1.5) + 
        geom_point(aes(shape=Species, color=Tissue, size=Species, stroke=2.25)) + 
        scale_shape_manual(values=c(16, 17, 0, 18, 15, 2, 6))  + 
        theme(text = element_text(size=22.5)) + 
        guides(colour = guide_legend(override.aes = list(size=4))) + 
        guides(shape = guide_legend(override.aes = list(size = c(5.0, 4.5, 3.5, 6.75, 5.0, 3.15, 3.15)))) + 
        scale_size_manual(values=c(5.0, 4.5, 3.5, 6.75, 5.0, 3.15, 3.15)) + 
        # shapes = filled round, filled rect, empty square, filled square_rot, filled square, empty rect, inverted empty rect
        scale_color_manual(values=c('#5850a3','#8591c7', '#00994f', '#95b73a','#fad819', '#de6daf','#f2a529', '#f23d29')) + 
        scale_fill_manual(values=c('#5850a3','#8591c7', '#00994f', '#95b73a','#fad819', '#de6daf','#f2a529', '#f23d29')) + 
        labs(x = x_lab, y = y_lab) + 
        theme(axis.line = element_line(colour = "black"), 
        	panel.grid.major = element_blank(), 
        	panel.grid.minor = element_blank(), 
        	panel.border = element_rect(colour = "black", fill=NA, size=1), 
        	panel.background = element_blank(), 
        	axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 4)), 
        	axis.title.x = element_text(margin = margin(t = 12.75, r = 0, b = 4, l = 0)), 
        	axis.text.x = element_text(size=21.25, angle=0, margin = margin(t = 5)), 
        	axis.text.y = element_text(size=21.25, angle=0, margin = margin(r = 5)), 
        	axis.ticks.length=unit(0.35, "cm"), 
        	axis.ticks = element_line(colour = "black", size = 0.7), 
        	plot.margin=unit(c(0.5,1,0.5,0.5),"cm"))

        ggsave(file = file.path(out_dir, "output", "plots", fname), plot = plot,
            width = 10.25, height = 8, dpi = 300, units = c("in"), 
            limitsize = FALSE)
    }


    if ((dataset_id == "DevSeq") && (is.element("Brassicaceae", devseq_spec))) {
    	plotPCA(data=DevSeq_pca_1_2_w_stamen, pc_var1=DevSeq_pc1_var_w_stamen, pc_var2=DevSeq_pc2_var_w_stamen, 
    		set="pc1_2", spec="Brassicaeae") 

    } else if ((dataset_id == "DevSeq") && (is.element("all", devseq_spec))) { 
    	plotPCA(data=DevSeq_pca_1_2_w_stamen, pc_var1=DevSeq_pc1_var_w_stamen, pc_var2=DevSeq_pc2_var_w_stamen, 
    		set="pc1_2", spec="all") 
    }



	# Make PCA plots for supplement (different font/line/shape sizes)
    plotPCA <- function(data, pc_var1, pc_var2, set=c("pc1_2","pc2_3"), spec=c("Brassicaceae", "all")) {

        fname <- sprintf('%s.png', paste(deparse(substitute(data)), "pca", expr_estimation, sep="_"))
        
        if(set == "pc1_2") {
        	x_coord <- "PC1 ("
        	y_coord <- "PC2 ("
        } else if(set == "pc2_3") {
        	x_coord <- "PC2 ("
        	y_coord <- "PC3 ("
        }

        x_lab <- paste(x_coord, pc_var1, "%)", sep="")
        y_lab <- paste(y_coord, pc_var2, "%)", sep="")

        plot <- ggplot(data, aes(x = PC1, y = PC2, colour=Tissue, fill = Tissue)) + 
        stat_bag(prop = 0.95, size=1.5) + 
        geom_point(aes(shape=Species, color=Tissue, size=Species, stroke=2.25)) + 
        scale_shape_manual(values=c(16, 17, 0, 18, 15, 2, 6))  + 
        theme(text = element_text(size=22.5)) + 
        guides(colour = guide_legend(override.aes = list(size=4))) + 
        guides(shape = guide_legend(override.aes = list(size = c(5.0, 4.5, 3.5, 6.75, 5.0, 3.15, 3.15)))) + 
        scale_size_manual(values=c(5.0, 4.5, 3.5, 6.75, 5.0, 3.15, 3.15)) + 
        # shapes = filled round, filled rect, empty square, filled square_rot, filled square, empty rect, inverted empty rect
        scale_color_manual(values=c('#5850a3','#8591c7', '#00994f', '#95b73a','#fad819', '#f2a529', '#f23d29', '#de6daf')) + 
        scale_fill_manual(values=c('#5850a3','#8591c7', '#00994f', '#95b73a','#fad819', '#f2a529', '#f23d29', '#de6daf')) + 
        labs(x = x_lab, y = y_lab) + 
        theme(axis.line = element_line(colour = "black"), 
        	panel.grid.major = element_blank(), 
        	panel.grid.minor = element_blank(), 
        	panel.border = element_rect(colour = "black", fill=NA, size=1), 
        	panel.background = element_blank(), 
        	axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 4)), 
        	axis.title.x = element_text(margin = margin(t = 12.75, r = 0, b = 4, l = 0)), 
        	axis.text.x = element_text(size=21.25, angle=0, margin = margin(t = 5)), 
        	axis.text.y = element_text(size=21.25, angle=0, margin = margin(r = 5)), 
        	axis.ticks.length=unit(0.35, "cm"), 
        	axis.ticks = element_line(colour = "black", size = 0.7), 
        	plot.margin=unit(c(0.5,1,0.5,0.5),"cm"))

        ggsave(file = file.path(out_dir, "output", "plots", fname), plot = plot,
            width = 10.25, height = 8, dpi = 300, units = c("in"), 
            limitsize = FALSE)
    }


    if ((dataset_id == "DevSeq") && (is.element("Brassicaceae", devseq_spec))) {
    	plotPCA(data=DevSeq_pca_2_3_w_stamen, pc_var1=DevSeq_pc2_var_w_stamen, pc_var2=DevSeq_pc3_var_w_stamen, 
    		set="pc2_3", spec="Brassicaeae") 

    } else if ((dataset_id == "DevSeq") && (is.element("all", devseq_spec))) { 
    	plotPCA(data=DevSeq_pca_2_3_w_stamen, pc_var1=DevSeq_pc2_var_w_stamen, pc_var2=DevSeq_pc3_var_w_stamen, 
    		set="pc2_3", spec="all") 
    }










}


makeCompAnylsis(dataset="DevSeq", expr_estimation="counts", coefficient="pearson", spec="Brassicaeae")
makeCompAnylsis(dataset="DevSeq", expr_estimation="counts", coefficient="spearman", spec="Brassicaeae")

makeCompAnylsis(dataset="DevSeq", expr_estimation="counts", coefficient="pearson", spec="all")
makeCompAnylsis(dataset="DevSeq", expr_estimation="counts", coefficient="spearman", spec="all")

makeCompAnylsis(dataset="DevSeq", expr_estimation="TPM", coefficient="pearson", spec="Brassicaeae")
makeCompAnylsis(dataset="DevSeq", expr_estimation="TPM", coefficient="spearman", spec="Brassicaeae")

makeCompAnylsis(dataset="DevSeq", expr_estimation="TPM", coefficient="pearson", spec="all")
makeCompAnylsis(dataset="DevSeq", expr_estimation="TPM", coefficient="spearman", spec="all")

makeCompAnylsis(dataset="Brawand", expr_estimation="counts", coefficient="pearson")
makeCompAnylsis(dataset="Brawand", expr_estimation="counts", coefficient="spearman")







#---------- Check clustering results of correlation plot with other implementations #-----------

# if (!require(corrplot)) install.packages('corrplot')
# library(corrplot)

# Make corrplots with corrplot function
# See SO https://stackoverflow.com/questions/45896231/r-corrplot-with-clustering-default-dissimilarity-measure-for-correlation-matrix
   # png(height = 5000, width = 5000, pointsize = 20, file = file.path(out_dir, "output", "plots", "DevSeq_TPM_hclust_avg_corrplot.png"))
   # par(lwd = 15) # dendrogram line width
   # corrplot(x,
        # method = "color", 
        # order = "hclust", 
        # hclust.method = "average", 
        # col=pal(100), 
        # addrect = NULL, 
        # rect.lwd = 0,
        # number.cex = 4.25
        # )
  # dev.off()


# Make corrplots with base R heatmap function
# This gives exact same reults as using heatmap.2 function
   # png(height = 5000, width = 5000, pointsize = 20, file = file.path(out_dir, "output", "plots", "DevSeq_TPM_pearson_dist_baseR_heatmap.png"))
   # par(lwd = 15) # dendrogram line width
   # heatmap(x,
        # col = pal(200),
        # margins = c(20, 20),
        # key.title = NULL,
        # distfun = function(c) get_dist(x, method = "pearson"), 
        # hclustfun = function(x) hclust(x, method = "average")
        # )
   # dev.off()


