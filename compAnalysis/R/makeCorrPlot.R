# Prepare Brawand and DevSeq comparative expression data
# Thresholds: 0.5 TPM (since there are no ERCC spike-ins in Brawand data)
# Data input: Brawand and DevSeq TPM expression tables of all samples


#------------------- Load packages, set directories and read sample tables ---------------------


# Install and load packages
if (!require(dplyr)) install.packages('dplyr')
library(dplyr)
if (!require(data.table)) install.packages('data.table')
library(data.table)
if (!require(plyr)) install.packages('plyr')
library(plyr)
if (!require(corrplot)) install.packages('corrplot')
library(corrplot)
if (!require(gplots)) install.packages('gplots')
library(gplots)
if (!require(factoextra)) install.packages('factoextra')
library(factoextra)
if (!require(dendextend)) install.packages('dendextend')
library(dendextend)


# Set file path and input files
in_dir <- "/Volumes/User/Shared/Christoph_manuscript/DevSeq_paper/Analysis/Analysis_2019/A_thaliana_gene_exression_map/20200527_CS_comparative/data"
out_dir <- "/Volumes/User/Shared/Christoph_manuscript/DevSeq_paper/Analysis/Analysis_2019/A_thaliana_gene_exression_map/20200527_CS_comparative"




# Read expression data
makeCorrPlot <- function(dataset = c("Brawand", "DevSeq"), expr_estimation = c("TPM", "counts"), 
	coefficient = c("pearson", "spearman")) {
	
	# Show error message if no species is chosen
    if (missing(dataset))
   
       stop(
       "Please choose one of the available data sets: 
	   'Brawand', 'DevSeq'",
	   call. = TRUE
       )

    # Show error message if unknown data set is chosen 
    if (!is.element(dataset, c("Brawand", "DevSeq")))
   
       stop(
       "Please choose one of the available data sets: 
	   'Brawand', 'DevSeq'",
	   call. = TRUE
       )

   	# Show error message if no species is chosen
    if (missing(expr_estimation))
   
       stop(
       "Please choose one of the available expression estimations: 
	   'TPM', 'counts'",
	   call. = TRUE
       )

    # Show error message if unknown expression estimation is chosen
    if (!is.element(expr_estimation, c("TPM", "counts")))
   
       stop(
       "Please choose one of the available expression estimations: 
	   'TPM', 'counts'",
	   call. = TRUE
       )

    # Show error message if no correlation coefficient is chosen
    if (missing(coefficient))
   
       stop(
       "Please choose one of the available correlation coefficients: 
	   'pearson', 'spearman'",
	   call. = TRUE
       )

    # Show error message if unknown correlation coefficient is chosen
    if (!is.element(coefficient, c("pearson", "spearman")))
   
       stop(
       "Please choose one of the available correlation coefficients: 
	   'pearson', 'spearman'",
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

    } else if ((is.element("DevSeq", dataset)) && (is.element("TPM", expr_estimation))) {
		genesExpr = file.path(in_dir, "Expression_data", "protein_gene_tpm_DESeq2_norm.csv")
		dataset_id <- "DevSeq"

    } else if ((is.element("DevSeq", dataset)) && (is.element("counts", expr_estimation))) {
		genesExpr = file.path(in_dir, "Expression_data", "rlog_counts_DESeq_norm.csv")
		dataset_id <- "DevSeq"
    }


  # Define simplified Brawand and DevSeq column names
  if (is.element("DevSeq", dataset)) {
    col_names <- rep(c("root", "hypocotyl", "leaf", "veg_apex", "inf_apex", 
      "flower", "carpel", "stamen"), each=3)
    replicate_tag_samples <- rep(c(".1",".2",".3"), times=8)
	  col_names <- paste0(col_names,replicate_tag_samples)
	  col_names <- rep(col_names, times=7)
    spec_names <- rep(c("_AT", "_AL", "_CR", "_ES", "_TH", 
      "_MT", "_BD"), each=24)
    col_names <- paste0(col_names, spec_names)
	}


	# Read expression data
	if (is.element("DevSeq", dataset)) {
	x <- read.table(genesExpr, sep=",", dec=".", skip = 1, header=FALSE, stringsAsFactors=FALSE)
	colnames(x)[1] <- "gene_id"
	} else if (is.element("Brawand", dataset)) {
	x <- read.table(genesExpr, sep=",", dec=".", header=TRUE, stringsAsFactors=FALSE)
	colnames(x)[1] <- "gene_id"
	}


    # Stop function here to allow specific analysis of a single data set
    # return_list <- list("dataset_id" = dataset_id, "expr_estimation" = expr_estimation, "x" = x, "coefficient" = coefficient, "col_names" = col_names)
    # return(return_list)
    # }
    # return_objects <- makeCorrPlot(dataset="DevSeq", expr_estimation="TPM", coefficient="pearson") # read in DevSeq expression data
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




#----------------------------- Calculate correlation and make plot -----------------------------


    # Create "plots" folder in /out_dir/output/plots
    if (!dir.exists(file.path(out_dir, "output", "plots"))) 
        dir.create(file.path(out_dir, "output", "plots"), recursive = TRUE)


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
    steps <- c("#600303", "#a80606", "#d20808", "yellow", "lightgoldenrodyellow")
    pal <- color.palette(steps, c(20, 20, 33, 4), space = "rgb")

    # Set filename
    fname <- sprintf('%s.png', paste(dataset_id, expr_estimation, coefficient, sep="_"))


    # Define column and row colors for color bars based on sample names and experiment
    # Build distance matrix & dendrogram then get dendrogram leaf colors to create color vector
    if (dataset_id == "DevSeq") {

      exp_col <- c(AT="gray1", AL="cornsilk3", CR="floralwhite", ES="wheat4", TH="#f1da5b", 
      	MT="#9f0000", BD="darkorchid4")

      species_col <- c(roo="#52428c", hyp="#808dc2", mes="#808dc2", lea="#00994f", veg="#95b73a", 
        inf="#eed410", spi="#eed410", flo="#de6daf", sta="#f23d29", mat="#a63126", car="#f2a529") # last two letters of experiment string

    } else if (dataset_id == "Brawand") {

      species_col <- c(se="palegreen4", fi="palegreen4", ap="orange2", hy="darkorchid4", fl="gold1", ro="orchid4", 
              le="olivedrab3", co="olivedrab3", ca="olivedrab3")

      organ_col <- c(GE="navajowhite2", eq="maroon4") # last two letters of experiment string
    }


    x[is.na(x)] <- 0 # replaces NAs by 0

    # Compute correlation and build distance matrix
    if (is.element("pearson", coefficient) && is.element("counts", expr_estimation)) {
        x <- cor(x[, 2:ncol(x)], method = "pearson")
        x_dist <- get_dist(x, stand = FALSE, method = "pearson")
        # correlation matrix does not need to be transposed for get_dist method (since ncol=nrow)

    } else if (is.element("pearson", coefficient) && is.element("TPM", expr_estimation)) {
        x[,2:ncol(x)] <- log2(x[,2:ncol(x)] + 1)
        x <- cor(x[, 2:ncol(x)], method = "pearson")
        x_dist <- get_dist(x, stand = FALSE, method = "pearson")

    } else if (is.element("spearman", coefficient)) {
        x <- cor(x[, 2:ncol(x)], method = "spearman")
        x_dist <- get_dist(x, stand = FALSE, method = "spearman")
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
  
    # Make corrplots
    png(height = 5000, width = 5000, pointsize = 20, file = file.path(out_dir, "output", "plots", fname))
    par(lwd = 15) # dendrogram line width
    getRowOrder = heatmap.2(x,
         revC = F,
         ColSideColors = col_cols, 
         RowSideColors = row_cols,
         distfun = function(c) get_dist(x, stand = FALSE, method = "pearson"), 
         hclustfun = function(x) hclust(x_dist, method = "average"))

    # Get order of rows and rearrange "row_cols" vector
    # fixes gplots heatmap.2 RowSideColors bug (colorbar does not reverse when revC=T)
    ordinary_order = getRowOrder$rowInd
    reversal = cbind(ordinary_order, rev(ordinary_order))
    rev_col = row_cols[reversal[,2]]; rev_col = rev_col[order(reversal[,1])];

    # Create heatmap with reversed RowSideColors
    heatmap.2(x, 
         revC = T,
         ColSideColors = col_cols, 
         RowSideColors = rev_col, 
         density.info = "none",
         trace = "none",
         col = pal(800),
         Colv=TRUE, 
         cexRow = 2,
         cexCol = 2,
         margins = c(60, 60),
         key.par = list(cex = 2.75),
         lwid = c(0.5,4,17.5), # column width
         lhei = c(0.5,4,17.5), # column height
         offsetRow = 1,
         offsetCol = 1,
         key.xlab = NA,
         key.title = NULL,
         distfun = function(c) get_dist(x, stand = FALSE, method = "pearson"), 
         hclustfun = function(x) hclust(x_dist, method = "average"))

    dev.off()
}

makeCorrPlot(dataset="DevSeq", expr_estimation="counts", coefficient="pearson")
makeCorrPlot(dataset="DevSeq", expr_estimation="counts", coefficient="spearman")

makeCorrPlot(dataset="DevSeq", expr_estimation="TPM", coefficient="pearson")
makeCorrPlot(dataset="DevSeq", expr_estimation="TPM", coefficient="spearman")










#---------- Check clustering results of correlation plot with other implementations #-----------


# Make corrplots with corrplot function
# See SO https://stackoverflow.com/questions/45896231/r-corrplot-with-clustering-default-dissimilarity-measure-for-correlation-matrix
    png(height = 5000, width = 5000, pointsize = 20, file = file.path(out_dir, "output", "plots", "DevSeq_TPM_hclust_avg_corrplot.png"))
    par(lwd = 15) # dendrogram line width
    corrplot(x,
         method = "color", 
         order = "hclust", 
         hclust.method = "average", 
         col=pal(100), 
         addrect = NULL, 
         rect.lwd = 0,
         number.cex = 4.25
         )
    dev.off()


# Make corrplots with base R heatmap function
# This gives exact same reults as using heatmap.2 function
    png(height = 5000, width = 5000, pointsize = 20, file = file.path(out_dir, "output", "plots", "DevSeq_TPM_pearson_dist_baseR_heatmap.png"))
    par(lwd = 15) # dendrogram line width
    heatmap(x,
         col = pal(200),
         margins = c(20, 20),
         key.title = NULL,
         distfun = function(c) get_dist(x, method = "pearson"), 
         hclustfun = function(x) hclust(x, method = "average")
         )
    dev.off()


