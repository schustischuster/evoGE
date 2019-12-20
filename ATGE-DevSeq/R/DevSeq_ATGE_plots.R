
# Vizualize the results of the DevSeq-ATGE comparative analysis using the relative expression
# tables generated in the previous step.


#---------------------------- Read sample tables and prepare data -----------------------------


# Set input files
devseq_re_vs_atge_re_input <- "DevSeq_RE_vs_ATGE_RE.csv"
devseq_log2_re_vs_atge_log2_re_input <- "DevSeq_log2_RE_vs_ATGE_log2_RE.csv"


# Read data tables
message("Reading ATGE-DevSeq cor data tables")
devseq_re_vs_atge_re <- read.table(file=file.path(in_dir, devseq_re_vs_atge_re_input), sep=";", dec=".", header=TRUE, stringsAsFactors = FALSE)
devseq_log2_re_vs_atge_log2_re <- read.table(file=file.path(in_dir, devseq_log2_re_vs_atge_log2_re_input), sep=";", dec=".", header=TRUE, stringsAsFactors = FALSE)


# Prepare data
# Add experiment information to data tables
experiment <- gsub(".*_","",devseq_re_vs_atge_re[,3])
devseq_re_vs_atge_re <- cbind(as.data.frame(experiment),devseq_re_vs_atge_re)
experiment_log <- gsub(".*_","",devseq_log2_re_vs_atge_log2_re[,3])
devseq_log2_re_vs_atge_log2_re <- cbind(as.data.frame(experiment),devseq_log2_re_vs_atge_log2_re)

# Seperate ATGE and DevSeq data
atge_re <- subset(devseq_re_vs_atge_re, experiment == "ATGE")
devseq_re <- subset(devseq_re_vs_atge_re, experiment == "DevSeq")
atge_log2_re <- subset(devseq_log2_re_vs_atge_log2_re, experiment_log == "ATGE")
devseq_log2_re <- subset(devseq_log2_re_vs_atge_log2_re, experiment_log == "DevSeq")


# Set ATGE and DevSeq sample names
atge_names <- c('experiment','gene_id','symbol','gene_id_exp','Spearman','Pearson','root_7d_ATGE',
    'root_17d_ATGE','root_21d_ATGE','hypocotyl_7d_ATGE','second_internode_24d_ATGE',
    'cotyledons_7d_ATGE','leaf_1+2_7d_ATGE','leaf_petiole_17d_ATGE','leaf_5+6_17d_ATGE',
    'leaves_senescing_35d_ATGE','apex_vegetative_7d_ATGE','apex_vegetative_14d_ATGE',
    'apex_inflorescence_21d_ATGE','flower_stg9_21d+_ATGE','flower_stg10/11_21d+_ATGE',
    'flower_stg12_21d+_ATGE','flower_stg15_21d+_ATGE','flower_stg12_sepals_ATGE',
    'flower_stg15_sepals_ATGE','flower_stg12_petals_ATGE','flower_stg15_petals_ATGE',
    'flower_stg12_stamens_ATGE','flower_stg15_stamens_ATGE','flower_stg12_carpels_ATGE',
    'flower_stg15_carpels_ATGE','flower_stg16_siliques_ATGE')

devseq_names <- c('experiment','gene_id','symbol','gene_id_exp','Spearman','Pearson','root_7d_DevSeq',
      'root_14d_DevSeq','root_21d_DevSeq','hypocotyl_10d_DevSeq','first_internode_28d_DevSeq',
      'cotyledons_7d_DevSeq','leaf_1+2_7d_DevSeq','leaf_petiole_10d_DevSeq','leaf_5+6_17d_DevSeq',
      'leaves_senescing_35d_DevSeq','apex_vegetative_7d_DevSeq','apex_vegetative_14d_DevSeq',
      'apex_inflorescence_21d_DevSeq','flower_stg9_21d+_DevSeq','flower_stg10/11_21d+_DevSeq',
      'flower_stg12_21d+_DevSeq','flower_stg15_21d+_DevSeq','flower_stg12_sepals_DevSeq',
      'flower_stg15_sepals_DevSeq','flower_stg12_petals_DevSeq','flower_stg15_petals_DevSeq',
      'flower_stg12_stamens_DevSeq','flower_stg15_stamens_DevSeq','flower_stg12_carpels_DevSeq',
      'flower_stg15_carpels_DevSeq','flower_stg16_siliques_DevSeq')

names(atge_re) <- atge_names
names(devseq_re) <- devseq_names
names(atge_log2_re) <- atge_names
names(devseq_log2_re) <- devseq_names


# Get pairwise gene correlation values
spearman_RE <- atge_re[, 5]
pearson_RE <- atge_re[, 6]
pearson_log2_RE <- atge_log2_re[, 6]


# Merge ATGE and DevSeq data frames and compute ATGE-DevSeq sample correlations
atge_devseq_re <- merge(atge_re[, c(2,7:ncol(atge_re))], devseq_re[, c(2,7:ncol(devseq_re))], by="gene_id")
atge_devseq_re_log <- merge(
    atge_log2_re[, c(2,7:ncol(atge_log2_re))], devseq_log2_re[, c(2,7:ncol(devseq_log2_re))], by="gene_id")


# Compute correlations between ATGE and DevSeq samples
atge_devseq_pearson <- diag(cor(atge_re[-1:-6], devseq_re[-1:-6], method="pearson"))
atge_devseq_log_pearson <- diag(cor(atge_log2_re[-1:-6], devseq_log2_re[-1:-6], method="pearson"))
atge_devseq_spearman <- diag(cor(atge_re[-1:-6], devseq_re[-1:-6], method="spearman"))


# Store results in /out_dir/output/plots
if (!dir.exists(file.path(out_dir, "output", "plots"))) 
  dir.create(file.path(out_dir, "output", "plots"), recursive = TRUE)
message("Storing plots in: ", file.path("output", "plots"))




#-------------------------------------- Generate plots -----------------------------------------


# Boxplot of pairwise ATGE-DevSeq gene correlations
plot_Gene_Corr <- function(data1, data2, data3) {

  png(file = file.path(out_dir, "output", "plots", "atge_devseq_gene_corr.png"), 
    width = 4250, height = 3620, res = 825)
  par(mar = c(4.5, 5, 3, 1.5))
  boxplot(data1, data2, data3, 
    ylim = c(-1, 1),
    names = c("Spearman", "Pearson", "Pearson_log"), 
    yaxt='n', 
    cex.lab = 1.2, 
    las = 1,
    cex.axis = 1.2, #adapt size of axis labels
    ylab = "Correlation coefficient", 
    col = c("#d8a900", "#00bc1f", "#00c094"), 
    boxwex = 0.71, 
    lwd = 1.4, 
    whisklty = 1, 
    at = c(1,2,3), 
    pars = list(outcol = "gray50"),
    notch = FALSE
    )
    title("Pairwise gene correlation", adj = 0.5, line = 1.4, font.main = 1, cex.main = 1.3)
    box(lwd = 1.4)
    axis(side=2, lwd = 1.4, las = 2)
    par(xpd=TRUE)
  dev.off()
}



# Boxplot of ATGE-DevSeq sample correlations
plot_Sample_Corr <- function(data1, data2, data3) {

  png(file = file.path(out_dir, "output", "plots", "atge_devseq_sample_corr.png"), 
    width = 4250, height = 3620, res = 825)
  par(mar = c(4.5, 5, 3, 1.5))
  boxplot(data1, data2, data3, 
    ylim = c(0, 1),
    names = c("Spearman", "Pearson", "Pearson_log"), 
    yaxt='n', 
    cex.lab = 1.2, 
    las = 1,
    cex.axis = 1.2, #adapt size of axis labels
    ylab = "Correlation coefficient", 
    col = c("#d8a900", "#00bc1f", "#00c094"), 
    boxwex = 0.71, 
    lwd = 1.4, 
    whisklty = 1, 
    at = c(1,2,3), 
    pars = list(outcol = "gray50"),
    notch = FALSE
    )
    title("Pairwise sample correlation", adj = 0.5, line = 1.4, font.main = 1, cex.main = 1.3)
    box(lwd = 1.4)
    axis(side=2, lwd = 1.4, las = 2)
    par(xpd=TRUE)
  dev.off()
}



# Generate correlation heatmap of merged ATGE_DevSeq data
makeCorrplot <- function(x, coefficient = c("pearson", "spearman")) {

    # Show error message if no scaling is chosen
    if (missing(coefficient))
   
       stop(
           "Please choose one of the following coefficients: 
           'pearson', 'spearman'",
           call. = TRUE
        )

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
    dfname <- deparse(substitute(x))
    coefficient_tag <- match.arg(coefficient)
    fname <- file.path(out_dir, "output", "plots",
              sprintf('%s.png', paste(dfname, coefficient_tag, sep="_"))
              )

    # Define column and row colors for color bars based on sample names and experiment
    # Build distance matrix & dendrogram then get dendrogram leaf colors to create color vector

    species_col <- c(se="palegreen4", fi="palegreen4", ap="orange2", hy="darkorchid4", fl="gold1", ro="orchid4", 
      le="olivedrab3", co="olivedrab3", ca="olivedrab3")
    
    exp_col <- c(GE="navajowhite2", eq="maroon4") # last two letters of experiment string

    x[is.na(x)] <- 0 # replaces NAs by 0
    df_t <- t(x[, 2:ncol(x)]) # transposes data frame

    # Compute correlation and build distance matrix
    if (is.element(coefficient, c("pearson"))) {
        x <- cor(x[, 2:ncol(x)], method = "pearson")
        df_t_dist.mat <- get_dist(df_t, stand = FALSE, method = "pearson")

    } else if (is.element(coefficient, c("spearman"))) {
        x <- cor(x[, 2:ncol(x)], method = "spearman")
        df_t_dist.mat <- get_dist(df_t, stand = FALSE, method = "spearman")
    }   

    df_clust.res <- hclust(df_t_dist.mat, method = "average") # agglomerate clustering

    # Get dendrogram w/ species colors
    col_dend <- dendrapply(as.dendrogram(df_clust.res), function(y){
        
        if (is.leaf(y)){
            dend_col <- species_col[substr(attr(y, "label"), 1, 2)]
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
    png(height = 5000, width = 5000, pointsize = 20, file = fname)
    par(lwd = 15) # dendrogram line width
    heatmap.2(x,
         density.info = "none",
         trace = "none",
         col = pal(800),
         ColSideColors = col_cols, 
         RowSideColors = row_cols,
         cexRow = 4.5,
         cexCol = 4.5,
         margins = c(60, 60),
         key.par = list(cex = 2.75),
         lwid = c(0.5,4,17.5), # column width
         lhei = c(0.5,4,17.5), # column height
         offsetRow = 1,
         offsetCol = 1,
         key.xlab = NA,
         key.title = NULL,
         hclustfun = function(x) hclust(x, method="average"))
    dev.off()
}



# Generate hclust dendrogram using relative expression data
makeDendrogram <- function(x, coefficient = c("pearson", "spearman")) {

    # Show error message if no scaling is chosen
    if (missing(coefficient))
   
       stop(
           "Please choose one of the following coefficients: 
           'pearson', 'spearman'",
           call. = TRUE
        )

    df_t <- t(x[, 7:ncol(x)]) # transposes data frame so rows become columns and vice versa
    df_t[is.na(df_t)] <- 0 # replaces NAs by 0

    # Set filename
    dfname <- deparse(substitute(x))
    coefficient_tag <- match.arg(coefficient)
    fname <- file.path(out_dir, "output", "plots",
              sprintf('%s_dend.png', paste(dfname, coefficient_tag, sep="_"))
              )

    # Define colors based on sample name
    label_col <- c(se="gray10", fi="gray10", ap="blue", hy="purple", fl="#e40000", ro="gray47", le="green3", 
      co="green3", ca="green3")

    # Build distance matrix
    if (is.element(coefficient, c("pearson"))) {
        df_t_dist.mat <- get_dist(df_t, stand = FALSE, method = "pearson")

    } else if (is.element(coefficient, c("spearman"))) {
      df_t_dist.mat <- get_dist(df_t, stand = FALSE, method = "spearman")
    } 

    df_clust.res <- hclust(df_t_dist.mat, method = "average") # agglomerate clustering using average linkage
  
    df_dend <- dendrapply(as.dendrogram(df_clust.res), function(n){
    
    if (is.leaf(n)){
      dend_col <- label_col[substr(attr(n,"label"),1,2)]
      attr(n, "nodePar") <- list(pch = NA, lab.col = dend_col) # to define label color
      attr(n, "edgePar") <- list(col = dend_col) # to color branch
    }
    return(n)
    })

    # make branch colors extend to last common node
    brc_col <- label_col[substr(colnames(x[, 7:ncol(x)]),1,2)]
    brc_col <- brc_col[order.dendrogram(df_dend)]
    brc_col <- factor(brc_col, unique(brc_col))

    png(height = 1195, width = 1200, pointsize = 10.94, file = fname)
    par(mar = c(14.5, 4, 4, 1.5), lwd = 8.5, cex = 3, cex.axis = 1)
    df_dend = color_branches(df_dend, clusters = as.numeric(brc_col), col = levels(brc_col))
    if (dfname == "atge_re") {
      df_dend <- flip_leaves(df_dend, c(10), c(5))
      df_dend <- flip_leaves(df_dend, c(23), c(17))
    }
    plot(df_dend)
    axis(side = 2, lwd = 3.5)
    dev.off()
}

