
# Vizualize the results of the DevSeq-ATGE comparative analysis using the relative expression
# tables generated in the previous step.


#---------------------------- Read sample tables and prepare data -----------------------------


# Set input files
devseq_vs_atge_input <- "devseq_vs_atge.csv"
devseq_log2_re_vs_atge_log2_re_input <- "DevSeq_log2_RE_vs_ATGE_log2_RE.csv"


# Read data tables
message("Reading ATGE-DevSeq cor data tables")
devseq_vs_atge <- read.table(file=file.path(in_dir, devseq_vs_atge_input), sep=";", dec=".", header=TRUE, stringsAsFactors = FALSE)
devseq_log2_re_vs_atge_log2_re <- read.table(file=file.path(in_dir, devseq_log2_re_vs_atge_log2_re_input), sep=";", dec=".", header=TRUE, stringsAsFactors = FALSE)


# Prepare data
# Add experiment information to data tables
experiment <- gsub(".*_","",devseq_vs_atge[,3])
devseq_vs_atge <- cbind(as.data.frame(experiment),devseq_vs_atge)
experiment_log <- gsub(".*_","",devseq_log2_re_vs_atge_log2_re[,3])
devseq_log2_re_vs_atge_log2_re <- cbind(as.data.frame(experiment),devseq_log2_re_vs_atge_log2_re)

# Seperate ATGE and DevSeq data
atge <- subset(devseq_vs_atge, experiment == "ATGE")
devseq <- subset(devseq_vs_atge, experiment == "DevSeq")
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

names(atge) <- atge_names
names(devseq) <- devseq_names
names(atge_log2_re) <- atge_names
names(devseq_log2_re) <- devseq_names

# log2-transform devseq TPM and ATGE gcRMA data
devseq[,7:ncol(devseq)] <- log2(devseq[,7:ncol(devseq)] + 1)
atge[,7:ncol(atge)] <- log2(atge[,7:ncol(atge)] + 1)


# Get pairwise gene correlation values
gene_cor <- rbind(data.frame(Cor_value=atge_log2_re$Pearson), data.frame(Cor_value=atge_log2_re$Spearman))
cor_g_method <- data.frame(Method = rep(c("Pearson", "Spearman"), each=nrow(atge_log2_re)))
cor_g_type <- data.frame(Cor_type = rep(c("Genes"), nrow(gene_cor)))
gene_cor <- cbind(cor_g_type, cor_g_method, gene_cor)


# Merge ATGE_log2_RE and DevSeq_log2_RE data frames for correlation heatmap
atge_devseq_re_log <- merge(
    atge_log2_re[, c(2,7:ncol(atge_log2_re))], devseq_log2_re[, c(2,7:ncol(devseq_log2_re))], by="gene_id")


# Compute correlations between ATGE and DevSeq samples
atge_devseq_log_pearson <- diag(cor(atge_log2_re[-1:-6], devseq_log2_re[-1:-6], method="pearson"))
atge_devseq_log_spearman <- diag(cor(atge_log2_re[-1:-6], devseq_log2_re[-1:-6], method="spearman"))
sample_cor <- rbind(data.frame(Cor_value=atge_devseq_log_pearson), data.frame(Cor_value=atge_devseq_log_spearman))
cor_s_method <- data.frame(Method = rep(c("Pearson", "Spearman"), each=length(atge_devseq_log_pearson)))
cor_stype <- data.frame(Cor_type = rep(c("Samples"), nrow(sample_cor)))
sample_cor <- cbind(cor_stype, cor_s_method, sample_cor)
gene_sample_cor <- rbind(gene_cor, sample_cor)


# Store results in /out_dir/output/plots
if (!dir.exists(file.path(out_dir, "output", "plots"))) 
  dir.create(file.path(out_dir, "output", "plots"), recursive = TRUE)
message("Storing plots in: ", file.path("output", "plots"))




#-------------------------------------- Generate plots -----------------------------------------


# Reorder factors for correct facet order
gene_sample_cor$Cor_type <- factor(gene_sample_cor$Cor_type, levels = unique(
    gene_sample_cor$Cor_type))
gene_sample_cor$Method <- factor(gene_sample_cor$Method, levels = unique(
    gene_sample_cor$Method))


# Boxplot of pairwise ATGE-DevSeq gene and sample correlations
plotCor <- function(data) {

    x_labels <- c("Pearson" = "Pearson", "Spearman" = "Spearman")

    fname <- sprintf('%s.jpg', paste("Parwise correlations"))

    p <- ggplot(data=data, aes(x=Method, y=Cor_value)) + 
    stat_boxplot(geom ='errorbar', width = 0, size=1.125, color="black") + 
    geom_boxplot(width = 0.75, size=1.125, fatten = 2.8, color="black", outlier.shape = 1, 
        outlier.size = 2, outlier.stroke = 1.5, outlier.color = "gray55", alpha = 1.0,
        fill=c("goldenrod2", "#00bc1f", "goldenrod2", "#00bc1f")) + 
    scale_y_continuous(expand = c(0.059, 0)) + 
    scale_x_discrete(labels = x_labels)

    q <- p + theme_classic() + xlab("") + ylab("Pairwise correlation") + 
    theme(text=element_text(size = 16), 
            strip.text = element_text(size = 18.45), 
            strip.text.x = element_text(margin = margin(0.359, 0, 0.359, 0, "cm")), 
            strip.background = element_rect(colour = 'black', fill = NA, size = 2.0), 
            axis.ticks.length = unit(0.4, "cm"), 
            axis.ticks = element_line(colour = "black", size = 0.9), 
            axis.line = element_line(colour = 'black', size = 0.9), 
            plot.margin = unit(c(0.55, 0.75, 0.925, 2.425),"cm"), 
            axis.title.y = element_text(size=19.5, margin = margin(t = 0, r = 5.5, b = 0, l = 12.5), colour="black", 
                face = "plain"), 
            axis.title.x = element_text(size=19.5, margin = margin(t = 7, r = 0, b = 0, l = 0), colour="black", 
                face = "plain"), 
            axis.text.x = element_text(size=19.5, angle=0, margin = margin(t = 3.5, b = 0), colour="grey5"), 
            axis.text.y = element_text(size=16.75, angle=0, margin = margin(l = 2.5, r = 2.5), colour="grey5"), 
            panel.spacing = unit(0.5, "cm"), 
            panel.grid.major = element_blank(),
            panel.grid.minor.x = element_blank(), 
            panel.grid.minor.y = element_blank()) 

        q <- q + facet_wrap(~ Cor_type, scales = "free_y", nrow = 1)

        ggsave(file = file.path(out_dir, "output", "plots", fname), plot = q, 
            width = 10.0, height = 7.0, dpi = 300, units = c("in"), limitsize = FALSE) 
    }

    plotCor(data = gene_sample_cor)




# Generate correlation heatmap of merged ATGE_DevSeq data
makeCorrplot <- function(exp_data, coefficient = c("pearson", "spearman"), 
    clustm = c("average", "complete")) {

    # Show error message if no coefficient is chosen
    if (missing(coefficient))
   
       stop(
           "Please choose one of the following coefficients: 
           'pearson', 'spearman'",
           call. = TRUE
        )

    # Show error message if no clustering method is chosen
    if (missing(clustm))
   
       stop(
           "Please choose one of the following hclust algorithms: 
           'average', 'complete'",
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
    steps <- c("#800e0e", "#b41414", "#e11b1b", "#ffff1f", "#fcfce4")
    pal <- color.palette(steps, c(10, 24, 39, 4), space = "rgb")

    # Set filename
    dfname <- deparse(substitute(exp_data))
    fname <- sprintf('%s.png', paste(dfname, coefficient, sep="_"))

    # Define column and row colors for color bars based on sample names and experiment
    # Build distance matrix & dendrogram then get dendrogram leaf colors to create color vector

    species_col <- c(se="#009700", fi="#009700", ap="#ff9100", hy="#1d2f55", fl="#e40000", ro="#3d62b4", 
      le="#00d200", co="#00d200", ca="#00d200")
    
    exp_col <- c(GE="#d9b800", eq="#b82e90") # last two letters of experiment string

    exp_data[is.na(exp_data)] <- 0 # replaces NAs by 0
    x_df <- exp_data[, 2:ncol(exp_data)]

    # Compute correlation and build distance matrix
    x_cor <- cor(x_df, method = coefficient)
    df_t_dist.mat <- as.dist(sqrt(1/2*(1-x_cor)))

    df_clust.res <- hclust(df_t_dist.mat, method = clustm) # agglomerate clustering

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


    # Re-order (rotate) some of the clusters to make the colorbar colors more distinguishable 
    if (is.element("pearson", coefficient) && is.element("average", clustm)) {

      dend_order = dendextend::rotate(as.dendrogram(df_clust.res),c(1:2,13:23,3:12,24:32,34,33,47,46,48:52,35:45))

    } else if (is.element("pearson", coefficient) && is.element("complete", clustm)) {

      dend_order = dendextend::rotate(as.dendrogram(df_clust.res),c(1:4,7:8,5:6,13:15,10:12,9,18:19,16:17,20:52))

    } else dend_order = as.dendrogram(df_clust.res)


    # Make corrplots
        png(height = 3500, width = 3500, pointsize = 20, file = file.path(out_dir, "output", "plots", fname))
        par(lwd = 17.5) # dendrogram line width
        getRowOrder = heatmap.2(x_cor,
            revC = F,
            ColSideColors = col_cols, 
            RowSideColors = row_cols,
            distfun = function(c) as.dist(sqrt(1/2*(1-c))), 
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
            cexRow = 1.4,
            cexCol = 1.4,
            margins = c(20,20),
            key = FALSE,
            lwid = c(0.2,2.3,28.5), # column width
            lhei = c(0.2,2.3,28.5), # column height
            offsetRow = 1,
            offsetCol = 1,
            key.xlab = NA,
            key.title = NULL,
            distfun = function(c) as.dist(sqrt(1/2*(1-c))), 
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

}



# Generate hclust dendrogram using relative expression data
makeDendrogram <- function(x, coefficient = c("pearson", "spearman"), 
    clustm = c("average", "complete")) {

    # Show error message if no coefficient is chosen
    if (missing(coefficient))
   
       stop(
           "Please choose one of the following coefficients: 
           'pearson', 'spearman'",
           call. = TRUE
        )

    # Show error message if no clustering method is chosen
    if (missing(clustm))
   
       stop(
           "Please choose one of the following hclust algorithms: 
           'average', 'complete'",
           call. = TRUE
        )

    # Set filename
    dfname <- deparse(substitute(x))
    fname <- file.path(out_dir, "output", "plots",
              sprintf('%s_dend.png', paste(dfname, coefficient, sep="_"))
              )

    # Define colors based on sample name
    label_col <- c(se="#009700", fi="#009700", ap="#ff9100", hy="#1d2f55", fl="#e40000", ro="#3d62b4", 
      le="#00d200", co="#00d200", ca="#00d200")

    x_df <- x[, 7:ncol(x)]
    x_df[is.na(x_df)] <- 0 # replaces NAs by 0

    # Compute correlation and build distance matrix
    x_cor <- cor(x_df, method = coefficient)
    df_t_dist.mat <- as.dist(sqrt(1/2*(1-x_cor)))

    df_clust.res <- hclust(df_t_dist.mat, method = clustm) # agglomerate clustering

    # Re-order (rotate) some of the clusters to make the colors more distinguishable 
    if (dfname == "atge") {

      df_dend = dendextend::rotate(as.dendrogram(df_clust.res),c(1,3,2,5,4,6:11,14:15,12:13,16:24,26,25))

    } else if (dfname == "devseq") {

      df_dend = dendextend::rotate(as.dendrogram(df_clust.res),c(1:5,7,6,8:24,26,25))

    } else df_dend = as.dendrogram(df_clust.res)
  
    df_dend <- dendrapply(df_dend, function(n){
    
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
    par(mar = c(17.31, 4, 3.5, 1.5), mgp = c(3, 0.525, 0), lwd = 8.8, cex = 3, cex.axis = 0.85)
    df_dend = color_branches(df_dend, clusters = as.numeric(brc_col), col = levels(brc_col))
    plot(df_dend)
    axis(side = 2, lwd = 3.5)
    dev.off()
}



# Make facet plots of RE values for some example genes
# List of genes for plotting
genelist <- c("WUS", "REV", "AP1", "AG", "LFY", "PLT1", "FLC", "PIN1")


plotRE <- function(exp_data, genelist) {

    data <- exp_data[exp_data$symbol %in% genelist,]


    # Define specific notation
    set_scientific <- function(l) {

        if (l < 0.01) {

            # turn in to character string in scientific notation
            l <- format(l, scientific = TRUE)
            # quote the part before the exponent to keep all the digits
            l <- gsub("^(.*)e", "'\\1'e", l)
            # turn the 'e+' into plotmath format
            l <- gsub("e", "%*%10^", l)
            # return this as an expression
            parse(text=l)

        } else return(l)
    }


    # get cor p value
    pwdata <- split(data[,7:ncol(data)], rep(1:(nrow(data)/2), each = 2))
    lst <- setNames(vector('list', length(pwdata)), 1:length(pwdata))
    for(i in 1:length(pwdata)) {
        lst[[i]] <- cor.test(as.numeric(as.character(unlist(pwdata[[i]][1,]))), 
            as.numeric(as.character(unlist(pwdata[[i]][2,]))), method = "pearson")$p.value
    }
    p_value <- c(do.call(rbind, lst))

    fname <- sprintf('%s.jpg', paste("pairwise_re_values", sep = "_"))

    samples <- data.frame(samples=rep(colnames(data[7:ncol(data)]), nrow(data)))
    experiment <- data.frame(experiment=rep(rep(data$experiment[1:2], 
        each=length(7:ncol(data))), nrow(data)/2))
    symbol <- data.frame(symbol=rep(data$symbol, each=length(7:ncol(data))))
    gene_id <- data.frame(gene_id=rep(data$gene_id, each=length(7:ncol(data))))
    Pearson <- data.frame(Pearson=rep(data$Pearson, each=length(7:ncol(data))))
    RE <- data.frame(RE=as.vector(t(data[,7:ncol(data)])))

    exp_df <- cbind(symbol, gene_id, experiment, samples, RE, Pearson)
    exp_df$samples <- factor(exp_df$samples, levels=unique(exp_df$samples))
    exp_df_ATGE <- subset(exp_df, experiment == "ATGE")
    exp_df_DevSeq <- subset(exp_df, experiment == "DevSeq")

    x_labels <- c("root_whole_root_7d" = "rt", "root_whole_root_14d.17d" = "", 
        "root_whole_root_21d" = "", "hypocotyl_10d.7d" = "", "X1st_internode_24d.28d" = "st", 
        "cotyledons_7d" = "", "leaf_1_2_7d" = "", "leaf_1_2_petiole_10d.17d" = "lf", 
        "leaf_5_6_17d" = "", "leaves_senescing_35d" = "", "apex_vegetative_7d" = "", 
        "apex_vegetative_14d" = "ap", "apex_inflorescence_21d" = "", "flower_stg9_21d." = "", 
        "flower_stg10_11_21d." = "fl", "flower_stg12_21d." = "", "flower_stg15_21d." = "", 
        "flower_stg12_sepals_21d." = "", "flower_stg15_sepals_21d." = "flo", "flower_stg12_petals_21d." = "", 
        "flower_stg15_petals_21d." = "", "flower_stg12_stamens_21d." = "", "flower_stg15_stamens_21d." = "", 
        "flower_stg12_carpels_21d." = "", "flower_stg15_carpels_21d." = "ft", "fruit_stg16_siliques_28d." = "")

    line_df <- data.frame(
        x = c(0.5, 3.5, 4.5, 5.5, 10.5, 13.5, 24.5),
        y = -0.135,
        xend = c(3.5, 4.5, 5.5, 10.5, 13.5, 24.5, 26.5),
        yend = -0.135,
        colour = c("#3d62b4", "#243968", "#009700", "#00d200", "#ff9100", "#e40000", "#b10000")
        )

    peacor = unique(data$Pearson) #cor value label
    tlabel <- paste0("r = ", round(peacor, 2))

    if (length(tlabel) == 8) {
        xcoor <- c(4.15, 4.15, 4.15, 22.75, 4.15, 22.75, 4.15, 4.15)
        ycoor <- c(0.98, 0.98, 0.98, 0.98, 0.98, 0.98, 0.98, 0.98)
    } else {
        xcoor <- c(rep(4.5, length(tlabel)))
        ycoor <- c(rep(0.975, length(tlabel)))
    }

    dat_text <- data.frame(
        label = tlabel,
        symbol   = unique(data$symbol),
        x     = xcoor,
        y     = ycoor
        )

    plabel <- unlist(lapply(p_value, function(x) {
        paste0("italic('P =')~", set_scientific(signif(x, 1)))
        } #p-value of cor test label
        ))

    if (length(plabel) == 8) {
        xpcoor <- c(5.9, 5.9, 5.9, 21.235, 5.9, 21.235, 5.9, 5.9)
        ypcoor <- c(0.86, -0.01, 0.86, 0.86, 0.86, 0.86, -0.01, 0.86)
    } else {
        xpcoor <- c(rep(5.2, length(plabel)))
        ypcoor <- c(rep(0.86, length(plabel)))
    }

    pdat_text <- data.frame(
        label = plabel,
        symbol   = unique(data$symbol),
        x     = xpcoor,
        y     = ypcoor
        )


    p <- ggplot() + 
        geom_point(data=exp_df_ATGE, aes(x=samples, y=RE), shape=17, color="#ccad00", stroke=3.5) + 
        geom_line(data=exp_df_ATGE, aes(x=samples, y=RE, group=1), color="#ccad00", linetype="solid", lwd=1.4) + 
        geom_point(data=exp_df_DevSeq, aes(x=samples, y=RE), shape=16, color="#b82e90", stroke=3.5) + 
        geom_line(data=exp_df_DevSeq, aes(x=samples, y=RE, group=1), color="#b82e90", linetype="21", lwd=1.4) + 
        scale_y_continuous(expand = c(0.075, 0), limits=c(-0.16, 1.01), breaks = seq(0, 1, len = 6)) + 
        scale_x_discrete(expand = c(0.025, 0), labels = x_labels) + 
        geom_segment(data = line_df, mapping = aes(x=x, y=y, xend=xend, yend=yend, color=colour), size=5, inherit.aes = FALSE) + 
        scale_colour_identity() + 
        geom_text(data = dat_text,
            mapping = aes(x = x, y = y, label = label), size=6.2) + 
        geom_text(data = pdat_text,
            mapping = aes(x = x, y = y, label = label), size=6.2, parse=TRUE) + 
        guides(shape = guide_legend(override.aes = list(stroke=1.5)))
        

        q <- p + theme_classic() + xlab("Organ samples") + ylab("Relative expression") + 
        theme(text=element_text(size = 16), 
            strip.text = element_text(size = 17.85), 
            strip.text.x = element_text(margin = margin(0.3485, 0, 0.3485, 0, "cm")), 
            strip.background = element_rect(colour = 'black', fill = NA, size = 1.95), 
            axis.ticks.length = unit(0.2, "cm"), 
            axis.ticks = element_line(colour = "black", size = 0.9), 
            axis.line = element_line(colour = 'black', size = 0.9), 
            plot.margin = unit(c(0.15, 0.2, 0.2, 0.4),"cm"), 
            axis.title.y = element_text(size=18.9, margin = margin(t = 0, r = 5.5, b = 0, l = 11.5), colour="black", 
                face = "plain"), 
            axis.title.x = element_text(size=18.9, margin = margin(t = 4.5, r = 0, b = 0, l = 0), colour="black", 
                face = "plain"), 
            axis.text.x = element_text(size=17.0, angle=0, margin = margin(t = 2.5, b = 4), colour="grey5"), 
            axis.text.y = element_text(size=16.2, angle=0, margin = margin(l = 2.5, r = 2.5), colour="grey5"), 
            panel.spacing = unit(0.45, "cm"), 
            panel.grid.major = element_blank(),
            panel.grid.minor.x = element_blank(), 
            panel.grid.minor.y = element_blank(), 
            legend.position = "bottom", 
            legend.title = element_blank(), 
            legend.text = element_text(size = 22), 
            legend.spacing.x = unit(0.5, 'cm'), 
            legend.key.size = unit(0.95, "cm"), 
            legend.background=element_blank()) 

        q <- q + facet_wrap(~ symbol, scales = "free_x", nrow = 2)

        ggsave(file = file.path(out_dir, "output", "plots", fname), plot = q, 
            width = 17.0, height = 8.5, dpi = 300, units = c("in"), limitsize = FALSE) 
    }



# Pairwise correlation plots (facets)
plotRE(exp_data = devseq_log2_re_vs_atge_log2_re, genelist = genelist)

# Correlation heatmap of merged ATGE-DevSeq data
makeCorrplot(exp_data=atge_devseq_re_log, coefficient="pearson", clustm="complete")

# hclust dendrogram of ATGE and DevSeq data
makeDendrogram(atge, coefficient = "pearson", clustm="complete")
makeDendrogram(devseq, coefficient = "pearson", clustm="complete")

