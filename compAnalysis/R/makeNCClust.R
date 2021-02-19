# Prepare DevSeq comparative lncRNAs ortholog gene expression data
# Thresholds: DevSeq 0.05 ERCC 
# Data input: DevSeq TPM expression tables of all samples


#-------------------------------------- Read data tables ---------------------------------------


makeNCClust <- function(expr_estimation = c("TPM", "counts"), coefficient = c("pearson", "spearman"), 
    devseq_organs = c("all", "subset"), transcripttype = c("coding", "non-coding") ) {
	

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
    if ((missing(devseq_organs)) || (!is.element(devseq_organs, c("all", "subset"))))
   
       stop(
       "Please choose one of the available DevSeq organ sets: 
       'all', 'subset'",
       call. = TRUE
       )

    # Show error message if no transcript type or unknown type is chosen
    if ((missing(transcripttype)) || (!is.element(transcripttype, c("coding", "non-coding"))))
   
       stop(
       "Please choose one of the available transcript types: 
       'coding', 'non-coding'",
       call. = TRUE
       )


    # Show startup message
    message("Reading data...")


	# Set expression input file
    if ((is.element("TPM", expr_estimation)) && (is.element("non-coding", transcripttype))) {
        genesExpr = file.path(in_dir, "Expression_data", "lnc_AT_brass_inter_tpm_mat_deseq_sample_names.csv")

    } else if ((is.element("counts", expr_estimation)) && (is.element("non-coding", transcripttype))) {
        genesExpr = file.path(in_dir, "Expression_data", "lnc_AT_brass_inter_count_mat_vsd_sample_names.csv")
    
    } else if ((is.element("TPM", expr_estimation)) && (is.element("coding", transcripttype))) {
        genesExpr = file.path(in_dir, "Expression_data", "AT_brass_inter_tpm_mat_deseq_sample_names.csv")

    } else if ((is.element("counts", expr_estimation)) && (is.element("coding", transcripttype))) {
        genesExpr = file.path(in_dir, "Expression_data", "AT_brass_inter_count_mat_vsd_sample_names.csv")
    }
    


	# Define simplified Brawand and DevSeq column names
    col_names <- rep(c("Root", "Hypocotyl", "Leaf", "veg_apex", "inf_apex", 
            "Flower", "Stamen", "Carpel"), each=12)
    replicate_tag_samples <- rep(c(".1",".2",".3"), times=4)
    col_names <- paste0(col_names,replicate_tag_samples)
    spec_names <- rep(c("_AT", "_AL", "_CR", "_ES"), each=3, times=8)
    col_names <- paste0(col_names, spec_names)
    col_names <- c("gene_id", col_names) 


	
	x <- read.table(genesExpr, sep=";", dec=".", header=TRUE, stringsAsFactors=FALSE)


    # Stop function here to allow specific analysis of a single data set
    # return_list <- list("expr_estimation" = expr_estimation, "x" = x, "coefficient" = coefficient, "col_names" = col_names, "devseq_organs" = devseq_organs, "transcripttype" = transcripttype)
    # return(return_list)
    # }
    # return_objects <- makeNCClust(expr_estimation="counts", coefficient="pearson", devseq_organs="all", transcripttype="non-coding") # read in DevSeq expression data
    # list2env(return_objects, envir = .GlobalEnv)

    if (transcripttype == "coding") x <- x[,1:97]


    # Set column names
    colnames(x) <- col_names




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
    steps <- c("#e14134", "#e14134", "#fab141", "#fff200", "#fcfcf7")

    pal <- color.palette(steps, c(30, 29, 21, 15), space = "rgb")

    # Set filename
    fname <- sprintf('%s.png', paste("Brassicaceae", transcripttype, expr_estimation, coefficient, sep="_"))

    if (devseq_organs == "subset") {

        fname <- sprintf('%s.png', paste("Brassicaceae", transcripttype, expr_estimation, coefficient, devseq_organs, sep="_"))
    }


    # Define column and row colors for color bars based on sample names and experiment
    # Build distance matrix & dendrogram then get dendrogram leaf colors to create color vector

    exp_col <- c(AT="gray1", AL="cornsilk3", CR="#fff7e8", ES="#95928c") # last two letters of sample name

    species_col <- c(Roo="#6a54a9", Hyp="#53b0db", Lea="#2c8654", veg="#96ba37", 
        inf="#fad819", Flo="#e075af", Sta="#ed311c", Mat="#a63126", Car="#f2a72f")


    x[is.na(x)] <- 0 # replaces NAs by 0

    # Remove ERCC spike-ins from data
    x <- x[!grepl("ERCC", x$gene_id),]

    x_df <- x

    # Remove apex and flower samples for hclust heatmap if subset is chosen 
    if (devseq_organs == "subset") {

        x_df <- x_df %>% select (-c(
            veg_apex.1_AT, veg_apex.2_AT, veg_apex.3_AT, veg_apex.1_AL, veg_apex.2_AL, 
            veg_apex.3_AL, veg_apex.1_CR, veg_apex.2_CR, veg_apex.3_CR, veg_apex.1_ES, veg_apex.2_ES, veg_apex.3_ES, 
            inf_apex.1_AT, inf_apex.2_AT, inf_apex.3_AT, inf_apex.1_AL, inf_apex.2_AL, 
            inf_apex.3_AL, inf_apex.1_CR, inf_apex.2_CR, inf_apex.3_CR, inf_apex.1_ES, inf_apex.2_ES, inf_apex.3_ES, 
            Flower.1_AT, Flower.2_AT, Flower.3_AT, Flower.1_AL, Flower.2_AL, 
            Flower.3_AL, Flower.1_CR, Flower.2_CR, Flower.3_CR, Flower.1_ES, Flower.2_ES, Flower.3_ES))
    }




#------------------------ Prepare DevSeq and Brawand data for corrplot -------------------------


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
    if (is.element("spearman", coefficient) && is.element("counts", expr_estimation) && 
        is.element("subset", devseq_organs) && is.element("non-coding", transcripttype)) {

        dend_order=dendextend::rotate(as.dendrogram(df_clust.res),c(1:6,10:12,7:9,13:24,31:36,28:30,25:27,43:48,40:42,37:39,55:60,52:54,49:51))

    } else if (is.element("spearman", coefficient) && is.element("counts", expr_estimation) && 
        is.element("subset", devseq_organs) && is.element("coding", transcripttype)) {

        dend_order=dendextend::rotate(as.dendrogram(df_clust.res),c(1:3,10:12,7:9,4:6,25:36,13:24,37:42,46:48,43:45,49:54,58:60,55:57))

    } else dend_order = TRUE 




#---------------------------- Make corrplot for DevSeq and Brawand -----------------------------



    # Make corrplots
    png(height = 3500, width = 3500, pointsize = 20, file = file.path(out_dir, "output", "plots", fname))
    par(lwd = 20) # dendrogram line width
    getRowOrder = heatmap.2(x_cor,
        revC = F,
        ColSideColors = col_cols, 
        RowSideColors = row_cols,
        distfun = function(c) as.dist(sqrt(1/2*(1-cor(x_df[, 2:ncol(x_df)], method = coefficient)))), 
        hclustfun = function(x) hclust(x_dist, method = "average"), 
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
        cexRow = 3.25,
        cexCol = 3.25,
        margins = c(22, 22),
        key = FALSE,
        lwid = c(0.2,2.3,28.5), # column width
        lhei = c(0.2,2.3,28.5), # column height
        offsetRow = 1,
        offsetCol = 1,
        key.xlab = NA,
        key.title = NULL,
        distfun = function(c) as.dist(sqrt(1/2*(1-cor(x_df[, 2:ncol(x_df)], method = coefficient)))), 
        hclustfun = function(x) hclust(x_dist, method = "average"), 
        # Order dendrogram in a way that it starts with distant species BD, MT, TH
        Rowv = dend_order, 
        Colv = "Rowv"
        )

    dev.off()

    # Save colorbar
    png(height = 1000, width = 850, pointsize = 20, file = file.path(out_dir, "output", "plots", "DevSeq_lnc_comp_colorbar.png"))
    par(cex = 7, mar = c(1,2.75,1,1), cex.lab = 1.25)
    x_c = 1
    y_c = seq(0,1,len = 100)
    z_c = matrix(1:100, nrow = 1)
    ylabel <- seq(0, 1, by = 0.5)
    image(x_c,y_c,z_c,col = pal(800), axes = FALSE, xlab = "", ylab = "", cex=10)
    axis(2, at = ylabel, las = 1, lwd = 15)
    dev.off()
    



#------------------------------------------ Make PCAs ------------------------------------------


  if ((transcripttype == "non-coding") && (devseq_organs == "all")) { 
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
    plotPCA <- function(data, pc_var1, pc_var2, set=c("pc1_2", "pc2_3", "pc1_3")) {

        fname <- sprintf('%s.png', paste(transcripttype, deparse(substitute(data)), "pca", expr_estimation, devseq_organs, data_scale, sep="_"))
        
        if (set == "pc1_2") {
            x_coord <- "PC1 ("
            y_coord <- "PC2 ("

        } else if(set == "pc2_3") {
            x_coord <- "PC2 ("
            y_coord <- "PC3 ("

        } else if(set == "pc1_3") {
            x_coord <- "PC1 ("
            y_coord <- "PC3 ("
        }

        if (data_scale == "scaled") {
            legend_pos <- "none"
            legend_spacing <- 0
            order_guide <- 0

        } else if (data_scale == "raw") {
            legend_pos <- c(0.5, 0.23125)
            legend_spacing <- 7.7
            order_guide <- 1
        }

        if (set == "pc2_3") legend_pos <- "none"
        if (set == "pc1_3") legend_pos <- "none"


        x_lab <- paste(x_coord, pc_var1, "%)", sep="")
        y_lab <- paste(y_coord, pc_var2, "%)", sep="")

        plot <- ggplot(data, aes(x = PC1, y = PC2, colour=Organ, fill = Organ)) + 
        stat_bag(prop = 0.95, size=1.5) + 
        geom_point(aes(shape=Species, color=Organ, size=Species, stroke=2.25)) + 
        scale_shape_manual(values = c(16, 17, 18, 15), labels = c("A.lyr.","A.thal.","C.rub.","E.sals.")) + 
        scale_x_continuous(expand = c(0.059, 0), limits = NULL) + 
        scale_y_continuous(expand = c(0.07, 0)) +
        guides(shape = guide_legend(override.aes = list(size = c(5.75, 5.5, 7.5, 5.5), stroke=2.75), 
               order = order_guide, ncol = 1), size = FALSE) + 
        guides(colour = guide_legend(override.aes = list(size=5, linetype = "blank", alpha=1))) + 
        scale_size_manual(values = c(5.75, 5.5, 7.5, 5.5)) + 
        # shapes = filled round, filled rect, empty square, filled square_rot, filled square, empty rect, inverted empty rect
        # colors = dark moderate violet, soft blue, dark green, moderate green, vivid yellow, orange, vivid red, soft pink
        scale_color_manual(values=c('#6a54a9','#53b0db', '#2c8654', '#96ba37','#fad819', '#f2a72f', '#ee412e', '#e075af'), 
            guide = "legend") + 
        scale_fill_manual(values=c('#6a54a9','#53b0db', '#2c8654', '#96ba37','#fad819', '#f2a72f', '#ee412e', '#e075af'), 
            guide = "legend") + 
        labs(x = x_lab, y = y_lab) + 
        theme(panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(), 
            panel.border = element_rect(colour = "black", fill=NA, size=2.3), 
            panel.background = element_blank(), 
            axis.title.y = element_text(size=25.75, margin = margin(t = 0, r = -3.5, b = 0, l = 18.0), colour="black"), 
            axis.title.x = element_text(size=25.75, margin = margin(t = 12.75, r = 0, b = 4.5, l = 0), colour="black"), 
            axis.text.x = element_text(size=24.5, angle=0, margin = margin(t = 5), colour="black"), 
            axis.text.y = element_text(size=24.5, angle=0, margin = margin(r = 5), colour="black"), 
            axis.ticks.length=unit(0.35, "cm"), 
            axis.ticks = element_line(colour = "black", size = 1.15), 
            legend.key = element_rect(colour = "transparent", fill = "white"), 
            legend.key.size = unit(1.65,"line"), # default is 1.2
            legend.text = element_text(size=24.5), 
            legend.background = element_rect(fill = "transparent"),
            legend.title = element_blank(),
            legend.spacing.x = unit(legend_spacing,'cm'), 
            legend.position = legend_pos,
            legend.box="horizontal", 
            legend.box.just = "bottom", 
            plot.margin = unit(c(0.55, 0.55, 0.5, 0.5),"cm"))

        ggsave(file = file.path(out_dir, "output", "plots", fname), plot = plot,
            width = 8.35, height = 8, dpi = 300, units = c("in"), 
            limitsize = FALSE)
    }


    if (data_scale == "raw") {
        plotPCA(data=DevSeq_pca_1_2_w_stamen, pc_var1=DevSeq_pc1_var_w_stamen, pc_var2=DevSeq_pc2_var_w_stamen, 
            set="pc1_2") 

        plotPCA(data=DevSeq_pca_2_3_w_stamen, pc_var1=DevSeq_pc2_var_w_stamen, pc_var2=DevSeq_pc3_var_w_stamen, 
            set="pc2_3") 

        plotPCA(data=DevSeq_pca_1_3_w_stamen, pc_var1=DevSeq_pc1_var_w_stamen, pc_var2=DevSeq_pc3_var_w_stamen, 
            set="pc1_3") 

    } else if (data_scale == "scaled") {
        plotPCA(data=DevSeq_pca_1_2_w_stamen, pc_var1=DevSeq_pc1_var_w_stamen, pc_var2=DevSeq_pc2_var_w_stamen, 
            set="pc1_2") 
    }

  }

  performPCA(data = x_df, data_scale = "raw", ntop = nrow(x_df))
  performPCA(data = x_df, data_scale = "scaled", ntop = nrow(x_df))

 }

}


makeCompAnylsis(dataset="DevSeq", expr_estimation="TPM", coefficient="pearson", devseq_spec="all", data_norm="inter-organ", devseq_organs="all")
makeCompAnylsis(dataset="DevSeq", expr_estimation="counts", coefficient="pearson", devseq_spec="Brassicaceae", data_norm="inter-organ", devseq_organs="all")
makeCompAnylsis(dataset="DevSeq", expr_estimation="counts", coefficient="pearson", devseq_spec="all", data_norm="inter-organ", devseq_organs="all")
makeCompAnylsis(dataset="DevSeq", expr_estimation="counts", coefficient="pearson", devseq_spec="all", data_norm="inter-organ", devseq_organs="subset")
makeCompAnylsis(dataset="Brawand", expr_estimation="counts", coefficient="pearson", data_norm="inter-organ")





