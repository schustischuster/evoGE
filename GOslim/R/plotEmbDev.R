# Visualize results of GO analysis 
# Input files: GO annotation files from TAIR, quantile expression gene lists of interest, 
# GOslim terms cntaining core orthologs that evolve at different rate than background



    # Heatmap showing expression of embryo development genes in Arabidopsis thaliana
    at_embr_dev <-atExpr[atExpr$gene_id %in% slim_ortho_ls[["embryo development"]]$V1,]


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

    at_embr_dev <- calculateAvgExpr(at_embr_dev)


    # Scale data to the unit interval
    scaleTPM <- function(x){(x-min(x))/(max(x)-min(x))}
    at_embr_dev_scd <- as.data.frame(t(apply(at_embr_dev[,2:ncol(at_embr_dev)], 1, scaleTPM)))


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
    steps <- c("#fae85a", "#f7ea40", "#fdc91c", "#ffa700", "#fe8300", "#f85b17", 
        "#ea2828", "#ea285a")

    pal <- color.palette(steps, c(2, 10, 11, 12, 13, 14, 5), space = "rgb")


    # Create heatmap with reversed RowSideColors
    png(height = 880, width = 1600, pointsize = 10, file = file.path(out_dir, "output", "plots", "at_embr_dev_scaled.png"))
    cc <- c(rep("#6a54a9",6), rep("#53b0db",4), rep("#2c8654",9), rep("#96ba37",3), rep("#fad819",3), rep("#e075af",4), rep("red3",8), rep("grey70",6))

    heatmap.2(as.matrix(at_embr_dev_scd), 
        density.info = "none",
        labRow = FALSE, 
        labCol = FALSE,
        dendrogram = "none", 
        col = pal(100), 
        scale = "none",
        trace = "none",
        lmat = rbind(c(0,0,0,0,0), c(0,5,0,4,0), c(0,3,0,2,0), c(0,0,0,1,0), c(0,0,0,0,0)), 
        lhei = c(0,2.5,5,0.28,0.1),
        lwid = c(0.1,2.4,0.25,5,0.5),
        key.par = list(cex = 2.8), 
        ColSideColors = cc, 
        margins = c(2, 2),
        key = TRUE,
        key.xlab = "",
        key.title = "",
        distfun = function(x) as.dist(sqrt(1/2*(1-cor(t(x))))),
        hclustfun = function(x) hclust(x, method = "average"),
        Rowv = TRUE, 
        Colv = FALSE
        )

    dev.off()


