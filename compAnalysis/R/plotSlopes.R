# Read regression slope data from log and LOESS model for both DevSeq and Brawand data
# Data input: Slope values are from individual organ regressions
# DevSeq = angiosperm data, Brawand11 = mammalian data (from 2011 paper), Brawand11 = re-
# analyzed mammalian data



#-------------------------------------- Read data tables ---------------------------------------


plotSlopes <- function() {

    # Show startup message
    message("Reading data...")


    DS_AT_Br_loess_slopes = file.path(out_dir, "output", "data", "DS_AT_Br_loess_slopes.txt")
    DS_AT_Br_nlm_slopes = file.path(out_dir, "output", "data", "DS_AT_Br_nlm_slopes.txt")


    # Read expression data
    DS_AT_Br_loess <- read.table(DS_AT_Br_loess_slopes, header=TRUE, sep="\t", dec=".", stringsAsFactors=FALSE)
    DS_AT_Br_nlm <- read.table(DS_AT_Br_nlm_slopes, header=TRUE, sep="\t", dec=".", stringsAsFactors=FALSE)

    DS_AT_Br_loess <- DS_AT_Br_loess[-c(15:16),]




#--------------------- Prepare data and define color palette for corrplot ---------------------


    # Create "plots" folder in /out_dir/output/plots
    if (!dir.exists(file.path(out_dir, "output", "plots"))) 
        dir.create(file.path(out_dir, "output", "plots"), recursive = TRUE)

    # Show message
    message("Starting analysis and generate plots...")


    getStats <- function(df) {

        Welch_test <- sapply(df[-1], function(x) {

            t_test <- t.test(x[1:8], x[9:14])$p.value
            return(t_test)
        })

        Wilcox_test <- sapply(df[-1], function(x) {

            rs_test <- wilcox.test(x[1:8], x[9:14])$p.value
            return(rs_test)
        })

        test_stats <- as.data.frame(rbind(Welch_test, Wilcox_test))
        sample <- data.frame(sample = rownames(test_stats))
        stats_out <- cbind(sample, test_stats)
        return(stats_out)
    }


    DS_AT_Br_loess <- as.data.frame(rbind(DS_AT_Br_loess, getStats(DS_AT_Br_loess)))
    rownames(DS_AT_Br_loess) <- NULL
    DS_AT_Br_nlm <- as.data.frame(rbind(DS_AT_Br_nlm, getStats(DS_AT_Br_nlm)))
    rownames(DS_AT_Br_nlm) <- NULL


    Angiosperms.AT <- as.data.frame(rep("Angiosperms.AT", 16))
    colnames(Angiosperms.AT) <- "Species"
    Mammals <- as.data.frame(rep("Mammals", 24))
    colnames(Mammals) <- "Species"
    Mammals.RA <- as.data.frame(rep("Mammals.RA", 24))
    colnames(Mammals.RA) <- "Species"



    formatTableDS <- function(df, dist = c("sOU", "Pearson"), spec = c("AT", "AL")) {

        df_name <- paste(deparse(substitute(df)))

        regr_tag <- sub('.*\\_', '', df_name)

        if (dist == "sOU") {

            sel <- df[1:8, 1:2]

        } else if (dist == "Pearson") {

            sel <- df[1:8, c(1,4)]
        }

        colnames(sel)[2] <- "Slope"
        reg_name <- paste(dist, regr_tag, sep="_")
        Regression <- as.data.frame(rep(reg_name, nrow(sel)))
        colnames(Regression) <- "Regression"
        ID <- paste(spec, reg_name, sep="_")
        ID <- as.data.frame(rep(ID, nrow(sel)))
        colnames(ID) <- "ID"

        df_out <- cbind(sel, Regression, ID)
        return(df_out)

    }


    DS_AT_sOU_loess <- formatTableDS(DS_AT_Br_loess, dist = "sOU", spec = "AT")
    DS_AT_Pea_loess <- formatTableDS(DS_AT_Br_loess, dist = "Pearson", spec = "AT")
    DS_AT_sOU_nlm <- formatTableDS(DS_AT_Br_nlm, dist = "sOU", spec = "AT")
    DS_AT_Pea_nlm <- formatTableDS(DS_AT_Br_nlm, dist = "Pearson", spec = "AT")
    DS_AT <- rbind(DS_AT_sOU_loess, DS_AT_Pea_loess, DS_AT_sOU_nlm, DS_AT_Pea_nlm)
    DS_AT <- cbind(DS_AT, Angiosperms.AT)



    formatTableBr <- function(df, dist = c("sOU", "Pearson"), study = c("Br", "Br11")) {

        df_name <- paste(deparse(substitute(df)))

        regr_tag <- sub('.*\\_', '', df_name)

        if (dist == "sOU" && study == "Br11") {

            sel <- df[9:14, 1:2]

        } else if (dist == "sOU" && study == "Br")  {

            sel <- df[9:14, c(1,3)]

        } else if (dist == "Pearson" && study == "Br11")  {

            sel <- df[9:14, c(1,4)]

        } else if (dist == "Pearson" && study == "Br")  {

            sel <- df[9:14, c(1,5)]
        }

        colnames(sel)[2] <- "Slope"
        reg_name <- paste(dist, regr_tag, sep="_")
        Regression <- as.data.frame(rep(reg_name, nrow(sel)))
        colnames(Regression) <- "Regression"
        ID <- paste(study, reg_name, sep="_")
        ID <- as.data.frame(rep(ID, nrow(sel)))
        colnames(ID) <- "ID"

        df_out <- cbind(sel, Regression, ID)
        return(df_out)

    }


    Br11_sOU_loess <- formatTableBr(DS_AT_Br_loess, dist = "sOU", study = "Br11")
    Br_sOU_loess <- formatTableBr(DS_AT_Br_loess, dist = "sOU", study = "Br")
    Br11_sOU_nlm <- formatTableBr(DS_AT_Br_nlm, dist = "sOU", study = "Br11")
    Br_sOU_nlm <- formatTableBr(DS_AT_Br_nlm, dist = "sOU", study = "Br")

    Br11_Pea_loess <- formatTableBr(DS_AT_Br_loess, dist = "Pearson", study = "Br11")
    Br_Pea_loess <- formatTableBr(DS_AT_Br_loess, dist = "Pearson", study = "Br")
    Br11_Pea_nlm <- formatTableBr(DS_AT_Br_nlm, dist = "Pearson", study = "Br11")
    Br_Pea_nlm <- formatTableBr(DS_AT_Br_nlm, dist = "Pearson", study = "Br")

    Br_2011 <- rbind(Br11_sOU_loess, Br11_sOU_nlm, Br11_Pea_loess, Br11_Pea_nlm)
    Br_2011 <- cbind(Br_2011, Mammals)

    Br_RA <- rbind(Br_sOU_loess, Br_sOU_nlm, Br_Pea_loess, Br_Pea_nlm)
    Br_RA <- cbind(Br_RA, Mammals.RA)



    # Generate final table containing all organ slopes and data sets
    DSBr_reg_slopes <- rbind(DS_AT, Br_2011, Br_RA)

    DSBr_reg_slopes$sample <- factor(DSBr_reg_slopes$sample, 
        levels = unique(DSBr_reg_slopes$sample))
    DSBr_reg_slopes$Slope <- as.numeric(DSBr_reg_slopes$Slope)
    DSBr_reg_slopes$Regression <- factor(DSBr_reg_slopes$Regression, 
        levels = unique(DSBr_reg_slopes$Regression))
    DSBr_reg_slopes$ID <- factor(DSBr_reg_slopes$ID, 
        levels = unique(DSBr_reg_slopes$ID))
    DSBr_reg_slopes$Species <- factor(DSBr_reg_slopes$Species, 
        levels = unique(DSBr_reg_slopes$Species))


    # Change names for facet strip
    DSBr_reg_slopes[,3] <- DSBr_reg_slopes[,3] %<>% 
        gsub("sOU_loess", "sOU-v LOESS", .) %>% 
        gsub("Pearson_loess", "Pearson LOESS", .) %>% 
        gsub("sOU_nlm", "sOU-v NE", .) %>% 
        gsub("Pearson_nlm", "Pearson NE", .)


    # Reorder factors for correct facet order
    DSBr_reg_slopes$Regression <- factor(DSBr_reg_slopes$Regression, levels=c(
        "sOU-v LOESS", "sOU-v NE", "Pearson LOESS", "Pearson NE"))


    # Generate dummy data to scale y_max in individual facets
    scale_layer <- rbind(DSBr_reg_slopes[1,], DSBr_reg_slopes[9,], DSBr_reg_slopes[17,], 
        DSBr_reg_slopes[25,])
    max_values <- DSBr_reg_slopes %>% group_by(Regression) %>% summarize(max.Slope = max(Slope))
    max_values$max.Slope <- (max_values$max.Slope) * 1.148
    scale_layer$Slope <- c(max_values$max.Slope[1], max_values$max.Slope[3], max_values$max.Slope[2], 
        max_values$max.Slope[4])
    scale_layer$Slope <- as.numeric(scale_layer$Slope)




    makeSlopePlot <- function(data, data2, dist_meth) {

        if (dist_meth == "pea") {
            x_labels <- c("AT_Pearson_loess" = "Angiosperms.AT", "AT_Pearson_nlm" = "Angiosperms.AT",  
            "Br11_Pearson_loess" = "Mammals.11", "Br11_Pearson_nlm" = "Mammals.11", 
            "Br_Pearson_loess" = "Mammals.re-an.", "Br_Pearson_nlm" = "Mammals.re-an.")

        } else if (dist_meth == "sOU_v") {
            x_labels <- c("AT_sOU_loess" = "Angiosperms.AT", "AT_sOU_nlm" = "Angiosperms.AT",   
            "Br11_sOU_loess" = "Mammals.11", "Br11_sOU_nlm" = "Mammals.11", 
            "Br_sOU_loess" = "Mammals.re-an.", "Br_sOU_nlm" = "Mammals.re-an.")
        }

        fname <- sprintf('%s.jpg', paste("Organ_regression_slopes", dist_meth, sep = "_"))
            

        p <- ggplot(data=data, aes(x=ID, y=Slope)) + 
        stat_boxplot(geom ='errorbar', width = 0, size=1.0, color="black") + 
        geom_boxplot(width = 0.75, size=1.05, fatten = 3.5, color="black", outlier.shape = NA, alpha = 0) + 
        geom_point(aes(shape = sample, color = Species, size = sample, stroke = 3.0)) + 
        scale_shape_manual(values = c(0, 8, 2, 5, 3, 4, 6, 1, 15, 16, 17, 18, 19, 10)) + 
        scale_size_manual(values = c(4.75, 4.75, 4.75, 4.75, 5.25, 5.25, 4.75, 5.25, 5.5, 5.5, 5.5, 8.8, 5.5, 5.5)) + 
        scale_color_manual(values=c('#728acb', 'red', 'red3'), 
            guide = "none") + 
        scale_fill_manual(values=c('#728acb', 'red', 'red3'), 
            guide = "legend") + 
        scale_y_continuous(expand = c(0.075, 0), labels = comma) + 
        scale_x_discrete(labels = x_labels) + 
        guides(shape = guide_legend(override.aes = list(stroke=1.5)))

        q <- p + theme_classic() + xlab("Data set") + ylab("Slope value") + 
        theme(text=element_text(size = 16), 
            strip.text = element_text(size = 23.75), 
            strip.text.x = element_text(margin = margin(0.44, 0, 0.44, 0, "cm")), 
            strip.background = element_rect(colour = 'black', fill = NA, size = 2.4), 
            axis.ticks.length = unit(0.325, "cm"), 
            axis.ticks = element_line(colour = "black", size = 1.15), 
            axis.line = element_line(colour = 'black', size = 1.15), 
            plot.margin = unit(c(0.55, 1.175, 1, 0.4),"cm"), 
            axis.title.y = element_text(size=24.6, margin = margin(t = 0, r = 13.0, b = 0, l = 10.8), colour="black", 
                face = "bold"), 
            axis.title.x = element_text(size=24.6, margin = margin(t = 0, r = 0, b = 0, l = 0), colour="black", 
                face = "bold"), 
            axis.text.x = element_text(size=22.5, angle=45, margin = margin(t = -61, b = 80), colour="grey5", 
                hjust = 0.99, vjust = 0.45), 
            axis.text.y = element_text(size=21.75, angle=0, margin = margin(l = 2.5, r = 2.5), colour="grey5"), 
            panel.spacing = unit(0.7, "cm"), 
            panel.grid.major = element_blank(),
            panel.grid.minor.x = element_blank(), 
            panel.grid.minor.y = element_blank(), 
            legend.position = "right", 
            legend.title = element_blank(), 
            legend.text = element_text(size = 22.5), 
            legend.spacing.x = unit(0.5, 'cm'), 
            legend.key.size = unit(1, "cm"), 
            legend.background=element_blank()) 

        q <- q + facet_wrap(~ Regression, scales = "free_x", nrow = 1) + 
          geom_blank(data=data2, aes(x=ID, y=Slope))

        ggsave(file = file.path(out_dir, "output", "plots", fname), plot = q, 
            width = 11, height = 8.5, dpi = 300, units = c("in"), limitsize = FALSE) 
    }

    makeSlopePlot(data = DSBr_reg_slopes[c(1:8,17:24,33:44,57:68),], data2 = scale_layer[c(1,3),], dist_meth = "sOU_v")
    makeSlopePlot(data = DSBr_reg_slopes[c(9:16,25:32,45:56,69:80),], data2 = scale_layer[c(2,4),], dist_meth = "pea")




    # Generate blank plot and add annotations (workaround to add Wilcox.test p-values)
    # This plot has same dimensions as "makeSlopePlot", but features complete transparency
    # with added text annotations

    df_blank <- data.frame()

    AT_Br11_sOU_loess <- c(paste("italic('P =')~", formatC(DS_AT_Br_loess[16,2], format = "e", digits = 0)))
    AT_Br_sOU_loess <- c(paste("italic('P =')~", formatC(DS_AT_Br_loess[16,3], format = "e", digits = 0)))
    AT_Br11_pea_loess <- c(paste("italic('P =')~", formatC(DS_AT_Br_loess[16,4], format = "e", digits = 0)))

    AT_Br11_sOU_nlm <- c(paste("italic('P =')~", formatC(DS_AT_Br_nlm[16,2], format = "e", digits = 0)))
    AT_Br_sOU_nlm <- c(paste("italic('P =')~", formatC(DS_AT_Br_nlm[16,3], format = "e", digits = 0)))
    AT_Br11_pea_nlm <- c(paste("italic('P =')~", formatC(DS_AT_Br_nlm[16,4], format = "e", digits = 0)))


    dat_text <- data.frame(
        label = c(AT_Br11_sOU_loess, AT_Br_sOU_loess, AT_Br11_sOU_nlm, AT_Br_sOU_nlm, 
            AT_Br11_pea_loess, AT_Br11_pea_nlm), # other p-values have same value
        x = c(0.481, 0.595, 0.5387, 1.802, 3.01, 4.22),
        y = c(8.825, 9.47, 4.538, 9.47, 9.225, 9.225)
    )

    options(scipen = -1) # This sets scientific notation below 1e-2

    p <- ggplot(df_blank) + geom_point() + xlim(0, 5) + ylim(0, 10) + theme_void() + 
    theme(
        panel.background = element_rect(fill = "transparent"), # bg of the panel
        plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
        panel.grid.major = element_blank(), # get rid of major grid
        panel.grid.minor = element_blank(), # get rid of minor grid
        legend.background = element_rect(fill = "transparent"), # get rid of legend bg
        legend.box.background = element_rect(fill = "transparent") # get rid of legend panel bg
        )

    q <- p + geom_text(
        data = dat_text,
        mapping = aes(x = x, y = y, label = format(label, scientific=TRUE)), 
        size = 8.25, parse = TRUE
    )

    png("NUL")

    ggsave(file = file.path(out_dir, "output", "plots", "p_value_annotation_layer.png"), plot = q,
        width = 28.5, height = 12.8, units = c("in"), dpi = 300, limitsize = FALSE, bg = "transparent")

   
}


plotSlopes()





