# Read regression slope data from log and LOESS model for both DevSeq and Brawand data
# Data input: Slope values are from individual organ regressions
# DevSeq = angiosperm data, Brawand11 = mammalian data (from 2011 paper), Brawand11 = re-
# analyzed mammalian data



#-------------------------------------- Read data tables ---------------------------------------


plotSlopes <- function() {

    # Show startup message
    message("Reading data...")


    DS_AT_Br_loess_slopes = file.path(out_dir, "output", "data", "DS_AT_Br_loess_slopes.txt")
    DS_AT_Br_log_slopes = file.path(out_dir, "output", "data", "DS_AT_Br_log_slopes.txt")
    DS_AL_loess_slopes = file.path(out_dir, "output", "data", "DevSeq_AL_loess_slopes.txt")
    DS_AL_log_slopes = file.path(out_dir, "output", "data", "DevSeq_AL_log_slopes.txt")


    # Read expression data
    DS_AT_Br_loess <- read.table(DS_AT_Br_loess_slopes, header=TRUE, sep="\t", dec=".", stringsAsFactors=FALSE)
    DS_AT_Br_log <- read.table(DS_AT_Br_log_slopes, header=TRUE, sep="\t", dec=".", stringsAsFactors=FALSE)
    DS_AL_loess <- read.table(DS_AL_loess_slopes, header=TRUE, sep="\t", dec=".", stringsAsFactors=FALSE)
    DS_AL_log <- read.table(DS_AL_log_slopes, header=TRUE, sep="\t", dec=".", stringsAsFactors=FALSE)

    DS_AT_Br_loess <- DS_AT_Br_loess[-c(15:16),]
    DS_AT_Br_log <- DS_AT_Br_log[-c(15:16),]




#--------------------- Prepare data and define color palette for corrplot ---------------------


    # Create "plots" folder in /out_dir/output/plots
    if (!dir.exists(file.path(out_dir, "output", "plots"))) 
        dir.create(file.path(out_dir, "output", "plots"), recursive = TRUE)

    # Show message
    message("Starting analysis and generate plots...")



    # Prepare data for ggplot2
    DS_AL_loess_dupl <- data.frame(sample = DS_AL_loess$sample, 
        DS_AL_Br11_sOU_loess = DS_AL_loess$DS_AL_sOU_loess, DS_AL_Br_sOU_loess = DS_AL_loess$DS_AL_sOU_loess, 
        DS_AL_Br11_pea_loess = DS_AL_loess$DS_AL_pea_loess, DS_AL_Br_pea_loess = DS_AL_loess$DS_AL_pea_loess)
    Br_loess <- DS_AT_Br_loess[9:nrow(DS_AT_Br_loess),]
    colnames(Br_loess) <- colnames(DS_AL_loess_dupl)
    DS_AL_Br_loess <- rbind(DS_AL_loess_dupl, Br_loess)

    DS_AL_log_dupl <- data.frame(sample = DS_AL_log$sample, 
        DS_AL_Br11_sOU_log = DS_AL_log$DS_AL_sOU_log, DS_AL_Br_sOU_log = DS_AL_log$DS_AL_sOU_log, 
        DS_AL_Br11_pea_log = DS_AL_log$DS_AL_pea_log, DS_AL_Br_pea_log = DS_AL_log$DS_AL_pea_log)
    Br_log <- DS_AT_Br_log[9:nrow(DS_AT_Br_log),]
    colnames(Br_log) <- colnames(DS_AL_log_dupl)
    DS_AL_Br_log <- rbind(DS_AL_log_dupl, Br_log)


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
    DS_AT_Br_log <- as.data.frame(rbind(DS_AT_Br_log, getStats(DS_AT_Br_log)))
    rownames(DS_AT_Br_log) <- NULL
    DS_AL_Br_loess <- as.data.frame(rbind(DS_AL_Br_loess, getStats(DS_AL_Br_loess)))
    rownames(DS_AL_Br_loess) <- NULL
    DS_AL_Br_log <- as.data.frame(rbind(DS_AL_Br_log, getStats(DS_AL_Br_log)))
    rownames(DS_AL_Br_log) <- NULL


    Angiosperms.AT <- as.data.frame(rep("Angiosperms.AT", 16))
    colnames(Angiosperms.AT) <- "Species"
    Angiosperms.AL <- as.data.frame(rep("Angiosperms.AL", 16))
    colnames(Angiosperms.AL) <- "Species"
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
    DS_AT_sOU_log <- formatTableDS(DS_AT_Br_log, dist = "sOU", spec = "AT")
    DS_AT_Pea_log <- formatTableDS(DS_AT_Br_log, dist = "Pearson", spec = "AT")
    DS_AT <- rbind(DS_AT_sOU_loess, DS_AT_Pea_loess, DS_AT_sOU_log, DS_AT_Pea_log)
    DS_AT <- cbind(DS_AT, Angiosperms.AT)

    DS_AL_sOU_loess <- formatTableDS(DS_AL_Br_loess, dist = "sOU", spec = "AL")
    DS_AL_Pea_loess <- formatTableDS(DS_AL_Br_loess, dist = "Pearson", spec = "AL")
    DS_AL_sOU_log <- formatTableDS(DS_AL_Br_log, dist = "sOU", spec = "AL")
    DS_AL_Pea_log <- formatTableDS(DS_AL_Br_log, dist = "Pearson", spec = "AL")
    DS_AL <- rbind(DS_AL_sOU_loess, DS_AL_Pea_loess, DS_AL_sOU_log, DS_AL_Pea_log)
    DS_AL <- cbind(DS_AL, Angiosperms.AL)



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
    Br11_sOU_log <- formatTableBr(DS_AT_Br_log, dist = "sOU", study = "Br11")
    Br_sOU_log <- formatTableBr(DS_AT_Br_log, dist = "sOU", study = "Br")

    Br11_Pea_loess <- formatTableBr(DS_AT_Br_loess, dist = "Pearson", study = "Br11")
    Br_Pea_loess <- formatTableBr(DS_AT_Br_loess, dist = "Pearson", study = "Br")
    Br11_Pea_log <- formatTableBr(DS_AT_Br_log, dist = "Pearson", study = "Br11")
    Br_Pea_log <- formatTableBr(DS_AT_Br_log, dist = "Pearson", study = "Br")

    Br_2011 <- rbind(Br11_sOU_loess, Br11_sOU_log, Br11_Pea_loess, Br11_Pea_log)
    Br_2011 <- cbind(Br_2011, Mammals)

    Br_RA <- rbind(Br_sOU_loess, Br_sOU_log, Br_Pea_loess, Br_Pea_log)
    Br_RA <- cbind(Br_RA, Mammals.RA)



    # Generate final table containing all organ slopes and data sets
    DSBr_reg_slopes <- rbind(DS_AT, DS_AL, Br_2011, Br_RA)

    DSBr_reg_slopes$sample <- factor(DSBr_reg_slopes$sample, 
        levels = unique(DSBr_reg_slopes$sample))
    DSBr_reg_slopes$Slope <- as.numeric(DSBr_reg_slopes$Slope)
    DSBr_reg_slopes$Regression <- factor(DSBr_reg_slopes$Regression, 
        levels = unique(DSBr_reg_slopes$Regression))
    DSBr_reg_slopes$ID <- factor(DSBr_reg_slopes$ID, 
        levels = unique(DSBr_reg_slopes$ID))
    DSBr_reg_slopes$Species <- factor(DSBr_reg_slopes$Species, 
        levels = unique(DSBr_reg_slopes$Species))




    makeSlopePlot <- function(data) {

        fname <- sprintf('%s.jpg', paste("Organ_regression_slopes"))
            

        p <- ggplot(data=data, aes(x=ID, y=Slope)) + 
        stat_boxplot(geom ='errorbar', width = 0.45, size=1.0, color="gray15") + 
        geom_boxplot(width = 0.75, size=1.0, color="gray15", outlier.shape = 21, 
                outlier.size = 2.5, outlier.stroke = 1.5, outlier.fill = NA, outlier.color="gray35") +  
        scale_y_continuous(expand = c(0.05, 0)) + 
        guides(color = guide_legend(ncol = 3))

        q <- p + theme_bw() + xlab("Divergence time from A.lyrata (Myr)") + ylab("Pearson's r w/ A.lyrata") + 
        theme(text=element_text(size=16), 
            axis.ticks.length=unit(0.35, "cm"), 
            axis.ticks = element_line(colour = "black", size = 0.7),  
            plot.margin = unit(c(0.55, 1.175, 0.5, 0.4),"cm"), 
            axis.title.y = element_text(size=25, margin = margin(t = 0, r = 15, b = 0, l = 11), colour="black"), 
            axis.title.x = element_text(size=25, margin = margin(t = 14.75, r = 0, b = 2, l = 0), colour="black"), 
            axis.text.x = element_text(size=21.25, angle=90, margin = margin(t = 5.5), colour="black"), 
            axis.text.y = element_text(size=21.25, angle=0, margin = margin(r = 5.5), colour="black"), 
            legend.box.background = element_rect(colour = "#d5d5d5", fill=NA, size=0.3), 
            panel.border = element_rect(colour = "black", fill=NA, size=1.75), 
            panel.grid.major = element_line(color="#d5d5d5"),
            panel.grid.minor.x = element_blank(), 
            panel.grid.minor.y = element_blank(), 
            legend.position = c(0.6425, 0.8615), 
            legend.title = element_blank(), 
            legend.text = element_text(size=21.5), 
            legend.spacing.x = unit(0.5, 'cm'), 
            legend.key.size = unit(0.9, "cm"), 
            legend.background=element_blank()) 

        q <- q + facet_wrap(~ Regression, scales = "free", nrow = 1)

        ggsave(file = file.path(out_dir, "output", "plots", fname), plot = q, 
            width = 30, height = 17, dpi = 300, units = c("in"), limitsize = FALSE) 
    }

    makeSlopePlot(data = DSBr_reg_slopes)
   
}


plotSlopes()





