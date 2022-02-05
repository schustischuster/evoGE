# Read regression slope data from log and LOESS model for both DevSeq and Brawand data
# Data input: Slope values are from individual organ regressions
# DevSeq = angiosperm data, Brawand11 = mammalian data (from 2011 paper), Brawand11 = re-
# analyzed mammalian data



#-------------------------------------- Read data tables ---------------------------------------


plotTreeLength <- function() {

    # Show startup message
    message("Reading data...")


    ATH_all_bs_tree_length = file.path(out_dir, "output", "data", "phyloBSTreeL_all_ATH.txt")
    ATH_Brassicaceae_bs_tree_length = file.path(out_dir, "output", "data", "phyloBSTreeL_Brassicaceae_ATH.txt")


    # Read expression data
    AT_all_tree_length <- read.table(ATH_all_bs_tree_length, header=TRUE, sep="\t", dec=".", stringsAsFactors=FALSE)
    AT_Brassicaceae_tree_length <- read.table(ATH_Brassicaceae_bs_tree_length, header=TRUE, sep="\t", dec=".", stringsAsFactors=FALSE)




#--------------------- Prepare data and define color palette for corrplot ---------------------


    # Create "plots" folder in /out_dir/output/plots
    if (!dir.exists(file.path(out_dir, "output", "plots"))) 
        dir.create(file.path(out_dir, "output", "plots"), recursive = TRUE)

    # Show message
    message("Starting analysis and generate plots...")



    # Prepare data for ggplot2
    dataset_angio <- rep(c("Angiosperms"), nrow(AT_all_tree_length))
    dataset_brass <- rep(c("Brassicaceae"), nrow(AT_Brassicaceae_tree_length))
    dataset <- data.frame(Order = c(dataset_angio, dataset_brass))

    bs_tree_length <- rbind(AT_all_tree_length, AT_Brassicaceae_tree_length)
    bs_tree_length <- cbind(bs_tree_length, dataset)
    colnames(bs_tree_length) <- c("Organ", "Total_tree_length", "Observed_tree_length", "Order")
    bs_tree_length$Organ <- factor(bs_tree_length$Organ)
    bs_tree_length$Order <- factor(bs_tree_length$Order)
    bs_tree_length$Total_tree_length <- as.numeric(bs_tree_length$Total_tree_length)
    bs_tree_length$Observed_tree_length <- as.numeric(bs_tree_length$Observed_tree_length)


    # Create data tables w/ and wo/ pollen
    bs_tree_length_wo_pollen <- bs_tree_length[bs_tree_length$Organ != "Pollen", ]
    bs_tree_length_pollen <- bs_tree_length[bs_tree_length$Organ == "Pollen", ]
    bs_tree_length_stamen <- bs_tree_length[bs_tree_length$Organ == "Stamen", ]
    bs_tree_length_stamen_pollen <- rbind(bs_tree_length_stamen, bs_tree_length_pollen)


    observed_values <- unique(bs_tree_length$Observed_tree_length)
    unique_organ <- rep(unique(bs_tree_length$Organ), 2)
    unique_order <- rep(unique(bs_tree_length$Order), each = 9)
    observed_tree_length <- data.frame(unique_organ, observed_values, unique_order)
    colnames(observed_tree_length) <- c("Organ", "Total_tree_length", "Order")


    # Create data tables w/ and wo/ pollen
    observed_tree_length_wo_pollen <- observed_tree_length[observed_tree_length$Organ != "Pollen", ]
    observed_tree_length_pollen <- observed_tree_length[observed_tree_length$Organ == "Pollen", ]
    observed_tree_length_stamen <- observed_tree_length[observed_tree_length$Organ == "Stamen", ]
    observed_tree_length_stamen_pollen <- rbind(observed_tree_length_stamen, observed_tree_length_pollen)

    bs_tree_length_wo_pollen$Organ <- factor(bs_tree_length_wo_pollen$Organ, levels = unique(bs_tree_length_wo_pollen$Organ))
    observed_tree_length_wo_pollen$Organ <- factor(observed_tree_length_wo_pollen$Organ, levels = unique(observed_tree_length_wo_pollen$Organ))
    bs_tree_length_stamen_pollen$Organ <- factor(bs_tree_length_stamen_pollen$Organ, levels = unique(bs_tree_length_stamen_pollen$Organ))
    observed_tree_length_stamen_pollen$Organ <- factor(observed_tree_length_stamen_pollen$Organ, levels = unique(observed_tree_length_stamen_pollen$Organ))


    plotTTL <- function(data, data2, spec_order) {

        if (spec_order == "AS") {

            y_breaks = c(1.3, 1.4, 1.5, 1.6, 1.7)
            y_labels = c(1.3, 1.4, 1.5, 1.6, 1.7)
            y_exp <- 0.025

        } else { 

            y_breaks = c(0.55, 0.6, 0.65, 0.7)
            y_labels = c("", 0.6, "", 0.7)
            y_exp <- 0.035
        }


        fname <- sprintf('%s.jpg', paste("Total_tree_length", spec_order, sep="_"))
            
        p <- ggplot(data=data, aes(x = Organ, y = Total_tree_length, group = Organ)) + 
        stat_boxplot(geom ='errorbar', width = 0.35, size=1.7, color="black") + 
        geom_boxplot(width = 0.85, size = 1.7, fatten = 1.35, color="black", outlier.shape = NA, 
            alpha = 1, aes(fill = Organ)) + 
        geom_point(data = data2, position = position_dodge(width=0.75), size = 5, col = "red2", 
            aes(x = Organ, y = Total_tree_length)) + 
        scale_fill_manual(values=c('#6a54a9','#53b0db', '#2c8654', '#96ba37', '#fad819','#e075af', '#f2a72f', '#ed311c')) + 
        scale_y_continuous(expand = c(y_exp, 0), breaks = y_breaks, labels = y_labels) + 
        scale_x_discrete(labels=c("Root" = "Root", "Hypocotyl" = "Hypocotyl", 
            "Leaf" = "Leaf", "veg_apex" = "Apex veg", "inf_apex" = "Apex inf", 
            "Flower" = "Flower", "Carpel" = "Carpel", "Stamen" = "Stamen"), 
            expand = c(0.035, 0))

        q <- p + theme_classic() + ylab("Total tree length") + 
        theme(text=element_text(size = 16), 
            strip.text = element_text(size = 36.25), 
            strip.text.x = element_text(margin = margin(0.6, 0, 0.6, 0, "cm")), 
            strip.background = element_rect(colour = 'black', fill = NA, size = 3.25), 
            axis.ticks.length = unit(0.45, "cm"), 
            axis.ticks = element_line(colour = "black", size = 1.7), 
            axis.line = element_line(colour = 'black', size = 1.7), 
            plot.margin = unit(c(0.5, 0.015, 0.5, 0.475),"cm"), 
            axis.title.y = element_text(size=36.25, margin = margin(t = 0, r = 18.5, b = 0, l = 2), colour="black", 
                face = "plain"), 
            axis.title.x = element_blank(), 
            axis.text.x = element_text(size=36.25, angle=45, margin = margin(t = -52, b = 42.5), colour="black", 
                hjust = 0.99, vjust = 0.5), 
            axis.text.y = element_text(size=36.25, angle=0, margin = margin(r = 5), colour="black"), 
            panel.spacing = unit(0.5, "cm"), 
            panel.grid.major = element_blank(),
            panel.grid.minor.x = element_blank(), 
            panel.grid.minor.y = element_blank(), 
            legend.position = "none", 
            legend.title = element_blank()) 

        ggsave(file = file.path(out_dir, "output", "plots", fname), plot = q, 
            width = 8.25, height = 12, dpi = 300, units = c("in"), limitsize = FALSE) 
    }

    plotTTL(data = bs_tree_length_wo_pollen[bs_tree_length_wo_pollen$Order == "Angiosperms", ], 
        data2 = observed_tree_length_wo_pollen[observed_tree_length_wo_pollen$Order == "Angiosperms", ], 
        spec_order = "AS")
    plotTTL(data = bs_tree_length_wo_pollen[bs_tree_length_wo_pollen$Order == "Brassicaceae", ], 
        data2 = observed_tree_length_wo_pollen[observed_tree_length_wo_pollen$Order == "Brassicaceae", ], 
        spec_order = "BR")



    plotTTL.SO <- function(data, data2, spec_order) {

        if (spec_order == "AS") {

            y_breaks = c(1.6, 1.8, 2.0, 2.2, 2.4)
            y_labels = c(1.6, 1.8, sprintf(2.0, fmt = '%#.1f'), 2.2, 2.4)

        } else { 

            y_breaks = c(0.75, 0.8, 0.85, 0.9, 0.95, 1.0)
            y_labels = c("", 0.8, "", 0.9, "", sprintf(1.0, fmt = '%#.1f'))
        }

        fname <- sprintf('%s.jpg', paste("Total_tree_length_pollen", spec_order, sep="_"))
            
        p <- ggplot(data=data, aes(x = Organ, y = Total_tree_length, group = Organ)) + 
        stat_boxplot(geom ='errorbar', width = 0.35, size=1.7, color="black") + 
        geom_boxplot(width = 0.9, size = 1.7, fatten = 1.35, color="black", outlier.shape = NA, 
            alpha = 1, aes(fill = Organ)) + 
        geom_point(data = data2, position = position_dodge(width=0.75), size = 5, col = "red2", 
            aes(x = Organ, y = Total_tree_length)) + 
        scale_fill_manual(values=c('#ed311c', '#a63126')) + 
        scale_y_continuous(expand = c(0.035, 0), breaks = y_breaks, labels = y_labels) + 
        scale_x_discrete(labels=c("Stamen" = "Stamen", "Pollen" = "Pollen"), expand = c(0.125, 0))

        q <- p + theme_classic() + ylab("Total tree length") + 
        theme(text=element_text(size = 16), 
            strip.text = element_text(size = 36.25), 
            strip.text.x = element_text(margin = margin(0.6, 0, 0.6, 0, "cm")), 
            strip.background = element_rect(colour = 'black', fill = NA, size = 3.25), 
            axis.ticks.length = unit(0.45, "cm"), 
            axis.ticks = element_line(colour = "black", size = 1.7), 
            axis.line = element_line(colour = 'black', size = 1.7), 
            plot.margin = unit(c(0.5, 0.5, 0.5, 10.55),"cm"), 
            axis.title.y = element_text(size=36.25, margin = margin(t = 0, r = 18.5, b = 0, l = 2), colour="black", 
                face = "plain"), 
            axis.title.x = element_blank(), 
            axis.text.x = element_text(size=36.25, angle=45, margin = margin(t=-44, b=59), colour="black", 
                hjust = 0.99, vjust = 0.5), 
            axis.text.y = element_text(size=36.25, angle=0, margin = margin(r = 5), colour="black"), 
            panel.spacing = unit(0.5, "cm"), 
            panel.grid.major = element_blank(),
            panel.grid.minor.x = element_blank(), 
            panel.grid.minor.y = element_blank(), 
            legend.position = "none", 
            legend.title = element_blank())  

        ggsave(file = file.path(out_dir, "output", "plots", fname), plot = q, 
            width = 8.25, height = 12, dpi = 300, units = c("in"), limitsize = FALSE) 
    }

    plotTTL.SO(data = bs_tree_length_stamen_pollen[bs_tree_length_stamen_pollen$Order == "Angiosperms", ], 
        data2 = observed_tree_length_stamen_pollen[observed_tree_length_stamen_pollen$Order == "Angiosperms", ], 
        spec_order = "AS")
    plotTTL.SO(data = bs_tree_length_stamen_pollen[bs_tree_length_stamen_pollen$Order == "Brassicaceae", ], 
        data2 = observed_tree_length_stamen_pollen[observed_tree_length_stamen_pollen$Order == "Brassicaceae", ], 
        spec_order = "BR")


   
}


plotTreeLength()





