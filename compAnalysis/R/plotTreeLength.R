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
    bs_tree_length$Total_tree_length <- as.numeric(bs_tree_length$Total_tree_length)
    bs_tree_length$Order <- factor(bs_tree_length$Order)


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
    observed_tree_length$Organ <- factor(observed_tree_length$Organ)
    observed_tree_length$Total_tree_length <- factor(observed_tree_length$Total_tree_length)
    observed_tree_length$Order <- factor(observed_tree_length$Order)


    # Create data tables w/ and wo/ pollen
    observed_tree_length_wo_pollen <- observed_tree_length[observed_tree_length$Organ != "Pollen", ]
    observed_tree_length_pollen <- observed_tree_length[observed_tree_length$Organ == "Pollen", ]
    observed_tree_length_stamen <- observed_tree_length[observed_tree_length$Organ == "Stamen", ]
    observed_tree_length_stamen_pollen <- rbind(observed_tree_length_stamen, observed_tree_length_pollen)
    observed_tree_length_wo_pollen$Total_tree_length <- as.numeric(as.character(observed_tree_length_wo_pollen$Total_tree_length))
    observed_tree_length_stamen_pollen$Total_tree_length <- as.numeric(as.character(observed_tree_length_stamen_pollen$Total_tree_length))

    bs_tree_length_wo_pollen$Organ <- factor(bs_tree_length_wo_pollen$Organ, levels = unique(bs_tree_length_wo_pollen$Organ))
    observed_tree_length_wo_pollen$Organ <- factor(observed_tree_length_wo_pollen$Organ, levels = unique(observed_tree_length_wo_pollen$Organ))
    bs_tree_length_stamen_pollen$Organ <- factor(bs_tree_length_stamen_pollen$Organ, levels = unique(bs_tree_length_stamen_pollen$Organ))
    observed_tree_length_stamen_pollen$Organ <- factor(observed_tree_length_stamen_pollen$Organ, levels = unique(observed_tree_length_stamen_pollen$Organ))

# Mat="#a63126"

    plotTTL <- function(data, data2) {

        fname <- sprintf('%s.jpg', paste("Total_tree_length_AT"))
            
        p <- ggplot(data=data, aes(x = Organ, y = Total_tree_length, group = Organ)) + 
        stat_boxplot(geom ='errorbar', width = 0.5, size=1.0, color="black") + 
        geom_boxplot(width = 0.85, size = 1.25, fatten = 2, color="black", outlier.shape = NA, 
            alpha = 1, aes(fill = Organ)) + 
        geom_point(data = data2, position = position_dodge(width=0.75), size = 5, col = "red2", 
            aes(x = Organ, y = Total_tree_length)) + 
        scale_fill_manual(values=c('#6a54a9','#53b0db', '#2c8654', '#96ba37', '#fad819','#e075af', '#f2a72f', '#ed311c')) + 
        scale_y_continuous(expand = c(0.025, 0)) + 
        scale_x_discrete(labels=c("Root" = "Root", "Hypocotyl" = "Hypocotyl", 
            "Leaf" = "Leaf", "veg_apex" = "Apex veg", "inf_apex" = "Apex inf", 
            "Flower" = "Flower", "Carpel" = "Carpel", "Stamen" = "Stamen"))

        q <- p + theme_classic() + xlab("Organ") + ylab("Total tree length") + 
        theme(text=element_text(size = 16), 
            strip.text = element_text(size = 29.5), 
            strip.text.x = element_text(margin = margin(0.58, 0, 0.58, 0, "cm")), 
            strip.background = element_rect(colour = 'black', fill = NA, size = 2.5), 
            axis.ticks.length = unit(0.45, "cm"), 
            axis.ticks = element_line(colour = "black", size = 1.25), 
            axis.line = element_line(colour = 'black', size = 1.25), 
            plot.margin = unit(c(0.85, 1.175, 2, 0.4),"cm"), 
            axis.title.y = element_text(size=30, margin = margin(t = 0, r = 15.2, b = 0, l = 10.8), colour="black", 
                face = "bold"), 
            axis.title.x = element_text(size=30, margin = margin(t = 10, r = 0, b = -15, l = 0), colour="black", 
                face = "bold"), 
            axis.text.x = element_text(size=29.5, angle=45, margin = margin(t = -48, b = 48), colour="black", 
                hjust = 0.99, vjust = 0.45), 
            axis.text.y = element_text(size=29, angle=0, margin = margin(r = 5.5), colour="black"), 
            panel.spacing = unit(0.5, "cm"), 
            panel.grid.major = element_blank(),
            panel.grid.minor.x = element_blank(), 
            panel.grid.minor.y = element_blank(), 
            legend.position = "none", 
            legend.title = element_blank()) 

        q <- q + facet_wrap(~ Order, scales = "free", nrow = 1)

        ggsave(file = file.path(out_dir, "output", "plots", fname), plot = q, 
            width = 15, height = 13.5, dpi = 300, units = c("in"), limitsize = FALSE) 
    }

    plotTTL(data = bs_tree_length_wo_pollen, data2 = observed_tree_length_wo_pollen)



    plotTTL.SO <- function(data, data2) {

        fname <- sprintf('%s.jpg', paste("Total_tree_length_AT_pollen"))
            
        p <- ggplot(data=data, aes(x = Organ, y = Total_tree_length, group = Organ)) + 
        stat_boxplot(geom ='errorbar', width = 0.5, size=1.0, color="black") + 
        geom_boxplot(width = 0.85, size = 1.25, fatten = 2, color="black", outlier.shape = NA, 
            alpha = 1, aes(fill = Organ)) + 
        geom_point(data = data2, position = position_dodge(width=0.75), size = 5, col = "red2", 
            aes(x = Organ, y = Total_tree_length)) + 
        scale_fill_manual(values=c('#ed311c', '#a63126')) + 
        scale_y_continuous(expand = c(0.025, 0)) + 
        scale_x_discrete(labels=c("Stamen" = "Stamen", "Pollen" = "Pollen"))

        q <- p + theme_classic() + xlab("Organ") + ylab("Total tree length") + 
        theme(text=element_text(size = 16), 
            strip.text = element_text(size = 29.5), 
            strip.text.x = element_text(margin = margin(0.58, 0, 0.58, 0, "cm")), 
            strip.background = element_rect(colour = 'black', fill = NA, size = 2.5), 
            axis.ticks.length = unit(0.45, "cm"), 
            axis.ticks = element_line(colour = "black", size = 1.25), 
            axis.line = element_line(colour = 'black', size = 1.25), 
            plot.margin = unit(c(0.85, 1.175, 2, 0.4),"cm"), 
            axis.title.y = element_text(size=30, margin = margin(t = 0, r = 15.2, b = 0, l = 10.8), colour="black", 
                face = "bold"), 
            axis.title.x = element_text(size=30, margin = margin(t = 10, r = 0, b = -15, l = 0), colour="black", 
                face = "bold"), 
            axis.text.x = element_text(size=29.5, angle=45, margin = margin(t = -37.5, b = 57.5), colour="black", 
                hjust = 0.99, vjust = 0.45), 
            axis.text.y = element_text(size=29, angle=0, margin = margin(r = 5.5), colour="black"), 
            panel.spacing = unit(0.5, "cm"), 
            panel.grid.major = element_blank(),
            panel.grid.minor.x = element_blank(), 
            panel.grid.minor.y = element_blank(), 
            legend.position = "none", 
            legend.title = element_blank()) 

        q <- q + facet_wrap(~ Order, scales = "free", nrow = 1)

        ggsave(file = file.path(out_dir, "output", "plots", fname), plot = q, 
            width = 10.5, height = 13.5, dpi = 300, units = c("in"), limitsize = FALSE) 
    }

    plotTTL.SO(data = bs_tree_length_stamen_pollen, data2 = observed_tree_length_stamen_pollen)


   
}


plotTreeLength()





