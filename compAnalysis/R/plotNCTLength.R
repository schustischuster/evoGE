# Plot total tree length of lncRNA and Brassicaceae protein-coding genes for SI
# Data input: Data tables w/ total tree length and bootstrap values


#-------------------------------------- Read data tables ---------------------------------------


plotNCTLength <- function() {

    # Show startup message
    message("Reading data...")


    Brass_nc_bs_tree_length_spe = file.path(out_dir, "output", "data", "Brass_non-coding_total_tree_length_spearman_counts.txt")
    Brass_cd_bs_tree_length_spe = file.path(out_dir, "output", "data", "Brass_coding_total_tree_length_spearman_counts.txt")


    # Read expression data
    Brass_nc_bs_tree_length <- read.table(Brass_nc_bs_tree_length_spe, header=TRUE, sep="\t", dec=".", stringsAsFactors=FALSE)
    Brass_cd_bs_tree_length <- read.table(Brass_cd_bs_tree_length_spe, header=TRUE, sep="\t", dec=".", stringsAsFactors=FALSE)




#--------------------- Prepare data and define color palette for corrplot ---------------------


    # Create "plots" folder in /out_dir/output/plots
    if (!dir.exists(file.path(out_dir, "output", "plots"))) 
        dir.create(file.path(out_dir, "output", "plots"), recursive = TRUE)

    # Show message
    message("Starting analysis and generate plots...")


    # Get SD based on bootstrapping analysis
    getSTD <- function(df) {

        n <- nrow(df)/9
        nr <- nrow(df)
        STD <- do.call(rbind, lapply(split(df, rep(1:ceiling(nr/n), each=n, length.out=nr)), 
            function(x) {

                std <- sd(x$tree_length, na.rm=TRUE)
                num <- nrow(x)
                std <- data.frame(rep(std, num))
                return(std)
            }))

        colnames(STD) <- "sd"

        return(STD)
    }

    Brass_nc_bs_tree_length_SD <- getSTD(Brass_nc_bs_tree_length)
    Brass_nc_bs_sd_tree_length <- cbind(Brass_nc_bs_tree_length, Brass_nc_bs_tree_length_SD)

    Brass_cd_bs_tree_length_SD <- getSTD(Brass_cd_bs_tree_length)
    Brass_cd_bs_sd_tree_length <- cbind(Brass_cd_bs_tree_length, Brass_cd_bs_tree_length_SD)



    # Prepare data for ggplot2
    Brass_nc_bs_sd_tree_length$organ <- factor(Brass_nc_bs_sd_tree_length$organ, levels = unique(Brass_nc_bs_sd_tree_length$organ))

    # Separate pollen samples of Brassicaceae coding gene table
    Brass_cd_bs_tree_length_wo_pollen <- dplyr::filter(Brass_cd_bs_sd_tree_length, !grepl('Pollen', organ))
    Brass_cd_bs_tree_length_w_pollen <- dplyr::filter(Brass_cd_bs_sd_tree_length, grepl('Stamen|Pollen', organ))

    Brass_cd_bs_tree_length_wo_pollen$organ <- factor(Brass_cd_bs_tree_length_wo_pollen$organ, 
        levels = unique(Brass_cd_bs_tree_length_wo_pollen$organ))
    Brass_cd_bs_tree_length_w_pollen$organ <- factor(Brass_cd_bs_tree_length_w_pollen$organ, 
        levels = unique(Brass_cd_bs_tree_length_w_pollen$organ))



    # Plot lncRNA ortholog total tree length
    plotTTL <- function(data) {

        obj_name <- deparse(substitute(data))

        if (obj_name == "Brass_nc_bs_tree_length") {

            level_order <- c('Root', 'Hypocotyl', 'Leaf', 'veg_apex', 'inf_apex', 'Flower', 
                'Carpel', 'Stamen', 'Pollen')
            org_colors <- c('#6a54a9', '#53b0db', '#2c8654', '#96ba37', '#fad819', '#e075af', 
                '#f2a72f', '#ed311c', '#a63126')
            org_labels <- c("Root" = "Root", "Hypocotyl" = "Hypocotyl", "Leaf" = "Leaf", 
                "veg_apex" = "Apex veg", "inf_apex" = "Apex inf", "Flower" = "Flower", 
                "Carpel" = "Carpel", "Stamen" = "Stamen", "Pollen" = "Pollen")
            plt_margin <- unit(c(1.5, 2.0, 1.31, 1.5), "cm")
            y_exp <- c(0.005, 0)
            x_exp <- c(0.035, 0)
            txt_x_mar <- margin(t = -38, b = 42.5)

        } else if (obj_name == "Brass_cd_bs_tree_length_wo_pollen") {

            level_order <- c('Root', 'Hypocotyl', 'Leaf', 'veg_apex', 'inf_apex', 'Flower', 
                'Carpel', 'Stamen')
            org_colors <- c('#6a54a9', '#53b0db', '#2c8654', '#96ba37', '#fad819', '#e075af', 
                '#f2a72f', '#ed311c')
            org_labels <- c("Root" = "Root", "Hypocotyl" = "Hypocotyl", "Leaf" = "Leaf", 
                "veg_apex" = "Apex veg", "inf_apex" = "Apex inf", "Flower" = "Flower", 
                "Carpel" = "Carpel", "Stamen" = "Stamen")
            plt_margin <- unit(c(1.5, 2.0, 1.31, 3.0), "cm")
            y_exp <- c(0.005, 0)
            x_exp <- c(0.04, 0)
            txt_x_mar <- margin(t = -38, b = 42.5)

        } else if (obj_name == "Brass_cd_bs_tree_length_w_pollen") {

            level_order <- c('Stamen', 'Pollen')
            org_colors <- c('#ed311c', '#a63126')
            org_labels <- c("Stamen" = "Stamen", "Pollen" = "Pollen")
            plt_margin <- unit(c(1.5, 2.0, 1.62, 12.5), "cm")
            y_exp <- c(0.005, 0)
            x_exp <- c(0.035, 0)
            txt_x_mar <- margin(t = -29.5, b = 42.5)

        }

        fname <- sprintf('%s.jpg', paste(obj_name))
            
        p <- ggplot(data=data, aes(x=factor(organ, level=level_order), y=tree_length, group=organ)) + 
        stat_boxplot(geom ='errorbar', width = 0.35, size=1.25, color="black") + 
        geom_boxplot(width = 0.8, size = 1.25, fatten = 1.35, color="black", outlier.shape = NA, 
            alpha = 1, aes(fill = organ)) + 
        geom_point(data = data, position = position_dodge(width=0.75), size = 5, col = "red2", 
            aes(x = organ, y = orig_tree_length)) + 
        scale_fill_manual(values = org_colors) + 
        scale_y_continuous(expand = y_exp) + 
        scale_x_discrete(labels = org_labels, expand = x_exp)

        q <- p + theme_classic() + ylab("Total tree length") + 
        theme(text=element_text(size = 16), 
            axis.ticks.length = unit(0.35, "cm"), 
            axis.ticks = element_line(colour = "black", size = 1.15), 
            axis.line = element_line(colour = 'black', size = 1.15), 
            plot.margin = plt_margin, 
            axis.title.y = element_text(size=27.0, margin = margin(t = 0, r = 17.0, b = 0, l = 3.5), colour="black", 
                face = "plain"), 
            axis.title.x = element_blank(), 
            axis.text.x = element_text(size=25.75, angle=45, margin = txt_x_mar, colour="black", 
                hjust = 1, vjust = 0.5), 
            axis.text.y = element_text(size=24.5, angle=0, margin = margin(r = 5), colour="black"), 
            panel.spacing = unit(0.5, "cm"), 
            panel.grid.major = element_blank(),
            panel.grid.minor.x = element_blank(), 
            panel.grid.minor.y = element_blank(), 
            legend.position = "none", 
            legend.title = element_blank()) 

        ggsave(file = file.path(out_dir, "output", "plots", fname), plot = q, 
            width = 10.0, height = 9.1, dpi = 300, units = c("in"), limitsize = FALSE) 
    }

    plotTTL(data = Brass_nc_bs_tree_length)
    plotTTL(data = Brass_cd_bs_tree_length_wo_pollen)
    plotTTL(data = Brass_cd_bs_tree_length_w_pollen)

   
}


plotNCTLength()





