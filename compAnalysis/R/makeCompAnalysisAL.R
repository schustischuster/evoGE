# Prepare Brawand and DevSeq AL comparative ortholog gene expression data
# Thresholds: DevSeq 0.05 ERCC; Brawand 0.5 TPM (no ERCC spike-ins available)
# Data input: Brawand and DevSeq TPM expression tables of all samples



#-------------------------------------- Read data tables ---------------------------------------


makeCompAnalysisAL <- function(expr_estimation = c("TPM", "counts"), coefficient = c("pearson", "spearman")) {


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


    # Show startup message
    message("Reading data...")


    if (is.element("TPM", expr_estimation)) {
        genesExpr = file.path(in_dir, "Expression_data", "AL_core_inter_tpm_mat_deseq_sample_names.csv")

    } else if (is.element("counts", expr_estimation)) {
        genesExpr = file.path(in_dir, "Expression_data", "AL_core_inter_count_mat_vsd_sample_names.csv")
    }


	  # Define simplified DevSeq column names
    col_names <- rep(c("Root", "Hypocotyl", "Leaf", "veg_apex", "inf_apex", 
        "Flower", "Stamen", "Carpel", "Pollen"), each=21)
    replicate_tag_samples <- rep(c(".1",".2",".3"), times=9)
    col_names <- paste0(col_names,replicate_tag_samples)
    spec_names <- rep(c("_AL", "_AT", "_CR", "_ES", "_TH", "_MT", "_BD"), each=3)
    spec_names <- rep(spec_names, times=9)
    col_names <- paste0(col_names, spec_names)
    col_names <- c("gene_id", col_names)


	  # Read expression data
	  x <- read.table(genesExpr, sep=";", dec=".", header=TRUE, stringsAsFactors=FALSE)


    # Read Brawand11 metric Pearson distance data and merge with AL pea dist data
    div_rates_list <- c("compDivRates", "compDivRates11", "compSouVDivRates", "compSouVDivRates11")

    for(i in 1:length(div_rates_list)){

      path_to_br11_div_rates <- file.path(out_dir, "output", "data", paste0(div_rates_list[i],".txt"))

      assign(div_rates_list[i], read.table(path_to_br11_div_rates, header=TRUE, sep="\t", dec=".", 
        stringsAsFactors=FALSE))
    }


    # Stop function here to allow specific analysis of a single data set
    # return_list <- list("expr_estimation" = expr_estimation, "x" = x, "coefficient" = coefficient, "col_names" = col_names, 
    # "compDivRates" = compDivRates, "compDivRates11" = compDivRates11, "compSouVDivRates" = compSouVDivRates, "compSouVDivRates11" = compSouVDivRates11)
    # return(return_list)
    # }
    # return_objects <- makeCompAnalysisAL(expr_estimation="TPM", coefficient="pearson") # read in DevSeq expression data
    # list2env(return_objects, envir = .GlobalEnv)


    # set column names
    colnames(x) <- col_names




#--------------------- Prepare data and define color palette for corrplot ---------------------


    # Create "plots" folder in /out_dir/output/plots
    if (!dir.exists(file.path(out_dir, "output", "plots"))) 
        dir.create(file.path(out_dir, "output", "plots"), recursive = TRUE)

    # Show message
    message("Starting analysis and generate plots...")


    x[is.na(x)] <- 0 # replaces NAs by 0
    # Remove ERCC spike-ins from data
    x <- x[!grepl("ERCC", x$gene_id),]




#---------------- Get gene expression divergence rates for ATH/AL vs species X -----------------


   # Use pearson correlation, intra-organ normalization and TPM
   # Use previously merged replicates of DevSeq data including pollen sampless

   if (expr_estimation == "TPM") {

      getOrganCor <- function(df, organ, coefficient, expr_estimation) {

         # log-transform data if TPM and Pearson are chosen
         if ((coefficient == "pearson") && (expr_estimation == "TPM")) {
            df <- log2(df + 1)
         }

         df_cor <- cor(df, method=coefficient)
         df_cor <- df_cor[4:nrow(df_cor), 1:3]

         getError <- function(cor_data) {
            std <- sd(cor_data)
            n_value <- length(cor_data)
            error <- qt(0.995, df = n_value-1) * std/sqrt(n_value)
            return(error)
         }

         sp1 <- mean(c(df_cor[1:3,]))
         sp1_li <- sp1 - getError(c(df_cor[1:3,]))
         sp1_ri <- sp1 + getError(c(df_cor[1:3,]))

         sp2 <- mean(c(df_cor[4:6,]))
         sp2_li <- sp2 - getError(c(df_cor[4:6,]))
         sp2_ri <- sp2 + getError(c(df_cor[4:6,]))

         sp3 <- mean(c(df_cor[7:9,]))
         sp3_li <- sp3 - getError(c(df_cor[7:9,]))
         sp3_ri <- sp3 + getError(c(df_cor[7:9,]))

         sp4 <- mean(c(df_cor[10:12,]))
         sp4_li <- sp4 - getError(c(df_cor[10:12,]))
         sp4_ri <- sp4 + getError(c(df_cor[10:12,]))

         sp5 <- mean(c(df_cor[13:15,]))
         sp5_li <- sp5 - getError(c(df_cor[13:15,]))
         sp5_ri <- sp5 + getError(c(df_cor[13:15,]))

         sp6 <- mean(c(df_cor[16:18,]))
         sp6_li <- sp6 - getError(c(df_cor[16:18,]))
         sp6_ri <- sp6 + getError(c(df_cor[16:18,]))

         df_cor_avg <- rbind(sp1, sp2, sp3, sp4, sp5, sp6)
         colnames(df_cor_avg) <- organ
         lower <- rbind(sp1_li, sp2_li, sp3_li, sp4_li, sp5_li, sp6_li)
         colnames(lower) <- "lower"
         upper <- rbind(sp1_ri, sp2_ri, sp3_ri, sp4_ri, sp5_ri, sp6_ri)
         colnames(upper) <- "upper"

         getRowNames = function(x,n){ substring(x,nchar(x)-n+1) }
         row_names_repl <- getRowNames(rownames(df_cor),2)
         rnames_div_rates <- unique(row_names_repl)
         rownames(df_cor_avg) <- rnames_div_rates
         df_cor_avg <- cbind(df_cor_avg, lower, upper)

         return(df_cor_avg)

      }

      root_div <- getOrganCor(df=x[,2:22], organ="Root  ", coefficient=coefficient, expr_estimation=expr_estimation)
      hypocotyl_div <- getOrganCor(df=x[,23:43], organ="Hypocotyl  ", coefficient=coefficient, expr_estimation=expr_estimation)
      leaf_div <- getOrganCor(df=x[,44:64], organ="Leaf  ", coefficient=coefficient, expr_estimation=expr_estimation)
      veg_apex_div <- getOrganCor(df=x[,65:85], organ="Apex veg  ", coefficient=coefficient, expr_estimation=expr_estimation)
      inf_apex_div <- getOrganCor(df=x[,86:106], organ="Apex inf  ", coefficient=coefficient, expr_estimation=expr_estimation)
      flower_div <- getOrganCor(df=x[,107:127], organ="Flower  ", coefficient=coefficient, expr_estimation=expr_estimation)
      stamen_div <- getOrganCor(df=x[,128:148], organ="Stamen  ", coefficient=coefficient, expr_estimation=expr_estimation)
      carpel_div <- getOrganCor(df=x[,149:169], organ="Carpel  ", coefficient=coefficient, expr_estimation=expr_estimation)
      pollen_div <- getOrganCor(df=x[,170:190], organ="Pollen  ", coefficient=coefficient, expr_estimation=expr_estimation)


      # Reshape data table for ggplot
      # divergence times are estimated taxon pair times from TimeTree
      # http://www.timetree.org/
      div_times <- rep(c(7.1, 9.4, 25.6, 46, 106, 160), times=9)
      comp_organ <- rep(c(colnames(root_div)[1], colnames(hypocotyl_div)[1], colnames(leaf_div)[1], 
        colnames(veg_apex_div)[1], colnames(inf_apex_div)[1], colnames(flower_div)[1], 
        colnames(stamen_div)[1], colnames(carpel_div)[1], colnames(pollen_div)[1]), each=6)
      comp_spec <- c(rownames(root_div), rownames(hypocotyl_div), rownames(leaf_div), rownames(veg_apex_div), 
        rownames(inf_apex_div), rownames(flower_div), rownames(stamen_div), rownames(carpel_div), 
        rownames(pollen_div))

      DevSeq_AL_GE_div <- rbind(root_div, hypocotyl_div, leaf_div, veg_apex_div, inf_apex_div, 
        flower_div, stamen_div, carpel_div, pollen_div)
      rownames(DevSeq_AL_GE_div) <- NULL
      colnames(DevSeq_AL_GE_div) <- c("correlation", "lower", "upper")

      DevSeq_AL_div <- data.frame(cbind(comp_spec, comp_organ, div_times, DevSeq_AL_GE_div), 
        stringsAsFactors=FALSE)

      DevSeq_AL_div$div_times <- as.numeric(DevSeq_AL_div$div_times)
      DevSeq_AL_div$correlation <- as.numeric(DevSeq_AL_div$correlation)
      DevSeq_AL_div$lower <- as.numeric(DevSeq_AL_div$lower)
      DevSeq_AL_div$upper <- as.numeric(DevSeq_AL_div$upper)
      
      # Change order of organs in df
      DevSeq_AL_div <- DevSeq_AL_div[c(7:12,37:42,31:36,1:6,19:30,43:48,13:18,49:54),]
      DevSeq_AL_div_rates <- DevSeq_AL_div[1:54,]
      DevSeq_AL_div_rates_pollen <- DevSeq_AL_div[49:54,]
      DevSeq_AL_div_rates$comp_organ <- factor(DevSeq_AL_div_rates$comp_organ, 
        levels = unique(DevSeq_AL_div_rates$comp_organ))
      DevSeq_AL_div_rates_pollen$comp_organ <- factor(DevSeq_AL_div_rates_pollen$comp_organ, 
        levels = unique(DevSeq_AL_div_rates_pollen$comp_organ))

      

      # Make GE divergence plot
      makeGEDivPlot <- function(data1, data2, plot_title, coefficient) {

        fname <- sprintf('%s.jpg', paste("GE_divergence_rates_AL", coefficient, expr_estimation, "ext", sep="_"))
        plot_wdt <- 9.5 # condenced plot width for suppl
        plot_hdt <- 6.75 # condenced plot width for suppl
        legend_x_pos <- 0.638
        legend_y_pos <- 0.862
        legend_key_s <- 0.9
        linewd <- 2.6
        brass_label <- "Brassiceae"
        text_x1_pos <- 22
        text_x4_pos <- 157
        text_y_pos <- 0.46525
        x_poll_pos <- 145.8
        y_poll_pos <- 0.81525
        leg_ln_s1_x <- 129.45
        leg_ln_s1_xend <- 131.7
        leg_ln_s2_x <- 133.15
        leg_ln_s2_xend <- 135.4
        leg_ln_s_y <- 0.8155
        leg_box_bd <- 0
        pan_boarder <- 1.8
        axis_txt_size <- 21.75
        axis_ticks_s <- 0.95
        title_face <- "bold"

        p <- ggplot(data=data1, aes(x=div_times, y=correlation, group=comp_organ, colour=comp_organ)) + 
        geom_ribbon(aes(ymin = data1$lower, ymax = data1$upper, fill= comp_organ), alpha = 0.25, 
            linetype = 0, show.legend = FALSE) + 
        scale_fill_manual(values = c("Hypocotyl  "="#53b0db", "Stamen  "="#ee412e", "Flower  "="#e075af", 
                "Root  "="#6a54a9", "Apex veg  "="#96ba37", "Apex inf  "="#fad819", "Carpel  "="#f2a72f", 
                "Leaf  "="#2c8654", "Pollen  "="#a63126")) + 
        geom_line(size = linewd) +  
        scale_x_continuous(limits = c(5.5,161.5), expand = c(0.02,0), breaks = c(7,9,25,46,106,160)) + 
        scale_y_continuous(limits = c(0.4425, 0.907), expand = c(0.02, 0)) + 
        scale_color_manual(values = c("#53b0db", "#ee412e", "#e075af", "#6a54a9", "#96ba37", "#fad819", 
            "#f2a72f", "#2c8654", "#a63126"), 
            # organ order: hypocotyl/stamen/flower/root/veg_apex/inf_apex/carpel/leaf
            breaks=c("Root  ", "Hypocotyl  ", "Leaf  ", "Apex veg  ", "Apex inf  ", "Flower  ", 
                "Stamen  ", "Carpel  ", "Pollen  ")) + 
        geom_line(aes(x=div_times, y=correlation), data=data2, color = "white", lty = "solid", 
            lwd = linewd) + # pollen
        geom_line(aes(x=div_times, y=correlation), data=data2, color = "#a63126", lty = "22", 
            lwd = linewd) + # pollen
        geom_segment(x=156.75, xend=156.75, y=0.4475, yend=0.4775, color="white", size=12.5) + 
        annotate("text", x=text_x1_pos, y=text_y_pos, label= brass_label, size=8) + 
        annotate("text", x=46, y=text_y_pos, label= "TH", size=8) + 
        annotate("text", x=106, y=text_y_pos, label= "MT", size=8) + 
        annotate("text", x=text_x4_pos, y=text_y_pos, label= "BD", size=8) + 
        geom_segment(x=7, xend=7, y=0.425, yend=0.4475, color="black", size=axis_ticks_s-0.05) + 
        geom_segment(x=9, xend=9, y=0.425, yend=0.4475, color="black", size=axis_ticks_s-0.05) + 
        geom_segment(x=25, xend=25, y=0.425, yend=0.4475, color="black", size=axis_ticks_s-0.05) + 
        geom_segment(x=46, xend=46, y=0.425, yend=0.4475, color="black", size=axis_ticks_s-0.05) + 
        geom_segment(x=106, xend=106, y=0.425, yend=0.4475, color="black", size=axis_ticks_s-0.05) + 
        geom_segment(x=160, xend=160, y=0.425, yend=0.4475, color="black", size=axis_ticks_s-0.05) + 
        guides(color = guide_legend(ncol = 3))

        q <- p + theme_bw() + xlab("Divergence time from A.lyrata (Myr)") + ylab("Pearson's r w/ A.lyrata") + 
        theme(text=element_text(size=16), 
            axis.ticks.length=unit(0.35, "cm"), 
            axis.ticks = element_line(colour = "black", size = axis_ticks_s),  
            plot.margin = unit(c(0.55, 1.175, 0.5, 0.4),"cm"), 
            axis.title.y = element_text(size=24.5, margin = margin(t = 0, r = 15, b = 0, l = 11), colour="black", 
                face = title_face), 
            axis.title.x = element_text(size=24.5, margin = margin(t = 15.75, r = 0, b = 1, l = 0), colour="black", 
                face = title_face), 
            axis.text.x = element_text(size=axis_txt_size, angle=0, margin = margin(t = 5.5), colour="black"), 
            axis.text.y = element_text(size=axis_txt_size, angle=0, margin = margin(r = 5.5), colour="black"), 
            legend.box.background = element_rect(colour = "#d5d5d5", fill=NA, size=leg_box_bd), 
            panel.border = element_rect(colour = "black", fill=NA, size=pan_boarder), 
            panel.grid.major = element_blank(),
            panel.grid.minor.x = element_blank(), 
            panel.grid.minor.y = element_blank(), 
            legend.position = c(legend_x_pos, legend_y_pos), 
            legend.title = element_blank(), 
            legend.text = element_text(size=21.5), 
            legend.spacing.x = unit(0.5, 'cm'), 
            legend.key.size = unit(legend_key_s, "cm"), 
            legend.background=element_blank()) 

        ggsave(file = file.path(out_dir, "output", "plots", fname), plot = q, 
            width = plot_wdt, height = plot_hdt, dpi = 300, units = c("in"), limitsize = FALSE) 
      }

      makeGEDivPlot(data1 = DevSeq_AL_div_rates, data2 = DevSeq_AL_div_rates_pollen, 
        coefficient = coefficient)




#---- Apply non-linear regression to sOU and pearson dist expression data and compare slopes -----

# Non-linear regression using negative exponential law fit: pairwise expression differences
# between species saturate with evolutionary time in a power law relationship
# Fits assumption of OU model underlying stabilizing GE selection as a decelarated process


      nl_model <- function(a, b, c, x){

        y = a * exp(c * x) + b * (1 - exp(c * x))
        return(y)
      }
      # a + b defines maximum y value
      # a defines intercept


      x_DS_grid <- seq(7.1, 160, length = 200)  ## prediction grid
      x_Br_grid <- seq(6.7, 159, length = 200)  ## prediction grid

      # Compute data points for DevSeq_AL_pearson_dist based on model
      # First try to manually find rough parameters, then use nls to fine tune
      m <- nls(correlation ~ a * exp(div_times * c) + b * (1-(exp(div_times * c))), start = list(
        a = 0.3, b = 0.5, c = -0.01), data = compDivRates11_AL[1:6,])
      # m # get the optimized parameters


      # Get fit for data from compDivRates11 (AT)
      DS_AT_pea_dist_root_nl <- as.data.frame(do.call(rbind, lapply(x_DS_grid, 
        nl_model, a = 0.22559, b = 0.43214, c = -0.02308))) # compDivRates11[1:6, ]
      DS_AT_pea_dist_hypo_nl <- as.data.frame(do.call(rbind, lapply(x_DS_grid, 
        nl_model, a = 0.23702, b = 0.48283, c = -0.01718))) # compDivRates11[7:12, ]
      DS_AT_pea_dist_leaf_nl <- as.data.frame(do.call(rbind, lapply(x_DS_grid, 
        nl_model, a = 0.21827, b = 0.46071, c = -0.01101))) # compDivRates11[13:18, ]
      DS_AT_pea_dist_apex_veg_nl <- as.data.frame(do.call(rbind, lapply(x_DS_grid, 
        nl_model, a = 0.21392, b = 0.42240, c = -0.02446))) # compDivRates11[19:24, ]
      DS_AT_pea_dist_apex_inf_nl <- as.data.frame(do.call(rbind, lapply(x_DS_grid, 
        nl_model, a = 0.21304, b = 0.41412, c = -0.02389))) # compDivRates11[25:30, ]
      DS_AT_pea_dist_flower_nl <- as.data.frame(do.call(rbind, lapply(x_DS_grid, 
        nl_model, a = 0.21738, b = 0.43576, c = -0.02739))) # compDivRates11[31:36, ]
      DS_AT_pea_dist_stamen_nl <- as.data.frame(do.call(rbind, lapply(x_DS_grid, 
        nl_model, a = 0.23374, b = 0.47328, c = -0.02325))) # compDivRates11[37:42, ]
      DS_AT_pea_dist_carpel_nl <- as.data.frame(do.call(rbind, lapply(x_DS_grid, 
        nl_model, a = 0.2123, b = 0.4342, c = -0.0208))) # compDivRates11[43:48, ]

      DS_AT_pea_dist_nl_list <- list(DS_AT_pea_dist_root_nl=DS_AT_pea_dist_root_nl,
        DS_AT_pea_dist_hypo_nl=DS_AT_pea_dist_hypo_nl, DS_AT_pea_dist_leaf_nl=DS_AT_pea_dist_leaf_nl, 
        DS_AT_pea_dist_apex_veg_nl=DS_AT_pea_dist_apex_veg_nl, DS_AT_pea_dist_apex_inf_nl=DS_AT_pea_dist_apex_inf_nl, 
        DS_AT_pea_dist_flower_nl=DS_AT_pea_dist_flower_nl, DS_AT_pea_dist_stamen_nl=DS_AT_pea_dist_stamen_nl, 
        DS_AT_pea_dist_carpel_nl=DS_AT_pea_dist_carpel_nl)


      # Get fit for data from compSouVDivRates11 (AT)
      DS_AT_sOU_v_root_nl <- as.data.frame(do.call(rbind, lapply(x_DS_grid, nl_model, 
        a = 0.13873, b = 1.19261, c = -0.01683))) # compSouVDivRates11[1:6, ]
      DS_AT_sOU_v_hypo_nl <- as.data.frame(do.call(rbind, lapply(x_DS_grid, nl_model, 
        a = 0.233979, b = 4.567626, c = -0.002489))) # compSouVDivRates11[7:12, ]
      DS_AT_sOU_v_leaf_nl <- as.data.frame(do.call(rbind, lapply(x_DS_grid, nl_model, 
        a = 0.178, b = 8.698, c = -0.000758))) # compSouVDivRates11[13:18, ]
      DS_AT_sOU_v_apex_veg_nl <- as.data.frame(do.call(rbind, lapply(x_DS_grid, nl_model, 
        a = 0.13038, b = 1.11839, c = -0.01728))) # compSouVDivRates11[19:24, ]
      DS_AT_sOU_v_apex_inf_nl <- as.data.frame(do.call(rbind, lapply(x_DS_grid, nl_model, 
        a = 0.14499, b = 1.07888, c = -0.01609))) # compSouVDivRates11[25:30, ]
      DS_AT_sOU_v_flower_nl <- as.data.frame(do.call(rbind, lapply(x_DS_grid, nl_model, 
        a = 0.15837, b = 1.34561, c = -0.01447))) # compSouVDivRates11[31:36, ]
      DS_AT_sOU_v_stamen_nl <- as.data.frame(do.call(rbind, lapply(x_DS_grid, nl_model, 
        a = 0.13809, b = 1.61684, c = -0.01508))) # compSouVDivRates11[37:42, ]
      DS_AT_sOU_v_carpel_nl <- as.data.frame(do.call(rbind, lapply(x_DS_grid, nl_model, 
        a = 0.165474, b = 1.489164, c = -0.008982))) # compSouVDivRates11[43:48, ]

      DS_AT_sOU_v_nl_list <- list(DS_AT_sOU_v_root_nl=DS_AT_sOU_v_root_nl,
        DS_AT_sOU_v_hypo_nl=DS_AT_sOU_v_hypo_nl, DS_AT_sOU_v_leaf_nl=DS_AT_sOU_v_leaf_nl, 
        DS_AT_sOU_v_apex_veg_nl=DS_AT_sOU_v_apex_veg_nl, DS_AT_sOU_v_apex_inf_nl=DS_AT_sOU_v_apex_inf_nl, 
        DS_AT_sOU_v_flower_nl=DS_AT_sOU_v_flower_nl, DS_AT_sOU_v_stamen_nl=DS_AT_sOU_v_stamen_nl, 
        DS_AT_sOU_v_carpel_nl=DS_AT_sOU_v_carpel_nl)


      # Get fit for data from Brawand11_sou_v_div_rates
      Br11_sOU_v_brain_nl <- as.data.frame(do.call(rbind, lapply(x_Br_grid, nl_model, 
        a = 0.084866, b = 0.745355, c = -0.004761))) # compSouVDivRates11[49:54, ]
      Br11_sOU_v_cereb_nl <- as.data.frame(do.call(rbind, lapply(x_Br_grid, nl_model, 
        a = 0.103590, b = 0.710710, c = -0.007094))) # compSouVDivRates11[55:60, ]
      Br11_sOU_v_heart_nl <- as.data.frame(do.call(rbind, lapply(x_Br_grid, nl_model, 
        a = 0.159070, b = 0.759801, c = -0.005688))) # compSouVDivRates11[61:66, ]
      Br11_sOU_v_kidney_nl <- as.data.frame(do.call(rbind, lapply(x_Br_grid, nl_model, 
        a = 0.07539, b = 0.60388, c = -0.01536))) # compSouVDivRates11[67:72, ]
      Br11_sOU_v_liver_nl <- as.data.frame(do.call(rbind, lapply(x_Br_grid, nl_model, 
        a = 0.11775, b = 0.54527, c = -0.01522))) # compSouVDivRates11[73:78, ]
      Br11_sOU_v_testis_nl <- as.data.frame(do.call(rbind, lapply(x_Br_grid, nl_model, 
        a = 0.268853, b = 2.245665, c = -0.002769))) # compSouVDivRates11[79:83, ]

      Br11_sOU_v_nl_list <- list(Br11_sOU_v_brain_nl=Br11_sOU_v_brain_nl,
        Br11_sOU_v_cereb_nl=Br11_sOU_v_cereb_nl, Br11_sOU_v_heart_nl=Br11_sOU_v_heart_nl, 
        Br11_sOU_v_kidney_nl=Br11_sOU_v_kidney_nl, Br11_sOU_v_liver_nl=Br11_sOU_v_liver_nl, 
        Br11_sOU_v_testis_nl=Br11_sOU_v_testis_nl)


      # Get fit for data from Brawand11_div_rates
      Br11_pea_dist_brain_nl <- as.data.frame(do.call(rbind, lapply(x_Br_grid, nl_model, 
        a = 0.206154, b = 0.423120, c = -0.004325))) # compDivRates11[49:54, ]
      Br11_pea_dist_cereb_nl <- as.data.frame(do.call(rbind, lapply(x_Br_grid, nl_model, 
        a = 0.19167, b = 0.35432, c = -0.01104))) # compDivRates11[55:60, ]
      Br11_pea_dist_heart_nl <- as.data.frame(do.call(rbind, lapply(x_Br_grid, nl_model, 
        a = 0.23327, b = 0.35153, c = -0.01156))) # compDivRates11[61:66, ]
      Br11_pea_dist_kidney_nl <- as.data.frame(do.call(rbind, lapply(x_Br_grid, nl_model, 
        a = 0.18597, b = 0.34359, c = -0.02366))) # compDivRates11[67:72, ]
      Br11_pea_dist_liver_nl <- as.data.frame(do.call(rbind, lapply(x_Br_grid, nl_model, 
        a = 0.20957, b = 0.33178, c = -0.02351))) # compDivRates11[73:78, ]
      Br11_pea_dist_testis_nl <- as.data.frame(do.call(rbind, lapply(x_Br_grid, nl_model, 
        a = 0.263177, b = 0.466058, c = -0.008119))) # compDivRates11[79:83, ]

      Br11_pea_dist_nl_list <- list(Br11_pea_dist_brain_nl=Br11_pea_dist_brain_nl,
        Br11_pea_dist_cereb_nl=Br11_pea_dist_cereb_nl, Br11_pea_dist_heart_nl=Br11_pea_dist_heart_nl, 
        Br11_pea_dist_kidney_nl=Br11_pea_dist_kidney_nl, Br11_pea_dist_liver_nl=Br11_pea_dist_liver_nl, 
        Br11_pea_dist_testis_nl=Br11_pea_dist_testis_nl)


      # Get fit for data from compDivRates (re-analyzed Brawand data)
      Br_pea_dist_brain_nl <- as.data.frame(do.call(rbind, lapply(x_Br_grid, nl_model, 
        a = 0.27502, b = 0.36373, c = -0.02263))) # compDivRates[49:54, ]
      Br_pea_dist_cereb_nl <- as.data.frame(do.call(rbind, lapply(x_Br_grid, nl_model, 
        a = 0.27552, b = 0.38217, c = -0.01613))) # compDivRates[55:60, ]
      Br_pea_dist_heart_nl <- as.data.frame(do.call(rbind, lapply(x_Br_grid, nl_model, 
        a = 0.28637, b = 0.37893, c = -0.01287))) # compDivRates[61:66, ]
      Br_pea_dist_kidney_nl <- as.data.frame(do.call(rbind, lapply(x_Br_grid, nl_model, 
        a = 0.27297, b = 0.38036, c = -0.02766))) # compDivRates[67:72, ]
      Br_pea_dist_liver_nl <- as.data.frame(do.call(rbind, lapply(x_Br_grid, nl_model, 
        a = 0.27202, b = 0.36434, c = -0.02531))) # compDivRates[73:78, ]
      Br_pea_dist_testis_nl <- as.data.frame(do.call(rbind, lapply(x_Br_grid, nl_model, 
        a = 0.31522, b = 0.46063, c = -0.01118))) # compDivRates[79:83, ]

      Br_pea_dist_nl_list <- list(Br_pea_dist_brain_nl=Br_pea_dist_brain_nl,
        Br_pea_dist_cereb_nl=Br_pea_dist_cereb_nl, Br_pea_dist_heart_nl=Br_pea_dist_heart_nl, 
        Br_pea_dist_kidney_nl=Br_pea_dist_kidney_nl, Br_pea_dist_liver_nl=Br_pea_dist_liver_nl, 
        Br_pea_dist_testis_nl=Br_pea_dist_testis_nl)


      # Get fit for data from compSouVDivRates (re-analyzed Brawand data)
      Br_sOU_v_brain_nl <- as.data.frame(do.call(rbind, lapply(x_Br_grid, nl_model, 
        a = 0.23473, b = 0.65771, c = -0.02006))) # compSouVDivRates[49:54, ]
      Br_sOU_v_cereb_nl <- as.data.frame(do.call(rbind, lapply(x_Br_grid, nl_model, 
        a = 0.323369, b = 0.922935, c = -0.008036))) # compSouVDivRates[55:60, ]
      Br_sOU_v_heart_nl <- as.data.frame(do.call(rbind, lapply(x_Br_grid, nl_model, 
        a = 0.33516, b = 0.85098, c = -0.00747))) # compSouVDivRates[61:66, ]
      Br_sOU_v_kidney_nl <- as.data.frame(do.call(rbind, lapply(x_Br_grid, nl_model, 
        a = 0.30916, b = 0.81354, c = -0.01609))) # compSouVDivRates[67:72, ]
      Br_sOU_v_liver_nl <- as.data.frame(do.call(rbind, lapply(x_Br_grid, nl_model, 
        a = 0.29498, b = 0.69386, c = -0.01699))) # compSouVDivRates[73:78, ]
      Br_sOU_v_testis_nl <- as.data.frame(do.call(rbind, lapply(x_Br_grid, nl_model, 
        a = 0.435973, b = 2.422398, c = -0.003065))) # compSouVDivRates[79:83, ]

      Br_sOU_v_nl_list <- list(Br_sOU_v_brain_nl=Br_sOU_v_brain_nl,
        Br_sOU_v_cereb_nl=Br_sOU_v_cereb_nl, Br_sOU_v_heart_nl=Br_sOU_v_heart_nl, 
        Br_sOU_v_kidney_nl=Br_sOU_v_kidney_nl, Br_sOU_v_liver_nl=Br_sOU_v_liver_nl, 
        Br_sOU_v_testis_nl=Br_sOU_v_testis_nl)



      # Create df with regression coordinates, comp_organ, div_times, correlation, and dataset
      formatNLM.table <- function(x) {

        fname <- deparse(substitute(x))

        dataset_id <- gsub( "_.*$", "", fname)

        fname_r <- substr(fname,1,nchar(fname)-3)
        organ <- sub('.*\\_', '', fname_r)

        if (dataset_id == "DS") {

          dataset <- "Angiosperms "
          div_times <- data.frame(div_times = x_DS_grid)  ## times from prediction grid

        } else {

          dataset <- "Mammals"
          div_times <- data.frame(div_times = x_Br_grid)  ## times from prediction grid
        }

        dataset <- rep(dataset, 200)
        dataset <- data.frame(dataset = dataset)
        organ <- rep(organ, 200)
        comp_organ <- data.frame(comp_organ = organ)

        nlm_table <- cbind(comp_organ, div_times, x, dataset)
        colnames(nlm_table)[3] <- "correlation"

        return(nlm_table)
      }


      # Create tables of nlm slope values for all organs with comp_organ, div_times and dataset
      # columns for plotting non-linear regressions for individual organs
      DS_AL_pea_dist_root_nl_df <- as.data.frame(formatNLM.table(DS_AL_pea_dist_root_nl))
      DS_AL_pea_dist_hypo_nl_df <- as.data.frame(formatNLM.table(DS_AL_pea_dist_hypo_nl))
      DS_AL_pea_dist_leaf_nl_df <- as.data.frame(formatNLM.table(DS_AL_pea_dist_leaf_nl))
      DS_AL_pea_dist_apex_veg_nl_df <- as.data.frame(formatNLM.table(DS_AL_pea_dist_apex_veg_nl))
      DS_AL_pea_dist_apex_inf_nl_df <- as.data.frame(formatNLM.table(DS_AL_pea_dist_apex_inf_nl))
      DS_AL_pea_dist_flower_nl_df <- as.data.frame(formatNLM.table(DS_AL_pea_dist_flower_nl))
      DS_AL_pea_dist_stamen_nl_df <- as.data.frame(formatNLM.table(DS_AL_pea_dist_stamen_nl))
      DS_AL_pea_dist_carpel_nl_df <- as.data.frame(formatNLM.table(DS_AL_pea_dist_carpel_nl))

      DS_AL_pea_dist_nlm_coord <- rbind(DS_AL_pea_dist_root_nl_df, DS_AL_pea_dist_hypo_nl_df, 
        DS_AL_pea_dist_leaf_nl_df, DS_AL_pea_dist_apex_veg_nl_df, DS_AL_pea_dist_apex_inf_nl_df, 
        DS_AL_pea_dist_flower_nl_df, DS_AL_pea_dist_stamen_nl_df, DS_AL_pea_dist_carpel_nl_df)

      DS_AL_pea_dist_nlm_sp <- DS_AL_pea_dist_nlm_coord
      DS_AL_pea_dist_nlm_sp$dataset <- rep(c("Angiosperms.AL "), nrow(DS_AL_pea_dist_nlm_sp))


      DS_AT_pea_dist_root_nl_df <- as.data.frame(formatNLM.table(DS_AT_pea_dist_root_nl))
      DS_AT_pea_dist_hypo_nl_df <- as.data.frame(formatNLM.table(DS_AT_pea_dist_hypo_nl))
      DS_AT_pea_dist_leaf_nl_df <- as.data.frame(formatNLM.table(DS_AT_pea_dist_leaf_nl))
      DS_AT_pea_dist_apex_veg_nl_df <- as.data.frame(formatNLM.table(DS_AT_pea_dist_apex_veg_nl))
      DS_AT_pea_dist_apex_inf_nl_df <- as.data.frame(formatNLM.table(DS_AT_pea_dist_apex_inf_nl))
      DS_AT_pea_dist_flower_nl_df <- as.data.frame(formatNLM.table(DS_AT_pea_dist_flower_nl))
      DS_AT_pea_dist_stamen_nl_df <- as.data.frame(formatNLM.table(DS_AT_pea_dist_stamen_nl))
      DS_AT_pea_dist_carpel_nl_df <- as.data.frame(formatNLM.table(DS_AT_pea_dist_carpel_nl))

      DS_AT_pea_dist_nlm_coord <- rbind(DS_AT_pea_dist_root_nl_df, DS_AT_pea_dist_hypo_nl_df, 
        DS_AT_pea_dist_leaf_nl_df, DS_AT_pea_dist_apex_veg_nl_df, DS_AT_pea_dist_apex_inf_nl_df, 
        DS_AT_pea_dist_flower_nl_df, DS_AT_pea_dist_stamen_nl_df, DS_AT_pea_dist_carpel_nl_df)

      DS_AT_pea_dist_nlm_sp <- DS_AT_pea_dist_nlm_coord
      DS_AT_pea_dist_nlm_sp$dataset <- rep(c("Angiosperms.AT "), nrow(DS_AT_pea_dist_nlm_sp))


      # Get final table for organ regression plot
      nlmPea_coor11_AT_AL <- rbind(DS_AT_pea_dist_nlm_sp, DS_AL_pea_dist_nlm_sp)


      # Make pea dist table for DevSeq-AT and DevSeq.AL organs
      compDivRates11_AT_sp <- compDivRates11[1:48,]
      compDivRates11_AL_sp <- compDivRates11_AL[1:48,]
      pea_dist_organs <- c("root", "hypo", "leaf", "veg", "inf", "flower", "stamen", "carpel")
      pea_dist_organs <- rep(pea_dist_organs, each=6)
      pea_dist_organs <- as.data.frame(pea_dist_organs)
      compDivRates11_AT_sp$comp_organ <- pea_dist_organs
      compDivRates11_AL_sp$comp_organ <- pea_dist_organs
      compDivRates11_AT_sp$dataset <- rep(c("Angiosperms.AT "), nrow(compDivRates11_AT_sp))
      compDivRates11_AT_sp$comp_organ <- factor(unlist(compDivRates11_AT_sp$comp_organ))
      compDivRates11_AL_sp$dataset <- rep(c("Angiosperms.AL "), nrow(compDivRates11_AL_sp))
      compDivRates11_AL_sp$comp_organ <- factor(unlist(compDivRates11_AL_sp$comp_organ))
      pea_dist_organs_AT_AL <- rbind(compDivRates11_AT_sp, compDivRates11_AL_sp)
      pea_dist_organs_AT_AL$comp_spec <- factor(pea_dist_organs_AT_AL$comp_spec)
      pea_dist_organs_AT_AL$dataset <- factor(pea_dist_organs_AT_AL$dataset)


      # Get mean value of slopes
      DS_AL_sOU_v_nl_mean <- rowMeans(do.call(cbind, DS_AL_sOU_v_nl_list))
      DS_AT_sOU_v_nl_mean <- rowMeans(do.call(cbind, DS_AT_sOU_v_nl_list))
      Br11_sOU_v_nl_mean <- rowMeans(do.call(cbind, Br11_sOU_v_nl_list))
      Br_sOU_v_nl_mean <- rowMeans(do.call(cbind, Br_sOU_v_nl_list))


      # Get cumulative slope values
      getCum.NL.Slope <- function(x){

        fname <- deparse(substitute(x))

        dataset_id <- gsub( "_.*$", "", fname)

        if (dataset_id == "DS") {

          x_grid <- x_DS_grid

        } else x_grid <- x_Br_grid

        x <- as.numeric(x)
        slopes = diff(x)/diff(x_grid)
        slopes_cumsum <- cumsum(as.data.frame(slopes)[,1])
        slopes_cumsum <- as.numeric(as.data.frame(slopes_cumsum)[,1])

        return(slopes_cumsum)
      }

      DS_AL_sOU_v_nl_mean_cumslopes <- getCum.NL.Slope(DS_AL_sOU_v_nl_mean)
      DS_AT_sOU_v_nl_mean_cumslopes <- getCum.NL.Slope(DS_AT_sOU_v_nl_mean)
      Br11_sOU_v_nl_mean_cumslopes <- getCum.NL.Slope(Br11_sOU_v_nl_mean)
      Br_sOU_v_nl_mean_cumslopes <- getCum.NL.Slope(Br_sOU_v_nl_mean)


      # Compute slope
      getDS.NL.Slope <- function(x){

        fname <- deparse(substitute(x))

        dataset_id <- gsub( "_.*$", "", fname)

        if (dataset_id == "DS") {

          x_grid <- x_DS_grid

        } else x_grid <- x_Br_grid

        x <- as.numeric(x[,1])
        slopes = diff(x)/diff(x_grid)
        slopes_avg <- mean(slopes)
        slopes_avg <- as.numeric(as.data.frame(slopes_avg))
        
        return(slopes_avg)

      }

      DevSeq_AL_pea_dist_nl_slopes <- as.data.frame(do.call(rbind, lapply(DS_AL_pea_dist_nl_list, getDS.NL.Slope)))
      DevSeq_AL_sOU_v_nl_slopes <- as.data.frame(do.call(rbind, lapply(DS_AL_sOU_v_nl_list, getDS.NL.Slope)))
      DevSeq_AT_pea_dist_nl_slopes <- as.data.frame(do.call(rbind, lapply(DS_AT_pea_dist_nl_list, getDS.NL.Slope)))
      DevSeq_AT_sOU_v_nl_slopes <- as.data.frame(do.call(rbind, lapply(DS_AT_sOU_v_nl_list, getDS.NL.Slope)))
      Br11_sOU_v_nl_slopes <- as.data.frame(do.call(rbind, lapply(Br11_sOU_v_nl_list, getDS.NL.Slope)))
      Br11_pea_dist_nl_slopes <- as.data.frame(do.call(rbind, lapply(Br11_pea_dist_nl_list, getDS.NL.Slope)))
      Br_pea_dist_nl_slopes <- as.data.frame(do.call(rbind, lapply(Br_pea_dist_nl_list, getDS.NL.Slope)))
      Br_sOU_v_nl_slopes <- as.data.frame(do.call(rbind, lapply(Br_sOU_v_nl_list, getDS.NL.Slope)))


      # Add organ id column
      formatNL.Slope <- function(x, data_set, regr_mod) {

        names(x) <- regr_mod

        if ((regr_mod == "DS_AL_pea_nlm") || (regr_mod == "DS_AL_sOU_nlm") 
          || (regr_mod == "DS_AT_pea_nlm") || (regr_mod == "DS_AT_sOU_nlm")) {

          organ <- data.frame(sample = c("Root", "Hypocotyl", "Leaf", "Apex_veg", "Apex_inf", 
            "Flower", "Stamen", "Carpel"))

        } else {

          organ <- data.frame(sample = c("Brain", "Cerebellum", "Heart", "Kidney", "Liver", 
            "Testis"))
        }

        slopes_df <- cbind(organ, x)
        rownames(slopes_df) <- NULL

        return(slopes_df)

      }


      DevSeq_AL_pea_nlm_slopes <- formatNL.Slope(DevSeq_AL_pea_dist_nl_slopes, regr_mod = "DS_AL_pea_nlm")
      DevSeq_AL_sOU_nlm_slopes <- formatNL.Slope(DevSeq_AL_sOU_v_nl_slopes, regr_mod = "DS_AL_sOU_nlm")
      DevSeq_AL_nlm_slopes <- cbind(DevSeq_AL_sOU_nlm_slopes, DevSeq_AL_sOU_nlm_slopes[-1], 
        DevSeq_AL_pea_nlm_slopes[-1], DevSeq_AL_pea_nlm_slopes[-1])
      colnames(DevSeq_AL_nlm_slopes) <- c("sample", "DS_AL_Br11_sOU_nlm", "DS_AL_Br_sOU_nlm", 
        "DS_AL_Br11_pea_nlm", "DS_AL_Br_pea_nlm")

      DevSeq_AT_pea_nlm_slopes <- formatNL.Slope(DevSeq_AT_pea_dist_nl_slopes, regr_mod = "DS_AT_pea_nlm")
      DevSeq_AT_sOU_nlm_slopes <- formatNL.Slope(DevSeq_AT_sOU_v_nl_slopes, regr_mod = "DS_AT_sOU_nlm")
      DevSeq_AT_nlm_slopes <- cbind(DevSeq_AT_sOU_nlm_slopes, DevSeq_AT_sOU_nlm_slopes[-1], 
        DevSeq_AT_pea_nlm_slopes[-1], DevSeq_AT_pea_nlm_slopes[-1])
      colnames(DevSeq_AT_nlm_slopes) <- c("sample", "DS_AT_Br11_sOU_nlm", "DS_AT_Br_sOU_nlm", 
        "DS_AT_Br11_pea_nlm", "DS_AT_Br_pea_nlm")

      Br11_pea_nlm_slopes <- formatNL.Slope(Br11_pea_dist_nl_slopes, regr_mod = "Br11_pea_nlm")
      Br11_sOU_nlm_slopes <- formatNL.Slope(Br11_sOU_v_nl_slopes, regr_mod = "Br11_sOU_nlm")

      Br_pea_nlm_slopes <- formatNL.Slope(Br_pea_dist_nl_slopes, regr_mod = "Br_pea_nlm")
      Br_sOU_nlm_slopes <- formatNL.Slope(Br_sOU_v_nl_slopes, regr_mod = "Br_sOU_nlm")
      Br_nlm_slopes <- cbind(Br11_sOU_nlm_slopes, Br_sOU_nlm_slopes[-1], Br11_pea_nlm_slopes[-1], 
        Br_pea_nlm_slopes[-1])

      # Generate final data table for DevSeq and Brawand nlm slopes
      Br_nlm_slopes_AL <- Br_nlm_slopes
      colnames(Br_nlm_slopes_AL) <- c("sample", "DS_AL_Br11_sOU_nlm", "DS_AL_Br_sOU_nlm", 
        "DS_AL_Br11_pea_nlm", "DS_AL_Br_pea_nlm")
      Br_nlm_slopes_AT <- Br_nlm_slopes
      colnames(Br_nlm_slopes_AT) <- c("sample", "DS_AT_Br11_sOU_nlm", "DS_AT_Br_sOU_nlm", 
        "DS_AT_Br11_pea_nlm", "DS_AT_Br_pea_nlm")

      DS_AL_Br_nlm_slopes <- rbind(DevSeq_AL_nlm_slopes, Br_nlm_slopes_AL)
      DS_AT_Br_nlm_slopes <- rbind(DevSeq_AT_nlm_slopes, Br_nlm_slopes_AT)




#---- Apply LOESS regression to sOU and pearson dist expression data and compare slopes -----


      getLOESS.Coord <- function(organ_data) {

        comp_organ <- unique(organ_data$comp_organ)

        temp <- loess.smooth(organ_data$div_times, organ_data$correlation, span = 1.18, 
            degree = 2, family="symmetric", evaluation = 200)

        # Obtain coordinates of the smooth curve
        div_times <- as.numeric(unlist(temp$x))
        correlation <- as.numeric(unlist(temp$y))
        comp_organ <- rep(comp_organ, 200)

        experiment <- unique(organ_data$dataset)
        experiment <- rep(experiment, 200)

        loess_out <- data.frame(comp_organ = comp_organ, div_times = div_times, correlation = correlation, 
          dataset = experiment)

        return(loess_out)
      }


      # Set up lists containing sOU expression distances
      devseqSouV_organ_lst <- list(compSouVDivRates11_AL[1:6,], compSouVDivRates11_AL[7:12,], compSouVDivRates11_AL[13:18,], 
        compSouVDivRates11_AL[19:24,], compSouVDivRates11_AL[25:30,], compSouVDivRates11_AL[31:36,], compSouVDivRates11_AL[37:42,], 
        compSouVDivRates11_AL[43:48,])

      brawandSouV11_organ_lst <- list(compSouVDivRates11_AL[49:54,], compSouVDivRates11_AL[55:60,], compSouVDivRates11_AL[61:66,], 
        compSouVDivRates11_AL[67:72,],compSouVDivRates11_AL[73:78,], compSouVDivRates11_AL[79:83,])

      devseq_AT_SouV_organ_lst <- list(compSouVDivRates11[1:6,], compSouVDivRates11[7:12,], compSouVDivRates11[13:18,], 
        compSouVDivRates11[19:24,], compSouVDivRates11[25:30,], compSouVDivRates11[31:36,], compSouVDivRates11[37:42,], 
        compSouVDivRates11[43:48,])

      brawandSouV_organ_lst <- list(compSouVDivRates[49:54,], compSouVDivRates[55:60,], compSouVDivRates[61:66,], 
        compSouVDivRates[67:72,],compSouVDivRates[73:78,], compSouVDivRates[79:83,])


      # Set up lists containing metric pearson expression distances
      devseq_organ_lst <- list(DevSeq_AL_div_rates[1:6,], DevSeq_AL_div_rates[7:12,], DevSeq_AL_div_rates[13:18,], 
        DevSeq_AL_div_rates[19:24,], DevSeq_AL_div_rates[25:30,], DevSeq_AL_div_rates[31:36,], DevSeq_AL_div_rates[37:42,], 
        DevSeq_AL_div_rates[43:48,])


      # Get LOESS coordinates for DevSeq and Brawand data
      DevSeqSouV_loess_coord <- as.data.frame(do.call(rbind, lapply(devseqSouV_organ_lst, getLOESS.Coord)))
      Brawand11SouV_loess_coord <- as.data.frame(do.call(rbind, lapply(brawandSouV11_organ_lst, getLOESS.Coord)))
      loessSouV_coor11_AL <- rbind(DevSeqSouV_loess_coord, Brawand11SouV_loess_coord)


      # Format loess df for ggplot2
      formatLOESS.DF <- function(df) {

        df$div_times <- as.numeric(df$div_times)
        df$correlation <- as.numeric(df$correlation)
        df$dataset <- as.factor(df$dataset)

        return(df)
      }

      loessSouV_coor11_AL <- formatLOESS.DF(loessSouV_coor11_AL)



      # Get cumulative slope values
      getCum.LOESS.Slope <- function(x, span = span, degree = degree, family = family){

        temp <- loess.smooth(x$div_times, x$correlation, span = span, 
          degree = degree, family = family, evaluation = 160)

        # Get slope values
        slopes = diff(temp$y)/diff(temp$x)

        slopes_cumsum <- cumsum(as.data.frame(slopes)[,1])
        slopes_cumsum <- as.numeric(as.data.frame(slopes_cumsum)[,1])

        return(slopes_cumsum)
      }


      # Get LOESS slope mean and CI
      getLoessStats <- function(loess_data) {

        nsamples <- ncol(loess_data)

        loess_mean <- rowMeans(loess_data)

        sd_out <- by(loess_data, 1:nrow(loess_data), function(row) sd <- sd(row))
        sd_out <- sd_out[1:length(sd_out)]
        sd_out <- as.data.frame(sd_out)
        sd_out <- as.numeric(sd_out[,1])

        error <- qnorm(0.975) * sd_out/sqrt(nsamples)

        mean_li <- loess_mean - error
        mean_ri <- loess_mean + error

        loess_stat <- data.frame(correlation = loess_mean, li = mean_li, ri = mean_ri)

        return(loess_stat)
      }


      # Get mean value of slopes
      DevSeqSouV_AL_loess_mean <- do.call(cbind, lapply(devseqSouV_organ_lst, getCum.LOESS.Slope, 
        span = 1, degree = 2, family = "gaussian"))
      DevSeqSouV_AL_loess_mean <- getLoessStats(DevSeqSouV_AL_loess_mean)

      DevSeqSouV_AT_loess_mean <- do.call(cbind, lapply(devseq_AT_SouV_organ_lst, getCum.LOESS.Slope, 
        span = 1, degree = 2, family = "gaussian"))
      DevSeqSouV_AT_loess_mean <- getLoessStats(DevSeqSouV_AT_loess_mean)

      brawandSouV11_loess_mean <- do.call(cbind, lapply(brawandSouV11_organ_lst, getCum.LOESS.Slope, 
        span = 1, degree = 2, family = "gaussian"))
      brawandSouV11_loess_mean <- getLoessStats(brawandSouV11_loess_mean)

      brawandSouV_loess_mean <- do.call(cbind, lapply(brawandSouV_organ_lst, getCum.LOESS.Slope, 
        span = 1, degree = 2, family = "gaussian"))
      brawandSouV_loess_mean <- getLoessStats(brawandSouV_loess_mean)


      # Prepare data for ggplot2
      x_DS_grid_cum <- seq(7.1, 160, length = 159)  ## prediction grid
      x_Br_grid_cum <- seq(6.7, 159, length = 159)  ## prediction grid

      div_times <- c(x_DS_grid_cum, x_DS_grid_cum, x_Br_grid_cum, x_Br_grid_cum)
      div_times <- as.data.frame(div_times)
      colnames(div_times) <- "div_times"
      dataset <- data.frame(rep(c("Angiosperms.AT", "Angiosperms.AL", "Mammals.11", "Mammals.re-an."), each=159))
      colnames(dataset) <- "dataset"
      correlation <- rbind(DevSeqSouV_AT_loess_mean, DevSeqSouV_AL_loess_mean, brawandSouV11_loess_mean, 
        brawandSouV_loess_mean)
      sOU_v_loess_cum_slopes <- data.frame(div_times, correlation, dataset)

      sOU_v_loess_cum_slopes$div_times <- as.numeric(sOU_v_loess_cum_slopes$div_times)
      sOU_v_loess_cum_slopes$correlation <- as.numeric(sOU_v_loess_cum_slopes$correlation)
      sOU_v_loess_cum_slopes$li <- as.numeric(sOU_v_loess_cum_slopes$li)
      sOU_v_loess_cum_slopes$ri <- as.numeric(sOU_v_loess_cum_slopes$ri)
      sOU_v_loess_cum_slopes$dataset <- factor(sOU_v_loess_cum_slopes$dataset)



      getLOESS.Slopes <- function(organ_data) {

        comp_organ <- unique(organ_data$comp_organ)

        temp <- loess.smooth(organ_data$div_times, organ_data$correlation, span = 1.18, 
          degree = 2, family="symmetric", evaluation = 200)

        # Get slope values
        slopes = diff(temp$y)/diff(temp$x)
        slopes_avg <- mean(slopes)
        slopes_avg <- as.numeric(as.data.frame(slopes_avg))

        return(slopes_avg)
      }


      # Get LOESS slopes for DevSeq and Brawand data
      # For sOU expression distances
      DevSeqSouV_AL_loess_slopes <- as.data.frame(do.call(rbind, lapply(devseqSouV_organ_lst, getLOESS.Slopes)))
      Brawand11SouV_loess_slopes <- as.data.frame(do.call(rbind, lapply(brawandSouV11_organ_lst, getLOESS.Slopes)))
      sOU_loess_DevSeq_AL_Br11_wilcox <- wilcox.test(as.numeric(unlist(DevSeqSouV_AL_loess_slopes)), as.numeric(unlist(Brawand11SouV_loess_slopes)))$p.value

      # For metric pearson expression distances
      DevSeq_AL_loess_slopes <- as.data.frame(do.call(rbind, lapply(devseq_organ_lst, getLOESS.Slopes)))


      # Write slope values to csv file
      DevSeq_AL_loess_slopes <- cbind(DevSeqSouV_AL_loess_slopes, DevSeq_AL_loess_slopes)
      colnames(DevSeq_AL_loess_slopes) <- c("DS_AL_sOU_loess", "DS_AL_pea_loess")
      ds_organs <- data.frame(sample = c("Root", "Hypocotyl", "Leaf", "Apex_veg", "Apex_inf", "Flower", "Stamen", "Carpel"))
      DevSeq_AL_loess_slopes <- cbind(ds_organs, DevSeq_AL_loess_slopes)



      # Show message
      message("Writing data tables...")

      slopes_regr_list <- list(DevSeq_AL_loess_slopes = DevSeq_AL_loess_slopes, 
        DevSeq_AL_log_slopes = DevSeq_AL_log_slopes, DS_AL_Br_nlm_slopes = DS_AL_Br_nlm_slopes, 
        DS_AT_Br_nlm_slopes = DS_AT_Br_nlm_slopes)

      for(i in names(slopes_regr_list)) { 
        write.table(slopes_regr_list[[i]], file = file.path(out_dir, "output", "data", paste0(i,".txt")), 
          sep="\t", col.names=TRUE, row.names = FALSE, dec=".", quote = FALSE)
      }


      # Create p-value containing test strings for plots
      sOU_loess_DevSeq_AL_Br11_slope_p <- paste("P =", formatC(sOU_loess_DevSeq_AL_Br11_wilcox, format="e", digits=0))

      # Change dataset ID for shorter legend
      compSouVDivRates11_AL$dataset <- c(rep('Angiosperms.AL ', 48), rep('Mammals.11', 35))
      loessSouV_coor11_AL$dataset <- c(rep('Angiosperms.AL ', 1600), rep('Mammals.11', 1200))



   # Make sOU GE divergence plot showing individual organ regressions for SI
   makeOrgRegPlot <- function(data1, data2, coefficient, expr_estimation, p_value, pos) {

      fname <- sprintf('%s.jpg', paste("compSouVDivRates11_AL_loess", expr_estimation, pos, sep="_"))
      y_min <- 0.055
      y_max <- 1.7
      col_breaks <- c("Angiosperms.AL ", "Mammals.11")
      fill_breaks <- c("Angiosperms.AL ", "Mammals.11")
      y_breaks <- c(0.2,0.4,0.6,0.8,1,1.2,1.4,1.6)
      y_title <- "Expression distance"

      if (pos == "main") {

            plot_wdt <- 12.535
            plot_hdt <- 8
            legend_x_pos <- 0.2825
            legend_y_pos <- 0.914
            linewd <- 3
            point_size <- 5.75
            txt_x_pos <- 16.25
            txt_y_pos <- 1.283
            pan_boarder <- 1.75
            axis_txt_size <- 21.25
            axis_ticks_s <- 0.7
            title_face <- "plain"

      } else if (pos == "ext") {

            plot_wdt <- 9.5 # condenced plot width for suppl
            plot_hdt <- 6.75 # condenced plot width for suppl
            legend_x_pos <- 0.368
            legend_y_pos <- 0.926
            linewd <- 2.5
            point_size <- 4.75
            txt_x_pos <- 16.25
            txt_y_pos <- 1.283
            pan_boarder <- 1.8
            axis_txt_size <- 21.75
            axis_ticks_s <- 0.95
            title_face <- "bold"
      }

      ds_col <- rep(c("#728acb"), 48)
      bw_col <- rep(c("red"), 35)
      col_scale <- c("#728acb", "red")
      fill_scale <- c("#728acb", "red")
      shape_scale <- c(16, 17)

      fill_col <- c(as.character(ds_col), as.character(bw_col))

      p <- ggplot() 
      p <- p + geom_line(size = linewd, data = data1, aes(x = div_times, y = correlation, group = comp_organ, 
        colour = dataset, linetype = dataset))
      p <- p + geom_point(size = point_size, data = data2, aes(x = div_times, y = correlation, group = comp_organ, 
        colour = dataset, shape = comp_organ, stroke = 2.7))
      p <- p + scale_x_continuous(limits = c(0,162), expand = c(0.02,0), breaks = c(0,20,40,60,80,100,120,140,160)) + 
      scale_y_continuous(limits = c(y_min, y_max), expand = c(0.02, 0), breaks = y_breaks) + 
      scale_color_manual(values = col_scale, breaks = col_breaks) + 
      scale_fill_manual(values = fill_scale, breaks = fill_breaks) + 
      scale_shape_manual(values = c(0, 8, 2, 5, 3, 4, 6, 1, 15, 16, 17, 18, 19, 10)) + 
      annotate("text", x = txt_x_pos, y = txt_y_pos, label = p_value, size = 8) + 
      guides(color = guide_legend(ncol=2, keywidth = 0.4, keyheight = 0.4, default.unit = "inch"), 
        shape = FALSE)

      q <- p + theme_bw() + xlab("Divergence time (Myr)") + ylab(y_title) + 
      theme(text=element_text(size=16), 
        axis.ticks.length=unit(0.35, "cm"), 
        axis.ticks = element_line(colour = "black", size = axis_ticks_s),  
        plot.margin = unit(c(0.55, 1.175, 0.5, 0.4),"cm"), 
        axis.title.y = element_text(size=24.5, margin = margin(t = 0, r = 15, b = 0, l = 11), colour="black", 
            face = title_face), 
        axis.title.x = element_text(size=24.5, margin = margin(t = 15.75, r = 0, b = 1, l = 0), colour="black", 
            face = title_face), 
        axis.text.x = element_text(size=axis_txt_size, angle=0, margin = margin(t = 5.5), colour="black"), 
        axis.text.y = element_text(size=axis_txt_size, angle=0, margin = margin(r = 5.5), colour="black"), 
        legend.box.background = element_rect(colour = NA, fill= "white" , size=1.0), 
        panel.border = element_rect(colour = "black", fill=NA, size=pan_boarder), 
        panel.grid.major = element_blank(),
        panel.grid.minor.x = element_blank(), 
        panel.grid.minor.y = element_blank(), 
        legend.position = c(legend_x_pos, legend_y_pos), 
        legend.title = element_blank(), 
        legend.text = element_text(size=22), 
        legend.spacing.x = unit(0.5, 'cm'), 
        legend.key.size = unit(0.95, "cm"), 
        legend.background=element_blank()) 

        ggsave(file = file.path(out_dir, "output", "plots", fname), plot = q, 
            width = plot_wdt, height = plot_hdt, dpi = 300, units = c("in"), limitsize = FALSE) 
  }

  makeOrgRegPlot(data1 = loessSouV_coor11_AL, data2 = compSouVDivRates11_AL, coefficient = coefficient, 
    expr_estimation = expr_estimation, pos = "ext", p_value = "")

  makeOrgRegPlot(data1 = loessSouV_coor11_AL, data2 = compSouVDivRates11_AL, coefficient = coefficient, 
    expr_estimation = expr_estimation, pos = "main", p_value = sOU_loess_DevSeq_AL_Br11_slope_p)



   # Make sOU GE divergence plot showing cumulative mean slope values for SI
   plotCumSlopes <- function(data, coefficient, expr_estimation) {


      fname <- sprintf('%s.jpg', paste("compSouV_loess_cum_slopes", coefficient, expr_estimation, sep="_"))
        
      col_breaks <- c("Angiosperms.AT", "Angiosperms.AL", "Mammals.11", "Mammals.re-an.")
      y_breaks <- c(0,0.2,0.4,0.6,0.8,1,1.2,1.4)

      col_scale <- c('#728acb', '#81bcef', 'red', 'red3')
      fill_scale <- c('#728acb', '#81bcef', 'red', 'red3')


      p <- ggplot(data=data, aes(x = div_times, y = correlation, group = dataset)) + 
      geom_ribbon(aes(ymin = data$li, ymax = data$ri, fill = dataset), alpha = 0.088, 
            linetype = 0, show.legend = FALSE) + 
      geom_line(size = 2.6, data = data, aes(x = div_times, y = correlation, group = dataset, 
        colour = dataset)) + 
      geom_line(data = data[1:159,], aes(x = div_times, y = correlation), color = '#81bcef', size = 2.5) +
      geom_line(data = data[160:318,], aes(x = div_times, y = correlation), color = '#728acb', size = 2.5) + 
      scale_x_continuous(limits = c(0,162), expand = c(0.02,0), breaks = c(0,20,40,60,80,100,120,140,160)) + 
      scale_y_continuous(limits = c(-0.01, 1.1675), expand = c(0.02, 0), breaks = y_breaks) + 
      scale_color_manual(values = col_scale, breaks = col_breaks) + 
      scale_fill_manual(values = fill_scale) + 
      scale_size(range = c(0.5, 12)) + 
      guides(color = guide_legend(ncol=1, keywidth = 0.4, keyheight = 0.4, default.unit = "inch"))

      q <- p + theme_bw() + xlab("Divergence time (Myr)") + ylab("Cumulative mean slope value") + 
      theme(text=element_text(size=16), 
        axis.ticks.length=unit(0.35, "cm"), 
        axis.ticks = element_line(colour = "black", size = 0.95),  
        plot.margin = unit(c(0.55, 1.175, 0.5, 0.4),"cm"), 
        axis.title.y = element_text(size=24.5, margin = margin(t = 0, r = 15, b = 0, l = 11), colour="black", 
            face = "bold"), 
        axis.title.x = element_text(size=24.5, margin = margin(t = 15.75, r = 0, b = 1, l = 0), colour="black", 
            face = "bold"), 
        axis.text.x = element_text(size=21.75, angle=0, margin = margin(t = 5.5), colour="black"), 
        axis.text.y = element_text(size=21.75, angle=0, margin = margin(r = 5.5), colour="black"), 
        legend.box.background = element_rect(colour = NA, fill= "white" , size=1.0), 
        panel.border = element_rect(colour = "black", fill=NA, size=1.8), 
        panel.grid.major = element_blank(),
        panel.grid.minor.x = element_blank(), 
        panel.grid.minor.y = element_blank(), 
        legend.position = c(0.22, 0.815), 
        legend.title = element_blank(), 
        legend.text = element_text(size=22), 
        legend.spacing.x = unit(0.5, 'cm'), 
        legend.key.size = unit(0.95, "cm"), 
        legend.background=element_blank()) 

        ggsave(file = file.path(out_dir, "output", "plots", fname), plot = q, 
            width = 9.5, height = 6.75, dpi = 300, units = c("in"), limitsize = FALSE) 
  }

  plotCumSlopes(data = sOU_v_loess_cum_slopes, coefficient = coefficient, expr_estimation = expr_estimation)



  # Plot AT and AL nlm's for pea dist
  # Change organ names for facet strip
  nlmPea_coor11_AT_AL_df <- nlmPea_coor11_AT_AL
  nlmPea_coor11_AT_AL_df$comp_organ <- nlmPea_coor11_AT_AL_df$comp_organ %<>% 
        gsub("root", "Root", .) %>% 
        gsub("hypo", "Hypocotyl", .) %>% 
        gsub("leaf", "Leaf", .) %>% 
        gsub("veg", "Apex veg", .) %>% 
        gsub("inf", "Apex inf", .) %>% 
        gsub("flower", "Flower", .) %>% 
        gsub("stamen", "Stamen", .) %>% 
        gsub("carpel", "Carpel", .)
  nlmPea_coor11_AT_AL_df$comp_organ <- factor(nlmPea_coor11_AT_AL_df$comp_organ, 
    levels=c("Root", "Hypocotyl", "Leaf", "Apex veg", "Apex inf", "Flower", "Stamen", "Carpel"))

  pea_dist_organs_AT_AL_df <- pea_dist_organs_AT_AL
  pea_dist_organs_AT_AL_df$comp_organ <- pea_dist_organs_AT_AL_df$comp_organ %<>% 
        gsub("root", "Root", .) %>% 
        gsub("hypo", "Hypocotyl", .) %>% 
        gsub("leaf", "Leaf", .) %>% 
        gsub("veg", "Apex veg", .) %>% 
        gsub("inf", "Apex inf", .) %>% 
        gsub("flower", "Flower", .) %>% 
        gsub("stamen", "Stamen", .) %>% 
        gsub("carpel", "Carpel", .)
  pea_dist_organs_AT_AL_df$comp_organ <- factor(pea_dist_organs_AT_AL_df$comp_organ, 
    levels=c("Root", "Hypocotyl", "Leaf", "Apex veg", "Apex inf", "Flower", "Stamen", "Carpel"))


  plotAT.AL.pea.NLM <- function(data, data2) {

    fname <- sprintf('%s.jpg', paste("AT_AL_pea_nlm_regression_slopes"))

    p <- ggplot(data=data, color = dataset, aes(x=div_times, y=correlation)) + 
    geom_line(size = 2.5, data = data, aes(x = div_times, y = correlation, group = dataset, 
        colour = dataset)) + 
    geom_point(data=data2, aes(shape = dataset, color = dataset, stroke = 7.5)) + 
    scale_color_manual(values=c('#728acb','#81bcef'), guide = "legend") + 
    scale_fill_manual(values=c('#728acb','#81bcef'), guide = "legend") + 
    scale_y_continuous(expand = c(0.15, 0), labels = comma) + 
    scale_x_continuous(expand = c(0.1, 0), breaks=c(0, 50, 100, 150)) + 
    guides(shape = guide_legend(override.aes = list(stroke = 7.5)))

    q <- p + theme_classic() + xlab("Divergence time (Myr)") + ylab("Pearson distance") + 
    theme(text=element_text(size = 16), 
      strip.text = element_text(size = 23.75), 
      strip.text.x = element_text(margin = margin(0.44, 0, 0.44, 0, "cm")), 
      strip.background = element_rect(colour = 'black', fill = NA, size = 2.2), 
      axis.ticks.length = unit(0.35, "cm"), 
      axis.ticks = element_line(colour = "black", size = 0.95), 
      axis.line = element_line(colour = 'black', size = 0.95), 
      plot.margin = unit(c(1, 1, 2.5, 0.4),"cm"), 
      axis.title.y = element_text(size=24.6, margin = margin(t = 0, r = 15.2, b = 0, l = 10.8), 
        colour="black", face = "bold"), 
      axis.title.x = element_text(size=24.6, margin = margin(t = 9.25, r = 0, b = 7.5, l = 0), 
        colour="black", face = "bold"), 
      axis.text.x = element_text(size=21.75, margin = margin(t = 5.5, b = 8), colour="black"), 
      axis.text.y = element_text(size=21.75, angle=0, margin = margin(l = 2.5, r = 7), colour="black"), 
      panel.spacing = unit(0.55, "cm"), 
      panel.grid.major = element_blank(),
      panel.grid.minor.x = element_blank(), 
      panel.grid.minor.y = element_blank(), 
      legend.position = c(0.125, -0.225), 
      legend.direction = "horizontal", 
      legend.title = element_blank(), 
      legend.text = element_text(size = 22.5), 
      legend.spacing.x = unit(0.5, 'cm'), 
      legend.key.size = unit(1.2, "cm"), 
      legend.background=element_blank()) 

    q <- q + facet_wrap(~ comp_organ, nrow = 1)

    ggsave(file = file.path(out_dir, "output", "plots", fname), plot = q, 
      width = 28.5, height = 6.5, dpi = 300, units = c("in"), limitsize = FALSE) 
  }

  plotAT.AL.pea.NLM(data = nlmPea_coor11_AT_AL_df, data2 = pea_dist_organs_AT_AL_df)


  }
   
}


makeCompAnalysisAL(expr_estimation="TPM", coefficient="pearson")





