# Prepare Brawand and DevSeq AL comparative ortholog gene expression data
# Thresholds: DevSeq 0.05 ERCC; Brawand 0.5 TPM (no ERCC spike-ins available)
# Data input: Brawand and DevSeq TPM expression tables of all samples



#-------------------------------------- Read data tables ---------------------------------------


makeCompAnylsisAL <- function(expr_estimation = c("TPM", "counts"), coefficient = c("pearson", "spearman"), 
    data_norm = c("intra-organ", "inter-organ")) {


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

    # Show error message if expression estimation or unknown expression estimation is chosen
    if ((missing(data_norm)) || (!is.element(data_norm, c("intra-organ", "inter-organ"))))
   
       stop(
       "Please choose one of the available data_norm data normalizations: 
       'intra-organ', 'inter-organ'",
       call. = TRUE
       )


    # Show startup message
    message("Reading data...")


    if (is.element("TPM", expr_estimation) && is.element("intra-organ", data_norm)) {
        genesExpr = file.path(in_dir, "Expression_data", "AL_core_intra_tpm_mat_deseq_sample_names.csv")

    } else if (is.element("TPM", expr_estimation) && is.element("inter-organ", data_norm)) {
        genesExpr = file.path(in_dir, "Expression_data", "AL_core_inter_tpm_mat_deseq_sample_names.csv")

    } else if (is.element("counts", expr_estimation) && is.element("intra-organ", data_norm)) {
        genesExpr = file.path(in_dir, "Expression_data", "AL_core_intra_count_mat_vsd_sample_names.csv")

    } else if (is.element("counts", expr_estimation) && is.element("inter-organ", data_norm)) {
        genesExpr = file.path(in_dir, "Expression_data", "AL_core_inter_count_mat_vsd_sample_names.csv")
    }


    # Get data normalization method
    if (is.element("intra-organ", data_norm)) {
        data_norm <- "intra-organ"

    } else if (is.element("inter-organ", data_norm)) {
        data_norm <- "inter-organ"
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


    # Stop function here to allow specific analysis of a single data set
    # return_list <- list("expr_estimation" = expr_estimation, "x" = x, "coefficient" = coefficient, "col_names" = col_names, "data_norm" = data_norm)
    # return(return_list)
    # }
    # return_objects <- makeCompAnylsisAL(expr_estimation="TPM", coefficient="pearson", data_norm="inter-organ") # read in DevSeq expression data
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

         # Reshape cor data frame to one column
         df_cor_rs <- data.frame(newcol = c(t(df_cor)), stringsAsFactors=FALSE)

         getError <- function(cor_data) {
            sd <- sd(cor_data)
            error <- qt(0.975, df = 9-1) * sd/sqrt(9)
         }

         sp1 <- mean(df_cor_rs[1:9,])
         sp1_li <- sp1 - getError(df_cor_rs[1:9,])
         sp1_ri <- sp1 + getError(df_cor_rs[1:9,])

         sp2 <- mean(df_cor_rs[10:18,])
         sp2_li <- sp2 - getError(df_cor_rs[10:18,])
         sp2_ri <- sp2 + getError(df_cor_rs[10:18,])

         sp3 <- mean(df_cor_rs[19:27,])
         sp3_li <- sp3 - getError(df_cor_rs[19:27,])
         sp3_ri <- sp3 + getError(df_cor_rs[19:27,])

         sp4 <- mean(df_cor_rs[28:36,])
         sp4_li <- sp4 - getError(df_cor_rs[28:36,])
         sp4_ri <- sp4 + getError(df_cor_rs[28:36,])

         sp5 <- mean(df_cor_rs[37:45,])
         sp5_li <- sp5 - getError(df_cor_rs[37:45,])
         sp5_ri <- sp5 + getError(df_cor_rs[37:45,])

         sp6 <- mean(df_cor_rs[46:54,])
         sp6_li <- sp6 - getError(df_cor_rs[46:54,])
         sp6_ri <- sp6 + getError(df_cor_rs[46:54,])

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

      DevSeq_GE_div <- rbind(root_div, hypocotyl_div, leaf_div, veg_apex_div, inf_apex_div, 
        flower_div, stamen_div, carpel_div, pollen_div)
      rownames(DevSeq_GE_div) <- NULL
      colnames(DevSeq_GE_div) <- c("correlation", "lower", "upper")

      DevSeq_div_rates <- data.frame(cbind(comp_spec, comp_organ, div_times, DevSeq_GE_div), 
        stringsAsFactors=FALSE)

      DevSeq_div_rates$div_times <- as.numeric(DevSeq_div_rates$div_times)
      DevSeq_div_rates$correlation <- as.numeric(DevSeq_div_rates$correlation)
      DevSeq_div_rates$lower <- as.numeric(DevSeq_div_rates$lower)
      DevSeq_div_rates$upper <- as.numeric(DevSeq_div_rates$upper)
      
      # Change order of organs in df
      DevSeq_div_rates <- DevSeq_div_rates[c(7:12,37:42,31:36,1:6,19:30,43:48,13:18,49:54),]
      DevSeq_div_rates_wo_pollen <- DevSeq_div_rates[1:48,]
      DevSeq_div_rates_pollen <- DevSeq_div_rates[49:54,]
      DevSeq_div_rates_wo_pollen$comp_organ <- factor(DevSeq_div_rates_wo_pollen$comp_organ, 
        levels = unique(DevSeq_div_rates_wo_pollen$comp_organ))
      DevSeq_div_rates_pollen$comp_organ <- factor(DevSeq_div_rates_pollen$comp_organ, 
        levels = unique(DevSeq_div_rates_pollen$comp_organ))

      

      # Make GE divergence plot
      makeGEDivPlot <- function(data1, data2, plot_title, coefficient, pos) {

        if (pos == "main") {

            fname <- sprintf('%s.jpg', paste("GE_divergence_rates_AL", coefficient, expr_estimation, pos, sep="_"))
            plot_wdt <- 12.535
            plot_hdt <- 8
            legend_x_pos <- 0.723
            legend_y_pos <- 0.88
            legend_key_s <- 0.95
            linewd <- 3.1
            brass_label <- "Brassicaceae"
            text_x1_pos <- 19.5
            text_x4_pos <- 158.6
            text_y_pos <- 0.4815
            x_poll_pos <- 147.47
            y_poll_pos <- 0.83575
            leg_ln_s1_x <- 135.55
            leg_ln_s1_xend <- 137.25
            leg_ln_s2_x <- 138.3
            leg_ln_s2_xend <- 139.975
            leg_ln_s_y <- 0.83575
            leg_box_bd <- 1.0

        } else {

            fname <- sprintf('%s.jpg', paste("GE_divergence_rates_AL", coefficient, expr_estimation, pos, sep="_"))
            plot_wdt <- 9.5 # condenced plot width for suppl
            plot_hdt <- 6.75 # condenced plot width for suppl
            legend_x_pos <- 0.6425
            legend_y_pos <- 0.8615
            legend_key_s <- 0.9
            linewd <- 2.85
            brass_label <- "Brassiceae"
            text_x1_pos <- 22
            text_x4_pos <- 157
            text_y_pos <- 0.4845
            x_poll_pos <- 145.5
            y_poll_pos <- 0.8223
            leg_ln_s1_x <- 129.3
            leg_ln_s1_xend <- 131.475
            leg_ln_s2_x <- 132.95
            leg_ln_s2_xend <- 135.15
            leg_ln_s_y <- 0.822
            leg_box_bd <- 0
        }

        

        p <- ggplot(data=data1, aes(x=div_times, y=correlation, group=comp_organ, colour=comp_organ)) + 
        geom_ribbon(aes(ymin = data1$lower, ymax = data1$upper, fill= comp_organ), alpha = 0.25, 
            linetype = 0, show.legend = FALSE) + 
        scale_fill_manual(values = c("Hypocotyl  "="#53b0db", "Stamen  "="#ee412e", "Flower  "="#e075af", 
                "Root  "="#6a54a9", "Apex veg  "="#96ba37", "Apex inf  "="#fad819", "Carpel  "="#f2a72f", 
                "Leaf  "="#2c8654")) + 
        geom_line(size = linewd) +  
        scale_x_continuous(limits = c(7,160), expand = c(0.02,0), breaks = c(7,9,25,46,106,160)) + 
        scale_y_continuous(limits = c(0.4675, 0.91), expand = c(0.02, 0)) + 
        scale_color_manual(values = c("#53b0db", "#ee412e", "#e075af", "#6a54a9", "#96ba37", "#fad819", 
            "#f2a72f", "#2c8654"), 
            # organ order: hypocotyl/stamen/flower/root/veg_apex/inf_apex/carpel/leaf
            breaks=c("Root  ", "Hypocotyl  ", "Leaf  ", "Apex veg  ", "Apex inf  ", "Flower  ", 
                "Stamen  ", "Carpel  ")) + 
        geom_line(aes(x=div_times, y=correlation), data=data2, color = "#a63126", lty = "22", 
            lwd = linewd) + # pollen
        annotate("text", x=x_poll_pos, y=y_poll_pos, label= "Pollen", size=7.56) + 
        geom_segment(x=leg_ln_s1_x, xend=leg_ln_s1_xend, y=leg_ln_s_y, yend=leg_ln_s_y, color="#a63126", size=linewd) + 
        geom_segment(x=leg_ln_s2_x, xend=leg_ln_s2_xend, y=leg_ln_s_y, yend=leg_ln_s_y, color="#a63126", size=linewd) + 
        geom_segment(x=157.5, xend=157.5, y=0.461, yend=0.49, color="white", size=12.5) + 
        annotate("text", x=text_x1_pos, y=text_y_pos, label= brass_label, size=8) + 
        annotate("text", x=46, y=text_y_pos, label= "TH", size=8) + 
        annotate("text", x=106, y=text_y_pos, label= "MT", size=8) + 
        annotate("text", x=text_x4_pos, y=text_y_pos, label= "BD", size=8) + 
        geom_segment(x=7, xend=7, y=0.435, yend=0.468, color="black", size=0.7) + 
        geom_segment(x=9, xend=9, y=0.435, yend=0.468, color="black", size=0.7) + 
        geom_segment(x=25, xend=25, y=0.435, yend=0.468, color="black", size=0.7) + 
        geom_segment(x=46, xend=46, y=0.435, yend=0.468, color="black", size=0.7) + 
        geom_segment(x=106, xend=106, y=0.435, yend=0.468, color="black", size=0.7) + 
        geom_segment(x=160, xend=160, y=0.435, yend=0.468, color="black", size=0.7) + 
        guides(color = guide_legend(ncol = 3))

        q <- p + theme_bw() + xlab("Divergence time from A.lyrata (Myr)") + ylab("Pearson's r w/ A.lyrata") + 
        theme(text=element_text(size=16), 
            axis.ticks.length=unit(0.35, "cm"), 
            axis.ticks = element_line(colour = "black", size = 0.7),  
            plot.margin = unit(c(0.55, 1.175, 0.5, 0.4),"cm"), 
            axis.title.y = element_text(size=25, margin = margin(t = 0, r = 15, b = 0, l = 11), colour="black"), 
            axis.title.x = element_text(size=25, margin = margin(t = 14.75, r = 0, b = 2, l = 0), colour="black"), 
            axis.text.x = element_text(size=21.25, angle=0, margin = margin(t = 5.5), colour="black"), 
            axis.text.y = element_text(size=21.25, angle=0, margin = margin(r = 5.5), colour="black"), 
            legend.box.background = element_rect(colour = "#d5d5d5", fill=NA, size=leg_box_bd), 
            panel.border = element_rect(colour = "black", fill=NA, size=1.75), 
            panel.grid.major = element_line(color="#d5d5d5"),
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

      makeGEDivPlot(data1 = DevSeq_div_rates_wo_pollen, data2 = DevSeq_div_rates_pollen, 
        coefficient = coefficient, pos = "main")
      makeGEDivPlot(data1 = DevSeq_div_rates_wo_pollen, data2 = DevSeq_div_rates_pollen, 
        coefficient = coefficient, pos = "ext")
   }
   
}


makeCompAnylsisAL(expr_estimation="TPM", coefficient="pearson", data_norm="inter-organ")
makeCompAnylsisAL(expr_estimation="counts", coefficient="pearson", data_norm="inter-organ")





