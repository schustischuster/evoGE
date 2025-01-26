# Fit non-linear models (NLMs) to angiosperm and mammalian organ transcriptome data 
# Prepare Brawand and DevSeq AL comparative ortholog gene expression data
# Thresholds: DevSeq 0.05 ERCC; Brawand 0.5 TPM (no ERCC spike-ins available)
# Data input: Brawand and DevSeq TPM expression tables of all samples



#-------------------------------------- Read data tables ---------------------------------------


getNLMs <- function(expr_estimation = c("TPM", "counts"), coefficient = c("pearson", "spearman")) {


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
    # return_objects <- getNLMs(expr_estimation="TPM", coefficient="pearson") # read in DevSeq expression data
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




#---------------- Get gene expression divergence rates for AL vs species X -----------------


   # Use pearson correlation, intra-organ normalization and TPM
   # Plot AL vs species X GE divergence rates for SI

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

        if (packageVersion("gplots") <  "3.0.0.2") {
            cvalues = c("#53b0db", "#ee412e", "#e075af", "#6a54a9", "#96ba37", "#fad819", 
            "#f2a72f", "#2c8654", "#a63126")
        } else {
            cvalues = c("#6a54a9", "#53b0db", "#2c8654", "#96ba37", "#fad819", "#e075af", 
            "#ee412e", "#f2a72f", "#a63126")
        }

        p <- ggplot(data=data1, aes(x=div_times, y=correlation, group=comp_organ, colour=comp_organ)) + 
        geom_ribbon(aes(ymin = data1$lower, ymax = data1$upper, fill= comp_organ), alpha = 0.25, 
            linetype = 0, show.legend = FALSE) + 
        scale_fill_manual(values = c("Hypocotyl  "="#53b0db", "Stamen  "="#ee412e", "Flower  "="#e075af", 
                "Root  "="#6a54a9", "Apex veg  "="#96ba37", "Apex inf  "="#fad819", "Carpel  "="#f2a72f", 
                "Leaf  "="#2c8654", "Pollen  "="#a63126")) + 
        geom_line(size = 2.9) +  
        scale_x_continuous(limits = c(5.5,161.5), expand = c(0.02,0), breaks = c(7,9,25,46,106,160), 
          labels = c( "7 ", " 9", 25, 46, 106, 160)) + 
        scale_y_continuous(limits = c(0.4425, 0.907), expand = c(0.02, 0)) + 
        scale_color_manual(values = cvalues, 
            # organ order: hypocotyl/stamen/flower/root/veg_apex/inf_apex/carpel/leaf
            breaks=c("Root  ", "Hypocotyl  ", "Leaf  ", "Apex veg  ", "Apex inf  ", "Flower  ", 
                "Stamen  ", "Carpel  ", "Pollen  ")) + 
        geom_line(aes(x=div_times, y=correlation), data=data2, color = "white", lty = "solid", 
            lwd = 2.9) + # pollen
        geom_line(aes(x=div_times, y=correlation), data=data2, color = "#a63126", lty = "22", 
            lwd = 2.9) + # pollen
        geom_segment(x=156.75, xend=156.75, y=0.4475, yend=0.4775, color="white", size=12.5) + 
        annotate("text", x=20.5, y=0.46525, label= "Brassiceae", size=8) + 
        annotate("text", x=46, y=0.46525, label= "TH", size=8) + 
        annotate("text", x=106, y=0.46525, label= "MT", size=8) + 
        annotate("text", x=157, y=0.46525, label= "BD", size=8) + 
        geom_segment(x=7, xend=7, y=0.425, yend=0.4475, color="black", size=1.1) + 
        geom_segment(x=9, xend=9, y=0.425, yend=0.4475, color="black", size=1.1) + 
        geom_segment(x=25, xend=25, y=0.425, yend=0.4475, color="black", size=1.1) + 
        geom_segment(x=46, xend=46, y=0.425, yend=0.4475, color="black", size=1.1) + 
        geom_segment(x=106, xend=106, y=0.425, yend=0.4475, color="black", size=1.1) + 
        geom_segment(x=160, xend=160, y=0.425, yend=0.4475, color="black", size=1.1) + 
        guides(color = guide_legend(ncol = 3))

        q <- p + theme_bw() + xlab("Divergence time from A.lyrata (Myr)") + ylab("Pearson's r w/ A.lyrata") + 
        theme(text=element_text(size=16), 
            axis.ticks.length=unit(0.325, "cm"), 
            axis.ticks = element_line(colour = "black", size = 1.15),  
            plot.margin = unit(c(0.55, 0.46, 4.26, 1.27),"cm"), 
            axis.title.y = element_text(size=24.6, margin = margin(t = 0, r = 13, b = 0, l = 10.8), colour="black", 
                face = "bold"), 
            axis.title.x = element_text(size=24.6, margin = margin(t = 13.75, r = 0, b = 1, l = 0), colour="black", 
                face = "bold"), 
            axis.text.x = element_text(size=21.75, angle=0, margin = margin(t = 4.5), colour="grey5"), 
            axis.text.y = element_text(size=21.75, angle=0, margin = margin(r = 2.5), colour="grey5"), 
            legend.box.background = element_rect(colour = "#d5d5d5", fill=NA, size=0), 
            panel.border = element_rect(colour = "black", fill=NA, size=2.3), 
            panel.grid.major = element_blank(),
            panel.grid.minor.x = element_blank(), 
            panel.grid.minor.y = element_blank(), 
            legend.position = c(0.671, 0.865), 
            legend.title = element_blank(), 
            legend.text = element_text(size=22), 
            legend.spacing.x = unit(0.5, 'cm'), 
            legend.key.size = unit(0.9, "cm"), 
            legend.background=element_blank()) 

        ggsave(file = file.path(out_dir, "output", "plots", fname), plot = q, 
            width = 10.5, height = 8.5, dpi = 300, units = c("in"), limitsize = FALSE) 
      }

      makeGEDivPlot(data1 = DevSeq_AL_div_rates, data2 = DevSeq_AL_div_rates_pollen, 
        coefficient = coefficient)
    }




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
        a = 0.3, b = 0.5, c = -0.01), data = compDivRates11[1:6,])
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
      # columns for plotting non-linear regressions for individual organs (facets)
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

      DS_AT_sOU_v_root_nl_df <- as.data.frame(formatNLM.table(DS_AT_sOU_v_root_nl))
      DS_AT_sOU_v_hypo_nl_df <- as.data.frame(formatNLM.table(DS_AT_sOU_v_hypo_nl))
      DS_AT_sOU_v_leaf_nl_df <- as.data.frame(formatNLM.table(DS_AT_sOU_v_leaf_nl))
      DS_AT_sOU_v_apex_veg_nl_df <- as.data.frame(formatNLM.table(DS_AT_sOU_v_apex_veg_nl))
      DS_AT_sOU_v_apex_inf_nl_df <- as.data.frame(formatNLM.table(DS_AT_sOU_v_apex_inf_nl))
      DS_AT_sOU_v_flower_nl_df <- as.data.frame(formatNLM.table(DS_AT_sOU_v_flower_nl))
      DS_AT_sOU_v_stamen_nl_df <- as.data.frame(formatNLM.table(DS_AT_sOU_v_stamen_nl))
      DS_AT_sOU_v_carpel_nl_df <- as.data.frame(formatNLM.table(DS_AT_sOU_v_carpel_nl))

      DS_AT_sOU_v_nlm_coord <- rbind(DS_AT_sOU_v_root_nl_df, DS_AT_sOU_v_hypo_nl_df, 
        DS_AT_sOU_v_leaf_nl_df, DS_AT_sOU_v_apex_veg_nl_df, DS_AT_sOU_v_apex_inf_nl_df, 
        DS_AT_sOU_v_flower_nl_df, DS_AT_sOU_v_stamen_nl_df, DS_AT_sOU_v_carpel_nl_df)


      Br11_pea_dist_brain_nl_df <- as.data.frame(formatNLM.table(Br11_pea_dist_brain_nl))
      Br11_pea_dist_cereb_nl_df <- as.data.frame(formatNLM.table(Br11_pea_dist_cereb_nl))
      Br11_pea_dist_heart_nl_df <- as.data.frame(formatNLM.table(Br11_pea_dist_heart_nl))
      Br11_pea_dist_kidney_nl_df <- as.data.frame(formatNLM.table(Br11_pea_dist_kidney_nl))
      Br11_pea_dist_liver_nl_df <- as.data.frame(formatNLM.table(Br11_pea_dist_liver_nl))
      Br11_pea_dist_testis_nl_df <- as.data.frame(formatNLM.table(Br11_pea_dist_testis_nl))

      Br11_pea_dist_nlm_coord <- rbind(Br11_pea_dist_brain_nl_df, Br11_pea_dist_cereb_nl_df, 
        Br11_pea_dist_heart_nl_df, Br11_pea_dist_kidney_nl_df, Br11_pea_dist_liver_nl_df, 
        Br11_pea_dist_testis_nl_df)

      Br11_sOU_v_brain_nl_df <- as.data.frame(formatNLM.table(Br11_sOU_v_brain_nl))
      Br11_sOU_v_cereb_nl_df <- as.data.frame(formatNLM.table(Br11_sOU_v_cereb_nl))
      Br11_sOU_v_heart_nl_df <- as.data.frame(formatNLM.table(Br11_sOU_v_heart_nl))
      Br11_sOU_v_kidney_nl_df <- as.data.frame(formatNLM.table(Br11_sOU_v_kidney_nl))
      Br11_sOU_v_liver_nl_df <- as.data.frame(formatNLM.table(Br11_sOU_v_liver_nl))
      Br11_sOU_v_testis_nl_df <- as.data.frame(formatNLM.table(Br11_sOU_v_testis_nl))

      Br11_sOU_v_nlm_coord <- rbind(Br11_sOU_v_brain_nl_df, Br11_sOU_v_cereb_nl_df, 
        Br11_sOU_v_heart_nl_df, Br11_sOU_v_kidney_nl_df, Br11_sOU_v_liver_nl_df, 
        Br11_sOU_v_testis_nl_df)




#-- Get cumulative slope values of mean organ regressions for pea and sOU non-linear models --


      # Define factors to get cumsum of DS slopes to 1
      cum_fact_pea <- 4.225954
      cum_fact_sOU <- 0.7859427


      # Get cumulative slope values for non-linear model
      getCum.NL.Slope <- function(x, cfact){

        fname <- deparse(substitute(x))

        dataset_id <- gsub( "_.*$", "", fname)

        if (dataset_id == "DS") {

          x_grid <- x_DS_grid

        } else x_grid <- x_Br_grid

        x <- as.numeric(x[,1])
        slopes = diff(x)/diff(x_grid)
        slopes_cumsum <- cumsum(as.data.frame(slopes)[,1])
        slopes_cumsum <- as.numeric(as.data.frame(slopes_cumsum)[,1])

        slopes_cumsum <- slopes_cumsum * cfact

        return(slopes_cumsum)
      }


      # Get cumsums for individual organs
      DS_AT_pea_nl_cumslopes <- do.call(cbind, lapply(DS_AT_pea_dist_nl_list, getCum.NL.Slope, 
        cfact=cum_fact_pea))
      Br11_pea_nl_cumslopes <- do.call(cbind, lapply(Br11_pea_dist_nl_list, getCum.NL.Slope, 
        cfact=cum_fact_pea))
      Br_pea_nl_cumslopes <- do.call(cbind, lapply(Br_pea_dist_nl_list, getCum.NL.Slope, 
        cfact=cum_fact_pea))
      DS_AT_sOU_v_nl_cumslopes <- do.call(cbind, lapply(DS_AT_sOU_v_nl_list, getCum.NL.Slope, 
        cfact=cum_fact_sOU))
      Br11_sOU_v_nl_cumslopes <- do.call(cbind, lapply(Br11_sOU_v_nl_list, getCum.NL.Slope, 
        cfact=cum_fact_sOU))
      Br_sOU_v_nl_cumslopes <- do.call(cbind, lapply(Br_sOU_v_nl_list, getCum.NL.Slope, 
        cfact=cum_fact_sOU))


      # Get stats
      getLoessStats <- function(loess_data) {

        nsamples <- ncol(loess_data)

        loess_mean <- rowMeans(loess_data)

        sd_out <- by(loess_data, 1:nrow(loess_data), function(row) sd <- sd(row))
        sd_out <- sd_out[1:length(sd_out)]
        sd_out <- as.data.frame(sd_out)
        sd_out <- as.numeric(sd_out[,1])

        error <- qt(0.975, df = nsamples-1) * sd_out/sqrt(nsamples)

        mean_li <- loess_mean - error
        mean_ri <- loess_mean + error

        loess_stat <- data.frame(correlation = loess_mean, li = mean_li, ri = mean_ri)

        return(loess_stat)

      }

      DS_AT_pea_nl_mean_cumslopes <- getLoessStats(DS_AT_pea_nl_cumslopes)
      Br11_pea_nl_mean_cumslopes <- getLoessStats(Br11_pea_nl_cumslopes)
      Br_pea_nl_mean_cumslopes <- getLoessStats(Br_pea_nl_cumslopes)
      DS_AT_sOU_v_nl_mean_cumslopes <- getLoessStats(DS_AT_sOU_v_nl_cumslopes)
      Br11_sOU_v_nl_mean_cumslopes <- getLoessStats(Br11_sOU_v_nl_cumslopes)
      Br_sOU_v_nl_mean_cumslopes <- getLoessStats(Br_sOU_v_nl_cumslopes)


      # Prepare data for ggplot2
      x_DS_grid_cum <- seq(7.1, 160, length = 199)  ## prediction grid
      x_Br_grid_cum <- seq(6.7, 159, length = 199)  ## prediction grid

      div_times <- as.numeric(data.frame(div_times=c(x_DS_grid_cum, x_Br_grid_cum, 
        x_Br_grid_cum))[,1])
      dataset <- factor(data.frame(dataset=rep(c("Angiosperms", "Mammals.11", "Mammals.ra"), 
        each=199))[,1])

      DS_Br11_Br_pea_nlm_cum <- rbind(DS_AT_pea_nl_mean_cumslopes, Br11_pea_nl_mean_cumslopes, 
        Br_pea_nl_mean_cumslopes)
      DS_Br11_Br_sOU_nlm_cum <- rbind(DS_AT_sOU_v_nl_mean_cumslopes, Br11_sOU_v_nl_mean_cumslopes, 
        Br_sOU_v_nl_mean_cumslopes)
      DS_Br11_Br_pea_nlm_cum <- data.frame(div_times, DS_Br11_Br_pea_nlm_cum, dataset)
      DS_Br11_Br_sOU_nlm_cum <- data.frame(div_times, DS_Br11_Br_sOU_nlm_cum, dataset)


      
      # Make sOU GE divergence plot showing cumulative mean slope values for SI
      plotsOUPeaCumSlopes <- function(data, coefficient, expr_estimation) {

        fname <- sprintf('%s.jpg', paste(deparse(substitute(data))))

        col_breaks <- factor(c("Angiosperms", "Mammals.11", "Mammals.ra"), 
            levels=c("Angiosperms", "Mammals.11", "Mammals.ra"))
        y_breaks <- c(0,0.2,0.4,0.6,0.8,1,1.2,1.4)

        col_scale <- c('#728acb', 'red', 'red3')
        fill_scale <- c('#728acb', 'red', 'red3')

        p <- ggplot(data=data, aes(x = div_times, y = correlation, group = dataset)) + 
        geom_ribbon(aes(ymin = data$li, ymax = data$ri, fill = dataset), alpha = 0.125, 
          linetype = 0, show.legend = FALSE) + 
        geom_line(size = 2.75, data = data, aes(x = div_times, y = correlation, group = dataset, 
          colour = dataset)) + 
        geom_line(data = data[1:159,], aes(x = div_times, y = correlation), color = '#728acb', size = 2.5) +
        scale_x_continuous(limits = c(0,161.25), expand = c(0.02,0), breaks = c(0,25,50,75,100,125,150)) + 
        scale_y_continuous(limits = c(-0.035, 1.3585), expand = c(0.01, 0), breaks = y_breaks) + 
        scale_color_manual(values = col_scale, breaks = col_breaks) + 
        scale_fill_manual(values = fill_scale) + 
        scale_size(range = c(0.5, 12)) + 
        guides(color = guide_legend(ncol=2, keywidth = 0.4, keyheight = 0.4, default.unit = "inch", 
          title = ""))

        q <- p + theme_bw() + xlab("Divergence time (Myr)") + ylab("Cumulative mean \n slope value") + 
        theme(text=element_text(size = 16), 
          panel.border = element_rect(colour = "white", fill=NA, size = 2.3), 
          axis.line = element_line(colour = 'black', size = 1.15), 
          axis.ticks.length = unit(0.325, "cm"), 
          axis.ticks = element_line(colour = "black", size = 1.15), 
          plot.margin = unit(c(1, 0.932, 4.8575, 0.025),"cm"), 
          axis.title.y = element_text(size=24.6, margin = margin(t = 0, r = 14, b = 0, l = 0), 
            colour="black", face = "bold"), 
          axis.title.x = element_text(size=24.6, margin = margin(t = 9.25, r = 0, b = 7.5, l = 0), 
            colour="black", face = "bold"), 
          axis.text.x = element_text(size=21.75, margin = margin(t = 3.5, b = 8), colour="grey5"), 
          axis.text.y = element_text(size=21.75, angle=0, margin = margin(l = 2.5, r = 2.5), colour="grey5"), 
          panel.spacing = unit(0.15, "cm"), 
          panel.grid.major = element_blank(),
          panel.grid.minor.x = element_blank(), 
          panel.grid.minor.y = element_blank(), 
          legend.position = c(0.507, 0.856), 
          legend.title = element_blank(), 
          legend.text = element_text(size=21.75), 
          legend.spacing.x = unit(0.5, 'cm'), 
          legend.key.size = unit(0.95, "cm"), 
          legend.background=element_blank()) 

        ggsave(file = file.path(out_dir, "output", "plots", fname), plot = q, 
          width = 6.5, height = 6.5, dpi = 300, units = c("in"), limitsize = FALSE)
      }

      plotsOUPeaCumSlopes(data = DS_Br11_Br_pea_nlm_cum)
      plotsOUPeaCumSlopes(data = DS_Br11_Br_sOU_nlm_cum)




#---- Compute mean slope value of non-linear regression for all DS and Br/Br11 organs -----


      # Compute mean slope value for each organ for non-linear regression
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

      DevSeq_AT_pea_dist_nl_slopes <- as.data.frame(do.call(rbind, lapply(DS_AT_pea_dist_nl_list, getDS.NL.Slope)))
      DevSeq_AT_sOU_v_nl_slopes <- as.data.frame(do.call(rbind, lapply(DS_AT_sOU_v_nl_list, getDS.NL.Slope)))
      Br11_sOU_v_nl_slopes <- as.data.frame(do.call(rbind, lapply(Br11_sOU_v_nl_list, getDS.NL.Slope)))
      Br11_pea_dist_nl_slopes <- as.data.frame(do.call(rbind, lapply(Br11_pea_dist_nl_list, getDS.NL.Slope)))
      Br_pea_dist_nl_slopes <- as.data.frame(do.call(rbind, lapply(Br_pea_dist_nl_list, getDS.NL.Slope)))
      Br_sOU_v_nl_slopes <- as.data.frame(do.call(rbind, lapply(Br_sOU_v_nl_list, getDS.NL.Slope)))


      # Add organ id column
      formatNL.Slope <- function(x, regr_mod) {

        names(x) <- regr_mod

        if ((regr_mod == "DS_AT_pea_nlm") || (regr_mod == "DS_AT_sOU_nlm")) {

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
      Br_nlm_slopes_AT <- Br_nlm_slopes
      colnames(Br_nlm_slopes_AT) <- c("sample", "DS_AT_Br11_sOU_nlm", "DS_AT_Br_sOU_nlm", 
        "DS_AT_Br11_pea_nlm", "DS_AT_Br_pea_nlm")
      DS_AT_Br_nlm_slopes <- rbind(DevSeq_AT_nlm_slopes, Br_nlm_slopes_AT)


      # Save Angiosperm and mammalian nlm slope values
      message("Writing data tables...")

      slopes_regr_list <- list(DS_AT_Br_nlm_slopes = DS_AT_Br_nlm_slopes)

      for(i in names(slopes_regr_list)) { 
        write.table(slopes_regr_list[[i]], file = file.path(out_dir, "output", "data", paste0(i,".txt")), 
          sep="\t", col.names=TRUE, row.names = FALSE, dec=".", quote = FALSE)
      }




#--- Apply LOESS to sOU and pea dist expression data and compute cumulative slope values ----


      # Set up lists containing sOU expression distances
      brawandSouV11_organ_lst <- list(compSouVDivRates11[49:54,], compSouVDivRates11[55:60,], compSouVDivRates11[61:66,], 
        compSouVDivRates11[67:72,],compSouVDivRates11[73:78,], compSouVDivRates11[79:83,])

      devseq_AT_SouV_organ_lst <- list(compSouVDivRates11[1:6,], compSouVDivRates11[7:12,], compSouVDivRates11[13:18,], 
        compSouVDivRates11[19:24,], compSouVDivRates11[25:30,], compSouVDivRates11[31:36,], compSouVDivRates11[37:42,], 
        compSouVDivRates11[43:48,])

      brawandSouV_organ_lst <- list(compSouVDivRates[49:54,], compSouVDivRates[55:60,], compSouVDivRates[61:66,], 
        compSouVDivRates[67:72,],compSouVDivRates[73:78,], compSouVDivRates[79:83,])


      # Get cumulative slope values for LOESS regression
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

        error <- qt(0.975, df = nsamples-1) * sd_out/sqrt(nsamples)

        mean_li <- loess_mean - error
        mean_ri <- loess_mean + error

        loess_stat <- data.frame(correlation = loess_mean, li = mean_li, ri = mean_ri)

        return(loess_stat)
      }


      # Get mean value of slopes
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

      div_times <- c(x_DS_grid_cum, x_Br_grid_cum, x_Br_grid_cum)
      div_times <- as.data.frame(div_times)
      colnames(div_times) <- "div_times"
      dataset <- data.frame(rep(c("Angiosperms", "Mammals.11", "Mammals.ra"), each=159))
      colnames(dataset) <- "dataset"
      correlation <- rbind(DevSeqSouV_AT_loess_mean, brawandSouV11_loess_mean, 
        brawandSouV_loess_mean)
      sOU_v_loess_cum_slopes <- data.frame(div_times, correlation, dataset)

      sOU_v_loess_cum_slopes$div_times <- as.numeric(sOU_v_loess_cum_slopes$div_times)
      sOU_v_loess_cum_slopes$correlation <- as.numeric(sOU_v_loess_cum_slopes$correlation)
      sOU_v_loess_cum_slopes$li <- as.numeric(sOU_v_loess_cum_slopes$li)
      sOU_v_loess_cum_slopes$ri <- as.numeric(sOU_v_loess_cum_slopes$ri)
      sOU_v_loess_cum_slopes$dataset <- factor(sOU_v_loess_cum_slopes$dataset)




   # Make sOU GE divergence plot showing cumulative mean slope values for SI
   plotCumSlopes <- function(data, coefficient, expr_estimation) {


      fname <- sprintf('%s.jpg', paste("compSouV_loess_cum_slopes", coefficient, expr_estimation, sep="_"))
        
      col_breaks <- c("Angiosperms", "Mammals.11", "Mammals.ra")
      y_breaks <- c(0,0.2,0.4,0.6,0.8,1,1.2,1.4)

      col_scale <- c('#728acb', 'red', 'red3')
      fill_scale <- c('#728acb', 'red', 'red3')


      p <- ggplot(data=data, aes(x = div_times, y = correlation, group = dataset)) + 
      geom_ribbon(aes(ymin = data$li, ymax = data$ri, fill = dataset), alpha = 0.088, 
            linetype = 0, show.legend = FALSE) + 
      geom_line(size = 2.75, data = data, aes(x = div_times, y = correlation, group = dataset, 
        colour = dataset)) + 
      geom_line(data = data[1:159,], aes(x = div_times, y = correlation), color = '#728acb', size = 2.5) +
      scale_x_continuous(limits = c(0,162), expand = c(0.02,0), breaks = c(0,20,40,60,80,100,120,140,160)) + 
      scale_y_continuous(limits = c(-0.035, 1.22), expand = c(0.01, 0), breaks = y_breaks) + 
      scale_color_manual(values = col_scale, breaks = col_breaks) + 
      scale_fill_manual(values = fill_scale) + 
      scale_size(range = c(0.5, 12)) + 
      guides(color = guide_legend(ncol=1, keywidth = 0.4, keyheight = 0.4, default.unit = "inch", 
        title = ""))

      q <- p + theme_bw() + xlab("Divergence time (Myr)") + ylab("Cumulative mean slope value") + 
      theme(text=element_text(size=16), 
        axis.ticks.length=unit(0.35, "cm"), 
        axis.ticks = element_line(colour = "black", size = 1.25),  
        plot.margin = unit(c(0.55, 1.175, 1.55, 0.4),"cm"), 
        axis.title.y = element_text(size=24.5, margin = margin(t = 0, r = 15, b = 0, l = 11), colour="black", 
            face = "bold"), 
        axis.title.x = element_text(size=24.5, margin = margin(t = 15.75, r = 0, b = 1, l = 0), colour="black", 
            face = "bold"), 
        axis.text.x = element_text(size=21.75, angle=0, margin = margin(t = 4), colour="black"), 
        axis.text.y = element_text(size=21.75, angle=0, margin = margin(r = 4), colour="black"), 
        legend.box.background = element_rect(colour = NA, fill= "white" , size=1.2), 
        panel.border = element_rect(colour = "black", fill=NA, size=2.4), 
        panel.grid.major = element_blank(),
        panel.grid.minor.x = element_blank(), 
        panel.grid.minor.y = element_blank(), 
        legend.position = c(0.18, 0.84), 
        legend.title = element_blank(), 
        legend.text = element_text(size=22), 
        legend.spacing.x = unit(0.5, 'cm'), 
        legend.key.size = unit(0.95, "cm"), 
        legend.background=element_blank()) 

        ggsave(file = file.path(out_dir, "output", "plots", fname), plot = q, 
            width = 9.5, height = 6.75, dpi = 300, units = c("in"), limitsize = FALSE) 
  }

  plotCumSlopes(data = sOU_v_loess_cum_slopes, coefficient = coefficient, expr_estimation = expr_estimation)




#-------- Plot Angiosperm and Mammalian nlm's for metric pearson and sOU-v distance ---------


  # Change organ names for facet strip

  formDsDistMod <- function(dist_df) {

     dist_df$comp_organ <- dist_df$comp_organ %<>% 
        gsub("root", "Root", .) %>% 
        gsub("hypo", "Hypocotyl", .) %>% 
        gsub("leaf", "Leaf", .) %>% 
        gsub("veg", "Apex veg", .) %>% 
        gsub("inf", "Apex inf", .) %>% 
        gsub("flower", "Flower", .) %>% 
        gsub("stamen", "Stamen", .) %>% 
        gsub("carpel", "Carpel", .)

     dist_df$comp_organ <- factor(dist_df$comp_organ, 
        levels=c("Root","Hypocotyl","Leaf","Apex veg","Apex inf","Flower","Stamen","Carpel"))

     return(dist_df)
  }

  DS_AT_pea_dist_nlm_coord_df <- formDsDistMod(DS_AT_pea_dist_nlm_coord)
  DS_AT_sOU_v_nlm_coord_df <- formDsDistMod(DS_AT_sOU_v_nlm_coord)


  formBrDistMod <- function(dist_df) {

     dist_df$comp_organ <- dist_df$comp_organ %<>% 
        gsub("brain", "Brain", .) %>% 
        gsub("cereb", "Cerebellum", .) %>% 
        gsub("heart", "Heart", .) %>% 
        gsub("kidney", "Kidney", .) %>% 
        gsub("liver", "Liver", .) %>% 
        gsub("testis", "Testis", .)

     dist_df$comp_organ <- factor(dist_df$comp_organ, 
        levels=c("Brain", "Cerebellum", "Heart", "Kidney", "Liver", "Testis"))

     return(dist_df)
  }

  Br11_pea_dist_nlm_coord_df <- formBrDistMod(Br11_pea_dist_nlm_coord)
  Br11_sOU_v_nlm_coord_df <- formBrDistMod(Br11_sOU_v_nlm_coord)


  # Format pea and sOU-v dist tables for DS and Br11 data
  pea_dist_organs_DS_df <- compDivRates11[1:48,]
  pea_dist_organs_Br11_df <- compDivRates11[49:nrow(compDivRates11),]
  sOU_v_dist_organs_DS_df <- compSouVDivRates11[1:48,]
  sOU_v_dist_organs_Br11_df <- compSouVDivRates11[49:nrow(compSouVDivRates11),]
  sOU_v_dist_organs_DS_df$comp_organ <- sOU_v_dist_organs_DS_df$comp_organ %<>% 
        gsub("vegApex", "Apex veg", .) %>% 
        gsub("infApex", "Apex inf", .)
  sOU_v_dist_organs_Br11_df$comp_organ <- sOU_v_dist_organs_Br11_df$comp_organ %<>% 
        gsub("br", "Brain", .) %>% 
        gsub("cb", "Cerebellum", .) %>% 
        gsub("ht", "Heart", .) %>% 
        gsub("kd", "Kidney", .) %>% 
        gsub("lv", "Liver", .) %>% 
        gsub("ts", "Testis", .)
  pea_dist_organs_DS_df$comp_organ <- factor(pea_dist_organs_DS_df$comp_organ, 
    levels=c("Root", "Hypocotyl", "Leaf", "Apex veg", "Apex inf", "Flower", "Stamen", "Carpel"))
  sOU_v_dist_organs_DS_df$comp_organ <- factor(sOU_v_dist_organs_DS_df$comp_organ, 
    levels=c("Root", "Hypocotyl", "Leaf", "Apex veg", "Apex inf", "Flower", "Stamen", "Carpel"))
  pea_dist_organs_Br11_df$comp_organ <- factor(pea_dist_organs_Br11_df$comp_organ, 
    levels=c("Brain", "Cerebellum", "Heart", "Kidney", "Liver", "Testis"))
  sOU_v_dist_organs_Br11_df$comp_organ <- factor(sOU_v_dist_organs_Br11_df$comp_organ, 
    levels=c("Brain", "Cerebellum", "Heart", "Kidney", "Liver", "Testis"))



  plotDS.Br.pea.NLM <- function(data, data2, data_set) {

    fname <- sprintf('%s.jpg', paste(data_set, "nlm_regression_slopes", sep="_"))

    data2$error[is.na(data2$error)] <- 0

    lower <- as.numeric(data.frame(lower=data2$correlation-data2$error)[,1])
    upper <- as.numeric(data.frame(lower=data2$correlation+data2$error)[,1])

    data2 <- cbind(data2, lower, upper)

    if ((data_set == "Angiosperm_pea") || (data_set == "Angiosperm_sOU")) {

      spec_col <- '#728acb'
      spec_label <- "Angiosperms"
      spec_shape <- 16
      margin_r <- 1

    } else { 

      spec_col <- 'red'
      spec_label <- "Mammals.11"
      spec_shape <- 17
      margin_r <- 18.397

    }

    # Define y axis breaks and labels
    if (data_set == "Angiosperm_pea") {
      y_title <- "Pearson distance"
      marg_r <- 15.2
      marg_l <- 10.8
      p_space <- 0.15
    } else if (data_set == "Angiosperm_sOU") {
      y_title <- "Expression distance"
      marg_r <- 13.0
      marg_l <- 23.0
      p_space <- 0.575
    } else if (data_set == "Mammalian_pea") {
      y_title <- "Pearson distance"
      marg_r <- 15.2
      marg_l <- 10.8
      p_space <- 0.15
    } else if (data_set == "Mammalian_sOU") {
      y_title <- "Expression distance"
      marg_r <- 13.0
      marg_l <- 23.0
      p_space <- 0.575
    }

    p <- ggplot(data=data, color = dataset, aes(x=div_times, y=correlation)) + 
    geom_line(size = 2.9, data = data, aes(x = div_times, y = correlation, group = dataset, 
        colour = spec_col)) + 
    geom_point(data=data2, shape = spec_shape, aes(color = spec_col, stroke = 7.85)) + 
    geom_pointrange(data = data2, shape = spec_shape, mapping = aes(x=div_times, y=correlation, 
      ymin = lower, ymax = upper), 
      fatten = 0.1, size = 1.75, color=spec_col) + 
    scale_color_manual(labels=c(spec_label) , values=c(spec_col), guide = "legend") + 
    scale_fill_manual(values=c(spec_col), guide = "legend") + 
    scale_y_continuous(expand = c(0.1, 0), breaks = pretty_breaks()) + 
    scale_x_continuous(expand = c(0.1, 0), breaks=c(0, 50, 100, 150)) + 
    guides(shape = guide_legend(override.aes = list(stroke = 7.75)))

    q <- p + theme_classic() + xlab("Divergence time (Myr)") + ylab(y_title) + 
    theme(text=element_text(size = 16), 
      strip.text = element_text(size = 23.75), 
      strip.text.x = element_text(margin = margin(0.44, 0, 0.44, 0, "cm")), 
      strip.background = element_rect(colour = 'black', fill = NA, size = 2.4), 
      axis.ticks.length = unit(0.325, "cm"), 
      axis.ticks = element_line(colour = "black", size = 1.15), 
      axis.line = element_line(colour = 'black', size = 1.15), 
      plot.margin = unit(c(1, margin_r, 4.85, 0.4),"cm"), 
      axis.title.y = element_text(size=24.6, margin = margin(t = 0, r = marg_r, b = 0, l = marg_l), 
        colour="black", face = "bold"), 
      axis.title.x = element_text(size=24.6, margin = margin(t = 9.25, r = 0, b = 7.5, l = 0), 
        colour="black", face = "bold"), 
      axis.text.x = element_text(size=21.75, margin = margin(t = 3.5, b = 8), colour="grey5"), 
      axis.text.y = element_text(size=21.75, angle=0, margin = margin(l = 2.5, r = 2.5), colour="grey5"), 
      panel.spacing = unit(p_space, "cm"), 
      panel.grid.major = element_blank(),
      panel.grid.minor.x = element_blank(), 
      panel.grid.minor.y = element_blank(), 
      legend.position = "none") 

    q <- q + facet_wrap(~ comp_organ, nrow = 1, scales = "free")

    ggsave(file = file.path(out_dir, "output", "plots", fname), plot = q, 
      width = 28.5, height = 6.5, dpi = 300, units = c("in"), limitsize = FALSE) 
  }

  plotDS.Br.pea.NLM(data = DS_AT_pea_dist_nlm_coord_df, data2 = pea_dist_organs_DS_df, data_set = "Angiosperm_pea")
  plotDS.Br.pea.NLM(data = DS_AT_sOU_v_nlm_coord_df, data2 = sOU_v_dist_organs_DS_df, data_set = "Angiosperm_sOU")
  plotDS.Br.pea.NLM(data = Br11_pea_dist_nlm_coord_df, data2 = pea_dist_organs_Br11_df, data_set = "Mammalian_pea")
  plotDS.Br.pea.NLM(data = Br11_sOU_v_nlm_coord_df, data2 = sOU_v_dist_organs_Br11_df, data_set = "Mammalian_sOU")


  }
   
}








