# Prepare Brawand and DevSeq AL comparative ortholog gene expression data
# Thresholds: DevSeq 0.05 ERCC; Brawand 0.5 TPM (no ERCC spike-ins available)
# Data input: Brawand and DevSeq TPM expression tables of all samples



#-------------------------------------- Read data tables ---------------------------------------


makeCompAnylsisAL <- function(expr_estimation = c("TPM", "counts"), coefficient = c("pearson", "spearman")) {


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


    # Stop function here to allow specific analysis of a single data set
    # return_list <- list("expr_estimation" = expr_estimation, "x" = x, "coefficient" = coefficient, "col_names" = col_names)
    # return(return_list)
    # }
    # return_objects <- makeCompAnylsisAL(expr_estimation="TPM", coefficient="pearson") # read in DevSeq expression data
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




#-------------------- Read taxa objects for DevSeq AL data with replicates --------------------


    if (is.element("TPM", expr_estimation)) {

   
        # Construc taxa object
        x_DS_AL_taxa_objects = TEconstruct(ExpValueFP = file.path(out_dir, 
            "output", "data", 'x_DS_AL_taxobj_input.txt'), taxa = "all", subtaxa = 'all')

        x_Br2011_taxa_objects = TEconstruct(ExpValueFP = file.path(out_dir, 
            "output", "data", 'x_Br2011_taxobj_input.txt'), taxa = "all", subtaxa = 'all')


        DevSeq_AL_organ_list <- list("Root", "Hypocotyl", "Leaf", "vegApex", "infApex", "Flower", 
            "Stamen", "Carpel", "Pollen")

        Brawand2011_organ_list <- list("br", "cb", "ht", "kd", "lv", "ts")


        # Function to apply extended OU model with dynamic expression optimum ("variable-µ method")
        getExtOU <- function(organ, taxa_obj, samples) {

            sou_v_out <- expdist(taxa_obj, taxa = "all",
                subtaxa = organ,
                method = "sou_v")

            sou_v_pi <- sou_v_out$pi
            sou_v_distance <- as.data.frame(sou_v_out$distance)
            sou_v_distance_div <- as.data.frame(sou_v_distance[,1])
            sou_v_distance_div <- rbind(sou_v_distance_div, sou_v_pi)

            spec_id <- c(sub("\\_.*", "", rownames(sou_v_distance)), "pi")
            organ_id <- unique(sub(".*_", "", rownames(sou_v_distance)))

            if ((organ_id == "ts") && (samples == "sel")) {

                ppy_testis <- "NA"
                sou_v_distance_div <- as.data.frame(c(sou_v_distance_div[1:3,], ppy_testis, 
                    sou_v_distance_div[4:7,]), stringsAsFactors = FALSE)
                spec_id <- c("hsa", "ppa", "ggo", "ppy", "mml", "mmu", "mdo", "pi")
            
            } 

            if ((organ_id == "ts") && (samples == "all")) {

                ppy_testis <- "NA"
                sou_v_distance_div <- as.data.frame(c(sou_v_distance_div[1:3,], ppy_testis, 
                    sou_v_distance_div[4:9,]), stringsAsFactors = FALSE)
                spec_id <- c("hsa", "ppa", "ggo", "ppy", "mml", "mmu", "mdo", "oan", "gga", "pi")

            }

            if (organ_id == "testis") {

                Orangutan_testis <- "NA"
                sou_v_distance_div <- as.data.frame(c(sou_v_distance_div[1:3,], Orangutan_testis, 
                    sou_v_distance_div[4:7,]), stringsAsFactors = FALSE)
                spec_id <- c("Human", "Bonobo", "Gorilla", "Orangutan", "Macaque", "Mouse", 
                    "Opossum", "pi")
            }

            rownames(sou_v_distance_div) <- spec_id
            colnames(sou_v_distance_div) <- organ_id

            return(sou_v_distance_div)
        }


        DS_AL_sou_v <- as.data.frame(do.call(cbind, lapply(DevSeq_AL_organ_list, getExtOU,
            taxa_obj = x_DS_AL_taxa_objects)))
        rows_to_remove_DS_AL <- "ALY"
        DS_AL_sou_v <- DS_AL_sou_v[!(row.names(DS_AL_sou_v) %in% rows_to_remove_DS_AL), ]


        Br2011_sou_v <- as.data.frame(do.call(cbind, lapply(Brawand2011_organ_list, getExtOU,
            taxa_obj = x_Br2011_taxa_objects, samples = "sel")))
        rows_to_remove_Br2011 <- "hsa"
        Br2011_sou_v <- Br2011_sou_v[!(row.names(Br2011_sou_v) %in% rows_to_remove_Br2011), ]
        Br2011_sou_v$ts <- suppressWarnings(as.numeric(Br2011_sou_v$ts))



        # Reshape data for ggplot2
        DevSeq_AL_sou_v_div <- as.data.frame(c(DS_AL_sou_v[1:6, 1], DS_AL_sou_v[1:6, 2], 
          DS_AL_sou_v[1:6, 3], DS_AL_sou_v[1:6, 4], DS_AL_sou_v[1:6, 5], 
          DS_AL_sou_v[1:6, 6], DS_AL_sou_v[1:6, 7], DS_AL_sou_v[1:6, 8]))
        colnames(DevSeq_AL_sou_v_div) <- "correlation"

        Brawand2011_sou_v_div <- as.data.frame(c(Br2011_sou_v[1:6, 1], Br2011_sou_v[1:6, 2], 
            Br2011_sou_v[1:6, 3], Br2011_sou_v[1:6, 4], Br2011_sou_v[1:6, 5], Br2011_sou_v[1:6, 6]))
        colnames(Brawand2011_sou_v_div) <- "correlation"

    }




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
      DevSeq_AL_div_rates_wo_pollen <- DevSeq_AL_div[1:48,]
      DevSeq_AL_div_rates_pollen <- DevSeq_AL_div[49:54,]
      DevSeq_AL_div_rates_wo_pollen$comp_organ <- factor(DevSeq_AL_div_rates_wo_pollen$comp_organ, 
        levels = unique(DevSeq_AL_div_rates_wo_pollen$comp_organ))
      DevSeq_AL_div_rates_pollen$comp_organ <- factor(DevSeq_AL_div_rates_pollen$comp_organ, 
        levels = unique(DevSeq_AL_div_rates_pollen$comp_organ))

      

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
            pan_boarder <- 1.75
            axis_txt_size <- 21.25
            axis_ticks_s <- 0.7
            title_face <- "plain"

        } else {

            fname <- sprintf('%s.jpg', paste("GE_divergence_rates_AL", coefficient, expr_estimation, pos, sep="_"))
            plot_wdt <- 9.5 # condenced plot width for suppl
            plot_hdt <- 6.75 # condenced plot width for suppl
            legend_x_pos <- 0.638
            legend_y_pos <- 0.862
            legend_key_s <- 0.9
            linewd <- 2.7
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
            pan_boarder <- 1.8
            axis_txt_size <- 21.75
            axis_ticks_s <- 0.95
            title_face <- "bold"
        }

        

        p <- ggplot(data=data1, aes(x=div_times, y=correlation, group=comp_organ, colour=comp_organ)) + 
        geom_ribbon(aes(ymin = data1$lower, ymax = data1$upper, fill= comp_organ), alpha = 0.25, 
            linetype = 0, show.legend = FALSE) + 
        scale_fill_manual(values = c("Hypocotyl  "="#53b0db", "Stamen  "="#ee412e", "Flower  "="#e075af", 
                "Root  "="#6a54a9", "Apex veg  "="#96ba37", "Apex inf  "="#fad819", "Carpel  "="#f2a72f", 
                "Leaf  "="#2c8654")) + 
        geom_line(size = linewd) +  
        scale_x_continuous(limits = c(5.5,161.5), expand = c(0.02,0), breaks = c(7,9,25,46,106,160)) + 
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
        geom_segment(x=7, xend=7, y=0.435, yend=0.468, color="black", size=axis_ticks_s) + 
        geom_segment(x=9, xend=9, y=0.435, yend=0.468, color="black", size=axis_ticks_s) + 
        geom_segment(x=25, xend=25, y=0.435, yend=0.468, color="black", size=axis_ticks_s) + 
        geom_segment(x=46, xend=46, y=0.435, yend=0.468, color="black", size=axis_ticks_s) + 
        geom_segment(x=106, xend=106, y=0.435, yend=0.468, color="black", size=axis_ticks_s) + 
        geom_segment(x=160, xend=160, y=0.435, yend=0.468, color="black", size=axis_ticks_s) + 
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

      makeGEDivPlot(data1 = DevSeq_AL_div_rates_wo_pollen, data2 = DevSeq_AL_div_rates_pollen, 
        coefficient = coefficient, pos = "main")
      makeGEDivPlot(data1 = DevSeq_AL_div_rates_wo_pollen, data2 = DevSeq_AL_div_rates_pollen, 
        coefficient = coefficient, pos = "ext")




#------- Get pearson dist and reshape pearson dist and sOU data for regression analysis --------


      # Use pearson correlation, inter-organ normalization and TPM for ms

      getDSOrganCor <- function(df, organ, coefficient) {

        # log-transform data if TPM and Pearson are chosen
        if ((coefficient == "pearson") && (expr_estimation == "TPM")) {
            df <- log2(df + 1)
        }

        df_cor <- sqrt(1 - cor(df, method=coefficient))
        df_cor <- df_cor[4:nrow(df_cor), 1:3]

        # Reshape cor data frame to one column
        df_cor_rs <- data.frame(newcol = c(t(df_cor)), stringsAsFactors=FALSE)

        sp1 <- mean(df_cor_rs[1:9,])
        sp2 <- mean(df_cor_rs[10:18,])
        sp3 <- mean(df_cor_rs[19:27,])
        sp4 <- mean(df_cor_rs[28:36,])
        sp5 <- mean(df_cor_rs[37:45,])
        sp6 <- mean(df_cor_rs[46:54,])

        df_cor_avg <- rbind(sp1, sp2, sp3, sp4, sp5, sp6)
        colnames(df_cor_avg) <- organ

        getRowNames = function(x,n){ substring(x,nchar(x)-n+1) }
        row_names_repl <- getRowNames(rownames(df_cor),2)
        rnames_div_rates <- unique(row_names_repl)
        rownames(df_cor_avg) <- rnames_div_rates

        return(df_cor_avg)

      }

      root_div <- getDSOrganCor(df=x[,2:22], organ="Root", coefficient=coefficient)
      hypocotyl_div <- getDSOrganCor(df=x[,23:43], organ="Hypocotyl", coefficient=coefficient)
      leaf_div <- getDSOrganCor(df=x[,44:64], organ="Leaf", coefficient=coefficient)
      veg_apex_div <- getDSOrganCor(df=x[,65:85], organ="Apex veg", coefficient=coefficient)
      inf_apex_div <- getDSOrganCor(df=x[,86:106], organ="Apex inf", coefficient=coefficient)
      flower_div <- getDSOrganCor(df=x[,107:127], organ="Flower", coefficient=coefficient)
      stamen_div <- getDSOrganCor(df=x[,128:148], organ="Stamen", coefficient=coefficient)
      carpel_div <- getDSOrganCor(df=x[,149:169], organ="Carpel", coefficient=coefficient)

      DevSeq_AL_organ_cor <- cbind(root_div, hypocotyl_div, leaf_div, veg_apex_div, inf_apex_div, 
        flower_div, stamen_div, carpel_div)


      # Reshape data table for ggplot
      # divergence times are estimated taxon pair times from TimeTree
      # http://www.timetree.org/
      div_times <- rep(c(7.1, 9.4, 25.6, 46, 106, 160), times=8)
      comp_organ <- rep(colnames(DevSeq_AL_organ_cor), each=6)
      comp_spec <- rep(rownames(DevSeq_AL_organ_cor), times=8)
      dataset <- rep("Angiosperms ", 48)

      DevSeq_AL_GE_dist <- rbind(root_div, hypocotyl_div, leaf_div, veg_apex_div, inf_apex_div, 
        flower_div, stamen_div, carpel_div)
      rownames(DevSeq_AL_GE_dist) <- NULL
      colnames(DevSeq_AL_GE_dist) <- "correlation"

      DevSeq_AL_div_rates <- data.frame(cbind(comp_spec, comp_organ, div_times, DevSeq_AL_GE_dist, dataset), 
        stringsAsFactors=FALSE)

      DevSeq_AL_div_rates$div_times <- as.numeric(DevSeq_AL_div_rates$div_times)
      DevSeq_AL_div_rates$correlation <- as.numeric(DevSeq_AL_div_rates$correlation)

      DevSeq_AL_div_rates$comp_organ <- factor(DevSeq_AL_div_rates$comp_organ, 
        levels = unique(DevSeq_AL_div_rates$comp_organ))


      # Reshape DevSeq AL sOU expression data
      DevSeq_AL_sou_v_div_rates <- data.frame(cbind(comp_spec, comp_organ, div_times, DevSeq_AL_sou_v_div, dataset), 
        stringsAsFactors=FALSE)

      DevSeq_AL_sou_v_div_rates$div_times <- as.numeric(DevSeq_AL_sou_v_div_rates$div_times)
      DevSeq_AL_sou_v_div_rates$correlation <- as.numeric(DevSeq_AL_sou_v_div_rates$correlation)

      DevSeq_AL_sou_v_div_rates$comp_organ <- factor(DevSeq_AL_sou_v_div_rates$comp_organ, 
        levels = unique(DevSeq_AL_sou_v_div_rates$comp_organ))



      # Reshape Brawand sOU expression data
      div_times <- rep(c(6.7, 9.1, 15.8, 29.4, 90, 159), times=6)
      comp_organ <- rep(colnames(Br2011_sou_v), each=6)
      comp_spec <- rep(rownames(Br2011_sou_v[1:6,]), times=6)
      dataset <- rep("Mammals", 36)

      Brawand11_sou_v_div_rates <- data.frame(cbind(comp_spec, comp_organ, div_times, 
      Brawand2011_sou_v_div, dataset), stringsAsFactors=FALSE)

      Brawand11_sou_v_div_rates$div_times <- as.numeric(Brawand11_sou_v_div_rates$div_times)
      Brawand11_sou_v_div_rates$correlation <- as.numeric(Brawand11_sou_v_div_rates$correlation)

      # Remove ppy testis sample (has NA value)
      Brawand11_sou_v_div_rates <- Brawand11_sou_v_div_rates[c(-33),]

      Brawand11_sou_v_div_rates$comp_organ <- factor(Brawand11_sou_v_div_rates$comp_organ, 
        levels = unique(Brawand11_sou_v_div_rates$comp_organ))


      # Combine DevSeq and Brawand 2011 GE divergence data
      compSouVDivRates11_AL <- rbind(DevSeq_AL_sou_v_div_rates, Brawand11_sou_v_div_rates)




#------------------------------ Get slopes of DevSeq AL log model ------------------------------


      # Retrieve slope value from individual organ regressions and compute p value
      getLogReg <- function(corrdata) {

        model_var <- correlation ~ log(div_times) * comp_organ #linear regression with log-x-transform

        model_lmp = lm(model_var, data = corrdata)

        div_trend <- lstrends(model_lmp, "comp_organ", var = "div_times") # get slopes for regressions
        div_trend_df <- summary(div_trend)

        slopes <- as.data.frame(div_trend_df[,2])
        colnames(slopes) <- "log_reg"
        rownames(slopes) <- div_trend_df[,1]

        return(slopes)
      }


      # Get organ slopes and p-values
      p_values_compDivRates_AL_io <- getLogReg(DevSeq_AL_div_rates)
      colnames(p_values_compDivRates_AL_io) <- "DS_AL_pea_log"
      p_values_compSouVDivRates_AL_io <- getLogReg(DevSeq_AL_sou_v_div_rates)
      colnames(p_values_compSouVDivRates_AL_io) <- "DS_AL_sOU_log"

      sample <- rownames(p_values_compDivRates_AL_io)
      sample <- gsub(" ", "_", sample)
      DevSeq_AL_log_slopes <- cbind(sample, p_values_compSouVDivRates_AL_io, p_values_compDivRates_AL_io)
      rownames(DevSeq_AL_log_slopes) <- NULL




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

      devseqSouV_organ_lst_sel <- list(compSouVDivRates11_AL[1:6,], compSouVDivRates11_AL[7:11,], compSouVDivRates11_AL[13:18,], 
        compSouVDivRates11_AL[19:24,], compSouVDivRates11_AL[25:30,], compSouVDivRates11_AL[31:36,], compSouVDivRates11_AL[37:42,], 
        compSouVDivRates11_AL[43:48,])

      brawandSouV11_organ_lst <- list(compSouVDivRates11_AL[49:54,], compSouVDivRates11_AL[55:60,], compSouVDivRates11_AL[61:66,], 
        compSouVDivRates11_AL[67:72,],compSouVDivRates11_AL[73:78,], compSouVDivRates11_AL[79:83,])


      # Set up lists containing metric pearson expression distances
      devseq_organ_lst <- list(DevSeq_AL_div_rates[1:6,], DevSeq_AL_div_rates[7:12,], DevSeq_AL_div_rates[13:18,], 
        DevSeq_AL_div_rates[19:24,], DevSeq_AL_div_rates[25:30,], DevSeq_AL_div_rates[31:36,], DevSeq_AL_div_rates[37:42,], 
        DevSeq_AL_div_rates[43:48,])

      devseq_organ_lst_sel <- list(DevSeq_AL_div_rates[1:6,], DevSeq_AL_div_rates[7:11,], DevSeq_AL_div_rates[13:18,], 
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



      getLOESS.Slopes <- function(organ_data, data_set) {

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
      DevSeqSouV_AL_loess_slopes <- as.data.frame(do.call(rbind, lapply(devseqSouV_organ_lst, getLOESS.Slopes, data_set="complete")))
      DevSeqSouV_sel_AL_loess_slopes <- as.data.frame(do.call(rbind, lapply(devseqSouV_organ_lst_sel, getLOESS.Slopes, data_set="selected"))) ## hypocotyl slope is 0.00779 instead 0.00758 if BD is left out
      Brawand11SouV_loess_slopes <- as.data.frame(do.call(rbind, lapply(brawandSouV11_organ_lst, getLOESS.Slopes, data_set="complete")))
      sOU_loess_DevSeq_AL_Br11_wilcox <- wilcox.test(as.numeric(unlist(DevSeqSouV_AL_loess_slopes)), as.numeric(unlist(Brawand11SouV_loess_slopes)))$p.value

      # For metric pearson expression distances
      DevSeq_AL_loess_slopes <- as.data.frame(do.call(rbind, lapply(devseq_organ_lst, getLOESS.Slopes, data_set="complete")))


      # Write slope values to csv file
      DevSeq_AL_loess_slopes <- cbind(DevSeqSouV_AL_loess_slopes, DevSeq_AL_loess_slopes)
      colnames(DevSeq_AL_loess_slopes) <- c("DS_AL_sOU_loess", "DS_AL_pea_loess")
      ds_organs <- data.frame(sample = c("Root", "Hypocotyl", "Leaf", "Apex_veg", "Apex_inf", "Flower", "Stamen", "Carpel"))
      DevSeq_AL_loess_slopes <- cbind(ds_organs, DevSeq_AL_loess_slopes)



      # Show message
      message("Writing data tables...")

      write.table(DevSeq_AL_loess_slopes, 
        file=file.path(out_dir, "output", "data", "DevSeq_AL_loess_slopes.txt"), sep="\t", 
        col.names=TRUE, row.names=FALSE, dec=".", quote = FALSE)

      write.table(DevSeq_AL_log_slopes, 
        file=file.path(out_dir, "output", "data", "DevSeq_AL_log_slopes.txt"), sep="\t", 
        col.names=TRUE, row.names=FALSE, dec=".", quote = FALSE)


      # Create p-value containing test strings for plots
    sOU_loess_DevSeq_AL_Br11_slope_p <- paste("P =", formatC(sOU_loess_DevSeq_AL_Br11_wilcox, format="e", digits=0))



   # Make sOU GE divergence plot showing individual organ regressions for SI
   makeOrgRegPlot <- function(data1, data2, coefficient, expr_estimation, p_value, pos) {

      fname <- sprintf('%s.jpg', paste("compSouVDivRates11_AL_loess", expr_estimation, pos, sep="_"))

      ymin <- 0.05
      ymax <- 1.45

      if (pos == "main") {

            plot_wdt <- 12.535
            plot_hdt <- 8
            legend_x_pos <- 0.241
            legend_y_pos <- 0.914
            linewd <- 3
            point_size <- 5.75
            txt_x_pos <- 16.25
            txt_y_pos <- 1.283
            pan_boarder <- 1.75
            axis_txt_size <- 21.25
            axis_ticks_s <- 0.7
            title_face <- "plain"

      } else {

            plot_wdt <- 9.5 # condenced plot width for suppl
            plot_hdt <- 6.75 # condenced plot width for suppl
            legend_x_pos <- 0.317
            legend_y_pos <- 0.9
            linewd <- 2.5
            point_size <- 5
            txt_x_pos <- 16.25
            txt_y_pos <- 1.283
            pan_boarder <- 1.8
            axis_txt_size <- 21.75
            axis_ticks_s <- 0.95
            title_face <- "bold"
      }

      ds_col <- rep(c("#798dc4"), 48)
      bw_col <- rep(c("red"), 35)
      col_scale <- c("#798dc4", "red")
      fill_scale <- c("#798dc4", "red")
      col_breaks <- c("Angiosperms ", "Mammals")
      fill_breaks <- c("Angiosperms ", "Mammals")
      shape_scale <- c(16, 17)

      fill_col <- c(as.character(ds_col), as.character(bw_col))

      p <- ggplot() 
      p <- p + geom_line(size = linewd, data = data1, aes(x = div_times, y = correlation, group = comp_organ, 
        colour = dataset, linetype = dataset))
      p <- p + geom_point(size = point_size, data = data2, aes(x = div_times, y = correlation, group = comp_organ, 
        colour = dataset, shape = dataset))
      p <- p + scale_x_continuous(limits = c(0,162), expand = c(0.02,0), breaks = c(0,20,40,60,80,100,120,140,160)) + 
      scale_y_continuous(limits = c(0.055, 1.5275), expand = c(0.02, 0), breaks = c(0.2,0.4,0.6,0.8,1,1.2,1.4)) + 
      scale_color_manual(values = col_scale, breaks = col_breaks) + 
      scale_fill_manual(values = fill_scale, breaks = fill_breaks) + 
      scale_shape_manual(values = shape_scale) + 
      scale_size(range = c(0.5, 12)) + 
      annotate("text", x = txt_x_pos, y = txt_y_pos, label = p_value, size = 8) + 
      guides(color = guide_legend(ncol=2, keywidth = 0.4, keyheight = 0.4, default.unit = "inch"))

      q <- p + theme_bw() + xlab("Divergence time (Myr)") + ylab("Expression distance") + 
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
        legend.box.background = element_rect(colour = "#d5d5d5", fill= "white" , size=1.0), 
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



   }
   
}


makeCompAnylsisAL(expr_estimation="TPM", coefficient="pearson")





