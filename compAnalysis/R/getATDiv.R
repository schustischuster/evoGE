# Prepare Brawand and DevSeq comparative expression data
# Thresholds: 0.5 TPM (since there are no ERCC spike-ins in Brawand data)
# Data input: Brawand and DevSeq TPM expression tables of all samples



#-------------------------------------- Read data tables ---------------------------------------


getATDiv <- function(coefficient = c("pearson", "spearman"), expr_estimation = c("TPM", "counts")) {
	

    # Show error message if no correlation or unknown correlation coefficient is chosen
    if ((missing(coefficient)) | (!is.element(coefficient, c("pearson", "spearman"))))
   
       stop(
       "Please choose one of the available correlation coefficients: 
	   'pearson', 'spearman'",
	   call. = TRUE
       )

    # Show error message if expression estimation or unknown expression estimation is chosen
    if ((missing(expr_estimation)) | (!is.element(expr_estimation, c("TPM", "counts"))))
   
       stop(
       "Please choose one of the available expression estimations: 
       'TPM', 'counts'",
       call. = TRUE
       )


    # Show startup message
    message("Reading data...")


    # Set file path for input files
    if (is.element(expr_estimation, "TPM")) {
        
        genesExprDS = file.path(in_dir, "Expression_data", "inter_organ_tpm_mat_deseq_sample_names_all.csv")
        genesExprBr = file.path(in_dir, "Expression_data", "TPM_Brawand_norm_inter_organ.csv")

    } else if (is.element(expr_estimation, "counts")) {
        
        genesExprDS = file.path(in_dir, "Expression_data", "inter_organ_count_mat_vsd_sample_names_all.csv")
        genesExprBr = file.path(in_dir, "Expression_data", "count_Brawand_norm_inter_organ.csv")

    }

    
    # Set colnames
    col_namesDS <- rep(c("Root", "Hypocotyl", "Leaf", "veg_apex", "inf_apex", 
        "Flower", "Stamen", "Carpel"), each=21)
    replicate_tag_samples_DS <- rep(c(".1",".2",".3"), times=8)
    col_namesDS <- paste0(col_namesDS, replicate_tag_samples_DS)
    spec_namesDS <- rep(c("_AT", "_AL", "_CR", "_ES", "_TH", "_MT", "_BD"), each=3)
    spec_namesDS <- rep(spec_namesDS, times=8)
    col_namesDS <- paste0(col_namesDS, spec_namesDS)


    # Read DevSeq table
	x_DS <- read.table(genesExprDS, sep=";", dec=".", skip = 1, header=FALSE, stringsAsFactors=FALSE)
    
    # Remove later on once expression table has gene_id column
    ID_repl <- as.data.frame(seq(1:nrow(x_DS)))
    colnames(ID_repl) <- "gene_id"
    x_DS <- cbind(ID_repl, x_DS)

    # set column names
    colnames(x_DS)[2:ncol(x_DS)] <- col_namesDS


	# Read Brawand table and set colnames
	x_Br <- read.table(genesExprBr, sep=";", dec=".", header=TRUE, stringsAsFactors=FALSE)

    # Remove later on once expression table has gene_id column
    ID_repl <- as.data.frame(seq(1:nrow(x_Br)))
    colnames(ID_repl) <- "gene_id"
    x_Br <- cbind(ID_repl, x_Br)

	# colnames(x_Br)[1] <- "gene_id"

    # Remove biological replicates that show log2 TPM Pearson's r < 0.85
    x_Br <- x_Br %>% select (-c(Human_brain__.1., Human_brain__.2..1, Human_heart__.1., Human_kidney__.2.))
	


    # Stop function here to allow specific analysis of a single data set
    # return_list <- list("x_DS" = x_DS, "x_Br" = x_Br, "coefficient" = coefficient, "expr_estimation" = expr_estimation)
    # return(return_list)
    # }
    # return_objects <- getATDiv(coefficient="pearson", expr_estimation=TPM) # read in Comparative expression data
    # list2env(return_objects, envir = .GlobalEnv)

    


#--------------------- Calculate correlation and prepare data for corrplot ---------------------


    # Create "plots" folder in /out_dir/output/plots
    if (!dir.exists(file.path(out_dir, "output", "plots"))) 
        dir.create(file.path(out_dir, "output", "plots"), recursive = TRUE)

    # Show message
    message("Starting analysis and generate plots...")



    x_Br[is.na(x_Br)] <- 0 # replaces NAs by 0
    x_DS[is.na(x_DS)] <- 0 # replaces NAs by 0


    if (is.element("pearson", coefficient) && is.element("TPM", expr_estimation)) {
            
        x_Br[,2:ncol(x_Br)] <- log2(x_Br[,2:ncol(x_Br)] + 1)
        x_DS[,2:ncol(x_DS)] <- log2(x_DS[,2:ncol(x_DS)] + 1)

    }




#------------------- Apply 0.5 TPM threshold to both DevSeq and Brawand data -------------------


# Need to implement expression data thresholding here!!!
# Use 0.5 TPM in at least 2 of 3 replicates; check condition for Brawand
# (number of biological replicates in Brawand data set is different in various organs)




#---------------- Get gene expression divergence rates for ATH/AL vs species X -----------------


    # Use pearson correlation, inter-organ normalization and TPM for ms

    getDSOrganCor <- function(df, organ, coefficient) {

        df_cor <- cor(df, method=coefficient)
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

    root_div <- getDSOrganCor(df=x_DS[,2:22], organ="Root", coefficient=coefficient)
    hypocotyl_div <- getDSOrganCor(df=x_DS[,23:43], organ="Hypocotyl", coefficient=coefficient)
    leaf_div <- getDSOrganCor(df=x_DS[,44:64], organ="Leaf", coefficient=coefficient)
    veg_apex_div <- getDSOrganCor(df=x_DS[,65:85], organ="Apex veg", coefficient=coefficient)
    inf_apex_div <- getDSOrganCor(df=x_DS[,86:106], organ="Apex inf", coefficient=coefficient)
    flower_div <- getDSOrganCor(df=x_DS[,107:127], organ="Flower", coefficient=coefficient)
    stamen_div <- getDSOrganCor(df=x_DS[,128:148], organ="Stamen", coefficient=coefficient)
    carpel_div <- getDSOrganCor(df=x_DS[,149:169], organ="Carpel", coefficient=coefficient)

    DevSeq_organ_cor <- cbind(root_div, hypocotyl_div, leaf_div, veg_apex_div, inf_apex_div, 
        flower_div, stamen_div, carpel_div)


    # Reshape data table for ggplot
    # divergence times are estimated taxon pair times from TimeTree
    # http://www.timetree.org/
    div_times <- rep(c(7.1, 9.4, 25.6, 46, 106, 160), times=8)
    comp_organ <- rep(colnames(DevSeq_organ_cor), each=6)
    comp_spec <- rep(rownames(DevSeq_organ_cor), times=8)
    dataset <- rep("Angiosperms", 48)

    DevSeq_GE_div <- rbind(root_div, hypocotyl_div, leaf_div, veg_apex_div, inf_apex_div, 
        flower_div, stamen_div, carpel_div)
    rownames(DevSeq_GE_div) <- NULL
    colnames(DevSeq_GE_div) <- "correlation"

    DevSeq_div_rates <- data.frame(cbind(comp_spec, comp_organ, div_times, DevSeq_GE_div, dataset), 
        stringsAsFactors=FALSE)

    DevSeq_div_rates$div_times <- as.numeric(DevSeq_div_rates$div_times)
    DevSeq_div_rates$correlation <- as.numeric(DevSeq_div_rates$correlation)
      
    # Remove Brachypodium mesocotyl data point
    DevSeq_div_rates <- DevSeq_div_rates[-12,]

    DevSeq_div_rates$comp_organ <- factor(DevSeq_div_rates$comp_organ, 
        levels = unique(DevSeq_div_rates$comp_organ))




#---------------- Get gene expression divergence rates for Human vs species X ------------------


    # Use pearson correlation, inter-organ normalization and TPM for ms

    getBrBrainCor <- function(df, organ, coefficient) {

        df_cor <- cor(df, method=coefficient)
        df_cor <- df_cor[5:nrow(df_cor), 1:4]

        # Reshape cor data frame to one column
        df_cor_rs <- data.frame(newcol = c(t(df_cor)), stringsAsFactors=FALSE)

        Ppa <- mean(df_cor_rs[1:12,]) # bonobo
        Pan <- mean(df_cor_rs[13:36,]) # chimp
        Ggo <- mean(df_cor_rs[37:44,]) # gorilla
        Ppy <- mean(df_cor_rs[45:52,]) # orangutan
        Mml <- mean(df_cor_rs[53:64,]) # macaque
        Mmu <- mean(df_cor_rs[65:76,]) # mouse
        Mdo <- mean(df_cor_rs[77:84,]) # opossum

        df_cor_avg <- rbind(Ppa, Ggo, Ppy, Mml, Mmu, Mdo)
        colnames(df_cor_avg) <- organ

        return(df_cor_avg)

    }

    brain_div <- getBrBrainCor(df=x_Br[,2:26], organ="Brain", coefficient=coefficient)


    getBrCerebCor <- function(df, organ, coefficient) {

        df_cor <- cor(df, method=coefficient)
        df_cor <- df_cor[3:nrow(df_cor), 1:2]

        # Reshape cor data frame to one column
        df_cor_rs <- data.frame(newcol = c(t(df_cor)), stringsAsFactors=FALSE)

        Ppa <- mean(df_cor_rs[1:4,]) # bonobo
        Pan <- mean(df_cor_rs[5:8,]) # chimp
        Ggo <- mean(df_cor_rs[9:12,]) # gorilla
        Ppy <- mean(df_cor_rs[13:14,]) # orangutan
        Mml <- mean(df_cor_rs[15:18,]) # macaque
        Mmu <- mean(df_cor_rs[19:24,]) # mouse
        Mdo <- mean(df_cor_rs[25:28,]) # opossum

        df_cor_avg <- rbind(Ppa, Ggo, Ppy, Mml, Mmu, Mdo)
        colnames(df_cor_avg) <- organ

        return(df_cor_avg)

    }

    cereb_div <- getBrCerebCor(df=x_Br[,27:42], organ="Cerebellum", coefficient=coefficient)


    getBrHtKdLvCor <- function(df, organ, coefficient) {

        df_cor <- cor(df, method=coefficient)
        df_cor <- df_cor[3:nrow(df_cor), 1:2]

        # Reshape cor data frame to one column
        df_cor_rs <- data.frame(newcol = c(t(df_cor)), stringsAsFactors=FALSE)

        Ppa <- mean(df_cor_rs[1:4,]) # bonobo
        Pan <- mean(df_cor_rs[5:8,]) # chimp
        Ggo <- mean(df_cor_rs[9:12,]) # gorilla
        Ppy <- mean(df_cor_rs[13:16,]) # orangutan
        Mml <- mean(df_cor_rs[17:20,]) # macaque
        Mmu <- mean(df_cor_rs[21:26,]) # mouse
        Mdo <- mean(df_cor_rs[27:30,]) # opossum

        df_cor_avg <- rbind(Ppa, Ggo, Ppy, Mml, Mmu, Mdo)
        colnames(df_cor_avg) <- organ

        return(df_cor_avg)

    }

    heart_div <- getBrHtKdLvCor(df=x_Br[,43:59], organ="Heart", coefficient=coefficient)
    kidney_div <- getBrHtKdLvCor(df=x_Br[,60:76], organ="Kidney", coefficient=coefficient)
    liver_div <- getBrHtKdLvCor(df=x_Br[,77:93], organ="Liver", coefficient=coefficient)


    getBrTestisCor <- function(df, organ, coefficient) {

        df_cor <- cor(df, method=coefficient)
        df_cor <- df_cor[3:nrow(df_cor), 1:2]

        # Reshape cor data frame to one column
        df_cor_rs <- data.frame(newcol = c(t(df_cor)), stringsAsFactors=FALSE)

        Ppa <- mean(df_cor_rs[1:2,]) # bonobo
        Pan <- mean(df_cor_rs[3:4,]) # chimp
        Ggo <- mean(df_cor_rs[5:6,]) # gorilla
        Ppy <- NA # orangutan: no data available
        Mml <- mean(df_cor_rs[7:10,]) # macaque
        Mmu <- mean(df_cor_rs[11:14,]) # mouse
        Mdo <- mean(df_cor_rs[15:18,]) # opossum

        df_cor_avg <- rbind(Ppa, Ggo, Ppy, Mml, Mmu, Mdo)
        colnames(df_cor_avg) <- organ

        return(df_cor_avg)

    }

    testis_div <- getBrTestisCor(df=x_Br[,94:104], organ="Testis", coefficient=coefficient)

    Brawand_organ_cor <- cbind(brain_div, cereb_div, heart_div, kidney_div, liver_div, testis_div)


    # Reshape data table for ggplot
    # divergence times are estimated taxon pair times from TimeTree
    # http://www.timetree.org/
    div_times <- rep(c(6.7, 9.1, 15.8, 29.4, 90, 159), times=6)
    comp_organ <- rep(colnames(Brawand_organ_cor), each=6)
    comp_spec <- rep(rownames(Brawand_organ_cor), times=6)
    dataset <- rep("Mammals", 36)

    Brawand_GE_div <- rbind(brain_div, cereb_div, heart_div, kidney_div, liver_div, testis_div)
    rownames(Brawand_GE_div) <- NULL
    colnames(Brawand_GE_div) <- "correlation"

    Brawand_div_rates <- data.frame(cbind(comp_spec, comp_organ, div_times, Brawand_GE_div, dataset), 
        stringsAsFactors=FALSE)

    Brawand_div_rates$div_times <- as.numeric(Brawand_div_rates$div_times)
    Brawand_div_rates$correlation <- as.numeric(Brawand_div_rates$correlation)

    # Remove Orangutan testis (missing data) and Opossum kidney (replicate corr < 0.85) data
    Brawand_div_rates <- Brawand_div_rates[c(-24,-33),]

    Brawand_div_rates$comp_organ <- factor(Brawand_div_rates$comp_organ, 
        levels = unique(Brawand_div_rates$comp_organ))


    # Combine DevSeq and Brawand GE divergence data
    compDivRates <- rbind(DevSeq_div_rates, Brawand_div_rates)




#--------- Make gene expression divergence rates plot of mammalian and angiosperm data ---------
      

      # Make GE divergence plot
      makeGEDivPlot <- function(data1, coefficient, expr_estimation) {

        fname <- sprintf('%s.jpg', paste("comp_divergence_rates", coefficient, expr_estimation, sep="_"))

        p <- ggplot(data=data1, aes(x = div_times, y = correlation, group = dataset, colour = dataset)) + 
        geom_point(size = 4) + 
        geom_smooth(method ='lm', size = 2) + 
        scale_x_continuous(limits = c(0,160), expand = c(0.02,0), breaks = c(0,20,40,60,80,100,120,140,160)) + 
        scale_y_continuous(limits = c(0.565, 0.91), expand = c(0.02, 0)) + 
        scale_color_manual(values = c("#8591c7", "red"), 
            # organ order: hypocotyl/stamen/flower/root/veg_apex/inf_apex/carpel/leaf
            breaks=c("Angiosperms", "Mammals")) + 
        guides(color = guide_legend(ncol = 3))

        q <- p + theme_bw() + xlab("Divergence time (Myr)") + ylab("Pearson's r") + 
        theme(text=element_text(size=16), 
            axis.ticks.length=unit(0.35, "cm"), 
            axis.ticks = element_line(colour = "black", size = 0.7),  
            plot.margin = unit(c(0.55, 1.1, 0.5, 0.4),"cm"), 
            axis.title.y = element_text(size=25, margin = margin(t = 0, r = 17, b = 0, l = 9)), 
            axis.title.x = element_text(size=25, margin = margin(t = 14.75, r = 0, b = 2, l = 0)), 
            axis.text.x = element_text(size=21.25, angle=0, margin = margin(t = 5.5)), 
            axis.text.y = element_text(size=21.25, angle=0, margin = margin(r = 5.5)), 
            legend.box.background = element_rect(colour = "#d5d5d5", fill=NA, size=1.0), 
            panel.border = element_rect(colour = "black", fill=NA, size=1.75), 
            panel.grid.major = element_line(color="#d5d5d5"),
            panel.grid.minor.x = element_blank(), 
            panel.grid.minor.y = element_blank(), 
            legend.position = c(0.723, 0.92), 
            legend.title = element_blank(), 
            legend.text = element_text(size=21.5), 
            legend.spacing.x = unit(0.5, 'cm'), 
            legend.key.size = unit(0.95, "cm"), 
            legend.background=element_blank()) 

        ggsave(file = file.path(out_dir, "output", "plots", fname), plot = q, 
            width = 12.535, height = 8, dpi = 300, units = c("in"), limitsize = FALSE) 
      }

      makeGEDivPlot(data1 = compDivRates, coefficient = coefficient, expr_estimation = expr_estimation)

}


getATDiv(coefficient = "pearson", expr_estimation = "TPM")
getATDiv(coefficient = "spearman", expr_estimation = "TPM")

getATDiv(coefficient = "pearson", expr_estimation = "counts")
getATDiv(coefficient = "spearman", expr_estimation = "counts")



