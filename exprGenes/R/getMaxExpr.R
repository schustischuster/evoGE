# Count number of genes with maximum expression in each organ in Arabidopsis thaliana for all organs,
# and additionally for only the comparative organs; do the same for all the other species
# The comparative organs are: root, hypocotyl, leaf, apex veg, apex inf, flower, stamen, carpel
# Data input: Normalized expr expression data containing protein-coding genes, NATs and lincRNAs
# Classification: (1) All protein-coding genes, (2) 7003 core orthologous protein-coding genes, 
# (3) all NATs, (4) all lincRNAs, (5) orthologous lncRNAs (Brassicaceae)
# Input sample tables should have the following format:
# DEVSEQ_SAMPLE_REPLICATES(between 27 and 132 depending on species), rownames = gene_id


#------------------- Load packages, set directories and read sample tables ---------------------


# Define function to get expressed genes

getMaxExpr <- function(species = c("AT", "all"), ...) {
   
   # Show error message if no species is chosen
   if (missing(species))

   stop(
      "Please choose one of the available species: 
      'AT', 'all'",
      call. = TRUE
      )

   species_id <- species


   # Show startup message
   message("Reading data...")

   # Set file path to expression data
   pathAT = file.path(in_dir, "Expression_data", "AT_genes_inter_norm_count_mat_vsd_sample_names.csv")
   pathAL = file.path(in_dir, "Expression_data", "AL_genes_inter_norm_count_mat_vsd_sample_names.csv")
   pathCR = file.path(in_dir, "Expression_data", "CR_genes_inter_norm_count_mat_vsd_sample_names.csv")
   pathES = file.path(in_dir, "Expression_data", "ES_genes_inter_norm_count_mat_vsd_sample_names.csv")
   pathTH = file.path(in_dir, "Expression_data", "TH_genes_inter_norm_count_mat_vsd_sample_names.csv")
   pathMT = file.path(in_dir, "Expression_data", "MT_genes_inter_norm_count_mat_vsd_sample_names.csv")
   pathBD = file.path(in_dir, "Expression_data", "BD_genes_inter_norm_count_mat_vsd_sample_names.csv")

   # Ortholog tables
   pathCore = file.path(in_dir, "Expression_data", "AT_core_inter_tpm_mat_deseq_sample_names.csv")
   pathPcBrass = file.path(in_dir, "Expression_data", "AT_brass_inter_tpm_mat_deseq_sample_names.csv")
   pathNcBrass = file.path(in_dir, "Expression_data", "lnc_AT_brass_inter_tpm_mat_deseq_sample_names.csv")

   # Tables containing raw expression values and biotype annotation
   pathAT_compl <- file.path(in_dir, "Expression_data", "AT_genes_complete_table_tpm_sample_names.csv")
   pathAL_compl <- file.path(in_dir, "Expression_data", "AL_genes_complete_table_tpm_sample_names.csv")
   pathCR_compl <- file.path(in_dir, "Expression_data", "CR_genes_complete_table_tpm_sample_names.csv")
   pathES_compl <- file.path(in_dir, "Expression_data", "ES_genes_complete_table_tpm_sample_names.csv")
   pathTH_compl <- file.path(in_dir, "Expression_data", "TH_genes_complete_table_tpm_sample_names.csv")
   pathMT_compl <- file.path(in_dir, "Expression_data", "MT_genes_complete_table_tpm_sample_names.csv")
   pathBD_compl <- file.path(in_dir, "Expression_data", "BD_genes_complete_table_tpm_sample_names.csv")


   # Set up list of expression tables
   if (species == "all") {

      expr_table_ls <- list(AT_expr = pathAT, AL_expr = pathAL, CR_expr = pathCR, ES_expr = pathES, 
         TH_expr = pathTH, MT_expr = pathMT, BD_expr = pathBD, Core_expr = pathCore, Brass_pc_expr = pathPcBrass, 
         Brass_nc_expr = pathNcBrass, AT_expr_compl = pathAT_compl, AL_expr_compl = pathAL_compl, 
         CR_expr_compl = pathCR_compl, ES_expr_compl = pathES_compl, TH_expr_compl = pathTH_compl, 
         MT_expr_compl = pathMT_compl, BD_expr_compl = pathBD_compl)

   } else if (species == "AT") {

      expr_table_ls <- list(AT_expr = pathAT, Core_expr = pathCore, Brass_pc_expr = pathPcBrass, 
         Brass_nc_expr = pathNcBrass, AT_expr_compl = pathAT_compl)
   }


   # Read expression data
   expr_tables <- lapply(expr_table_ls, function(x) {
      read.table(x, sep=";", dec=".", header=TRUE, stringsAsFactors=FALSE)
      })


   # Stop function here to run tests
   # return_list <- list("species_id" = species_id, "expr_tables" = expr_tables)
   # return(return_list)
   # }
   # return_objects <- getMaxExpr(species="AT") # read in expr expression data
   # list2env(return_objects, envir = .GlobalEnv)


   list2env(expr_tables, envir = .GlobalEnv)

   # Show message
   message("Starting analysis...")




   # ---------------------- Perform analysis for A.thaliana across all organs ----------------------


   if (species_id == "AT") {


      # Get replicate expression
      calculateAvgExpr <- function(df) {

         # Split data frame by sample replicates into a list and get rowMeans for each subset
         averaged_replicates <- do.call(cbind, lapply(split.default(df, 
            rep(seq_along(df), 
               each = 3, 
               length.out = ncol(df))
            ), rowMeans)
         )

         averaged_replicates <- as.data.frame(averaged_replicates)

         return(averaged_replicates)
      }

      AT_expr_repl <- calculateAvgExpr(AT_expr)


      # Get replicate names
      getReplNames <- function(n) {

         df_names <- names(n)
         df_names <- unique(substring(df_names, 1, nchar(df_names)-4))
         df_names <- gsub('[.]', '', df_names)

         return(df_names)
      }

      repl_names <- getReplNames(AT_expr)

      colnames(AT_expr_repl) <- repl_names


      # Extract protein-coding core ortholog and Brassicaceae lncRNA ortholog gene IDs
      core_ids <- sub("\\:.*", "", Core_expr[,1])
      core_ids <- core_ids[!grepl("ERCC", core_ids)]


      # Extract all A.thaliana lncRNAs and get orthologs
      AT_lncRNA_ids <- subset(AT_expr_compl, subset = biotype %in% c("lnc_exonic_antisense", 
         "lnc_intronic_antisense", "lnc_intergenic"))[,1]
      core_brass_ids <- sub("\\:.*", "", Brass_pc_expr[,1])
      core_lnc_ids <- sub("\\:.*", "", Brass_nc_expr[,1])


      # Get number of genes with maximum expression for each organ
      # "Merge" dev stages for each organ by calculating average number of genes w/ max expression

      getNumGenes <- function(df, scripttype = c("coding", "lncRNA"), c_level = c("all", "brass", "core")) {

         if ((c_level == "all") && (scripttype == "coding")) {

            df_at1 <- df[rownames(df) %like% "AT1G", ]
            df_at2 <- df[rownames(df) %like% "AT2G", ]
            df_at3 <- df[rownames(df) %like% "AT3G", ]
            df_at4 <- df[rownames(df) %like% "AT4G", ]
            df_at5 <- df[rownames(df) %like% "AT5G", ]

            df <- rbind(df_at1, df_at2, df_at3, df_at4, df_at5)

         } else if ((c_level == "all") && (scripttype == "lncRNA")) {

            df <- df[rownames(df) %in% AT_lncRNA_ids, ]
         
         # Reduce data to core orthologs if core set is chosen

         } else if ((c_level == "core") && (scripttype == "coding")) {

            df <- df[rownames(df) %in% core_ids,]

         } else if ((c_level == "brass") && (scripttype == "coding")) {

            df <- df[rownames(df) %in% core_brass_ids,]
         
         } else if ((c_level == "core") && (scripttype == "lncRNA")) {

            df <- df[rownames(df) %in% core_lnc_ids, ]
         }


         # Create group names
         groups <- c(rep("Root", 6), "Hypocotyl" ,rep("Stem", 3), rep("Leaf", 9), rep("Apex", 6), 
            rep("Flower", 4), rep("Sepals", 2), rep("Petals", 2), rep("Stamen", 2), rep("Carpel", 2), 
            rep("Fruit", 3), rep("Seed", 3))

         df$max <- colnames(df)[apply(df, 1, which.max)]

         x_max <- df$max
         sample_count <- do.call(rbind, lapply(colnames(df)[1:(length(df)-1)], function(x) {
            length(grep(x, x_max))}))

         # Calculate average maximum expression per organ group
         # This is the Weighted Arithmetic Mean (mean of organ weighted by number of stages)
         avg_n_max <- c(rep(mean(sample_count[1:6]),6), sample_count[7], rep(mean(sample_count[8:10]),3), 
            rep(mean(sample_count[11:19]),9), rep(mean(sample_count[20:25]),6), rep(mean(sample_count[26:29]),4), 
            rep(mean(sample_count[30:31]),2), rep(mean(sample_count[32:33]),2), rep(mean(sample_count[34:35]),2), 
            rep(mean(sample_count[36:37]),2), rep(mean(sample_count[38:40]),3), rep(mean(sample_count[41:43]),3))

         df_out <- data.frame(
            biotype = rep(scripttype), 
            conservation = rep(c_level),
            organ = colnames(df)[1:(length(df)-1)], 
            group = groups, 
            count = sample_count, 
            average_count = avg_n_max
            )

         total_avg <- sum(unique(avg_n_max[1:6]), avg_n_max[7], unique(avg_n_max[8:10]), unique(avg_n_max[11:19]), 
            unique(avg_n_max[20:25]), unique(avg_n_max[26:29]), unique(avg_n_max[30:31]), unique(avg_n_max[32:33]), 
            unique(avg_n_max[34:35]), unique(avg_n_max[36:37]), unique(avg_n_max[38:40]), unique(avg_n_max[41:43]))

         df_out$average_perc <- c(rep(mean(sample_count[1:6])/total_avg,6), avg_n_max[7]/total_avg, 
            rep(mean(sample_count[8:10])/total_avg,3), rep(mean(sample_count[11:19])/total_avg,9), 
            rep(mean(sample_count[20:25])/total_avg,6), rep(mean(sample_count[26:29])/total_avg,4), 
            rep(mean(sample_count[30:31])/total_avg,2), rep(mean(sample_count[32:33])/total_avg,2), 
            rep(mean(sample_count[34:35])/total_avg,2), rep(mean(sample_count[36:37])/total_avg,2), 
            rep(mean(sample_count[38:40])/total_avg,3), rep(mean(sample_count[41:43])/total_avg,3))

         return(df_out)

      }

      cd_all <- getNumGenes(AT_expr_repl, scripttype = "coding", c_level = "all")
      cd_brass <- getNumGenes(AT_expr_repl, scripttype = "coding", c_level = "brass")
      cd_core <- getNumGenes(AT_expr_repl, scripttype = "coding", c_level = "core")
      nc_all <- getNumGenes(AT_expr_repl, scripttype = "lncRNA", c_level = "all")
      nc_core <- getNumGenes(AT_expr_repl, scripttype = "lncRNA", c_level = "core")

      at_stats <- rbind(cd_all, cd_brass, cd_core, nc_all, nc_core)



      # Show message
      message("Writing output...")

      # Set filename
      fname_max_expr <- sprintf('%s.csv', paste(species_id, "max_expr_stats", sep = "_"))

      # Write final data tables to csv files and store them in /out_dir/output/max_expr
      if (!dir.exists(file.path(out_dir, "output", "max_expr"))) 
      dir.create(file.path(out_dir, "output", "max_expr"), recursive = TRUE)

      write.table(at_stats, file = file.path(out_dir, "output", "max_expr", fname_max_expr), 
         sep=";", dec=".", row.names = FALSE, col.names = TRUE)



      # Prepare data for plotting
      at_stats_red <- at_stats[c(1,7,8,11,20,26,30,32,34,36,38,41,87,93,94,97,106,112,116,
         118,120,122,124,127,130,136,137,140,149,155,159,161,163,165,167,170,173,179,180,
         183,192,198,202,204,206,208,210,213),]

      at_stats_pc <- at_stats_red[at_stats_red$biotype == "coding",]
      at_stats_nc <- at_stats_red[at_stats_red$biotype == "lncRNA",]


      
      # Generate plots
      plotMaxExprAT <- function(data, biotype) {

         if (biotype == "coding") {

            p_title <- "Protein-coding genes"
         
         } else if (biotype == "lncRNA") {

            p_title <- "lncRNAs"         
         }

         fname <- sprintf('%s.jpg', paste(deparse(substitute(data)), sep="_"))

         data$biotype <- factor(data$biotype, levels = unique(data$biotype))
         data$conservation <- factor(data$conservation, levels = unique(data$conservation))
         data$group <- factor(data$group, levels = unique(data$group))

         p <- ggplot(data, aes(x = conservation, y = average_perc, color = organ)) + 
         geom_point(size = 5.0, position = position_dodge(width = 0.75), aes(color = group)) + 
         geom_linerange(size = 0.75, aes(x = conservation, ymin = 0, ymax = average_perc, 
            colour = group), position = position_dodge(width = 0.75)) + 
         scale_x_discrete(expand = c(0.05, 0)) + 
         scale_y_continuous(limits = c(0, 0.4), expand = c(0, 0))


         q <- p + 
         # scale_color_manual(values=c("#b2b2b2","#e8a215","#f0d737","#069870","#0770ab","#4fb6f0","#ea6965")) + 
         # Uses a slightly modified colorblind-friendly palette from Wong (Nature Methods, 2011)
         theme_classic() + 
         xlab("") + ylab("Percentage") + ggtitle("") + 
         theme(text = element_text(size = 23.5), 
            panel.grid.major = element_line(colour = "white"), 
            panel.grid.minor = element_line(colour = "white"),  
            axis.ticks.length = unit(.2, "cm"),
            axis.ticks = element_line(colour = "gray15", size = 0.7),
            axis.title.x = element_text(colour = "black", size = 21.55, 
               margin = margin(t = 12.5, r = 0, b = 50.2, l = 0)),  
            axis.title.y = element_text(colour = "black", size = 21.5, 
               margin = margin(t = 0, r = 5.8, b = 0, l = 1.5)), 
            axis.text.x = element_text(colour = "black", margin = margin(t = 3.5, r = 0, b = 1.6, l = 0), size = 20.5), 
            axis.text.y = element_text(colour = "black", margin = margin(t = 0, r = 3.25, b = 0, l = 4), size = 18.55), 
            plot.title = element_text(colour = "black", size = 23.5, 
               margin = margin(t = 35.8, r = 0, b = 11.5, l = 0), hjust = 0.5), 
            plot.margin = unit(c(0, 0.5, 0, 1), "points"),
            legend.position = "right",
            legend.title = element_text(colour = "black", size = 20, face = "bold"),
            legend.text = element_text(size = 20), 
            legend.background = element_rect(fill = NA))

         ggsave(file = file.path(out_dir, "output", "plots", fname), plot = q,
            scale = 1, width = 10.25, height = 7.3, units = c("in"), 
            dpi = 600, limitsize = FALSE)
      }

      plotMaxExprAT(data = at_stats_pc, biotype = "coding")
      plotMaxExprAT(data = at_stats_nc, biotype = "lncRNA")





      # Generate stacked bar charts
      plotMaxExpAT <- function(data, biotype) {

         if (biotype == "coding") {

            p_title <- "Protein-coding"

            p_title_mar <- margin(t = 10, r = 0, b = 10, l = 0)

            x_text_mar <- margin(t = -8.5, r = 0, b = 1.6, l = 0)

            x_labels = c(
            "all" = expression(atop(NA, atop(textstyle('All'), textstyle('Genes')))), 
            "core" = expression(atop(NA, atop(textstyle('Core'), textstyle('Orthologs')))))
         
         } else if (biotype == "lncRNA") {

            p_title <- "lncRNAs"

            p_title_mar <- margin(t = 10, r = 0, b = 12.45, l = 0)

            x_text_mar <- margin(t = -8.35, r = 0, b = 1.6, l = 0)

            x_labels = c(
            "all" = expression(atop(NA, atop(textstyle('All'), textstyle('lncRNAs')))), 
            "core" = expression(atop(NA, atop(textstyle('Brassicaceae'), textstyle('Orthologs')))))       
         }

         fname <- sprintf('%s.jpg', paste(deparse(substitute(data)), "bc", sep="_"))

         data$conservation <- factor(data$conservation, levels = unique(data$conservation))
         data$group <- factor(data$group, levels = c("Stamen", "Carpel", 
            "Petals", "Sepals", "Seed", "Fruit", "Flower", "Apex", "Leaf", "Stem", "Hypocotyl", "Root"))

         p <- ggplot(data, aes(fill = group, x = conservation, y = average_perc)) + 
         geom_bar(position = "stack", stat = "identity", width = 0.725) + 
         scale_x_discrete(expand = c(0.07, 0), labels = x_labels) + 
         scale_y_continuous(expand = c(0, 0), breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1), 
            limits = c(0, 1.0))

         q <- p + 
         scale_fill_manual(values = c("Root" = "#5d4a95", "Hypocotyl" = "#53b0db", 
            "Stem" = "#0c703d", "Leaf" = "#0a9955", "Apex" = "#f0d737", "Flower" = "#e075af", 
            "Sepals" = "#84cd6a", "Petals" = "#ead1c7", "Stamen" = "#ee412e", 
            "Carpel" = "#f2a72f", "Fruit" = "#b54185", "Seed" = "#e9a3b3")) + 
         # Uses a slightly modified colorblind-friendly palette from Wong (Nature Methods, 2011)
         theme_classic() + 
         guides(fill = guide_legend(override.aes = list(size = 7.5))) + 
         xlab("") + ylab("Fraction") + ggtitle(p_title) + 
         theme(text = element_text(size = 22.5), 
            panel.grid.major = element_line(colour = "white"), 
            panel.grid.minor = element_line(colour = "white"),  
            axis.ticks.length = unit(0.26, "cm"), 
            axis.line = element_line(colour = "black", size = 0.8), 
            axis.ticks = element_line(colour = "black", size = 0.8),
            axis.title.x = element_text(colour = "black", size = 21.0, 
               margin = margin(t = 12.5, r = 0, b = 10, l = 0)),  
            axis.title.y = element_text(colour = "black", size = 21.0, 
               margin = margin(t = 0, r = 5.0, b = 0, l = 1.5)), 
            axis.text.x = element_text(colour = "black", margin = x_text_mar, size = 19.5), 
            axis.text.y = element_text(colour = "black", margin = margin(t = 0, r = 3.25, b = 0, l = 4), size = 18.7), 
            plot.title = element_text(colour = "black", size = 21.0, 
               margin = p_title_mar, hjust = 0.5), 
            plot.margin = unit(c(74, 350, 74, 10), "points"),
            legend.position = "right",
            legend.box.margin = margin(-4, 0, 0, -14), 
            legend.title = element_blank(),
            legend.text = element_text(size = 19.5), 
            legend.background = element_rect(fill = NA))

         ggsave(file = file.path(out_dir, "output", "plots", fname), plot = q,
            scale = 1, width = 10.25, height = 7.3, units = c("in"), 
            dpi = 600, limitsize = FALSE)
      }

      plotMaxExpAT(data = at_stats_pc, biotype = "coding")
      plotMaxExpAT(data = at_stats_nc, biotype = "lncRNA")


   }


   # ----------------- Perform analysis for all species across comparative organs ------------------


   else if (species_id == "all") {

      expr_table_ls <- list(AT_expr = AT_expr, AL_expr = AL_expr, CR_expr = CR_expr, ES_expr = ES_expr, 
         TH_expr = TH_expr, MT_expr = MT_expr, BD_expr = BD_expr)



      # Calculate threshold and get expressed genes for AT
      # Need to do this because less samples are used for analysis than in normalized tables
      # This will filter out genes that are expressed below threshold

      # First select comparative organs
      AT_expr_compl <- dplyr::select(AT_expr_compl, c("id", "biotype", "source", "root_whole_root_5d_1", 
         "root_whole_root_5d_2", "root_whole_root_5d_3", "hypocotyl_10d_1", "hypocotyl_10d_2", 
         "hypocotyl_10d_3", "leaf_1.2_10d_1", "leaf_1.2_10d_2", "leaf_1.2_10d_3", "apex_vegetative_7d_1", 
         "apex_vegetative_7d_2", "apex_vegetative_7d_3", "apex_inflorescence_21d_1", "apex_inflorescence_21d_2", 
         "apex_inflorescence_21d_3", "flower_stg12_21d._1", "flower_stg12_21d._2", "flower_stg12_21d._3", 
         "flower_stg12_stamens_21d._1", "flower_stg12_stamens_21d._2", "flower_stg12_stamens_21d._3", 
         "flower_early_stg12_carpels_21d._1", "flower_early_stg12_carpels_21d._2", "flower_early_stg12_carpels_21d._3"))


      getTH <- function(x) {

         # Extract ERCC data
         ERCC <- x[x$id %like% "ERCC", ]

         # Get sample-specific expr threshold
         ERCC_cutoff <- cbind(
            data.frame(id = "expr_cutoff"), data.frame(biotype = "<NA>"), data.frame(source = "<NA>"), 
            as.data.frame(t(data.frame(
               expr_cutoff = apply(ERCC[,4:ncol(ERCC)], 2, function(i)quantile(i[i>0], 0.05)))
            ))
         )

         # Bind sample-specific threshold to expression table
         all_genes_expr_cutoff <- rbind(x, ERCC_cutoff)


         applyThreshold <- function(express_df) {

            # Add keys to data frame
            key <- seq(1, nrow(express_df), 1)
            express_kdf <- cbind(as.data.frame(key),express_df)

            # Replace all values with "0" that are below sample threshold (either 0 or ERCC)
            getSampleTH <- function(df) {

            # Split data frame by sample replicates into a list then apply threshold for each subset

               th_replicates <- do.call(cbind, lapply(split.default(df[5:ncol(df)], #adjust columns
                        rep(seq_along(df), each = 1, length.out = ncol(df)-4)),
                        function(x) {
                           x[x <= x[nrow(df),], ] <- 0;
                           x
                        }
                     ))

               # Bind key/id/prt_id/symbol/biotype/source columns to thresholded data frame
               th_replicates <- cbind(df[1:4], th_replicates)

               return(th_replicates)
            }

            df_exp <- getSampleTH(express_kdf)

            # Define threshold function
            # This function will remove all rows that do not show expression in at least two of three
            # replicates in at least one sample type
            getThreshold <- function(df) {

            # Split data frame by sample replicates into a list then apply threshold for each subset

               th_replicates <- do.call(cbind, lapply(split.default(df[5:ncol(df)], #adjust columns
                        rep(seq_along(df), each = 3, length.out = ncol(df)-4)), #adjust columns
                        function(x) {
                           x[rowSums(x > 0) < 2, ] <- 0; 
                           x
                        }
                     ))

               # Bind key/id/prt_id/symbol/biotype/source columns to thresholded data frame
               th_replicates <- cbind(df[1:4], th_replicates)

               # Remove all rows that only contain "0"
               th_replicates <- th_replicates[which(rowSums(th_replicates[,-1:-4, drop = FALSE] > 0) > 0),]

               return(th_replicates)
            }

            # Apply threshold to data and extract keys ("key")
            keys_data_repl <- getThreshold(df_exp)
            keys_data <- keys_data_repl[,1:2]
            names(keys_data) <- c("key","ID")

            # Generate thresholded data frame based on keys
            th_df <- merge(keys_data, express_kdf, by="key")
            th_df <- th_df[order(th_df$key),]
            th_df <- th_df[-1:-2]

            return(th_df)
         }

         return_objects_th <- applyThreshold(all_genes_expr_cutoff)

         return(return_objects_th)

      }

      AT_expr_compl <- getTH(AT_expr_compl)

      
      # Get replicate expression
      calculateAvgExpr <- function(df) {

         # Split data frame by sample replicates into a list and get rowMeans for each subset
         averaged_replicates <- do.call(cbind, lapply(split.default(df, 
            rep(seq_along(df), 
               each = 3, 
               length.out = ncol(df))
            ), rowMeans)
         )

         averaged_replicates <- as.data.frame(averaged_replicates)

         # Get replicate names
         getReplNames <- function(n) {

            df_names <- names(n)
            df_names <- unique(substring(df_names, 1, nchar(df_names)-4))
            df_names <- gsub('[.]', '', df_names)

            return(df_names)
         }

         repl_names <- getReplNames(df)

         colnames(averaged_replicates) <- repl_names

         return(averaged_replicates)
      }

      expr_table_repl_ls <- lapply(expr_table_ls, calculateAvgExpr)


      # Select comparative organs for AT and AL
      expr_table_repl_ls$AT_expr <- dplyr::select(expr_table_repl_ls$AT_expr, c("root_whole_root_5d", 
         "hypocotyl_10d", "leaf_12_7d", "apex_vegetative_7d", "apex_inflorescence_21d", 
         "flower_stg12_21d", "flower_stg12_stamens_21d", "flower_early_stg12_carpels_21d"))

      expr_table_repl_ls$AL_expr <- expr_table_repl_ls$AL_expr[, -which(names(expr_table_repl_ls$AL_expr) %in% c(
         "flower_stg11_stamens_8w10w25d", "flower_early_stg12_stamens_8w10w23d", 
         "flower_late_stg12_stamens_8w10w21d"))]


      # Add dataset identifier to list elements
      spec_exp_names <- lapply(seq_along(expr_table_repl_ls), function(i) { 
         paste(names(expr_table_repl_ls)[[i]])
      })

      for (i in seq_along(expr_table_repl_ls)) {
         expr_table_repl_ls[[i]]$dataset <- rep(spec_exp_names[i], nrow(expr_table_repl_ls[[i]]))
      }



      # Get number of genes with maximum expression for each organ and species
      # For AT, use list of lncRNAs since some have AT identifier, some have "lnc..."
      # For other species, lncRNA genes can be identified by "lnc..." in id

      getNumGenes <- function(df, scripttype = c("coding", "lncRNA"), c_level = c("all", "brass", "core")) {

         spec_id <- unique(sub("\\_.*", "", df$dataset))
         df <- within(df, rm(dataset))

         Core_expr <- Core_expr[!grepl("ERCC", Core_expr[,1]),]

         AT_lncRNA_ids <- subset(AT_expr_compl, subset = biotype %in% c("lnc_exonic_antisense", 
               "lnc_intronic_antisense", "lnc_intergenic"))[,1]

         # Get protein-coding core IDs for non-AT species

         if (spec_id == "AL") {

            core_ids <- as.data.frame(sapply(Core_expr[,1], function(x) unlist(strsplit(x, "\\:"))[2]))
            core_brass_ids <- as.data.frame(sapply(Brass_pc_expr[,1], function(x) unlist(strsplit(x, "\\:"))[2]))
            core_lnc_ids <- as.data.frame(sapply(Brass_nc_expr[,1], function(x) unlist(strsplit(x, "\\:"))[2]))

         } else if (spec_id == "CR") {

            core_ids <- as.data.frame(sapply(Core_expr[,1], function(x) unlist(strsplit(x, "\\:"))[3]))
            core_brass_ids <- as.data.frame(sapply(Brass_pc_expr[,1], function(x) unlist(strsplit(x, "\\:"))[3]))
            core_lnc_ids <- as.data.frame(sapply(Brass_nc_expr[,1], function(x) unlist(strsplit(x, "\\:"))[3]))

         } else if (spec_id == "ES") {

            core_ids <- as.data.frame(sapply(Core_expr[,1], function(x) unlist(strsplit(x, "\\:"))[4]))
            core_brass_ids <- as.data.frame(sapply(Brass_pc_expr[,1], function(x) unlist(strsplit(x, "\\:"))[4]))
            core_lnc_ids <- as.data.frame(sapply(Brass_nc_expr[,1], function(x) unlist(strsplit(x, "\\:"))[4]))

         } else if (spec_id == "TH") {

            core_ids <- as.data.frame(sapply(Core_expr[,1], function(x) unlist(strsplit(x, "\\:"))[5]))

         } else if (spec_id == "MT") {

            core_ids <- as.data.frame(sapply(Core_expr[,1], function(x) unlist(strsplit(x, "\\:"))[6]))

         } else if (spec_id == "BD") {

            core_ids <- as.data.frame(sapply(Core_expr[,1], function(x) unlist(strsplit(x, "\\:"))[7]))
         }

         # ---------------------------- Process AT data ----------------------------

         if ((spec_id == "AT") && (c_level == "all") && (scripttype == "coding")) {

            df_at1 <- df[rownames(df) %like% "AT1G", ]
            df_at2 <- df[rownames(df) %like% "AT2G", ]
            df_at3 <- df[rownames(df) %like% "AT3G", ]
            df_at4 <- df[rownames(df) %like% "AT4G", ]
            df_at5 <- df[rownames(df) %like% "AT5G", ]

            df <- rbind(df_at1, df_at2, df_at3, df_at4, df_at5)

         } else if ((spec_id == "AT") && (c_level == "all") && (scripttype == "lncRNA")) {

            df <- df[rownames(df) %in% AT_lncRNA_ids, ]

         # Reduce data to core orthologs if core set is chosen

         } else if ((spec_id == "AT") && (c_level == "core") && (scripttype == "coding")) {

            # Extract AT protein-coding core IDs
            core_ids <- sub("\\:.*", "", Core_expr[,1])

            df <- df[rownames(df) %in% core_ids,]

         } else if ((spec_id == "AT") && (c_level == "brass") && (scripttype == "coding")) {

            # Get Brassicaceae protein-coding orthologs

            core_brass_ids <- sub("\\:.*", "", Brass_pc_expr[,1])

            df <- df[rownames(df) %in% core_brass_ids, ]
         
         } else if ((spec_id == "AT") && (c_level == "brass") && (scripttype == "lncRNA")) {

            # Extract AT lncRNA IDs and get Brassicaceae core lncRNAs

            core_lnc_ids <- sub("\\:.*", "", Brass_nc_expr[,1])

            df <- df[rownames(df) %in% core_lnc_ids, ]

         
         # -------------------- Process data from other species --------------------

         } else if ((spec_id != "AT") && (c_level == "all") && (scripttype == "coding")) {

            # Negate like function
            `%!like%` = Negate(`%like%`)

            df <- df[rownames(df) %!like% "ERCC|lnc|LTR", ]

         } else if ((spec_id != "AT") && (c_level == "all") && (scripttype == "lncRNA")) {

            df <- df[rownames(df) %like% "lnc", ]

         # Reduce data to Brassicaceae or core orthologs if brass/core set is chosen
            
         } else if ((spec_id != "AT") && (c_level == "core") && (scripttype == "coding")) {

            df <- df[rownames(df) %in% as.character(core_ids[,1]),]

         } else if ((spec_id != "AT") && (c_level == "brass") && (scripttype == "coding")) {

            df <- df[rownames(df) %in% as.character(core_brass_ids[,1]),]
         
         } else if ((spec_id != "AT") && (c_level == "brass") && (scripttype == "lncRNA")) {

            df <- df[rownames(df) %in% as.character(core_lnc_ids[,1]), ]
         }


         # Create group names
         groups <- c("Root", "Hypocotyl", "Leaf", "Apex_veg", "Apex_inf", "Flower", "Stamen", "Carpel")

         df$max <- colnames(df)[apply(df, 1, which.max)]

         x_max <- df$max
         sample_count <- do.call(rbind, lapply(colnames(df)[1:(length(df)-1)], function(x) {
            length(grep(x, x_max))}))
         perc <- sample_count/sum(sample_count)


         df_out <- data.frame(
            species = spec_id, 
            biotype = rep(scripttype), 
            conservation = rep(c_level),
            organ = colnames(df)[1:(length(df)-1)], 
            group = groups, 
            count = sample_count, 
            fraction = perc
            )

         return(df_out)

      }

      cd_all <- do.call(rbind, lapply(expr_table_repl_ls, getNumGenes, scripttype = "coding", 
         c_level = "all"))
      cd_brass <- do.call(rbind, lapply(expr_table_repl_ls[1:4], getNumGenes, scripttype = "coding", 
         c_level = "brass"))
      cd_core <- do.call(rbind, lapply(expr_table_repl_ls, getNumGenes, scripttype = "coding", 
         c_level = "core"))
      nc_all <- do.call(rbind, lapply(expr_table_repl_ls, getNumGenes, scripttype = "lncRNA", 
         c_level = "all"))
      nc_brass <- do.call(rbind, lapply(expr_table_repl_ls[1:4], getNumGenes, scripttype = "lncRNA", 
         c_level = "core")) # Ortholog lncRNA dataset is limited to Brassicaceae

      pc_stats <- rbind(cd_all, cd_core)
      nc_stats <- rbind(nc_all, nc_brass)



      # Generate plots
      plotMaxExprOS <- function(data, biotype) {

         if (biotype == "coding") {

            p_title <- "Organ with highest expression of protein-coding genes"

            y_scale <- c(0, 0.3125)

            plt_mar <- c(0.5, 1.75, 1.75, 1.75)
         
         } else if (biotype == "lncRNA") {

            p_title <- "Organ with highest expression of lncRNAs" 

            y_scale <- c(0, 0.44)

            plt_mar <- c(0.5, 1.75, 1.75, 1.75)
         
         } else if (biotype == "BrlncRNA") {

            p_title <- "Organ with highest expression of lncRNAs (n = 307)" 

            y_scale <- c(0, 0.34)

            plt_mar <- c(0.5, 30.52, 1.75, 1.75)
         }

         fname <- sprintf('%s.jpg', paste(deparse(substitute(data)), sep="_"))

         x_lab <- c(Root = "Rt", Hypocotyl = "Hc", Leaf = "Lf", Apex_veg = "Av", 
            Apex_inf = "Ai", Flower = "Fl", Stamen = "St", Carpel = "Ca")

         data$biotype <- factor(data$biotype, levels = unique(data$biotype))
         data$conservation <- factor(data$conservation, levels = unique(data$conservation))
         data$group <- factor(data$group, levels = unique(data$group))
         data$species <- factor(data$species, levels = unique(data$species))

         p <- ggplot(data, aes(x = group, y = fraction, color = group)) + 
         geom_segment(aes(y = 0, yend = fraction, xend = group), size = 2.5, colour = "grey77") + 
         geom_point(size = 7.7, position = position_dodge(width = 0.75), aes(color = group)) +
         scale_x_discrete(expand = c(0.05, 0), labels = x_lab) + 
         scale_y_continuous(limits = y_scale, expand = c(0, 0))
         q <- p + 
         scale_color_manual(values = c("Root" = "#6a54a9", "Hypocotyl" = "#53b0db", 
            "Leaf" = "#0a9955", "Apex_veg" = "#96ba37", "Apex_inf" = "#f0d737", 
            "Flower" = "#e075af", "Stamen" = "#ed311c", "Carpel" = "#f2a72f")) + 
         # Uses a slightly modified colorblind-friendly palette from Wong (Nature Methods, 2011)
         theme_classic() + 
         xlab("") + ylab("Fraction") + ggtitle(p_title) + 
         theme(text = element_text(size = 23.5), 
            strip.text = element_text(size = 24.5, face = "plain"), 
                strip.text.x = element_text(margin = margin(0.45, 0, 0.45, 0, "cm")), 
                strip.background = element_rect(colour = 'black', fill = NA, size = 3.125), 
                axis.ticks.length = unit(0.25, "cm"), 
                axis.ticks = element_line(colour = "black", size = 1.5), 
                axis.line = element_line(colour = 'black', size = 1.5), 
                plot.margin = unit(plt_mar, "cm"), 
                axis.title.y = element_text(size = 25, margin = margin(t = 0, r = 8, b = 0, l = 1), 
                    colour = "black", face = "bold"), 
                axis.title.x = element_text(size = 25, margin = margin(t = 6.5, r = 0, b = 5.75, l = 0), 
                    colour = "black", face = "bold"), 
                axis.text.x = element_text(size=21.5, margin = margin(t = 4, b = 7.75), colour = "grey35", 
                    angle = 0, vjust = 1, hjust = 0.5), 
                axis.text.y = element_text(size = 21.5, angle = 0, margin = margin(l = 0.75, r = 1.5), colour = "grey35"), 
                plot.title = element_text(size = 25.5, margin = margin(t = 5.5, b = 15.3), face = "plain"), 
                panel.spacing = unit(0.55, "cm"), 
                panel.grid.major = element_blank(),
                panel.grid.minor.x = element_blank(), 
                panel.grid.minor.y = element_blank(),  
                legend.position ="none")

         q <- q + facet_wrap(~ factor(species, levels = c("AT", "AL", "CR", "ES", "TH", "MT", "BD")) , nrow = 1, scales = "free_x")

            ggsave(file = file.path(out_dir, "output", "plots", fname), plot = q, 
                width = 28.5, height = 6.5, dpi = 300, units = c("in"), limitsize = FALSE)
      }

      plotMaxExprOS(data = cd_all, biotype = "coding")
      plotMaxExprOS(data = nc_all, biotype = "lncRNA")
      plotMaxExprOS(data = nc_brass, biotype = "BrlncRNA")



      # Show message
      message("Writing output...")

      # Set filename
      fname_max_pc <- sprintf('%s.csv', paste(species_id, "max_expr_pc_stats", sep = "_"))
      fname_max_nc <- sprintf('%s.csv', paste(species_id, "max_expr_nc_stats", sep = "_"))

      # Write final data tables to csv files and store them in /out_dir/output/max_expr
      if (!dir.exists(file.path(out_dir, "output", "max_expr"))) 
      dir.create(file.path(out_dir, "output", "max_expr"), recursive = TRUE)

      write.table(pc_stats, file = file.path(out_dir, "output", "max_expr", fname_max_pc), 
         sep=";", dec=".", row.names = FALSE, col.names = TRUE)

      write.table(nc_stats, file = file.path(out_dir, "output", "max_expr", fname_max_nc), 
         sep=";", dec=".", row.names = FALSE, col.names = TRUE)




      #-- Analyse distribution of max expression for coding+non-coding genes across species --


      expr_table_ls_br <- expr_table_ls[c(1:4)]

      expr_table_repl_ls <- lapply(expr_table_ls_br, calculateAvgExpr)


      # Select comparative organs for AT and AL
      expr_table_repl_ls$AT_expr <- dplyr::select(expr_table_repl_ls$AT_expr, c("root_whole_root_5d", 
         "hypocotyl_10d", "leaf_12_7d", "apex_vegetative_7d", "apex_inflorescence_21d", 
         "flower_stg12_21d", "flower_stg12_stamens_21d", "flower_early_stg12_carpels_21d"))


      # Add dataset identifier to list elements
      spec_exp_names <- lapply(seq_along(expr_table_repl_ls), function(i) { 
         paste(names(expr_table_repl_ls)[[i]])
      })

      for (i in seq_along(expr_table_repl_ls)) {
         expr_table_repl_ls[[i]]$dataset <- rep(spec_exp_names[i], nrow(expr_table_repl_ls[[i]]))
      }


      
      # Get maximum expression values for coding and non-coding genes for Brassicaceae species

      getMaxExprDist <- function(df, scripttype = c("coding", "lncRNA"), c_level = c("all", "non-core" , "core")) {

         spec_id <- unique(sub("\\_.*", "", df$dataset))
         df <- within(df, rm(dataset))

         Core_expr <- Core_expr[!grepl("ERCC", Core_expr[,1]),]

         `%nin%` = Negate(`%in%`)


         # Get protein-coding and non-coding core IDs for non-AT species

         if (spec_id == "AT") {

            compl_table <- AT_expr_compl
            core_ids <- as.data.frame(sapply(Core_expr[,1], function(x) unlist(strsplit(x, "\\:"))[1]))
            core_lnc_ids <- as.data.frame(sapply(Brass_nc_expr[,1], function(x) unlist(strsplit(x, "\\:"))[1]))

         } else if (spec_id == "AL") {

            compl_table <- AL_expr_compl
            core_ids <- as.data.frame(sapply(Core_expr[,1], function(x) unlist(strsplit(x, "\\:"))[2]))
            core_lnc_ids <- as.data.frame(sapply(Brass_nc_expr[,1], function(x) unlist(strsplit(x, "\\:"))[2]))

         } else if (spec_id == "CR") {

            compl_table <- CR_expr_compl
            core_ids <- as.data.frame(sapply(Core_expr[,1], function(x) unlist(strsplit(x, "\\:"))[3]))
            core_lnc_ids <- as.data.frame(sapply(Brass_nc_expr[,1], function(x) unlist(strsplit(x, "\\:"))[3]))

         } else if (spec_id == "ES") {

            compl_table <- ES_expr_compl
            core_ids <- as.data.frame(sapply(Core_expr[,1], function(x) unlist(strsplit(x, "\\:"))[4]))
            core_lnc_ids <- as.data.frame(sapply(Brass_nc_expr[,1], function(x) unlist(strsplit(x, "\\:"))[4]))

         }

         # -------------------------- Process protein-coding data --------------------------

         if ((c_level == "all") && (scripttype == "coding")) {

            df <- df[rownames(df) %in% compl_table[compl_table$biotype == "protein_coding",]$id, ]


         # ------------------------------ Process lncRNA data ------------------------------

         } else if ((c_level == "all") && (scripttype == "lncRNA")) {

            df <- df[rownames(df) %in% subset(compl_table, subset = biotype %in% c(
               "lnc_exonic_antisense", "lnc_intronic_antisense", "lnc_intergenic"))$id, ]


         # ----------------- Process non-core ortholog protein-coding data -----------------

         } else if ((c_level == "non-core") && (scripttype == "coding")) {

            df <- df[rownames(df) %in% compl_table[compl_table$biotype == "protein_coding",]$id, ]
            df <- df[rownames(df) %nin% as.character(core_ids[,1]),]


         # --------------------- Process non-core ortholog lncRNA data ---------------------
         
         } else if ((c_level == "non-core") && (scripttype == "lncRNA")) {

            df <- df[rownames(df) %in% subset(compl_table, subset = biotype %in% c(
               "lnc_exonic_antisense", "lnc_intronic_antisense", "lnc_intergenic"))$id, ]
            df <- df[rownames(df) %nin% as.character(core_lnc_ids[,1]), ]


         # ------------------- Process core ortholog protein-coding data -------------------

         } else if ((c_level == "core") && (scripttype == "coding")) {

            df <- df[rownames(df) %in% as.character(core_ids[,1]),]


         # ----------------------- Process core ortholog lncRNA data -----------------------
         
         } else if ((c_level == "core") && (scripttype == "lncRNA")) {

            df <- df[rownames(df) %in% as.character(core_lnc_ids[,1]), ]
         }


         # Get max expression value per gene
         df$max <- apply(df, 1, max)

         df_out <- data.frame(
            species = rep(spec_id), 
            biotype = rep(scripttype), 
            conservation = rep(c_level), 
            class = rep(paste(scripttype, c_level, sep = "_")),
            max_expr = df$max
            )

         return(df_out)

      }

      cd_expr_dist_all <- do.call(rbind, lapply(expr_table_repl_ls, getMaxExprDist, scripttype = "coding", 
         c_level = "all"))
      cd_expr_dist_n_core <- do.call(rbind, lapply(expr_table_repl_ls, getMaxExprDist, scripttype = "coding", 
         c_level = "non-core"))
      cd_expr_dist_core <- do.call(rbind, lapply(expr_table_repl_ls, getMaxExprDist, scripttype = "coding", 
         c_level = "core"))
      nc_expr_dist_all <- do.call(rbind, lapply(expr_table_repl_ls, getMaxExprDist, scripttype = "lncRNA", 
         c_level = "all"))
      nc_expr_dist_n_core <- do.call(rbind, lapply(expr_table_repl_ls, getMaxExprDist, scripttype = "lncRNA", 
         c_level = "non-core"))
      nc_expr_dist_brass <- do.call(rbind, lapply(expr_table_repl_ls[1:4], getMaxExprDist, scripttype = "lncRNA", 
         c_level = "core")) # Ortholog lncRNA dataset is limited to Brassicaceae


      # Combine data
      max_expr_dist <- rbind(cd_expr_dist_all, cd_expr_dist_n_core, cd_expr_dist_core, nc_expr_dist_all, 
         nc_expr_dist_n_core, nc_expr_dist_brass)




      # Generate plots
      plotMaxExprDist <- function(data) {

         #if (biotype == "coding") {

           # p_title <- "Organ with highest expression of protein-coding genes"

           y_scale <- c(-0.5, 13.75)

           plt_mar <- c(0.1, 2.75, 1.5, 1.75)
         
        # } else if (biotype == "lncRNA") {


         fname <- sprintf('%s.jpg', paste(deparse(substitute(data)), sep="_"))

         x_lab <- c(Root = "Rt", Hypocotyl = "Hc", Leaf = "Lf", Apex_veg = "Av", 
            Apex_inf = "Ai", Flower = "Fl", Stamen = "St", Carpel = "Ca")

         data$class <- factor(data$class, levels = unique(data$class))
         data$conservation <- factor(data$conservation, levels = unique(data$conservation))
         data$species <- factor(data$species, levels = unique(data$species))

         p <- ggplot(data, aes(x = class, y = max_expr, color = class)) + 
         geom_boxplot(aes(fill = class), colour = "black", width = 0.85) +
         scale_x_discrete(expand = c(0.055, 0)) + 
         scale_y_continuous(limits = y_scale, expand = c(0, 0), breaks = c(0,2.5,5,7.5,10,12.5))
         q <- p + 
         #scale_color_manual(values = c("Root" = "#6a54a9", "Hypocotyl" = "#53b0db", 
         #   "Leaf" = "#0a9955", "Apex_veg" = "#96ba37")) + 
         # Uses a slightly modified colorblind-friendly palette from Wong (Nature Methods, 2011)
         theme_classic() + 
         xlab("") + ylab("Maximum expression") + ggtitle("") + 
         theme(text = element_text(size = 23.5), 
            strip.text = element_text(size = 24.5, face = "plain"), 
                strip.text.x = element_text(margin = margin(0.25, 0, 0.25, 0, "cm")), 
                strip.background = element_rect(colour = 'white', fill = NA, size = 0.1), 
                axis.ticks.length = unit(0.25, "cm"), 
                axis.ticks = element_line(colour = "black", size = 1.5), 
                axis.line = element_line(colour = 'black', size = 1.5), 
                plot.margin = unit(plt_mar, "cm"), 
                axis.title.y = element_text(size = 25, margin = margin(t = 0, r = 8, b = 0, l = 1), 
                    colour = "black", face = "bold"), 
                axis.title.x = element_text(size = 25, margin = margin(t = 6.5, r = 0, b = 5.75, l = 0), 
                    colour = "black", face = "bold"), 
                axis.text.x = element_text(size=21.5, margin = margin(t = 4, b = 7.75), colour = "grey35", 
                    angle = 0, vjust = 1, hjust = 0.5), 
                axis.text.y = element_text(size = 21.5, angle = 0, margin = margin(l = 0.75, r = 1.5), colour = "grey35"), 
                panel.spacing = unit(0.8, "cm"), 
                panel.grid.major = element_blank(),
                panel.grid.minor.x = element_blank(), 
                panel.grid.minor.y = element_blank(),  
                legend.position ="none")

         q <- q + facet_wrap(~ factor(species, levels = c("AT", "AL", "CR", "ES")) , nrow = 1, scales = "free_x")

            ggsave(file = file.path(out_dir, "output", "plots", fname), plot = q, 
                width = 20, height = 6.5, dpi = 300, units = c("in"), limitsize = FALSE)
      }

      plotMaxExprDist(data = max_expr_dist)


      # Wilcoxon rank sum test with continuity correction (all genes vs core genes)
      # W = 1729000, p-value < 2.2e-16 for comparison cd/lnc for all species




   }



}





