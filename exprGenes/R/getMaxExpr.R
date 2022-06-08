# Count number of genes with maximum expression in each organ in Arabidopsis thaliana for all organs,
# and additionally for only the comparative organs; do the same for all the other species
# The comparative organs are: root, hypocotyl, leaf, apex veg, apex inf, flower, stamen, carpel
# Data input: Normalized TPM expression data containing protein-coding genes, NATs and lincRNAs
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
   pathAT = file.path(in_dir, "Expression_data", "AT_genes_inter_norm_tpm_mat_deseq_sample_names.csv")
   pathAL = file.path(in_dir, "Expression_data", "AL_genes_inter_norm_tpm_mat_deseq_sample_names.csv")
   pathCR = file.path(in_dir, "Expression_data", "CR_genes_inter_norm_tpm_mat_deseq_sample_names.csv")
   pathES = file.path(in_dir, "Expression_data", "ES_genes_inter_norm_tpm_mat_deseq_sample_names.csv")
   pathTH = file.path(in_dir, "Expression_data", "TH_genes_inter_norm_tpm_mat_deseq_sample_names.csv")
   pathMT = file.path(in_dir, "Expression_data", "MT_genes_inter_norm_tpm_mat_deseq_sample_names.csv")
   pathBD = file.path(in_dir, "Expression_data", "BD_genes_inter_norm_tpm_mat_deseq_sample_names.csv")

   # Ortholog tables
   pathCore = file.path(in_dir, "Expression_data", "AT_core_inter_tpm_mat_deseq_sample_names.csv")
   pathPcBrass = file.path(in_dir, "Expression_data", "AT_brass_inter_tpm_mat_deseq_sample_names.csv")
   pathNcBrass = file.path(in_dir, "Expression_data", "lnc_AT_brass_inter_tpm_mat_deseq_sample_names.csv")

   # A.thaliana table containing raw expression values and biotype annotation
   pathAT_compl <- file.path(in_dir, "Expression_data", "AT_genes_complete_table_tpm_sample_names.csv")


   # Set up list of expression tables
   if (species == "all") {

      expr_table_ls <- list(AT_tpm = pathAT, AL_tpm = pathAL, CR_tpm = pathCR, ES_tpm = pathES, 
         TH_tpm = pathTH, MT_tpm = pathMT, BD_tpm = pathBD, Core_tpm = pathCore, Brass_pc_tpm = pathPcBrass, 
         Brass_nc_tpm = pathNcBrass, AT_tpm_compl = pathAT_compl)

   } else if (species == "AT") {

      expr_table_ls <- list(AT_tpm = pathAT, Core_tpm = pathCore, Brass_pc_tpm = pathPcBrass, 
         Brass_nc_tpm = pathNcBrass, AT_tpm_compl = pathAT_compl)
   }


   # Read expression data
   expr_tables <- lapply(expr_table_ls, function(x) {
      read.table(x, sep=";", dec=".", header=TRUE, stringsAsFactors=FALSE)
      })


   # Stop function here to run tests
   # return_list <- list("species_id" = species_id, "expr_tables" = expr_tables)
   # return(return_list)
   # }
   # return_objects <- getMaxExpr(species="AT") # read in TPM expression data
   # list2env(return_objects, envir = .GlobalEnv)


   list2env(expr_tables, envir = .GlobalEnv)

   # Show message
   message("Starting analysis...")




   # ---------------------- Perform analysis for A.thaliana across all organs ----------------------


   if (species_id == "AT") {

      # log-transform data
      AT_tpm_log <- log2(AT_tpm + 1)


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

      AT_tpm_repl <- calculateAvgExpr(AT_tpm_log)


      # Get replicate names
      getReplNames <- function(n) {

         df_names <- names(n)
         df_names <- unique(substring(df_names, 1, nchar(df_names)-4))
         df_names <- gsub('[.]', '', df_names)

         return(df_names)
      }

      repl_names <- getReplNames(AT_tpm_log)

      colnames(AT_tpm_repl) <- repl_names


      # Extract protein-coding core ortholog and Brassicaceae lncRNA ortholog gene IDs
      core_ids <- sub("\\:.*", "", Core_tpm[,1])
      core_ids <- core_ids[!grepl("ERCC", core_ids)]


      # Extract all A.thaliana lncRNAs and get orthologs
      AT_lncRNA_ids <- subset(AT_tpm_compl, subset = biotype %in% c("lnc_exonic_antisense", 
         "lnc_intronic_antisense", "lnc_intergenic"))[,1]
      core_brass_ids <- sub("\\:.*", "", Brass_pc_tpm[,1])
      core_lnc_ids <- sub("\\:.*", "", Brass_nc_tpm[,1])


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

         return(df_out)

      }

      cd_all <- getNumGenes(AT_tpm_repl, scripttype = "coding", c_level = "all")
      cd_brass <- getNumGenes(AT_tpm_repl, scripttype = "coding", c_level = "brass")
      cd_core <- getNumGenes(AT_tpm_repl, scripttype = "coding", c_level = "core")
      nc_all <- getNumGenes(AT_tpm_repl, scripttype = "lncRNA", c_level = "all")
      nc_core <- getNumGenes(AT_tpm_repl, scripttype = "lncRNA", c_level = "core")

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


   }


   # ----------------- Perform analysis for all species across comparative organs ------------------


   else if (species_id == "all") {

      tpm_table_ls <- list(AT_tpm = AT_tpm, AL_tpm = AL_tpm, CR_tpm = CR_tpm, ES_tpm = ES_tpm, 
         TH_tpm = TH_tpm, MT_tpm = MT_tpm, BD_tpm = BD_tpm)

      # log-transform data
      tpm_table_ls <- lapply(tpm_table_ls, function(t) log2(t + 1))

      
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

      tpm_table_repl_ls <- lapply(tpm_table_ls, calculateAvgExpr)


      # Select comparative organs for AT and AL
      tpm_table_repl_ls$AT_tpm <- dplyr::select(tpm_table_repl_ls$AT_tpm, c("root_whole_root_5d", 
         "hypocotyl_10d", "leaf_12_7d", "apex_vegetative_7d", "apex_inflorescence_21d", 
         "flower_stg12_21d", "flower_stg12_stamens_21d", "flower_early_stg12_carpels_21d"))

      tpm_table_repl_ls$AL_tpm <- tpm_table_repl_ls$AL_tpm[, -which(names(tpm_table_repl_ls$AL_tpm) %in% c(
         "flower_stg11_stamens_8w10w25d", "flower_early_stg12_stamens_8w10w23d", 
         "flower_late_stg12_stamens_8w10w21d"))]


      # Add dataset identifier to list elements
      spec_exp_names <- lapply(seq_along(tpm_table_repl_ls), function(i) { 
         paste(names(tpm_table_repl_ls)[[i]])
      })

      for (i in seq_along(tpm_table_repl_ls)) {
         tpm_table_repl_ls[[i]]$dataset <- rep(spec_exp_names[i], nrow(tpm_table_repl_ls[[i]]))
      }



      # Get number of genes with maximum expression for each organ and species
      # For AT, use list of lncRNAs since some have AT identifier, some have "lnc..."
      # For other species, lncRNA genes can be identified by "lnc..." in id

      getNumGenes <- function(df, scripttype = c("coding", "lncRNA"), c_level = c("all", "brass", "core")) {

         spec_id <- unique(sub("\\_.*", "", df$dataset))
         df <- within(df, rm(dataset))

         Core_tpm <- Core_tpm[!grepl("ERCC", Core_tpm[,1]),]

         AT_lncRNA_ids <- subset(AT_tpm_compl, subset = biotype %in% c("lnc_exonic_antisense", 
               "lnc_intronic_antisense", "lnc_intergenic"))[,1]

         # Get protein-coding core IDs for non-AT species

         if (spec_id == "AL") {

            core_ids <- as.data.frame(sapply(Core_tpm[,1], function(x) unlist(strsplit(x, "\\:"))[2]))
            core_brass_ids <- as.data.frame(sapply(Brass_pc_tpm[,1], function(x) unlist(strsplit(x, "\\:"))[2]))
            core_lnc_ids <- as.data.frame(sapply(Brass_nc_tpm[,1], function(x) unlist(strsplit(x, "\\:"))[2]))

         } else if (spec_id == "CR") {

            core_ids <- as.data.frame(sapply(Core_tpm[,1], function(x) unlist(strsplit(x, "\\:"))[3]))
            core_brass_ids <- as.data.frame(sapply(Brass_pc_tpm[,1], function(x) unlist(strsplit(x, "\\:"))[3]))
            core_lnc_ids <- as.data.frame(sapply(Brass_nc_tpm[,1], function(x) unlist(strsplit(x, "\\:"))[3]))

         } else if (spec_id == "ES") {

            core_ids <- as.data.frame(sapply(Core_tpm[,1], function(x) unlist(strsplit(x, "\\:"))[4]))
            core_brass_ids <- as.data.frame(sapply(Brass_pc_tpm[,1], function(x) unlist(strsplit(x, "\\:"))[4]))
            core_lnc_ids <- as.data.frame(sapply(Brass_nc_tpm[,1], function(x) unlist(strsplit(x, "\\:"))[4]))

         } else if (spec_id == "TH") {

            core_ids <- as.data.frame(sapply(Core_tpm[,1], function(x) unlist(strsplit(x, "\\:"))[5]))

         } else if (spec_id == "MT") {

            core_ids <- as.data.frame(sapply(Core_tpm[,1], function(x) unlist(strsplit(x, "\\:"))[6]))

         } else if (spec_id == "BD") {

            core_ids <- as.data.frame(sapply(Core_tpm[,1], function(x) unlist(strsplit(x, "\\:"))[7]))
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
            core_ids <- sub("\\:.*", "", Core_tpm[,1])

            df <- df[rownames(df) %in% core_ids,]

         } else if ((spec_id == "AT") && (c_level == "brass") && (scripttype == "coding")) {

            # Get Brassicaceae protein-coding orthologs

            core_brass_ids <- sub("\\:.*", "", Brass_pc_tpm[,1])

            df <- df[rownames(df) %in% core_brass_ids, ]
         
         } else if ((spec_id == "AT") && (c_level == "brass") && (scripttype == "lncRNA")) {

            # Extract AT lncRNA IDs and get Brassicaceae core lncRNAs

            core_lnc_ids <- sub("\\:.*", "", Brass_nc_tpm[,1])

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

      cd_all <- do.call(rbind, lapply(tpm_table_repl_ls, getNumGenes, scripttype = "coding", 
         c_level = "all"))
      cd_brass <- do.call(rbind, lapply(tpm_table_repl_ls[1:4], getNumGenes, scripttype = "coding", 
         c_level = "brass"))
      cd_core <- do.call(rbind, lapply(tpm_table_repl_ls, getNumGenes, scripttype = "coding", 
         c_level = "core"))
      nc_all <- do.call(rbind, lapply(tpm_table_repl_ls, getNumGenes, scripttype = "lncRNA", 
         c_level = "all"))
      nc_brass <- do.call(rbind, lapply(tpm_table_repl_ls[1:4], getNumGenes, scripttype = "lncRNA", 
         c_level = "core")) # Ortholog lncRNA dataset is limited to Brassicaceae

      at_stats <- rbind(cd_all, cd_brass, cd_core, nc_all, nc_brass)



      # Show message
      message("Writing output...")

      # Set filename
      fname_max_expr <- sprintf('%s.csv', paste(species_id, "max_expr_stats", sep = "_"))

      # Write final data tables to csv files and store them in /out_dir/output/max_expr
      if (!dir.exists(file.path(out_dir, "output", "max_expr"))) 
      dir.create(file.path(out_dir, "output", "max_expr"), recursive = TRUE)

      write.table(at_stats, file = file.path(out_dir, "output", "max_expr", fname_max_expr), 
         sep=";", dec=".", row.names = FALSE, col.names = TRUE)


   }



}



