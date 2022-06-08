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

   	expr_table_ls <- list(AT_tpm = pathAT, Core_tpm = pathCore, Brass_nc_tpm = pathNcBrass, AT_tpm_compl = pathAT_compl)
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


   	# Extract all A.thaliana lncRNAs
   	AT_lncRNA_ids <- subset(AT_tpm_compl, subset = biotype %in% c("lnc_exonic_antisense", 
   		"lnc_intronic_antisense", "lnc_intergenic"))[,1]
   	core_lnc_ids <- sub("\\:.*", "", Brass_nc_tpm[,1])


   	# Get number of genes with maximum expression for each organ
   	# "Merge" dev stages for each organ by calculating average number of genes w/ max expression

   	getNumGenes <- function(df, scripttype = c("coding", "lncRNA"), c_level = c("all", "core")) {

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
   	cd_core <- getNumGenes(AT_tpm_repl, scripttype = "coding", c_level = "core")
   	nc_all <- getNumGenes(AT_tpm_repl, scripttype = "lncRNA", c_level = "all")
   	nc_core <- getNumGenes(AT_tpm_repl, scripttype = "lncRNA", c_level = "core")

   	at_stats <- rbind(cd_all, cd_core, nc_all, nc_core)



   	# Show message
   	message("Writing output...")

   	# Set filename
   	fname_max_expr <- sprintf('%s.csv', paste(species_id, "max_expr_stats", sep = "_"))

   	# Write final data tables to csv files and store them in /out_dir/output/data_tables
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

















	all_genes_tpm <- read.table(genesTPM, sep=";", dec=".", header=TRUE, stringsAsFactors=FALSE)
	colnames(all_genes_tpm)[1] <- "gene_id"
	all_genes_counts <- read.table(genesCounts, sep=";", dec=".", header=TRUE, stringsAsFactors=FALSE)
	genes_ptmt <- read.table(genesPtMt, sep=";", dec=".", header=TRUE, stringsAsFactors=FALSE)


	# Format expression data and rename pollen samples
    if ((is.element("ATH", species)) && (is.element("comparative", experiment))) {

		all_genes_tpm <- dplyr::select(all_genes_tpm, c(
			gene_id, biotype, source, info, 
			root_whole_root_5d_.1.,
			root_whole_root_5d_.2.,
			root_whole_root_5d_.3.,
			hypocotyl_10d_.1.,
			hypocotyl_10d_.2.,
			hypocotyl_10d_.3.,
			leaf_1.2_7d_.1.,
			leaf_1.2_7d_.2.,
			leaf_1.2_7d_.3.,
			apex_vegetative_7d_.1.,
			apex_vegetative_7d_.2.,
			apex_vegetative_7d_.3.,
			apex_inflorescence_21d_.1.,
			apex_inflorescence_21d_.2.,
			apex_inflorescence_21d_.3.,
			flower_stg12_21d._.1.,
			flower_stg12_21d._.2.,
			flower_stg12_21d._.3.,
			flower_stg12_stamens_21d._.1.,
			flower_stg12_stamens_21d._.2.,
			flower_stg12_stamens_21d._.3.,
			flowers_mature_pollen_28d_.1.,
			flowers_mature_pollen_28d_.2.,
			flowers_mature_pollen_28d_.3.,
			flower_early_stg12_carpels_21d._.1.,
			flower_early_stg12_carpels_21d._.2.,
			flower_early_stg12_carpels_21d._.3.)) #tibble w/o pollen samples

		species_id <- "ATH_comparative_samples"


    } else if ((is.element("AL", species)) && (is.element("comparative", experiment))) {

		all_genes_tpm <- dplyr::select(all_genes_tpm, -c(
			flower_stg11_stamens_8w.10w.25d_.1., 
			flower_stg11_stamens_8w.10w.25d_.2., 
			flower_stg11_stamens_8w.10w.25d_.3.,
			flower_early_stg12_stamens_8w.10w.23d_.1.,
			flower_early_stg12_stamens_8w.10w.23d_.2.,
			flower_early_stg12_stamens_8w.10w.23d_.3.,
			flower_late_stg12_stamens_8w.10w.21d_.1.,
			flower_late_stg12_stamens_8w.10w.21d_.2.,
			flower_late_stg12_stamens_8w.10w.21d_.3.)) #tibble w/o pollen samples

    }



    # Stop function here to allow specific analysis of a single species
    # return_list <- list("species_id" = species_id, "experiment" = experiment, "all_genes_tpm" = all_genes_tpm, "all_genes_counts" = all_genes_counts, "threshold" = threshold, "genes_ptmt" = genes_ptmt)
    # return(return_list)
    # }
    # return_objects <- getExprGenes(species="ATH", experiment="single-species", 0.05) # read in GTF and expression data for A.thaliana
    # list2env(return_objects, envir = .GlobalEnv)



#--------------------------- Get and apply sample-specific threshold  --------------------------


    # Show message
    message("Starting analysis...")


    # Remove info column
    all_genes_tpm <- dplyr::select(all_genes_tpm, -c(info))
    

	# Extract ERCC data
	ERCC <- all_genes_tpm[all_genes_tpm$gene_id %like% "ERCC", ]


	if (threshold > 0) {

	   # Get sample-specific TPM threshold
	   ERCC_cutoff <- cbind(
		   data.frame(gene_id="TPM_cutoff"), data.frame(biotype="<NA>"), data.frame(source="<NA>"), 
		   as.data.frame(t(data.frame(
			TPM_cutoff = apply(ERCC[,4:ncol(ERCC)], 2, function(i)quantile(i[i>0], threshold)))
		   ))
	   )

	   # Bind sample-specific threshold to expression table
	   all_genes_tpm_cutoff <- rbind(all_genes_tpm, ERCC_cutoff)
	}


	# Set static threshold of 0.5 TPM for "0 ERCC threshold"
	if (threshold == 0) {
		th_values <- data.frame(t(rep(0.5, ncol(ERCC)-3)))
		names(th_values) <- colnames(all_genes_tpm)[4:ncol(all_genes_tpm)]
		ERCC_cutoff <- data.frame(gene_id="TPM_cutoff", biotype="NA", source="NA", th_values)
		all_genes_tpm_cutoff <- rbind(all_genes_tpm, ERCC_cutoff)
	}




	#* test data
	#* df <- data.frame (ID  = c('data1', 'data2', 'data3', 'data4', 'data5', 'data6', 'data7', 'TPM_cutoff'), 
	#* 	biotype  = c('coding', 'coding', 'coding', 'coding', 'coding', 'coding', 'coding', 'NA'), 
	#* 	source  = c('DevSeq', 'DevSeq', 'DevSeq', 'DevSeq', 'DevSeq', 'Araport', 'DevSeq', 'NA'), 
	#*     sample1 = c(2, 7, 1, 18, 3, 0.1, 0, 3),
	#*     sample2 = c(4, 0, 3, 17, 16, 0.3, 0, 5),
	#*     sample3 = c(3, 4, 2, 11, 2, 0.2, 0, 2),
	#*     sample4 = c(22, 14, 9, 11, 35, 0, 0, 11),
	#*     sample5 = c(10, 8, 5, 8, 22, 0, 0, 7),
	#*     sample6 = c(17, 11, 6, 9, 11, 0.1, 0, 8))




	# This is a threshold function that can be applied to expression tables
	# Settings: TPM > ERCC spike-in threshold in at least 2 of 3 replicates in at least one sample
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
		keys_data_repl <- keys_data_repl[-1]

		# Create list of return objects
		return_th_list <- list("express_data_th"=th_df, "express_data_th_replicates"=keys_data_repl)
		return(return_th_list)
	}

    return_objects_th <- applyThreshold(all_genes_tpm_cutoff)
    list2env(return_objects_th, envir = .GlobalEnv)
    #* there are two return objects: 
    #* (1) express_data_th = expression data that only contains genes that show in at least one 
       #* sample type in at least two out of three replicates an expression above the individual 
       #* ERCC spike-in threshold (at defined level of either 0.01/0.05/0.1) or above 0 (if
       #* threshold=0)
    #* (2) express_data_th_replicates = same file as above, but if expression value of replicate
       #* is below the individual ERCC spike-in threshold (or 0), the expression values of the
       #* replicates will be replaced by "0", and if two out of three replicates of one sample
       #* type get "0", the third replicate value also gets "0" (therefore sum of all three
       #* replicates will be zero, which defines a gene that is not expressed [above threshold
       #* level])




#----------- Merge replicates and retrieve number of expressed genes for each sample -----------


    calculateAvgExpr <- function(df) {

	# Split data frame by sample replicates into a list
	# then get rowMeans for each subset and bind averaged data to gene_id/biotype/source column
	
	averaged_replicates <- do.call(cbind, lapply(split.default(df[4:ncol(df)], 
			rep(seq_along(df), 
			each = 3, 
			length.out=ncol(df)-3)
			), rowMeans)
		)

		averaged_replicates <- cbind(df[1:3], averaged_replicates)
		
		return(averaged_replicates)
	}


	express_data_th_avg <- calculateAvgExpr(express_data_th_replicates)

	
	# Set replicate colnames
	repl_names <- unique(substr(names(express_data_th_replicates[-1:-3]), 1, nchar(names(express_data_th_replicates[-1:-3]))-2))
	repl_names <- gsub("^[^.]*.", "", repl_names)
	repl_names <- substr(repl_names, 1, nchar(repl_names)-2)
	repl_names <- gsub('^\\.|\\.$', '', repl_names)

	colnames(express_data_th_avg) <- c("gene_id", "biotype", "source", repl_names)


	# Check which biotypes are in df
	express_data_th_avg$biotype[!(duplicated(express_data_th_avg$biotype))]

	
	# Remove Pt and Mt genes
	express_data_th_avg <- express_data_th_avg[!(express_data_th_avg$gene_id %in% genes_ptmt$gene_id),]


	# Get number of genes expressed in each sample type
	protein_coding_subset <- subset(express_data_th_avg, biotype=="protein_coding")
	lnc_intergenic_subset <- subset(express_data_th_avg, biotype=="lnc_intergenic")
	lnc_antisense_subset <- express_data_th_avg[express_data_th_avg$biotype %like% "antisense", ]
	LTR_subset <- subset(express_data_th_avg, biotype=="LTR_retrotransposon")

	
	# Double-check that no chloroplast and mito genes from coding gene list
	`%nlike%` = Negate(`%like%`)

	if ((is.element("ATH", species_id))) { 

		protein_coding_subset <- protein_coding_subset[protein_coding_subset$gene_id %nlike% "ATMG", ]
		protein_coding_subset <- protein_coding_subset[protein_coding_subset$gene_id %nlike% "ATCG", ]
	}


	expr_protein_coding <- colSums(protein_coding_subset > 0)
	expr_lnc_antisense <- colSums(lnc_antisense_subset > 0)
	expr_lnc_intergenic <- colSums(lnc_intergenic_subset > 0)
	expr_LTR <- colSums(LTR_subset > 0)


	# Create final list of expressed genes per organ/ sample type
	expressed_genes_th_avg <- rbind(expr_protein_coding, expr_lnc_antisense, expr_lnc_intergenic, expr_LTR)
	expressed_genes_th_avg <- subset(expressed_genes_th_avg, select = -c(biotype, source))
	colnames(expressed_genes_th_avg)[1] <- "total_expressed"
	biotype <- substring(rownames(expressed_genes_th_avg), 6)
	expressed_genes_th_avg <- as.data.frame(cbind(biotype, expressed_genes_th_avg))




#--------------- Calculate correlation for biological replicates of each sample ----------------


	# Remove ERCC spike-ins from VST count thresholded expression table
	`%nlike%` = Negate(`%like%`)
	express_data_th_counts <- all_genes_counts[rownames(all_genes_counts) %nlike% "ERCC", ]


	# Get replicate correlation
	replCorr <- function(df, coefficient=c("spearman", "pearson")) {

		replicate_corr <- do.call(cbind, lapply(split.default(df[1:ncol(df)], #adjust columns
								rep(seq_along(df), each = 3, length.out = ncol(df))), #adjust columns
									function(x) {
									repl_corr <- data.frame(cbind(
										cor(x[,1],x[,2],method=coefficient),
										cor(x[,1],x[,3],method=coefficient),
										cor(x[,2],x[,3],method=coefficient)));
									repl_corr
								}
							))

		return(replicate_corr)
	}


	# Get replicate correlations based on normalized VST counts
	repl_corr_df_counts <- replCorr(express_data_th_counts, coefficient="pearson")
	colnames(repl_corr_df_counts) <- colnames(express_data_th_counts)




#-------- Generate replicate expression lists with thresholded genes for each biotype ---------


	# For ERCC-thresholded VST count data
	protein_coding_th_repl_counts <- subset(all_genes_counts, rownames(all_genes_counts) %in% protein_coding_subset$gene_id)
	lnc_antisense_th_repl_counts <- subset(all_genes_counts, rownames(all_genes_counts) %in% lnc_antisense_subset$gene_id)
	lnc_intergenic_th_repl_counts <- subset(all_genes_counts, rownames(all_genes_counts) %in% lnc_intergenic_subset$gene_id)
	LTR_th_repl_counts <- subset(all_genes_counts, rownames(all_genes_counts) %in% LTR_subset$gene_id)
	th_genes_counts <- rbind(protein_coding_th_repl_counts, lnc_antisense_th_repl_counts, lnc_intergenic_th_repl_counts, 
		LTR_th_repl_counts)




#--------------------------------------- Write csv file  ---------------------------------------


	# Show message
    message("Writing output...")


	# Set filename
    fname_expr_genes <- sprintf('%s.csv', paste(species_id, "expr_genes", threshold, sep="_"))
    fname_repl_corr_counts <- sprintf('%s.csv', paste(species_id, "repl_corr_counts", threshold, sep="_"))
    fname_th_genes_repl_counts <- sprintf('%s.csv', paste(species_id, "th_genes_repl_counts", threshold, sep="_"))


	# Write final data tables to csv files and store them in /out_dir/output/data_tables
	if (!dir.exists(file.path(out_dir, "output", "expr_genes"))) 
		dir.create(file.path(out_dir, "output", "expr_genes"), recursive = TRUE)

	write.table(expressed_genes_th_avg, file=file.path(out_dir, "output", "expr_genes", fname_expr_genes), 
		sep=";", dec=".", row.names=FALSE, col.names=TRUE)
	write.table(repl_corr_df_counts, file=file.path(out_dir, "output", "expr_genes", fname_repl_corr_counts), 
		sep=";", dec=".", row.names=FALSE, col.names=TRUE)

	# Only write thresholded expression tables for 0.5 TPM threshold to file
	if (threshold == 0) {
		write.table(th_genes_counts, file=file.path(out_dir, "output", "expr_genes", fname_th_genes_repl_counts), 
			sep=";", dec=".", row.names=TRUE, col.names=TRUE)
	}

}



thresholds <- list(0, 0.01, 0.05, 0.1) # ERCC threshold values are 0 (a fixed TPM threshold of 0.5)
# or perc of expressed spike-ins for 0.01/0.05/0.1

lapply(thresholds, getExprGenes, species = "ATH", experiment = "single-species")
lapply(thresholds, getExprGenes, species = "AL", experiment = "comparative")
lapply(thresholds, getExprGenes, species = "CR", experiment = "comparative")
lapply(thresholds, getExprGenes, species = "ES", experiment = "comparative")
lapply(thresholds, getExprGenes, species = "TH", experiment = "comparative")
lapply(thresholds, getExprGenes, species = "MT", experiment = "comparative")
lapply(thresholds, getExprGenes, species = "BD", experiment = "comparative")


