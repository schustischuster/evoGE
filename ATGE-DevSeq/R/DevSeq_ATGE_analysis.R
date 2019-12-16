
# Perform a comparative analysis between the AtGenExpress Development genome array and the Arabidopsis
# DevSeq RNA-Seq data sets. The final DevSeq/ATGE expression tables generated are used as database
# input for the DevSeq-ATGE comparative analysis application @ www.devseqplant.org


#------------------- Load packages, download files and read sample tables ---------------------


# Install and load the following R packages
if (!require(dplyr)) install.packages('dplyr')
library(dplyr)
if (!require(gplots)) install.packages('gplots')
library(gplots)
if (!require(factoextra)) install.packages('factoextra')
library(factoextra)
if (!require(dendextend)) install.packages('dendextend')
library(dendextend)


# Create 'DevSeq-ATGE' and 'data' folders in working directory
dir.create(file.path("DevSeq-ATGE","data"), recursive = TRUE)
setwd(file.path("DevSeq-ATGE", "data"))


# Download csv files containing DevSeq and ATGE expression data
download.file("https://raw.githubusercontent.com/schustischuster/DevSeq-ATGE/master/data/AtGE_dev_gcRMA_linear.csv?token=AFTF4YHGWCBMZL3JT4X3BCS5YVUFW", "ATGE_input_file.csv")
download.file("https://raw.githubusercontent.com/schustischuster/DevSeq-ATGE/master/data/No_TE_genes_tpm_sample_ids_w_gene_symbol_20190626.csv?token=AFTF4YHKKL6PU3IEPTWCCDK5YVUF4", "DevSeq_input_file.csv")

# Download sample description, replicate information and Araport symbol list
download.file("https://raw.githubusercontent.com/schustischuster/DevSeq-ATGE/master/data/Comp_sample_ids_descriptions.csv?token=AFTF4YBVW64NXXPWPOBYMN25YVUGA", "Samples_input_file.csv")
download.file("https://raw.githubusercontent.com/schustischuster/DevSeq-ATGE/master/data/Comp_sample_ids_replicates.csv?token=AFTF4YEEZX553F63UKMVPQC5YVUGG", "Samples_replicates_input_file.csv")
download.file("https://raw.githubusercontent.com/schustischuster/DevSeq-ATGE/master/data/Gene_IDs_ATH_names_wo_dupl.csv?token=AFTF4YFTRFPVAV4J3UTWQ325YVUGK", "Gene_names_input_file.csv")


# Read data tables
message("Reading data tables: input files, sample information, symbol list")
atge_all_samples <- read.table(file="ATGE_input_file.csv", sep=";", dec=".", header=TRUE, stringsAsFactors = FALSE)
devseq_all_samples <- read.table(file="DevSeq_input_file.csv", sep=";", dec=".", header=TRUE, stringsAsFactors = FALSE)
comparative_samples <- read.table(file="Samples_input_file.csv", sep=";", dec=".", header=TRUE, stringsAsFactors = FALSE)
replicate_samples <- read.table(file="Samples_replicates_input_file.csv", sep=";", dec=".", header=TRUE, stringsAsFactors = FALSE)
gene_names <- read.table(file="Gene_names_input_file.csv", sep=";", dec=".", header=TRUE, stringsAsFactors = FALSE)




#--------------- Select a subset of expression data for comparative analysis -----------------

message("Starting analysis...")

# Prepare replicate list and extract sample names 
replicate_samples <- rbind(c("id", "AGI.code", "NA"), replicate_samples)
names(replicate_samples) <- NULL

naOmitList <- function(y) { 
		   return(y[!sapply(y, function(x) all(is.na(x)))]) 
		}
		
devseq_replicates  <- naOmitList(replicate_samples[,1])
atge_replicates <- replicate_samples[,2]


# Select comparative samples from input tables and arrange them in specific order
devseq_samples_replicates <- devseq_all_samples[, which(names(devseq_all_samples) %in% devseq_replicates)]
devseq_samples_replicates <- devseq_samples_replicates[devseq_replicates]

atge_samples_replicates <- atge_all_samples[, which(names(atge_all_samples) %in% atge_replicates)]
atge_samples_replicates <- atge_samples_replicates[atge_replicates]

colnames(atge_samples_replicates)[1] <- "gene_id"
		    



#--------------------------- Apply TPM threshold to DevSeq samples ----------------------------


# This is a threshold function that can be applied to expression tables
applyThreshold <- function(df, threshold) {
  
  	#* Add an error if radius < 0
  	if (threshold < 0)
    	stop(
        "'threshold' must be >= 0",
	   	call. = TRUE
    	)

    # Apply an average threshold function (average TPM > 0.5)
	df_average_threshold <- df[which(rowMeans(df[,-1, drop = FALSE]) > 0.5),]

	# Add keys to data frame
	key <- seq(1, nrow(df), 1)
	df <- cbind(as.data.frame(key),df)

	# Define replicate threshold function
	getThreshold <- function(df) {

		# Split data frame by sample replicates into a list then apply threshold for each subset
	
		th_replicates <- do.call(cbind, lapply(split.default(df[3:ncol(df)], #adjust columns
							rep(seq_along(df), each = 3, length.out = ncol(df)-2)), #adjust columns
							function(x) {
								x[rowSums(x > threshold) < 2, ] <- 0; 
								x
							}
						))

		# Bind key/id columns to thresholded data frame
		th_replicates <- cbind(df[1:2], th_replicates)

		# Remove all rows that only contain "0"
		th_replicates <- th_replicates[which(rowSums(th_replicates[,-1:-2, drop = FALSE] > 0) > 0),]

		return(th_replicates)
	}

	# Apply threshold to data and extract keys ("key")
	keys_data <- getThreshold(df)
	keys_data <- keys_data[,1:2]
	names(keys_data) <- c("key","ID")

	# Generate thresholded data frame based on keys
	th_df <- merge(keys_data, df, by="key")
	th_df <- th_df[order(th_df$key),]
	th_df <- th_df[-1:-2]

	# Combine average and replicate threshold data
	th_df <- rbind(th_df,
	subset(df_average_threshold, !(id %in% th_df$id)))

	return(th_df)
}


# Apply threshold (TPM > 2 in two of three replicates OR average expression > 0.5) to devseq data 
devseq_replicates_threshold <- applyThreshold(devseq_samples_replicates, 2)




#--------------------- Calculate average expression and remove replicates ---------------------


calculateAvgExpr <- function(df) {

	# Split data frame by sample replicates into a list
	# then get rowMeans for each subset, simplify output and bind to gene_id column
	
	averaged_replicates <- data.frame(df[1],

		sapply(split.default(df[2:ncol(df)], 
			rep(seq_along(df), 
			each = 3, 
			length.out=ncol(df)-1)
			), rowMeans)
		)
		
	return(averaged_replicates)
}

# ---- 1.ATGE
# Merge replicates and add column names
atge_samples <- calculateAvgExpr(atge_samples_replicates)

colnames(atge_samples) <- unique(gsub("_[A-C]", "", names(atge_samples_replicates)))

# Merge and remove ATGE_98 and ATGE_99 replicates and rearrange columns
atge_samples <- transform(atge_samples, ATGE_98_99 = rowMeans(atge_samples[,c("ATGE_98", "ATGE_99")]))
atge_samples <- subset(atge_samples, select=-c(ATGE_98, ATGE_99))
atge_samples = atge_samples %>% select(gene_id, ATGE_3, ATGE_9, ATGE_98_99, everything())

# ---- 2.DevSeq
# Merge replicates and add column names
devseq_samples <- calculateAvgExpr(devseq_replicates_threshold)

colnames(devseq_samples) <- c("gene_id", comparative_samples[,1])

# Filter for genes present in both data sets
devseq_samples_w_atge_ids <- devseq_samples %>% filter(gene_id %in% atge_samples$gene_id)
atge_samples_w_devseq_ids <- atge_samples %>% filter(gene_id %in% devseq_samples$gene_id)

# Remove input files from memory
rm(atge_all_samples, devseq_all_samples)




#---------------------- Calculate Relative Expression (RE) and log_RE ------------------------


# The following wrapper function calculates relative expression (RE) and logRE
# it scales all expression values to the range between 0 and 1

getRE <- function(x, scale = c("none", "log2")) { 

    # Show error message if no scaling is chosen
    if (missing(scale))
   
       stop(
           "Please choose one of the available scalings: 
	   'none', 'log2'",
	   call. = TRUE
        )

    reWrapper <- function(x) {
    
        data.frame(x[1], 
	    as.data.frame(t(apply(x[,c(2:ncol(x))], 1, # transposes matrix output from apply function
	    
	        normalize <- function(x) {
		
		    # log2 scaling function
		    if (is.element(scale, c("log2"))) {
		        x <- log2(x + 1)
		    } 
		    
		    # calculate relative expression
		    calculateRE <- function(x) {
		        (x - min(x, na.rm = TRUE)) / 
  			(max(x, na.rm = TRUE) - min(x, na.rm = TRUE))
		    }
		    
		    x <- calculateRE(x)
		    return(x)
                }
            )))
	)
    }
    x <- reWrapper(x)

    if (is.element(scale, c("none"))) {

        colnames(x)[-(1)] <- paste0(colnames(x)[-(1)], "_RE")
	return(x)

    } else if (is.element(scale, c("log2"))) {

	colnames(x)[-(1)] <- paste0(colnames(x)[-(1)], "_log_RE")
	return(x)
    }	
}

# Apply getRE function
atge_re <- getRE(atge_samples_w_devseq_ids, scale = "none")
atge_log2_re <- getRE(atge_samples_w_devseq_ids, scale = "log2")
devseq_re <- getRE(devseq_samples_w_atge_ids, scale = "none")
devseq_log2_re <- getRE(devseq_samples_w_atge_ids, scale = "log2")




#--------------------------- Add gene symbols to expression tables ----------------------------


# Create a table of gene_ids and corresponding gene symbols
all_atge_devseq_genes <- as.data.frame(devseq_re[, 1])
names(all_atge_devseq_genes) <- c("gene_id")

gene_names_id <- as.data.frame(gene_names[, 1])
names(gene_names_id) <- c("gene_id")

atge_devseq_genes_wo_symbol <- anti_join(all_atge_devseq_genes, gene_names_id, by = "gene_id")
atge_devseq_genes_wo_symbol$symbol = atge_devseq_genes_wo_symbol$gene_id
gene_names_all_genes <- rbind(gene_names, atge_devseq_genes_wo_symbol)

# Add gene symbols
re_df_list <- list(atge_re = atge_re, atge_log2_re = atge_log2_re, devseq_re = devseq_re,
		   devseq_log2_re = devseq_log2_re)

re_sbl_list <- lapply(re_df_list, function(df) {
		    df <- merge(df, gene_names_all_genes[, c("gene_id", "symbol")], by="gene_id")
		    df = df %>% select(gene_id, symbol, everything())
		    df <- df[order(df$gene_id),]
		})

list2env(re_sbl_list, envir = .GlobalEnv)




#------------------------------------ Compute correlations -------------------------------------


getCor <- function(df1, df2) {

	# startup message
	message("Computing correlation...")

	df1$Spearman <- sapply(1:nrow(df1), function(i) 
	    cor(as.numeric(df1[i, 3:28]), as.numeric(df2[i, 3:28]), method=c("spearman")))

	df1$Pearson <- sapply(1:nrow(df1), function(i) 
	    cor(as.numeric(df1[i, 3:28]), as.numeric(df2[i, 3:28]), method=c("pearson")))

	return(df1)
    }

atge_re <- getCor(atge_re, devseq_re)
atge_log2_re <- getCor(atge_log2_re, devseq_log2_re)
devseq_re <- getCor(devseq_re, atge_re)
devseq_log2_re <- getCor(devseq_log2_re, atge_log2_re)




#-------------- Generate the final expression tables and save them to "data_tables" folder -------------


# Combine ATGE_RE/DevSeq_RE and ATGE_log2_RE/DevSeq_log2_RE tables and order by gene_id
# Change header of data frames
# Note: ATGE 2nd_internode is 1st_internode based on global expression profile analysis,
# therefore it is re-named 1st_internode_24d in the description of the final expression tables 

descriptions <- c(
		"gene_id", 
		"symbol",
		"gene_id_exp",
		"Spearman",
		"Pearson", 
		"root_whole_root_7d",
		"root_whole_root_14d.17d",
		"root_whole_root_21d",
		"hypocotyl_10d.7d",
		"1st_internode_24d.28d",
		"cotyledons_7d", 
		"leaf_1_2_7d",
		"leaf_1_2_petiole_10d.17d",
		"leaf_5_6_17d", 
		"leaves_senescing_35d",
		"apex_vegetative_7d", 
		"apex_vegetative_14d",
		"apex_inflorescence_21d",
		"flower_stg9_21d+", 
		"flower_stg10_11_21d+",
		"flower_stg12_21d+",
		"flower_stg15_21d+", 
		"flower_stg12_sepals_21d+",
		"flower_stg15_sepals_21d+",
		"flower_stg12_petals_21d+", 
		"flower_stg15_petals_21d+",
		"flower_stg12_stamens_21d+",
		"flower_stg15_stamens_21d+",
		"flower_stg12_carpels_21d+",
		"flower_stg15_carpels_21d+",
		"fruit_stg16_siliques_28d+"
		)


re_cor_list <- list(atge_re = atge_re, atge_log2_re = atge_log2_re, 
		    devseq_re = devseq_re, devseq_log2_re = devseq_log2_re)

# Add experimental tags

addExpTag <- function(x) {

    exp_tag <- sub("_.*", "", colnames(x[3])) # extract ATGE/DevSeq from colname
    exp_tag <- paste0("_", exp_tag) # add linker
    x$gene_id_exp <- x$gene_id # duplicate gene_id column
    x$gene_id_exp <- paste0(x$gene_id_exp, exp_tag) # attach tag to gene_id
    return(x)
}

re_cortg_list <- lapply(re_cor_list, addExpTag)

# Reorder data frames and add descriptions
re_fin_list <- lapply(re_cortg_list, function(df) { 
			df = df %>% select(
			gene_id, 
			symbol, 
			gene_id_exp, 
			Spearman, 
			Pearson, 
			everything())
			colnames(df) <- descriptions
			return(df)
		    })

list2env(re_fin_list, envir = .GlobalEnv)

# Combine data frames and order by gene_id
devseq_re_vs_atge_re <- rbind(devseq_re, atge_re)
devseq_re_vs_atge_re <- devseq_re_vs_atge_re[order(devseq_re_vs_atge_re$gene_id),]

devseq_log2_re_vs_atge_log2_re <- rbind(devseq_log2_re, atge_log2_re)
devseq_log2_re_vs_atge_log2_re <- devseq_log2_re_vs_atge_log2_re[order(devseq_log2_re_vs_atge_log2_re$gene_id),]

# Write final data tables to csv files
setwd("..")
dir.create(file.path("output","data_tables"), recursive = TRUE)
setwd(file.path("output", "data_tables"))
message("Storing results in: ", file.path("output", "data_tables"))
write.table(devseq_re_vs_atge_re, file="DevSeq_RE_vs_ATGE_RE.csv", sep=";", dec=".", row.names=FALSE, col.names=TRUE)
write.table(devseq_log2_re_vs_atge_log2_re, file="DevSeq_log2_RE_vs_ATGE_log2_RE.csv", sep=";", dec=".", row.names=FALSE, col.names=TRUE)


