# Summarize statistics of DevSeq samples
# Data input: mapping statistics of all samples and species


# Install and load packages
if (!require(dplyr)) install.packages('dplyr')
library(dplyr)


# Set file path and input files
in_dir <- "/Volumes/User/Shared/Christoph_manuscript/DevSeq_paper/Analysis/Analysis_2019/A_thaliana_gene_exression_map/20200401_CS_exprGenes/data/Mapping_statistics"
out_dir <- "/Volumes/User/Shared/Christoph_manuscript/DevSeq_paper/Analysis/Analysis_2019/A_thaliana_gene_exression_map/20200401_CS_exprGenes"


# Format data statistic tables
getStats <- function() {

	# Read all csv files in input file path
	readTable <- function(path, pattern = "*.tsv") {
    files = list.files(path, pattern, full.names = TRUE)
    lapply(files, function(x) read.table(x, sep="\t", dec=".", header = TRUE, stringsAsFactors = FALSE))
	}

	stats_tables <- readTable(in_dir)
	stats_table_list <- as.character(list.files(in_dir, pattern = "*.tsv"))
	stats_table_names <- gsub('\\_mapping_stats.tsv$', '', stats_table_list)

	# Change data frame names in list
	names(stats_tables) <- stats_table_names
	list2env(stats_tables, envir = .GlobalEnv)

	# Load sample information
	samples <- read.table(file=file.path(in_dir, "Plant samples for profiling_final_list.csv"), 
	sep=";", dec=".", header=TRUE, stringsAsFactors = FALSE)
	samples_repl <- samples[rep(row.names(samples), samples$Replicates), 1:10]
	comp_sample_names <- as.data.frame(samples_repl[277:303,6])
	names(comp_sample_names) <- "Comparative_Sample"

	# Prepare data tables for plotting
	# ATH_all data
	ATH_data <- A_thaliana[-5]
	names(ATH_data)[1] <- "Sample"
	species_tag = as.data.frame(rep(c("ATH"),each=nrow(ATH_data)))
	names(species_tag) <- "Species"
	ATH_stats <- cbind(species_tag, ATH_data, samples_repl[1:132,c(2,5,7,6)])
	rownames(ATH_stats) <- c()

	# ATH_comp data
	ATH_comp_stats <- subset(ATH_stats, Comparative_Sample!="NA")

	# Non-ATH data
	addSpecTag <- function(species=c("AL", "CR", "ES", "TH", "MT", "BD"), data_stats) {
		
		species_name = as.data.frame(rep(species, each = nrow(data_stats)))
		names(species_name) <- "Species"

		df <- data.frame(species_name, data_stats)
		df <- df[-6]
		names(df)[2] <- "Sample"
		return(df)
	}

	AL <- addSpecTag(species = "AL", data_stats = A_lyrata)
	CR <- addSpecTag(species = "CR", data_stats = C_rubella)
	ES <- addSpecTag(species = "ES", data_stats = E_salsugineum)
	TH <- addSpecTag(species = "TH", data_stats = T_hassleriana)
	MT <- addSpecTag(species = "MT", data_stats = M_truncatula)
	BD <- addSpecTag(species = "BD", data_stats = B_distachyon)

	non_ATH_samples <- rbind(AL, CR, ES, TH, MT, BD)
	non_ATH_stats <- cbind(non_ATH_samples, samples_repl[133:nrow(samples_repl),c(2,5,7,6)])
	rownames(non_ATH_stats) <- c()
	
	# Comparative data
	comp_stats <- rbind(ATH_comp_stats, non_ATH_stats)
	comp_stats <- comp_stats[-c(52:54,55:57,61:63),]
	rownames(comp_stats) <- c()

	# Create list of stat tables
	stat_list <- list(ATH_stats = ATH_stats, non_ATH_stats = non_ATH_stats, comp_stats = comp_stats)

	# Write final data tables to csv files and store them in /out_dir/output/data_tables
	if (!dir.exists(file.path(out_dir, "output", "mapping_statistics"))) 
		dir.create(file.path(out_dir, "output", "mapping_statistics"), recursive = TRUE)

	filepath <- file.path(out_dir, "output", "mapping_statistics", "")
	for (df in names(stat_list)) write.table(stat_list[[df]], file = paste0((filepath), df, ".csv"), 
		sep=";", dec=".", row.names = FALSE, col.names = TRUE)

}

# Execute getStats function
getStats()

