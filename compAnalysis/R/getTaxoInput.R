# Create data tables that have an input format which is compatible with TreeExp2
# Data tables contain sample replicates with expression values in FPKM (original Brawand 2011 data)
# or TPM expression metric (DevSeq and re-analyzed Brawand data)
# First column contains gene IDs or ortholog numbers



getTaxoInput <- function() {

    DS_table = file.path(in_dir, "Expression_data", "AT_core_inter_tpm_mat_deseq_sample_names.csv")
    Br_table = file.path(in_dir, "Expression_data", "Brawand_inter_tpm_mat_deseq_sample_names_0_5_threshold.csv")
    Br2011_table = file.path(in_dir, "Expression_data", "Brawand_Supplementary_Data1", "NormalizedRPKM_ConstitutiveAlignedExons_Amniote1to1Orthologues.txt")


    # Show startup message
    message("Reading data...")


    # Read DevSeq table
    x_DS_tbj <- read.table(DS_table, sep=";", dec=".", header=TRUE, stringsAsFactors=FALSE)


    # Read Brawand table and set colnames
    x_Br_tbj <- read.table(Br_table, sep=";", dec=".", header=TRUE, stringsAsFactors=FALSE)

    # Remove later in case expression table has gene_id column
    ID_repl_Br <- as.data.frame(seq(1:nrow(x_Br_tbj)))
    colnames(ID_repl_Br) <- "Gene"
    x_Br_tbj <- cbind(ID_repl_Br, x_Br_tbj)


    # Read original Brawand expression table from Nature 2011
    x_Br2011_tbj <- read.table(Br2011_table, sep="\t", dec=".", header=TRUE, stringsAsFactors=FALSE)

    # Remove biological replicates that show log2 RPMK Pearson's r < 0.85
    # remove ortholog names, two hsa samples with low cor and platypus and chicken data 
    x_Br2011_tbj <- x_Br2011_tbj %>% select (-c(hsa:gga, hsa.br.M.4, hsa.br.M.5))
    ID_repl <- as.data.frame(seq(1:nrow(x_Br2011_tbj)))
    colnames(ID_repl) <- "Gene"
    x_Br2011_tbj <- cbind(ID_repl, x_Br2011_tbj)




    #--------------------- Calculate correlation and prepare data for corrplot ---------------------


    # Create "plots" folder in /out_dir/output/plots
    if (!dir.exists(file.path(out_dir, "output", "data"))) 
        dir.create(file.path(out_dir, "output", "data"), recursive = TRUE)

    # Show message
    message("Writing data tables...")


    x_Br_tbj[is.na(x_Br_tbj)] <- 0 # replaces NAs by 0
    x_DS_tbj[is.na(x_DS_tbj)] <- 0 # replaces NAs by 0
    x_Br2011_tbj[is.na(x_Br2011_tbj)] <- 0 # replaces NAs by 0




    #-------------- Create taxa objects for DevSeq and Brawand with replicates  ----------------

    substrRight <- function(x, n){ 
        substr(x, nchar(x)-n+1, nchar(x))
    }


    # Generate taxa object for re-analyzed Brawand mammalian expression data
    x_Br_repl_colnames <- colnames(x_Br_tbj[2:ncol(x_Br_tbj)])
    x_Br_repl_del <- gsub('.{2}$', '', x_Br_repl_colnames)
    x_Br_repl_replname <- substrRight(x_Br_repl_colnames, 1)
    x_Br_taxa_object_names <- c("Gene", paste0(x_Br_repl_del, x_Br_repl_replname))
    colnames(x_Br_tbj) <- x_Br_taxa_object_names


    # Remove Pan troglodytes (Chimp) samples
    x_Br_tbj <- subset(x_Br_tbj, select = -c(Chimp_brain_M4, Chimp_brain_M5, Chimp_brain_M3, 
    	Chimp_brain_M2, Chimp_brain_F1, Chimp_brain_M1, Chimp_cerebellum_F1, Chimp_cerebellum_M1, 
    	Chimp_heart_F1, Chimp_heart_M1, Chimp_kidney_F1, Chimp_kidney_M1, Chimp_liver_F1, 
        Chimp_liver_M1, Chimp_testis_M1))

    write.table(x_Br_tbj, 
    	file=file.path(out_dir, "output", "data", "x_Br_taxobj_input.txt"), sep="\t", 
        col.names=TRUE, row.names=FALSE, dec=".", quote = FALSE)



    # Generate taxa object for original Brawand 2011 mammalian expression data
    x_Br2011_repl_colnames <- colnames(x_Br2011_tbj[2:ncol(x_Br2011_tbj)])
    x_Br2011_repl_taxa <- substr(x_Br2011_repl_colnames, start = 1, stop = 3)
    x_Br2011_repl_subtaxa <- substr(x_Br2011_repl_colnames, start = 5, stop = 6)
    x_Br2011_repl_replname <- substrRight(x_Br2011_repl_colnames, 3)
    x_Br2011_repl_replname <- gsub("[.]","" , x_Br2011_repl_replname ,ignore.case = TRUE)
    x_Br2011_taxa_object_names <- c("Gene", paste(x_Br2011_repl_taxa, x_Br2011_repl_subtaxa, 
        x_Br2011_repl_replname, sep="_"))
    colnames(x_Br2011_tbj) <- x_Br2011_taxa_object_names


    # Remove Pan troglodytes (chimp) samples
    x_Br2011_sel_tbj <- subset(x_Br2011_tbj, select = -c(ptr_br_M1, ptr_br_M3, ptr_br_M5, 
    	ptr_br_M2, ptr_br_M4, ptr_br_F1, ptr_cb_M1, ptr_cb_F1, ptr_ht_M1, ptr_ht_F1, 
        ptr_kd_M1, ptr_kd_F1, ptr_lv_M1, ptr_lv_F1, ptr_ts_M1, oan_br_M1, oan_br_F1, 
        oan_cb_M1, oan_cb_F1, oan_ht_M1, oan_ht_F1, oan_kd_M1, oan_kd_F1, oan_lv_M1, 
        oan_lv_F1, oan_ts_M3, oan_ts_M2, oan_ts_M1, gga_br_M1, gga_br_F1, gga_cb_M1, 
        gga_cb_F1, gga_ht_M1, gga_ht_F1, gga_kd_M1, gga_kd_F1, gga_lv_M1, gga_lv_F1, 
        gga_ts_M1, gga_ts_M2))


    x_Br2011_all_tbj <- subset(x_Br2011_tbj, select = -c(ptr_br_M1, ptr_br_M3, ptr_br_M5, 
        ptr_br_M2, ptr_br_M4, ptr_br_F1, ptr_cb_M1, ptr_cb_F1, ptr_ht_M1, ptr_ht_F1, 
        ptr_kd_M1, ptr_kd_F1, ptr_lv_M1, ptr_lv_F1, ptr_ts_M1))


    write.table(x_Br2011_sel_tbj, 
    	file=file.path(out_dir, "output", "data", "x_Br2011_taxobj_input.txt"), sep="\t", 
        col.names=TRUE, row.names=FALSE, dec=".", quote = FALSE)

    write.table(x_Br2011_all_tbj, 
        file=file.path(out_dir, "output", "data", "x_Br2011_all_taxobj_input.txt"), sep="\t", 
        col.names=TRUE, row.names=FALSE, dec=".", quote = FALSE)



    # Generate taxa object for DevSeq angiosperm gene expression data
    col_namesDS <- rep(c("Root", "Hypocotyl", "Leaf", "vegApex", "infApex", "Flower", "Stamen", 
        "Carpel", "Pollen"), each=21)
    replicate_tag_samples <- rep(c("r1","r2","r3"), times=9)
    col_namesDS <- paste(col_namesDS, replicate_tag_samples, sep="_")
    spec_namesDS <- rep(c("ATH", "ALY", "CRU", "ESA", "THA", "MTR", "BDY"), each=3)
    spec_namesDS <- rep(spec_namesDS, times=9)
    col_namesDS <- paste(spec_namesDS, col_namesDS, sep="_")
    col_namesDS <- c("Gene", col_namesDS)
    colnames(x_DS_tbj) <- col_namesDS

    write.table(x_DS_tbj, 
    	file=file.path(out_dir, "output", "data", "x_DS_taxobj_input.txt"), sep="\t", 
        col.names=TRUE, row.names=FALSE, dec=".", quote = FALSE)

}



