# Compute coefficient of variation (CV) for 7003 angiosperm core orthologous genes
# Classifiy genes into evolutionarily stable and variable genes based on average CV across organs
# Match stable and variable genes based on their average gene expression level, using caliper matching
# Similar procedure as described in Berthelot et al., Nat Ecol Evol (2018)
# Check proportion of stable/variable genes for GOslim terms of interest  
# GOslim categories retrieved from TAIR version 20211101


getCV <- function(aspect = c("biological_process", "molecular_function"), estimate = c("VST", "TPM"), 
                  pFDR, ...) {

    # Show error message if no/unknown GO aspect is chosen
    if ((missing(aspect)) || (!is.element(aspect, c("biological_process", "molecular_function"))))

        stop("Please choose one of the available aspects: 
            'biological_process', 'molecular_function'",
            call. = TRUE
            )

    # Show error message if no sample_size for GO term size is chosen
    if ((missing(estimate)) || (!is.element(estimate, c("VST", "TPM"))))

        stop("Please choose one of the available expression estimates: 
            'VST', 'TPM'",
            call. = TRUE
            )

    # Show error message if no sample_size for GO term size is chosen
    if ((missing(pFDR)) || (pFDR > 1))

        stop("Please choose p-value cutoff",
            call. = TRUE
            )

    # Set file path for input files
	GOSLIM = file.path(in_dir, "ATH_GO_GOSLIM.txt")
	GOCAT = file.path(in_dir, "TAIR_GO_slim_categories.txt")

    if (estimate == "VST" ) {

        orthoEst = file.path(in_dir, "AT_core_inter_count_mat_vsd_sample_names.csv")

    } else { 

        orthoEst = file.path(in_dir, "AT_core_inter_TPM_mat_deseq_sample_names.csv")
    }

	filenames <- list.files(file.path(out_dir, "output", "data"), pattern="*.txt", full.names=TRUE)
    filenames <- filenames[!filenames %in% "./GOslim/output/data/cor_bsv_traject_1000.txt"]
    names_txt <- gsub('.*/', '', filenames)
    names <- gsub("\\..*","", names_txt)


    ldf <- lapply(filenames, function(x) {
        df <- read.table(x, sep="\t", dec=".", quote = "", header=TRUE, stringsAsFactors=FALSE)
        return(df)
    })

    names(ldf) <- names
    res <- lapply(ldf, summary) # check if file input is OK

    GOSLIM <- read.table(GOSLIM, sep="\t", dec=".", quote = "", header=FALSE, skip=4, fill = TRUE, stringsAsFactors=FALSE)
    GOCAT <- read.table(GOCAT, sep="\t", dec=".", header=TRUE, skip=7, fill = TRUE, stringsAsFactors=FALSE)
    orthoEst <- read.table(orthoEst, sep=";", dec=".", header=TRUE, stringsAsFactors=FALSE)


    # return_list <- list("ldf" = ldf, "orthoEst" = orthoEst, "GOSLIM" = GOSLIM, "GOCAT" = GOCAT, "aspect" = aspect, "estimate" = estimate, "pFDR" = pFDR, "res" = res)
    # return(return_list)
    # }
    # return_objects <- getCV(aspect = "biological_process", estimate = "VST", pFDR = 0.01)
    # list2env(return_objects, envir = .GlobalEnv)

    # Show message
    message("Starting analysis...")

    list2env(ldf, envir = .GlobalEnv)



    #------------------- Calculate average core ortholog expression and CV --------------------


    # Prepare angiosperm ortholog data
    orthoExpr <- data.frame(gene_id=sub("\\:.*", "", orthoEst[,1]),orthoEst[,2:ncol(orthoEst)])
    
    if (estimate == "TPM") { 
        orthoExpr[,2:ncol(orthoExpr)] <- log2(orthoExpr[,2:ncol(orthoExpr)] + 1)
    }

    orthoExpr <- orthoExpr[!grepl("ERCC", orthoExpr$gene_id),]


    # Negate dplyr %in%
    `%!in%` = Negate(`%in%`)

    # Remove pollen samples
    orthoExpr <- orthoExpr %>% select (-c(
        A.thaliana_flowers_mature_pollen_28d_.2.:B.distachyon_flowers_mature_pollen_32d_.1.))


    calculateAvgExpr <- function(df) {

        # Split data frame by sample replicates into a list
        # then get rowMeans for each subset and bind averaged data to gene_id column

            averaged_replicates <- do.call(cbind, lapply(split.default(df[2:ncol(df)], 
                rep(seq_along(df), 
                each = 3, 
                length.out=ncol(df)-1)
                ), rowMeans)
              )

            averaged_replicates <- cbind(df[1], averaged_replicates)
        
            return(averaged_replicates)
    }

    x_avg <- calculateAvgExpr(orthoExpr)

    DevSeq_col_names <- rep(c("Root", "Hypocotyl", "Leaf", "veg_apex", "inf_apex", "Flower", 
    	"Stamen", "Carpel"), each=7)
    DevSeq_spec_names <- rep(c("_AT", "_AL", "_CR", "_ES", "_TH", "_MT", "_BD"), times=8)
    repl_names <- paste0(DevSeq_col_names, DevSeq_spec_names)

    colnames(x_avg)[2:ncol(x_avg)] <- repl_names


    # Compute average expression and CV for each organ
    calculateAvgCV <- function(df) {

            averaged_spec <- do.call(cbind, lapply(split.default(df[2:ncol(df)], 
                rep(seq_along(df), 
                each = 7, 
                length.out=ncol(df)-1)
                ), rowMeans)
              )

            base_averaged <- rowMeans(df[2:ncol(df)])
   

            RowSD <- function(x) {
                sqrt(rowSums((x - rowMeans(x))^2)/(dim(x)[2] - 1))
            }

            averaged_sd <- do.call(cbind, lapply(split.default(df[2:ncol(df)], 
                rep(seq_along(df), 
                each = 7, 
                length.out=ncol(df)-1)
                ), RowSD)
              )

            organCV <- averaged_sd/averaged_spec

            CV_averaged <- rowMeans(organCV)

            names_organ <- unique(sub("\\_.*", "", colnames(df)[2:ncol(df)]))
            avg_names <- paste("CV", names_organ, sep="_")
            colnames(organCV) <- avg_names

            averaged_spec <- cbind(df[1], organCV, CV_averaged, base_averaged)
        
            return(averaged_spec)
    }

    spec_CV <- calculateAvgCV(x_avg)






