# Compute coefficient of variation (CV) for 7003 angiosperm core orthologous genes
# Classifiy genes into evolutionarily stable and variable genes based on average CV across organs
# Match stable and variable genes based on their average gene expression level, using caliper matching
# Similar procedure as described in Berthelot et al., Nat Ecol Evol (2018)
# Check proportion of stable/variable genes for GOslim terms of interest  
# GOslim categories retrieved from TAIR version 20211101


getCV <- function(aspect = c("biological_process", "molecular_function"), estimate = c("VST", "TPM"), 
    sample_size, ...) {

    # Show error message if no/unknown GO aspect is chosen
    if ((missing(aspect)) || (!is.element(aspect, c("biological_process", "molecular_function"))))

        stop("Please choose one of the available aspects: 
            'biological_process', 'molecular_function'",
            call. = TRUE
            )

    # Show error message if no expression estimate is chosen
    if ((missing(estimate)) || (!is.element(estimate, c("VST", "TPM"))))

        stop("Please choose one of the available expression estimates: 
            'VST', 'TPM'",
            call. = TRUE
            )

    # Show error message if no sample_size for GO term size is chosen
    if ((missing(sample_size)) || (sample_size < 1))

        stop("Please choose one of the available aspects",
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


    # return_list <- list("ldf" = ldf, "orthoEst" = orthoEst, "GOSLIM" = GOSLIM, "GOCAT" = GOCAT, "aspect" = aspect, "estimate" = estimate, "sample_size" = sample_size, "res" = res)
    # return(return_list)
    # }
    # return_objects <- getCV(aspect = "biological_process", estimate = "VST", sample_size = 412)
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


    # Order coeff var table by base expression value
    spec_CV <- spec_CV[order(spec_CV$CV_averaged),]

    spec_CV$sign <- c(rep(1, 3501), rep(0, 3502))


    # Perform nearest distance matching with caliper option
    # Check balance statistics to find optimal caliper
    message("Matching genes...")

    calpr <- 0.5

        matchSample <- function(x) {

            success <- FALSE
            while (!success) {

                # Create background gene set
                match_res <- matchit(sign ~ base_averaged, x, method = "nearest", caliper = calpr, 
                    std.caliper = TRUE, distance = "logit", replace = FALSE, m.order = "data", 
                    ratio = 1)
                match_res_m <- match_res$match.matrix

                # Extract standard mean difference from matchIt summary data
                comp <- as.data.frame(summary(match_res, standardize = TRUE)["sum.matched"])
                stmdif <- abs(comp[1,3])
                varR <- abs(comp[1,4])

                calpr <- calpr-0.01

                # check for success
                success <- ((stmdif <= 0.01) && (varR >= 0.99))
            }

            return(match_res_m)
        }

    match_res_df <- matchSample(spec_CV)

    match_res_sep <- data.frame(treat = rownames(match_res_df), control = match_res_df)
    match_res_sep <- match_res_sep[!is.na(match_res_sep$control),]

    stable_genes <- subset(x_avg, rownames(x_avg) %in% match_res_sep$treat)
    variable_genes <- subset(x_avg, rownames(x_avg) %in% match_res_sep$control)

    stable_genes <- merge(stable_genes, spec_CV)
    variable_genes <- merge(variable_genes, spec_CV)


    message("Process GOslim terms...")







