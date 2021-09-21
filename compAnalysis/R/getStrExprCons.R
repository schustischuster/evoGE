# Extract GOslim terms with aspect of interest (biological process, molecular function) 
# from TAIR list downloaded on 19th April 2021
# GOSLIM table contains the GOslim terms for 28553 A.thaliana genes with "AT" identifier
# Create GOSLIM lists of 7003 angiosperm orthologous genes

library(dplyr)
library(MatchIt)
library(ggplot2)
library(scales)

getStrExprCons <- function(nquant, ...) {


    # Show error message if no quartile number is chosen
    if ((missing(nquant)) || (nquant < 1))

        stop("Please choose a nquant value greater 1",
            call. = TRUE
            )

	# Set file path for input files
	orthoTPM = file.path(in_dir, "Expression_data", "AT_core_inter_tpm_mat_deseq_sample_names.csv")

	orthoTPM <- read.table(orthoTPM, sep=";", dec=".", header=TRUE, stringsAsFactors=FALSE)


    # return_list <- list("orthoTPM" = orthoTPM, "nquant" = nquant)
    # return(return_list)
    # }
    # return_objects <- getStrExprCons(nquant = 350)
    # list2env(return_objects, envir = .GlobalEnv)

    # Show message
    message("Starting analysis...")



    #------------ Combine DevSeq core ortholog expression tables with GOslim data -------------


    # Prepare angiosperm ortholog data
    orthoExpr <- data.frame(gene_id=sub("\\:.*", "", orthoTPM[,1]),orthoTPM[,2:ncol(orthoTPM)])
    orthoExpr[,2:ncol(orthoExpr)] <- log2(orthoExpr[,2:ncol(orthoExpr)] + 1)
    orthoExpr <- orthoExpr[!grepl("ERCC", orthoExpr$gene_id),]


    # Negate dplyr %in%
    `%!in%` = Negate(`%in%`)

    # Remove pollen samples
    orthoExpr <- orthoExpr %>% select (-c(
        A.thaliana_flowers_mature_pollen_28d_.2.:B.distachyon_flowers_mature_pollen_32d_.1.))

    # Get quartiles
    quant_num <- round(nrow(orthoExpr)/nquant)


    orthoExpr$base_averaged <- rowMeans(orthoExpr[2:ncol(orthoExpr)])
    orthoExpr$quartile <- ntile(orthoExpr$base_averaged, quant_num) # dplyr function to get n quartiles

    quart_ls <- split(orthoExpr, f = orthoExpr$quartile)




    #---------------- Get gene expression divergence rates for ATH/AL vs species X -----------------


    # Show message
    message("Calculating expression distances...")


        # Use pearson correlation, inter-organ normalization and TPM for ms

        getDSOrganCor <- function(df, organ) {

            # Select rows for each organ
            if (organ == "Root"){
                df <- df[,2:22]
            } else if (organ == "Hypocotyl"){
                df <- df[,23:43]
            } else if (organ == "Leaf"){
                df <- df[,44:64]
            } else if (organ == "Apex_veg"){
                df <- df[,65:85]
            } else if (organ == "Apex_inf"){
                df <- df[,86:106]
            } else if (organ == "Flower"){
                df <- df[,107:127]
            } else if (organ == "Stamen"){
                df <- df[,128:148]
            } else if (organ == "Carpel"){
                df <- df[,149:169]
            }

            df_cor <- sqrt(1/2*(1 - cor(df, method="pearson")))

            replcor <- function(xrepl){

                xrepl <- c(xrepl)[xrepl>0][!(duplicated(c(xrepl)[xrepl>0]))]
                return(xrepl)
            }

            sp0_repl <- c(mean(replcor(df_cor[1:3,1:3])), mean(replcor(df_cor[4:6,4:6])), 
                mean(replcor(df_cor[7:9,7:9])), mean(replcor(df_cor[10:12,10:12])), 
                mean(replcor(df_cor[13:15,13:15])), mean(replcor(df_cor[16:18,16:18])), 
                mean(replcor(df_cor[19:21,19:21]))) # AT-AT

            # Get mean
            df_cor <- as.data.frame(df_cor, stringsAsFactors=FALSE)


            avgRepl <- function(x_df) {

                getRepl <- function(x) {

                    split.default(x, 
                        rep(seq_along(x), 
                            each = 3, 
                            length.out=ncol(x)
                            )
                        )
                }

                repl_lst <- getRepl(x_df)
                repl_sum <- lapply(repl_lst, sum)
                repl_mean <- as.numeric(unlist(repl_sum))/9

                return(repl_mean)
            }

            sp1_repl <- avgRepl(df_cor[4:6,1:3]) # AT-AL
            sp2_repl <- avgRepl(df_cor[7:9,1:6]) # AT-CR AL-CR
            sp3_repl <- avgRepl(df_cor[10:12,1:9]) # AT-ES AL-ES CR-ES
            sp4_repl <- avgRepl(df_cor[13:15,1:12]) # AT-TH AL-TH CR-TH ES-TH
            sp5_repl <- avgRepl(df_cor[16:18,1:15]) # AT-MT AL-MT CR-MT ES-MT TH-MT
            sp6_repl <- avgRepl(df_cor[19:21,1:18]) # AT-BD AL-BD CR-BD ES-BD TH-BD MT-BD

            getError <- function(cor_data) {
                std <- sd(cor_data, na.rm=TRUE)
                num <- length(cor_data)
                error <- std/sqrt(num)
                return(error)
            } # Use this function to replace sd if reqired

            df_cor_error <- data.frame(error = c(rep(as.numeric(c(sd(sp0_repl))),length(sp0_repl)),
                    as.numeric(c(sd(sp1_repl))), rep(as.numeric(c(sd(sp2_repl))),length(sp2_repl)), 
                    rep(as.numeric(c(sd(sp3_repl))),length(sp3_repl)), rep(as.numeric(c(sd(sp4_repl))),length(sp4_repl)), 
                    rep(as.numeric(c(sd(sp5_repl))),length(sp5_repl)), rep(as.numeric(c(sd(sp6_repl))),length(sp6_repl))))

            df_cor_avg <- data.frame(correlation = c(sp0_repl, sp1_repl, sp2_repl, sp3_repl, sp4_repl, sp5_repl, sp6_repl))
            div_tag <- data.frame(clade = c(rep("T0", length(sp0_repl)), "T1", rep("T2", length(sp2_repl)), rep("T3", length(sp3_repl)), 
                rep("T4", length(sp4_repl)), rep("T5", length(sp5_repl)), rep("T6", length(sp6_repl))))
            organ_id <- data.frame(comp_organ = rep(organ, nrow(df_cor_avg)))
            div_times <- data.frame(div_times = c(rep(0, length(sp0_repl)), 7.1, rep(9.4, length(sp2_repl)), rep(25.6, length(sp3_repl)), 
                rep(46, length(sp4_repl)), rep(106, length(sp5_repl)), rep(160, length(sp6_repl))))
            dataset <- data.frame(dataset = rep("Angiosperms ", nrow(df_cor_avg)))
            df_cor_avg <- cbind(div_tag, organ_id, div_times, df_cor_avg, df_cor_error, dataset)

            return(df_cor_avg)

        }

        root_divc <- lapply(quart_ls, getDSOrganCor, organ="Root")
        hypocotyl_divc <- lapply(quart_ls, getDSOrganCor, organ="Hypocotyl")
        leaf_divc <- lapply(quart_ls, getDSOrganCor, organ="Leaf")
        veg_apex_divc <- lapply(quart_ls, getDSOrganCor, organ="Apex_veg")
        inf_apex_divc <- lapply(quart_ls, getDSOrganCor, organ="Apex_inf")
        flower_divc <- lapply(quart_ls, getDSOrganCor, organ="Flower")
        stamen_divc <- lapply(quart_ls, getDSOrganCor, organ="Stamen")
        carpel_divc <- lapply(quart_ls, getDSOrganCor, organ="Carpel")


    

        #---- Apply non-linear regression to sOU and pearson dist expression data and compare slopes -----

        # Non-linear regression using negative exponential law fit: pairwise expression differences
        # between species saturate with evolutionary time in a power law relationship
        # Fits assumption of OU model underlying stabilizing GE selection as a decelarated process

        getNLEstimates <- function(corrdata) {

            comp_organ <- unique(corrdata$comp_organ)


            nl_model <- function(a, b, c, x){

                y = a * exp(c * x) + b * (1 - exp(c * x))
                return(y)
            }
            # a + b defines maximum y value
            # a defines intercept


            x_DS_grid <- seq(0, 160, length = 200)  ## prediction grid

            cor_0 <- corrdata$correlation[corrdata$clade=="T0"]

            weights <- c(rep(3.5,7), 0.5, rep(1,2), rep(1.5,3), rep(2,4), rep(2.5,5), rep(3,6))

            # Compute data points for DevSeq_AL_pearson_dist based on model
            # First try to manually find rough parameters, then use nls to fine tune
            mcoeff <- nls(correlation ~ a * exp(div_times * c) + b * (1-(exp(div_times * c))), 
                start = list(a = 0.01, b = 0.5, c = -0.01), data = corrdata, control = list(maxiter = 500), 
                weights=weights)
            coeff <- as.data.frame(summary(mcoeff)["coefficients"])

            model_expr_dist <- data.frame(y = do.call(rbind, lapply(x_DS_grid, nl_model, 
                a = coeff["a",1], b = coeff["b",1], c = coeff["c",1])))

            model_coord <- data.frame(x = x_DS_grid, model_expr_dist)

            # Get slope values
            slopes = diff(model_coord$y)/diff(model_coord$x)
            slopes_avg <- mean(slopes)

            nlm_coord <- data.frame(model_coord, organ=rep(comp_organ, nrow(model_coord)), 
                nlm_slope=rep(slopes_avg, nrow(model_coord)))

            return(nlm_coord)

        }


        # Get mean nlm regression slopes for all goslim control sets
        rtc_nlm_slopes <- data.frame(do.call(rbind, lapply(root_divc, getNLEstimates)))
        hcc_nlm_slopes <- data.frame(do.call(rbind, lapply(hypocotyl_divc, getNLEstimates)))
        lfc_nlm_slopes <- data.frame(do.call(rbind, lapply(leaf_divc, getNLEstimates)))
        avc_nlm_slopes <- data.frame(do.call(rbind, lapply(veg_apex_divc, getNLEstimates)))
        aic_nlm_slopes <- data.frame(do.call(rbind, lapply(inf_apex_divc, getNLEstimates)))
        flc_nlm_slopes <- data.frame(do.call(rbind, lapply(flower_divc, getNLEstimates)))
        stc_nlm_slopes <- data.frame(do.call(rbind, lapply(stamen_divc, getNLEstimates)))
        clc_nlm_slopes <- data.frame(do.call(rbind, lapply(carpel_divc, getNLEstimates)))


        # Check if 1st and last quantile evolve at significant different rate than the rest
        p.qrt1 <- wilcox.test(c(unique(rtc_nlm_slopes[,4])[1], unique(hcc_nlm_slopes[,4])[1], unique(lfc_nlm_slopes[,4])[1],
            unique(avc_nlm_slopes[,4])[1], unique(aic_nlm_slopes[,4])[1], unique(flc_nlm_slopes[,4])[1],
            unique(stc_nlm_slopes[,4])[1], unique(clc_nlm_slopes[,4])[1]), 
            c(unique(rtc_nlm_slopes[,4])[2:19], unique(hcc_nlm_slopes[,4])[2:19], unique(lfc_nlm_slopes[,4])[2:19],
            unique(avc_nlm_slopes[,4])[2:19], unique(aic_nlm_slopes[,4])[2:19], unique(flc_nlm_slopes[,4])[2:19],
            unique(stc_nlm_slopes[,4])[2:19], unique(clc_nlm_slopes[,4])[2:19]))$p.value
        # p = 0.0006853317

        p.qrt20 <- wilcox.test(c(unique(rtc_nlm_slopes[,4])[20], unique(hcc_nlm_slopes[,4])[20], unique(lfc_nlm_slopes[,4])[20],
            unique(avc_nlm_slopes[,4])[20], unique(aic_nlm_slopes[,4])[20], unique(flc_nlm_slopes[,4])[20],
            unique(stc_nlm_slopes[,4])[20], unique(clc_nlm_slopes[,4])[20]), 
            c(unique(rtc_nlm_slopes[,4])[2:19], unique(hcc_nlm_slopes[,4])[2:19], unique(lfc_nlm_slopes[,4])[2:19],
            unique(avc_nlm_slopes[,4])[2:19], unique(aic_nlm_slopes[,4])[2:19], unique(flc_nlm_slopes[,4])[2:19],
            unique(stc_nlm_slopes[,4])[2:19], unique(clc_nlm_slopes[,4])[2:19]))$p.value
        # p = 0.04540067



        # Set up list containing all organ goslim control sets
        control_nlm_slp_lst <- list(rtc_nlm_slopes, hcc_nlm_slopes, lfc_nlm_slopes, 
            avc_nlm_slopes, aic_nlm_slopes, flc_nlm_slopes, stc_nlm_slopes, clc_nlm_slopes)





}
