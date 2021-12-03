# Estimate the stability of Pearson expression correlations using Monte Carlo simulations 
# Similar procedure as in Schönbrodt and Perugini, 2013


estimatePOS <- function(nbootstrap, coswidth, bss, ...) {


	# Show error message if no sample size for nbootstrap is chosen
	if ((missing(nbootstrap)) || (nbootstrap < 1))

		stop("Please choose number of bootstraps",
			call. = TRUE
			)

    # Show error message if no sample size for coswidth is chosen
    if ((missing(coswidth)) || (coswidth >= 1 | coswidth < 0))

        stop("Please choose corridor width between 0 and 1",
            call. = TRUE
            )

    # Show error message if no sample size for bootstrap support is chosen
    if ((missing(bss)) || (bss  >= 1 | bss < 0))

        stop("Please choose bootstrap support between 0 and 1",
            call. = TRUE
            )


	# Set file path for input files
	orthoexp = file.path(in_dir, "AT_core_inter_tpm_mat_deseq_sample_names.csv")

	orthoexp <- read.table(orthoexp, sep=";", dec=".", header=TRUE, stringsAsFactors=FALSE)


    # return_list <- list("orthoexp" = orthoexp, "nbootstrap" = nbootstrap, "coswidth" = coswidth, "bss" = bss)
    # return(return_list)
    # }
    # return_objects <- estimatePOS(nbootstrap = 1000, coswidth = 0.1, bss = 0.8)
    # list2env(return_objects, envir = .GlobalEnv)

    # Show message
    message("Starting analysis...")



    #------------ Combine DevSeq core ortholog expression tables with GOslim data -------------


    # Prepare angiosperm ortholog data
    orthoexp <- data.frame(gene_id=sub("\\:.*", "", orthoexp[,1]),orthoexp[,2:ncol(orthoexp)])
    orthoexp[,2:ncol(orthoexp)] <- log2(orthoexp[,2:ncol(orthoexp)] + 1)
    orthoexp <- orthoexp[!grepl("ERCC", orthoexp$gene_id),]


    # Negate dplyr %in%
    `%!in%` = Negate(`%in%`)

    # Remove pollen samples
    orthoexp <- orthoexp %>% select (-c(
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

    x_avg <- calculateAvgExpr(orthoexp)

    DevSeq_col_names <- rep(c("Root", "Hypocotyl", "Leaf", "veg_apex", "inf_apex", "Flower", 
            "Stamen", "Carpel"), each=7)
    DevSeq_spec_names <- rep(c("_AT", "_AL", "_CR", "_ES", "_TH", "_MT", "_BD"), times=8)
    repl_names <- paste0(DevSeq_col_names, DevSeq_spec_names)

    colnames(x_avg)[2:ncol(x_avg)] <- repl_names




    # Real correlations
    rr08 <- round(cor(x_avg[,2:ncol(x_avg)])[c("Flower_AL"), c("Flower_ES")], 2) # r=0.8
    rr07 <- round(cor(x_avg[,2:ncol(x_avg)])[c("veg_apex_TH"), c("veg_apex_CR")], 2) # r=0.7
    rr06 <- round(cor(x_avg[,2:ncol(x_avg)])[c("Carpel_BD"), c("Carpel_TH")], 2) # r=0.6
    rr05 <- round(cor(x_avg[,2:ncol(x_avg)])[c("Flower_AT"), c("Hypocotyl_BD")], 2) # r=0.5

    x_avg_sel <- x_avg %>% select(c(Hypocotyl_BD, veg_apex_CR, veg_apex_TH, Flower_AT, Flower_AL, 
        Flower_ES, Carpel_TH, Carpel_BD))



    # Define a list of sample sizes for simulation
    sample_size_ls <- seq(20, 1000, by = 1)


    getCorBsv <- function(i){

        # Sample expression values for n genes
        do.Bts <- function(z){

            data_random <- z[sample(nrow(z), max(sample_size_ls), replace = FALSE), ]

            bsrepl <- sapply(sample_size_ls,
                function(x) data_random[1:x,], simplify = FALSE)

            return(bsrepl)

        }

        bts_samples <- do.Bts(i)

        getCor <- function(df) {

            pea <- cor(df)
            pea <- as.data.frame(pea)

            cor08 <- pea[c("Flower_AL"), c("Flower_ES")]
            cor07 <- pea[c("veg_apex_TH"), c("veg_apex_CR")]
            cor06 <- pea[c("Carpel_BD"), c("Carpel_TH")]
            cor05 <- pea[c("Flower_AT"), c("Hypocotyl_BD")]

            bscor <- data.frame(c08=cor08, c07=cor07, c06=cor06, c05=cor05)
        }

        sampleCor <- data.frame(do.call(cbind, lapply(bts_samples, getCor)))
        cor08_df = sampleCor[,seq(1, ncol(sampleCor), 4)]
        cor07_df = sampleCor[,seq(2, ncol(sampleCor), 4)]
        cor06_df = sampleCor[,seq(3, ncol(sampleCor), 4)]
        cor05_df = sampleCor[,seq(4, ncol(sampleCor), 4)]
        sampleCor <- t(cbind(cor05_df, cor06_df, cor07_df, cor08_df))

        return(sampleCor)
    }


    # Repeat sampling function n times
    set.seed(250) 
    # Note: Seed will give different results in R versions < 3.6.0 because algorithm for 
    # random-number generator underlying sample() function was changed. See post on SO 
    # https://stackoverflow.com/questions/47199415/is-set-seed-consistent-over-different-versions-of-r-and-ubuntu/56381613#56381613
    # and GitHub https://github.com/wch/r-source/blob/7f6cc784523dfa69087958633f7deb309d9c8718/doc/NEWS.Rd#L150-L161:
    cor_bsv_ls <- replicate(nbootstrap, getCorBsv(x_avg_sel), simplify = FALSE)
    cor_bsv <- as.data.frame(do.call("cbind", cor_bsv_ls))




    # Estimate the corridor of stability (COS)
    getCOS <- function(corv){

        # Fisher's z transformation
        z = 0.5*(log(1+corv)-log(1-corv))

        corridor_uw <- z+coswidth
        corridor_lw <- z-coswidth

        # Back-transform upper and lower boundaries to correlation metric (inverse Fisher's z transformation)
        btfz <- function(i) {

            rho <- (exp(2*i) - 1)/(exp(2*i)+1)
            return(rho)
        }

        ui <- btfz(corridor_uw)
        li <- btfz(corridor_lw)

        cdr <- data.frame(cor=corv, ui=ui, li=li)

        return(cdr)
    }

    corv_ls <- c(rr05, rr06, rr07, rr08)

    cos <- data.frame(do.call(rbind, lapply(corv_ls, getCOS)))



    # Check when trajectory stays permanently within COS (POS)
    getPOS <- function(x) {

        max_size <- nrow(x)/4

        c05d <- cor_bsv[1:max_size,]
        c06d <- cor_bsv[(max_size+1):(max_size*2),]
        c07d <- cor_bsv[((max_size*2)+1):(max_size*3),]
        c08d <- cor_bsv[((max_size*3)+1):(max_size*4),]

        getColPOS <- function(df){

            corc <- unique(gsub("\\..*","", rownames(df)))

            if (corc == "c05"){
                cth  <- cos[cos$cor == 0.5, ] 
            } else if (corc == "c06"){
                cth  <- cos[cos$cor == 0.6, ]
            } else if (corc == "c07"){
                cth  <- cos[cos$cor == 0.7, ]
            } else if (corc == "c08"){
                cth  <- cos[cos$cor == 0.79, ]
            }

            value <- suppressWarnings(apply(df, 2, function(i) {

                i <- as.data.frame(i)

                iinv <- i[nrow(i):1,]

                maxv <- min(which(iinv > cth$ui))
                minv <- min(which(iinv < cth$li))

                maxvi <- nrow(i)-maxv+1
                minvi <- nrow(i)-minv+1

                cvalue <- c(abs(maxvi), abs(minvi))
                pos <- max(cvalue[is.finite(cvalue)])

                if (abs(pos)=="Inf") pos <- 0

                return(pos)
            }))

            return(value)
        }

        pos05 <- getColPOS(df=c05d)
        pos06 <- getColPOS(df=c06d)
        pos07 <- getColPOS(df=c07d)
        pos08 <- getColPOS(df=c08d)

        pos_df <- data.frame(pos05=pos05, pos06=pos06, pos07=pos07, pos08=pos08)

        return(pos_df)
    }

    cor_pos <- getPOS(cor_bsv)



    # Plot trajectories
    prepTrajectories <- function(x){

        value <- unlist(x)
        size <- rep(sample_size_ls, times=ncol(x))
        trajectory <- rep(colnames(x), each=nrow(x))

        table <- data.frame(trajectory=trajectory , size=size , value=value)

        return(table)
    }

    maxs_size <- nrow(cor_bsv)/4

    traject05_df <- prepTrajectories(cor_bsv[1:maxs_size,])
    traject06_df <- prepTrajectories(cor_bsv[(maxs_size+1):(maxs_size*2),])
    traject07_df <- prepTrajectories(cor_bsv[((maxs_size*2)+1):(maxs_size*3),])
    traject08_df <- prepTrajectories(cor_bsv[((maxs_size*3)+1):(maxs_size*4),])


    plotTrajectories <- function(df, loc) {

        fname <- sprintf('%s.jpg', paste(deparse(substitute(df)), "COS", sep="_"))

        if (deparse(substitute(df)) == "traject05_df"){
            ui <- cos[cos$cor == 0.5, 2]
            li <- cos[cos$cor == 0.5, 3]
            tc <- 0.5
            th <- quantile(cor_pos$pos05, probs = loc)

        } else if (deparse(substitute(df)) == "traject06_df"){
            ui <- cos[cos$cor == 0.6, 2]
            li <- cos[cos$cor == 0.6, 3]
            tc <- 0.6
            th <- quantile(cor_pos$pos06, probs = loc)

        } else if (deparse(substitute(df)) == "traject07_df"){
            ui <- cos[cos$cor == 0.7, 2]
            li <- cos[cos$cor == 0.7, 3]
            tc <- 0.7
            th <- quantile(cor_pos$pos07, probs = loc)

        } else if (deparse(substitute(df)) == "traject08_df"){
            ui <- cos[cos$cor == 0.79, 2]
            li <- cos[cos$cor == 0.79, 3]
            tc <- 0.79
            th <- quantile(cor_pos$pos08, probs = loc)
        }

        if (deparse(substitute(df)) == "traject06_df" && nbootstrap == 1000 && coswidth == 0.1 && bss == 0.8){

            possegpos <- 0.0585
            th09 <- quantile(cor_pos$pos06, probs = 0.9)
            th095 <- quantile(cor_pos$pos06, probs = 0.95)
            pos09_txt <- c("POS (CL 90%)")
            pos09_ntxt <- paste0("n=", round(quantile(cor_pos$pos06, probs = 0.9)))
            pos09segpos <- 0.253
            y09pos <- 0.253
            y09posntxt <- 0.1765
            pos095_txt <- c("POS (CL 95%)")
            pos095_ntxt <- paste0("n=", round(quantile(cor_pos$pos06, probs = 0.95)))
            pos095segpos <- 0.446
            y095pos <- 0.446
            y095posntxt <- 0.3695
        
        } else {

            possegpos <- -0.2
            th09 <- -0.2
            th095 <- -0.2
            pos09_txt <- ""
            pos09_ntxt <- ""
            pos09segpos <- -0.2
            y09pos <- -0.2
            y09posntxt <- -0.2
            pos095_txt <- ""
            pos095_ntxt <- ""
            pos095segpos <- -0.2
            y095pos <- -0.2
            y095posntxt <- -0.2
        }

        pos_txt <- c("POS (CL 80%)")
        pos_ntxt <- paste0("n=", round(quantile(cor_pos$pos06, probs = loc)))
        cos_txt <- paste0("(w=", coswidth, ")")
        cor_txt <- paste0("actual correlation (r=", tc, ")")

        p <- ggplot(df, aes(size, value, group = trajectory)) +
        geom_line(data = df, aes(size, value),
            color = "#e26887", alpha = 0.50, size = 0.25) + 
        scale_x_continuous(limits = c(0, 955), breaks = c(0, 200, 400, 600, 800), expand = c(0.007, 0)) + 
        scale_y_continuous(limits = c(-0.15, 0.99), breaks = c(0, 0.2, 0.4, 0.6, 0.8), expand = c(0, 0)) + 
        geom_hline(yintercept = ui, linetype = 2, size = 1.145) + 
        geom_hline(yintercept = li, linetype = 2, size = 1.145) + 
        geom_hline(yintercept = tc, size = 1.115) + 
        geom_vline(xintercept = th, col="grey35", size = 1.15) + 
        annotate("text", x=557, y=0.0585, label=pos_txt, size=7.75, col="grey35") + 
        annotate("text", x=504.0, y=-0.018, label=pos_ntxt, size=7.75, col="grey35") +
        geom_segment(aes(x = 454, y = possegpos, xend = 424, yend = possegpos), arrow = arrow(length = unit(0.5, "cm")), 
            size=1.1, col="grey35") + 
        annotate("text", x=710, y=y09pos, label=pos09_txt, size=7.75, col="grey50") + 
        annotate("text", x=656.9, y=y09posntxt, label=pos09_ntxt, size=7.75, col="grey50") +
        geom_segment(aes(x = 606.9, y = pos09segpos, xend = 576.9, yend = pos09segpos), arrow = arrow(length = unit(0.5, "cm")), 
            size=1.1, col="grey50") + 
        annotate("text", x=853.5, y=y095pos, label=pos095_txt, size=7.75, col="grey65") + 
        annotate("text", x=800.4, y=y095posntxt, label=pos095_ntxt, size=7.75, col="grey65") +
        geom_segment(aes(x = 750.5, y = pos095segpos, xend = 720.5, yend = pos095segpos), arrow = arrow(length = unit(0.5, "cm")), 
            size=1.1, col="grey65") +
        geom_segment(aes(x=71.05, xend=83.25, y=0.045, yend=0.045), colour = "black", show.legend = FALSE, size = 1.1) + 
        geom_segment(aes(x=95.3, xend=107.5, y=0.045, yend=0.045), colour = "black", show.legend = FALSE, size = 1.1) + 
        geom_segment(aes(x=71.05, xend=83.25, y=0.0715, yend=0.0715), colour = "black", show.legend = FALSE, size = 1.1) + 
        geom_segment(aes(x=95.3, xend=107.5, y=0.0715, yend=0.0715), colour = "black", show.legend = FALSE, size = 1.1) + 
        geom_segment(aes(x=th09, xend=th09, y=0.205, yend=0.945), colour = "grey50", show.legend = FALSE, size = 1.15) + 
        geom_segment(aes(x=th095, xend=th095, y=0.4, yend=0.8), colour = "grey65", show.legend = FALSE, size = 1.15) + 
        annotate("text", x=239, y=0.0585, label= "Corridor of Stability", size=7.75) + 
        annotate("text", x=167, y=-0.018, label= cos_txt, size=7.75) + 
        annotate("text", x=792.75, y=0.92, label= cor_txt, size=7.75) + 
        geom_segment(aes(x = 792.75, y = 0.875, xend = 792.75, yend = 0.6125), arrow = arrow(length = unit(0.5, "cm")), 
            size=1.1, col="black") +
        labs(x = "Sample size", y = "Correlation") +
        theme(panel.background = element_blank(), 
            axis.ticks.length = unit(0.29, "cm"), 
            axis.ticks = element_line(colour = "black", size = 1.25), 
            axis.line = element_line(colour = 'black', size = 1.25), 
            plot.margin = unit(c(0.2, 0.1, 0, 0),"cm"), 
            axis.title.y = element_text(size=22.75, margin = margin(t = 0, r = 7.0, b = 0, l = 10), 
                colour="black", face = "bold"), 
            axis.title.x = element_text(size=22.75, margin = margin(t = 4.0, r = 0, b = 7.15, l = 0), 
                colour="black", face = "bold"), 
            axis.text.x = element_text(size=18.8, margin = margin(t = 2.5, b = 8), colour="grey20"), 
            axis.text.y = element_text(size=18.8, angle=0, margin = margin(l = 4.1, r = 1.5), colour="grey20")
        )

        suppressWarnings(ggsave(file = file.path(out_dir, "output", "plots", fname), plot = p, 
               width = 11.5, height = 6.5, dpi = 300, units = c("in"), limitsize = FALSE))

    }

    plotTrajectories(traject06_df, loc = bss)



    # Write data to file
    # Show message
    message("Writing data tables...")

    # Create "data" folder in /out_dir/output
    if (!dir.exists(file.path(out_dir, "output", "data"))) 
        dir.create(file.path(out_dir, "output", "data"), recursive = TRUE)

    mcs_out_list <- list(cor_bsv_traject_ = cor_bsv)

    for(i in names(mcs_out_list)){
        write.table(mcs_out_list[[i]], file=file.path(out_dir, "output", "data", paste0(i, nbootstrap, ".txt")), 
            sep="\t", col.names=TRUE, row.names=TRUE, dec=".", quote = FALSE)
    }


}

estimatePOS(nbootstrap = 1000, coswidth = 0.1, bss = 0.8)



