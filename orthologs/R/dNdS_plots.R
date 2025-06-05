# This script loads and plots dNdS data statistics for angiosperms and mammals
# Input files are either angiosperm or mammalian dNdS data tables


#------------------------------------ Read data tables --------------------------------------

plotdNdS <- function(clade = c("Angiosperms", "Mammals")) {

	# Read all csv files in input file path
	readTable <- function(path, pattern = "*.csv") {

		files = list.files(path, pattern, full.names = TRUE)
		lapply(files, function(x) read.table(x, sep=";", dec=".", header=TRUE, stringsAsFactors=FALSE))
	}

	if (is.element("Angiosperms", clade)) {

		AT_AL_dNdS <- read.table(file.path(in_dir, "map_q=Athaliana.fa_s=Alyrata.fa_orthodetec=RBH_eval=1E-5.csv"), sep=";", dec=".", header=TRUE, stringsAsFactors=FALSE)
		AT_CR_dNdS <- read.table(file.path(in_dir, "map_q=Athaliana.fa_s=Crubella.fa_orthodetec=RBH_eval=1E-5.csv"), sep=";", dec=".", header=TRUE, stringsAsFactors=FALSE)
		AT_ES_dNdS <- read.table(file.path(in_dir, "map_q=Athaliana.fa_s=Esalsugineum.fa_orthodetec=RBH_eval=1E-5.csv"), sep=";", dec=".", header=TRUE, stringsAsFactors=FALSE)
		AT_TH_dNdS <- read.table(file.path(in_dir, "map_q=Athaliana.fa_s=Thalophila.fa_orthodetec=RBH_eval=1E-5.csv"), sep=";", dec=".", header=TRUE, stringsAsFactors=FALSE)
		AT_MT_dNdS <- read.table(file.path(in_dir, "map_q=Athaliana.fa_s=Mtruncatula.fa_orthodetec=RBH_eval=1E-5.csv"), sep=";", dec=".", header=TRUE, stringsAsFactors=FALSE)
		AT_BD_dNdS <- read.table(file.path(in_dir, "map_q=Athaliana.fa_s=Bdistachyon.fa_orthodetec=RBH_eval=1E-5.csv"), sep=";", dec=".", header=TRUE, stringsAsFactors=FALSE)

	} else if(is.element("Mammals", clade)) {

		hs_bonobo_dNdS <- read.table(file.path(in_dir, "map_q=Human.fa_s=Bonobo.fa_orthodetec=RBH_eval=1E-5.csv"), sep=";", dec=".", header=TRUE, stringsAsFactors=FALSE)
		hs_chimp_dNdS <- read.table(file.path(in_dir, "map_q=Human.fa_s=Chimp.fa_orthodetec=RBH_eval=1E-5.csv"), sep=";", dec=".", header=TRUE, stringsAsFactors=FALSE)
		hs_gorilla_dNdS <- read.table(file.path(in_dir, "map_q=Human.fa_s=Gorilla.fa_orthodetec=RBH_eval=1E-5.csv"), sep=";", dec=".", header=TRUE, stringsAsFactors=FALSE)
		hs_macaque_dNdS <- read.table(file.path(in_dir, "map_q=Human.fa_s=Macaque.fa_orthodetec=RBH_eval=1E-5.csv"), sep=";", dec=".", header=TRUE, stringsAsFactors=FALSE)
		hs_mouse_dNdS <- read.table(file.path(in_dir, "map_q=Human.fa_s=Mouse.fa_orthodetec=RBH_eval=1E-5.csv"), sep=";", dec=".", header=TRUE, stringsAsFactors=FALSE)
		hs_oppossum_dNdS <- read.table(file.path(in_dir, "map_q=Human.fa_s=Opossum.fa_orthodetec=RBH_eval=1E-5.csv"), sep=";", dec=".", header=TRUE, stringsAsFactors=FALSE)
		hs_orangutan_dNdS <- read.table(file.path(in_dir, "map_q=Human.fa_s=Orangutan.fa_orthodetec=RBH_eval=1E-5.csv"), sep=";", dec=".", header=TRUE, stringsAsFactors=FALSE)
	}


	# Stop function here to allow specific analysis of a clade
    # if (is.element("Angiosperms", clade)) {
    #    return_list <- list("clade" = clade, "AT_AL_dNdS" = AT_AL_dNdS, "AT_CR_dNdS" = AT_CR_dNdS, "AT_ES_dNdS" = AT_ES_dNdS, "AT_TH_dNdS" = AT_TH_dNdS, "AT_MT_dNdS" = AT_MT_dNdS, "AT_BD_dNdS" = AT_BD_dNdS)
    #    return(return_list)
    # } else if (is.element("Mammals", clade)) {
    #    return_list <- list("clade" = clade, "hs_ppa_dNdS" = hs_bonobo_dNdS, "hs_ptr_dNdS" = hs_chimp_dNdS, "hs_ggo_dNdS" = hs_gorilla_dNdS, "hs_ppy_dNdS" = hs_orangutan_dNdS, "hs_mml_dNdS" = hs_macaque_dNdS, "hs_mmu_dNdS" = hs_mouse_dNdS, "hs_mdo_dNdS" = hs_oppossum_dNdS)
    #    return(return_list)
    # }
    # }
    # return_objects <- plotdNdS(clade = "Angiosperms") # read in dNdS data
    # list2env(return_objects, envir = .GlobalEnv)


    # Create "plots" folder in /out_dir

    if (!dir.exists(file.path(out_dir, "dNdS_plots"))) 
    dir.create(file.path(out_dir, "dNdS_plots"), recursive = TRUE)



    #----------------------------------- Plotting dNdS data -------------------------------------



    # Prepare data for ggplot
    preparedNdS <- function(x) {

    	spec_name <- x[1,2]

    	if (clade == "Angiosperms") {

    		if (substr(spec_name, start = 1, stop = 2) == "AL") {
    			spec <- "AL"

			} else if (substr(spec_name, start = 1, stop = 2) == "Ca") {
				spec <- "CR"

			} else if (substr(spec_name, start = 1, stop = 2) == "Th") {
				spec <- "ES"

			} else if (substr(spec_name, start = 1, stop = 2) == "XM") {
				spec <- "TH"

			} else if (substr(spec_name, start = 1, stop = 2) == "Me") {
				spec <- "MT"

			} else if (substr(spec_name, start = 1, stop = 2) == "Br") {
				spec <- "BD"

		    }

		} else if (clade == "Mammals") {

			if (substr(spec_name, start = 1, stop = 6) == "ENSPTR") {
				spec <- "Ptr"

			} else if (substr(spec_name, start = 1, stop = 6) == "ENSPPA") {
				spec <- "Ppa"

			} else if (substr(spec_name, start = 1, stop = 6) == "ENSGGO") {
				spec <- "Ggo"

			} else if (substr(spec_name, start = 1, stop = 6) == "ENSPPY") {
				spec <- "Ppy"

			} else if (substr(spec_name, start = 1, stop = 6) == "ENSMMU") {
				spec <- "Mml"

			} else if (substr(spec_name, start = 1, stop = 6) == "ENSMUS") {
				spec <- "Mmu"

			} else if (substr(spec_name, start = 1, stop = 6) == "ENSMOD") {
				spec <- "Mdo"

		    }
		}

		x_df <- data.frame(species = rep(spec, nrow(x)), dNdS = x[,"dNdS"])

		return(x_df)
	}


	if (clade == "Angiosperms") {

		dNdS_ls <- list("AT_AL_dNdS" = AT_AL_dNdS, "AT_CR_dNdS" = AT_CR_dNdS, "AT_ES_dNdS" = AT_ES_dNdS, "AT_TH_dNdS" = AT_TH_dNdS, "AT_MT_dNdS" = AT_MT_dNdS, "AT_BD_dNdS" = AT_BD_dNdS)

	} else if (clade == "Mammals") {

		dNdS_ls <- list("hs_ptr_dNdS" = hs_ptr_dNdS, "hs_ggo_dNdS" = hs_ggo_dNdS, "hs_ppy_dNdS" = hs_ppy_dNdS, "hs_mml_dNdS" = hs_mml_dNdS, "hs_mmu_dNdS" = hs_mmu_dNdS, "hs_mdo_dNdS" = hs_mdo_dNdS)
	}

	dNdS_dist <- do.call(rbind, lapply(dNdS_ls, preparedNdS))
	dNdS_dist <- na.omit(dNdS_dist)
	# for mammals, bonobo data (ppa) not shown (same evolutionary distance as chimp)




	# Make replicate correlation plot
	makePlotdNdS <- function(data, clade) {

		fname <- sprintf('%s.pdf', paste(deparse(substitute(data)), clade, sep = "_"))

		data$species <- factor(data$species, levels = unique(data$species))


		if (clade == "Angiosperms") {

			axs_title_y_mar <- margin(t = 0, r = 35, b = 0, l = 1.5)

		} else if (clade == "Mammals") {

			axs_title_y_mar <- margin(t = 0, r = 1.8, b = 0, l = 1.5)

		}

	p <- ggplot(data, aes(x = dNdS, y = species, fill = species, height = after_stat(density))) + 
	     geom_density_ridges(aes(fill = species), scale = 2.0, size = 2.0, alpha = 0.7, stat = "density", trim = TRUE)

	q <- p + scale_fill_manual(values = c("#e8a215", "#f0d737", "#069870", "#0770ab", "#4fb6f0", "#ea6965")) + 
	# Uses a slightly modified colorblind-friendly palette from Wong (Nature Methods, 2011)
	theme_minimal() + 
	xlab("dN/dS") + ylab("Species") + 
	scale_x_continuous(limits = c(0, 2.55), breaks = c(0, 0.5, 1, 1.5, 2, 2.5)) + 
	scale_y_discrete(limits = rev, expand = c(0.05, 0.1)) +
	theme(text = element_text(size = 23.5), 
  		panel.grid.major = element_line(colour = "white"), 
  		panel.grid.minor = element_line(colour = "white"),  
  		axis.ticks.length = unit(.55, "cm"),
  		axis.ticks = element_line(colour = "gray10", linewidth = 2.0),
  		axis.line = element_line(colour = "black", linewidth = 2.0),
  		axis.title.x = element_text(colour = "black", size = 47, 
  			margin = margin(t = 12.8, r = 0, b = 10, l = 0)),  
  		axis.title.y = element_text(colour = "black", size = 47, 
  			margin = axs_title_y_mar), 
  		axis.text.x = element_text(colour = "grey20", margin = margin(t = 8.8, r = 0, b = 1.0, l = 0), size = 43), 
  		axis.text.y = element_text(colour = "grey20", margin = margin(t = 0, r = 5.75, b = 0, l = 4.2), size = 43), 
  		plot.title = element_text(colour = "black", size = 22.85, 
  			margin = margin(t = 36.25, r = 0, b = 15.5, l = 0), hjust = 0.5), 
  		plot.margin = unit(c(17, 2.5, 25, 5), "points"),
		legend.position = "none",
  		panel.border = element_blank())

  	ggsave(file = file.path(out_dir, "dNdS_plots", fname), plot = q,
		scale = 1, width = 10.7, height = 7.3, units = c("in"))
}

makePlotdNdS(data = dNdS_dist, clade)
# for angiosperms, 6 data points above dNdS rate of 2.55 not shown
# for mammals, 610 data points above dNdS rate of 2.55 not shown






}


